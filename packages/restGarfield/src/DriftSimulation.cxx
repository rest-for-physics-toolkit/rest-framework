#include "DriftSimulation.h"

#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "FieldMap.h"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/MediumGas.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "TFile.h"
#include "TRestDetectorDriftResponse.h"
#include "TRestDetectorMedium.h"
#include "TRestDetectorReadout.h"
#include "TRestDetectorReadoutManager.h"
#include "TRestDetectorTransportManager.h"
#include "TRestVersion.h"
#include "TRestYAMLUtils.h"
#include "TransportSimulation.h"

using ROOT::Math::XYZVector;

namespace {
// Garfield units: cm, ns. REST units: mm, us
constexpr double kCmPerMm = 0.1;
constexpr double kUsPerNs = 1e-3;
/// Distance (mm) of an end point to the boundary plane to count it as arrived (Garfield ends the drift at
/// the edge of the area within 1e-8 cm)
constexpr double kArrivalTolerance = 1e-3;

double Component(const XYZVector& vector, int axis) {
    return axis == 0 ? vector.X() : (axis == 1 ? vector.Y() : vector.Z());
}

/// Limits (cm) of the mesh of a field map. ComponentGrid::GetMesh gives the numbers of nodes as std::size_t
/// in the versions after Garfield++ 4.0 (e.g. that of the LCG views) and as unsigned int in 4.0: the first
/// overload exists only for the former
template <typename Grid>
auto GetMeshLimits(const Grid& grid, double* low, double* high, int)
    -> decltype(grid.GetMesh(std::declval<std::size_t&>(), std::declval<std::size_t&>(),
                             std::declval<std::size_t&>(), low[0], high[0], low[1], high[1], low[2], high[2]),
                void()) {
    std::size_t nx = 0, ny = 0, nz = 0;
    grid.GetMesh(nx, ny, nz, low[0], high[0], low[1], high[1], low[2], high[2]);
}

template <typename Grid>
void GetMeshLimits(const Grid& grid, double* low, double* high, long) {
    unsigned int nx = 0, ny = 0, nz = 0;
    grid.GetMesh(nx, ny, nz, low[0], high[0], low[1], high[1], low[2], high[2]);
}

/// Writes the response of a node to a pipe: primaries, arrived, ranges and counts
bool WriteNode(FILE* pipe, const TRestDetectorDriftResponse::Node& node) {
    const int counts[2] = {node.fPrimaries, node.fArrived};
    const double range[6] = {node.fUMin, node.fUMax, node.fVMin, node.fVMax, node.fTMin, node.fTMax};
    return fwrite(counts, sizeof(int), 2, pipe) == 2 && fwrite(range, sizeof(double), 6, pipe) == 6 &&
           fwrite(node.fPositionCounts.data(), sizeof(int), node.fPositionCounts.size(), pipe) ==
               node.fPositionCounts.size() &&
           fwrite(node.fTimeCounts.data(), sizeof(int), node.fTimeCounts.size(), pipe) ==
               node.fTimeCounts.size();
}

/// Reads the response of a node written by WriteNode
bool ReadNode(FILE* pipe, TRestDetectorDriftResponse::Node& node, size_t positionSize, size_t timeSize) {
    int counts[2];
    double range[6];
    if (fread(counts, sizeof(int), 2, pipe) != 2 || fread(range, sizeof(double), 6, pipe) != 6) return false;
    node.fPrimaries = counts[0];
    node.fArrived = counts[1];
    node.fUMin = range[0];
    node.fUMax = range[1];
    node.fVMin = range[2];
    node.fVMax = range[3];
    node.fTMin = range[4];
    node.fTMax = range[5];
    node.fPositionCounts.resize(positionSize);
    node.fTimeCounts.resize(timeSize);
    return fread(node.fPositionCounts.data(), sizeof(int), positionSize, pipe) == positionSize &&
           fread(node.fTimeCounts.data(), sizeof(int), timeSize, pipe) == timeSize;
}
}  // namespace

void SetGridOnReadout(TRestDetectorDriftResponse& response, const TRestDetectorReadout& readout) {
    if (response.HasRange(0) && response.HasRange(1) && response.HasRange(2)) return;
    const std::string error = "restGarfield: drift response " + response.fResponseFile + ": ";
    const auto [low, high] = readout.GetChannelsBoundingBox();
    if (std::isnan(low.X()) || std::isnan(high.X())) {
        throw std::runtime_error(error + "the readout has no channels: give xRange, yRange and zRange");
    }
    const XYZVector normal = readout.GetReadoutNormal((low + high) / 2);
    if (normal.R() == 0) throw std::runtime_error(error + "the readout has no normal");
    for (int axis = 0; axis < 3; axis++) {
        if (response.HasRange(axis)) continue;
        if (std::abs(Component(normal, axis)) > 1e-6) {
            throw std::runtime_error(error + "the range along the readout normal (" + "xyz"[axis] +
                                     "Range) must be given, and the normal must be along x, y or z");
        }
        response.SetRange(axis, Component(low, axis) - response.fMargin,
                          Component(high, axis) + response.fMargin);
    }
}

void ComputeDriftResponse(const TRestDetectorMedium& medium, const TRestDetectorTransportTable& table,
                          const TRestDetectorReadout& readout, TRestDetectorDriftResponse& response,
                          bool verbose, int jobs) {
    const std::string error = "restGarfield: drift response " + response.fResponseFile + ": ";
    if (!response.HasRange(0) || !response.HasRange(1) || !response.HasRange(2)) {
        throw std::runtime_error(error + "the ranges of the grid are not set (SetGridOnReadout)");
    }
    auto gas = CreateGasFromTable(medium, table, medium.GetPressure());
    const std::array<double, 3> offset = {response.fFieldMapOffset[0], response.fFieldMapOffset[1],
                                          response.fFieldMapOffset[2]};
    auto map = LoadCSTFieldMap(response.fFieldMap, offset);
    map->SetMedium(gas.get());
    double low[3], high[3];  // cm
    GetMeshLimits(*map, low, high, 0);

    // The readout normal at the centre of the grid, along an axis: the drift region is cut at the boundary
    const XYZVector centre((response.fXRange.first + response.fXRange.second) / 2,
                           (response.fYRange.first + response.fYRange.second) / 2,
                           (response.fZRange.first + response.fZRange.second) / 2);
    const XYZVector normal = readout.GetReadoutNormal(centre);
    if (normal.R() == 0) throw std::runtime_error(error + "the readout has no distance to it");
    int axis = -1;
    for (int a = 0; a < 3; a++) {
        if (std::abs(Component(normal, a)) > 1 - 1e-6) axis = a;
    }
    if (axis < 0) throw std::runtime_error(error + "the normal of the readout must be along x, y or z");
    const double boundary = response.fBoundary;
    const XYZVector onBoundary = centre - (readout.GetDistanceToReadout(centre) - boundary) * normal;
    const double plane = Component(onBoundary, axis) * kCmPerMm;
    if (Component(normal, axis) > 0) {
        low[axis] = std::max(low[axis], plane);
    } else {
        high[axis] = std::min(high[axis], plane);
    }
    if (!(low[axis] < high[axis])) {
        throw std::runtime_error(error + "the field map does not reach the drift region above the boundary");
    }

    Garfield::Sensor sensor;
    sensor.AddComponent(map.get());
    sensor.SetArea(low[0], low[1], low[2], high[0], high[1], high[2]);
    Garfield::AvalancheMC drift;
    drift.SetSensor(&sensor);
    drift.SetDistanceSteps(response.fStepSize * kCmPerMm);
    // Without attachment in the table (e.g. Ar + isobutane), AvalancheMC does not integrate it along each
    // drift line: faster, and without its message (ComputeAlphaEta, status -5) for the lines whose steps
    // cross a cell out of the gas near a wall
    const auto electron = table.fCarriers.find("electron");
    if (electron == table.fCarriers.end() ||
        std::all_of(electron->second.fAttachment.begin(), electron->second.fAttachment.end(),
                    [](double attachment) { return attachment <= 0; })) {
        drift.DisableAttachment();
    }

    response.fMedium = medium.GetName();
    response.fPressure = medium.GetPressure();
    response.fTemperature = medium.GetTemperature();
    response.fTransportTable = table.GetKey();
    // Not a structured binding: captured by the lambdas below (C++17)
    const auto axes = TRestDetectorDriftResponse::PlaneAxes(normal);
    const XYZVector u = axes.first, v = axes.second;
    const size_t nodes = response.GetNumberOfNodes();
    response.fNodes.assign(nodes, TRestDetectorDriftResponse::Node());

    // The response of a node: its electrons drifted
    auto computeNode = [&](size_t index) {
        const XYZVector start = response.GetNodePosition(index);
        const double distance = readout.GetDistanceToReadout(start);
        const XYZVector startOnBoundary = start - (distance - boundary) * normal;
        std::vector<std::array<double, 3>> arrivals;
        const double x[3] = {start.X() * kCmPerMm, start.Y() * kCmPerMm, start.Z() * kCmPerMm};
        // Start points out of the gas of the map (e.g. the corners of a grid around a cylinder, nodes with
        // the field 0) are not drifted: survival 0, without the message of AvalancheMC for each electron
        const bool inArea = x[0] >= low[0] && x[0] <= high[0] && x[1] >= low[1] && x[1] <= high[1] &&
                            x[2] >= low[2] && x[2] <= high[2] && map->GetMedium(x[0], x[1], x[2]) != nullptr;
        if (distance <= boundary + kArrivalTolerance) {
            // On or below the boundary plane: the drift belongs to the amplification
            arrivals.assign(response.fElectrons, std::array<double, 3>{0, 0, 0});
        } else if (inArea) {
            for (int e = 0; e < response.fElectrons; e++) {
                if (!drift.DriftElectron(x[0], x[1], x[2], 0)) continue;
                if (drift.GetNumberOfElectronEndpoints() == 0) continue;
                double x0, y0, z0, t0, x1, y1, z1, t1;
                int status = 0;
                drift.GetElectronEndpoint(0, x0, y0, z0, t0, x1, y1, z1, t1, status);
                const XYZVector end(x1 / kCmPerMm, y1 / kCmPerMm, z1 / kCmPerMm);
                if (std::abs(readout.GetDistanceToReadout(end) - boundary) > kArrivalTolerance) continue;
                const XYZVector displacement = end - startOnBoundary;
                arrivals.push_back({displacement.Dot(u), displacement.Dot(v), t1 * kUsPerNs});
            }
        }
        response.SetNode(index, response.fElectrons, arrivals);
    };
    // The nodes of a process (part of parts: every parts-th node, the drift lengths mixed), with its progress
    auto computePart = [&](size_t part, size_t parts) {
        const size_t partNodes = (nodes - part + parts - 1) / parts;
        size_t done = 0, reported = 0;
        for (size_t index = part; index < nodes; index += parts) {
            computeNode(index);
            done++;
            if (verbose && done * 10 / partNodes > reported) {
                reported = done * 10 / partNodes;
                RESTInfo << "restGarfield: drift response, " << done << " of " << partNodes << " nodes"
                         << (parts > 1 ? " (process " + std::to_string(part + 1) + ")" : "") << RESTendl;
            }
        }
    };

    const size_t parts = std::max<size_t>(1, std::min<size_t>(std::max(jobs, 1), nodes));
    if (parts == 1) {
        computePart(0, 1);
        return;
    }
    // In parallel (Garfield++ is not thread safe): forked processes, each with its own random seed (from the
    // generator of Garfield++, drawn once here), send their nodes through a pipe
    const auto seed = static_cast<unsigned int>(Garfield::RndmUniform() * 1e9) + 1;
    struct Child {
        pid_t pid;
        FILE* pipe;
    };
    std::vector<Child> children;
    for (size_t part = 0; part < parts; part++) {
        int fd[2];
        if (pipe(fd) != 0) throw std::runtime_error(error + "cannot create a pipe");
        std::cout.flush();
        std::fflush(nullptr);
        const pid_t pid = fork();
        if (pid < 0) throw std::runtime_error(error + "cannot create a process");
        if (pid == 0) {
            close(fd[0]);
            int status = 1;
            try {
                Garfield::randomEngine.Seed(seed + part);
                computePart(part, parts);
                FILE* output = fdopen(fd[1], "wb");
                bool ok = output != nullptr;
                for (size_t index = part; ok && index < nodes; index += parts) {
                    ok = WriteNode(output, response.fNodes[index]);
                }
                if (output && fclose(output) == 0 && ok) status = 0;
            } catch (const std::exception& exception) {
                std::cerr << "restGarfield: " << exception.what() << std::endl;
            }
            std::cout.flush();
            // Without the exit handlers of the parent (ROOT)
            _exit(status);
        }
        close(fd[1]);
        children.push_back({pid, fdopen(fd[0], "rb")});
    }
    RESTInfo << "restGarfield: " << nodes << " nodes of the drift response in " << parts << " processes"
             << RESTendl;
    bool ok = true;
    const size_t positionSize = static_cast<size_t>(response.fPositionBins) * response.fPositionBins;
    for (size_t part = 0; part < parts; part++) {
        auto& child = children[part];
        for (size_t index = part; ok && index < nodes; index += parts) {
            ok = child.pipe && ReadNode(child.pipe, response.fNodes[index], positionSize, response.fTimeBins);
        }
        if (child.pipe) fclose(child.pipe);
        int status = 0;
        waitpid(child.pid, &status, 0);
        ok &= WIFEXITED(status) && WEXITSTATUS(status) == 0;
    }
    if (!ok) throw std::runtime_error(error + "a process computing the nodes failed");
}

void RunDriftResponse(const TRestDetectorMedium& medium, TRestDetectorDriftResponse& response, bool verbose,
                      int jobs) {
    const std::string error = "restGarfield: drift response " + response.fResponseFile + ": ";
    const std::string fileName = response.fResponseFile;
    if (!response.fOverwrite && std::filesystem::exists(fileName)) {
        RESTInfo << "restGarfield: drift response " << fileName
                 << " already exists (set overwrite: true to compute it again)" << RESTendl;
        return;
    }
    if (!medium.IsGas()) throw std::runtime_error(error + "the medium must be a gas");

    // The readout of the readouts file
    if (response.fReadoutsFile.empty()) throw std::runtime_error(error + "readouts (a YAML file) is needed");
    const YAML::Node config = TRestYAMLUtils::OpenConfigFile(response.fReadoutsFile);
    const auto section = TRestYAMLUtils::GetMetadataClass(config, "TRestDetectorReadoutManager");
    if (!section.second) {
        throw std::runtime_error(error + "no TRestDetectorReadoutManager in " + response.fReadoutsFile);
    }
    TRestDetectorReadoutManager readouts(section.first, section.second);
    readouts.LoadReadout();
    std::string readoutName = response.fReadoutName;
    const TRestDetectorReadout* readout = readouts.SelectReadout(readoutName, "restGarfield drift response");

    // The transport table of the medium file
    const auto mediumFile = medium.GetMediumFile();
    std::unique_ptr<TRestDetectorTransportTable> table;
    {
        const MediumFileLock lock(mediumFile);
        std::unique_ptr<TFile> file(TFile::Open(mediumFile.string().c_str(), "READ"));
        if (!file || file->IsZombie()) throw std::runtime_error(error + "cannot open " + mediumFile.string());
        std::string key = response.fTransportTable;
        if (key.empty()) {
            const YAML::Node tables = TRestDetectorTransportTable::GetTables(file.get(), medium.GetName());
            if (tables.size() != 1) {
                throw std::runtime_error(error + "the medium file has " + std::to_string(tables.size()) +
                                         " transport tables of " + medium.GetName() +
                                         ": give transportTable (1 needed)");
            }
            key = tables[0]["key"].as<std::string>();
        }
        table = std::make_unique<TRestDetectorTransportTable>(
            TRestDetectorTransportTable::Read(file.get(), medium.GetName(), key));
    }

    SetGridOnReadout(response, *readout);
    RESTInfo << "restGarfield: computing the drift response " << fileName << " ("
             << response.GetNumberOfNodes() << " nodes, " << response.fElectrons << " electrons each)"
             << RESTendl;
    const std::pair<TRestWithUnits, TRestWithUnits>* ranges[3] = {&response.fXRange, &response.fYRange,
                                                                  &response.fZRange};
    for (int axis = 0; axis < 3; axis++) {
        RESTInfo << "    " << "xyz"[axis] << ": [" << ranges[axis]->first << ", " << ranges[axis]->second
                 << "] mm, " << response.fGridNodes[axis] << " nodes" << RESTendl;
    }
    const auto start = std::chrono::steady_clock::now();
    ComputeDriftResponse(medium, *table, *readout, response, verbose, jobs);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    // The time of an electron in one process (the start points out of the gas included), to plan a grid
    const double electrons = static_cast<double>(response.GetNumberOfNodes()) * response.fElectrons;
    RESTInfo << "restGarfield: drift response computed in " << seconds << " s with " << std::max(jobs, 1)
             << " processes, " << seconds * std::max(jobs, 1) / electrons * 1e3 << " ms per electron"
             << RESTendl;

    YAML::Node generator;
    generator["program"] = "Garfield++ AvalancheMC (field map)";
    const std::time_t now = std::time(nullptr);
    char date[32];
    std::strftime(date, sizeof(date), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    generator["date"] = date;
    generator["code"] = REST_Version::GetSummary();
    const std::filesystem::path path(fileName);
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
    if (!file || file->IsZombie()) throw std::runtime_error(error + "cannot create the file");
    response.Write(file.get(), generator);
    RESTInfo << "restGarfield: drift response written to " << fileName << RESTendl;
}
