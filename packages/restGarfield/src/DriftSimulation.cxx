#include "DriftSimulation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "FieldMap.h"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/MediumGas.hh"
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
}  // namespace

void ComputeDriftResponse(const TRestDetectorMedium& medium, const TRestDetectorTransportTable& table,
                          const TRestDetectorReadout& readout, TRestDetectorDriftResponse& response,
                          bool verbose) {
    const std::string error = "restGarfield: drift response " + response.fResponseFile + ": ";
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

    response.fMedium = medium.GetName();
    response.fPressure = medium.GetPressure();
    response.fTemperature = medium.GetTemperature();
    response.fTransportTable = table.GetKey();
    const auto [u, v] = TRestDetectorDriftResponse::PlaneAxes(normal);
    const size_t nodes = response.GetNumberOfNodes();
    response.fNodes.assign(nodes, TRestDetectorDriftResponse::Node());
    size_t reported = 0;
    for (size_t index = 0; index < nodes; index++) {
        const XYZVector start = response.GetNodePosition(index);
        const double distance = readout.GetDistanceToReadout(start);
        const XYZVector startOnBoundary = start - (distance - boundary) * normal;
        std::vector<std::array<double, 3>> arrivals;
        const double x[3] = {start.X() * kCmPerMm, start.Y() * kCmPerMm, start.Z() * kCmPerMm};
        const bool inArea = x[0] >= low[0] && x[0] <= high[0] && x[1] >= low[1] && x[1] <= high[1] &&
                            x[2] >= low[2] && x[2] <= high[2];
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
        if (verbose && (index + 1) * 10 / nodes > reported) {
            reported = (index + 1) * 10 / nodes;
            RESTInfo << "restGarfield: drift response, " << index + 1 << " of " << nodes << " nodes"
                     << RESTendl;
        }
    }
}

void RunDriftResponse(const TRestDetectorMedium& medium, TRestDetectorDriftResponse& response, bool verbose) {
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

    RESTInfo << "restGarfield: computing the drift response " << fileName << " ("
             << response.GetNumberOfNodes() << " nodes, " << response.fElectrons << " electrons each)"
             << RESTendl;
    ComputeDriftResponse(medium, *table, *readout, response, verbose);

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
