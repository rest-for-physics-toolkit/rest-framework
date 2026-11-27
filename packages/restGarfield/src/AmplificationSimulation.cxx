#include "AmplificationSimulation.h"

#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "FieldMap.h"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "TFile.h"
#include "TRestDetectorAmplificationResponse.h"
#include "TRestDetectorMedium.h"
#include "TRestVersion.h"
#include "TransportSimulation.h"

namespace {
// Garfield units: cm, ns, eV. REST units: mm, us
constexpr double kCmPerMm = 0.1;
constexpr double kUsPerNs = 1e-3;
/// Initial energy of the primary electron (eV)
constexpr double kInitialEnergy = 0.1;
/// Cells on each side of the cell along a periodic axis inside the area of the sensor (the avalanche spreads
/// over a few of them at most)
constexpr double kPeriodicCells = 10;

/// Upper limit (eV) of the energy table of Magboltz for the avalanches: above the energies reached in the
/// gap (about 40-50 eV at ~60 kV/cm), so that Magboltz does not extend it during the avalanches
constexpr double kMaxElectronEnergy = 100;

using Losses = TRestDetectorAmplificationResponse::Losses;
using Avalanche = TRestDetectorAmplificationResponse::Avalanche;

/// The losses as an array, for the pipes
std::array<int, 5> LossCounts(const Losses& losses) {
    return {losses.fOut, losses.fConductor, losses.fDielectric, losses.fAttached, losses.fOther};
}

/// Writes the result of a process to a pipe: the avalanches of its primaries and the losses
bool WriteResult(FILE* pipe, const Losses& losses, const std::vector<Avalanche>& avalanches) {
    const int count = avalanches.size();
    if (fwrite(&count, sizeof(int), 1, pipe) != 1) return false;
    const auto lost = LossCounts(losses);
    if (fwrite(lost.data(), sizeof(int), lost.size(), pipe) != lost.size()) return false;
    for (const auto& avalanche : avalanches) {
        const double values[5] = {avalanche.fX, avalanche.fY, avalanche.fSize, avalanche.fTime,
                                  avalanche.fShared};
        if (fwrite(values, sizeof(double), 5, pipe) != 5) return false;
    }
    return true;
}

/// Reads the result written by WriteResult, adding it
bool ReadResult(FILE* pipe, Losses& losses, std::vector<Avalanche>& avalanches) {
    int count = 0;
    if (fread(&count, sizeof(int), 1, pipe) != 1) return false;
    std::array<int, 5> lost;
    if (fread(lost.data(), sizeof(int), lost.size(), pipe) != lost.size()) return false;
    losses.fOut += lost[0];
    losses.fConductor += lost[1];
    losses.fDielectric += lost[2];
    losses.fAttached += lost[3];
    losses.fOther += lost[4];
    for (int i = 0; i < count; i++) {
        double values[5];
        if (fread(values, sizeof(double), 5, pipe) != 5) return false;
        avalanches.push_back({values[0], values[1], values[2], values[3], values[4]});
    }
    return true;
}
}  // namespace

void ComputeAmplificationResponse(const TRestDetectorMedium& medium,
                                  TRestDetectorAmplificationResponse& response, bool verbose, int jobs) {
    const std::string error = "restGarfield: amplification response " + response.fResponseFile + ": ";
    if (!medium.IsGas()) throw std::runtime_error(error + "the medium must be a gas");
    auto gas = CreateGas(medium, medium.GetPressure(), medium.GetTemperature());
    gas->SetMaxElectronEnergy(kMaxElectronEnergy);
    const double penning = response.fUsePenning ? response.fPenningRate : 0;
    if (penning > 0 && !gas->EnablePenningTransfer(penning, 0.)) {
        throw std::runtime_error(error + "cannot enable the Penning transfer for " + medium.GetName());
    }

    // The field maps of the cell, the first one taking precedence, repeated along the periodic axes
    std::vector<CSTFieldMapFiles> files;
    for (size_t i = 0; i < response.fFieldMaps.size(); i++) {
        files.push_back(
            {response.fFieldMaps[i], i < response.fDielectricMaps.size() ? response.fDielectricMaps[i] : ""});
    }
    const std::array<double, 3> offset = {response.fFieldMapOffset[0], response.fFieldMapOffset[1],
                                          response.fFieldMapOffset[2]};
    auto map = LoadCSTFieldMaps(files, offset);
    map->SetMedium(gas.get());
    std::array<bool, 3> periodic = {false, false, false};
    for (const char axis : response.fPeriodicity) periodic[axis - 'x'] = true;
    if (periodic[0]) response.fMirror ? map->EnableMirrorPeriodicityX() : map->EnablePeriodicityX();
    if (periodic[1]) response.fMirror ? map->EnableMirrorPeriodicityY() : map->EnablePeriodicityY();
    if (periodic[2]) response.fMirror ? map->EnableMirrorPeriodicityZ() : map->EnablePeriodicityZ();
    double low[3], high[3];  // cm
    map->GetElementaryCell(low[0], low[1], low[2], high[0], high[1], high[2]);
    const double start = response.fStart * kCmPerMm, anode = response.fAnode * kCmPerMm;
    if (start <= low[2] || start > high[2] || anode < low[2] || anode >= start) {
        throw std::runtime_error(error + "start and anode must be inside the field maps along z (" +
                                 std::to_string(low[2] / kCmPerMm) + " to " +
                                 std::to_string(high[2] / kCmPerMm) + " mm), start above anode");
    }

    Garfield::Sensor sensor;
    sensor.AddComponent(map.get());
    double areaLow[3], areaHigh[3];
    for (int a = 0; a < 3; a++) {
        const double margin = periodic[a] ? kPeriodicCells * (high[a] - low[a]) : 0;
        areaLow[a] = low[a] - margin;
        areaHigh[a] = high[a] + margin;
    }
    sensor.SetArea(areaLow[0], areaLow[1], areaLow[2], areaHigh[0], areaHigh[1], areaHigh[2]);
    Garfield::AvalancheMicroscopic avalanche;
    avalanche.SetSensor(&sensor);

    response.fMedium = medium.GetName();
    response.fPressure = medium.GetPressure();
    response.fTemperature = medium.GetTemperature();
    response.fPenning = penning;
    // The cell of the map of the response, x and y in mm
    response.fCellLow = {low[0] / kCmPerMm, low[1] / kCmPerMm};
    response.fCellHigh = {high[0] / kCmPerMm, high[1] / kCmPerMm};

    // The avalanches of a process (part of parts), with their progress: start point, size, time and
    // electrons in another cell of each primary, and where the lost ones ended
    auto computePart = [&](size_t part, size_t parts, Losses& losses, std::vector<Avalanche>& avalanches) {
        const int electrons = response.fElectrons;
        const int partElectrons = (electrons - static_cast<int>(part) + static_cast<int>(parts) - 1) / parts;
        int reported = 0, primaries = 0;
        for (int e = static_cast<int>(part); e < electrons; e += parts) {
            const double x = low[0] + Garfield::RndmUniform() * (high[0] - low[0]);
            const double y = low[1] + Garfield::RndmUniform() * (high[1] - low[1]);
            avalanche.AvalancheElectron(x, y, start, 0, kInitialEnergy);
            primaries++;
            int size = 0, shared = 0;
            double time = 0;
            for (unsigned int i = 0; i < avalanche.GetNumberOfElectronEndpoints(); i++) {
                double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1;
                int status = 0;
                avalanche.GetElectronEndpoint(i, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
                if (z1 > anode) continue;
                size++;
                time += t1;
                // On the pad of another cell (periodic copies)
                if (x1 < low[0] || x1 > high[0] || y1 < low[1] || y1 > high[1]) shared++;
            }
            avalanches.push_back({x / kCmPerMm, y / kCmPerMm, static_cast<double>(size),
                                  size > 0 ? time / size * kUsPerNs : 0, static_cast<double>(shared)});
            if (size > 0) {
                // Collected: nothing else
            } else if (avalanche.GetNumberOfElectronEndpoints() > 0) {
                // Lost: where its first electron that stopped ended
                double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1;
                int status = 0;
                avalanche.GetElectronEndpoint(0, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
                const CSTRegion region = map->GetRegionNear(x1, y1, z1);
                if (status == Garfield::StatusAttached) {
                    losses.fAttached++;
                } else if (status == Garfield::StatusLeftDriftArea || region == CSTRegion::kOutside) {
                    losses.fOut++;
                } else if (region == CSTRegion::kDielectric) {
                    losses.fDielectric++;
                } else if (region == CSTRegion::kConductor) {
                    losses.fConductor++;
                } else {
                    losses.fOther++;
                }
            } else {
                losses.fOther++;
            }
            if (verbose && primaries * 10 / std::max(1, partElectrons) > reported) {
                reported = primaries * 10 / std::max(1, partElectrons);
                RESTInfo << "restGarfield: amplification response, " << primaries << " of " << partElectrons
                         << " avalanches" << (parts > 1 ? " (process " + std::to_string(part + 1) + ")" : "")
                         << RESTendl;
            }
        }
    };

    Losses losses;
    std::vector<Avalanche> avalanches;
    const size_t parts = std::max<size_t>(1, std::min<size_t>(std::max(jobs, 1), response.fElectrons));
    if (parts == 1) {
        computePart(0, 1, losses, avalanches);
        response.SetResult(avalanches, losses);
        return;
    }
    // In parallel (Garfield++ is not thread safe): forked processes, each with its own random seed (from the
    // generator of Garfield++, drawn once here), send their avalanches through a pipe
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
                Losses partLosses;
                std::vector<Avalanche> partAvalanches;
                computePart(part, parts, partLosses, partAvalanches);
                FILE* output = fdopen(fd[1], "wb");
                const bool ok = output != nullptr && WriteResult(output, partLosses, partAvalanches);
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
    RESTInfo << "restGarfield: " << response.fElectrons << " avalanches of the amplification response in "
             << parts << " processes" << RESTendl;
    bool ok = true;
    for (auto& child : children) {
        ok = ok && child.pipe && ReadResult(child.pipe, losses, avalanches);
        if (child.pipe) fclose(child.pipe);
        int status = 0;
        waitpid(child.pid, &status, 0);
        ok &= WIFEXITED(status) && WEXITSTATUS(status) == 0;
    }
    if (!ok) throw std::runtime_error(error + "a process computing the avalanches failed");
    response.SetResult(avalanches, losses);
}

void RunAmplificationResponse(const TRestDetectorMedium& medium, TRestDetectorAmplificationResponse& response,
                              bool verbose, int jobs) {
    const std::string error = "restGarfield: amplification response " + response.fResponseFile + ": ";
    const std::string fileName = response.fResponseFile;
    if (!response.fOverwrite && std::filesystem::exists(fileName)) {
        RESTInfo << "restGarfield: amplification response " << fileName
                 << " already exists (set overwrite: true to compute it again)" << RESTendl;
        return;
    }
    RESTInfo << "restGarfield: computing the amplification response " << fileName << " ("
             << response.fElectrons << " avalanches)" << RESTendl;
    const auto start = std::chrono::steady_clock::now();
    ComputeAmplificationResponse(medium, response, verbose, jobs);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    RESTInfo << "restGarfield: amplification response computed in " << seconds << " s with "
             << std::max(jobs, 1) << " processes, " << seconds * std::max(jobs, 1) / response.fElectrons
             << " s per avalanche; transparency " << response.GetTransparency() << ", mean gain "
             << response.fMeanSize << ", theta " << response.GetTheta() << RESTendl;
    const auto& lost = response.fLosses;
    RESTInfo << "restGarfield: primaries lost " << response.fPrimaries - response.fCollected
             << ": out of the cell " << lost.fOut << ", on a conductor " << lost.fConductor
             << ", on a dielectric " << lost.fDielectric << ", attached " << lost.fAttached << ", other "
             << lost.fOther << RESTendl;
    if (response.fAnodeElectrons > 0) {
        RESTInfo
            << "restGarfield: electrons at the anode out of the cell of their primary (charge sharing, not "
               "applied): "
            << response.fSharedElectrons / response.fAnodeElectrons << RESTendl;
    }
    if (response.HasCellMap()) {
        double minimum = 1, maximum = 0;
        for (size_t bin = 0; bin < response.fCellPrimaries.size(); bin++) {
            if (response.fCellPrimaries[bin] <= 0) continue;
            minimum = std::min(minimum, response.GetCellTransparency(bin));
            maximum = std::max(maximum, response.GetCellTransparency(bin));
        }
        RESTInfo << "restGarfield: map of the cell, " << response.fCellBins << " x " << response.fCellBins
                 << " bins, " << static_cast<double>(response.fPrimaries) / response.fCellPrimaries.size()
                 << " primaries per bin, transparency from " << minimum << " to " << maximum << RESTendl;
    }

    YAML::Node generator;
    generator["program"] = "Garfield++ AvalancheMicroscopic (field maps of a unit cell)";
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
    RESTInfo << "restGarfield: amplification response written to " << fileName << RESTendl;
}
