/// \class TRestDetectorAmplificationResponse
/// Amplification of the drifted electrons in a non-uniform field (the field map of a unit cell of the
/// amplification structure, e.g. a microbulk Micromegas exported from CST Studio Suite) as a response, for
/// special simulations; the amplification in a uniform field across the gap (parallel plates) is the gain
/// table of the medium (TRestDetectorGainTable). Both are used by TRestDetectorAvalancheProcess.
///
/// restGarfield computes it (the `amplification` key of a TRestDetectorTransportManager, after the drift
/// responses): it loads the field maps in Garfield++ (restGarfield LoadCSTFieldMaps: the first map takes
/// precedence where they overlap; the nodes with the field 0, the conductors, and the dielectric nodes of
/// the D maps stop the electrons), repeats the cell along the `periodicity` axes (simple or `mirror`) and
/// simulates the microscopic avalanches (AvalancheMicroscopic, the Garfield gas of the composition of the
/// medium at its pressure and temperature, Magboltz cross sections, with the Penning transfer if
/// usePenning) of `electrons` primary electrons started at random over the cell (x and y) on the plane z =
/// `start`, in the frame of the maps moved by fieldMapOffset. The electrons drift towards -z; an electron of
/// the avalanche ending below z = `anode` is collected by the anode. A primary whose avalanche brings no
/// electron to the anode is lost (on the mesh, on a dielectric, attached); the transparency is the fraction
/// of the primaries with an avalanche at the anode. Where the lost primaries ended is counted (Losses, from
/// the first electron of the avalanche that stopped): out of the maps or the cell (e.g. back above the start
/// plane), on a conductor (the nodes with the field 0 next to its end point, e.g. the mesh), on a dielectric
/// (the dielectric nodes next to it, e.g. the kapton), attached, or other. `start` is the boundary plane of
/// the drift response (TRestDetectorDriftResponse `boundary`, measured from the mesh): the drift response
/// ends there and the amplification starts there, so the electrons lost on the mesh are counted once.
///
/// Response: the number of primaries and of collected avalanches, the histogram of the size of the
/// collected avalanches (electrons at the anode, sizeBins over [0, the largest size]) with its mean and
/// variance, and the histogram of the time from the start plane to the anode (the mean time of the electrons
/// of the avalanche at the anode, timeBins over the range of the times). The map of the cell (cellBins x
/// cellBins bins over x and y of the cell, by the start point of the primaries): the transparency and the
/// mean size of each bin, e.g. the losses over the solid mesh around the holes. The electrons at the anode
/// out of the cell where their primary started (on the pad of another cell) are counted, to check that the
/// charge stays on the pad of its cell.
/// TRestDetectorAvalancheProcess with `responseFile` uses it for each drifted electron: it passes with the
/// transparency, its gain is sampled from the size histogram and its time is delayed by a time sampled from
/// the time histogram. With the map, the transparency is that of the bin of the electron in the cell and the
/// gain is scaled by the mean size of the bin over the mean of all; the cell is the volume of the readout
/// containing the electron (a pixel, TRestDetectorReadout::GetPixelPosition): its x and y from the lower
/// corner, as fractions of its size, are those of the cell (so the cell must be centred on the pixel, with
/// its size).
///
/// Parameters (`amplification` key): fieldMaps (required), dielectricMaps ([]: none; one per field map in
/// the same order, "" for a map without), fieldMapOffset ([0, 0, 0] mm), periodicity ("xy"), mirror (false),
/// start and anode (z in mm, start > anode), electrons (1000), usePenning (false), penningRate (0), sizeBins
/// (100), timeBins (50), cellBins (10; 0: no map), responseFile and overwrite (false: an existing
/// responseFile is not computed again). About 100 primaries per bin of the map give its transparency to
/// +-0.04: electrons ~ 100 cellBins^2.
/// One response per operating point (the voltages of the CST model: mesh voltage and drift field).
///
///    amplification:
///      fieldMaps: ["pixel_mesh320_E_fine.txt", "pixel_mesh320_E_coarse.txt"]
///      dielectricMaps: ["pixel_mesh320_D_fine.txt"]
///      periodicity: "xy"
///      start: 1 mm
///      anode: -54 um
///      electrons: 2000
///      responseFile: "amplification_320V.root"
///
/// Response file: the metadata `AmplificationResponse` in RESTMetadataStore (YAML: `configuration`, the
/// `amplification` keys; medium, pressure, temperature, penning, primaries, collected, losses (out,
/// conductor, dielectric, attached, other), meanSize, sizeVariance, sizeMax, sizeCounts, tMin, tMax,
/// timeCounts, cellLow, cellHigh, cellPrimaries, cellCollected, cellSizeSum, anodeElectrons,
/// sharedElectrons and `generator`, the program, date and code that computed it), the histograms TH1D
/// `AmplificationSize` and `AmplificationTime` (us) and the map of the cell, TH2D
/// `AmplificationTransparencyMap` and `AmplificationGainMap` (mean size), to draw them.
///
/// \author JuanAn García <juanangp@unizar.es>

#include "TRestDetectorAmplificationResponse.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"

static const bool TRestDetectorAmplificationResponse_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using A = TRestDetectorAmplificationResponse;
    reg.RegisterField<A>("fieldMaps", &A::fFieldMaps);
    reg.RegisterField<A>("dielectricMaps", &A::fDielectricMaps);
    reg.RegisterField<A>("fieldMapOffset", &A::fFieldMapOffset);
    reg.RegisterField<A>("periodicity", &A::fPeriodicity);
    reg.RegisterField<A>("mirror", &A::fMirror);
    reg.RegisterField<A>("start", &A::fStart);
    reg.RegisterField<A>("anode", &A::fAnode);
    reg.RegisterField<A>("electrons", &A::fElectrons);
    reg.RegisterField<A>("usePenning", &A::fUsePenning);
    reg.RegisterField<A>("penningRate", &A::fPenningRate);
    reg.RegisterField<A>("sizeBins", &A::fSizeBins);
    reg.RegisterField<A>("timeBins", &A::fTimeBins);
    reg.RegisterField<A>("cellBins", &A::fCellBins);
    reg.RegisterField<A>("responseFile", &A::fResponseFile);
    reg.RegisterField<A>("overwrite", &A::fOverwrite);
    return true;
}();

namespace {
/// Bin of a value in [min, max] with n bins (0 if the range is empty)
int Bin(double value, double min, double max, int n) {
    if (!(max > min)) return 0;
    const int bin = static_cast<int>(std::floor((value - min) / (max - min) * n));
    return std::clamp(bin, 0, n - 1);
}

/// A value sampled from counts over [min, max]: a bin with its counts, uniform within it
double Sample(const std::vector<int>& counts, double min, double max, TRandom& random) {
    long total = 0;
    for (const int count : counts) total += count;
    if (total <= 0 || !(max > min)) return min;
    long target = static_cast<long>(random.Uniform() * total);
    size_t bin = 0;
    for (; bin + 1 < counts.size(); bin++) {
        if (target < counts[bin]) break;
        target -= counts[bin];
    }
    const double width = (max - min) / counts.size();
    return min + (bin + random.Uniform()) * width;
}

/// A histogram of counts over [min, max] to draw it
void WriteHistogram(const std::string& name, const std::string& title, const std::vector<int>& counts,
                    double min, double max) {
    TH1D histogram(name.c_str(), title.c_str(), counts.size(), min, max > min ? max : min + 1);
    histogram.SetDirectory(nullptr);
    for (size_t bin = 0; bin < counts.size(); bin++) histogram.SetBinContent(bin + 1, counts[bin]);
    histogram.Write(name.c_str(), TObject::kOverwrite);
}
}  // namespace

TRestDetectorAmplificationResponse::TRestDetectorAmplificationResponse() : TRestMetadata() {
    fName = "TRestDetectorAmplificationResponse";
}

TRestDetectorAmplificationResponse::TRestDetectorAmplificationResponse(const std::string& instanceName,
                                                                       const YAML::Node& node)
    : TRestMetadata(instanceName, node) {
    LoadConfig();
}

void TRestDetectorAmplificationResponse::LoadConfig() {
    UpdateParamsFromYAML<TRestDetectorAmplificationResponse>(fNode);
    const std::string prefix = "TRestDetectorAmplificationResponse '" + fName + "': ";
    if (fFieldMaps.empty()) throw std::runtime_error(prefix + "fieldMaps is needed");
    if (fDielectricMaps.size() > fFieldMaps.size()) {
        throw std::runtime_error(prefix + "dielectricMaps has more files than fieldMaps");
    }
    if (fResponseFile.empty()) throw std::runtime_error(prefix + "responseFile is needed");
    if (fFieldMapOffset.size() != 3) throw std::runtime_error(prefix + "fieldMapOffset must have 3 values");
    if (fPeriodicity.find_first_not_of("xyz") != std::string::npos) {
        throw std::runtime_error(prefix + "periodicity must be the axes of the cell (e.g. \"xy\")");
    }
    if (!(fStart > fAnode)) throw std::runtime_error(prefix + "start must be above anode (start > anode)");
    if (fElectrons <= 0) throw std::runtime_error(prefix + "electrons must be > 0");
    if (fSizeBins <= 0 || fTimeBins <= 0)
        throw std::runtime_error(prefix + "sizeBins and timeBins must be > 0");
    if (fCellBins < 0) throw std::runtime_error(prefix + "cellBins must be >= 0");
    if (fPenningRate < 0 || fPenningRate > 1)
        throw std::runtime_error(prefix + "penningRate must be in [0, 1]");
    UpdateYAMLFromParams<TRestDetectorAmplificationResponse>(fNode);
}

void TRestDetectorAmplificationResponse::SetResult(const std::vector<Avalanche>& avalanches,
                                                   const Losses& losses) {
    fPrimaries = static_cast<int>(avalanches.size());
    fCollected = 0;
    fLosses = losses;
    fSizeCounts.assign(fSizeBins, 0);
    fTimeCounts.assign(fTimeBins, 0);
    fMeanSize = fSizeVariance = fSizeMax = fTMin = fTMax = 0;
    fAnodeElectrons = fSharedElectrons = 0;
    const bool map = fCellBins > 0 && fCellHigh[0] > fCellLow[0] && fCellHigh[1] > fCellLow[1];
    const size_t cells = map ? static_cast<size_t>(fCellBins) * fCellBins : 0;
    fCellPrimaries.assign(cells, 0);
    fCellCollected.assign(cells, 0);
    fCellSizeSum.assign(cells, 0);

    double sum = 0, sum2 = 0;
    for (const auto& avalanche : avalanches) {
        const int cell = map ? GetCellBin((avalanche.fX - fCellLow[0]) / (fCellHigh[0] - fCellLow[0]),
                                          (avalanche.fY - fCellLow[1]) / (fCellHigh[1] - fCellLow[1]))
                             : -1;
        if (map) fCellPrimaries[cell]++;
        if (!(avalanche.fSize > 0)) continue;
        if (map) {
            fCellCollected[cell]++;
            fCellSizeSum[cell] += avalanche.fSize;
        }
        fTMin = fCollected == 0 ? avalanche.fTime : std::min(fTMin, avalanche.fTime);
        fTMax = fCollected == 0 ? avalanche.fTime : std::max(fTMax, avalanche.fTime);
        fCollected++;
        sum += avalanche.fSize;
        sum2 += avalanche.fSize * avalanche.fSize;
        fAnodeElectrons += avalanche.fSize;
        fSharedElectrons += avalanche.fShared;
        fSizeMax = std::max(fSizeMax, avalanche.fSize);
    }
    if (fCollected == 0) return;
    const double n = fCollected;
    fMeanSize = sum / n;
    fSizeVariance = n > 1 ? std::max(0.0, (sum2 - n * fMeanSize * fMeanSize) / (n - 1)) : 0;
    for (const auto& avalanche : avalanches) {
        if (!(avalanche.fSize > 0)) continue;
        fSizeCounts[Bin(avalanche.fSize, 0, fSizeMax, fSizeBins)]++;
        fTimeCounts[Bin(avalanche.fTime, fTMin, fTMax, fTimeBins)]++;
    }
}

double TRestDetectorAmplificationResponse::GetTransparency() const {
    return fPrimaries > 0 ? static_cast<double>(fCollected) / fPrimaries : 0;
}

double TRestDetectorAmplificationResponse::GetTheta() const {
    return fSizeVariance > 0 ? fMeanSize * fMeanSize / fSizeVariance - 1 : 0;
}

double TRestDetectorAmplificationResponse::SampleSize(TRandom& random) const {
    return Sample(fSizeCounts, 0, fSizeMax, random);
}

double TRestDetectorAmplificationResponse::SampleTime(TRandom& random) const {
    return Sample(fTimeCounts, fTMin, fTMax, random);
}

int TRestDetectorAmplificationResponse::GetCellBin(double fractionX, double fractionY) const {
    if (fCellBins <= 0) return 0;
    return Bin(fractionX, 0, 1, fCellBins) + fCellBins * Bin(fractionY, 0, 1, fCellBins);
}

double TRestDetectorAmplificationResponse::GetCellTransparency(int bin) const {
    if (bin < 0 || static_cast<size_t>(bin) >= fCellPrimaries.size() || fCellPrimaries[bin] <= 0) {
        return GetTransparency();
    }
    return static_cast<double>(fCellCollected[bin]) / fCellPrimaries[bin];
}

double TRestDetectorAmplificationResponse::GetCellGainFactor(int bin) const {
    if (bin < 0 || static_cast<size_t>(bin) >= fCellCollected.size() || fCellCollected[bin] <= 0 ||
        !(fMeanSize > 0)) {
        return 1;
    }
    return fCellSizeSum[bin] / fCellCollected[bin] / fMeanSize;
}

void TRestDetectorAmplificationResponse::Write(TFile* file, const YAML::Node& generator) const {
    if (fSizeCounts.size() != static_cast<size_t>(fSizeBins) ||
        fTimeCounts.size() != static_cast<size_t>(fTimeBins)) {
        throw std::runtime_error("TRestDetectorAmplificationResponse: the response has not the bins set");
    }
    YAML::Node metadata;
    // The configuration as read (LoadConfig writes the values in REST units)
    metadata["configuration"] = YAML::Clone(fNode);
    metadata["medium"] = fMedium;
    metadata["pressure"] = fPressure;
    metadata["temperature"] = fTemperature;
    metadata["penning"] = fPenning;
    metadata["primaries"] = fPrimaries;
    metadata["collected"] = fCollected;
    metadata["losses"]["out"] = fLosses.fOut;
    metadata["losses"]["conductor"] = fLosses.fConductor;
    metadata["losses"]["dielectric"] = fLosses.fDielectric;
    metadata["losses"]["attached"] = fLosses.fAttached;
    metadata["losses"]["other"] = fLosses.fOther;
    metadata["meanSize"] = fMeanSize;
    metadata["sizeVariance"] = fSizeVariance;
    metadata["sizeMax"] = fSizeMax;
    metadata["sizeCounts"] = fSizeCounts;
    metadata["tMin"] = fTMin;
    metadata["tMax"] = fTMax;
    metadata["timeCounts"] = fTimeCounts;
    metadata["cellLow"] = std::vector<double>(fCellLow.begin(), fCellLow.end());
    metadata["cellHigh"] = std::vector<double>(fCellHigh.begin(), fCellHigh.end());
    metadata["cellPrimaries"] = fCellPrimaries;
    metadata["cellCollected"] = fCellCollected;
    metadata["cellSizeSum"] = fCellSizeSum;
    metadata["anodeElectrons"] = fAnodeElectrons;
    metadata["sharedElectrons"] = fSharedElectrons;
    metadata["generator"] = YAML::Clone(generator);
    TRestMetadata::WriteMetadata(file, kName, metadata);

    file->cd();
    WriteHistogram("AmplificationSize", "Avalanche size;electrons at the anode;avalanches", fSizeCounts, 0,
                   fSizeMax);
    WriteHistogram("AmplificationTime", "Time to the anode;time (#mus);avalanches", fTimeCounts, fTMin,
                   fTMax);
    if (!HasCellMap()) return;
    TH2D transparency("AmplificationTransparencyMap", "Transparency by start point;x (mm);y (mm)", fCellBins,
                      fCellLow[0], fCellHigh[0], fCellBins, fCellLow[1], fCellHigh[1]);
    TH2D gain("AmplificationGainMap", "Mean avalanche size by start point;x (mm);y (mm)", fCellBins,
              fCellLow[0], fCellHigh[0], fCellBins, fCellLow[1], fCellHigh[1]);
    for (auto* histogram : {&transparency, &gain}) histogram->SetDirectory(nullptr);
    for (int binY = 0; binY < fCellBins; binY++) {
        for (int binX = 0; binX < fCellBins; binX++) {
            const int bin = binX + fCellBins * binY;
            if (fCellPrimaries[bin] > 0) {
                transparency.SetBinContent(binX + 1, binY + 1, GetCellTransparency(bin));
            }
            if (fCellCollected[bin] > 0) {
                gain.SetBinContent(binX + 1, binY + 1, GetCellGainFactor(bin) * fMeanSize);
            }
        }
    }
    transparency.Write(transparency.GetName(), TObject::kOverwrite);
    gain.Write(gain.GetName(), TObject::kOverwrite);
}

TRestDetectorAmplificationResponse TRestDetectorAmplificationResponse::Read(TFile* file) {
    const std::string error =
        "TRestDetectorAmplificationResponse: file " + std::string(file->GetName()) + ": ";
    const YAML::Node metadata = TRestMetadata::ReadMetadata(file, kName);
    if (!metadata || !metadata.IsMap() || !metadata["configuration"]) {
        throw std::runtime_error(error + "no " + kName + " metadata");
    }
    TRestDetectorAmplificationResponse response("amplification", YAML::Clone(metadata["configuration"]));
    response.fMedium = metadata["medium"].as<std::string>("");
    response.fPressure = metadata["pressure"].as<double>(0);
    response.fTemperature = metadata["temperature"].as<double>(0);
    response.fPenning = metadata["penning"].as<double>(0);
    response.fPrimaries = metadata["primaries"].as<int>(0);
    response.fCollected = metadata["collected"].as<int>(0);
    if (const YAML::Node losses = metadata["losses"]) {
        response.fLosses.fOut = losses["out"].as<int>(0);
        response.fLosses.fConductor = losses["conductor"].as<int>(0);
        response.fLosses.fDielectric = losses["dielectric"].as<int>(0);
        response.fLosses.fAttached = losses["attached"].as<int>(0);
        response.fLosses.fOther = losses["other"].as<int>(0);
    }
    response.fMeanSize = metadata["meanSize"].as<double>(0);
    response.fSizeVariance = metadata["sizeVariance"].as<double>(0);
    response.fSizeMax = metadata["sizeMax"].as<double>(0);
    response.fSizeCounts = metadata["sizeCounts"].as<std::vector<int>>(std::vector<int>());
    response.fTMin = metadata["tMin"].as<double>(0);
    response.fTMax = metadata["tMax"].as<double>(0);
    response.fTimeCounts = metadata["timeCounts"].as<std::vector<int>>(std::vector<int>());
    // The map of the cell (none in the responses computed before it)
    const auto low = metadata["cellLow"].as<std::vector<double>>(std::vector<double>());
    const auto high = metadata["cellHigh"].as<std::vector<double>>(std::vector<double>());
    if (low.size() == 2 && high.size() == 2) {
        response.fCellLow = {low[0], low[1]};
        response.fCellHigh = {high[0], high[1]};
    }
    response.fCellPrimaries = metadata["cellPrimaries"].as<std::vector<int>>(std::vector<int>());
    response.fCellCollected = metadata["cellCollected"].as<std::vector<int>>(std::vector<int>());
    response.fCellSizeSum = metadata["cellSizeSum"].as<std::vector<double>>(std::vector<double>());
    response.fAnodeElectrons = metadata["anodeElectrons"].as<double>(0);
    response.fSharedElectrons = metadata["sharedElectrons"].as<double>(0);
    const size_t cells =
        response.fCellPrimaries.empty() ? 0 : static_cast<size_t>(response.fCellBins) * response.fCellBins;
    if (response.fSizeCounts.size() != static_cast<size_t>(response.fSizeBins) ||
        response.fTimeCounts.size() != static_cast<size_t>(response.fTimeBins) ||
        response.fCollected > response.fPrimaries || response.fCellPrimaries.size() != cells ||
        response.fCellCollected.size() != cells || response.fCellSizeSum.size() != cells) {
        throw std::runtime_error(error + "the response is not consistent (bins or counts)");
    }
    return response;
}
