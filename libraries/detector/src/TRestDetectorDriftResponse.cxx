/// \class TRestDetectorDriftResponse
/// Drift of the electrons in a non-uniform electric field (a field map exported from CST Studio Suite) as a
/// response per start point, for special simulations: TRestDetectorHitsDriftProcess with `responseFile`,
/// instead of the uniform drift field of the medium.
///
/// restGarfield computes it (the `drift` key of a TRestDetectorTransportManager, after the transport and gain
/// tables): it loads the field map in Garfield++ (ComponentGrid), builds the Garfield gas from a transport
/// table of the medium file (without Magboltz) and drifts `electrons` primary electrons (AvalancheMC, steps
/// of `stepSize`, with diffusion, and attachment if the table has it) from each node of a regular grid of
/// start points in the detector frame (xRange, yRange, zRange, gridStep). The drift ends at the boundary
/// plane, at the distance
/// `boundary` from the readout (TRestDetectorReadout::GetDistanceToReadout), above the mesh, where the
/// amplification starts; the drift region is cut there (the normal of the readout must be along an axis).
/// An electron reaches the boundary if it ends on that plane; ending elsewhere (cathode, walls, out of the
/// map) or attached, it is lost. Each electron drifts independently, so the arrivals from a point depend only
/// on the field and the gas: the response replaces the drift of every electron of every event, up to the
/// spacing of the grid.
///
/// Response of each node: the number of primary electrons and of those that arrive, and the histograms of the
/// arrival position, the displacement (u, v) from the node projected on the boundary plane along the axes of
/// the readout plane (PlaneAxes; e.g. normal +z: u = x, v = y), positionBins x positionBins over the range of
/// the arrivals of the node, and of the drift time (timeBins). TRestDetectorHitsDriftProcess uses it
/// for each electron of a hit: it survives with the fraction arrived / primaries interpolated between the 8
/// neighbouring nodes (GetSurvival), and its displacement and time are sampled from one of those nodes,
/// chosen with its interpolation weight times its survival (SampleArrival); the displacement is added to the
/// hit projected on the readout, exact for a uniform field.
///
/// Parameters (`drift` key): fieldMap, fieldMapOffset ([0, 0, 0] mm), readouts (YAML file with the
/// TRestDetectorReadoutManager), readout (empty: its data readout), boundary (0 mm), transportTable (key of
/// the table of the medium file; empty: the only one), xRange, yRange, zRange (mm; a range not given along an
/// axis of the readout plane is the extent of the readout channels, TRestDetectorReadout::
/// GetChannelsBoundingBox, plus `margin` on each side; the range along the normal is needed), margin (5 mm),
/// gridStep (required: the largest spacing of the nodes along x, y and z; each range is divided into the
/// fewest equal steps not larger than it, at least 2 nodes per axis), electrons (1000, at most 65535),
/// stepSize (1 mm), positionBins (20), timeBins (20), responseFile and overwrite (false: an existing
/// responseFile is not computed again). The pressure and temperature are those of the medium. The time of the
/// computation grows with the nodes, the electrons and the drift length over stepSize. The response file
/// keeps the ranges used (those of the readout included).
///
///    drift:
///      fieldMap: "drift_field.txt"
///      fieldMapOffset: [0 mm, 0 mm, 0 mm]
///      readouts: "readouts.yaml"
///      boundary: 0.5 mm
///      zRange: [0 mm, 190 mm]          # x and y: the readout channels plus margin
///      margin: 5 mm
///      gridStep: [4 mm, 4 mm, 10 mm]
///      electrons: 1000
///      responseFile: "drift_response.root"
///
/// Response file: the metadata `DriftResponse` in RESTMetadataStore (YAML: `configuration`, the `drift` keys;
/// medium, pressure, temperature, transportTable, positionBins, timeBins and `generator`, the program, date
/// and code that computed it) and the tree `DriftResponse`, an entry per node (x fastest, then y, then z)
/// with primaries, arrived, range (float[6]: u min, u max, v min, v max, t min, t max; mm and us),
/// positionCounts (UShort_t[positionBins * positionBins], u bin + positionBins * v bin) and timeCounts
/// (UShort_t[timeBins]); hence electrons <= 65535. The file is compressed with LZMA. To draw the response,
/// each entry also has the position of the node (x, y, z, mm) and its summary (GetNodeSummary, floats):
/// survival, meanU, meanV (mean displacement on the readout plane, mm: the distortion of the field), meanT
/// (mean drift time, us), sigmaU, sigmaV (transverse spread, mm) and sigmaT (us), e.g.
/// `DriftResponse->Draw("meanU:x", "z == 90")`; restGarfield macro REST_Garfield_DrawDriftResponse draws
/// them as maps on the readout plane.
///
/// \author JuanAn García <juanangp@unizar.es>

#include "TRestDetectorDriftResponse.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include "Compression.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TRandom.h"
#include "TTree.h"

using ROOT::Math::XYZVector;

static const bool TRestDetectorDriftResponse_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using D = TRestDetectorDriftResponse;
    reg.RegisterField<D>("fieldMap", &D::fFieldMap);
    reg.RegisterField<D>("fieldMapOffset", &D::fFieldMapOffset);
    reg.RegisterField<D>("readouts", &D::fReadoutsFile);
    reg.RegisterField<D>("readout", &D::fReadoutName);
    reg.RegisterField<D>("boundary", &D::fBoundary);
    reg.RegisterField<D>("transportTable", &D::fTransportTable);
    reg.RegisterField<D>("xRange", &D::fXRange);
    reg.RegisterField<D>("yRange", &D::fYRange);
    reg.RegisterField<D>("zRange", &D::fZRange);
    reg.RegisterField<D>("gridStep", &D::fGridStep);
    reg.RegisterField<D>("margin", &D::fMargin);
    reg.RegisterField<D>("electrons", &D::fElectrons);
    reg.RegisterField<D>("stepSize", &D::fStepSize);
    reg.RegisterField<D>("positionBins", &D::fPositionBins);
    reg.RegisterField<D>("timeBins", &D::fTimeBins);
    reg.RegisterField<D>("responseFile", &D::fResponseFile);
    reg.RegisterField<D>("overwrite", &D::fOverwrite);
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
}  // namespace

TRestDetectorDriftResponse::TRestDetectorDriftResponse() : TRestMetadata() {
    fName = "TRestDetectorDriftResponse";
}

TRestDetectorDriftResponse::TRestDetectorDriftResponse(const std::string& instanceName,
                                                       const YAML::Node& node)
    : TRestMetadata(instanceName, node) {
    LoadConfig();
}

void TRestDetectorDriftResponse::LoadConfig() {
    UpdateParamsFromYAML<TRestDetectorDriftResponse>(fNode);
    const std::string prefix = "TRestDetectorDriftResponse '" + fName + "': ";
    if (fFieldMap.empty()) throw std::runtime_error(prefix + "fieldMap is needed");
    if (fResponseFile.empty()) throw std::runtime_error(prefix + "responseFile is needed");
    if (fFieldMapOffset.size() != 3) throw std::runtime_error(prefix + "fieldMapOffset must have 3 values");
    for (int axis = 0; axis < 3; axis++) {
        const auto& range = Range(axis);
        if (!HasRange(axis) && !(range.first == 0 && range.second == 0)) {
            throw std::runtime_error(prefix + "xRange, yRange and zRange must be [min, max] with min < max");
        }
    }
    if (fGridStep.size() != 3 || !(*std::min_element(fGridStep.begin(), fGridStep.end()) > 0)) {
        throw std::runtime_error(prefix + "gridStep must be 3 lengths > 0");
    }
    if (fMargin < 0) throw std::runtime_error(prefix + "margin must be >= 0");
    // The counts of a node are stored as 16-bit integers
    if (fElectrons <= 0 || fElectrons > 65535) {
        throw std::runtime_error(prefix + "electrons must be > 0 and <= 65535");
    }
    if (fStepSize <= 0) throw std::runtime_error(prefix + "stepSize must be > 0");
    if (fPositionBins <= 0 || fTimeBins <= 0) {
        throw std::runtime_error(prefix + "positionBins and timeBins must be > 0");
    }
    if (fBoundary < 0) throw std::runtime_error(prefix + "boundary must be >= 0");
    UpdateGrid();
    UpdateYAMLFromParams<TRestDetectorDriftResponse>(fNode);
}

std::pair<TRestWithUnits, TRestWithUnits>& TRestDetectorDriftResponse::Range(int axis) {
    return axis == 0 ? fXRange : (axis == 1 ? fYRange : fZRange);
}

const std::pair<TRestWithUnits, TRestWithUnits>& TRestDetectorDriftResponse::Range(int axis) const {
    return axis == 0 ? fXRange : (axis == 1 ? fYRange : fZRange);
}

bool TRestDetectorDriftResponse::HasRange(int axis) const { return Range(axis).first < Range(axis).second; }

void TRestDetectorDriftResponse::UpdateGrid() {
    for (int axis = 0; axis < 3; axis++) {
        fGridNodes[axis] = 0;
        if (!HasRange(axis) || fGridStep.size() != 3) continue;
        // The fewest equal steps not larger than gridStep (a relative tolerance for an exact division)
        const double steps = (Range(axis).second - Range(axis).first) / fGridStep[axis];
        fGridNodes[axis] = std::max(2, static_cast<int>(std::ceil(steps * (1 - 1e-9))) + 1);
    }
    fNodes.clear();
}

void TRestDetectorDriftResponse::SetRange(int axis, double min, double max) {
    if (!(min < max)) {
        throw std::runtime_error("TRestDetectorDriftResponse '" + fName + "': a range must have min < max");
    }
    Range(axis) = {min, max};
    UpdateGrid();
    UpdateYAMLFromParams<TRestDetectorDriftResponse>(fNode);
}

size_t TRestDetectorDriftResponse::GetNumberOfNodes() const {
    return static_cast<size_t>(fGridNodes[0]) * fGridNodes[1] * fGridNodes[2];
}

size_t TRestDetectorDriftResponse::NodeIndex(int i, int j, int k) const {
    return i + static_cast<size_t>(fGridNodes[0]) * (j + static_cast<size_t>(fGridNodes[1]) * k);
}

XYZVector TRestDetectorDriftResponse::GetNodePosition(size_t index) const {
    const int i = index % fGridNodes[0];
    const int j = (index / fGridNodes[0]) % fGridNodes[1];
    const int k = index / (static_cast<size_t>(fGridNodes[0]) * fGridNodes[1]);
    auto coordinate = [](const std::pair<TRestWithUnits, TRestWithUnits>& range, int n, int node) {
        return range.first + (range.second - range.first) * node / (n - 1);
    };
    return {coordinate(fXRange, fGridNodes[0], i), coordinate(fYRange, fGridNodes[1], j),
            coordinate(fZRange, fGridNodes[2], k)};
}

bool TRestDetectorDriftResponse::Contains(const XYZVector& position) const {
    std::array<size_t, 8> nodes;
    std::array<double, 8> weights;
    return Neighbours(position, nodes, weights);
}

std::pair<XYZVector, XYZVector> TRestDetectorDriftResponse::PlaneAxes(const XYZVector& normal) {
    const XYZVector n = normal.Unit();
    const XYZVector axis = std::abs(n.X()) > 0.9 ? XYZVector(0, 1, 0) : XYZVector(1, 0, 0);
    const XYZVector u = (axis - axis.Dot(n) * n).Unit();
    return {u, n.Cross(u)};
}

void TRestDetectorDriftResponse::SetNode(size_t index, int primaries,
                                         const std::vector<std::array<double, 3>>& arrivals) {
    if (fNodes.size() != GetNumberOfNodes()) fNodes.assign(GetNumberOfNodes(), Node());
    Node& node = fNodes.at(index);
    node.fPrimaries = primaries;
    node.fArrived = arrivals.size();
    node.fPositionCounts.assign(static_cast<size_t>(fPositionBins) * fPositionBins, 0);
    node.fTimeCounts.assign(fTimeBins, 0);
    if (arrivals.empty()) return;
    node.fUMin = node.fUMax = arrivals[0][0];
    node.fVMin = node.fVMax = arrivals[0][1];
    node.fTMin = node.fTMax = arrivals[0][2];
    for (const auto& [u, v, t] : arrivals) {
        node.fUMin = std::min(node.fUMin, u);
        node.fUMax = std::max(node.fUMax, u);
        node.fVMin = std::min(node.fVMin, v);
        node.fVMax = std::max(node.fVMax, v);
        node.fTMin = std::min(node.fTMin, t);
        node.fTMax = std::max(node.fTMax, t);
    }
    for (const auto& [u, v, t] : arrivals) {
        const int uBin = Bin(u, node.fUMin, node.fUMax, fPositionBins);
        const int vBin = Bin(v, node.fVMin, node.fVMax, fPositionBins);
        node.fPositionCounts[uBin + static_cast<size_t>(fPositionBins) * vBin]++;
        node.fTimeCounts[Bin(t, node.fTMin, node.fTMax, fTimeBins)]++;
    }
}

bool TRestDetectorDriftResponse::Neighbours(const XYZVector& position, std::array<size_t, 8>& nodes,
                                            std::array<double, 8>& weights) const {
    if (GetNumberOfNodes() == 0 || fNodes.size() != GetNumberOfNodes()) return false;
    const double coordinates[3] = {position.X(), position.Y(), position.Z()};
    const std::pair<TRestWithUnits, TRestWithUnits>* ranges[3] = {&fXRange, &fYRange, &fZRange};
    int first[3];
    double fraction[3];
    for (int axis = 0; axis < 3; axis++) {
        const int n = fGridNodes[axis];
        const double min = ranges[axis]->first, max = ranges[axis]->second;
        const double f = (coordinates[axis] - min) / (max - min) * (n - 1);
        constexpr double tolerance = 1e-9;
        if (!(f >= -tolerance && f <= n - 1 + tolerance)) return false;
        first[axis] = std::clamp(static_cast<int>(std::floor(f)), 0, n - 2);
        fraction[axis] = std::clamp(f - first[axis], 0.0, 1.0);
    }
    for (int corner = 0; corner < 8; corner++) {
        const int di = corner & 1, dj = (corner >> 1) & 1, dk = (corner >> 2) & 1;
        nodes[corner] = NodeIndex(first[0] + di, first[1] + dj, first[2] + dk);
        weights[corner] = (di ? fraction[0] : 1 - fraction[0]) * (dj ? fraction[1] : 1 - fraction[1]) *
                          (dk ? fraction[2] : 1 - fraction[2]);
    }
    return true;
}

TRestDetectorDriftResponse::Summary TRestDetectorDriftResponse::GetNodeSummary(size_t index) const {
    Summary summary;
    const Node& node = fNodes.at(index);
    if (node.fPrimaries > 0) summary.fSurvival = static_cast<double>(node.fArrived) / node.fPrimaries;
    // Mean and standard deviation of counts over [min, max] (bin centres)
    auto moments = [](const std::vector<int>& counts, size_t stride, size_t bins, double min, double max,
                      double& mean, double& sigma) {
        std::vector<double> projected(bins, 0);
        for (size_t i = 0; i < counts.size(); i++) projected[(i / stride) % bins] += counts[i];
        double total = 0, sum = 0, sum2 = 0;
        for (size_t bin = 0; bin < bins; bin++) {
            const double centre = min + (bin + 0.5) * (max - min) / bins;
            total += projected[bin];
            sum += projected[bin] * centre;
            sum2 += projected[bin] * centre * centre;
        }
        if (total <= 0) return;
        mean = sum / total;
        sigma = std::sqrt(std::max(0.0, sum2 / total - mean * mean));
    };
    if (node.fArrived <= 0) return summary;
    const size_t positionBins = fPositionBins, timeBins = fTimeBins;
    moments(node.fPositionCounts, 1, positionBins, node.fUMin, node.fUMax, summary.fMeanU, summary.fSigmaU);
    moments(node.fPositionCounts, positionBins, positionBins, node.fVMin, node.fVMax, summary.fMeanV,
            summary.fSigmaV);
    moments(node.fTimeCounts, 1, timeBins, node.fTMin, node.fTMax, summary.fMeanT, summary.fSigmaT);
    return summary;
}

double TRestDetectorDriftResponse::GetSurvival(const XYZVector& position) const {
    std::array<size_t, 8> nodes;
    std::array<double, 8> weights;
    if (!Neighbours(position, nodes, weights)) return 0;
    double survival = 0;
    for (int corner = 0; corner < 8; corner++) {
        const Node& node = fNodes[nodes[corner]];
        if (node.fPrimaries > 0) survival += weights[corner] * node.fArrived / node.fPrimaries;
    }
    return survival;
}

bool TRestDetectorDriftResponse::SampleArrival(const XYZVector& position, TRandom& random, double& u,
                                               double& v, double& time) const {
    std::array<size_t, 8> nodes;
    std::array<double, 8> weights;
    if (!Neighbours(position, nodes, weights)) return false;
    // Weight of each node: its interpolation weight times its survival
    double total = 0;
    for (int corner = 0; corner < 8; corner++) {
        const Node& node = fNodes[nodes[corner]];
        weights[corner] *= node.fPrimaries > 0 ? static_cast<double>(node.fArrived) / node.fPrimaries : 0;
        total += weights[corner];
    }
    if (total <= 0) return false;
    double target = random.Uniform() * total;
    int corner = 0;
    for (; corner < 7; corner++) {
        if (target < weights[corner]) break;
        target -= weights[corner];
    }
    // The last corners can have no weight (rounding): the last one with arrivals
    while (corner > 0 && weights[corner] <= 0) corner--;
    const Node& node = fNodes[nodes[corner]];

    // Position bin from the 2D counts, then uniform within it
    long count = 0;
    for (const int c : node.fPositionCounts) count += c;
    long bin = 0;
    if (count > 0) {
        long t = static_cast<long>(random.Uniform() * count);
        for (; bin + 1 < static_cast<long>(node.fPositionCounts.size()); bin++) {
            if (t < node.fPositionCounts[bin]) break;
            t -= node.fPositionCounts[bin];
        }
    }
    const int uBin = bin % fPositionBins, vBin = bin / fPositionBins;
    const double uWidth = (node.fUMax - node.fUMin) / fPositionBins;
    const double vWidth = (node.fVMax - node.fVMin) / fPositionBins;
    u = node.fUMin + (uBin + random.Uniform()) * uWidth;
    v = node.fVMin + (vBin + random.Uniform()) * vWidth;
    time = Sample(node.fTimeCounts, node.fTMin, node.fTMax, random);
    return true;
}

void TRestDetectorDriftResponse::Write(TFile* file, const YAML::Node& generator) const {
    if (fNodes.size() != GetNumberOfNodes()) {
        throw std::runtime_error("TRestDetectorDriftResponse: the response has not a value per node");
    }
    YAML::Node metadata;
    // The configuration as read (LoadConfig writes the values in REST units)
    metadata["configuration"] = YAML::Clone(fNode);
    metadata["medium"] = fMedium;
    metadata["pressure"] = fPressure;
    metadata["temperature"] = fTemperature;
    metadata["transportTable"] = fTransportTable;
    metadata["positionBins"] = fPositionBins;
    metadata["timeBins"] = fTimeBins;
    metadata["generator"] = YAML::Clone(generator);
    TRestMetadata::WriteMetadata(file, kName, metadata);

    // LZMA: the counts compress about 20 % better than with the default algorithm (the response file is its
    // own); the counts as 16-bit integers (electrons <= 65535) and the ranges as float: 120 MB were 70 MB
    file->SetCompressionSettings(ROOT::CompressionSettings(ROOT::RCompressionSetting::EAlgorithm::kLZMA, 8));
    file->cd();
    auto tree = std::make_unique<TTree>(kName.c_str(), "Drift response per node of the grid");
    const size_t positionSize = static_cast<size_t>(fPositionBins) * fPositionBins;
    const size_t timeSize = fTimeBins;
    int primaries = 0, arrived = 0;
    float range[6];
    std::vector<UShort_t> positionCounts(positionSize), timeCounts(timeSize);
    tree->Branch("primaries", &primaries, "primaries/I");
    tree->Branch("arrived", &arrived, "arrived/I");
    tree->Branch("range", range, "range[6]/F");
    tree->Branch("positionCounts", positionCounts.data(),
                 ("positionCounts[" + std::to_string(positionSize) + "]/s").c_str());
    tree->Branch("timeCounts", timeCounts.data(), ("timeCounts[" + std::to_string(timeSize) + "]/s").c_str());
    // To draw the response (not read back): the node and its summary, e.g. Draw("meanU:x", "z == 90")
    float x = 0, y = 0, z = 0;
    float summary[7];
    tree->Branch("x", &x, "x/F");
    tree->Branch("y", &y, "y/F");
    tree->Branch("z", &z, "z/F");
    const char* names[7] = {"survival", "meanU", "meanV", "meanT", "sigmaU", "sigmaV", "sigmaT"};
    for (int i = 0; i < 7; i++) tree->Branch(names[i], &summary[i], (std::string(names[i]) + "/F").c_str());
    for (size_t index = 0; index < fNodes.size(); index++) {
        const Node& node = fNodes[index];
        if (node.fPositionCounts.size() != positionSize || node.fTimeCounts.size() != timeSize) {
            throw std::runtime_error("TRestDetectorDriftResponse: node " + std::to_string(index) +
                                     " has not the bins of the response");
        }
        const XYZVector position = GetNodePosition(index);
        x = position.X();
        y = position.Y();
        z = position.Z();
        const Summary values = GetNodeSummary(index);
        const double summaryValues[7] = {values.fSurvival, values.fMeanU,  values.fMeanV, values.fMeanT,
                                         values.fSigmaU,   values.fSigmaV, values.fSigmaT};
        std::copy(summaryValues, summaryValues + 7, summary);
        primaries = node.fPrimaries;
        arrived = node.fArrived;
        const double ranges[6] = {node.fUMin, node.fUMax, node.fVMin, node.fVMax, node.fTMin, node.fTMax};
        std::copy(ranges, ranges + 6, range);
        std::copy(node.fPositionCounts.begin(), node.fPositionCounts.end(), positionCounts.begin());
        std::copy(node.fTimeCounts.begin(), node.fTimeCounts.end(), timeCounts.begin());
        tree->Fill();
    }
    tree->Write(kName.c_str(), TObject::kOverwrite);
}

TRestDetectorDriftResponse TRestDetectorDriftResponse::Read(TFile* file) {
    const std::string error = "TRestDetectorDriftResponse: file " + std::string(file->GetName()) + ": ";
    const YAML::Node metadata = TRestMetadata::ReadMetadata(file, kName);
    if (!metadata || !metadata.IsMap() || !metadata["configuration"]) {
        throw std::runtime_error(error + "no " + kName + " metadata");
    }
    TRestDetectorDriftResponse response("drift", YAML::Clone(metadata["configuration"]));
    response.fMedium = metadata["medium"].as<std::string>("");
    response.fPressure = metadata["pressure"].as<double>(0);
    response.fTemperature = metadata["temperature"].as<double>(0);
    response.fTransportTable = metadata["transportTable"].as<std::string>("");

    auto* tree = file->Get<TTree>(kName.c_str());
    if (!tree) throw std::runtime_error(error + "no " + kName + " tree");
    if (static_cast<size_t>(tree->GetEntries()) != response.GetNumberOfNodes()) {
        throw std::runtime_error(error + "the tree has not an entry per node of the grid");
    }
    // The counts: fixed arrays of 16-bit integers with the bins of the response (files of the earlier format,
    // vectors of int, must be computed again)
    const size_t positionSize = static_cast<size_t>(response.fPositionBins) * response.fPositionBins;
    const size_t timeSize = response.fTimeBins;
    auto checkLeaf = [&](const char* name, const char* type, size_t length) {
        const TLeaf* leaf = tree->GetLeaf(name);
        if (!leaf || std::string(leaf->GetTypeName()) != type ||
            static_cast<size_t>(leaf->GetLen()) != length) {
            throw std::runtime_error(error + "the branch " + name + " is not " + type + "[" +
                                     std::to_string(length) +
                                     "] (a file of an earlier format: compute it again)");
        }
    };
    checkLeaf("range", "Float_t", 6);
    checkLeaf("positionCounts", "UShort_t", positionSize);
    checkLeaf("timeCounts", "UShort_t", timeSize);
    int primaries = 0, arrived = 0;
    float range[6];
    std::vector<UShort_t> positionCounts(positionSize), timeCounts(timeSize);
    tree->SetBranchStatus("*", false);
    for (const char* name : {"primaries", "arrived", "range", "positionCounts", "timeCounts"}) {
        tree->SetBranchStatus(name, true);
    }
    tree->SetBranchAddress("primaries", &primaries);
    tree->SetBranchAddress("arrived", &arrived);
    tree->SetBranchAddress("range", range);
    tree->SetBranchAddress("positionCounts", positionCounts.data());
    tree->SetBranchAddress("timeCounts", timeCounts.data());
    response.fNodes.resize(response.GetNumberOfNodes());
    for (Long64_t entry = 0; entry < tree->GetEntries(); entry++) {
        tree->GetEntry(entry);
        Node& node = response.fNodes[entry];
        node.fPrimaries = primaries;
        node.fArrived = arrived;
        node.fUMin = range[0];
        node.fUMax = range[1];
        node.fVMin = range[2];
        node.fVMax = range[3];
        node.fTMin = range[4];
        node.fTMax = range[5];
        node.fPositionCounts.assign(positionCounts.begin(), positionCounts.end());
        node.fTimeCounts.assign(timeCounts.begin(), timeCounts.end());
    }
    tree->ResetBranchAddresses();
    tree->SetBranchStatus("*", true);
    return response;
}
