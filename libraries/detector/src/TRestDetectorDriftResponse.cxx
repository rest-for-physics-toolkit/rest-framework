/// \class TRestDetectorDriftResponse
/// Drift of the electrons in a non-uniform electric field (a field map exported from CST Studio Suite) as a
/// response per start point, for special simulations; the drift in a uniform field (the driftField of the
/// medium) is TRestDetectorHitsDriftProcess.
///
/// restGarfield computes it (the `drift` key of a TRestDetectorTransportManager, after the transport and gain
/// tables): it loads the field map in Garfield++ (ComponentGrid), builds the Garfield gas from a transport
/// table of the medium file (without Magboltz) and drifts `electrons` primary electrons (AvalancheMC, steps
/// of `stepSize`, with diffusion and attachment) from each node of a regular grid of start points in the
/// detector frame (xRange, yRange, zRange, gridNodes). The drift ends at the boundary plane, at the distance
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
/// the arrivals of the node, and of the drift time (timeBins). TRestDetectorHitsDriftResponseProcess uses it
/// for each electron of a hit: it survives with the fraction arrived / primaries interpolated between the 8
/// neighbouring nodes (GetSurvival), and its displacement and time are sampled from one of those nodes,
/// chosen with its interpolation weight times its survival (SampleArrival); the displacement is added to the
/// hit projected on the readout, exact for a uniform field.
///
/// Parameters (`drift` key): fieldMap, fieldMapOffset ([0, 0, 0] mm), readouts (YAML file with the
/// TRestDetectorReadoutManager), readout (empty: its data readout), boundary (0 mm), transportTable (key of
/// the table of the medium file; empty: the only one), xRange, yRange, zRange (mm), gridNodes ([2, 2, 2], at
/// least 2 per axis), electrons (1000), stepSize (1 mm), positionBins (20), timeBins (20), responseFile and
/// overwrite (false: an existing responseFile is not computed again). The pressure and temperature are those
/// of the medium. The time of the computation grows with the nodes, the electrons and the drift length over
/// stepSize.
///
///    drift:
///      fieldMap: "drift_field.txt"
///      fieldMapOffset: [0 mm, 0 mm, 0 mm]
///      readouts: "readouts.yaml"
///      boundary: 0.5 mm
///      xRange: [-100 mm, 100 mm]
///      yRange: [-100 mm, 100 mm]
///      zRange: [0 mm, 190 mm]
///      gridNodes: [21, 21, 20]
///      electrons: 1000
///      responseFile: "drift_response.root"
///
/// Response file: the metadata `DriftResponse` in RESTMetadataStore (YAML: `configuration`, the `drift` keys;
/// medium, pressure, temperature, transportTable, positionBins, timeBins and `generator`, the program, date
/// and code that computed it) and the tree `DriftResponse`, an entry per node (x fastest, then y, then z)
/// with primaries, arrived, range (u min, u max, v min, v max, t min, t max; mm and us), positionCounts (u
/// bin + positionBins * v bin) and timeCounts.
///
/// \author JuanAn García <juanangp@unizar.es>

#include "TRestDetectorDriftResponse.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

#include "TFile.h"
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
    reg.RegisterField<D>("gridNodes", &D::fGridNodes);
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
    for (const auto* range : {&fXRange, &fYRange, &fZRange}) {
        if (!(range->first < range->second)) {
            throw std::runtime_error(prefix + "xRange, yRange and zRange must be [min, max] with min < max");
        }
    }
    if (fGridNodes.size() != 3 || *std::min_element(fGridNodes.begin(), fGridNodes.end()) < 2) {
        throw std::runtime_error(prefix + "gridNodes must be 3 values >= 2");
    }
    if (fElectrons <= 0) throw std::runtime_error(prefix + "electrons must be > 0");
    if (fStepSize <= 0) throw std::runtime_error(prefix + "stepSize must be > 0");
    if (fPositionBins <= 0 || fTimeBins <= 0) {
        throw std::runtime_error(prefix + "positionBins and timeBins must be > 0");
    }
    if (fBoundary < 0) throw std::runtime_error(prefix + "boundary must be >= 0");
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
    if (fNodes.size() != GetNumberOfNodes()) return false;
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

    file->cd();
    auto tree = std::make_unique<TTree>(kName.c_str(), "Drift response per node of the grid");
    int primaries = 0, arrived = 0;
    double range[6];
    std::vector<int> positionCounts, timeCounts;
    tree->Branch("primaries", &primaries, "primaries/I");
    tree->Branch("arrived", &arrived, "arrived/I");
    tree->Branch("range", range, "range[6]/D");
    tree->Branch("positionCounts", &positionCounts);
    tree->Branch("timeCounts", &timeCounts);
    for (const Node& node : fNodes) {
        primaries = node.fPrimaries;
        arrived = node.fArrived;
        const double values[6] = {node.fUMin, node.fUMax, node.fVMin, node.fVMax, node.fTMin, node.fTMax};
        std::copy(values, values + 6, range);
        positionCounts = node.fPositionCounts;
        timeCounts = node.fTimeCounts;
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
    int primaries = 0, arrived = 0;
    double range[6];
    std::vector<int>* positionCounts = nullptr;
    std::vector<int>* timeCounts = nullptr;
    tree->SetBranchAddress("primaries", &primaries);
    tree->SetBranchAddress("arrived", &arrived);
    tree->SetBranchAddress("range", range);
    tree->SetBranchAddress("positionCounts", &positionCounts);
    tree->SetBranchAddress("timeCounts", &timeCounts);
    const size_t positionSize = static_cast<size_t>(response.fPositionBins) * response.fPositionBins;
    response.fNodes.resize(response.GetNumberOfNodes());
    for (Long64_t entry = 0; entry < tree->GetEntries(); entry++) {
        tree->GetEntry(entry);
        if (!positionCounts || !timeCounts || positionCounts->size() != positionSize ||
            timeCounts->size() != static_cast<size_t>(response.fTimeBins)) {
            throw std::runtime_error(error + "node " + std::to_string(entry) +
                                     " has not the bins of the grid");
        }
        Node& node = response.fNodes[entry];
        node.fPrimaries = primaries;
        node.fArrived = arrived;
        node.fUMin = range[0];
        node.fUMax = range[1];
        node.fVMin = range[2];
        node.fVMax = range[3];
        node.fTMin = range[4];
        node.fTMax = range[5];
        node.fPositionCounts = *positionCounts;
        node.fTimeCounts = *timeCounts;
    }
    tree->ResetBranchAddresses();
    delete positionCounts;
    delete timeCounts;
    return response;
}
