#include "FieldMap.h"

#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Garfield/ComponentGrid.hh"
#include "TRestSystemOfUnits.h"

namespace {
// Garfield units: cm, V/cm. REST units: mm, V/mm
constexpr double kCmPerMm = 0.1;
constexpr double kVcmPerVmm = 10.;
/// Vacuum permittivity (C/(V m))
constexpr double kEpsilon0 = 8.8541878128e-12;
/// Relative permittivity above which a node is in a dielectric (1 in the gas; e.g. 3.5 in kapton)
constexpr double kDielectricPermittivity = 1.5;

/// Whether a line (without leading blanks) is a node: it starts with a number (not a line of dashes)
bool IsNode(const std::string& line) {
    if (line.empty()) return false;
    if (line.compare(0, 2, "--") == 0) return false;
    return std::string("0123456789+-.").find(line[0]) != std::string::npos;
}

/// ComponentGrid with nodes out of the gas: in a cell with one of them there is no medium (status -5), the
/// drift ends there. Checked here: the flags of ComponentGrid::LoadElectricField are kept only for 2D maps in
/// Garfield++ 4.0. Each node has its region (CSTRegion: gas, conductor or out of the solved volume,
/// dielectric)
class CSTFieldGrid : public Garfield::ComponentGrid {
   public:
    using Garfield::ComponentGrid::ElectricField;
    using Garfield::ComponentGrid::GetMedium;

    /// \param nodes Nodes along x, y, z; \param low, high Limits of the mesh (cm); \param regions The region
    /// of each node (CSTRegion), x fastest, then y, then z
    CSTFieldGrid(const std::array<unsigned int, 3>& nodes, const std::array<double, 3>& low,
                 const std::array<double, 3>& high, std::vector<unsigned char> regions)
        : fNodes(nodes), fLow(low), fHigh(high), fRegions(std::move(regions)) {}

    void ElectricField(const double x, const double y, const double z, double& ex, double& ey, double& ez,
                       double& v, Garfield::Medium*& m, int& status) override {
        Garfield::ComponentGrid::ElectricField(x, y, z, ex, ey, ez, v, m, status);
        if (status == 0 && !IsActive(x, y, z)) {
            m = nullptr;
            status = -5;
        }
    }
    void ElectricField(const double x, const double y, const double z, double& ex, double& ey, double& ez,
                       Garfield::Medium*& m, int& status) override {
        Garfield::ComponentGrid::ElectricField(x, y, z, ex, ey, ez, m, status);
        if (status == 0 && !IsActive(x, y, z)) {
            m = nullptr;
            status = -5;
        }
    }
    Garfield::Medium* GetMedium(const double x, const double y, const double z) override {
        return IsActive(x, y, z) ? Garfield::ComponentGrid::GetMedium(x, y, z) : nullptr;
    }

    /// The region of the nodes within `reach` cells of the cell of a point (cm) inside the mesh: dielectric
    /// if one of them is, else conductor, else gas
    CSTRegion RegionNear(double x, double y, double z, unsigned int reach) const {
        const double point[3] = {x, y, z};
        unsigned int first[3], last[3];
        for (int a = 0; a < 3; a++) {
            const double step = fNodes[a] > 1 ? (fHigh[a] - fLow[a]) / (fNodes[a] - 1) : 1;
            const double f = std::clamp((point[a] - fLow[a]) / step, 0.0, static_cast<double>(fNodes[a] - 1));
            const auto cell = static_cast<unsigned int>(std::floor(f));
            first[a] = cell > reach ? cell - reach : 0;
            last[a] = std::min(cell + 1 + reach, fNodes[a] - 1);
        }
        unsigned char region = static_cast<unsigned char>(CSTRegion::kGas);
        for (unsigned int k = first[2]; k <= last[2]; k++) {
            for (unsigned int j = first[1]; j <= last[1]; j++) {
                for (unsigned int i = first[0]; i <= last[0]; i++) {
                    const unsigned char node =
                        fRegions[(static_cast<size_t>(k) * fNodes[1] + j) * fNodes[0] + i];
                    if (node == static_cast<unsigned char>(CSTRegion::kDielectric))
                        return CSTRegion::kDielectric;
                    region = std::max(region, node);
                }
            }
        }
        return static_cast<CSTRegion>(region);
    }

   private:
    /// Whether the 8 nodes of the cell of a point (cm) are in the gas; out of the mesh, ComponentGrid decides
    bool IsActive(double x, double y, double z) const {
        const double point[3] = {x, y, z};
        unsigned int first[3], last[3];
        for (int a = 0; a < 3; a++) {
            if (point[a] < fLow[a] || point[a] > fHigh[a]) return true;
            const double step = fNodes[a] > 1 ? (fHigh[a] - fLow[a]) / (fNodes[a] - 1) : 1;
            const auto cell = static_cast<unsigned int>(std::floor((point[a] - fLow[a]) / step));
            first[a] = std::min(cell, fNodes[a] - 1);
            last[a] = std::min(first[a] + 1, fNodes[a] - 1);
        }
        for (const unsigned int k : {first[2], last[2]}) {
            for (const unsigned int j : {first[1], last[1]}) {
                for (const unsigned int i : {first[0], last[0]}) {
                    if (fRegions[(static_cast<size_t>(k) * fNodes[1] + j) * fNodes[0] + i] !=
                        static_cast<unsigned char>(CSTRegion::kGas)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    std::array<unsigned int, 3> fNodes;
    std::array<double, 3> fLow, fHigh;
    std::vector<unsigned char> fRegions;
};

/// Temporary file removed when destroyed
class TemporaryFile {
   public:
    TemporaryFile() {
        const auto directory = std::filesystem::temp_directory_path();
        std::string pattern = (directory / "restGarfield_fieldmap_XXXXXX").string();
        const int descriptor = mkstemp(pattern.data());
        if (descriptor < 0) throw std::runtime_error("restGarfield: cannot create a temporary file");
        close(descriptor);
        fPath = pattern;
    }
    ~TemporaryFile() { std::filesystem::remove(fPath); }
    TemporaryFile(const TemporaryFile&) = delete;
    TemporaryFile& operator=(const TemporaryFile&) = delete;
    const std::string& Path() const { return fPath; }

   private:
    std::string fPath;
};

/// Units of the columns of a CST export: of the coordinates and of the field
struct CSTUnits {
    std::string fLength, fField;
};

/// Reads a CST ASCII export: its header (the units between brackets of x y z and of the 3 field
/// components, before the nodes) and calls onNode(units, x, value) for each node, in the units of the file
template <typename Callback>
void ReadCSTExport(const std::string& fileName, const std::string& error, Callback&& onNode) {
    std::ifstream input(fileName);
    if (!input) throw std::runtime_error(error + "cannot be read");
    CSTUnits units;
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(input, line)) {
        lineNumber++;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos) continue;
        line = line.substr(first);
        if (!IsNode(line)) {
            // The header: the units of the columns between brackets, x y z and the 3 components
            std::vector<std::string> columns;
            const std::regex brackets(R"(\[([^\]]*)\])");
            for (auto it = std::sregex_iterator(line.begin(), line.end(), brackets);
                 it != std::sregex_iterator(); ++it) {
                columns.push_back((*it)[1]);
            }
            if (columns.empty()) continue;
            if (columns.size() < 6 || columns[0] != columns[1] || columns[0] != columns[2] ||
                columns[3] != columns[4] || columns[3] != columns[5]) {
                throw std::runtime_error(error + "line " + std::to_string(lineNumber) +
                                         ": the header must give the units of x, y, z and of the 3 field "
                                         "components");
            }
            units = {columns[0], columns[3]};
            continue;
        }
        if (units.fLength.empty()) {
            throw std::runtime_error(error + "no header with the units before the nodes");
        }
        std::istringstream values(line);
        double x[3], value[3];
        values >> x[0] >> x[1] >> x[2] >> value[0] >> value[1] >> value[2];
        if (values.fail()) {
            throw std::runtime_error(error + "line " + std::to_string(lineNumber) +
                                     ": x y z and 3 field components expected");
        }
        onNode(units, x, value);
    }
}

/// A unit of REST_Units in mm or V/mm, with the error of the file if unknown
double ParseUnit(const std::string& unit, const std::string& error) {
    try {
        return REST_Units::ParseUnit(unit);
    } catch (const std::exception& exception) {
        throw std::runtime_error(error + "unknown units in the header: " + exception.what());
    }
}

/// Reduces a coordinate to the cell [low, high] of a periodic axis (as ComponentGrid::Reduce)
double Reduce(double x, double low, double high, bool periodic, bool mirror, bool& mirrored) {
    const double length = high - low;
    if (!(length > 0) || (!periodic && !mirror)) return x;
    double reduced = low + std::fmod(x - low, length);
    if (reduced < low) reduced += length;
    if (mirror) {
        const long cells = std::lround(std::floor((x - low) / length));
        if (cells % 2 != 0) {
            reduced = low + high - reduced;
            mirrored = true;
        }
    }
    return reduced;
}
}  // namespace

std::unique_ptr<Garfield::ComponentGrid> LoadCSTFieldMap(const std::string& fileName,
                                                         const std::array<double, 3>& offset,
                                                         const std::string& dielectricFile) {
    const std::string error = "restGarfield: field map " + fileName + ": ";

    // The nodes in Garfield units, moved to the detector frame, written to a file for ComponentGrid (it only
    // reads files)
    TemporaryFile grid;
    std::ofstream output(grid.Path());
    output.precision(12);
    double lengthFactor = 0, fieldFactor = 0;  // to cm and V/cm
    // The nodes (cm), whether they are in the gas and the magnitude of their field (V/cm), for the mesh, the
    // inactive nodes and the dielectrics
    std::vector<std::array<double, 3>> positions;
    std::vector<bool> inGas;
    std::vector<float> magnitudes;
    ReadCSTExport(fileName, error, [&](const CSTUnits& units, const double* x, const double* e) {
        if (lengthFactor == 0) {
            lengthFactor = ParseUnit(units.fLength, error) * kCmPerMm;
            fieldFactor = ParseUnit(units.fField, error) * kVcmPerVmm;
        }
        std::array<double, 3> position;
        for (int i = 0; i < 3; i++) position[i] = x[i] * lengthFactor + offset[i] * kCmPerMm;
        for (int i = 0; i < 3; i++) output << " " << position[i];
        for (int i = 0; i < 3; i++) output << " " << e[i] * fieldFactor;
        output << "\n";
        positions.push_back(position);
        // In the gas unless the field is exactly 0: CST writes 0 out of the solved volume (the corners of
        // the box around a cylinder) and inside the conductors
        inGas.push_back(e[0] != 0 || e[1] != 0 || e[2] != 0);
        magnitudes.push_back(std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]) * fieldFactor);
    });
    output.close();
    const size_t nodes = positions.size();
    if (nodes == 0) throw std::runtime_error(error + "no nodes");

    // The mesh: the coordinates of the nodes along each axis, equally spaced, every node given once
    const std::string regular = "a regular grid is needed (every node once, equal steps)";
    std::array<unsigned int, 3> n;
    std::array<double, 3> low, high, step;
    for (int a = 0; a < 3; a++) {
        std::set<double> values;
        for (const auto& position : positions) values.insert(position[a]);
        n[a] = static_cast<unsigned int>(values.size());
        low[a] = *values.begin();
        high[a] = *values.rbegin();
        step[a] = n[a] > 1 ? (high[a] - low[a]) / (n[a] - 1) : 0;
        size_t index = 0;
        for (const double value : values) {
            if (std::abs(value - (low[a] + index * step[a])) > 1e-6 * std::max(1., std::abs(step[a]))) {
                throw std::runtime_error(error + regular);
            }
            index++;
        }
        if (step[a] == 0) step[a] = 1;
    }
    if (static_cast<size_t>(n[0]) * n[1] * n[2] != nodes) throw std::runtime_error(error + regular);
    std::vector<unsigned char> regions(nodes, static_cast<unsigned char>(CSTRegion::kConductor));
    std::vector<bool> given(nodes, false);
    std::vector<float> magnitude(nodes, 0);
    for (size_t node = 0; node < nodes; node++) {
        size_t index = 0;
        for (int a = 2; a >= 0; a--) {
            index = index * n[a] + static_cast<size_t>(std::lround((positions[node][a] - low[a]) / step[a]));
        }
        if (given[index]) throw std::runtime_error(error + regular);
        given[index] = true;
        if (inGas[node]) regions[index] = static_cast<unsigned char>(CSTRegion::kGas);
        magnitude[index] = magnitudes[node];
    }
    positions = {};
    inGas = {};
    magnitudes = {};

    // The dielectrics: nodes where |D| > 1.5 eps0 |E|
    if (!dielectricFile.empty()) {
        const std::string dielectricError =
            "restGarfield: dielectric map " + dielectricFile + " of " + fileName + ": ";
        double dielectricLength = 0;  // to cm
        ReadCSTExport(
            dielectricFile, dielectricError, [&](const CSTUnits& units, const double* x, const double* d) {
                if (dielectricLength == 0) {
                    dielectricLength = ParseUnit(units.fLength, dielectricError) * kCmPerMm;
                    if (units.fField != "C/m^2" && units.fField != "As/m^2") {
                        throw std::runtime_error(dielectricError +
                                                 "the displacement field D in C/m^2 is needed");
                    }
                }
                size_t index = 0;
                for (int a = 2; a >= 0; a--) {
                    const double position = x[a] * dielectricLength + offset[a] * kCmPerMm;
                    const double f = (position - low[a]) / step[a];
                    const long i = std::lround(f);
                    if (i < 0 || i >= static_cast<long>(n[a]) || std::abs(f - i) > 1e-3) {
                        throw std::runtime_error(dielectricError + "a node is not a node of the field map");
                    }
                    index = index * n[a] + static_cast<size_t>(i);
                }
                // |D| / (eps0 |E|): D in C/m^2, E in V/cm (100 V/m)
                const double displacement = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                if (displacement > kDielectricPermittivity * kEpsilon0 * magnitude[index] * 100) {
                    regions[index] = static_cast<unsigned char>(CSTRegion::kDielectric);
                }
            });
    }

    // The mesh set before loading: ComponentGrid sizes its own flags with it (deduced from the file, they
    // were sized 1 x 1 x 1 in Garfield++ 4.0)
    auto component = std::make_unique<CSTFieldGrid>(n, low, high, std::move(regions));
    if (!component->SetMesh(n[0], n[1], n[2], low[0], high[0], low[1], high[1], low[2], high[2]) ||
        !component->LoadElectricField(grid.Path(), "XYZ", false, false)) {
        throw std::runtime_error(error + "rejected by Garfield++ ComponentGrid (" + regular + ")");
    }
    return component;
}

CSTFieldMap::CSTFieldMap() : Garfield::Component("CSTFieldMap") {}

CSTFieldMap::~CSTFieldMap() = default;

void CSTFieldMap::AddMap(std::unique_ptr<Garfield::ComponentGrid> map) {
    std::array<double, 6> limits;
    if (!map || !map->GetElementaryCell(limits[0], limits[2], limits[4], limits[1], limits[3], limits[5])) {
        throw std::runtime_error("restGarfield: a field map without mesh");
    }
    if (fMedium) map->SetMedium(fMedium);
    for (int a = 0; a < 3; a++) {
        fLow[a] = fMaps.empty() ? limits[2 * a] : std::min(fLow[a], limits[2 * a]);
        fHigh[a] = fMaps.empty() ? limits[2 * a + 1] : std::max(fHigh[a], limits[2 * a + 1]);
    }
    fMaps.push_back(std::move(map));
    fLimits.push_back(limits);
    m_ready = true;
}

void CSTFieldMap::SetMedium(Garfield::Medium* medium) {
    fMedium = medium;
    for (auto& map : fMaps) map->SetMedium(medium);
}

int CSTFieldMap::Locate(double x, double y, double z, std::array<double, 3>& cell,
                        std::array<bool, 3>& mirrored) const {
    const double point[3] = {x, y, z};
    for (int a = 0; a < 3; a++) {
        mirrored[a] = false;
        cell[a] = Reduce(point[a], fLow[a], fHigh[a], m_periodic[a], m_mirrorPeriodic[a], mirrored[a]);
    }
    // The first map containing the point; else, within the rounding of the point (e.g. 0.1 mm in cm is above
    // a face at 0.01 cm), the first map with the point moved onto its face (ComponentGrid checks its limits
    // exactly: a point just above them has no medium)
    constexpr double tolerance = 1e-9;  // cm
    for (const double margin : {0.0, tolerance}) {
        for (size_t map = 0; map < fMaps.size(); map++) {
            const auto& limits = fLimits[map];
            bool inside = true;
            for (int a = 0; a < 3 && inside; a++) {
                inside = cell[a] >= limits[2 * a] - margin && cell[a] <= limits[2 * a + 1] + margin;
            }
            if (!inside) continue;
            for (int a = 0; a < 3; a++) cell[a] = std::clamp(cell[a], limits[2 * a], limits[2 * a + 1]);
            return static_cast<int>(map);
        }
    }
    return -1;
}

void CSTFieldMap::ElectricField(const double x, const double y, const double z, double& ex, double& ey,
                                double& ez, Garfield::Medium*& m, int& status) {
    double v = 0;
    ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

void CSTFieldMap::ElectricField(const double x, const double y, const double z, double& ex, double& ey,
                                double& ez, double& v, Garfield::Medium*& m, int& status) {
    std::array<double, 3> cell;
    std::array<bool, 3> mirrored;
    const int map = Locate(x, y, z, cell, mirrored);
    if (map < 0) {
        ex = ey = ez = v = 0;
        m = nullptr;
        status = -6;
        return;
    }
    fMaps[map]->ElectricField(cell[0], cell[1], cell[2], ex, ey, ez, v, m, status);
    if (mirrored[0]) ex = -ex;
    if (mirrored[1]) ey = -ey;
    if (mirrored[2]) ez = -ez;
}

CSTRegion CSTFieldMap::GetRegionNear(double x, double y, double z) const {
    std::array<double, 3> cell;
    std::array<bool, 3> mirrored;
    const int map = Locate(x, y, z, cell, mirrored);
    if (map < 0) return CSTRegion::kOutside;
    const auto* grid = dynamic_cast<const CSTFieldGrid*>(fMaps[map].get());
    return grid ? grid->RegionNear(cell[0], cell[1], cell[2], 1) : CSTRegion::kGas;
}

Garfield::Medium* CSTFieldMap::GetMedium(const double x, const double y, const double z) {
    std::array<double, 3> cell;
    std::array<bool, 3> mirrored;
    const int map = Locate(x, y, z, cell, mirrored);
    return map < 0 ? nullptr : fMaps[map]->GetMedium(cell[0], cell[1], cell[2]);
}

bool CSTFieldMap::GetVoltageRange(double& vmin, double& vmax) {
    vmin = vmax = 0;
    return false;
}

bool CSTFieldMap::GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax, double& ymax,
                                 double& zmax) {
    if (!m_ready) return false;
    double* low[3] = {&xmin, &ymin, &zmin};
    double* high[3] = {&xmax, &ymax, &zmax};
    for (int a = 0; a < 3; a++) {
        const bool periodic = m_periodic[a] || m_mirrorPeriodic[a];
        *low[a] = periodic ? -INFINITY : fLow[a];
        *high[a] = periodic ? INFINITY : fHigh[a];
    }
    return true;
}

bool CSTFieldMap::GetElementaryCell(double& xmin, double& ymin, double& zmin, double& xmax, double& ymax,
                                    double& zmax) {
    if (!m_ready) return false;
    xmin = fLow[0];
    ymin = fLow[1];
    zmin = fLow[2];
    xmax = fHigh[0];
    ymax = fHigh[1];
    zmax = fHigh[2];
    return true;
}

void CSTFieldMap::Reset() {
    fMaps.clear();
    fLimits.clear();
    m_ready = false;
}

void CSTFieldMap::UpdatePeriodicity() {
    for (int a = 0; a < 3; a++) {
        if (m_periodic[a] && m_mirrorPeriodic[a]) {
            throw std::runtime_error("restGarfield: simple and mirror periodicity along the same axis");
        }
    }
    if (m_axiallyPeriodic[0] || m_axiallyPeriodic[1] || m_axiallyPeriodic[2] || m_rotationSymmetric[0] ||
        m_rotationSymmetric[1] || m_rotationSymmetric[2]) {
        throw std::runtime_error("restGarfield: a field map supports only simple and mirror periodicity");
    }
}

std::unique_ptr<CSTFieldMap> LoadCSTFieldMaps(const std::vector<CSTFieldMapFiles>& maps,
                                              const std::array<double, 3>& offset) {
    if (maps.empty()) throw std::runtime_error("restGarfield: no field maps");
    auto component = std::make_unique<CSTFieldMap>();
    for (const auto& map : maps) component->AddMap(LoadCSTFieldMap(map.fField, offset, map.fDielectric));
    return component;
}
