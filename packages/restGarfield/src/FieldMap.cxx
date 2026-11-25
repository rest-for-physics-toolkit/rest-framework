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

/// Whether a line (without leading blanks) is a node: it starts with a number (not a line of dashes)
bool IsNode(const std::string& line) {
    if (line.empty()) return false;
    if (line.compare(0, 2, "--") == 0) return false;
    return std::string("0123456789+-.").find(line[0]) != std::string::npos;
}

/// ComponentGrid with nodes out of the gas: in a cell with one of them there is no medium (status -5), the
/// drift ends there. Checked here: the flags of ComponentGrid::LoadElectricField are kept only for 2D maps in
/// Garfield++ 4.0
class CSTFieldGrid : public Garfield::ComponentGrid {
   public:
    using Garfield::ComponentGrid::ElectricField;
    using Garfield::ComponentGrid::GetMedium;

    /// \param nodes Nodes along x, y, z; \param low, high Limits of the mesh (cm); \param active Whether
    /// each node is in the gas, x fastest, then y, then z
    CSTFieldGrid(const std::array<unsigned int, 3>& nodes, const std::array<double, 3>& low,
                 const std::array<double, 3>& high, std::vector<bool> active)
        : fNodes(nodes), fLow(low), fHigh(high), fActive(std::move(active)) {}

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
                    if (!fActive[(static_cast<size_t>(k) * fNodes[1] + j) * fNodes[0] + i]) return false;
                }
            }
        }
        return true;
    }

    std::array<unsigned int, 3> fNodes;
    std::array<double, 3> fLow, fHigh;
    std::vector<bool> fActive;
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
}  // namespace

std::unique_ptr<Garfield::ComponentGrid> LoadCSTFieldMap(const std::string& fileName,
                                                         const std::array<double, 3>& offset) {
    const std::string error = "restGarfield: field map " + fileName + ": ";
    std::ifstream input(fileName);
    if (!input) throw std::runtime_error(error + "cannot be read");

    // The nodes in Garfield units, moved to the detector frame, written to a file for ComponentGrid (it only
    // reads files)
    TemporaryFile grid;
    std::ofstream output(grid.Path());
    output.precision(12);
    double lengthFactor = 0, fieldFactor = 0;  // to cm and V/cm
    std::string line;
    size_t lineNumber = 0, nodes = 0;
    // The nodes (cm) and whether they are in the gas, for the mesh and the inactive nodes
    std::vector<std::array<double, 3>> positions;
    std::vector<bool> inGas;
    while (std::getline(input, line)) {
        lineNumber++;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos) continue;
        line = line.substr(first);
        if (!IsNode(line)) {
            // The header: the units of the columns between brackets, x y z Ex Ey Ez
            std::vector<std::string> units;
            const std::regex brackets(R"(\[([^\]]*)\])");
            for (auto it = std::sregex_iterator(line.begin(), line.end(), brackets);
                 it != std::sregex_iterator(); ++it) {
                units.push_back((*it)[1]);
            }
            if (units.empty()) continue;
            if (units.size() < 6 || units[0] != units[1] || units[0] != units[2] || units[3] != units[4] ||
                units[3] != units[5]) {
                throw std::runtime_error(error + "line " + std::to_string(lineNumber) +
                                         ": the header must give the units of x, y, z and Ex, Ey, Ez");
            }
            try {
                lengthFactor = REST_Units::ParseUnit(units[0]) * kCmPerMm;
                fieldFactor = REST_Units::ParseUnit(units[3]) * kVcmPerVmm;
            } catch (const std::exception& unit) {
                throw std::runtime_error(error + "unknown units in the header: " + unit.what());
            }
            continue;
        }
        if (lengthFactor == 0) throw std::runtime_error(error + "no header with the units before the nodes");
        std::istringstream values(line);
        double x[3], e[3];
        values >> x[0] >> x[1] >> x[2] >> e[0] >> e[1] >> e[2];
        if (values.fail()) {
            throw std::runtime_error(error + "line " + std::to_string(lineNumber) +
                                     ": x y z Ex Ey Ez expected");
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
        nodes++;
    }
    output.close();
    if (nodes == 0) throw std::runtime_error(error + "no nodes");

    // The mesh: the coordinates of the nodes along each axis, equally spaced, every node given once
    const std::string regular = "a regular grid is needed (every node once, equal steps)";
    std::array<unsigned int, 3> n;
    std::array<double, 3> low, high;
    for (int a = 0; a < 3; a++) {
        std::set<double> values;
        for (const auto& position : positions) values.insert(position[a]);
        n[a] = static_cast<unsigned int>(values.size());
        low[a] = *values.begin();
        high[a] = *values.rbegin();
        const double step = n[a] > 1 ? (high[a] - low[a]) / (n[a] - 1) : 0;
        size_t index = 0;
        for (const double value : values) {
            if (std::abs(value - (low[a] + index * step)) > 1e-6 * std::max(1., std::abs(step))) {
                throw std::runtime_error(error + regular);
            }
            index++;
        }
    }
    if (static_cast<size_t>(n[0]) * n[1] * n[2] != nodes) throw std::runtime_error(error + regular);
    std::vector<bool> active(nodes, false);
    std::vector<bool> given(nodes, false);
    for (size_t node = 0; node < nodes; node++) {
        size_t index = 0;
        for (int a = 2; a >= 0; a--) {
            const double step = n[a] > 1 ? (high[a] - low[a]) / (n[a] - 1) : 1;
            index = index * n[a] + static_cast<size_t>(std::lround((positions[node][a] - low[a]) / step));
        }
        if (given[index]) throw std::runtime_error(error + regular);
        given[index] = true;
        active[index] = inGas[node];
    }

    // The mesh set before loading: ComponentGrid sizes its own flags with it (deduced from the file, they
    // were sized 1 x 1 x 1 in Garfield++ 4.0)
    auto component = std::make_unique<CSTFieldGrid>(n, low, high, std::move(active));
    if (!component->SetMesh(n[0], n[1], n[2], low[0], high[0], low[1], high[1], low[2], high[2]) ||
        !component->LoadElectricField(grid.Path(), "XYZ", false, false)) {
        throw std::runtime_error(error + "rejected by Garfield++ ComponentGrid (" + regular + ")");
    }
    return component;
}
