#include "FieldMap.h"

#include <unistd.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
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
        for (int i = 0; i < 3; i++) output << " " << (x[i] * lengthFactor + offset[i] * kCmPerMm);
        for (int i = 0; i < 3; i++) output << " " << e[i] * fieldFactor;
        output << "\n";
        nodes++;
    }
    output.close();
    if (nodes == 0) throw std::runtime_error(error + "no nodes");

    auto component = std::make_unique<Garfield::ComponentGrid>();
    if (!component->LoadElectricField(grid.Path(), "XYZ", false, false)) {
        throw std::runtime_error(error + "rejected by Garfield++ ComponentGrid (a regular grid is needed)");
    }
    return component;
}
