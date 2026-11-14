#include "TRestParticles.h"

#include <cctype>
#include <map>
#include <stdexcept>

namespace {

/// Geant4 names and PDG encodings of the usual particles and light nuclei
const std::map<std::string, int>& ParticleTable() {
    static const std::map<std::string, int> table = {
        {"gamma", 22},
        {"e-", 11},
        {"e+", -11},
        {"mu-", 13},
        {"mu+", -13},
        {"tau-", 15},
        {"tau+", -15},
        {"nu_e", 12},
        {"anti_nu_e", -12},
        {"nu_mu", 14},
        {"anti_nu_mu", -14},
        {"nu_tau", 16},
        {"anti_nu_tau", -16},
        {"pi+", 211},
        {"pi-", -211},
        {"pi0", 111},
        {"kaon+", 321},
        {"kaon-", -321},
        {"kaon0", 311},
        {"anti_kaon0", -311},
        {"kaon0L", 130},
        {"kaon0S", 310},
        {"proton", 2212},
        {"anti_proton", -2212},
        {"neutron", 2112},
        {"anti_neutron", -2112},
        {"lambda", 3122},
        {"deuteron", 1000010020},
        {"triton", 1000010030},
        {"He3", 1000020030},
        {"alpha", 1000020040},
    };
    return table;
}

/// Element symbols, index Z - 1
const char* const kElements[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
constexpr int kNElements = sizeof(kElements) / sizeof(kElements[0]);

/// PDG encoding of a nucleus name (Ar40, Ar40[1460.859]); 0 if it is not one
int IonID(const std::string& name) {
    size_t n = 0;
    while (n < name.size() && std::isalpha(static_cast<unsigned char>(name[n]))) n++;
    const std::string symbol = name.substr(0, n);
    size_t m = n;
    while (m < name.size() && std::isdigit(static_cast<unsigned char>(name[m]))) m++;
    if (n == 0 || m == n || m - n > 3) return 0;
    int level = 0;
    if (m < name.size()) {
        // Geant4 excitation energy in keV: "[0.000]" is the ground state, any other an excited one
        if (name[m] != '[' || name.back() != ']') return 0;
        try {
            size_t used = 0;
            const std::string energy = name.substr(m + 1, name.size() - m - 2);
            const double value = std::stod(energy, &used);
            if (used != energy.size() || value < 0) return 0;
            if (value > 0) level = 9;
        } catch (const std::exception&) {
            return 0;
        }
    }
    const int a = std::stoi(name.substr(n, m - n));
    for (int z = 1; z <= kNElements; z++) {
        if (symbol == kElements[z - 1]) {
            if (a < z) return 0;
            return 1000000000 + z * 10000 + a * 10 + level;
        }
    }
    return 0;
}

}  // namespace

namespace REST_Particles {

int GetParticleID(const std::string& name) {
    const auto& table = ParticleTable();
    const auto it = table.find(name);
    if (it != table.end()) return it->second;
    if (const int ion = IonID(name)) return ion;
    // A PDG encoding given as a number
    try {
        size_t used = 0;
        const int id = std::stoi(name, &used);
        if (used == name.size()) return id;
    } catch (const std::exception&) {
    }
    return 0;
}

std::string GetParticleName(int id) {
    if (id == 0) return "";
    for (const auto& [name, code] : ParticleTable()) {
        if (code == id) return name;
    }
    const int z = GetIonZ(id);
    if (z >= 1 && z <= kNElements) {
        const std::string name = kElements[z - 1] + std::to_string(GetIonA(id));
        return id % 10 == 0 ? name : name + "*";
    }
    return std::to_string(id);
}

bool Matches(int id, const std::string& name) {
    if (name == "ion") return IsIon(id);
    const int wanted = GetParticleID(name);
    if (wanted == 0) throw std::invalid_argument("REST_Particles: unknown particle '" + name + "'");
    if (IsIon(wanted)) return IsIon(id) && GetIonZ(id) == GetIonZ(wanted) && GetIonA(id) == GetIonA(wanted);
    return id == wanted;
}

}  // namespace REST_Particles
