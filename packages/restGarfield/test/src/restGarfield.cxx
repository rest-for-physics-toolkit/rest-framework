/// Unit tests of restGarfield (GoogleTest, BUILD_TESTING): the Garfield++ gas of a medium (CreateGas), and
/// the transport and gain tables computed with small grids (a few field nodes, one Magboltz collision, few
/// avalanches at a low field), written to a medium file of the work directory and read back with the
/// detector library; a CST field map (LoadCSTFieldMap) and several with a dielectric and periodicity
/// (LoadCSTFieldMaps), the Garfield gas of a transport table (CreateGasFromTable), the drift response in a
/// uniform field map (ComputeDriftResponse) and the amplification response in a periodic cell
/// (RunAmplificationResponse).

#include <TFile.h>
#include <TRestDetectorAmplificationResponse.h>
#include <TRestDetectorDriftResponse.h>
#include <TRestDetectorMedium.h>
#include <TRestDetectorTransportManager.h>
#include <TRestMicromegasReadout.h>
#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "AmplificationSimulation.h"
#include "DriftSimulation.h"
#include "FieldMap.h"
#include "GainSimulation.h"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/MediumGas.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "TransportSimulation.h"

namespace fs = std::filesystem;

namespace {

/// Medium file of the current test in the work directory, removed (ctest runs the tests in parallel)
std::string MediumFile(const std::string& name) {
    const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
    const fs::path path = fs::path(TESTS_WORK_DIR) / (std::string(test->name()) + "_" + name + ".root");
    fs::remove(path);
    return path.string();
}

/// A medium from its definition, with its medium file
TRestDetectorMedium MakeMedium(const std::string& definition, const std::string& mediumFile) {
    YAML::Node node = YAML::Load(definition);
    node["mediumFile"] = mediumFile;
    return TRestDetectorMedium("medium", node);
}

const std::string kArIso = "{name: ArIso, state: gas, components: {Ar: 99, iC4H10: 1}}";

/// A Micromegas readout with its channels in x in [-3, 4] mm, y in [-2, 5] mm (SetGridOnReadout)
class ChannelsBoxReadout : public TRestMicromegasReadout {
   public:
    using TRestMicromegasReadout::TRestMicromegasReadout;
    std::pair<ROOT::Math::XYZVector, ROOT::Math::XYZVector> GetChannelsBoundingBox() const override {
        return {{-3, -2, 0}, {4, 5, 0.1}};
    }
};

}  // namespace

TEST(RestGarfield, GasOfTheMedium) {
    const auto medium = MakeMedium(kArIso, MediumFile("ArIso"));
    const auto gas = CreateGas(medium, 1.4, 300);
    ASSERT_NE(gas, nullptr);
    EXPECT_NEAR(gas->GetPressure(), 1.4 * 750.062, 1e-3);  // Torr
    EXPECT_NEAR(gas->GetTemperature(), 300, 1e-9);

    // Errors: not a gas, no components, a component unknown to Garfield
    EXPECT_THROW(
        CreateGas(MakeMedium("{name: Si, state: solid, components: {Si: 1}}", MediumFile("Si")), 1, 293.15),
        std::runtime_error);
    EXPECT_THROW(CreateGas(MakeMedium("{name: Empty, state: gas}", MediumFile("Empty")), 1, 293.15),
                 std::runtime_error);
    const auto unknown = MakeMedium("{name: X, state: gas, components: {Unobtainium: 1}}", MediumFile("X"));
    EXPECT_THROW(CreateGas(unknown, 1, 293.15), std::runtime_error);
}

TEST(RestGarfield, SiliconTransportTable) {
    const std::string fileName = MediumFile("Si");
    const auto medium = MakeMedium("{name: Si, state: solid, components: {Si: 1}}", fileName);
    TRestDetectorTransportTable simulation("transport",
                                           YAML::Load("{fieldRange: [100 V/cm, 1000 V/cm], fieldNodes: 3}"));
    RunTransportSimulation(medium, simulation, false);

    // One table with electrons and holes, not scaled with the pressure (a solid)
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    ASSERT_TRUE(file && !file->IsZombie());
    const YAML::Node tables = TRestDetectorTransportTable::GetTables(file.get(), "Si");
    ASSERT_EQ(tables.size(), 1u);
    const std::string key = tables[0]["key"].as<std::string>();
    const auto table = TRestDetectorTransportTable::Read(file.get(), "Si", key);
    EXPECT_EQ(table.fFields.size(), 3u);
    EXPECT_FALSE(table.fPressureScaling);
    ASSERT_TRUE(table.HasCarrier("electron"));
    ASSERT_TRUE(table.HasCarrier("hole"));
    // Low field mobilities of about 1400 and 450 cm2/V/s: 1.4 and 0.45 mm/us at 100 V/cm
    EXPECT_NEAR(table.GetDriftVelocity(10, 1), 1.4, 0.4);
    EXPECT_NEAR(table.GetDriftVelocity(10, 1, "hole"), 0.45, 0.2);
    EXPECT_GT(table.GetDriftVelocity(100, 1), table.GetDriftVelocity(10, 1));

    // The medium reads the table of its file
    const auto withTables = MakeMedium("{name: Si, state: solid, charge: {driftField: 500 V/cm}}", fileName);
    const auto& charge = withTables.GetChargeResponse("test");
    EXPECT_TRUE(charge.HasTransportTables());
    EXPECT_NEAR(charge.GetDriftVelocity(), table.GetDriftVelocity(50, 1), 1e-9);

    // usePenning is only for gases
    TRestDetectorTransportTable penning(
        "transport", YAML::Load("{fieldRange: [100 V/cm, 1000 V/cm], usePenning: true, penningRate: 0.3}"));
    EXPECT_THROW(RunTransportSimulation(medium, penning, false), std::runtime_error);
}

// Magboltz in full (about a minute): a test suite named *Long, skipped with ctest -E Long
TEST(RestGarfieldLong, GasTransportTableInParallel) {
    const std::string fileName = MediumFile("ArIso");
    const auto medium = MakeMedium(kArIso, fileName);
    TRestDetectorTransportTable simulation(
        "transport", YAML::Load("{fieldRange: [100 V/cm, 1000 V/cm], fieldNodes: 4, nCollisions: 1}"));
    // Two processes of two field nodes each
    RunTransportSimulation(medium, simulation, false, 2);

    std::string key;
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
        ASSERT_TRUE(file && !file->IsZombie());
        const YAML::Node tables = TRestDetectorTransportTable::GetTables(file.get(), "ArIso");
        ASSERT_EQ(tables.size(), 1u);
        key = tables[0]["key"].as<std::string>();
        EXPECT_EQ(key, "1bar_293.15K_E100-1000Vcm");
        const auto table = TRestDetectorTransportTable::Read(file.get(), "ArIso", key);
        ASSERT_EQ(table.fFields.size(), 4u);
        EXPECT_TRUE(table.fPressureScaling);
        const auto& electron = table.GetCarrier("electron");
        for (size_t i = 0; i < 4; i++) {
            EXPECT_GT(electron.fDriftVelocity[i], 0) << i;
            EXPECT_GT(electron.fLongitudinalDiffusion[i], 0) << i;
            EXPECT_GT(electron.fTransversalDiffusion[i], 0) << i;
        }
    }

    // A table already in the file is not computed again (the file is not changed)
    const auto modified = fs::last_write_time(fileName);
    RunTransportSimulation(medium, simulation, false, 2);
    EXPECT_EQ(fs::last_write_time(fileName), modified);
}

TEST(RestGarfield, GainTableOfFewAvalanches) {
    const std::string fileName = MediumFile("ArIso");
    const auto medium = MakeMedium(kArIso, fileName);
    // Garfield++ seeds its generator from the time (TRandom3(0)): a fixed seed for a reproducible test
    Garfield::randomEngine.Seed(1);
    // A low field (small gain): fast avalanches
    TRestDetectorGainTable simulation(
        "gain", YAML::Load("{gap: 128 um, reducedFieldRange: [10 kV/cm/bar, 12 kV/cm/bar], fieldNodes: 2, "
                           "avalanches: 20}"));
    RunGainSimulation(medium, simulation, false);

    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    ASSERT_TRUE(file && !file->IsZombie());
    const YAML::Node tables = TRestDetectorGainTable::GetTables(file.get(), "ArIso");
    ASSERT_EQ(tables.size(), 1u);
    const std::string key = tables[0]["key"].as<std::string>();
    EXPECT_EQ(key, "1bar_293.15K_gap128um_E10000-12000Vcm");
    EXPECT_EQ(tables[0]["avalanches"].as<int>(), 20);
    const auto table = TRestDetectorGainTable::Read(file.get(), "ArIso", key);
    ASSERT_EQ(table.fGain.size(), 2u);
    EXPECT_NEAR(table.fTableGap, 0.128, 1e-9);
    for (size_t i = 0; i < 2; i++) {
        // A gain of about one: the primary electron, started at the mesh, may diffuse back out of the gap
        EXPECT_GT(table.fGain[i], 0.5) << i;
        EXPECT_TRUE(std::isfinite(table.fTheta[i])) << i;
    }
    // The avalanche sizes of each field node
    EXPECT_NE(file->Get(("amplification/ArIso/" + key + "/ArIso_" + key + "_size_0").c_str()), nullptr);

    // Only gases
    const auto solid = MakeMedium("{name: Si, state: solid, components: {Si: 1}}", MediumFile("Si"));
    EXPECT_THROW(RunGainSimulation(solid, simulation, false), std::runtime_error);
}

TEST(RestGarfield, CSTFieldMap) {
    // A CST ASCII export (header with the units, dashes, x fastest): 3 x 3 x 3 nodes in mm, Ez in V/m growing
    // with z (1e5 V/m at z = 0, 2e5 V/m at z = 10 mm); the field is 0 at x = y = 10 mm (out of the solved
    // volume of CST, e.g. a corner of the box around a cylinder)
    const fs::path cstFile = fs::path(TESTS_WORK_DIR) / "CSTFieldMap_efield.txt";
    {
        std::ofstream file(cstFile);
        file << "           x [mm]           y [mm]           z [mm]         Ex [V/m]         Ey [V/m]"
                "         Ez [V/m]\n";
        file << std::string(102, '-') << "\n";
        for (const double z : {0., 5., 10.}) {
            for (const double y : {-10., 0., 10.}) {
                for (const double x : {-10., 0., 10.}) {
                    const double ez = x == 10 && y == 10 ? 0 : 1e5 * (1 + z / 10);
                    file << x << " " << y << " " << z << " 0 0 " << ez << "\n";
                }
            }
        }
    }
    // The origin of the CST frame at (1, 2, 3) mm in the detector frame
    const auto grid = LoadCSTFieldMap(cstFile.string(), {1, 2, 3});
    Garfield::MediumGas gas;
    grid->SetMedium(&gas);
    unsigned int nx = 0, ny = 0, nz = 0;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    ASSERT_TRUE(grid->GetMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax));
    EXPECT_EQ(nx, 3u);
    EXPECT_EQ(ny, 3u);
    EXPECT_EQ(nz, 3u);
    // cm
    EXPECT_NEAR(xmin, -0.9, 1e-9);
    EXPECT_NEAR(xmax, 1.1, 1e-9);
    EXPECT_NEAR(ymin, -0.8, 1e-9);
    EXPECT_NEAR(zmin, 0.3, 1e-9);
    EXPECT_NEAR(zmax, 1.3, 1e-9);

    // Between nodes, at z = 2.5 mm of the CST frame: 1.25e5 V/m = 1250 V/cm
    double ex = 0, ey = 0, ez = 0;
    Garfield::Medium* medium = nullptr;
    int status = 0;
    grid->ElectricField(0.1 + 0.05, 0.2 - 0.3, 0.3 + 0.25, ex, ey, ez, medium, status);
    EXPECT_EQ(status, 0);
    EXPECT_NEAR(ex, 0, 1e-9);
    EXPECT_NEAR(ey, 0, 1e-9);
    EXPECT_NEAR(ez, 1250, 1e-6);
    // Out of the map: z = 15 mm in the CST frame
    grid->ElectricField(0.1, 0.2, 0.3 + 1.5, ex, ey, ez, medium, status);
    EXPECT_EQ(status, -11);
    // A cell with a node of field 0 is out of the gas: no medium, at x = y = 5 mm of the CST frame
    grid->ElectricField(0.1 + 0.5, 0.2 + 0.5, 0.3 + 0.25, ex, ey, ez, medium, status);
    EXPECT_EQ(status, -5);
    EXPECT_EQ(medium, nullptr);

    // Errors: no file, no header with the units, unknown units, a wrong node line, a node missing
    EXPECT_THROW(LoadCSTFieldMap((fs::path(TESTS_WORK_DIR) / "CSTFieldMap_none.txt").string()),
                 std::runtime_error);
    auto write = [](const fs::path& path, const std::string& content) {
        std::ofstream file(path);
        file << content;
        return path.string();
    };
    const fs::path bad = fs::path(TESTS_WORK_DIR) / "CSTFieldMap_bad.txt";
    EXPECT_THROW(LoadCSTFieldMap(write(bad, "0 0 0 0 0 1\n")), std::runtime_error);
    const std::string unknown = "x [parsec] y [parsec] z [parsec] Ex [V/m] Ey [V/m] Ez [V/m]\n0 0 0 0 0 1\n";
    EXPECT_THROW(LoadCSTFieldMap(write(bad, unknown)), std::runtime_error);
    const std::string wrongNode = "x [mm] y [mm] z [mm] Ex [V/m] Ey [V/m] Ez [V/m]\n0 0 zero 0 0 1\n";
    EXPECT_THROW(LoadCSTFieldMap(write(bad, wrongNode)), std::runtime_error);
    const std::string missing =
        "x [mm] y [mm] z [mm] Ex [V/m] Ey [V/m] Ez [V/m]\n0 0 0 0 0 1\n1 0 0 0 0 1\n"
        "0 1 0 0 0 1\n";
    EXPECT_THROW(LoadCSTFieldMap(write(bad, missing)), std::runtime_error);
}

TEST(RestGarfield, CSTFieldMapsDielectricAndPeriodicity) {
    // A fine map (x, y in [-1, 1] mm, z in [0, 2] mm, steps of 1 mm; Ex = 1e4 V/m per mm of x, Ez = 1e5 V/m)
    // with the displacement field D of its layers z = 0 and 1 mm (eps0 E in the gas, 3.5 eps0 E at the node
    // x = y = 1 mm, z = 0: a dielectric), and a coarse map above it (z in [2, 6] mm, Ez = 2e5 V/m)
    const double eps0 = 8.8541878128e-12;
    auto write = [](const std::string& name, const std::string& unit, const std::vector<double>& zs,
                    auto value) {
        const fs::path path = fs::path(TESTS_WORK_DIR) / ("CSTFieldMaps_" + name + ".txt");
        std::ofstream file(path);
        file << "x [mm] y [mm] z [mm] X [" << unit << "] Y [" << unit << "] Z [" << unit << "]\n";
        file << std::string(102, '-') << "\n";
        for (const double z : zs) {
            for (const double y : {-1., 0., 1.}) {
                for (const double x : {-1., 0., 1.}) {
                    const auto [vx, vy, vz] = value(x, y, z);
                    file << x << " " << y << " " << z << " " << vx << " " << vy << " " << vz << "\n";
                }
            }
        }
        return path.string();
    };
    auto fineField = [](double x, double, double) { return std::array<double, 3>{1e4 * x, 0, 1e5}; };
    const std::string fine = write("fine", "V/m", {0, 1, 2}, fineField);
    const std::string dielectric = write("D", "C/m^2", {0, 1}, [&](double x, double y, double z) {
        const double permittivity = x == 1 && y == 1 && z == 0 ? 3.5 : 1;
        const auto e = fineField(x, y, z);
        return std::array<double, 3>{permittivity * eps0 * e[0], 0, permittivity * eps0 * e[2]};
    });
    const std::string coarse = write("coarse", "V/m", {2, 4, 6},
                                     [](double, double, double) { return std::array<double, 3>{0, 0, 2e5}; });

    const auto map = LoadCSTFieldMaps({{fine, dielectric}, {coarse, ""}});
    EXPECT_EQ(map->GetNumberOfMaps(), 2u);
    Garfield::MediumGas gas;
    map->SetMedium(&gas);
    double ex = 0, ey = 0, ez = 0;
    Garfield::Medium* medium = nullptr;
    int status = 0;
    // The fine map below z = 2 mm and on it (it takes precedence), the coarse one above (cm, V/cm)
    map->ElectricField(-0.05, -0.05, 0.15, ex, ey, ez, medium, status);
    EXPECT_EQ(status, 0);
    EXPECT_NEAR(ex, -50, 1e-6);
    EXPECT_NEAR(ez, 1000, 1e-6);
    map->ElectricField(-0.05, -0.05, 0.2, ex, ey, ez, medium, status);
    EXPECT_NEAR(ez, 1000, 1e-6);
    map->ElectricField(-0.05, -0.05, 0.3, ex, ey, ez, medium, status);
    EXPECT_EQ(status, 0);
    EXPECT_NEAR(ez, 2000, 1e-6);
    // Out of both maps
    map->ElectricField(0, 0, 0.7, ex, ey, ez, medium, status);
    EXPECT_EQ(status, -6);
    EXPECT_EQ(medium, nullptr);
    // The cell with the dielectric node has no medium; not periodic, x = 1.5 mm is out of the maps
    map->ElectricField(0.05, 0.05, 0.05, ex, ey, ez, medium, status);
    EXPECT_EQ(status, -5);
    EXPECT_EQ(map->GetMedium(0.05, 0.05, 0.05), nullptr);
    EXPECT_NE(map->GetMedium(-0.05, 0.05, 0.05), nullptr);
    EXPECT_EQ(map->GetMedium(0.15, -0.05, 0.05), nullptr);

    // The region around a point: the dielectric node within a cell, the gas of the coarse map, out of both
    EXPECT_EQ(map->GetRegionNear(0.095, 0.095, 0.005), CSTRegion::kDielectric);
    EXPECT_EQ(map->GetRegionNear(0, 0, 0.4), CSTRegion::kGas);
    EXPECT_EQ(map->GetRegionNear(0, 0, 0.7), CSTRegion::kOutside);

    // Periodic in x and y (period 2 mm): x = 2.5 mm is x = 0.5 mm
    map->EnablePeriodicityX();
    map->EnablePeriodicityY();
    map->ElectricField(0.25, -0.05, 0.05, ex, ey, ez, medium, status);
    EXPECT_EQ(status, 0);
    EXPECT_NEAR(ex, 50, 1e-6);
    map->ElectricField(0.25, 0.25, 0.05, ex, ey, ez, medium, status);
    EXPECT_EQ(status, -5);
    // Mirror periodic in x: x = 1.5 mm is x = 0.5 mm mirrored, Ex changes sign
    map->EnablePeriodicityX(false);
    map->EnableMirrorPeriodicityX();
    map->ElectricField(0.15, -0.05, 0.05, ex, ey, ez, medium, status);
    EXPECT_EQ(status, 0);
    EXPECT_NEAR(ex, -50, 1e-6);
    EXPECT_NEAR(ez, 1000, 1e-6);

    // Errors: a node of D out of the field map, D not in C/m^2, no maps
    const std::string outside = write("D_outside", "C/m^2", {0.5}, [](double, double, double) {
        return std::array<double, 3>{0, 0, 1e-6};
    });
    EXPECT_THROW(LoadCSTFieldMap(fine, {0, 0, 0}, outside), std::runtime_error);
    const std::string units =
        write("D_units", "V/m", {0}, [](double, double, double) { return std::array<double, 3>{0, 0, 1}; });
    EXPECT_THROW(LoadCSTFieldMap(fine, {0, 0, 0}, units), std::runtime_error);
    EXPECT_THROW(LoadCSTFieldMaps({}), std::runtime_error);
}

TEST(RestGarfield, GasFromTransportTable) {
    const auto medium = MakeMedium(kArIso, MediumFile("ArIso"));
    // A table at 1 bar (no Magboltz): fields 10, 50, 100 V/mm
    TRestDetectorTransportTable table;
    table.fMedium = "ArIso";
    table.fPressure = 1;
    table.fFields = {10, 50, 100};
    auto& electron = table.fCarriers["electron"];
    electron.fDriftVelocity = {10, 30, 50};                 // mm/us
    electron.fLongitudinalDiffusion = {0.03, 0.025, 0.02};  // sqrt(cm)
    electron.fTransversalDiffusion = {0.02, 0.015, 0.01};
    electron.fTownsend = {0, 1, 100};  // 1/cm
    electron.fAttachment = {0.1, 0.2, 0.3};

    // At 2 bar (E/p scaling): the nodes at 20, 100, 200 V/mm
    const auto gas = CreateGasFromTable(medium, table, 2);
    ASSERT_NE(gas, nullptr);
    double vx = 0, vy = 0, vz = 0, dl = 0, dt = 0, alpha = 0, eta = 0;
    // 100 V/mm = 1000 V/cm, the second node: the values of the table at E/p = 50 V/mm/bar
    ASSERT_TRUE(gas->ElectronVelocity(0, 0, 1000, 0, 0, 0, vx, vy, vz));
    EXPECT_NEAR(std::sqrt(vx * vx + vy * vy + vz * vz), 30e-4, 1e-9);  // cm/ns
    EXPECT_LT(vz, 0);                                                  // electrons against the field
    ASSERT_TRUE(gas->ElectronDiffusion(0, 0, 1000, 0, 0, 0, dl, dt));
    EXPECT_NEAR(dl, 0.025 / std::sqrt(2), 1e-9);
    EXPECT_NEAR(dt, 0.015 / std::sqrt(2), 1e-9);
    ASSERT_TRUE(gas->ElectronTownsend(0, 0, 1000, 0, 0, 0, alpha));
    EXPECT_NEAR(alpha, 2, 1e-9);
    ASSERT_TRUE(gas->ElectronAttachment(0, 0, 2000, 0, 0, 0, eta));
    EXPECT_NEAR(eta, 0.6, 1e-9);
    // A Townsend coefficient 0 stays 0
    ASSERT_TRUE(gas->ElectronTownsend(0, 0, 200, 0, 0, 0, alpha));
    EXPECT_NEAR(alpha, 0, 1e-12);

    // Errors: a table without electrons, a medium that is not a gas
    TRestDetectorTransportTable noElectrons = table;
    noElectrons.fCarriers.clear();
    EXPECT_THROW(CreateGasFromTable(medium, noElectrons, 1), std::runtime_error);
    const auto solid = MakeMedium("{name: Si, state: solid, components: {Si: 1}}", MediumFile("Si"));
    EXPECT_THROW(CreateGasFromTable(solid, table, 1), std::runtime_error);
}

TEST(RestGarfield, AmplificationResponseInAPeriodicCell) {
    const auto medium = MakeMedium(kArIso, MediumFile("ArIso"));
    Garfield::randomEngine.Seed(1);
    // A cell of 0.2 x 0.2 mm, periodic in x and y, with a uniform field of 10 kV/cm along +z (a small gain:
    // fast avalanches) above an anode at z = 0 (nodes with the field 0: the cells below z = 0.05 mm stop the
    // electrons)
    const fs::path cstFile = fs::path(TESTS_WORK_DIR) / "AmplificationResponse_efield.txt";
    {
        std::ofstream file(cstFile);
        file << "x [um] y [um] z [um] Ex [V/m] Ey [V/m] Ez [V/m]\n" << std::string(60, '-') << "\n";
        for (const double z : {0., 50., 100., 150., 200.}) {
            for (const double y : {-100., 0., 100.}) {
                for (const double x : {-100., 0., 100.}) {
                    file << x << " " << y << " " << z << " 0 0 " << (z == 0 ? 0 : 1e6) << "\n";
                }
            }
        }
    }
    const std::string responseFile = (fs::path(TESTS_WORK_DIR) / "AmplificationResponse.root").string();
    fs::remove(responseFile);
    const std::string config = "{fieldMaps: [\"" + cstFile.string() + "\"], responseFile: \"" + responseFile +
                               "\", start: 150 um, anode: 60 um, electrons: 10, sizeBins: 5, timeBins: 5}";
    TRestDetectorAmplificationResponse response("amplification", YAML::Load(config));
    RunAmplificationResponse(medium, response, false, 2);
    EXPECT_EQ(response.fMedium, "ArIso");
    EXPECT_EQ(response.fPrimaries, 10);
    // The primaries reach the anode (no mesh, no attachment) unless they diffuse out of the top of the map
    EXPECT_GE(response.fCollected, 5);
    EXPECT_GE(response.fMeanSize, 1);
    // Every lost primary counted somewhere
    const auto& lost = response.fLosses;
    EXPECT_EQ(lost.fOut + lost.fConductor + lost.fDielectric + lost.fAttached + lost.fOther,
              response.fPrimaries - response.fCollected);
    // About 0.1 mm at several cm/us
    EXPECT_GT(response.fTMin, 0);
    EXPECT_LT(response.fTMax, 0.1);

    std::unique_ptr<TFile> file(TFile::Open(responseFile.c_str(), "READ"));
    ASSERT_TRUE(file && !file->IsZombie());
    const auto read = TRestDetectorAmplificationResponse::Read(file.get());
    EXPECT_EQ(read.fCollected, response.fCollected);
    EXPECT_EQ(read.fSizeCounts, response.fSizeCounts);

    // The start plane above the maps
    TRestDetectorAmplificationResponse above(
        "amplification", YAML::Load("{fieldMaps: [\"" + cstFile.string() +
                                    "\"], responseFile: r.root, start: 1 mm, anode: 60 um, electrons: 1}"));
    EXPECT_THROW(ComputeAmplificationResponse(medium, above, false), std::runtime_error);
}

TEST(RestGarfield, DriftResponseInAUniformField) {
    const auto medium = MakeMedium(kArIso, MediumFile("ArIso"));
    // A uniform field of 100 V/mm along +z (electrons drift to -z) in x, y in [-10, 10] mm, z in [0, 20] mm
    const fs::path cstFile = fs::path(TESTS_WORK_DIR) / "DriftResponse_efield.txt";
    {
        std::ofstream file(cstFile);
        file << "x [mm] y [mm] z [mm] Ex [V/m] Ey [V/m] Ez [V/m]\n" << std::string(60, '-') << "\n";
        for (const double z : {0., 5., 10., 15., 20.}) {
            for (const double y : {-10., 0., 10.}) {
                for (const double x : {-10., 0., 10.}) file << x << " " << y << " " << z << " 0 0 1e5\n";
            }
        }
    }
    // Electrons at 10 mm/us without diffusion, attachment or multiplication
    TRestDetectorTransportTable table;
    table.fMedium = "ArIso";
    table.fPressure = medium.GetPressure();
    table.fFields = {50, 100, 200};  // V/mm
    auto& electron = table.fCarriers["electron"];
    electron.fDriftVelocity = {10, 10, 10};  // mm/us
    electron.fLongitudinalDiffusion = {0, 0, 0};
    electron.fTransversalDiffusion = {0, 0, 0};
    electron.fTownsend = {0, 0, 0};
    electron.fAttachment = {0, 0, 0};
    // A readout plane at z = 0 with the normal +z (no decoding needed for the distance)
    const TRestMicromegasReadout readout(
        "mm", YAML::Load("{readoutParameters: {nChannels: 4, pitch: 1 mm, thickness: 0.1 mm, "
                         "gasThickness: 50 mm, normal: [0, 0, 1]}}"));

    // Nodes at z = 0, 5, 10, 15 mm, boundary plane at z = 1 mm
    const std::string config = "{fieldMap: \"" + cstFile.string() +
                               "\", responseFile: r.root, boundary: 1 mm, xRange: [-5 mm, 5 mm], "
                               "yRange: [-5 mm, 5 mm], zRange: [0 mm, 15 mm], "
                               "gridStep: [10 mm, 10 mm, 5 mm], electrons: 20, stepSize: 1 mm, "
                               "positionBins: 4, timeBins: 4}";
    TRestDetectorDriftResponse response("drift", YAML::Load(config));
    ComputeDriftResponse(medium, table, readout, response, false);
    EXPECT_EQ(response.fMedium, "ArIso");
    ASSERT_EQ(response.fNodes.size(), 16u);
    for (size_t index = 0; index < response.fNodes.size(); index++) {
        const auto& node = response.fNodes[index];
        const double z = response.GetNodePosition(index).Z();
        // Every electron reaches the boundary plane straight down, at (z - 1 mm) / v (0 below the plane)
        EXPECT_EQ(node.fPrimaries, 20) << "z = " << z;
        EXPECT_EQ(node.fArrived, 20) << "z = " << z;
        EXPECT_NEAR(node.fUMin, 0, 1e-6);
        EXPECT_NEAR(node.fUMax, 0, 1e-6);
        EXPECT_NEAR(node.fVMin, 0, 1e-6);
        EXPECT_NEAR(node.fVMax, 0, 1e-6);
        const double time = z > 1 ? (z - 1) / 10 : 0;  // us
        EXPECT_NEAR(node.fTMin, time, 1e-3) << "z = " << z;
        EXPECT_NEAR(node.fTMax, time, 1e-3) << "z = " << z;
    }
    EXPECT_NEAR(response.GetSurvival({1, 2, 7}), 1, 1e-12);

    // In 3 processes (every third node each): the same response (no diffusion)
    TRestDetectorDriftResponse parallel("drift", YAML::Load(config));
    ComputeDriftResponse(medium, table, readout, parallel, false, 3);
    ASSERT_EQ(parallel.fNodes.size(), response.fNodes.size());
    for (size_t index = 0; index < parallel.fNodes.size(); index++) {
        EXPECT_EQ(parallel.fNodes[index].fArrived, response.fNodes[index].fArrived) << index;
        EXPECT_NEAR(parallel.fNodes[index].fTMin, response.fNodes[index].fTMin, 1e-9) << index;
        EXPECT_EQ(parallel.fNodes[index].fTimeCounts, response.fNodes[index].fTimeCounts) << index;
    }

    // The boundary plane above the field map: no drift region
    TRestDetectorDriftResponse high("drift", YAML::Load(config));
    high.fBoundary = 30;
    EXPECT_THROW(ComputeDriftResponse(medium, table, readout, high, false), std::runtime_error);

    // The grid on the readout: x and y not given, the extent of the channels (x in [-3, 4] mm, y in [-2, 5]
    // mm) plus the margin; the range along the normal (z) is needed
    const ChannelsBoxReadout channels(
        "mm", YAML::Load("{readoutParameters: {nChannels: 4, pitch: 1 mm, thickness: 0.1 mm, "
                         "gasThickness: 50 mm, normal: [0, 0, 1]}}"));
    const std::string onReadout = "{fieldMap: \"" + cstFile.string() +
                                  "\", responseFile: r.root, boundary: 1 mm, zRange: [0 mm, 15 mm], "
                                  "margin: 1 mm, gridStep: [3 mm, 3 mm, 5 mm], electrons: 1}";
    TRestDetectorDriftResponse grid("drift", YAML::Load(onReadout));
    EXPECT_THROW(ComputeDriftResponse(medium, table, channels, grid, false), std::runtime_error);
    SetGridOnReadout(grid, channels);
    EXPECT_DOUBLE_EQ(grid.fXRange.first, -4);
    EXPECT_DOUBLE_EQ(grid.fXRange.second, 5);
    EXPECT_DOUBLE_EQ(grid.fYRange.first, -3);
    EXPECT_DOUBLE_EQ(grid.fYRange.second, 6);
    EXPECT_EQ(grid.GetNumberOfNodes(), 4u * 4u * 4u);  // 9 mm in 3 steps of 3 mm, 15 mm in 3 of 5 mm
    TRestDetectorDriftResponse noZ("drift", YAML::Load("{fieldMap: f.txt, responseFile: r.root, "
                                                       "gridStep: [3 mm, 3 mm, 5 mm]}"));
    EXPECT_THROW(SetGridOnReadout(noZ, channels), std::runtime_error);
}
