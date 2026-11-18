/// Unit tests of restGarfield (GoogleTest, BUILD_TESTING): the Garfield++ gas of a medium (CreateGas), and
/// the transport and gain tables computed with small grids (a few field nodes, one Magboltz collision, few
/// avalanches at a low field), written to a medium file of the work directory and read back with the
/// detector library.

#include <TFile.h>
#include <TRestDetectorMedium.h>
#include <TRestDetectorTransportManager.h>
#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

#include "GainSimulation.h"
#include "Garfield/MediumMagboltz.hh"
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

TEST(RestGarfield, GasTransportTableInParallel) {
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
        // At least the primary electron reaches the anode (no attachment)
        EXPECT_GE(table.fGain[i], 1) << i;
        EXPECT_TRUE(std::isfinite(table.fTheta[i])) << i;
    }
    // The avalanche sizes of each field node
    EXPECT_NE(file->Get(("amplification/ArIso/" + key + "/ArIso_" + key + "_size_0").c_str()), nullptr);

    // Only gases
    const auto solid = MakeMedium("{name: Si, state: solid, components: {Si: 1}}", MediumFile("Si"));
    EXPECT_THROW(RunGainSimulation(solid, simulation, false), std::runtime_error);
}
