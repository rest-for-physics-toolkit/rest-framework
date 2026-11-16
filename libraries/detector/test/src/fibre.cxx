/// Unit tests of the fibre readout (TRestFibreReadout) and of the light processes: the primary scintillation
/// (TRestDetectorScintillationProcess, light yield of the medium) and the collection of the light by the
/// channels and its photoelectrons (TRestDetectorLightCollectionProcess, models localCapture and attenuation,
/// quantum efficiency, a hit per photoelectron with the scintillation time profile of its particle) and their
/// signals (TRestDetectorHitsToSignalProcess). The readout is a box of 100 mm with two families of fibres of
/// pitch 20 mm, built in memory; the attenuation model uses a minimal readout of the test (a bar read at
/// x = 0). The processes get the readout and the medium from a TRestManager created from a YAML string.

#include <TFile.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TRestDetectorHitsEvent.h>
#include <TRestDetectorHitsToSignalProcess.h>
#include <TRestDetectorLightCollectionProcess.h>
#include <TRestDetectorMedium.h>
#include <TRestDetectorReadoutManager.h>
#include <TRestDetectorScintillationProcess.h>
#include <TRestDetectorSignalEvent.h>
#include <TRestDetectorSignalToHitsProcess.h>
#include <TRestFibreReadout.h>
#include <TRestManager.h>
#include <TRestRun.h>
#include <TTree.h>
#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

using ROOT::Math::XYZVector;

namespace {

/// Analysis tree in memory (no file) given to the processes through a TRestRun, to read their observables
struct AnalysisTree {
    TRestRun run;
    std::unique_ptr<TTree> tree;

    AnalysisTree() : tree(std::make_unique<TTree>("AnalysisTree", "AnalysisTree")) {
        tree->SetDirectory(nullptr);
        run.fOutputAnalysisTree = tree.get();
    }
    ~AnalysisTree() { run.fOutputAnalysisTree = nullptr; }

    template <typename T>
    T Get(const std::string& name) const {
        return run.GetObservable<T>(name);
    }
};

// --- Fibre readout of the tests: box of 100 mm centred at (0, 0, 50) (x, y in [-50, 50], z in [0, 100]).
// Fibres along Y (transverse Z and X): z = 10 + 20 i, x = -40 + 20 j, physical IDs 0-24 (i * 5 + j).
// Fibres along X (transverse Y and Z): y = -40 + 20 i, z = 20 + 20 j, physical IDs 25-49 (25 + i * 5 + j).
// Without decoding file: DAQ ID = physical ID ---
YAML::Node ReadoutNode() {
    return YAML::Load(
        "class: TRestFibreReadout\n"
        "name: fibreTest\n"
        "readoutParameters:\n"
        "  size: [100mm, 100mm, 100mm]\n"
        "  position: [0, 0, 50mm]\n"
        "  fibres:\n"
        "    - {axis: Y, pitch: 20mm, offset: [10mm, 10mm], volume: fibreYVolume}\n"
        "    - {axis: X, pitch: 20mm, offset: [10mm, 20mm], volume: fibreXVolume}\n");
}

std::unique_ptr<TRestFibreReadout> BuildReadout() {
    auto readout = std::make_unique<TRestFibreReadout>("fibreTest", ReadoutNode());
    readout->BuildGeometry();
    return readout;
}

/// A bar along X (x in [0, 100]) read at x = 0: one channel (0), distance to the readout x
class BarReadout : public TRestDetectorReadout {
   public:
    BarReadout() : TRestDetectorReadout() { fName = "bar"; }
    void BuildGeometry() override {}
    int GetChannelFromPosition(double x, double, double) const override {
        return x >= 0 && x <= 100 ? 0 : -1;
    }
    using TRestDetectorReadout::GetDistanceToReadout;
    double GetDistanceToReadout(int daqID, const XYZVector& position) const override {
        return daqID == 0 ? position.X() : std::numeric_limits<double>::quiet_NaN();
    }
};

/// A TRestManager with the fibre readout (fibreTest) and the bar (bar) in its TRestDetectorReadoutManager
/// and a scintillator medium with only `light` (no charge), plus the sections of `extra` (other media)
struct Detector {
    std::unique_ptr<TRestManager> manager;
    TRestDetectorReadout* readout = nullptr;
};

Detector MakeDetector(const std::string& light =
                          "{lightYield: 100, attenuationLength: 100 mm, "
                          "effectiveLightSpeed: 100 mm/us}",
                      const std::string& extra = "") {
    const YAML::Node node = YAML::Load(
        "class: TRestManager\n"
        "readouts:\n"
        "  class: TRestDetectorReadoutManager\n"
        "  readout: fibreTest\n"
        "scintillator:\n"
        "  class: TRestDetectorMedium\n"
        "  name: LS\n"
        "  state: liquid\n"
        "  light: " +
        light + "\n" + extra);
    Detector detector;
    detector.manager = std::make_unique<TRestManager>("manager", node);
    auto* readouts = dynamic_cast<TRestDetectorReadoutManager*>(
        detector.manager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!readouts) throw std::runtime_error("no TRestDetectorReadoutManager in the manager");
    auto readout = BuildReadout();
    detector.readout = readout.get();
    readouts->fReadoutMap["fibreTest"] = readout.release();  // owned by the readout manager
    readouts->fReadoutMap["bar"] = new BarReadout();
    return detector;
}

/// A process initialized with the manager of the detector and the analysis tree
template <typename Process>
std::unique_ptr<Process> MakeProcess(const Detector& detector, AnalysisTree& analysis,
                                     const std::string& name, const std::string& yaml) {
    auto process = std::make_unique<Process>(name, YAML::Load(yaml));
    process->SetManager(detector.manager.get());
    process->SetRunInfo(&analysis.run);
    process->Initialize();
    return process;
}

}  // namespace

// ============================================================================
// Fibre readout
// ============================================================================

TEST(FibreReadout, GridChannelsAndPositions) {
    const auto readout = BuildReadout();
    EXPECT_EQ(readout->GetNumberOfFibres(), 50);

    // Closest fibre: the one along Y at z = 30, x = 0 (physical ID 1 * 5 + 2)
    EXPECT_EQ(readout->GetChannelFromPosition(0.5, 33, 31), 7);
    // On a fibre along X (y = -20, z = 60: 25 + 1 * 5 + 2)
    EXPECT_EQ(readout->GetChannelFromPosition(10, -20, 60), 32);
    // Outside the box
    EXPECT_EQ(readout->GetChannelFromPosition(0, 0, -1), -1);
    EXPECT_EQ(readout->GetChannelFromPosition(0, 60, 50), -1);

    // A fibre measures its transverse coordinates, NaN along it
    const auto yFibre = readout->GetPositionFromChannel(7);
    EXPECT_NEAR(yFibre.X(), 0, 1e-9);
    EXPECT_TRUE(std::isnan(yFibre.Y()));
    EXPECT_NEAR(yFibre.Z(), 30, 1e-9);
    const auto xFibre = readout->GetPositionFromChannel(32);
    EXPECT_TRUE(std::isnan(xFibre.X()));
    EXPECT_NEAR(xFibre.Y(), -20, 1e-9);
    EXPECT_NEAR(xFibre.Z(), 60, 1e-9);
    EXPECT_TRUE(std::isnan(readout->GetPositionFromChannel(50).X()));  // not a fibre

    // Neighbours in the same row, none across rows or families
    EXPECT_EQ(readout->GetAdjacentChannels(7), std::make_pair(6, 8));
    EXPECT_EQ(readout->GetAdjacentChannels(5), std::make_pair(-1, 6));
    EXPECT_EQ(readout->GetAdjacentChannels(24), std::make_pair(23, -1));

    const auto [low, high] = readout->GetChannelsBoundingBox();
    EXPECT_NEAR(low.X(), -50, 1e-9);
    EXPECT_NEAR(high.Z(), 100, 1e-9);
}

TEST(FibreReadout, ChannelsInRange) {
    const auto readout = BuildReadout();

    // On the fibre along Y at z = 30, x = 0 (ID 7): it and the fibres along X at y = 0, z = 20 and z = 40
    // (IDs 35 and 36, 10 mm away) are within 15 mm; the others are 20 mm or more away
    const auto channels = readout->GetChannelsInRange({0, 0, 30}, 15);
    ASSERT_EQ(channels.size(), 3u);
    std::set<int> ids;
    for (const auto& channel : channels) {
        ids.insert(channel.daqID);
        EXPECT_NEAR(channel.distance, channel.daqID == 7 ? 0 : 10, 1e-9);
        // The position is on the fibre: it gives back the same channel
        const auto& p = channel.position;
        EXPECT_EQ(readout->GetChannelFromPosition(p.X(), p.Y(), p.Z()), channel.daqID);
    }
    EXPECT_EQ(ids, std::set<int>({7, 35, 36}));
    EXPECT_EQ(readout->GetChannelsInRange({0, 0, 30}, 5).size(), 1u);

    // Along X outside the box: the fibres along X do not reach it, those along Y are too far
    EXPECT_TRUE(readout->GetChannelsInRange({60, 0, 30}, 15).empty());

    // Read end ("+" by default): the fibre 7 along Y at y = +50, 50 mm along it from (0, 0, 30)
    const auto end = readout->GetReadoutPosition(7);
    EXPECT_NEAR(end.X(), 0, 1e-9);
    EXPECT_NEAR(end.Y(), 50, 1e-9);
    EXPECT_NEAR(end.Z(), 30, 1e-9);
    EXPECT_NEAR(readout->GetDistanceToReadout(7, {0, 0, 30}), 50, 1e-9);
    EXPECT_TRUE(std::isnan(readout->GetDistanceToReadout(XYZVector(0, 0, 30))));  // no drift
    EXPECT_TRUE(std::isnan(readout->GetDistanceToReadout(99, {0, 0, 30})));       // not a channel
}

TEST(FibreReadout, ReadAtBothEnds) {
    // The fibres along X read at both ends: "-" end IDs 25-49, "+" end IDs 50-74
    auto node = ReadoutNode();
    node["readoutParameters"]["fibres"][1]["readoutEnd"] = "both";
    TRestFibreReadout readout("fibreTest", node);
    readout.BuildGeometry();
    EXPECT_EQ(readout.GetNumberOfFibres(), 50);
    EXPECT_EQ(readout.GetNumberOfChannels(), 75);

    // The fibre along X at y = 0, z = 20: channel 35 at x = -50, channel 60 at x = +50
    EXPECT_NEAR(readout.GetReadoutPosition(35).X(), -50, 1e-9);
    EXPECT_NEAR(readout.GetReadoutPosition(60).X(), 50, 1e-9);
    EXPECT_NEAR(readout.GetDistanceToReadout(35, {10, 0, 20}), 60, 1e-9);
    EXPECT_NEAR(readout.GetDistanceToReadout(60, {10, 0, 20}), 40, 1e-9);
    // The channel of the end closer to the position; both ends measure the same coordinates
    EXPECT_EQ(readout.GetChannelFromPosition(10, 0, 20), 60);
    EXPECT_EQ(readout.GetChannelFromPosition(-10, 0, 20), 35);
    EXPECT_NEAR(readout.GetPositionFromChannel(60).Z(), 20, 1e-9);
    EXPECT_TRUE(std::isnan(readout.GetPositionFromChannel(60).X()));
    EXPECT_EQ(readout.GetAdjacentChannels(60), std::make_pair(-1, 61));

    // Around the fibre 7: it, and both ends of the fibres along X at z = 20 and z = 40, at their read ends
    const auto channels = readout.GetChannelsInRange({0, 0, 30}, 15);
    std::set<int> ids;
    for (const auto& channel : channels) {
        ids.insert(channel.daqID);
        const auto& p = channel.position;
        EXPECT_EQ(readout.GetChannelFromPosition(p.X(), p.Y(), p.Z()), channel.daqID);
    }
    EXPECT_EQ(ids, std::set<int>({7, 35, 36, 60, 61}));

    node["readoutParameters"]["fibres"][1]["readoutEnd"] = "middle";
    EXPECT_ANY_THROW(TRestFibreReadout("fibreTest", node));
}

TEST(FibreReadout, WrongParameters) {
    YAML::Node node = ReadoutNode();
    node["readoutParameters"]["fibres"][0]["axis"] = "W";
    EXPECT_ANY_THROW(TRestFibreReadout("fibreTest", node));

    node = ReadoutNode();
    node["readoutParameters"]["fibres"][1]["pitch"] = "-1 mm";
    EXPECT_ANY_THROW(TRestFibreReadout("fibreTest", node));
}

TEST(FibreReadout, GDMLPlacements) {
    const auto readout = BuildReadout();
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/fibres.xml";
    ASSERT_TRUE(readout->WriteGDMLPlacements(fileName));

    std::ifstream file(fileName);
    std::stringstream buffer;
    buffer << file.rdbuf();
    const std::string text = buffer.str();
    size_t count = 0;
    for (size_t at = text.find("<physvol"); at != std::string::npos; at = text.find("<physvol", at + 1)) {
        count++;
    }
    EXPECT_EQ(count, 50u);
    // The first fibre along Y (z = 10, x = -40), relative to the centre (0, 0, 50), rotated to Y
    EXPECT_NE(text.find("<physvol name=\"fibre_0\" copynumber=\"0\">"), std::string::npos);
    EXPECT_NE(text.find("name=\"fibre_0_position\" unit=\"mm\" x=\"-40\" y=\"0\" z=\"-40\""),
              std::string::npos);
    EXPECT_NE(text.find("name=\"fibre_0_rotation\" unit=\"deg\" x=\"90\" y=\"0\" z=\"0\""),
              std::string::npos);
    EXPECT_NE(text.find("<volumeref ref=\"fibreXVolume\"/>"), std::string::npos);
}

TEST(FibreReadout, ExportAndImportWithManager) {
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/fibre_test_readout.root";
    {
        const auto readout = BuildReadout();
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        ASSERT_TRUE(file && !file->IsZombie());
        ASSERT_TRUE(readout->Export(file.get(), "fibreTest", "test_decoding"));
        file->Close();
    }

    const YAML::Node node = YAML::Load(
        "readout: fibreTest\n"
        "readouts:\n"
        "  - {instanceName: fibreTest, decodingName: test_decoding, inputFileName: \"" +
        fileName + "\"}\n");
    TRestDetectorReadoutManager manager("readouts", node);
    manager.LoadReadout();
    std::string name;
    TRestDetectorReadout* readout = manager.SelectReadout(name, "test");
    ASSERT_NE(readout, nullptr);
    EXPECT_EQ(readout->GetClassName(), "TRestFibreReadout");
    EXPECT_EQ(readout->GetChannelFromPosition(0.5, 33, 31), 7);
    EXPECT_EQ(readout->GetChannelsInRange({0, 0, 30}, 15).size(), 3u);
}

// ============================================================================
// Light: scintillation, collection and photoelectrons
// ============================================================================

TEST(DetectorMedium, LightYield) {
    TRestDetectorMedium medium("LS", YAML::Load("{state: liquid, light: {lightYield: 12}}"));
    EXPECT_DOUBLE_EQ(medium.GetLightResponse("test").GetLightYield(), 12);
    EXPECT_EQ(medium.GetWLSResponse(), nullptr);
    EXPECT_THROW(medium.GetWLSResponse("test"), std::runtime_error);

    // Wavelength shifting fibres: attenuation, speed and re-emission along them
    TRestDetectorMedium fibre(
        "Y11", YAML::Load("{state: solid, wls: {attenuationLength: 3.5 m, lightSpeed: 170 mm/ns, "
                          "decayTime: 8 ns}}"));
    const auto& wls = fibre.GetWLSResponse("test");
    EXPECT_DOUBLE_EQ(wls.GetAttenuationLength(), 3500);
    EXPECT_DOUBLE_EQ(wls.GetLightSpeed(), 170000);
    EXPECT_DOUBLE_EQ(wls.GetDecayTime(), 0.008);
    TRandom3 random(2);
    double mean = 0;
    for (int i = 0; i < 10000; i++) mean += wls.GetReemissionTime(random) / 10000;
    EXPECT_NEAR(mean, 0.008, 5 * 0.008 / 100);
    EXPECT_ANY_THROW(TRestDetectorMedium("Y11", YAML::Load("{state: solid, wls: {decayTime: -1 ns}}")));
    EXPECT_ANY_THROW(TRestDetectorMedium("LS", YAML::Load("{state: liquid, light: {lightYield: -1}}")));
}

TEST(DetectorScintillationProcess, PhotonsPerHit) {
    auto detector = MakeDetector();
    AnalysisTree analysis;
    auto process = MakeProcess<TRestDetectorScintillationProcess>(detector, analysis, "scintillation",
                                                                  "{seed: 7, observables: [NPhotons]}");
    EXPECT_EQ(process->fMediumName, "LS");

    // 10 keV, 100 photons per keV: 1000 photons on average, at the same position and time
    TRestDetectorHitsEvent hits;
    hits.SetID(3);
    hits.AddHit(1, 2, 30, 10, 2);
    hits.AddHit(5, 5, 5, 0, 2);  // no energy: dropped
    TRestDetectorHitsEvent photons;
    ASSERT_TRUE(process->ProcessEvent(hits, photons));
    EXPECT_EQ(photons.GetID(), 3);
    ASSERT_EQ(photons.GetNumberOfHits(), 1u);
    EXPECT_NEAR(photons.GetEnergy(0), 1000, 5 * std::sqrt(1000.));
    EXPECT_EQ(photons.GetEnergy(0), std::round(photons.GetEnergy(0)));
    EXPECT_NEAR(photons.GetZ(0), 30, 1e-9);
    EXPECT_DOUBLE_EQ(photons.GetTime(0), 2);
    EXPECT_DOUBLE_EQ(analysis.Get<double>("scintillation_NPhotons"), photons.GetEnergy(0));

    // Same seed and event: the same photons
    TRestDetectorHitsEvent again;
    ASSERT_TRUE(process->ProcessEvent(hits, again));
    EXPECT_DOUBLE_EQ(again.GetEnergy(0), photons.GetEnergy(0));

    // Without lightYield
    auto noYield = MakeDetector("{attenuationLength: 1 m}");
    TRestDetectorScintillationProcess failing("scintillation", YAML::Load("{}"));
    failing.SetManager(noYield.manager.get());
    EXPECT_THROW(failing.Initialize(), std::runtime_error);
}

TEST(DetectorLightCollectionProcess, LocalCaptureByTheFibres) {
    auto detector = MakeDetector();
    AnalysisTree analysis;
    auto process = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection",
        "{model: localCapture, captureProbability: \"0.2+0*x\", maxDistance: 15mm, seed: 5, "
        "observables: [NPhotoelectrons, NHitsWithoutChannel]}");

    // 1000 photons on the fibre 7: the fibres 7, 35 and 36 capture 200 each on average
    TRestDetectorHitsEvent photons;
    photons.AddHit(0, 0, 30, 1000, 2);
    photons.AddHit(60, 0, 30, 1000, 2);  // no fibre in range
    TRestDetectorHitsEvent collected;
    ASSERT_TRUE(process->ProcessEvent(photons, collected));
    EXPECT_EQ(analysis.Get<int>("lightCollection_NHitsWithoutChannel"), 1);
    // A hit per photoelectron, of charge 1
    const double total = analysis.Get<double>("lightCollection_NPhotoelectrons");
    ASSERT_EQ(collected.GetNumberOfHits(), static_cast<size_t>(total));
    EXPECT_NEAR(collected.GetTotalEnergy(), total, 1e-6);
    EXPECT_NEAR(total, 600, 5 * std::sqrt(600.));
    std::set<int> channels;
    for (size_t n = 0; n < collected.GetNumberOfHits(); n++) {
        EXPECT_DOUBLE_EQ(collected.GetTime(n), 2);
        const auto p = collected.GetPosition(n);
        channels.insert(detector.readout->GetChannelFromPosition(p.X(), p.Y(), p.Z()));
    }
    EXPECT_EQ(channels, std::set<int>({7, 35, 36}));

    // Probabilities adding up to more than 1 around a position
    auto tooLarge = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection", "{captureProbability: \"0.5+0*x\", maxDistance: 15mm}");
    TRestDetectorHitsEvent output;
    EXPECT_THROW(tooLarge->ProcessEvent(photons, output), std::runtime_error);

    // WLS fibres: the photons captured by the fibre 7 travel 50 mm along it to its read end (y = +50),
    // exp(-50 / 100) of them arrive, 0.5 us later
    auto withFibres = MakeDetector("{lightYield: 100}",
                                   "fibres:\n"
                                   "  class: TRestDetectorMedium\n"
                                   "  name: Y11\n"
                                   "  state: solid\n"
                                   "  wls: {attenuationLength: 100 mm, lightSpeed: 100 mm/us}\n");
    auto guided = MakeProcess<TRestDetectorLightCollectionProcess>(
        withFibres, analysis, "lightCollection",
        "{captureProbability: \"1+0*x\", maxDistance: 5mm, wlsMedium: Y11, seed: 5}");
    TRestDetectorHitsEvent onFibre, guidedOutput;
    onFibre.AddHit(0, 0, 30, 10000, 2);
    ASSERT_TRUE(guided->ProcessEvent(onFibre, guidedOutput));
    const double expected = 10000 * std::exp(-0.5);
    EXPECT_NEAR(guidedOutput.GetNumberOfHits(), expected, 5 * std::sqrt(expected));
    EXPECT_NEAR(guidedOutput.GetTime(0), 2.5, 1e-6);
    EXPECT_NEAR(guidedOutput.GetY(0), 50, 1e-9);  // at the read end
    // wlsMedium only for localCapture, and it needs a wls response
    EXPECT_ANY_THROW(MakeProcess<TRestDetectorLightCollectionProcess>(
        withFibres, analysis, "lightCollection", "{model: attenuation, medium: LS, wlsMedium: Y11}"));
    EXPECT_ANY_THROW(MakeProcess<TRestDetectorLightCollectionProcess>(
        withFibres, analysis, "lightCollection",
        "{captureProbability: \"1+0*x\", maxDistance: 5mm, wlsMedium: LS}"));

    // Wrong configurations
    EXPECT_ANY_THROW(TRestDetectorLightCollectionProcess("c", YAML::Load("{model: unknown}")));
    EXPECT_ANY_THROW(TRestDetectorLightCollectionProcess("c", YAML::Load("{captureProbability: \"0.1\"}")));
    EXPECT_ANY_THROW(TRestDetectorLightCollectionProcess(
        "c", YAML::Load("{captureProbability: \"x*(\", maxDistance: 1}")));
}

TEST(DetectorLightCollectionProcess, CaptureProbabilityFromGraph) {
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/light_capture.root";
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        ASSERT_TRUE(file && !file->IsZombie());
        TGraph graph(2);
        graph.SetPoint(0, 0, 0.5);
        graph.SetPoint(1, 20, 0);
        graph.Write("capture");
        file->Close();
    }
    TRestDetectorLightCollectionProcess process(
        "lightCollection",
        YAML::Load("{captureProbability: \"" + fileName + ":capture\", maxDistance: 15 mm}"));
    EXPECT_NEAR(process.GetCaptureProbability(0), 0.5, 1e-9);
    EXPECT_NEAR(process.GetCaptureProbability(10), 0.25, 1e-9);
    EXPECT_DOUBLE_EQ(process.GetCaptureProbability(16), 0);  // beyond maxDistance

    EXPECT_ANY_THROW(TRestDetectorLightCollectionProcess(
        "c", YAML::Load("{captureProbability: \"" + fileName + ":missing\", maxDistance: 15 mm}")));
}

TEST(DetectorLightCollectionProcess, AttenuationToTheReadout) {
    auto detector = MakeDetector();
    AnalysisTree analysis;
    auto process = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection",
        "{model: attenuation, readout: bar, seed: 3, observables: [NPhotoelectrons, NHitsWithoutChannel]}");

    // 10000 photons at 50 mm from the readout: exp(-50 / 100) of them, 0.5 us later
    TRestDetectorHitsEvent photons;
    photons.AddHit(50, 0, 0, 10000, 1);
    photons.AddHit(200, 0, 0, 10000, 1);  // no channel
    TRestDetectorHitsEvent collected;
    ASSERT_TRUE(process->ProcessEvent(photons, collected));
    const double expected = 10000 * std::exp(-0.5);
    EXPECT_NEAR(collected.GetNumberOfHits(), expected, 5 * std::sqrt(expected));
    EXPECT_NEAR(collected.GetTime(0), 1.5, 1e-6);
    EXPECT_NEAR(collected.GetX(0), 50, 1e-9);
    EXPECT_EQ(analysis.Get<int>("lightCollection_NHitsWithoutChannel"), 1);

    // A fibre readout: the distance along the fibre of the hit to its read end (50 mm to y = +50)
    auto fibres = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection", "{model: attenuation, readout: fibreTest, seed: 3}");
    TRestDetectorHitsEvent onFibre;
    onFibre.AddHit(0, 0, 30, 10000, 0);
    TRestDetectorHitsEvent output;
    ASSERT_TRUE(fibres->ProcessEvent(onFibre, output));
    EXPECT_NEAR(output.GetNumberOfHits(), expected, 5 * std::sqrt(expected));
    EXPECT_NEAR(output.GetTime(0), 0.5, 1e-6);
}

TEST(DetectorLightCollectionProcess, QuantumEfficiency) {
    auto detector = MakeDetector();
    AnalysisTree analysis;
    // On the fibre 7 only (maxDistance 5 mm), all the photons captured: Binomial(10000, 0.4) photoelectrons
    auto process = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection",
        "{captureProbability: \"1+0*x\", maxDistance: 5mm, quantumEfficiency: 0.4, seed: 9, "
        "observables: [NPhotoelectrons]}");

    TRestDetectorHitsEvent photons;
    photons.AddHit(0, 0, 30, 10000, 1);
    TRestDetectorHitsEvent photoelectrons;
    ASSERT_TRUE(process->ProcessEvent(photons, photoelectrons));
    EXPECT_NEAR(photoelectrons.GetNumberOfHits(), 4000, 5 * std::sqrt(10000 * 0.4 * 0.6));
    EXPECT_DOUBLE_EQ(analysis.Get<double>("lightCollection_NPhotoelectrons"),
                     photoelectrons.GetNumberOfHits());

    auto blind = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection",
        "{captureProbability: \"1+0*x\", maxDistance: 5mm, quantumEfficiency: 0}");
    TRestDetectorHitsEvent none;
    EXPECT_FALSE(blind->ProcessEvent(photons, none));

    EXPECT_ANY_THROW(TRestDetectorLightCollectionProcess(
        "c", YAML::Load("{captureProbability: \"1\", maxDistance: 5mm, quantumEfficiency: 1.5}")));
}

TEST(DetectorLightProcesses, FromEnergyToHitsOfTheFibres) {
    auto detector = MakeDetector();
    AnalysisTree analysis;
    auto scintillation =
        MakeProcess<TRestDetectorScintillationProcess>(detector, analysis, "scintillation", "{seed: 3}");
    auto collection = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection",
        "{captureProbability: \"0.2+0*x\", maxDistance: 15mm, quantumEfficiency: 0.5, seed: 3}");
    // Generic processes: no medium and no drift
    auto hitsToSignal =
        MakeProcess<TRestDetectorHitsToSignalProcess>(detector, analysis, "hitsToSignal", "{}");
    auto signalToHits = MakeProcess<TRestDetectorSignalToHitsProcess>(detector, analysis, "signalToHits",
                                                                      "{method: integral}");

    TRestDetectorHitsEvent hits;
    hits.AddHit(0, 0, 30, 10, 4);
    TRestDetectorHitsEvent photons, photoelectrons;
    ASSERT_TRUE(scintillation->ProcessEvent(hits, photons));
    ASSERT_TRUE(collection->ProcessEvent(photons, photoelectrons));
    TRestDetectorSignalEvent signals;
    ASSERT_TRUE(hitsToSignal->ProcessEvent(photoelectrons, signals));
    ASSERT_EQ(signals.GetNumberOfSignals(), 3);
    for (int id : {7, 35, 36}) {
        ASSERT_TRUE(signals.HasSignalID(id));
        EXPECT_NEAR(signals.GetSignalByID(id).GetPoint(0).first, 4, 1e-9);
    }

    // Back to hits: an XZ hit for the fibre along Y, YZ hits for those along X
    TRestDetectorHitsEvent reconstructed;
    ASSERT_TRUE(signalToHits->ProcessEvent(signals, reconstructed));
    ASSERT_EQ(reconstructed.GetNumberOfHits(), 3u);
    int xz = 0, yz = 0;
    for (size_t n = 0; n < reconstructed.GetNumberOfHits(); n++) {
        if (reconstructed.GetType(n) == TRestHitsData::REST_HitType::XZ) xz++;
        if (reconstructed.GetType(n) == TRestHitsData::REST_HitType::YZ) yz++;
    }
    EXPECT_EQ(xz, 1);
    EXPECT_EQ(yz, 2);

    // The media are found by response: with a scintillator (light) and WLS fibres (wls) no name is needed
    auto withFibres = MakeDetector("{lightYield: 100}",
                                   "fibres:\n"
                                   "  class: TRestDetectorMedium\n"
                                   "  name: Y11\n"
                                   "  state: solid\n"
                                   "  wls: {attenuationLength: 100 mm, lightSpeed: 100 mm/us}\n");
    const auto& manager = *withFibres.manager;
    EXPECT_EQ(TRestDetectorMedium::Find(manager, "", "light", "test")->GetName(), "LS");
    EXPECT_EQ(TRestDetectorMedium::Find(manager, "", "wls", "test")->GetName(), "Y11");
    EXPECT_EQ(TRestDetectorMedium::Find(manager, "", "charge", "test"), nullptr);
    EXPECT_TRUE(TRestDetectorMedium::Find(manager, "Y11", "wls", "test")->HasResponse("wls"));
    EXPECT_ANY_THROW(TRestDetectorMedium::Find(manager, "Y11", "light", "test"));    // without light
    EXPECT_ANY_THROW(TRestDetectorMedium::Find(manager, "Other", "light", "test"));  // not found
    EXPECT_ANY_THROW(TRestDetectorMedium::Find(manager, "", "optics", "test"));      // unknown response
    EXPECT_NO_THROW(
        MakeProcess<TRestDetectorScintillationProcess>(withFibres, analysis, "scintillation", "{}"));
    EXPECT_NO_THROW(MakeProcess<TRestDetectorLightCollectionProcess>(
        withFibres, analysis, "lightCollection",
        "{captureProbability: \"0.2+0*x\", maxDistance: 15mm, wlsMedium: Y11}"));
    // Two media with light: the name is needed
    auto twoScintillators = MakeDetector("{lightYield: 100}",
                                         "other:\n"
                                         "  class: TRestDetectorMedium\n"
                                         "  name: BC408\n"
                                         "  state: solid\n"
                                         "  light: {lightYield: 10}\n");
    EXPECT_ANY_THROW(
        MakeProcess<TRestDetectorScintillationProcess>(twoScintillators, analysis, "scintillation", "{}"));
    EXPECT_NO_THROW(MakeProcess<TRestDetectorScintillationProcess>(twoScintillators, analysis,
                                                                   "scintillation", "{medium: BC408}"));
}

TEST(DetectorLightCollectionProcess, PhotoelectronTimesPerParticle) {
    // Default profile: 10 ns; alpha: 100 ns (times in us)
    auto detector = MakeDetector(
        "{lightYield: 100, scintillation: [{components: [{decayTime: 10 ns}]},"
        "{particles: [alpha], components: [{decayTime: 100 ns}]}]}");
    AnalysisTree analysis;
    // Every photon on a fibre within 5 mm reaches its SiPM and gives a photoelectron
    const std::string capture = "captureProbability: \"1+0*x\", maxDistance: 5mm";
    auto process = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection",
        "{" + capture + ", seed: 4, observables: [NPhotoelectrons, NHitsWithoutChannel]}");
    EXPECT_EQ(process->fMediumName, "LS");

    // 2000 photons of an electron on the fibre 7 (along Y, x = 0, z = 30) and of an alpha on the fibre 35
    // (along X, y = 0, z = 20), at t = 1 us
    TRestDetectorHitsEvent photons;
    photons.SetID(2);
    photons.AddHit(0, 0, 30, 2000, 1, TRestHitsData::REST_HitType::XYZ, 11);
    photons.AddHit(10, 0, 20, 2000, 1, TRestHitsData::REST_HitType::XYZ, 1000020040);
    photons.AddHit(0, 60, 50, 10, 1);  // outside the readout
    TRestDetectorHitsEvent photoelectrons;
    ASSERT_TRUE(process->ProcessEvent(photons, photoelectrons));
    EXPECT_EQ(photoelectrons.GetID(), 2);
    EXPECT_DOUBLE_EQ(analysis.Get<double>("lightCollection_NPhotoelectrons"), 4000);
    EXPECT_EQ(analysis.Get<int>("lightCollection_NHitsWithoutChannel"), 1);
    ASSERT_EQ(photoelectrons.GetNumberOfHits(), 4000u);

    // Mean time of the photoelectrons of each channel: 1 us + the decay time of the profile of its particle
    auto meanTimes = [&detector](const TRestDetectorHitsEvent& event) {
        std::map<int, std::pair<double, double>> sums;  // channel: (time * charge, charge)
        for (size_t n = 0; n < event.GetNumberOfHits(); n++) {
            const auto p = event.GetPosition(n);
            auto& [weighted, charge] = sums[detector.readout->GetChannelFromPosition(p.X(), p.Y(), p.Z())];
            EXPECT_GE(event.GetTime(n), 1);
            weighted += event.GetTime(n) * event.GetEnergy(n);
            charge += event.GetEnergy(n);
        }
        std::map<int, double> means;
        for (const auto& [channel, sum] : sums) means[channel] = sum.first / sum.second;
        return means;
    };
    auto means = meanTimes(photoelectrons);
    ASSERT_EQ(means.size(), 2u);
    EXPECT_NEAR(means[7], 1.010, 5 * 0.010 / std::sqrt(2000.));
    EXPECT_NEAR(means[35], 1.100, 5 * 0.100 / std::sqrt(2000.));
    for (size_t n = 0; n < photoelectrons.GetNumberOfHits(); n++) {
        EXPECT_DOUBLE_EQ(photoelectrons.GetEnergy(n), 1);
    }

    // Same seed and event: the same photoelectrons
    TRestDetectorHitsEvent again;
    ASSERT_TRUE(process->ProcessEvent(photons, again));
    ASSERT_EQ(again.GetNumberOfHits(), photoelectrons.GetNumberOfHits());
    EXPECT_DOUBLE_EQ(again.GetTime(10), photoelectrons.GetTime(10));

    // The signals (generic TRestDetectorHitsToSignalProcess): 2000 photoelectrons in each channel
    auto hitsToSignal = MakeProcess<TRestDetectorHitsToSignalProcess>(detector, analysis, "hitsToSignal",
                                                                      "{sampling: 50 ns}");
    TRestDetectorSignalEvent signals;
    ASSERT_TRUE(hitsToSignal->ProcessEvent(photoelectrons, signals));
    ASSERT_EQ(signals.GetNumberOfSignals(), 2);
    for (int id : {7, 35}) {
        ASSERT_TRUE(signals.HasSignalID(id));
        double charge = 0;
        for (int i = 0; i < signals.GetSignalByID(id).GetNPoints(); i++) {
            charge += signals.GetSignalByID(id).GetPoint(i).second;
        }
        EXPECT_NEAR(charge, 2000, 1e-6);
    }

    // WLS fibres: each photoelectron also delayed by the re-emission (8 ns more on average)
    auto withFibres = MakeDetector("{lightYield: 100, scintillation: [{components: [{decayTime: 10 ns}]}]}",
                                   "fibres:\n"
                                   "  class: TRestDetectorMedium\n"
                                   "  name: Y11\n"
                                   "  state: solid\n"
                                   "  wls: {decayTime: 8 ns}\n");
    auto shifted = MakeProcess<TRestDetectorLightCollectionProcess>(
        withFibres, analysis, "lightCollection", "{" + capture + ", wlsMedium: Y11, seed: 4}");
    TRestDetectorHitsEvent delayed;
    ASSERT_TRUE(shifted->ProcessEvent(photons, delayed));
    EXPECT_NEAR(meanTimes(delayed)[7], 1.018, 5 * 0.013 / std::sqrt(2000.));
    // A wlsMedium without wls
    TRestDetectorLightCollectionProcess noWLS("lightCollection",
                                              YAML::Load("{" + capture + ", wlsMedium: LS}"));
    noWLS.SetManager(withFibres.manager.get());
    EXPECT_THROW(noWLS.Initialize(), std::runtime_error);

    // Gain spread: the charge of each photoelectron Gauss(1, 0.3), 2000 of them on the fibre 7
    auto spread = MakeProcess<TRestDetectorLightCollectionProcess>(
        detector, analysis, "lightCollection", "{" + capture + ", seed: 4, gainSpread: 0.3}");
    TRestDetectorHitsEvent fluctuating;
    ASSERT_TRUE(spread->ProcessEvent(photons, fluctuating));
    double charge7 = 0;
    for (size_t n = 0; n < fluctuating.GetNumberOfHits(); n++) {
        const auto p = fluctuating.GetPosition(n);
        if (detector.readout->GetChannelFromPosition(p.X(), p.Y(), p.Z()) == 7) {
            charge7 += fluctuating.GetEnergy(n);
        }
    }
    EXPECT_NEAR(charge7, 2000, 5 * 0.3 * std::sqrt(2000.));
    EXPECT_NE(charge7, 2000);
    EXPECT_DOUBLE_EQ(analysis.Get<double>("lightCollection_NPhotoelectrons"), 4000);

    // Negative jitter or gain spread
    EXPECT_ANY_THROW(
        TRestDetectorLightCollectionProcess("c", YAML::Load("{" + capture + ", timeJitter: -1 ns}")));
    EXPECT_ANY_THROW(
        TRestDetectorLightCollectionProcess("c", YAML::Load("{" + capture + ", gainSpread: -1}")));
}
