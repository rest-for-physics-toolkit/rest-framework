/// Unit tests of the detector library (GoogleTest, BUILD_TESTING): the hits and signal events, a Micromegas
/// readout built in memory (with a decoding file written in the work directory) and its export/import, a
/// GDML readout of a small geometry, the medium with fixed transport parameters or a gain table, a drift
/// response built by hand, and the hits and signal processes fed with events built by hand.
/// The processes that need a readout or a medium get them from a TRestManager created from a YAML string,
/// with the readout added to its TRestDetectorReadoutManager; the observables are read from an in-memory
/// analysis tree given to the process through a TRestRun.

#include <TFile.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TRestDataSet.h>
#include <TRestDataSetGainMap.h>
#include <TRestDataSetTardis.h>
#include <TRestDetectorAmplificationResponse.h>
#include <TRestDetectorAvalancheProcess.h>
#include <TRestDetectorDriftResponse.h>
#include <TRestDetectorHitsAnalysisProcess.h>
#include <TRestDetectorHitsDriftProcess.h>
#include <TRestDetectorHitsEvent.h>
#include <TRestDetectorHitsSmearingProcess.h>
#include <TRestDetectorHitsTimeProjectionProcess.h>
#include <TRestDetectorHitsToSignalProcess.h>
#include <TRestDetectorHitsTransformationProcess.h>
#include <TRestDetectorMedium.h>
#include <TRestDetectorReadoutManager.h>
#include <TRestDetectorSignalEvent.h>
#include <TRestDetectorSignalPeaksFinderProcess.h>
#include <TRestDetectorSignalRecoveryProcess.h>
#include <TRestDetectorSignalShapingProcess.h>
#include <TRestDetectorSignalToHitsProcess.h>
#include <TRestDetectorTransportManager.h>
#include <TRestGDMLReadout.h>
#include <TRestManager.h>
#include <TRestMicromegasReadout.h>
#include <TRestRun.h>
#include <TTree.h>
#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using ROOT::Math::XYZVector;

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
using HitType = TRestHitsData::REST_HitType;

/// Analysis tree in memory (no file) given to the processes through a TRestRun, to read their observables
struct AnalysisTree {
    TRestRun run;
    std::unique_ptr<TTree> tree;

    AnalysisTree() : tree(std::make_unique<TTree>("AnalysisTree", "AnalysisTree")) {
        tree->SetDirectory(nullptr);
        run.fOutputAnalysisTree = tree.get();
    }
    ~AnalysisTree() { run.fOutputAnalysisTree = nullptr; }

    /// Current value of the observable processName_name
    template <typename T>
    T Get(const std::string& name) const {
        return run.GetObservable<T>(name);
    }
};

// --- Micromegas readout of the tests: 4 X strips and 4 Y strips of pitch 1 mm, plane at z = 0, normal +Z,
// drift volume of 50 mm. DAQ channel = 100 + physical ID (Y strips 0-3, X strips 4-7) ---
constexpr int kNChannels = 4;
constexpr double kPitch = 1;
constexpr double kThickness = 0.1;
constexpr double kGasThickness = 50;
// Offsets of the module (TRestMicromegasReadout::BuildGeometry)
constexpr double kOffsetX = -((kNChannels + 1) * kPitch - 0.5 * kPitch) / 2;
constexpr double kOffsetY = -((kNChannels + 1) * kPitch - 0.75 * kPitch) / 2;

int DaqID(int physicalID) { return 100 + physicalID; }
/// DAQ channel of the X strip n (physical ID nChannels + n) and of the Y strip n (physical ID n)
int XStrip(int n) { return DaqID(kNChannels + n); }
int YStrip(int n) { return DaqID(n); }
/// X of the X strip n and Y of the Y strip n (GetPositionFromChannel)
double XStripPosition(int n) { return n * kPitch - kPitch / 4 + kOffsetX; }
double YStripPosition(int n) { return (1 + n) * kPitch + kOffsetY; }
/// Centre of the pixel nPix of the X strip n and of the Y strip n, on the plane (z = 0)
XYZVector XStripPixel(int n, int nPix) { return {XStripPosition(n), (0.5 + nPix) * kPitch + kOffsetY, 0}; }
XYZVector YStripPixel(int n, int nPix) {
    return {kPitch / 4 + nPix * kPitch + kOffsetX, YStripPosition(n), 0};
}

/// Decoding file of the readout in the work directory, one per test (ctest runs them in parallel)
std::string WriteDecoding() {
    const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string path =
        std::string(TESTS_WORK_DIR) + "/readout_" + test->test_suite_name() + "_" + test->name() + ".dec";
    std::ofstream out(path);
    out << "# daqChannel physicalID\n";
    for (int physicalID = 0; physicalID < 2 * kNChannels; physicalID++) {
        out << DaqID(physicalID) << "\t" << physicalID << "\n";
    }
    return path;
}

YAML::Node ReadoutNode(const std::string& normal = "[0, 0, 1]") {
    return YAML::Load(
        "class: TRestMicromegasReadout\n"
        "name: mmTest\n"
        "decodingFile: \"" +
        WriteDecoding() +
        "\"\n"
        "readoutParameters:\n"
        "  nChannels: " +
        std::to_string(kNChannels) +
        "\n"
        "  pitch: 1 mm\n"
        "  thickness: 0.1 mm\n"
        "  gasThickness: 50 mm\n"
        "  normal: " +
        normal + "\n");
}

std::unique_ptr<TRestMicromegasReadout> BuildReadout(const std::string& normal = "[0, 0, 1]") {
    auto readout = std::make_unique<TRestMicromegasReadout>("mmTest", ReadoutNode(normal));
    readout->BuildGeometry();
    return readout;
}

/// Charge response with fixed transport parameters (no medium file)
const std::string kChargeNoDiffusion = "{wValue: 26.2 eV, fanoFactor: 0.2, driftVelocity: 10}";
const std::string kChargeDiffusion =
    "{wValue: 26.2 eV, fanoFactor: 0.2, driftVelocity: 10, longDiff: 0.03, transDiff: 0.04}";

/// A TRestManager with a TRestDetectorReadoutManager holding the readout of the tests and a medium
struct Detector {
    std::unique_ptr<TRestManager> manager;
    TRestDetectorReadoutManager* readouts = nullptr;
    TRestDetectorReadout* readout = nullptr;
};

Detector MakeDetector(const std::string& charge = kChargeNoDiffusion) {
    const YAML::Node node = YAML::Load(
        "class: TRestManager\n"
        "readouts:\n"
        "  class: TRestDetectorReadoutManager\n"
        "  readout: mmTest\n"
        "gas:\n"
        "  class: TRestDetectorMedium\n"
        "  name: ArIso\n"
        "  charge: " +
        charge + "\n");
    Detector detector;
    detector.manager = std::make_unique<TRestManager>("manager", node);
    detector.readouts = dynamic_cast<TRestDetectorReadoutManager*>(
        detector.manager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!detector.readouts) throw std::runtime_error("no TRestDetectorReadoutManager in the manager");
    auto readout = BuildReadout();
    detector.readout = readout.get();
    detector.readouts->fReadoutMap["mmTest"] = readout.release();  // owned by the readout manager
    return detector;
}

double TotalEnergy(const TRestDetectorHitsEvent& event) { return event.GetTotalEnergy(); }

}  // namespace

// ============================================================================
// Events
// ============================================================================

TEST(DetectorHitsEvent, AddHitsAndObservables) {
    TRestDetectorHitsEvent event;
    event.SetID(3);
    event.AddHit(1, 2, 3, 4, 0.5);
    event.AddHit(XYZVector(3, 2, 7), 4, 1.5);
    event.AddHit(5, kNaN, 10, 2, 0, HitType::XZ, 1000020040);  // alpha

    ASSERT_EQ(event.GetNumberOfHits(), 3u);
    EXPECT_EQ(event.GetParticleID(0), 0);  // unknown by default
    EXPECT_EQ(event.GetParticleID(2), 1000020040);
    EXPECT_DOUBLE_EQ(event.GetX(1), 3);
    EXPECT_DOUBLE_EQ(event.GetTime(1), 1.5);
    EXPECT_EQ(event.GetType(0), HitType::XYZ);
    EXPECT_EQ(event.GetType(2), HitType::XZ);
    EXPECT_TRUE(std::isnan(event.GetY(2)));
    EXPECT_DOUBLE_EQ(event.GetTotalEnergy(), 10);
    EXPECT_DOUBLE_EQ(event.GetMaximumHitEnergy(), 4);
    EXPECT_DOUBLE_EQ(event.GetMinimumHitEnergy(), 2);
    // Energy weighted means of the hits defining each coordinate
    EXPECT_NEAR(event.GetMeanPositionX(), (1 * 4 + 3 * 4 + 5 * 2) / 10., 1e-6);
    EXPECT_NEAR(event.GetMeanPositionY(), 2, 1e-6);
    EXPECT_NEAR(event.GetMeanPositionZ(), (3 * 4 + 7 * 4 + 10 * 2) / 10., 1e-6);
    EXPECT_EQ(event.GetXZHits()->GetNumberOfHits(), 1u);
    EXPECT_EQ(event.GetXYZHits()->GetNumberOfHits(), 2u);
    EXPECT_EQ(event.GetYZHits()->GetNumberOfHits(), 0u);

    TRestDetectorHitsEvent copy;
    copy.CopyFrom(&event);
    EXPECT_EQ(copy.GetID(), 3);
    ASSERT_EQ(copy.GetNumberOfHits(), 3u);
    EXPECT_DOUBLE_EQ(copy.GetZ(2), 10);
    EXPECT_EQ(copy.GetParticleID(2), 1000020040);

    event.Initialize();
    EXPECT_EQ(event.GetNumberOfHits(), 0u);
    EXPECT_EQ(event.GetID(), 0);
    EXPECT_EQ(copy.GetNumberOfHits(), 3u);
}

TEST(DetectorHitsEvent, SortAndEventSeed) {
    TRestDetectorHitsEvent event;
    event.AddHit(0, 0, 0, 1, 0, HitType::XYZ, 11);
    event.AddHit(0, 0, 1, 3, 0, HitType::XYZ, 13);
    event.AddHit(0, 0, 2, 2, 0, HitType::XYZ, 2212);
    event.Sort();  // descending energy by default, the particles with their hits
    EXPECT_DOUBLE_EQ(event.GetEnergy(0), 3);
    EXPECT_DOUBLE_EQ(event.GetEnergy(2), 1);
    EXPECT_EQ(event.GetParticleID(0), 13);
    EXPECT_EQ(event.GetParticleID(2), 11);
    EXPECT_TRUE(event.GetHits()->isSortedByEnergy());

    event.Sort([&event](int a, int b) { return event.GetZ(a) < event.GetZ(b); });
    EXPECT_DOUBLE_EQ(event.GetZ(0), 0);
    EXPECT_DOUBLE_EQ(event.GetZ(2), 2);
    EXPECT_EQ(event.GetParticleID(2), 2212);

    // The seed of an event: deterministic, never 0, different for other events
    TRestDetectorHitsEvent other;
    event.SetID(1);
    other.SetID(2);
    const auto seed = TRestEvent::EventSeed(17, event);
    EXPECT_EQ(seed, TRestEvent::EventSeed(17, event));
    EXPECT_NE(seed, 0u);
    EXPECT_NE(seed, TRestEvent::EventSeed(17, other));
    EXPECT_NE(seed, TRestEvent::EventSeed(18, event));
}

TEST(DetectorSignalEvent, AddRemoveAndAccess) {
    TRestDetectorSignalEvent event;
    event.AddSignal(10, {{0, 1}, {1, 2}, {2, 3}});
    event.AddSignal(20, std::vector<double>{0, 1}, std::vector<double>{5, 6});
    event.AddSignal(10, {{0, 100}});                                         // repeated ID: ignored
    event.AddSignal(30, std::vector<double>{0, 1}, std::vector<double>{1});  // size mismatch: ignored

    ASSERT_EQ(event.GetNumberOfSignals(), 2);
    EXPECT_TRUE(event.HasSignalID(20));
    EXPECT_FALSE(event.HasSignalID(30));
    EXPECT_EQ(event.GetSignal(0).GetSignalID(), 10);
    EXPECT_EQ(event.GetSignal(0).GetNPoints(), 3);
    EXPECT_EQ(event.GetSignalByID(20).GetNPoints(), 2);
    EXPECT_DOUBLE_EQ(event.GetSignalByID(20).GetPoint(1).second, 6);
    EXPECT_THROW(event.GetSignalByID(99), std::runtime_error);
    EXPECT_THROW(event.GetSignal(5), std::out_of_range);

    event.RemoveSignalByID(10);
    ASSERT_EQ(event.GetNumberOfSignals(), 1);
    EXPECT_EQ(event.GetSignal(0).GetSignalID(), 20);
    EXPECT_EQ(event.GetSignal(0).GetCharge(), (std::vector<double>{5, 6}));
    EXPECT_EQ(event.GetSignal(0).GetTime(), (std::vector<double>{0, 1}));

    TRestDetectorSignalEvent copy;
    copy.CopyFrom(&event);
    ASSERT_EQ(copy.GetNumberOfSignals(), 1);
    EXPECT_DOUBLE_EQ(copy.GetSignalByID(20).GetPoint(0).second, 5);
}

TEST(DetectorSignalEvent, SignalOperations) {
    TRestDetectorSignalEvent event;
    event.AddSignal(1, {{2, 20}, {0, 10}});
    event.AddSignal(2, {{0, 1}, {2, 3}});

    auto& first = event.GetSignalByID(1);
    // A point added to the first signal does not change the second one
    first.AddPoint(1, 15);
    EXPECT_EQ(first.GetNPoints(), 3);
    EXPECT_EQ(event.GetSignalByID(2).GetCharge(), (std::vector<double>{1, 3}));

    first.Sort();
    EXPECT_EQ(first.GetTime(), (std::vector<double>{0, 1, 2}));
    EXPECT_EQ(first.GetCharge(), (std::vector<double>{10, 15, 20}));
    EXPECT_EQ(first.GetTimeIndex(2), 2);
    EXPECT_EQ(first.GetTimeIndex(5), -1);

    first.IncreaseAmplitude(1, 5);  // existing time
    first.IncreaseAmplitude(3, 4);  // new point
    EXPECT_EQ(first.GetNPoints(), 4);
    EXPECT_DOUBLE_EQ(first.GetPoint(1).second, 20);
    EXPECT_DOUBLE_EQ(first.GetPoint(3).second, 4);

    first.MultiplySignalBy(2);
    first.AddOffset(-1);
    EXPECT_EQ(first.GetCharge(), (std::vector<double>{19, 39, 39, 7}));

    // Addition point by point, only with the same times
    auto& second = event.GetSignalByID(2);
    TRestDetectorSignalEvent other;
    other.AddSignal(5, {{0, 1}, {2, 1}});
    other.AddSignal(6, {{0, 1}, {1, 1}});
    second.SignalAddition(&other.GetSignalByID(5));
    EXPECT_EQ(second.GetCharge(), (std::vector<double>{2, 4}));
    second.SignalAddition(&other.GetSignalByID(6));  // different times: nothing added
    EXPECT_EQ(second.GetCharge(), (std::vector<double>{2, 4}));

    const TGraph graph = second.GetGraph();
    EXPECT_EQ(graph.GetN(), 2);
}

// ============================================================================
// Readout and medium
// ============================================================================

TEST(DetectorReadout, MicromegasPositionToChannel) {
    const auto readout = BuildReadout();
    ASSERT_NE(readout->GetGeoManager(), nullptr);

    for (int n = 0; n < kNChannels; n++) {
        for (int nPix : {0, kNChannels - 1}) {
            const auto x = XStripPixel(n, nPix);
            EXPECT_EQ(readout->GetChannelFromPosition(x.X(), x.Y(), x.Z()), XStrip(n))
                << "X strip " << n << ", pixel " << nPix;
            const auto y = YStripPixel(n, nPix);
            EXPECT_EQ(readout->GetChannelFromPosition(y.X(), y.Y(), y.Z()), YStrip(n))
                << "Y strip " << n << ", pixel " << nPix;
        }
    }
    // Outside the plane, and in the drift volume (no channel)
    EXPECT_EQ(readout->GetChannelFromPosition(100, 100, 0), -1);
    EXPECT_EQ(readout->GetChannelFromPosition(0, 0, 10), -1);

    // The position in the pixel (squares of pitch / sqrt(2) rotated 45 degrees): 0.1 mm from the centre
    // along a side of the pixel, also with the plane rotated by 30 degrees about its normal
    const double half = kPitch / std::sqrt(8.), d = 0.1;
    YAML::Node rotatedNode = ReadoutNode();
    rotatedNode["readoutParameters"]["globalRotation"] = "30 deg";
    TRestMicromegasReadout rotated("mmTest", rotatedNode);
    rotated.BuildGeometry();
    const double angle = 30 * std::acos(-1.) / 180;
    for (const bool rotate : {false, true}) {
        const TRestDetectorReadout& plane = rotate ? static_cast<TRestDetectorReadout&>(rotated) : *readout;
        const auto centre = YStripPixel(1, 2);
        const XYZVector inPlane = centre + XYZVector(d, d, 0) / std::sqrt(2.);
        const XYZVector point =
            rotate ? XYZVector(inPlane.X() * std::cos(angle) - inPlane.Y() * std::sin(angle),
                               inPlane.X() * std::sin(angle) + inPlane.Y() * std::cos(angle), 0)
                   : inPlane;
        TRestDetectorReadout::PixelPosition pixel;
        ASSERT_TRUE(plane.GetPixelPosition(point, pixel)) << "rotated " << rotate;
        EXPECT_EQ(pixel.daqID, YStrip(1));
        EXPECT_NEAR(std::abs(pixel.local.X()) + std::abs(pixel.local.Y()), d, 1e-9);
        EXPECT_NEAR(std::abs(pixel.local.X() * pixel.local.Y()), 0, 1e-9);
        EXPECT_NEAR(std::abs(pixel.local.Z()), 0, 1e-9);
        EXPECT_NEAR(pixel.halfSize.X(), half, 1e-9);
        EXPECT_NEAR(pixel.halfSize.Y(), half, 1e-9);
    }
    TRestDetectorReadout::PixelPosition pixel;
    EXPECT_FALSE(readout->GetPixelPosition(XYZVector(100, 100, 0), pixel));
    EXPECT_EQ(pixel.daqID, -1);
}

TEST(DetectorReadout, MicromegasChannelPositionsAndDecoding) {
    const auto readout = BuildReadout();

    const auto xStrip = readout->GetPositionFromChannel(XStrip(2));
    EXPECT_NEAR(xStrip.X(), XStripPosition(2), 1e-9);
    EXPECT_TRUE(std::isnan(xStrip.Y()));
    EXPECT_TRUE(std::isnan(xStrip.Z()));
    const auto yStrip = readout->GetPositionFromChannel(YStrip(1));
    EXPECT_NEAR(yStrip.Y(), YStripPosition(1), 1e-9);
    EXPECT_TRUE(std::isnan(yStrip.X()));
    EXPECT_TRUE(std::isnan(readout->GetPositionFromChannel(999).X()));

    EXPECT_TRUE(readout->IsChannelConnected(YStrip(0)));
    EXPECT_FALSE(readout->IsChannelConnected(999));
    EXPECT_EQ(readout->GetPhysicalID(XStrip(1)), kNChannels + 1);
    EXPECT_EQ(readout->GetPhysicalID(999), -1);
    EXPECT_NE(readout->GetDecodingAsString().find(std::to_string(XStrip(2)) + "\t" +
                                                  std::to_string(kNChannels + 2)),
              std::string::npos);

    // Neighbours in the same plane: the last Y strip and the first X strip are not adjacent
    EXPECT_EQ(readout->GetAdjacentChannels(YStrip(1)), std::make_pair(YStrip(0), YStrip(2)));
    EXPECT_EQ(readout->GetAdjacentChannels(YStrip(kNChannels - 1)).second, -1);
    EXPECT_EQ(readout->GetAdjacentChannels(XStrip(0)).first, -1);
    EXPECT_EQ(readout->GetAdjacentChannels(XStrip(0)).second, XStrip(1));
    EXPECT_EQ(readout->GetAdjacentChannels(999), std::make_pair(-1, -1));
}

TEST(DetectorReadout, MicromegasDriftGeometry) {
    const auto readout = BuildReadout();
    EXPECT_DOUBLE_EQ(readout->GetDistanceToReadout(XYZVector(1, 2, 30)), 30);
    EXPECT_DOUBLE_EQ(readout->GetDistanceToReadout(XYZVector(1, 2, -5)), -5);
    EXPECT_NEAR(readout->GetMaxDistanceToReadout(XYZVector(0, 0, 0)), kThickness / 2 + kGasThickness, 1e-9);
    // Per channel: NaN for a channel not in the decoding
    EXPECT_DOUBLE_EQ(readout->GetDistanceToReadout(XStrip(0), XYZVector(0, 0, 7)), 7);
    EXPECT_TRUE(std::isnan(readout->GetDistanceToReadout(999, XYZVector(0, 0, 7))));

    const auto normal = readout->GetReadoutNormal(XYZVector(1, 2, 30));
    EXPECT_NEAR(normal.Z(), 1, 1e-12);
    EXPECT_NEAR(normal.X(), 0, 1e-12);
    const auto projected = readout->ProjectToReadout(XYZVector(1, 2, 30));
    EXPECT_NEAR(projected.X(), 1, 1e-12);
    EXPECT_NEAR(projected.Y(), 2, 1e-12);
    EXPECT_NEAR(projected.Z(), 0, 1e-12);

    // The normal sets the drift side
    const auto inverted = BuildReadout("[0, 0, -1]");
    EXPECT_DOUBLE_EQ(inverted->GetDistanceToReadout(XYZVector(0, 0, -20)), 20);
    EXPECT_NEAR(inverted->GetReadoutNormal(XYZVector(0, 0, -20)).Z(), -1, 1e-12);

    // The normal must be along an axis
    EXPECT_THROW(TRestMicromegasReadout("bad", ReadoutNode("[1, 1, 0]")), std::runtime_error);

    // The X and Y pixels do not overlap (sampling check of the strips and the gas)
    EXPECT_EQ(readout->CheckOverlaps(), 0);
}

TEST(DetectorReadout, ExportAndImportWithManager) {
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/detector_test_readout.root";
    {
        const auto readout = BuildReadout();
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        ASSERT_TRUE(file && !file->IsZombie());
        ASSERT_TRUE(readout->Export(file.get(), "mmTest", "test_decoding"));
        file->Close();
    }

    const YAML::Node node = YAML::Load(
        "readout: mmTest\n"
        "readouts:\n"
        "  - {instanceName: mmTest, decodingName: test_decoding, inputFileName: \"" +
        fileName + "\"}\n");
    TRestDetectorReadoutManager manager("readouts", node);
    EXPECT_EQ(manager.fReadoutInfo.size(), 1u);
    manager.LoadReadout();
    ASSERT_TRUE(manager.HasReadout("mmTest"));
    EXPECT_EQ(manager.GetDataReadoutName(), "mmTest");

    std::string name;
    TRestDetectorReadout* readout = manager.SelectReadout(name, "test");
    ASSERT_NE(readout, nullptr);
    EXPECT_EQ(name, "mmTest");
    EXPECT_EQ(readout->GetClassName(), "TRestMicromegasReadout");
    const auto pixel = XStripPixel(1, 2);
    EXPECT_EQ(readout->GetChannelFromPosition(pixel.X(), pixel.Y(), pixel.Z()), XStrip(1));
    EXPECT_NEAR(readout->GetPositionFromChannel(YStrip(3)).Y(), YStripPosition(3), 1e-9);

    std::string missing = "other";
    EXPECT_THROW(manager.SelectReadout(missing, "test"), std::runtime_error);
}

TEST(DetectorReadout, AddDecodingToAFile) {
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/detector_test_add_decoding.root";
    {
        const auto readout = BuildReadout();
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        ASSERT_TRUE(readout->Export(file.get(), "mmTest", "nominal"));
    }
    // Another decoding (DAQ channel = 200 + physical ID), and one with a physical ID of no strip
    const std::string other = std::string(TESTS_WORK_DIR) + "/detector_test_other.dec";
    const std::string wrong = std::string(TESTS_WORK_DIR) + "/detector_test_wrong.dec";
    {
        std::ofstream out(other);
        for (int physicalID = 0; physicalID < 2 * kNChannels; physicalID++) {
            out << 200 + physicalID << "\t" << physicalID << "\n";
        }
        std::ofstream(wrong) << "300\t" << 2 * kNChannels << "\n";
    }
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "UPDATE"));
        EXPECT_TRUE(TRestDetectorReadout::AddDecoding(file.get(), "mmTest", "other", other));
        EXPECT_FALSE(TRestDetectorReadout::AddDecoding(file.get(), "mmTest", "wrong", wrong));
        EXPECT_FALSE(TRestDetectorReadout::AddDecoding(file.get(), "missing", "other", other));
    }
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    const auto names = TRestDetectorReadout::GetDecodingNames(file.get(), "mmTest");
    EXPECT_EQ(std::set<std::string>(names.begin(), names.end()), (std::set<std::string>{"nominal", "other"}));

    // The readout manager uses the added decoding (the geometry of the file)
    const YAML::Node node = YAML::Load(
        "readouts:\n  - {instanceName: mmTest, decodingName: other, inputFileName: \"" + fileName + "\"}\n");
    TRestDetectorReadoutManager manager("readouts", node);
    manager.LoadReadout();
    ASSERT_TRUE(manager.HasReadout("mmTest"));
    const auto pixel = XStripPixel(1, 2);
    EXPECT_EQ(manager.GetReadout("mmTest")->GetChannelFromPosition(pixel.X(), pixel.Y(), pixel.Z()),
              200 + kNChannels + 1);
}

TEST(DetectorReadoutManager, WritesOnlyTheDataReadout) {
    // Two readouts (the same geometry under two names, as the two planes of TREX-DM)
    const std::string input = std::string(TESTS_WORK_DIR) + "/detector_test_two_readouts.root";
    {
        const auto readout = BuildReadout();
        std::unique_ptr<TFile> file(TFile::Open(input.c_str(), "RECREATE"));
        ASSERT_TRUE(readout->Export(file.get(), "plane0", "nominal"));
        ASSERT_TRUE(readout->Export(file.get(), "plane1", "nominal"));
    }
    auto write = [&](const std::string& dataReadout, const std::string& output) {
        const std::string source = "decodingName: nominal, inputFileName: \"" + input + "\"";
        const YAML::Node node =
            YAML::Load("readout: \"" + dataReadout + "\"\nreadouts:\n" + "  - {instanceName: plane0, " +
                       source + "}\n" + "  - {instanceName: plane1, " + source + "}\n");
        TRestDetectorReadoutManager manager("readoutManager", node);
        manager.LoadReadout();
        EXPECT_TRUE(manager.HasReadout("plane0") && manager.HasReadout("plane1"));
        std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "RECREATE"));
        manager.WriteMetadata(file.get());
    };

    // With a data readout, only that one is stored
    const std::string onePlane = std::string(TESTS_WORK_DIR) + "/detector_test_data_readout.root";
    write("plane1", onePlane);
    {
        std::unique_ptr<TFile> file(TFile::Open(onePlane.c_str(), "READ"));
        EXPECT_TRUE(TRestMetadata::ReadMetadata(file.get(), "plane0").IsNull());
        EXPECT_TRUE(TRestMetadata::ReadMetadata(file.get(), "plane1").IsMap());
        const YAML::Node stored = TRestMetadata::ReadMetadata(file.get(), "readoutManager");
        ASSERT_TRUE(stored.IsMap() && stored["readouts"]);
        ASSERT_EQ(stored["readouts"].size(), 1u);
        EXPECT_EQ(stored["readouts"][0]["instanceName"].as<std::string>(), "plane1");
    }
    {
        std::unique_ptr<TFile> file(TFile::Open(onePlane.c_str(), "READ"));
        const YAML::Node stored = TRestMetadata::ReadMetadata(file.get(), "readoutManager");
        TRestDetectorReadoutManager manager("readoutManager", stored);
        manager.LoadFromFile(file.get());
        EXPECT_FALSE(manager.HasReadout("plane0"));
        EXPECT_TRUE(manager.HasReadout("plane1"));
    }

    // Without a data readout, all of them
    const std::string allPlanes = std::string(TESTS_WORK_DIR) + "/detector_test_all_readouts.root";
    write("", allPlanes);
    std::unique_ptr<TFile> file(TFile::Open(allPlanes.c_str(), "READ"));
    EXPECT_TRUE(TRestMetadata::ReadMetadata(file.get(), "plane0").IsMap());
    EXPECT_TRUE(TRestMetadata::ReadMetadata(file.get(), "plane1").IsMap());
}

namespace {
/// GDML of the TRestGDMLReadout tests (mm): two copies of a module, at x = -200 and at x = 200 rotated 180
/// degrees around X, each one with a scintillator (local z = -20) and its light guide (local z = 60)
std::string WriteVetoGeometry() {
    const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string path = std::string(TESTS_WORK_DIR) + "/veto_" + test->name() + ".gdml";
    std::ofstream out(path);
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           "<gdml>\n"
           "<define/>\n"
           "<materials>\n"
           "  <element name=\"Hydrogen\" formula=\"H\" Z=\"1\"><atom value=\"1.008\"/></element>\n"
           "  <material name=\"Plastic\" state=\"solid\"><D value=\"1.03\"/>"
           "<fraction n=\"1\" ref=\"Hydrogen\"/></material>\n"
           "</materials>\n"
           "<solids>\n"
           "  <box name=\"worldBox\" x=\"1000\" y=\"1000\" z=\"1000\" lunit=\"mm\"/>\n"
           "  <box name=\"moduleBox\" x=\"100\" y=\"100\" z=\"200\" lunit=\"mm\"/>\n"
           "  <box name=\"scintillatorBox\" x=\"50\" y=\"50\" z=\"100\" lunit=\"mm\"/>\n"
           "  <box name=\"guideBox\" x=\"50\" y=\"50\" z=\"20\" lunit=\"mm\"/>\n"
           "</solids>\n"
           "<structure>\n"
           "  <volume name=\"scintillatorVolume\"><materialref ref=\"Plastic\"/>"
           "<solidref ref=\"scintillatorBox\"/></volume>\n"
           "  <volume name=\"lightGuideVolume\"><materialref ref=\"Plastic\"/>"
           "<solidref ref=\"guideBox\"/></volume>\n"
           "  <volume name=\"moduleVolume\"><materialref ref=\"Plastic\"/><solidref ref=\"moduleBox\"/>\n"
           "    <physvol name=\"scintillator\"><volumeref ref=\"scintillatorVolume\"/>"
           "<position name=\"scintillatorPosition\" z=\"-20\" unit=\"mm\"/></physvol>\n"
           "    <physvol name=\"lightGuide\"><volumeref ref=\"lightGuideVolume\"/>"
           "<position name=\"guidePosition\" z=\"60\" unit=\"mm\"/></physvol>\n"
           "  </volume>\n"
           "  <volume name=\"World\"><materialref ref=\"Plastic\"/><solidref ref=\"worldBox\"/>\n"
           "    <physvol name=\"module0\"><volumeref ref=\"moduleVolume\"/>"
           "<position name=\"module0Position\" x=\"-200\" unit=\"mm\"/></physvol>\n"
           "    <physvol name=\"module1\"><volumeref ref=\"moduleVolume\"/>"
           "<position name=\"module1Position\" x=\"200\" unit=\"mm\"/>"
           "<rotation name=\"module1Rotation\" x=\"180\" unit=\"deg\"/></physvol>\n"
           "  </volume>\n"
           "</structure>\n"
           "<setup name=\"Default\" version=\"1.0\"><world ref=\"World\"/></setup>\n"
           "</gdml>\n";
    return path;
}

/// A TRestGDMLReadout of the veto geometry: the scintillators are the channels 10 + module (DAQ 300 and
/// 301), with readoutInfo extra keys appended to the entry
std::unique_ptr<TRestGDMLReadout> BuildVetoReadout(const std::string& readoutInfo) {
    const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string decoding = std::string(TESTS_WORK_DIR) + "/veto_" + test->name() + ".dec";
    std::ofstream(decoding) << "# daqChannel physicalID\n300\t10\n301\t11\n";
    auto readout = std::make_unique<TRestGDMLReadout>(
        "veto", YAML::Load("name: veto\n"
                           "decodingFile: \"" +
                           decoding +
                           "\"\n"
                           "readoutParameters:\n"
                           "  gdmlFileName: \"" +
                           WriteVetoGeometry() +
                           "\"\n"
                           "  readoutInfo:\n"
                           "    - {volume: \"module(\\\\d)/scintillator\", offset: 10" +
                           readoutInfo + "}\n"));
    readout->BuildGeometry();
    return readout;
}
}  // namespace

TEST(DetectorReadout, GDMLChannelsAndReadoutVolumes) {
    const auto readout = BuildVetoReadout(", readoutVolume: lightGuide");
    // Physical ID offset + the module number captured by the expression
    EXPECT_EQ(readout->GetPhysicalID(300), 10);
    EXPECT_EQ(readout->GetPhysicalID(301), 11);
    ASSERT_EQ(readout->fPhysicalIDToReadoutPath.size(), 2u);

    // Channel volumes and their readout volumes in the global frame (the second module rotated)
    const auto channel0 = readout->GetPositionFromChannel(300);
    const auto channel1 = readout->GetPositionFromChannel(301);
    EXPECT_NEAR((channel0 - XYZVector(-200, 0, -20)).R(), 0, 1e-9);
    EXPECT_NEAR((channel1 - XYZVector(200, 0, 20)).R(), 0, 1e-9);
    EXPECT_NEAR((readout->GetReadoutPosition(300) - XYZVector(-200, 0, 60)).R(), 0, 1e-9);
    EXPECT_NEAR((readout->GetReadoutPosition(301) - XYZVector(200, 0, -60)).R(), 0, 1e-9);
    EXPECT_NEAR(readout->GetDistanceToReadout(300, XYZVector(-200, 0, -20)), 80, 1e-9);
    EXPECT_NEAR(readout->GetDistanceToReadout(301, XYZVector(200, 0, 0)), 60, 1e-9);
    // The hits are not drifted
    EXPECT_TRUE(std::isnan(readout->GetDistanceToReadout(XYZVector(-200, 0, -20))));

    // Positions in the scintillators give their channels; outside them (light guide, world), none
    EXPECT_EQ(readout->GetChannelFromPosition(-200, 10, -30), 300);
    EXPECT_EQ(readout->GetChannelFromPosition(200, -10, 30), 301);
    EXPECT_EQ(readout->GetChannelFromPosition(-200, 0, 60), -1);
    EXPECT_EQ(readout->GetChannelFromPosition(0, 0, 0), -1);
}

TEST(DetectorReadout, GDMLWithoutReadoutVolume) {
    const auto readout = BuildVetoReadout("");
    EXPECT_EQ(readout->GetPhysicalID(301), 11);
    EXPECT_TRUE(readout->fPhysicalIDToReadoutPath.empty());
    EXPECT_TRUE(std::isnan(readout->GetReadoutPosition(300).X()));
    EXPECT_TRUE(std::isnan(readout->GetDistanceToReadout(300, XYZVector(-200, 0, -20))));

    // Errors: a readout volume not in the module, an invalid expression
    EXPECT_THROW(BuildVetoReadout(", readoutVolume: notAVolume"), std::runtime_error);
    EXPECT_THROW(BuildVetoReadout(", readoutVolume: \"(lightGuide\""), std::runtime_error);
}

TEST(DetectorMedium, FixedChargeResponse) {
    const YAML::Node node = YAML::Load(
        "name: ArIso\n"
        "pressure: 1.4 bar\n"
        "charge: {wValue: 26.2 eV, fanoFactor: 0.2, driftVelocity: 1.2 mm/us, longDiff: 0.03, "
        "transDiff: 0.04, attachment: 0.1}\n"
        "light: {attenuationLength: 3 m, effectiveLightSpeed: 150 mm/us}\n");
    TRestDetectorMedium medium("medium", node);
    EXPECT_EQ(medium.GetName(), "ArIso");
    EXPECT_TRUE(medium.IsGas());
    EXPECT_DOUBLE_EQ(medium.GetPressure(), 1.4);

    const auto& charge = medium.GetChargeResponse("test");
    EXPECT_NEAR(charge.GetWValue(), 0.0262, 1e-12);  // keV
    EXPECT_DOUBLE_EQ(charge.GetFanoFactor(), 0.2);
    EXPECT_FALSE(charge.HasTransportTables());
    EXPECT_DOUBLE_EQ(charge.GetDriftVelocity(), 1.2);
    EXPECT_DOUBLE_EQ(charge.GetLongitudinalDiffusion(), 0.03);
    EXPECT_DOUBLE_EQ(charge.GetTransversalDiffusion(), 0.04);
    EXPECT_DOUBLE_EQ(charge.GetAttachment(), 0.1);
    // Without tables the fixed values for any field, no Townsend and only electrons
    EXPECT_DOUBLE_EQ(charge.GetDriftVelocity(50.), 1.2);
    EXPECT_DOUBLE_EQ(charge.GetTownsend(50.), 0);
    EXPECT_THROW(charge.GetDriftVelocity(50., "ion"), std::runtime_error);
    EXPECT_THROW(charge.GetGainTable(0.05), std::runtime_error);

    const auto& light = medium.GetLightResponse("test");
    EXPECT_DOUBLE_EQ(light.GetAttenuationLength(), 3000);
    EXPECT_DOUBLE_EQ(light.GetEffectiveLightSpeed(), 150);
    // Without scintillation profiles the photons are emitted at once
    TRandom3 random(3);
    EXPECT_EQ(light.GetScintillationProfile(11), nullptr);
    EXPECT_DOUBLE_EQ(light.GetScintillationTime(11, random), 0);

    // Scintillation profiles per particle: default (fast) and nuclei (slow, no rise time)
    const TRestDetectorLightResponse scintillator(
        "light", YAML::Load("{lightYield: 10, scintillation: ["
                            "{components: [{fraction: 1, riseTime: 1 ns, decayTime: 2 ns}]},"
                            "{particles: [alpha, ion], components: [{fraction: 0.5, decayTime: 10 ns},"
                            "{fraction: 0.5, decayTime: 100 ns}]}]}"));
    const auto* fast = scintillator.GetScintillationProfile(11);
    const auto* slow = scintillator.GetScintillationProfile(1000180400);  // Ar40 recoil
    ASSERT_NE(fast, nullptr);
    ASSERT_NE(slow, nullptr);
    EXPECT_TRUE(fast->IsDefault());
    EXPECT_FALSE(slow->IsDefault());
    EXPECT_EQ(scintillator.GetScintillationProfile(0), fast);           // unknown particle (data)
    EXPECT_EQ(scintillator.GetScintillationProfile(1000020040), slow);  // alpha
    // Mean emission time: decay + rise (us)
    const int n = 20000;
    double fastMean = 0, slowMean = 0;
    for (int i = 0; i < n; i++) {
        const double t = scintillator.GetScintillationTime(13, random);
        ASSERT_GE(t, 0);
        fastMean += t / n;
        slowMean += scintillator.GetScintillationTime(1000020040, random) / n;
    }
    EXPECT_NEAR(fastMean, 0.003, 0.0002);
    EXPECT_NEAR(slowMean, 0.055, 0.003);
    // Errors: fractions not adding up to 1, unknown particle, two defaults, negative time
    EXPECT_THROW(TRestDetectorLightResponse(
                     "l", YAML::Load("{scintillation: [{components: [{fraction: 0.5, decayTime: 2 ns}]}]}")),
                 std::runtime_error);
    EXPECT_THROW(TRestDetectorLightResponse("l", YAML::Load("{scintillation: [{particles: [electron], "
                                                            "components: [{decayTime: 2 ns}]}]}")),
                 std::runtime_error);
    EXPECT_THROW(TRestDetectorLightResponse("l", YAML::Load("{scintillation: ["
                                                            "{components: [{decayTime: 2 ns}]},"
                                                            "{components: [{decayTime: 3 ns}]}]}")),
                 std::runtime_error);
    EXPECT_THROW(
        TRestDetectorLightResponse("l", YAML::Load("{scintillation: [{components: [{decayTime: -2 ns}]}]}")),
        std::runtime_error);

    // Without name, the name of the instance (as any metadata)
    EXPECT_EQ(TRestDetectorMedium("m", YAML::Load("{charge: {wValue: 26 eV}}")).GetName(), "m");
    // Errors: negative parameters, pressure of a solid
    EXPECT_THROW(TRestDetectorMedium("m", YAML::Load("{name: a, charge: {wValue: -1}}")), std::runtime_error);
    EXPECT_THROW(TRestDetectorMedium("m", YAML::Load("{name: a, state: solid, pressure: 2 bar}")),
                 std::runtime_error);
    TRestDetectorMedium noCharge("m", YAML::Load("{name: a}"));
    EXPECT_EQ(noCharge.GetChargeResponse(), nullptr);
    EXPECT_THROW(noCharge.GetChargeResponse("test"), std::runtime_error);
}

TEST(DetectorMedium, DiffuseElectron) {
    TRandom3 random(5);
    // Without diffusion and attachment the electron stays where it is
    XYZVector position(1, 2, 30);
    ASSERT_TRUE(
        TRestDetectorChargeResponse::DiffuseElectron(position, XYZVector(0, 0, 1), 30, 0, 0, 0, random));
    EXPECT_NEAR((position - XYZVector(1, 2, 30)).R(), 0, 1e-12);

    // Gaussian spread sigma = D sqrt(distance), D in sqrt(cm) and distance in mm
    const double distance = 40, longDiff = 0.03, transDiff = 0.05;
    const int n = 20000;
    double sumX2 = 0, sumZ = 0, sumZ2 = 0;
    for (int i = 0; i < n; i++) {
        XYZVector p(0, 0, distance);
        ASSERT_TRUE(TRestDetectorChargeResponse::DiffuseElectron(p, XYZVector(0, 0, 1), distance, longDiff,
                                                                 transDiff, 0, random));
        sumX2 += p.X() * p.X();
        sumZ += p.Z();
        sumZ2 += (p.Z() - distance) * (p.Z() - distance);
    }
    const double sigmaT = 10 * transDiff * std::sqrt(distance / 10);
    const double sigmaL = 10 * longDiff * std::sqrt(distance / 10);
    EXPECT_NEAR(std::sqrt(sumX2 / n), sigmaT, 0.03 * sigmaT);
    EXPECT_NEAR(std::sqrt(sumZ2 / n), sigmaL, 0.03 * sigmaL);
    EXPECT_NEAR(sumZ / n, distance, 0.02);

    // Attachment: survival exp(-attachment * distance(cm))
    int survived = 0;
    for (int i = 0; i < n; i++) {
        XYZVector p(0, 0, distance);
        survived +=
            TRestDetectorChargeResponse::DiffuseElectron(p, XYZVector(0, 0, 1), distance, 0, 0, 0.2, random);
    }
    EXPECT_NEAR(survived / double(n), std::exp(-0.2 * distance / 10), 0.015);
}

// ============================================================================
// Hits processes
// ============================================================================

TEST(DetectorHitsDriftProcess, ElectronsWithoutFluctuations) {
    auto detector = MakeDetector(kChargeNoDiffusion);
    TRestDetectorHitsDriftProcess drift("drift", YAML::Load("{fluctuations: none, maxHits: 1000, seed: 3}"));
    drift.SetManager(detector.manager.get());
    drift.Initialize();
    EXPECT_EQ(drift.fMediumName, "ArIso");

    TRestDetectorHitsEvent input, output;
    input.SetID(1);
    input.AddHit(1, 2, 20, 5.9, 0.5);
    input.AddHit(0, 0, -1, 1);  // behind the readout plane
    input.AddHit(0, 0, 60, 1);  // beyond the cathode
    ASSERT_TRUE(drift.ProcessEvent(input, output));

    EXPECT_EQ(output.GetID(), 1);
    EXPECT_DOUBLE_EQ(drift.fDriftVelocity, 10);  // from the medium
    // energy / W electrons, each one with the energy W, on the readout plane below the hit (no diffusion)
    // at its arrival time 0.5 + 20 / 10
    EXPECT_EQ(output.GetNumberOfHits(), (size_t)std::lround(5.9 / 0.0262));
    EXPECT_NEAR(TotalEnergy(output), 5.9, 1e-3);
    for (size_t n = 0; n < output.GetNumberOfHits(); n++) {
        EXPECT_NEAR(output.GetX(n), 1, 1e-5);
        EXPECT_NEAR(output.GetY(n), 2, 1e-5);
        EXPECT_NEAR(output.GetZ(n), 0, 1e-5);
        EXPECT_NEAR(output.GetTime(n), 2.5, 1e-6);
    }

    // maxHits: each electron carries the energy of several
    TRestDetectorHitsEvent big, bigOutput;
    big.AddHit(0, 0, 10, 100);
    ASSERT_TRUE(drift.ProcessEvent(big, bigOutput));
    EXPECT_EQ(bigOutput.GetNumberOfHits(), 1000u);
    EXPECT_NEAR(TotalEnergy(bigOutput), 100, 1e-2);

    // Only hits outside the drift volume: rejected
    TRestDetectorHitsEvent outside, outsideOutput;
    outside.AddHit(0, 0, -3, 5);
    EXPECT_FALSE(drift.ProcessEvent(outside, outsideOutput));
}

TEST(DetectorHitsDriftProcess, DiffusionIsSeeded) {
    auto detector = MakeDetector(kChargeDiffusion);
    const YAML::Node node = YAML::Load("{fluctuations: fano, seed: 11}");
    TRestDetectorHitsDriftProcess drift("drift", node);
    drift.SetManager(detector.manager.get());
    drift.Initialize();
    EXPECT_DOUBLE_EQ(drift.fLongitudinalDiffusion, 0.03);
    EXPECT_DOUBLE_EQ(drift.fTransversalDiffusion, 0.04);

    TRestDetectorHitsEvent input, first, second;
    input.SetID(7);
    input.AddHit(0, 0, 30, 5.9);
    ASSERT_TRUE(drift.ProcessEvent(input, first));
    ASSERT_TRUE(drift.ProcessEvent(input, second));

    // Same seed and event: the same electrons
    ASSERT_EQ(first.GetNumberOfHits(), second.GetNumberOfHits());
    for (size_t n = 0; n < first.GetNumberOfHits(); n++) {
        EXPECT_FLOAT_EQ(first.GetX(n), second.GetX(n));
        EXPECT_FLOAT_EQ(first.GetTime(n), second.GetTime(n));
    }
    // Fano fluctuation of the number of electrons and the diffusion around the hit: transversal in the
    // position on the plane, longitudinal in the arrival time (distance / 10, between the plane and the
    // cathode)
    EXPECT_NEAR(TotalEnergy(first), 5.9, 1.0);
    EXPECT_NEAR(first.GetMeanPositionX(), 0, 0.3);
    EXPECT_GT(first.GetSigmaX(), 0);
    double meanTime = 0;
    for (size_t n = 0; n < first.GetNumberOfHits(); n++) {
        EXPECT_NEAR(first.GetZ(n), 0, 1e-5);
        EXPECT_GE(first.GetTime(n), -1e-6);
        EXPECT_LE(first.GetTime(n), (kThickness / 2 + kGasThickness) / 10 + 1e-4);
        meanTime += first.GetTime(n) / first.GetNumberOfHits();
    }
    EXPECT_NEAR(meanTime, 3, 0.03);

    // Without diffusion: the whole hit on the readout plane below it at 30 / 10 us, with its charge
    auto attaching = MakeDetector(
        "{wValue: 26.2 eV, fanoFactor: 0.2, driftVelocity: 10, longDiff: 0.03, transDiff: 0.04, "
        "attachment: 0.1}");
    TRestDetectorHitsDriftProcess whole("drift", YAML::Load("{fluctuations: none, diffusion: false}"));
    whole.SetManager(attaching.manager.get());
    whole.Initialize();
    TRestDetectorHitsEvent projected;
    input.AddHit(0, 0, -1, 1);  // behind the readout plane: removed
    ASSERT_TRUE(whole.ProcessEvent(input, projected));
    ASSERT_EQ(projected.GetNumberOfHits(), 1u);
    EXPECT_NEAR(projected.GetX(0), 0, 1e-5);
    EXPECT_NEAR(projected.GetZ(0), 0, 1e-5);
    EXPECT_NEAR(projected.GetTime(0), 3, 1e-6);
    EXPECT_NEAR(projected.GetEnergy(0), 5.9 * std::exp(-0.1 * 3), 1e-4);  // mean attachment over 3 cm

    // Fano fluctuations need fanoFactor > 0
    auto noFano = MakeDetector("{wValue: 26.2 eV, driftVelocity: 10}");
    TRestDetectorHitsDriftProcess fano("drift", YAML::Load("{fluctuations: fano}"));
    fano.SetManager(noFano.manager.get());
    EXPECT_THROW(fano.Initialize(), std::runtime_error);
    TRestDetectorHitsDriftProcess noManager("drift", YAML::Load("{}"));
    EXPECT_THROW(noManager.Initialize(), std::runtime_error);
}

namespace {
/// Drift response of a grid x, y in [-10, 10] mm, z in [0, 50] mm, 2 x 2 x 2 nodes, 4 bins (no file)
TRestDetectorDriftResponse SmallDriftResponse(const std::string& responseFile) {
    return TRestDetectorDriftResponse(
        "drift", YAML::Load("{fieldMap: field.txt, responseFile: \"" + responseFile +
                            "\", xRange: [-10 mm, 10 mm], yRange: [-10 mm, 10 mm], zRange: [0 mm, 50 mm], "
                            "gridStep: [20 mm, 20 mm, 50 mm], positionBins: 4, timeBins: 4}"));
}

/// The 4 nodes of a layer k of the small response with `arrived` of 100 electrons at (u, v, t)
void SetLayer(TRestDetectorDriftResponse& response, int k, int arrived,
              const std::array<double, 3>& arrival) {
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            response.SetNode(response.NodeIndex(i, j, k), 100,
                             std::vector<std::array<double, 3>>(arrived, arrival));
        }
    }
}
}  // namespace

TEST(DetectorDriftResponse, WriteReadAndSample) {
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/DriftResponse_WriteReadAndSample.root";
    auto response = SmallDriftResponse(fileName);
    EXPECT_EQ(response.GetNumberOfNodes(), 8u);
    EXPECT_DOUBLE_EQ(response.GetNodePosition(response.NodeIndex(1, 0, 1)).X(), 10);
    EXPECT_DOUBLE_EQ(response.GetNodePosition(response.NodeIndex(1, 0, 1)).Z(), 50);
    // Readout plane axes: normal +z gives u = x, v = y
    const auto [u, v] = TRestDetectorDriftResponse::PlaneAxes({0, 0, 2});
    EXPECT_NEAR(u.X(), 1, 1e-12);
    EXPECT_NEAR(v.Y(), 1, 1e-12);

    // The histograms of a node over the range of its arrivals
    response.SetNode(0, 10, {{0, 0, 0}, {1, 1, 1}});
    const auto& node = response.fNodes[0];
    EXPECT_EQ(node.fArrived, 2);
    EXPECT_EQ(node.fPositionCounts[0], 1);
    EXPECT_EQ(node.fPositionCounts[3 + 4 * 3], 1);
    EXPECT_EQ(node.fTimeCounts[3], 1);
    // Its summary from the bin centres (0.125 and 0.875 mm, us): survival 2 / 10
    const auto summary = response.GetNodeSummary(0);
    EXPECT_NEAR(summary.fSurvival, 0.2, 1e-12);
    EXPECT_NEAR(summary.fMeanU, 0.5, 1e-12);
    EXPECT_NEAR(summary.fMeanV, 0.5, 1e-12);
    EXPECT_NEAR(summary.fMeanT, 0.5, 1e-12);
    EXPECT_NEAR(summary.fSigmaU, 0.375, 1e-12);
    EXPECT_NEAR(summary.fSigmaT, 0.375, 1e-12);

    // All the electrons of the readout side arrive (u = 1, v = -2 mm at 0 us), half of those of the cathode
    // side (at 5 us)
    SetLayer(response, 0, 100, {1, -2, 0});
    SetLayer(response, 1, 50, {1, -2, 5});
    response.fMedium = "ArIso";
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        response.Write(file.get(), YAML::Load("{program: test}"));
    }
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    const auto read = TRestDetectorDriftResponse::Read(file.get());
    EXPECT_EQ(read.fMedium, "ArIso");
    ASSERT_EQ(read.fNodes.size(), 8u);
    EXPECT_EQ(read.fNodes[read.NodeIndex(0, 1, 1)].fArrived, 50);
    const auto cathode = read.GetNodeSummary(read.NodeIndex(0, 1, 1));
    EXPECT_NEAR(cathode.fSurvival, 0.5, 1e-12);
    EXPECT_NEAR(cathode.fMeanU, 1, 1e-12);
    EXPECT_NEAR(cathode.fMeanV, -2, 1e-12);
    EXPECT_NEAR(cathode.fMeanT, 5, 1e-12);
    EXPECT_NEAR(cathode.fSigmaU, 0, 1e-12);
    // The summary in the tree, to draw the response
    auto* tree = file->Get<TTree>("DriftResponse");
    ASSERT_NE(tree, nullptr);
    EXPECT_EQ(tree->Draw("meanT", "z == 50 && x == -10 && y == 10", "goff"), 1);
    EXPECT_NEAR(tree->GetV1()[0], 5, 1e-12);

    EXPECT_TRUE(read.Contains({0, 0, 25}));
    EXPECT_FALSE(read.Contains({0, 0, 51}));
    EXPECT_FALSE(read.Contains({11, 0, 0}));
    EXPECT_NEAR(read.GetSurvival({0, 0, 0}), 1, 1e-12);
    EXPECT_NEAR(read.GetSurvival({3, -4, 50}), 0.5, 1e-12);
    EXPECT_NEAR(read.GetSurvival({0, 0, 25}), 0.75, 1e-12);
    EXPECT_EQ(read.GetSurvival({0, 0, 60}), 0);

    // At z = 25 mm a node of each layer, the cathode side with weight 0.5 * 0.5 / (0.5 + 0.25) = 1/3
    TRandom3 random(1);
    int late = 0;
    const int samples = 6000;
    for (int i = 0; i < samples; i++) {
        double du = 0, dv = 0, time = -1;
        ASSERT_TRUE(read.SampleArrival({0, 0, 25}, random, du, dv, time));
        EXPECT_DOUBLE_EQ(du, 1);
        EXPECT_DOUBLE_EQ(dv, -2);
        ASSERT_TRUE(time == 0 || time == 5) << time;
        if (time == 5) late++;
    }
    EXPECT_NEAR(static_cast<double>(late) / samples, 1. / 3, 0.03);
    double du, dv, time;
    EXPECT_FALSE(read.SampleArrival({0, 0, 60}, random, du, dv, time));

    // Configuration errors
    EXPECT_THROW(TRestDetectorDriftResponse("drift", YAML::Load("{fieldMap: f.txt, responseFile: r.root, "
                                                                "xRange: [1 mm, 0 mm], yRange: [0 mm, 1 mm], "
                                                                "zRange: [0 mm, 1 mm]}")),
                 std::runtime_error);
    const std::string twoSteps =
        "{fieldMap: f.txt, responseFile: r.root, xRange: [0 mm, 1 mm], yRange: [0 mm, 1 mm], "
        "zRange: [0 mm, 1 mm], gridStep: [1 mm, 1 mm]}";
    EXPECT_THROW(TRestDetectorDriftResponse("drift", YAML::Load(twoSteps)), std::runtime_error);
    const std::string manyElectrons =
        "{fieldMap: f.txt, responseFile: r.root, xRange: [0 mm, 1 mm], yRange: [0 mm, 1 mm], "
        "zRange: [0 mm, 1 mm], gridStep: [1 mm, 1 mm, 1 mm], electrons: 70000}";  // counts of 16 bits
    EXPECT_THROW(TRestDetectorDriftResponse("drift", YAML::Load(manyElectrons)), std::runtime_error);

    // gridStep: the fewest equal steps not larger than it (at least 2 nodes); x and y not given (the readout
    // area, set with SetRange): no nodes until they are set, and the file keeps them
    TRestDetectorDriftResponse onReadout(
        "drift", YAML::Load("{fieldMap: f.txt, responseFile: r.root, zRange: [0 mm, 10 mm], "
                            "gridStep: [3 mm, 30 mm, 4 mm]}"));
    EXPECT_FALSE(onReadout.HasRange(0));
    EXPECT_TRUE(onReadout.HasRange(2));
    EXPECT_EQ(onReadout.fGridNodes[2], 4);  // 10 mm in 3 steps of 3.33 mm
    EXPECT_EQ(onReadout.GetNumberOfNodes(), 0u);
    EXPECT_FALSE(onReadout.Contains({0, 0, 5}));
    onReadout.SetRange(0, -6, 6);
    onReadout.SetRange(1, -6, 6);
    EXPECT_EQ(onReadout.fGridNodes[0], 5);  // 12 mm in 4 steps of 3 mm (exact)
    EXPECT_EQ(onReadout.fGridNodes[1], 2);
    EXPECT_EQ(onReadout.GetNumberOfNodes(), 40u);
    EXPECT_DOUBLE_EQ(onReadout.GetNodePosition(onReadout.NodeIndex(1, 0, 0)).X(), -3);
    for (size_t index = 0; index < onReadout.GetNumberOfNodes(); index++) onReadout.SetNode(index, 1, {});
    const std::string onReadoutFile = std::string(TESTS_WORK_DIR) + "/DriftResponse_OnReadout.root";
    {
        std::unique_ptr<TFile> out(TFile::Open(onReadoutFile.c_str(), "RECREATE"));
        onReadout.Write(out.get(), YAML::Load("{program: test}"));
    }
    std::unique_ptr<TFile> in(TFile::Open(onReadoutFile.c_str(), "READ"));
    const auto readOnReadout = TRestDetectorDriftResponse::Read(in.get());
    EXPECT_EQ(readOnReadout.GetNumberOfNodes(), 40u);
    EXPECT_DOUBLE_EQ(readOnReadout.fXRange.first, -6);
    EXPECT_THROW(onReadout.SetRange(0, 1, 1), std::runtime_error);
}

TEST(DetectorHitsDriftProcess, ElectronsFromTheResponse) {
    auto detector = MakeDetector(kChargeNoDiffusion);
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/DriftResponse_Process.root";
    {
        // Every electron arrives displaced by (0.5, -0.25) mm (u = x, v = y for the normal +z) after 3 us
        auto response = SmallDriftResponse(fileName);
        SetLayer(response, 0, 100, {0.5, -0.25, 3});
        SetLayer(response, 1, 100, {0.5, -0.25, 3});
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        response.Write(file.get(), YAML::Load("{program: test}"));
    }
    TRestDetectorHitsDriftProcess drift(
        "drift", YAML::Load("{responseFile: \"" + fileName + "\", fluctuations: none, seed: 3}"));
    drift.SetManager(detector.manager.get());
    drift.Initialize();
    EXPECT_EQ(drift.fMediumName, "ArIso");

    TRestDetectorHitsEvent input, output;
    input.SetID(2);
    input.AddHit(1, 2, 20, 5.9, 0.5);
    input.AddHit(0, 0, 60, 1);  // out of the grid of the response
    ASSERT_TRUE(drift.ProcessEvent(input, output));
    EXPECT_EQ(output.GetID(), 2);
    // energy / W electrons on the readout plane, displaced, at 0.5 + 3 us
    EXPECT_EQ(output.GetNumberOfHits(), (size_t)std::lround(5.9 / 0.0262));
    EXPECT_NEAR(TotalEnergy(output), 5.9, 1e-3);
    for (size_t n = 0; n < output.GetNumberOfHits(); n++) {
        EXPECT_NEAR(output.GetX(n), 1.5, 1e-5);
        EXPECT_NEAR(output.GetY(n), 1.75, 1e-5);
        EXPECT_NEAR(output.GetZ(n), 0, 1e-5);
        EXPECT_NEAR(output.GetTime(n), 3.5, 1e-6);
    }

    // Only hits out of the grid: rejected
    TRestDetectorHitsEvent outside, outsideOutput;
    outside.AddHit(0, 0, 60, 5);
    EXPECT_FALSE(drift.ProcessEvent(outside, outsideOutput));

    // Errors: no response file, diffusion false with a response
    TRestDetectorHitsDriftProcess noFile("drift", YAML::Load("{responseFile: none.root}"));
    noFile.SetManager(detector.manager.get());
    EXPECT_ANY_THROW(noFile.Initialize());
    TRestDetectorHitsDriftProcess whole("drift",
                                        YAML::Load("{responseFile: \"" + fileName + "\", diffusion: false}"));
    whole.SetManager(detector.manager.get());
    EXPECT_ANY_THROW(whole.Initialize());
}

TEST(DetectorHitsSmearingProcess, OneGainPerEvent) {
    TRestDetectorHitsSmearingProcess smearing(
        "smearing", YAML::Load("{energyReference: 5.9 keV, resolutionReference: 15, seed: 42}"));
    smearing.Initialize();

    TRestDetectorHitsEvent input;
    input.AddHit(1, 2, 3, 2.0, 0.1);
    input.AddHit(4, 5, 6, 3.9, 0.2);

    // Deterministic with the seed and the same gain for all the hits; positions and times unchanged
    TRestDetectorHitsEvent first, second;
    input.SetID(1);
    ASSERT_TRUE(smearing.ProcessEvent(input, first));
    ASSERT_TRUE(smearing.ProcessEvent(input, second));
    ASSERT_EQ(first.GetNumberOfHits(), 2u);
    EXPECT_FLOAT_EQ(first.GetEnergy(0), second.GetEnergy(0));
    EXPECT_NEAR(first.GetEnergy(0) / 2.0, first.GetEnergy(1) / 3.9, 1e-5);
    EXPECT_FLOAT_EQ(first.GetX(1), 4);
    EXPECT_FLOAT_EQ(first.GetTime(1), 0.2);

    // Relative spread FWHM / 2.355 at the reference energy
    const int n = 4000;
    double sum = 0, sum2 = 0;
    for (int i = 0; i < n; i++) {
        TRestDetectorHitsEvent output;
        input.SetID(i);
        ASSERT_TRUE(smearing.ProcessEvent(input, output));
        const double gain = TotalEnergy(output) / 5.9;
        sum += gain;
        sum2 += gain * gain;
    }
    const double mean = sum / n;
    const double sigma = std::sqrt(sum2 / n - mean * mean);
    EXPECT_NEAR(mean, 1, 0.006);
    EXPECT_NEAR(sigma, 0.15 / 2.35482, 0.006);

    // No resolution: no smearing
    TRestDetectorHitsSmearingProcess none("smearing", YAML::Load("{resolutionReference: 0, seed: 1}"));
    none.Initialize();
    TRestDetectorHitsEvent output;
    ASSERT_TRUE(none.ProcessEvent(input, output));
    EXPECT_NEAR(TotalEnergy(output), 5.9, 1e-5);
}

TEST(DetectorHitsSmearingProcess, SubtractFano) {
    auto detector = MakeDetector(kChargeNoDiffusion);
    // Resolution below the Fano contribution at this energy: not smeared
    TRestDetectorHitsSmearingProcess smearing(
        "smearing",
        YAML::Load("{energyReference: 5.9 keV, resolutionReference: 0.5, subtractFano: true, seed: 3}"));
    smearing.SetManager(detector.manager.get());
    smearing.Initialize();
    EXPECT_EQ(smearing.fMediumName, "ArIso");

    TRestDetectorHitsEvent input, output;
    input.AddHit(0, 0, 0, 5.9);
    ASSERT_TRUE(smearing.ProcessEvent(input, output));
    EXPECT_NEAR(TotalEnergy(output), 5.9, 1e-5);

    TRestDetectorHitsSmearingProcess noManager("smearing", YAML::Load("{subtractFano: true}"));
    EXPECT_THROW(noManager.Initialize(), std::runtime_error);
}

TEST(DetectorHitsTransformationProcess, TranslationRotationMirror) {
    TRestDetectorHitsEvent input;
    input.AddHit(10, 0, 5, 2, 1);
    input.AddHit(1, kNaN, 3, 4, 2, HitType::XZ);  // not XYZ: only the energy changes

    {
        TRestDetectorHitsTransformationProcess process(
            "transform", YAML::Load("{energyFactor: 0.5, rotationAxis: [0, 0, 1], rotationAngle: 90 deg, "
                                    "translation: [100 mm, 0, 1 cm]}"));
        process.Initialize();
        TRestDetectorHitsEvent output;
        ASSERT_TRUE(process.ProcessEvent(input, output));
        ASSERT_EQ(output.GetNumberOfHits(), 2u);
        // (10, 0, 5) rotated 90 deg around Z: (0, 10, 5), then translated
        EXPECT_NEAR(output.GetX(0), 100, 1e-4);
        EXPECT_NEAR(output.GetY(0), 10, 1e-4);
        EXPECT_NEAR(output.GetZ(0), 15, 1e-4);
        EXPECT_FLOAT_EQ(output.GetEnergy(0), 1);
        EXPECT_FLOAT_EQ(output.GetTime(0), 1);
        EXPECT_FLOAT_EQ(output.GetX(1), 1);
        EXPECT_FLOAT_EQ(output.GetZ(1), 3);
        EXPECT_FLOAT_EQ(output.GetEnergy(1), 2);
        EXPECT_EQ(output.GetType(1), HitType::XZ);
    }
    {
        // Rotation about an axis through a centre, and a mirror on the plane z = 1
        TRestDetectorHitsTransformationProcess process(
            "transform", YAML::Load("{rotationCenter: [10, 0, 0], rotationAxis: [0, 0, 2], rotationAngle: "
                                    "180 deg, mirrorPosition: [0, 0, 1], mirrorNormal: [0, 0, 1]}"));
        process.Initialize();
        TRestDetectorHitsEvent output;
        ASSERT_TRUE(process.ProcessEvent(input, output));
        // Mirror first: z = 5 -> -3; the point is on the rotation axis
        EXPECT_NEAR(output.GetX(0), 10, 1e-4);
        EXPECT_NEAR(output.GetY(0), 0, 1e-4);
        EXPECT_NEAR(output.GetZ(0), -3, 1e-4);
    }

    TRestDetectorHitsTransformationProcess bad("transform",
                                               YAML::Load("{rotationAxis: [0, 0, 0], rotationAngle: 1}"));
    EXPECT_THROW(bad.Initialize(), std::runtime_error);
}

TEST(DetectorHitsAnalysisProcess, ProjectionObservables) {
    AnalysisTree analysis;
    TRestDetectorHitsAnalysisProcess process("hitsAna",
                                             YAML::Load("{projections: [XZ, YZ], observables: [all]}"));
    process.SetRunInfo(&analysis.run);
    process.Initialize();

    TRestDetectorHitsEvent input, output;
    input.AddHit(1, kNaN, 10, 2, 0, HitType::XZ);
    input.AddHit(3, kNaN, 10, 2, 0, HitType::XZ);
    input.AddHit(kNaN, -2, 20, 6, 0, HitType::YZ);
    ASSERT_TRUE(process.ProcessEvent(input, output));
    EXPECT_EQ(output.GetNumberOfHits(), 3u);  // copied unchanged

    EXPECT_EQ(analysis.Get<int>("hitsAna_NHits"), 3);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_Energy"), 10, 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_MaxHitEnergy"), 6, 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_MinHitEnergy"), 2, 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_MeanHitEnergy"), 10 / 3., 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_MeanX"), 2, 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_MeanY"), -2, 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_MeanZ"), 16, 1e-5);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_SigmaX"), 1, 1e-5);
    EXPECT_EQ(analysis.Get<int>("hitsAna_NHits_XZ"), 2);
    EXPECT_EQ(analysis.Get<int>("hitsAna_NHits_YZ"), 1);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_Energy_YZ"), 6, 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_NHitsBalance_XZ_YZ"), 1 / 3., 1e-9);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_EnergyBalance_XZ_YZ"), -0.2, 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_SigmaBalance_XZ_YZ"), 1, 1e-5);
    // Not registered without their parameters
    EXPECT_THROW(analysis.Get<double>("hitsAna_FiducialEnergy"), std::runtime_error);

    TRestDetectorHitsEvent empty, emptyOutput;
    EXPECT_FALSE(process.ProcessEvent(empty, emptyOutput));
}

TEST(DetectorHitsAnalysisProcess, FiducialAndErrors) {
    AnalysisTree analysis;
    TRestDetectorHitsAnalysisProcess process(
        "hitsAna", YAML::Load("{fiducialRadius: 5 mm, fiducialLength: 10 mm, observables: [FiducialEnergy, "
                              "MeanX, SigmaY]}"));
    process.SetRunInfo(&analysis.run);
    process.Initialize();

    TRestDetectorHitsEvent input, output;
    input.AddHit(1, 0, 0, 3);   // inside
    input.AddHit(10, 0, 0, 2);  // outside the radius
    input.AddHit(0, 1, 8, 4);   // outside the length
    ASSERT_TRUE(process.ProcessEvent(input, output));
    EXPECT_NEAR(analysis.Get<double>("hitsAna_FiducialEnergy"), 3, 1e-6);
    EXPECT_NEAR(analysis.Get<double>("hitsAna_MeanX"), (3 + 20) / 9., 1e-5);
    // Not in the observables of the configuration
    EXPECT_THROW(analysis.Get<int>("hitsAna_NHits"), std::runtime_error);

    TRestDetectorHitsAnalysisProcess unknown("hitsAna", YAML::Load("{observables: [Unknown]}"));
    unknown.SetRunInfo(&analysis.run);
    EXPECT_THROW(unknown.Initialize(), std::runtime_error);
    EXPECT_THROW(TRestDetectorHitsAnalysisProcess("hitsAna", YAML::Load("{fiducialRadius: -1}")),
                 std::runtime_error);
}

// ============================================================================
// Signal processes
// ============================================================================

TEST(DetectorHitsToSignalProcess, ChannelsAndTimes) {
    auto detector = MakeDetector(kChargeNoDiffusion);
    AnalysisTree analysis;
    TRestDetectorHitsToSignalProcess process("hitsToSignal", YAML::Load("{observables: [NHitsRejected]}"));
    process.SetManager(detector.manager.get());
    process.SetRunInfo(&analysis.run);
    process.Initialize();

    // Hits on the readout plane (e.g. drifted by TRestDetectorHitsDriftProcess): the channel at their
    // position, at their time (no drift)
    const auto xPixel = XStripPixel(2, 1);
    const auto yPixel = YStripPixel(1, 2);
    TRestDetectorHitsEvent hits;
    hits.SetID(5);
    hits.AddHit(xPixel.X(), xPixel.Y(), 0, 3, 4);
    hits.AddHit(xPixel.X(), xPixel.Y(), 0, 1, 4);  // same channel and time: added
    hits.AddHit(yPixel.X(), yPixel.Y(), 0, 2, 1);
    hits.AddHit(xPixel.X(), xPixel.Y(), 30, 1, 0);  // in the gas, not on a channel
    hits.AddHit(100, 100, 0, 1, 0);                 // no channel

    TRestDetectorSignalEvent signals;
    ASSERT_TRUE(process.ProcessEvent(hits, signals));
    EXPECT_EQ(signals.GetID(), 5);
    EXPECT_EQ(analysis.Get<int>("hitsToSignal_NHitsRejected"), 2);
    ASSERT_EQ(signals.GetNumberOfSignals(), 2);
    ASSERT_TRUE(signals.HasSignalID(XStrip(2)));
    ASSERT_TRUE(signals.HasSignalID(YStrip(1)));
    const auto& xSignal = signals.GetSignalByID(XStrip(2));
    ASSERT_EQ(xSignal.GetNPoints(), 1);
    EXPECT_NEAR(xSignal.GetPoint(0).first, 4, 1e-6);
    EXPECT_NEAR(xSignal.GetPoint(0).second, 4, 1e-6);
    EXPECT_NEAR(signals.GetSignalByID(YStrip(1)).GetPoint(0).first, 1, 1e-6);

    // Sampling: times rounded down to its multiples
    TRestDetectorHitsToSignalProcess sampled("hitsToSignal", YAML::Load("{sampling: 0.5}"));
    sampled.SetManager(detector.manager.get());
    sampled.Initialize();
    TRestDetectorHitsEvent one;
    one.AddHit(xPixel.X(), xPixel.Y(), 0, 1, 0.75);
    TRestDetectorSignalEvent sampledSignals;
    ASSERT_TRUE(sampled.ProcessEvent(one, sampledSignals));
    EXPECT_NEAR(sampledSignals.GetSignal(0).GetPoint(0).first, 0.5, 1e-9);

    // Only rejected hits: no signals. The drift parameters are not keys of this process
    TRestDetectorHitsEvent rejected;
    rejected.AddHit(100, 100, 0, 1, 0);
    TRestDetectorSignalEvent none;
    EXPECT_FALSE(process.ProcessEvent(rejected, none));
    EXPECT_ANY_THROW(MetadataClassRegistry::Instance().Create("TRestDetectorHitsToSignalProcess", "h",
                                                              YAML::Load("{driftVelocity: 4}")));
}

TEST(DetectorSignalShapingProcess, GaussianResponse) {
    TRestDetectorSignalShapingProcess process("shaping",
                                              YAML::Load("{shapingType: gaus, shapingTime: 1 us, gain: 2}"));
    process.Initialize();

    TRestDetectorSignalEvent input, output;
    input.AddSignal(7, {{100, 3}, {10, 5}});  // far apart: two groups
    ASSERT_TRUE(process.ProcessEvent(input, output));
    ASSERT_EQ(output.GetNumberOfSignals(), 1);
    const auto& signal = output.GetSignal(0);
    EXPECT_EQ(signal.GetSignalID(), 7);

    // Height gain * charge at the time of each point, no grid points between the groups
    double max10 = 0, max100 = 0, timeMax10 = 0;
    for (int i = 0; i < signal.GetNPoints(); i++) {
        const auto [time, value] = signal.GetPoint(i);
        EXPECT_FALSE(time > 15 && time < 95) << "point at " << time;
        if (time < 50 && value > max10) {
            max10 = value;
            timeMax10 = time;
        }
        if (time > 50) max100 = std::max(max100, value);
    }
    EXPECT_NEAR(max10, 10, 1e-3);
    EXPECT_NEAR(max100, 6, 1e-3);
    EXPECT_NEAR(timeMax10, 10, 0.051);
    // Sigma of the response: value at one shaping time from the maximum
    for (int i = 0; i < signal.GetNPoints(); i++) {
        const auto [time, value] = signal.GetPoint(i);
        if (std::abs(time - 11) < 1e-6) {
            EXPECT_NEAR(value, 10 * std::exp(-0.5), 1e-3);
        }
    }

    TRestDetectorSignalShapingProcess bad("shaping", YAML::Load("{shapingTime: 0}"));
    EXPECT_THROW(bad.Initialize(), std::runtime_error);
}

TEST(DetectorSignalPeaksFinderProcess, PeaksAndWindows) {
    AnalysisTree analysis;
    TRestDetectorSignalPeaksFinderProcess process(
        "peaks",
        YAML::Load("{threshold: 100, distance: 2 us, halfWindow: 2, maxGreaterEqual: 0, "
                   "refineMaximum: false, amplitudeHalfWindow: 0, window: 1 us, observables: [all]}"));
    process.SetRunInfo(&analysis.run);
    process.Initialize();

    // Gaussian pulses (sigma 0.5 us, sampled every 0.1 us): signal 1 at 5 and 20 us, signal 2 at 5.3 us
    auto pulses = [](const std::vector<double>& centers) {
        std::vector<std::pair<double, double>> points;
        for (int i = 0; i <= 300; i++) {
            const double t = i * 0.1;
            double value = 0;
            for (double c : centers) value += 1000 * std::exp(-0.5 * (t - c) * (t - c) / 0.25);
            points.emplace_back(t, value);
        }
        return points;
    };
    TRestDetectorSignalEvent input, output;
    input.AddSignal(1, pulses({5, 20}));
    input.AddSignal(2, pulses({5.3}));
    ASSERT_TRUE(process.ProcessEvent(input, output));
    EXPECT_EQ(output.GetNumberOfSignals(), 2);  // copied unchanged

    EXPECT_EQ(analysis.Get<int>("peaks_NPeaks"), 3);
    EXPECT_EQ(analysis.Get<int>("peaks_NPeaksSignals"), 2);
    EXPECT_EQ(analysis.Get<int>("peaks_NWindows"), 2);
    EXPECT_NEAR(analysis.Get<double>("peaks_PeaksAmplitudeSum"), 3000, 1);

    const auto times = analysis.Get<std::vector<double>>("peaks_PeaksTime");
    ASSERT_EQ(times.size(), 3u);
    EXPECT_NEAR(times[0], 5, 1e-6);
    EXPECT_NEAR(times[1], 5.3, 1e-6);
    EXPECT_NEAR(times[2], 20, 1e-6);
    EXPECT_EQ(analysis.Get<std::vector<int>>("peaks_PeaksSignalID"), (std::vector<int>{1, 2, 1}));
    EXPECT_EQ(analysis.Get<std::vector<int>>("peaks_PeaksWindowIndex"), (std::vector<int>{0, 0, 1}));
    EXPECT_EQ(analysis.Get<std::vector<int>>("peaks_WindowMultiplicity"), (std::vector<int>{2, 1}));
}

/// A process initialized with the manager of the detector
template <typename Process>
std::unique_ptr<Process> Initialized(const Detector& detector, const std::string& yaml) {
    auto process = std::make_unique<Process>("process", YAML::Load(yaml));
    process->SetManager(detector.manager.get());
    process->Initialize();
    return process;
}

TEST(DetectorSignalToHitsProcess, IntegralAndQCenter) {
    auto detector = MakeDetector(kChargeNoDiffusion);
    auto process = Initialized<TRestDetectorSignalToHitsProcess>(detector, "{method: integral}");
    // The drift is reconstructed afterwards, with the drift velocity of the medium
    auto reconstruction = Initialized<TRestDetectorHitsTimeProjectionProcess>(detector, "{}");
    EXPECT_DOUBLE_EQ(reconstruction->fDriftVelocity, 10);
    EXPECT_EQ(reconstruction->fMediumName, "ArIso");

    TRestDetectorSignalEvent signals;
    signals.SetID(9);
    signals.AddSignal(XStrip(2), {{1, 10}, {2, 50}, {3, 20}});
    signals.AddSignal(YStrip(1), {{0.5, 4}});
    signals.AddSignal(999, {{1, 10}});  // not in the readout: skipped

    TRestDetectorHitsEvent onReadout, hits;
    ASSERT_TRUE(process->ProcessEvent(signals, onReadout));
    EXPECT_EQ(onReadout.GetID(), 9);
    ASSERT_EQ(onReadout.GetNumberOfHits(), 2u);
    // Only the coordinates of the channels: X of the X strip, Y of the Y strip
    EXPECT_EQ(onReadout.GetType(0), HitType::X);
    EXPECT_TRUE(std::isnan(onReadout.GetZ(0)));
    EXPECT_EQ(onReadout.GetType(1), HitType::Y);
    ASSERT_TRUE(reconstruction->ProcessEvent(onReadout, hits));
    EXPECT_EQ(hits.GetID(), 9);
    ASSERT_EQ(hits.GetNumberOfHits(), 2u);
    // X strip: X of the strip, Z from the drift (time of the maximum * velocity)
    EXPECT_EQ(hits.GetType(0), HitType::XZ);
    EXPECT_NEAR(hits.GetX(0), XStripPosition(2), 1e-5);
    EXPECT_TRUE(std::isnan(hits.GetY(0)));
    EXPECT_NEAR(hits.GetZ(0), 20, 1e-5);
    EXPECT_NEAR(hits.GetEnergy(0), 80, 1e-4);
    EXPECT_NEAR(hits.GetTime(0), 2, 1e-6);
    EXPECT_EQ(hits.GetType(1), HitType::YZ);
    EXPECT_NEAR(hits.GetY(1), YStripPosition(1), 1e-5);
    EXPECT_NEAR(hits.GetZ(1), 5, 1e-5);

    auto qCenter = Initialized<TRestDetectorSignalToHitsProcess>(detector, "{method: qCenter}");
    auto slower = Initialized<TRestDetectorHitsTimeProjectionProcess>(detector, "{driftVelocity: 2}");
    EXPECT_DOUBLE_EQ(slower->fDriftVelocity, 2);  // given: the medium is not used
    TRestDetectorHitsEvent qOnReadout, qHits;
    ASSERT_TRUE(qCenter->ProcessEvent(signals, qOnReadout));
    ASSERT_TRUE(slower->ProcessEvent(qOnReadout, qHits));
    ASSERT_GE(qHits.GetNumberOfHits(), 1u);
    EXPECT_NEAR(qHits.GetTime(0), (10 + 100 + 60) / 80., 1e-5);
    EXPECT_NEAR(qHits.GetEnergy(0), 80 / 3., 1e-4);
    EXPECT_NEAR(qHits.GetZ(0), 2 * (10 + 100 + 60) / 80., 1e-4);

    EXPECT_THROW(TRestDetectorSignalToHitsProcess("s", YAML::Load("{method: unknownMethod}")),
                 std::exception);
    // The drift parameters are not keys of SignalToHits
    EXPECT_ANY_THROW(MetadataClassRegistry::Instance().Create("TRestDetectorSignalToHitsProcess", "s",
                                                              YAML::Load("{driftVelocity: 2}")));

    // Normal along -Y: the drift fills Y (distance to the readout = velocity * time), the hits keep their Z
    auto inverted = MakeDetector(kChargeNoDiffusion);
    auto readout = BuildReadout("[0, -1, 0]");
    const auto* alongYReadout = readout.get();
    delete inverted.readouts->fReadoutMap["mmTest"];
    inverted.readouts->fReadoutMap["mmTest"] = readout.release();
    auto alongY = Initialized<TRestDetectorHitsTimeProjectionProcess>(inverted, "{}");
    TRestDetectorHitsEvent xHit, xyHit;
    xHit.AddHit(3, 0, 7, 1, 2, HitType::XZ);
    ASSERT_TRUE(alongY->ProcessEvent(xHit, xyHit));
    EXPECT_EQ(xyHit.GetType(0), HitType::XYZ);
    EXPECT_NEAR(xyHit.GetX(0), 3, 1e-5);
    EXPECT_NEAR(xyHit.GetZ(0), 7, 1e-5);
    EXPECT_NEAR(alongYReadout->GetDistanceToReadout(xyHit.GetPosition(0)), 20, 1e-4);

    // Without driftVelocity the medium must exist
    TRestDetectorHitsTimeProjectionProcess unknownMedium("r", YAML::Load("{medium: Other}"));
    unknownMedium.SetManager(detector.manager.get());
    EXPECT_ANY_THROW(unknownMedium.Initialize());
}

TEST(DetectorSignalToHitsProcess, RoundTripWithHitsToSignal) {
    // Drift (one electron, no fluctuations nor diffusion), signals, hits and drift reconstruction
    auto detector = MakeDetector(kChargeNoDiffusion);
    auto drift = Initialized<TRestDetectorHitsDriftProcess>(detector, "{fluctuations: none, maxHits: 1}");
    auto toSignal = Initialized<TRestDetectorHitsToSignalProcess>(detector, "{}");
    auto toHits = Initialized<TRestDetectorSignalToHitsProcess>(detector, "{method: integral}");
    auto reconstruction = Initialized<TRestDetectorHitsTimeProjectionProcess>(detector, "{}");

    const auto pixel = YStripPixel(3, 0);
    TRestDetectorHitsEvent input, drifted, onReadout, output;
    input.AddHit(pixel.X(), pixel.Y(), 27, 5.9, 0);
    TRestDetectorSignalEvent signals;
    ASSERT_TRUE(drift->ProcessEvent(input, drifted));
    ASSERT_TRUE(toSignal->ProcessEvent(drifted, signals));
    ASSERT_TRUE(toHits->ProcessEvent(signals, onReadout));
    ASSERT_TRUE(reconstruction->ProcessEvent(onReadout, output));
    ASSERT_EQ(output.GetNumberOfHits(), 1u);
    EXPECT_NEAR(output.GetY(0), YStripPosition(3), 1e-5);
    EXPECT_NEAR(output.GetZ(0), 27, 1e-4);
    EXPECT_NEAR(output.GetEnergy(0), 5.9, 1e-5);
}

TEST(DetectorSignalRecoveryProcess, DeadChannelFromNeighbours) {
    auto detector = MakeDetector(kChargeNoDiffusion);
    TRestDetectorSignalRecoveryProcess process(
        "recovery", YAML::Load("{channelIds: [" + std::to_string(YStrip(1)) + "]}"));
    process.SetManager(detector.manager.get());
    process.Initialize();

    TRestDetectorSignalEvent input, output;
    input.AddSignal(YStrip(0), {{1, 100}, {2, 40}});
    input.AddSignal(YStrip(2), {{1, 60}});
    ASSERT_TRUE(process.ProcessEvent(input, output));

    // Half of the charge of each neighbour goes to the dead channel
    ASSERT_TRUE(output.HasSignalID(YStrip(1)));
    const auto& recovered = output.GetSignalByID(YStrip(1));
    ASSERT_EQ(recovered.GetNPoints(), 2);
    EXPECT_NEAR(recovered.GetPoint(0).second, 0.5 * 100 + 0.5 * 60, 1e-9);
    EXPECT_NEAR(recovered.GetPoint(1).second, 0.5 * 40, 1e-9);
    EXPECT_NEAR(output.GetSignalByID(YStrip(0)).GetPoint(0).second, 50, 1e-9);
    EXPECT_NEAR(output.GetSignalByID(YStrip(2)).GetPoint(0).second, 30, 1e-9);

    // The total charge is preserved
    double before = 0, after = 0;
    for (int s = 0; s < input.GetNumberOfSignals(); s++)
        for (double q : input.GetSignal(s).GetCharge()) before += q;
    for (int s = 0; s < output.GetNumberOfSignals(); s++)
        for (double q : output.GetSignal(s).GetCharge()) after += q;
    EXPECT_NEAR(after, before, 1e-9);
}

// A copied event (copy constructor and assignment) has its own storage: changing the original after the
// copy does not change the copy (its views do not point to the storage of the original)
TEST(DetectorHitsEvent, CopyHasItsOwnStorage) {
    TRestDetectorHitsEvent original;
    original.SetID(7);
    original.AddHit(1, 2, 3, 10);
    original.AddHit(4, 5, 6, 20);
    TRestDetectorHitsEvent copy(original);
    TRestDetectorHitsEvent assigned;
    assigned = original;
    original.Initialize();
    original.AddHit(100, 200, 300, 1000);
    for (const auto* event : {&copy, &assigned}) {
        EXPECT_EQ(event->GetID(), 7);
        ASSERT_EQ(event->GetNumberOfHits(), 2u);
        EXPECT_DOUBLE_EQ(event->GetX(0), 1);
        EXPECT_DOUBLE_EQ(event->GetX(1), 4);
        EXPECT_DOUBLE_EQ(event->GetEnergy(1), 20);
    }
}

TEST(DetectorSignalEvent, CopyHasItsOwnStorage) {
    TRestDetectorSignalEvent original;
    original.AddSignal(3, {{0., 1.}, {1., 2.}});
    original.GetSignal(0);  // the views are built before the copy
    TRestDetectorSignalEvent copy(original);
    TRestDetectorSignalEvent assigned;
    assigned = original;
    original.Initialize();
    original.AddSignal(9, {{0., 5.}});
    original.GetSignal(0);
    for (auto* event : {&copy, &assigned}) {
        ASSERT_EQ(event->GetNumberOfSignals(), 1);
        EXPECT_EQ(event->GetSignal(0).GetSignalID(), 3);
    }
}

TEST(DetectorTransportTable, SeveralMediaInAFile) {
    // Two media of a mixture (Ar99_iC4H10_1, Ar98_iC4H10_2) with tables of the same key in the same file
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/detector_test_media.root";
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        for (const auto& [name, velocity] : {std::pair<std::string, double>{"Ar99_iC4H10_1", 10.},
                                             std::pair<std::string, double>{"Ar98_iC4H10_2", 20.}}) {
            TRestDetectorTransportTable table;
            table.fMedium = name;
            table.fFields = {1., 10.};
            auto& electron = table.fCarriers["electron"];
            electron.fDriftVelocity = {velocity, velocity};
            electron.fLongitudinalDiffusion = electron.fTransversalDiffusion = {0.03, 0.03};
            electron.fTownsend = electron.fAttachment = {0, 0};
            table.Write(file.get(), YAML::Load("{name: " + name + ", state: gas}"), YAML::Load("{}"));
        }
    }
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    for (const auto& [name, velocity] : {std::pair<std::string, double>{"Ar99_iC4H10_1", 10.},
                                         std::pair<std::string, double>{"Ar98_iC4H10_2", 20.}}) {
        const YAML::Node tables = TRestDetectorTransportTable::GetTables(file.get(), name);
        ASSERT_EQ(tables.size(), 1u);
        const auto table =
            TRestDetectorTransportTable::Read(file.get(), name, tables[0]["key"].as<std::string>());
        EXPECT_DOUBLE_EQ(table.GetDriftVelocity(5., 1.), velocity);
    }
}

TEST(DetectorAvalancheProcess, PolyaGainFromTheTable) {
    // Gain table of a 128 um gap at 1 bar and 293.15 K: gain 100 -> 1000 and theta 0.5 -> 1.5 from 2000 to
    // 4000 V/mm
    const std::string fileName = std::string(TESTS_WORK_DIR) + "/detector_test_gain.root";
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        TRestDetectorGainTable table;
        table.fMedium = "ArIsoGain";
        table.fTableGap = 0.128;
        table.fFields = {2000., 4000.};
        table.fGain = {100., 1000.};
        table.fTheta = {0.5, 1.5};
        table.Write(file.get(), YAML::Load("{name: ArIsoGain, state: gas}"), YAML::Load("{}"), {});
    }
    const auto manager =
        std::make_unique<TRestManager>("manager", YAML::Load("class: TRestManager\n"
                                                             "gas:\n"
                                                             "  class: TRestDetectorMedium\n"
                                                             "  name: ArIsoGain\n"
                                                             "  mediumFile: \"" +
                                                             fileName +
                                                             "\"\n"
                                                             "  charge: {wValue: 26.2 eV}\n"));
    auto makeAvalanche = [&manager](const std::string& config) {
        auto avalanche = std::make_unique<TRestDetectorAvalancheProcess>("avalanche", YAML::Load(config));
        avalanche->SetManager(manager.get());
        return avalanche;
    };

    // 384 V / 128 um = 3000 V/mm: gain and theta interpolated (log G and theta linear in the field)
    auto avalanche = makeAvalanche("{amplificationVoltage: 384 V, gap: 128 um, seed: 5}");
    avalanche->Initialize();
    EXPECT_EQ(avalanche->fMediumName, "ArIsoGain");
    EXPECT_EQ(avalanche->fTable, "1bar_293.15K_gap128um_E20000-40000Vcm");
    EXPECT_NEAR(avalanche->fGain, std::sqrt(100. * 1000.), 1e-6);
    EXPECT_NEAR(avalanche->fTheta, 1, 1e-9);

    // One electron per hit: Polya gains, mean G and relative variance 1 / (theta + 1)
    const double wValue = 0.0262, gain = avalanche->fGain;
    const int n = 20000;
    double sum = 0, sum2 = 0;
    for (int i = 0; i < n; i++) {
        TRestDetectorHitsEvent input, output;
        input.SetID(i);
        input.AddHit(1, 2, 3, wValue, 0.5);
        ASSERT_TRUE(avalanche->ProcessEvent(input, output));
        ASSERT_EQ(output.GetNumberOfHits(), 1u);
        const double ratio = output.GetEnergy(0) / wValue / gain;
        sum += ratio;
        sum2 += ratio * ratio;
    }
    const double mean = sum / n;
    EXPECT_NEAR(mean, 1, 0.02);
    EXPECT_NEAR(sum2 / n - mean * mean, 0.5, 0.03);

    // A hit of 1000 electrons: the gain per electron of the sum, same position and time; seeded per event
    TRestDetectorHitsEvent input, first, second;
    input.SetID(3);
    input.AddHit(1, 2, 3, 1000 * wValue, 0.5);
    input.AddHit(0, 0, 0, 0);  // no electrons: removed
    ASSERT_TRUE(avalanche->ProcessEvent(input, first));
    ASSERT_TRUE(avalanche->ProcessEvent(input, second));
    ASSERT_EQ(first.GetNumberOfHits(), 1u);
    EXPECT_DOUBLE_EQ(first.GetEnergy(0), second.GetEnergy(0));
    EXPECT_NEAR(first.GetEnergy(0) / (1000 * wValue * gain), 1, 5 / std::sqrt(2000.));
    EXPECT_NEAR(first.GetX(0), 1, 1e-6);
    EXPECT_NEAR(first.GetTime(0), 0.5, 1e-9);

    // Errors: field out of the range of the table, no table with the gap, no voltage, no manager
    EXPECT_THROW(makeAvalanche("{amplificationVoltage: 600 V, gap: 128 um}")->Initialize(),
                 std::runtime_error);
    EXPECT_THROW(makeAvalanche("{amplificationVoltage: 300 V, gap: 100 um}")->Initialize(),
                 std::runtime_error);
    EXPECT_THROW(makeAvalanche("{gap: 128 um}")->Initialize(), std::runtime_error);
    TRestDetectorAvalancheProcess noManager("avalanche",
                                            YAML::Load("{amplificationVoltage: 384 V, gap: 128 um}"));
    EXPECT_THROW(noManager.Initialize(), std::runtime_error);
}

TEST(DetectorAvalancheProcess, GainFromTheAmplificationResponse) {
    // 4 primaries, 3 avalanches at the anode: sizes 100, 100 and 200 at 10, 10 and 30 ns
    const std::string config =
        "{fieldMaps: [cell.txt], responseFile: r.root, start: 1 mm, anode: -55 um, sizeBins: 10, "
        "timeBins: 4}";
    TRestDetectorAmplificationResponse response("amplification", YAML::Load(config));
    EXPECT_DOUBLE_EQ(response.fStart, 1);
    EXPECT_DOUBLE_EQ(response.fAnode, -0.055);
    EXPECT_EQ(response.fPeriodicity, "xy");
    response.fMedium = "ArIso";
    TRestDetectorAmplificationResponse::Losses losses;
    losses.fConductor = 1;  // the fourth primary, e.g. on the mesh
    // Without the cell (fCellLow = fCellHigh): no map
    response.SetResult({{0, 0, 100, 0.01, 0}, {0, 0, 100, 0.01, 0}, {0, 0, 200, 0.03, 0}, {0, 0, 0, 0, 0}},
                       losses);
    EXPECT_FALSE(response.HasCellMap());
    EXPECT_DOUBLE_EQ(response.fAnodeElectrons, 400);
    EXPECT_DOUBLE_EQ(response.GetTransparency(), 0.75);
    EXPECT_NEAR(response.fMeanSize, 400. / 3, 1e-9);
    EXPECT_NEAR(response.fSizeVariance, 10000. / 3, 1e-6);
    EXPECT_NEAR(response.GetTheta(), (400. / 3) * (400. / 3) / (10000. / 3) - 1, 1e-9);
    // Bins of 20 over [0, 200]: 100 in [100, 120), 200 in the last one; times in [10, 30] ns, 4 bins
    EXPECT_EQ(response.fSizeCounts[5], 2);
    EXPECT_EQ(response.fSizeCounts[9], 1);
    EXPECT_EQ(response.fTimeCounts[0], 2);
    EXPECT_EQ(response.fTimeCounts[3], 1);

    const std::string fileName = std::string(TESTS_WORK_DIR) + "/AmplificationResponse.root";
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        response.Write(file.get(), YAML::Load("{program: test}"));
        EXPECT_NE(file->Get("AmplificationSize"), nullptr);
    }
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
        const auto read = TRestDetectorAmplificationResponse::Read(file.get());
        EXPECT_EQ(read.fMedium, "ArIso");
        EXPECT_EQ(read.fPrimaries, 4);
        EXPECT_EQ(read.fCollected, 3);
        EXPECT_EQ(read.fLosses.fConductor, 1);
        EXPECT_EQ(read.fLosses.fDielectric, 0);
        EXPECT_EQ(read.fSizeCounts, response.fSizeCounts);
        EXPECT_EQ(read.fTimeCounts, response.fTimeCounts);
        EXPECT_DOUBLE_EQ(read.fSizeMax, 200);
        TRandom3 random(2);
        for (int i = 0; i < 100; i++) {
            const double size = read.SampleSize(random);
            EXPECT_TRUE((size >= 100 && size < 120) || (size >= 180 && size <= 200)) << size;
            const double time = read.SampleTime(random);
            EXPECT_TRUE((time >= 0.01 && time < 0.015) || (time >= 0.025 && time <= 0.03)) << time;
        }
    }

    // The process: one electron per hit, amplified with the transparency 0.75, its gain and time sampled
    auto detector = MakeDetector(kChargeNoDiffusion);
    TRestDetectorAvalancheProcess avalanche("avalanche",
                                            YAML::Load("{responseFile: \"" + fileName + "\", seed: 7}"));
    avalanche.SetManager(detector.manager.get());
    avalanche.Initialize();
    EXPECT_NEAR(avalanche.fGain, 400. / 3, 1e-9);
    EXPECT_DOUBLE_EQ(avalanche.fTransparency, 0.75);
    const double wValue = 0.0262;
    const int n = 4000;
    int amplified = 0;
    double sum = 0;
    for (int i = 0; i < n; i++) {
        TRestDetectorHitsEvent input, output;
        input.SetID(i);
        input.AddHit(1, 2, 0, wValue, 0.5);
        if (!avalanche.ProcessEvent(input, output)) continue;
        ASSERT_EQ(output.GetNumberOfHits(), 1u);
        amplified++;
        sum += output.GetEnergy(0) / wValue;
        EXPECT_NEAR(output.GetX(0), 1, 1e-6);
        EXPECT_GE(output.GetTime(0), 0.51 - 1e-9);
        EXPECT_LE(output.GetTime(0), 0.53 + 1e-9);
    }
    EXPECT_NEAR(static_cast<double>(amplified) / n, 0.75, 0.03);
    // The bin centres of the sampled sizes: (2 x 110 + 190) / 3
    EXPECT_NEAR(sum / amplified, 410. / 3, 3);

    // Errors: responseFile and a gap, no response file, a response without avalanches, the configuration
    auto initialize = [&detector](const std::string& node) {
        TRestDetectorAvalancheProcess process("avalanche", YAML::Load(node));
        process.SetManager(detector.manager.get());
        process.Initialize();
    };
    EXPECT_THROW(initialize("{responseFile: \"" + fileName + "\", gap: 128 um}"), std::runtime_error);
    EXPECT_ANY_THROW(initialize("{responseFile: none.root}"));
    const std::string emptyFile = std::string(TESTS_WORK_DIR) + "/AmplificationResponse_Empty.root";
    {
        response.SetResult(std::vector<TRestDetectorAmplificationResponse::Avalanche>(4));
        std::unique_ptr<TFile> file(TFile::Open(emptyFile.c_str(), "RECREATE"));
        response.Write(file.get(), YAML::Load("{program: test}"));
    }
    EXPECT_THROW(initialize("{responseFile: \"" + emptyFile + "\"}"), std::runtime_error);
    const std::string base = "fieldMaps: [cell.txt], responseFile: r.root";
    EXPECT_THROW(
        TRestDetectorAmplificationResponse("a", YAML::Load("{" + base + ", start: 0 mm, anode: 0 mm}")),
        std::runtime_error);
    EXPECT_THROW(TRestDetectorAmplificationResponse(
                     "a", YAML::Load("{" + base + ", start: 1 mm, dielectricMaps: [a.txt, b.txt]}")),
                 std::runtime_error);
    EXPECT_THROW(
        TRestDetectorAmplificationResponse("a", YAML::Load("{" + base + ", start: 1 mm, periodicity: xw}")),
        std::runtime_error);
    EXPECT_THROW(TRestDetectorAmplificationResponse("a", YAML::Load("{responseFile: r.root, start: 1 mm}")),
                 std::runtime_error);
    EXPECT_THROW(
        TRestDetectorAmplificationResponse("a", YAML::Load("{" + base + ", start: 1 mm, cellBins: -1}")),
        std::runtime_error);
}

TEST(DetectorAvalancheProcess, AmplificationDependingOnThePositionInThePixel) {
    // A cell of the size of the pixels of the test readout (pitch / sqrt(2)), 2 x 2 bins: the primaries
    // started in the half x < 0 are lost (e.g. on the solid mesh), those in x > 0 give 100 electrons in the
    // bin y < 0 and 300 in y > 0 (mean 200); 1 electron of 400 on the pad of another cell
    const double half = kPitch / std::sqrt(8.);
    TRestDetectorAmplificationResponse response(
        "amplification",
        YAML::Load("{fieldMaps: [cell.txt], responseFile: r.root, start: 1 mm, anode: -55 um, "
                   "sizeBins: 30, timeBins: 1, cellBins: 2}"));
    response.fMedium = "ArIso";
    response.fCellLow = {-half, -half};
    response.fCellHigh = {half, half};
    const double q = half / 2;
    response.SetResult({{-q, -q, 0, 0, 0}, {-q, q, 0, 0, 0}, {q, -q, 100, 0.01, 1}, {q, q, 300, 0.01, 0}});
    ASSERT_TRUE(response.HasCellMap());
    EXPECT_DOUBLE_EQ(response.GetTransparency(), 0.5);
    const int lostBin = response.GetCellBin(0.25, 0.75), lowBin = response.GetCellBin(0.75, 0.25),
              highBin = response.GetCellBin(0.75, 0.75);
    EXPECT_EQ(lostBin, 2);
    EXPECT_EQ(lowBin, 1);
    EXPECT_EQ(highBin, 3);
    EXPECT_DOUBLE_EQ(response.GetCellTransparency(lostBin), 0);
    EXPECT_DOUBLE_EQ(response.GetCellTransparency(highBin), 1);
    EXPECT_DOUBLE_EQ(response.GetCellGainFactor(lowBin), 0.5);
    EXPECT_DOUBLE_EQ(response.GetCellGainFactor(highBin), 1.5);
    EXPECT_DOUBLE_EQ(response.GetCellGainFactor(lostBin), 1);  // no avalanche: the mean
    EXPECT_DOUBLE_EQ(response.fSharedElectrons, 1);
    EXPECT_DOUBLE_EQ(response.fAnodeElectrons, 400);

    const std::string fileName = std::string(TESTS_WORK_DIR) + "/AmplificationResponse_Cell.root";
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        response.Write(file.get(), YAML::Load("{program: test}"));
        EXPECT_NE(file->Get("AmplificationGainMap"), nullptr);
    }
    {
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
        const auto read = TRestDetectorAmplificationResponse::Read(file.get());
        EXPECT_EQ(read.fCellBins, 2);
        EXPECT_EQ(read.fCellPrimaries, response.fCellPrimaries);
        EXPECT_EQ(read.fCellCollected, response.fCellCollected);
        EXPECT_EQ(read.fCellSizeSum, response.fCellSizeSum);
        EXPECT_DOUBLE_EQ(read.fCellHigh[0], half);
    }

    // The process: hits 0.1 mm from the centre of a pixel along its local x and 0.05 mm along its local y,
    // on each side of x
    auto detector = MakeDetector(kChargeNoDiffusion);
    TRestDetectorAvalancheProcess avalanche("avalanche",
                                            YAML::Load("{responseFile: \"" + fileName + "\", seed: 3}"));
    avalanche.SetManager(detector.manager.get());
    avalanche.Initialize();
    EXPECT_EQ(avalanche.fReadoutName, "mmTest");
    const auto centre = YStripPixel(1, 2);
    const XYZVector offset = XYZVector(0.1, 0.1, 0) / std::sqrt(2.);
    const XYZVector across = XYZVector(-0.05, 0.05, 0) / std::sqrt(2.);
    TRestDetectorReadout::PixelPosition pixel;
    ASSERT_TRUE(detector.readout->GetPixelPosition(centre + offset + across, pixel));
    ASSERT_NEAR(std::abs(pixel.local.X()), 0.1, 1e-9);
    ASSERT_NEAR(std::abs(pixel.local.Y()), 0.05, 1e-9);
    // The side of the pixel with local x > 0 (amplified) and the other one (lost)
    const bool positive = pixel.local.X() > 0;
    const XYZVector amplifiedSide = centre + across + (positive ? offset : -offset);
    const XYZVector lostSide = centre + across - (positive ? offset : -offset);
    ASSERT_TRUE(detector.readout->GetPixelPosition(amplifiedSide, pixel));
    const int amplifiedBin =
        response.GetCellBin(0.5 + pixel.local.X() / (2 * half), 0.5 + pixel.local.Y() / (2 * half));
    EXPECT_TRUE(amplifiedBin == lowBin || amplifiedBin == highBin);
    const double wValue = 0.0262;
    int amplified = 0, lost = 0;
    double sum = 0;
    const int n = 2000;
    for (int i = 0; i < n; i++) {
        TRestDetectorHitsEvent input, output;
        input.SetID(i);
        input.AddHit(amplifiedSide.X(), amplifiedSide.Y(), 0, wValue, 0.5);
        input.AddHit(lostSide.X(), lostSide.Y(), 0, wValue, 0.5);
        if (!avalanche.ProcessEvent(input, output)) continue;
        ASSERT_EQ(output.GetNumberOfHits(), 1u);
        // The positions of the hits are stored as float; the two sides are 0.14 mm apart along x
        if (std::abs(output.GetX(0) - amplifiedSide.X()) < 1e-4) {
            amplified++;
            sum += output.GetEnergy(0) / wValue;
        } else {
            lost++;
        }
    }
    EXPECT_EQ(amplified, n);
    EXPECT_EQ(lost, 0);
    // The sizes of the whole cell (105 or 295, bins of 10 over [0, 300]: mean 200) scaled by the mean of the
    // bin of the hit over the mean of all
    EXPECT_NEAR(sum / amplified, 200 * response.GetCellGainFactor(amplifiedBin), 15);

    // Out of the pixels: the response averaged over the cell (transparency 0.5)
    int outside = 0;
    for (int i = 0; i < n; i++) {
        TRestDetectorHitsEvent input, output;
        input.SetID(n + i);
        input.AddHit(100, 100, 0, wValue, 0.5);
        if (avalanche.ProcessEvent(input, output)) outside++;
    }
    EXPECT_NEAR(static_cast<double>(outside) / n, 0.5, 0.05);

    // A response with the map needs a readout
    auto noReadout = std::make_unique<TRestManager>(
        "manager", YAML::Load("class: TRestManager\ngas:\n  class: TRestDetectorMedium\n  name: ArIso\n"
                              "  charge: " +
                              kChargeNoDiffusion + "\n"));
    TRestDetectorAvalancheProcess withoutReadout("avalanche",
                                                 YAML::Load("{responseFile: \"" + fileName + "\"}"));
    withoutReadout.SetManager(noReadout.get());
    EXPECT_THROW(withoutReadout.Initialize(), std::runtime_error);
}

// ---------------------------------------------------------------------------
// TRestDataSetTardis
// ---------------------------------------------------------------------------

namespace {

/// An entry of TRestDataSetTardis::TimeDelays at seconds + nanoseconds with physical channels
TRestDataSetTardis::Entry TardisEntry(long long seconds, int nanoseconds, std::vector<int> channels = {}) {
    TRestDataSetTardis::Entry entry;
    entry.seconds = seconds;
    entry.nanoseconds = nanoseconds;
    entry.channels = std::move(channels);
    return entry;
}

/// Physical IDs of the readout of the tests: X strip n and Y strip n
int XStripID(int n) { return kNChannels + n; }
int YStripID(int n) { return n; }

/// A file of the current test in the work directory (ctest runs the tests in parallel)
std::string TardisPath(const std::string& name) {
    const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
    return std::string(TESTS_WORK_DIR) + "/tardis_" + test->name() + "_" + name;
}

/// A dataset of two runs with the timestamps out of order: eventID = entry, runIndex = entry % 2.
/// Chronological order: 4 (0 s), 3 (500 ns), 1 (999999999 ns), 0 and 2 (1 s, same timestamp).
std::string WriteTardisDataSet(const std::string& name) {
    const std::string path = TardisPath(name);
    TFile file(path.c_str(), "RECREATE");
    TTree tree("AnalysisTree", "AnalysisTree");
    int runIndex = 0, eventID = 0, subEventID = 0, nanos = 0;
    long long seconds = 0;
    tree.Branch("runIndex", &runIndex);
    tree.Branch("eventID", &eventID);
    tree.Branch("subEventID", &subEventID);
    tree.Branch("timeSeconds", &seconds);
    tree.Branch("timeNanoSecs", &nanos);
    const std::vector<long long> times = {1700000001, 1700000000, 1700000001, 1700000000, 1700000000};
    const std::vector<int> fractions = {0, 999999999, 0, 500, 0};
    for (size_t i = 0; i < times.size(); ++i) {
        runIndex = i % 2;
        eventID = i;
        seconds = times[i];
        nanos = fractions[i];
        tree.Fill();
    }
    tree.Write();
    TObjString runs(
        "- run: {runNumber: 1}\n  fileName: a.root\n  entries: 3\n  processing: 0\n"
        "- run: {runNumber: 2}\n  fileName: b.root\n  entries: 2\n  processing: 0\n");
    file.WriteTObject(&runs, "RESTDataSetInfo");
    return path;
}

/// eventID -> <column> of a dataset file
std::map<int, double> TardisDelays(const std::string& fileName, const std::string& column) {
    TFile file(fileName.c_str(), "READ");
    auto* tree = file.Get<TTree>("AnalysisTree");
    if (!tree) throw std::runtime_error("no AnalysisTree in " + fileName);
    int id = 0;
    double delay = 0;
    tree->SetBranchAddress("eventID", &id);
    tree->SetBranchAddress(column.c_str(), &delay);
    std::map<int, double> delays;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        delays[id] = delay;
    }
    return delays;
}

}  // namespace

TEST(DataSetTardis, TimeDelaysWithoutReadout) {
    const auto delays = TRestDataSetTardis::TimeDelays(
        {TardisEntry(10, 999999999), TardisEntry(11, 0), TardisEntry(11, 0), TardisEntry(13, 500)});
    ASSERT_EQ(delays.size(), 4u);
    EXPECT_TRUE(std::isnan(delays[0]));
    EXPECT_NEAR(delays[1], 1e-9, 1e-15);
    EXPECT_DOUBLE_EQ(delays[2], 0);
    EXPECT_NEAR(delays[3], 2 + 500e-9, 1e-12);
}

TEST(DataSetTardis, TimeDelaysSinceEventsOverAThreshold) {
    // Only the entries with canBePrevious (previousObservable over previousThreshold) are previous events
    auto entry = [](long long seconds, bool canBePrevious, std::vector<int> channels = {}) {
        auto result = TardisEntry(seconds, 0, std::move(channels));
        result.canBePrevious = canBePrevious;
        return result;
    };
    const auto delays = TRestDataSetTardis::TimeDelays(
        {entry(0, true), entry(1, false), entry(2, false), entry(5, true), entry(6, false)});
    ASSERT_EQ(delays.size(), 5u);
    EXPECT_TRUE(std::isnan(delays[0]));
    EXPECT_DOUBLE_EQ(delays[1], 1);
    EXPECT_DOUBLE_EQ(delays[2], 2);
    EXPECT_DOUBLE_EQ(delays[3], 5);
    EXPECT_DOUBLE_EQ(delays[4], 1);

    // With the readout: the previous event in the same place and over the threshold
    const auto readout = BuildReadout();
    const auto inPlace = TRestDataSetTardis::TimeDelays(
        {entry(0, true, {XStripID(0), YStripID(0)}), entry(1, false, {XStripID(0), YStripID(0)}),
         entry(2, false, {XStripID(1), YStripID(1)}), entry(3, true, {XStripID(0), YStripID(0)})},
        readout.get(), 1);
    EXPECT_DOUBLE_EQ(inPlace[1], 1);
    EXPECT_DOUBLE_EQ(inPlace[2], 2);  // the entry 0, not the entry 1 (under the threshold)
    EXPECT_DOUBLE_EQ(inPlace[3], 3);
}

TEST(DataSetTardis, TimeDelaysInTheSameStrips) {
    const auto readout = BuildReadout();
    const std::vector<TRestDataSetTardis::Entry> entries = {
        TardisEntry(0, 0, {XStripID(0), YStripID(0)}),
        TardisEntry(1, 0, {XStripID(3), YStripID(3)}),
        TardisEntry(2, 0, {XStripID(1), YStripID(1)}),  // next to the entry 0 in X and in Y
        TardisEntry(3, 0, {XStripID(1), YStripID(3)}),  // X of 0 and 2, Y of 1: none in both
        TardisEntry(4, 0, {XStripID(2)}),               // only X: next to the entry 3
        TardisEntry(5, 500000000, {YStripID(3)}),       // only Y: the entry 3 (not the X strip 0 next to it)
        TardisEntry(6, 0, {-1, 99}),                    // no channel of the readout
    };

    const auto delays = TRestDataSetTardis::TimeDelays(entries, readout.get(), 1);
    ASSERT_EQ(delays.size(), entries.size());
    EXPECT_TRUE(std::isnan(delays[0]));
    EXPECT_TRUE(std::isnan(delays[1]));
    EXPECT_DOUBLE_EQ(delays[2], 2);
    EXPECT_TRUE(std::isnan(delays[3]));
    EXPECT_DOUBLE_EQ(delays[4], 1);
    EXPECT_DOUBLE_EQ(delays[5], 2.5);
    EXPECT_TRUE(std::isnan(delays[6]));

    // Only the same channels
    const auto same = TRestDataSetTardis::TimeDelays(entries, readout.get(), 0);
    EXPECT_TRUE(std::isnan(same[2]));
    EXPECT_TRUE(std::isnan(same[4]));
    EXPECT_DOUBLE_EQ(same[5], 2.5);

    // Two neighbours: the entry 3 is next to the entry 2 (the same X strip, Y strips 3 and 1)
    const auto two = TRestDataSetTardis::TimeDelays(entries, readout.get(), 2);
    EXPECT_DOUBLE_EQ(two[3], 1);
}

TEST(DataSetTardis, ApplyAcrossRunsWithoutReordering) {
    const std::string input = WriteTardisDataSet("input.root");
    TRestDataSetTardis tardis("tardis", YAML::Load("{class: TRestDataSetTardis, verbose: warning}"));
    for (const bool multithreaded : {false, true}) {
        if (multithreaded) ROOT::EnableImplicitMT(2);
        TRestDataSet dataSet(input);
        tardis.ObtainParameters(dataSet, nullptr);
        tardis.Apply(dataSet);
        const std::string output = TardisPath(multithreaded ? "output_mt.root" : "output.root");
        dataSet.Write(output);
        if (multithreaded) ROOT::DisableImplicitMT();

        const auto delays = TardisDelays(output, "tardis_timeDelay");
        ASSERT_EQ(delays.size(), 5u);
        EXPECT_TRUE(std::isnan(delays.at(4)));
        EXPECT_NEAR(delays.at(3), 500e-9, 1e-15);
        EXPECT_NEAR(delays.at(1), 0.999999999 - 500e-9, 1e-15);
        // Same timestamp (both of run index 0): ordered by event, 0 first
        EXPECT_NEAR(delays.at(0), 1e-9, 1e-15);
        EXPECT_DOUBLE_EQ(delays.at(2), 0);

        // The dataset is not reordered (one thread)
        if (multithreaded) continue;
        TFile file(output.c_str(), "READ");
        auto* tree = file.Get<TTree>("AnalysisTree");
        ASSERT_NE(tree, nullptr);
        int id = 0;
        tree->SetBranchAddress("eventID", &id);
        for (int i = 0; i < 5; ++i) {
            tree->GetEntry(i);
            EXPECT_EQ(id, i);
        }
    }
}

TEST(DataSetTardis, CutsFilterTheDataSet) {
    const std::string input = WriteTardisDataSet("input.root");
    TRestManager manager("manager", YAML::Load("class: TRestManager\n"
                                               "even:\n"
                                               "  class: TRestCut\n"
                                               "  conditions: {even: \"runIndex == 0\"}\n"));
    TRestDataSetTardis tardis("tardis",
                              YAML::Load("{class: TRestDataSetTardis, verbose: warning, cuts: [even]}"));
    TRestDataSet dataSet(input);
    tardis.ObtainParameters(dataSet, &manager);
    tardis.Apply(dataSet);
    const std::string output = TardisPath("output.root");
    dataSet.Write(output);

    // Entries 0, 2 and 4 (run index 0): 4 at 0 s, then 0 and 2 at 1 s
    const auto delays = TardisDelays(output, "tardis_timeDelay");
    ASSERT_EQ(delays.size(), 3u);
    EXPECT_TRUE(std::isnan(delays.at(4)));
    EXPECT_DOUBLE_EQ(delays.at(0), 1);
    EXPECT_DOUBLE_EQ(delays.at(2), 0);
}

TEST(DataSetTardis, Errors) {
    // A reference is not used: cleared (warning)
    TRestDataSetTardis tardis("tardis", YAML::Load("{class: TRestDataSetTardis, verbose: error, "
                                                   "referenceFileName: other.root}"));
    EXPECT_TRUE(tardis.fReferenceFileName.empty());
    const YAML::Node negative = YAML::Load("{class: TRestDataSetTardis, channelTolerance: -1}");
    EXPECT_THROW(TRestDataSetTardis("tardis", negative), std::runtime_error);
    const YAML::Node noChannels =
        YAML::Load("{class: TRestDataSetTardis, channelObservable: rawToSignal_SignalsMaxAmplitude}");
    EXPECT_THROW(TRestDataSetTardis("tardis", noChannels), std::runtime_error);

    // Cuts without a manager
    TRestDataSetTardis withCuts("tardis", YAML::Load("{class: TRestDataSetTardis, cuts: [even]}"));
    TRestDataSet dataSet(WriteTardisDataSet("input.root"));
    EXPECT_THROW(withCuts.ObtainParameters(dataSet, nullptr), std::runtime_error);

    // Channels without a readout in the dataset
    TRestDataSetTardis withChannels(
        "tardis", YAML::Load("{class: TRestDataSetTardis, channels: rawToSignal_SignalsPhysicalID}"));
    EXPECT_THROW(withChannels.Apply(dataSet), std::runtime_error);
}

// ---------------------------------------------------------------------------
// TRestDataSetGainMap
// ---------------------------------------------------------------------------

namespace {

/// A dataset of one run with the columns charge, posX and posY: 6000 events in each half of the readout
/// (posX in (-10, 0) with 100 per keV, posX in (0, 10) with 200 per keV), 6 of every 7 of the peak of
/// 22.5 keV and the others of the peak of 8 keV (sigma 3 %), and 10 events without position (NaN)
std::string WriteGainMapDataSet() {
    const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string path = std::string(TESTS_WORK_DIR) + "/gainmap_" + test->name() + ".root";
    TFile file(path.c_str(), "RECREATE");
    TTree tree("AnalysisTree", "AnalysisTree");
    int runIndex = 0, eventID = 0;
    double charge = 0, posX = 0, posY = 0;
    tree.Branch("runIndex", &runIndex);
    tree.Branch("eventID", &eventID);
    tree.Branch("charge", &charge);
    tree.Branch("posX", &posX);
    tree.Branch("posY", &posY);
    TRandom3 random(17);
    for (int i = 0; i < 12010; i++) {
        eventID = i;
        const double energy = i % 7 == 0 ? 8 : 22.5;
        const bool right = i % 2 == 1;
        const double mean = energy * (right ? 200 : 100);
        charge = random.Gaus(mean, 0.03 * mean);
        posX = i >= 12000 ? kNaN : (right ? random.Uniform(0, 10) : random.Uniform(-10, 0));
        posY = random.Uniform(-10, 10);
        tree.Fill();
    }
    tree.Write();
    TObjString runs("- run: {runNumber: 1}\n  fileName: a.root\n  entries: 12010\n  processing: 0\n");
    file.WriteTObject(&runs, "RESTDataSetInfo");
    return path;
}

const std::string kGainMapConfig =
    "{class: TRestDataSetGainMap, verbose: error, observable: charge, spatialObservableX: posX, "
    "spatialObservableY: posY, splitX: [-10, 0, 10], splitY: [-10, 10], range: [0, 6000], nBins: 600, "
    "peaks: [{energy: 22.5 keV}, {energy: 8 keV}]}";

}  // namespace

TEST(DataSetGainMap, SegmentsCalibratedWithTheirGain) {
    const std::string path = WriteGainMapDataSet();
    TRestDataSet dataSet(path);
    TRestDataSetGainMap gainMap("gainMap", YAML::Load(kGainMapConfig));
    EXPECT_EQ(gainMap.GetNumberOfSegmentsX(), 2);
    EXPECT_EQ(gainMap.GetNumberOfSegmentsY(), 1);
    gainMap.ObtainParameters(dataSet, nullptr);

    // energy = slope * charge + intercept: 1 keV per 100 on the left, per 200 on the right
    EXPECT_EQ(gainMap.GetSegment(-5, 0), std::make_pair(0, 0));
    EXPECT_EQ(gainMap.GetSegment(5, 0), std::make_pair(1, 0));
    EXPECT_EQ(gainMap.GetSegment(20, 0), std::make_pair(1, 0));  // outside: the last segment
    EXPECT_NEAR(gainMap.GetSlope(-5, 0), 0.01, 1e-4);
    EXPECT_NEAR(gainMap.GetSlope(5, 0), 0.005, 5e-5);
    EXPECT_NEAR(gainMap.GetIntercept(-5, 0), 0, 0.2);
    EXPECT_NEAR(gainMap.GetIntercept(5, 0), 0, 0.2);
    // The whole readout: the windows from the maximum of its spectrum (the 22.5 keV peak on the left)
    EXPECT_NEAR(gainMap.fFullSlope, 0.01, 1e-4);

    // Spectra of the segments and linear fits of the two peaks
    ASSERT_NE(gainMap.GetSegmentSpectrum(1, 0), nullptr);
    EXPECT_NEAR(gainMap.GetSegmentSpectrum(1, 0)->Integral(), 6000, 1);
    ASSERT_NE(gainMap.GetSegmentLinearFit(0, 0), nullptr);
    EXPECT_EQ(gainMap.GetSegmentLinearFit(0, 0)->GetN(), 2);
    EXPECT_EQ(gainMap.GetSegmentSpectrum(2, 0), nullptr);

    // The columns: calibrated with the segment of the event (NaN without position) and the whole readout
    gainMap.Apply(dataSet);
    EXPECT_EQ(gainMap.GetEnergyColumn(), "gainMap_Energy");
    const auto charge = *dataSet.GetDataFrame().Take<double>("charge");
    const auto posX = *dataSet.GetDataFrame().Take<double>("posX");
    const auto energy = *dataSet.GetDataFrame().Take<double>("gainMap_Energy");
    const auto full = *dataSet.GetDataFrame().Take<double>("gainMap_EnergyNoSegmentation");
    ASSERT_EQ(energy.size(), 12010u);
    for (size_t i = 0; i < energy.size(); i++) {
        EXPECT_NEAR(full[i], gainMap.fFullSlope * charge[i] + gainMap.fFullIntercept, 1e-9);
        if (std::isnan(posX[i])) {
            EXPECT_TRUE(std::isnan(energy[i]));
            continue;
        }
        EXPECT_NEAR(energy[i], gainMap.GetSlope(posX[i], 0) * charge[i] + gainMap.GetIntercept(posX[i], 0),
                    1e-9);
    }

    // A refit of a peak of a segment keeps its calibration
    gainMap.RefitAt(5, 0, 22.5, {4000, 5000});
    EXPECT_NEAR(gainMap.GetSlope(5, 0), 0.005, 5e-5);

    // Stored in the dataset and read back with its objects
    gainMap.UpdateFile(path);
    TRestDataSetGainMap read;
    read.ReadFromFile(path, "gainMap");
    EXPECT_NEAR(read.GetSlope(5, 0), gainMap.GetSlope(5, 0), 1e-12);
    EXPECT_NEAR(read.fFullSlope, gainMap.fFullSlope, 1e-12);
    ASSERT_NE(read.GetSegmentSpectrum(0, 0), nullptr);
    EXPECT_NEAR(read.GetSegmentSpectrum(0, 0)->Integral(), 6000, 1);
}

TEST(DataSetGainMap, Errors) {
    YAML::Node noPeaks = YAML::Load(kGainMapConfig);
    noPeaks.remove("peaks");
    EXPECT_THROW(TRestDataSetGainMap("gainMap", noPeaks), std::runtime_error);
    YAML::Node badSplit = YAML::Load(kGainMapConfig);
    badSplit["splitX"] = YAML::Load("[0, 0]");
    EXPECT_THROW(TRestDataSetGainMap("gainMap", badSplit), std::runtime_error);
    YAML::Node noSpatial = YAML::Load(kGainMapConfig);
    noSpatial.remove("spatialObservableY");
    EXPECT_THROW(TRestDataSetGainMap("gainMap", noSpatial), std::runtime_error);

    // Segments without edges need the readout of the dataset (none here)
    YAML::Node segments = YAML::Load(kGainMapConfig);
    segments.remove("splitX");
    segments.remove("splitY");
    segments["segments"] = YAML::Load("[2, 2]");
    TRestDataSetGainMap withoutReadout("gainMap", segments);
    TRestDataSet dataSet(WriteGainMapDataSet());
    EXPECT_THROW(withoutReadout.ObtainParameters(dataSet, nullptr), std::runtime_error);
}
