/// Unit tests of restDAQ (GoogleTest, BUILD_TESTING), without electronics: a FEMINOS data frame built word
/// by word, split into events (TryExtractNextEvent, also when it arrives in two parts) and decoded into a
/// TRestRawSignalEvent (ParseEventFromWords), and the configuration of a run (RunConfig) read from a YAML
/// file written in the work directory.

#include <TRestRawSignalEvent.h>
#include <gtest/gtest.h>

#include <deque>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "FEMINOSPacket.h"
#include "RunConfig.h"

namespace fs = std::filesystem;

namespace {

/// Words of a FEMINOS data frame with one event: timestamp 0x000300020001, event count 0x00050004, and the
/// channel 5 of the chip 2 of the card 1 with the samples 100, 200 and 300 from the time bin 10
std::vector<uint16_t> FeminosFrame() {
    const uint16_t hit = PFX_CARD_CHIP_CHAN_HIT_IX | (1 << 9) | (2 << 7) | 5;
    return {PFX_START_OF_DFRAME,
            40,  // start of frame and its size
            PFX_START_OF_EVENT,
            0x0001,
            0x0002,
            0x0003,  // timestamp, 16 bits each
            0x0004,
            0x0005,  // event count
            hit,
            static_cast<uint16_t>(PFX_TIME_BIN_IX | 10),
            PFX_ADC_SAMPLE | 100,
            PFX_ADC_SAMPLE | 200,
            PFX_ADC_SAMPLE | 300,
            PFX_END_OF_EVENT,
            0x0000,  // end of event and its trailer
            PFX_END_OF_FRAME};
}

/// Physical channel of the FEMINOS channel: channel + 72 * chip + 288 * card
constexpr int kChannel = 5 + 2 * 72 + 288;

/// A restDAQ configuration in the work directory (FEMINOS, one FEM) with extra keys of the run
std::string WriteConfig(const std::string& runKeys) {
    const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
    const fs::path path = fs::path(TESTS_WORK_DIR) / (std::string(test->name()) + ".yaml");
    std::ofstream(path) << "daq:\n"
                           "  class: TRestRawDAQMetadata\n"
                           "  verbose: warning\n"
                           "  electronics: FEMINOS\n"
                           "  fems: [{id: 0, IP: \"192.168.10.2\"}]\n"
                           "  TCM: \"192.168.10.1\"\n"
                           "  nEvents: 100\n"
                           "  updateRate: 2 s\n"
                           "run:\n"
                           "  class: TRestRun\n"
                           "  experimentName: Test\n"
                           "  runTag: DUMMY\n"
                           "  runType: background\n"
                           "  mainDataPath: \""
                        << TESTS_WORK_DIR << "/data/\"\n"
                        << runKeys;
    return path.string();
}

}  // namespace

TEST(FEMINOSPacket, EventExtractedAndDecoded) {
    const auto frame = FeminosFrame();
    EXPECT_TRUE(FEMINOSPacket::isDataFrame(const_cast<uint16_t*>(frame.data())));
    EXPECT_FALSE(FEMINOSPacket::isMFrame(const_cast<uint16_t*>(frame.data())));

    // The first part of the frame (up to the samples): no complete event yet
    const size_t end = frame.size() - 3;  // the end of event, its trailer and the end of frame are missing
    std::deque<uint16_t> buffer(frame.begin(), frame.begin() + end);
    size_t index = 0;
    std::vector<uint16_t> event;
    EXPECT_FALSE(FEMINOSPacket::TryExtractNextEvent(buffer, index, event));
    EXPECT_EQ(buffer.size(), end);

    // The rest arrives: the event up to its trailer, the end of frame stays in the buffer
    buffer.insert(buffer.end(), frame.begin() + end, frame.end());
    ASSERT_TRUE(FEMINOSPacket::TryExtractNextEvent(buffer, index, event));
    EXPECT_EQ(event.size(), frame.size() - 1);
    ASSERT_EQ(buffer.size(), 1u);
    EXPECT_EQ(buffer.front(), PFX_END_OF_FRAME);
    EXPECT_EQ(index, 0u);

    TRestRawSignalEvent signals;
    uint64_t timestamp = 0;
    uint32_t count = 0;
    FEMINOSPacket::ParseEventFromWords(event, signals, timestamp, count);
    EXPECT_EQ(timestamp, 0x000300020001ULL);
    EXPECT_EQ(count, 0x00050004u);
    ASSERT_EQ(signals.GetNumberOfSignals(), 1);
    const auto& signal = signals.GetSignalByID(kChannel);
    ASSERT_EQ(signal.GetData().size(), 512u);
    EXPECT_EQ(signal.GetData()[9], 0);
    EXPECT_EQ(signal.GetData()[10], 100);
    EXPECT_EQ(signal.GetData()[11], 200);
    EXPECT_EQ(signal.GetData()[12], 300);
    EXPECT_EQ(signal.GetData()[13], 0);
}

TEST(RunConfig, RunAndDAQFromTheConfiguration) {
    RunConfig config(WriteConfig("  runNumber: 7\n"), {});
    EXPECT_EQ(config.electronics, "FEMINOS");
    EXPECT_EQ(config.nEvents, 100);
    EXPECT_EQ(config.updateRateTime, 2);
    EXPECT_EQ(config.verboseLevel.load(), RunConfig::Verbosity::REST_Warning);
    ASSERT_NE(config.run, nullptr);
    EXPECT_EQ(config.run->GetRunNumber(), 7);
    // Without outputFileName, the default one, in mainDataPath
    const std::string output = config.GetOutputFileName();
    EXPECT_EQ(output.rfind(std::string(TESTS_WORK_DIR) + "/data/Run00007_DUMMY_Test_background_", 0), 0u)
        << output;

    // runNumber auto: a new number in each run
    RunConfig automatic(WriteConfig("  runNumber: auto\n"), {});
    automatic.NewRun();
    const int first = automatic.run->GetRunNumber();
    automatic.NewRun();
    EXPECT_EQ(automatic.run->GetRunNumber(), first + 1);

    // Only the TCM
    RunConfig tcm(WriteConfig(""), {}, true);
    EXPECT_EQ(tcm.electronics, "TCM");
    EXPECT_EQ(tcm.run, nullptr);

    // Errors: unknown key of the run, no configuration file
    EXPECT_THROW(RunConfig(WriteConfig("  unknownKey: 1\n"), {}), std::exception);
    EXPECT_THROW(RunConfig(std::string(TESTS_WORK_DIR) + "/missing.yaml", {}), std::exception);
}
