/// Pipeline tests (GoogleTest, BUILD_TESTING), as the legacy pipeline validation: restManager runs the YAML
/// of the examples, as a user would, on the data of rest-common-data (downloaded by URL): real data cut to a
/// few hundred events and fixed restG4 simulations (the processing is tested, not Geant4; the random
/// processes of the examples have fixed seeds), each test in its own work directory, and the analysis tree
/// of each output is compared event by event with a reference of examples/test/reference: every numeric
/// column (scalars and vectors; strings are not compared), with a relative tolerance (looser for the results
/// of fits, which depend on the platform). With the environment variable REST_UPDATE_REFERENCE=1 the
/// references are written instead (after an intended change, committed with it).

#include <TFile.h>
#include <TLeaf.h>
#include <TRestGDMLReadout.h>
#include <TRestTools.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

/// Release of rest-common-data with the input files
const std::string kData =
    "https://github.com/rest-for-physics-toolkit/rest-common-data/releases/download/v1.1/";
/// Tolerance of the comparison: relative, or absolute for values close to 0
constexpr double kTolerance = 1e-5;
/// Results of fits (columns whose name contains "Gauss", e.g. trackAnalysis_MaxTrackGaussSigmaX): Minuit
/// converges to slightly different points on other platforms (libm, architecture), and a fit may fail on one
/// of them only. Their tolerance, and the fraction of their values allowed to differ
constexpr double kFitTolerance = 1e-2;
constexpr double kFitMismatches = 0.02;
const std::regex kFitColumn(".*Gauss.*");
/// Mismatches printed per column
constexpr int kMaxPrinted = 5;

/// Work directory of the current test (work/<suite>_<test>), created again for each test
fs::path WorkDirectory() {
    const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
    const fs::path directory =
        fs::path(TESTS_WORK_DIR) / (std::string(test->test_suite_name()) + "_" + test->name());
    static std::string created;
    if (created != directory.string()) {
        fs::remove_all(directory);
        fs::create_directories(directory);
        created = directory.string();
    }
    return directory;
}

/// A file of rest-common-data in the work directory with its own name (the run parameters are read from
/// it): downloaded once (cached) and linked
std::string Input(const std::string& name) {
    const fs::path link = WorkDirectory() / name;
    if (!fs::exists(link)) fs::create_symlink(TRestTools::DownloadFile(kData + name), link);
    return name;
}

/// restManager in the work directory with a configuration of the repository (path from its root), one
/// thread (the order of the events) and other options; true if it succeeds
bool RestManager(const std::string& config, const std::string& input, const std::string& options = "") {
    // Without display (ROOT in batch mode) and without terminal input: the plots with previewPlot do not
    // open a window nor wait for Enter (TRestAnalysisPlot)
    const std::string configPath = (fs::path(REPOSITORY_DIR) / config).string();
    const std::string line = "cd '" + WorkDirectory().string() + "' && env -u DISPLAY -u WAYLAND_DISPLAY " +
                             RESTMANAGER_EXECUTABLE + " --c '" + configPath + "' --i '" + input + "' --j 1 " +
                             options + " < /dev/null";
    std::cout << "Running: " << line << std::endl;
    return std::system(line.c_str()) == 0;
}

/// The veto readout of libraries/detector/examples/GDMLReadout.yaml in a file, as the macro
/// REST_Detector_GenerateGDMLReadout, built from the repository root (the paths of the YAML); true if done
bool GenerateVetoReadout(const std::string& fileName) {
    const fs::path previous = fs::current_path();
    fs::current_path(REPOSITORY_DIR);
    bool done = false;
    try {
        TRestGDMLReadout readout("libraries/detector/examples/GDMLReadout.yaml", "readout");
        readout.BuildGeometry();
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "RECREATE"));
        done = file && !file->IsZombie() && readout.Export(file.get(), readout.GetName(), "veto_decoding");
    } catch (const std::exception& error) {
        std::cerr << "GenerateVetoReadout: " << error.what() << std::endl;
    }
    fs::current_path(previous);
    return done;
}

/// The only file of the work directory whose name matches a regular expression (empty if none or several)
std::string Output(const std::string& pattern) {
    std::vector<std::string> found;
    for (const auto& entry : fs::directory_iterator(WorkDirectory())) {
        if (std::regex_match(entry.path().filename().string(), std::regex(pattern))) {
            found.push_back(entry.path().string());
        }
    }
    EXPECT_EQ(found.size(), 1u) << "files matching " << pattern;
    return found.size() == 1 ? found[0] : "";
}

/// Values of a numeric column of an analysis tree: the values of all the entries one after the other and the
/// number of values of each entry (1 for a scalar, the size for a vector)
struct Column {
    std::vector<double> values;
    std::vector<int> sizes;
};

/// Whether a branch holds numbers (a basic type or a std::vector of one), not strings or other objects
bool IsNumeric(TBranch* branch) {
    static const std::regex number(
        "(Double_t|Float_t|Int_t|UInt_t|Long64_t|ULong64_t|Long_t|ULong_t|Short_t|UShort_t|Bool_t|double|"
        "float|int|unsigned int|long|unsigned long|long long|unsigned long long|short|unsigned short|bool)");
    const std::string className = branch->GetClassName();
    if (!className.empty()) {
        std::smatch match;
        static const std::regex vector("vector<(.+)>");
        return std::regex_match(className, match, vector) && std::regex_match(match[1].str(), number);
    }
    auto* leaves = branch->GetListOfLeaves();
    if (!leaves || leaves->GetEntries() != 1) return false;
    return std::regex_match(std::string(static_cast<TLeaf*>(leaves->At(0))->GetTypeName()), number);
}

/// The numeric columns of the AnalysisTree of a file
std::map<std::string, Column> ReadAnalysisTree(const std::string& fileName) {
    std::map<std::string, Column> columns;
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        ADD_FAILURE() << "cannot open " << fileName;
        return columns;
    }
    auto* tree = file->Get<TTree>("AnalysisTree");
    if (!tree) {
        ADD_FAILURE() << "no AnalysisTree in " << fileName;
        return columns;
    }
    // An output without events would pass any comparison
    if (tree->GetEntries() == 0) {
        ADD_FAILURE() << "AnalysisTree without entries in " << fileName;
        return columns;
    }
    for (auto* object : *tree->GetListOfBranches()) {
        auto* branch = static_cast<TBranch*>(object);
        if (!IsNumeric(branch)) continue;
        TTreeFormula formula(branch->GetName(), branch->GetName(), tree);
        if (formula.GetNdim() == 0) continue;
        Column& column = columns[branch->GetName()];
        for (Long64_t entry = 0; entry < tree->GetEntries(); entry++) {
            tree->LoadTree(entry);
            const int size = formula.GetNdata();
            column.sizes.push_back(size);
            for (int i = 0; i < size; i++) column.values.push_back(formula.EvalInstance(i));
        }
    }
    return columns;
}

bool Equal(double value, double reference, double tolerance) {
    if (std::isnan(value) || std::isnan(reference)) return std::isnan(value) && std::isnan(reference);
    if (std::isinf(value) || std::isinf(reference)) return value == reference;
    return std::abs(value - reference) <= tolerance * std::max(1., std::abs(reference));
}

/// Compares the analysis tree of an output with the reference `name` (examples/test/reference/<name>.root),
/// or writes it with REST_UPDATE_REFERENCE=1
void CompareWithReference(const std::string& fileName, const std::string& name) {
    ASSERT_FALSE(fileName.empty());
    const auto columns = ReadAnalysisTree(fileName);
    ASSERT_FALSE(columns.empty()) << "no numeric column in the AnalysisTree of " << fileName;
    const fs::path referenceFile = fs::path(TESTS_SOURCE_DIR) / "reference" / (name + ".root");

    const char* update = std::getenv("REST_UPDATE_REFERENCE");
    if (update && std::string(update) == "1") {
        fs::create_directories(referenceFile.parent_path());
        // LZMA: the references are committed
        TFile file(referenceFile.c_str(), "RECREATE", "", 207);
        TTree tree("Reference", "Numeric columns of the AnalysisTree");
        std::string column;
        std::vector<double> values;
        std::vector<int> sizes;
        tree.Branch("column", &column);
        tree.Branch("values", &values);
        tree.Branch("sizes", &sizes);
        for (const auto& [columnName, data] : columns) {
            column = columnName;
            values = data.values;
            sizes = data.sizes;
            tree.Fill();
        }
        tree.Write();
        std::cout << "Reference written: " << referenceFile << " (" << columns.size() << " columns)"
                  << std::endl;
        return;
    }

    std::unique_ptr<TFile> file(TFile::Open(referenceFile.c_str(), "READ"));
    ASSERT_TRUE(file && !file->IsZombie())
        << "no reference " << referenceFile << ": write it with REST_UPDATE_REFERENCE=1";
    auto* tree = file->Get<TTree>("Reference");
    ASSERT_NE(tree, nullptr);
    std::string* column = nullptr;
    std::vector<double>* values = nullptr;
    std::vector<int>* sizes = nullptr;
    tree->SetBranchAddress("column", &column);
    tree->SetBranchAddress("values", &values);
    tree->SetBranchAddress("sizes", &sizes);
    std::set<std::string> referenced;
    for (Long64_t n = 0; n < tree->GetEntries(); n++) {
        tree->GetEntry(n);
        referenced.insert(*column);
        const auto found = columns.find(*column);
        if (found == columns.end()) {
            ADD_FAILURE() << *column << ": in the reference, not in the output";
            continue;
        }
        const Column& output = found->second;
        if (output.sizes != *sizes) {
            ADD_FAILURE() << *column << ": " << output.sizes.size() << " entries (reference " << sizes->size()
                          << ") or different number of values in an entry";
            continue;
        }
        const bool fit = std::regex_match(*column, kFitColumn);
        const double tolerance = fit ? kFitTolerance : kTolerance;
        // The mismatches allowed: a fraction of the values of a fit column, none otherwise
        const size_t allowed = fit ? static_cast<size_t>(kFitMismatches * values->size()) : 0;
        std::vector<std::string> mismatches;
        size_t index = 0;
        for (size_t entry = 0; entry < sizes->size(); entry++) {
            for (int i = 0; i < (*sizes)[entry]; i++, index++) {
                if (Equal(output.values[index], (*values)[index], tolerance)) continue;
                std::ostringstream text;
                text << *column << ", entry " << entry << ", value " << i << ": " << output.values[index]
                     << " (reference " << (*values)[index] << ")";
                mismatches.push_back(text.str());
            }
        }
        if (mismatches.size() <= allowed) {
            if (!mismatches.empty()) {
                std::cout << *column << ": " << mismatches.size() << " fit values differ (allowed " << allowed
                          << ")" << std::endl;
            }
            continue;
        }
        ADD_FAILURE() << *column << ": " << mismatches.size() << " values differ (allowed " << allowed << ")";
        for (size_t m = 0; m < mismatches.size() && m < static_cast<size_t>(kMaxPrinted); m++) {
            ADD_FAILURE() << mismatches[m];
        }
    }
    for (const auto& [columnName, data] : columns) {
        if (!referenced.count(columnName)) {
            ADD_FAILURE() << columnName << ": in the output, not in the reference";
        }
    }
}

}  // namespace

// IAXO-D1 data of the FEMDAQ: raw signals to tracks
TEST(Pipelines, IAXOFemDAQToTracks) {
    ASSERT_TRUE(
        RestManager("examples/FemDAQToTracks.yaml", Input("Run00229_Nominal_Test_Calibration_008.root")));
    CompareWithReference(Output("R00229_008_ANA_.*\\.root"), "FemDAQToTracks");
}

// TREX-DM data of the FEMDAQ: both planes to tracks, then the time since the previous event (Tardis) of the
// dataset of the plane 0 and its plots
TEST(Pipelines, TREXDMToTracksAndTardis) {
    const std::string input = Input("Run03318_Calibration109Cd_TREXDM_DAQ_001.root");
    for (const std::string plane : {"0", "1"}) {
        ASSERT_TRUE(RestManager("examples/TREXDMToTracks.yaml", input, "--set plane=" + plane));
        CompareWithReference(Output("Run03318_Calibration109Cd_TREXDM_Plane" + plane + "_001\\.root"),
                             "TREXDMToTracks_Plane" + plane);
    }
    ASSERT_TRUE(RestManager("examples/TREXDMTardis.yaml", "Run*_TREXDM_Plane0_*.root"));
    const std::string dataset = Output("Dataset_03318_Calibration109Cd_TREXDM_Plane0_Tardis\\.root");
    CompareWithReference(dataset, "TREXDMTardis_Plane0");
    ASSERT_FALSE(dataset.empty());
    EXPECT_TRUE(RestManager("examples/TREXDMTardisPlots.yaml", fs::path(dataset).filename().string()));
}

// .aqs file of the Feminos cards (TREX-DM, TCM mode): raw signals and their analysis
TEST(Pipelines, MultiFEMINOSToRaw) {
    ASSERT_TRUE(RestManager("libraries/raw/examples/MultiFEMINOSToRaw.yaml",
                            Input("R01928_tagTest_Vm_250_Vd_160_Pr_6_Gain_0x0_Shape_0xF_Clock_0x4-068.aqs")));
    ASSERT_TRUE(RestManager("libraries/raw/examples/RawSignalAnalysis.yaml", "RawSignals.root"));
    CompareWithReference(Output("R.*_RAWANA\\.root"), "MultiFEMINOSRawAnalysis");
}

// ROOT file of feminos-daq (TREX-DM): both planes to tracks
TEST(Pipelines, FeminosRootToTracks) {
    const std::string input =
        Input("R03137_Calibration37Ar_Vm_360_Vd_100_Pr_2_Gain_0x0_Shape_0xF_Clock_0x4.root");
    for (const std::string plane : {"0", "1"}) {
        ASSERT_TRUE(
            RestManager("libraries/raw/examples/FeminosRootToRaw.yaml", input, "--set plane=" + plane));
        CompareWithReference(Output("FeminosRoot_Plane" + plane + "\\.root"), "FeminosRoot_Plane" + plane);
    }
}

// Simulated cosmic neutrons in IAXO-D1 (restG4 example 13.IAXO, Neutrons.yaml): quenching, then the Geant4
// analysis and the hits of the gas and of the vetoes
TEST(Pipelines, IAXOSimulatedNeutrons) {
    ASSERT_TRUE(
        RestManager("libraries/geant4/examples/Quenching.yaml", Input("Run00010_CosmicNeutrons_IAXO.root")));
    const std::string quenching = "Run00010_CosmicNeutrons_IAXO_QUENCHING.root";
    CompareWithReference(Output(quenching), "Quenching");
    ASSERT_TRUE(RestManager("libraries/geant4/examples/Geant4Analysis.yaml", quenching));
    CompareWithReference(Output("Run00010_CosmicNeutrons_IAXO_G4ANA\\.root"), "Geant4Analysis");
    ASSERT_TRUE(RestManager("libraries/connectors/examples/Geant4ToHits.yaml", quenching));
    CompareWithReference(Output("Run00010_CosmicNeutrons_IAXO_HITS\\.root"), "Geant4ToHits");
}

// Simulated cosmic muons in IAXO-D1 (13.IAXO, Muons.yaml): raw signals of the vetoes, with their readout
// generated from the GDML of the example (GenerateVetoReadout)
TEST(Pipelines, IAXOSimulatedVetoes) {
    ASSERT_TRUE(GenerateVetoReadout((WorkDirectory() / "veto_readout_output.root").string()));
    ASSERT_TRUE(RestManager("libraries/connectors/examples/Geant4ToRawVeto.yaml",
                            Input("Run00001_CosmicMuons_IAXO.root")));
    CompareWithReference(Output("Run00001_CosmicMuons_IAXO_VETO_RAW\\.root"), "Geant4ToRawVeto");
}

// Simulated 241Am calibration of AlphaCAMM (Simulation.yaml of the AlphaCAMM repository): raw signals of the
// Micromegas, and the drifted hits
TEST(Pipelines, AlphaCAMMSimulatedCalibration) {
    const std::string input = Input("Run00001_Am241_calibration_AlphaCAMM.root");
    ASSERT_TRUE(RestManager("libraries/connectors/examples/Geant4ToRaw.yaml", input));
    CompareWithReference(Output("Run00001_Am241_calibration_AlphaCAMM_RAW\\.root"), "Geant4ToRaw");
    ASSERT_TRUE(RestManager("libraries/connectors/examples/Geant4ToDriftHits.yaml", input));
    CompareWithReference(Output("Run00001_Am241_calibration_AlphaCAMM_DRIFT\\.root"), "Geant4ToDriftHits");
}
