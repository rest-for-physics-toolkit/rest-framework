/// Pipeline tests (GoogleTest, BUILD_TESTING), as the legacy pipeline validation: restManager runs the YAML
/// of the examples, as a user would, on real data cut to a few hundred events (rest-common-data, downloaded
/// by URL), each test in its own work directory, and the analysis tree of each output is compared event by
/// event with a reference of examples/test/reference: every numeric column (scalars and vectors; strings are
/// not compared), with a relative tolerance. With the environment variable REST_UPDATE_REFERENCE=1 the
/// references are written instead (after an intended change, committed with it).

#include <TFile.h>
#include <TLeaf.h>
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
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

/// Release of rest-common-data with the input files
const std::string kData =
    "https://github.com/rest-for-physics-toolkit/rest-common-data/releases/download/v1.0/";
/// Tolerance of the comparison: relative, or absolute for values close to 0
constexpr double kTolerance = 1e-5;
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
    const std::string line = "cd '" + WorkDirectory().string() + "' && " + RESTMANAGER_EXECUTABLE + " --c '" +
                             (fs::path(REPOSITORY_DIR) / config).string() + "' --i '" + input + "' --j 1 " +
                             options;
    std::cout << "Running: " << line << std::endl;
    return std::system(line.c_str()) == 0;
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

bool Equal(double value, double reference) {
    if (std::isnan(value) || std::isnan(reference)) return std::isnan(value) && std::isnan(reference);
    if (std::isinf(value) || std::isinf(reference)) return value == reference;
    return std::abs(value - reference) <= kTolerance * std::max(1., std::abs(reference));
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
        int printed = 0;
        size_t index = 0;
        for (size_t entry = 0; entry < sizes->size(); entry++) {
            for (int i = 0; i < (*sizes)[entry]; i++, index++) {
                if (Equal(output.values[index], (*values)[index])) continue;
                if (printed++ < kMaxPrinted) {
                    ADD_FAILURE() << *column << ", entry " << entry << ", value " << i << ": "
                                  << output.values[index] << " (reference " << (*values)[index] << ")";
                }
            }
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
