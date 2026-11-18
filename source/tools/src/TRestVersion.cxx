#include "TRestVersion.h"

// Written at each build by cmake/RestVersion.cmake (build directory of source/tools)
#include "RestVersionInfo.h"

namespace REST_Version {

std::string GetVersion() { return REST_VERSION_STRING; }

std::string GetCommit() { return REST_GIT_COMMIT; }

std::string GetDescribe() { return REST_GIT_DESCRIBE; }

std::string GetCommitDate() { return REST_GIT_DATE; }

bool IsCommitKnown() { return GetCommit() != "unknown"; }

bool IsClean() { return REST_GIT_CLEAN == 1; }

bool IsCleanKnown() { return REST_GIT_CLEAN >= 0; }

std::string GetSummary() {
    std::string summary = "rest-for-physics-toolkit " + GetVersion() + " (commit " + GetCommit();
    if (!GetDescribe().empty() && GetDescribe() != GetCommit().substr(0, GetDescribe().size())) {
        summary += ", " + GetDescribe();
    }
    if (!GetCommitDate().empty()) summary += ", " + GetCommitDate().substr(0, 10);
    if (IsCleanKnown()) summary += IsClean() ? ", clean" : ", uncommitted changes";
    return summary + ")";
}

}  // namespace REST_Version
