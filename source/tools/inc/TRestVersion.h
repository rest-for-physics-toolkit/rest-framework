#ifndef TRESTVERSION_H
#define TRESTVERSION_H

#include <string>

/// \brief Version of the code of REST (rest-for-physics-toolkit), stored in the files it writes (the
/// codeVersion, codeCommit and codeClean of TRestRun and TRestDataSet) and printed by the executables
/// (`--version`):
///
/// - the version, the VERSION of the project in CMakeLists.txt (e.g. 3.0.0); a release is the tag vX.Y.Z of
///   the commit with that version;
/// - the commit of the build (git rev-parse HEAD), "unknown" without git (a build from a tarball or an image
///   can give it with the CMake option REST_GIT_COMMIT);
/// - whether the tracked files had uncommitted changes (the files not tracked by git do not count).
///
/// They are obtained at each build (cmake/RestVersion.cmake), not only when CMake is configured: a build
/// after a new commit or a change gives the new values.
///
/// \author JuanAn García <juanangp@unizar.es>
namespace REST_Version {

/// \brief Version of the code, e.g. "3.0.0".
std::string GetVersion();

/// \brief Commit of the build (40 hexadecimal characters), "unknown" if not known.
std::string GetCommit();

/// \brief The last tag and the commits since it (git describe --tags, e.g. v2.3.1-5-gabc1234), or the
/// abbreviated commit without tags; empty without git.
std::string GetDescribe();

/// \brief Date of the commit of the build (e.g. "2026-11-18 12:00:00 +0100"), empty without git.
std::string GetCommitDate();

/// \brief Whether the commit is known.
bool IsCommitKnown();

/// \brief Whether the tracked files were committed when the code was built (false also if unknown).
bool IsClean();

/// \brief Whether the clean state is known (built from a git checkout).
bool IsCleanKnown();

/// \brief One line with the version, the commit and its state, e.g.
/// "rest-for-physics-toolkit 3.0.0 (commit 963b345..., 2026-11-18, uncommitted changes)".
std::string GetSummary();

}  // namespace REST_Version

#endif
