#!/usr/bin/env bash
# Installer of rest-for-physics-toolkit (rest-framework): clones a release, builds and installs it with its
# dependencies, and writes <prefix>/setup.sh, the environment to source (dependencies and REST).
#
#   install.sh [--mode system|cvmfs|source] [--prefix <dir>] [--version <vX.Y.Z|master>] [options]
#
# Modes:
#   system  ROOT (required), Geant4 and Garfield++ (optional) already in the environment, e.g. a server; --env
#           sources a script loading them first (also in setup.sh). Default.
#   cvmfs   The LCG view of /cvmfs/sft.cern.ch given by --lcg (ROOT, Geant4, GSL...), e.g. lxplus or NAF.
#           Garfield++ GARFIELD_VERSION (that of the CI) is built on top of it: the view has another version
#           (LCG_109: a later commit), used instead with --lcg-garfield.
#   source  Everything with the versions of the rest-env image of the CI (the rest-docker Dockerfile): the
#           packages of the system checked (apt, dnf or brew; with sudo only if --install-packages, otherwise
#           the command for the administrator is printed), ROOT (the binary of root.cern of the CI on
#           Ubuntu 24.04, built elsewhere), Geant4, Garfield++ and CRY. Any Linux or macOS.
# In every mode CRY (cry generator of restG4) and yaml-cpp are built if missing.
#
# Layout of <prefix> (default ~/rest): src/rest-framework (the clone), build/, install/ (REST_PATH), deps/
# (the dependencies built), setup.sh and install.log. Run again with the same prefix it continues: built
# dependencies are kept, and the clone (which must have no local changes) is moved to the version asked.

# The whole script in a block, read by bash before running it: a change of the file while it runs (e.g. a git
# pull) does not reach the run (bash reads a script as it runs it)
{
set -euo pipefail

# Versions of the rest-env image (ARG of the rest-docker Dockerfile): keep them equal
ROOT_VERSION=6.36.04
GEANT4_VERSION=11.0.3
GARFIELD_VERSION=4.0
CRY_VERSION=1.7
YAML_CPP_VERSION=0.8.0
DEFAULT_LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt
REPOSITORY=https://github.com/rest-for-physics-toolkit/rest-framework.git

usage() {
    cat <<EOF
Usage: install.sh [options]

  --mode <system|cvmfs|source>  where the dependencies come from (default system)
  --prefix <dir>                installation directory (default ~/rest)
  --version <vX.Y.Z|master>     version of rest-framework (default the last release)
  --env <script>                (system) script loading ROOT, Geant4... sourced before building
  --lcg <view>                  (cvmfs) LCG view (default $DEFAULT_LCG)
  --lcg-garfield                (cvmfs) the Garfield++ of the view, not $GARFIELD_VERSION (that of the CI)
  --jobs <n>                    parallel jobs (default: CPUs, at most one per 2 GB of memory)
  --install-packages            (source) install the missing packages of the system with sudo
  --daq                         also restDAQ (femclient; needs readline)
  --no-g4                       without restG4 (and CRY)
  --no-garfield                 without restGarfield
  -h, --help                    this help
EOF
}

info() { echo -e "\033[1;34m[install]\033[0m $*"; }
fail() {
    echo -e "\033[1;31m[install] error:\033[0m $*" >&2
    exit 1
}

mode=system prefix="$HOME/rest" version="" env_script="" lcg="$DEFAULT_LCG" jobs="" daq=OFF g4=ON garfield=ON
install_packages=OFF lcg_garfield=OFF
while [[ $# -gt 0 ]]; do
    case "$1" in
    --mode) mode=${2:?}; shift ;;
    --prefix) prefix=${2:?}; shift ;;
    --version) version=${2:?}; shift ;;
    --env) env_script=${2:?}; shift ;;
    --lcg) lcg=${2:?}; shift ;;
    --jobs) jobs=${2:?}; shift ;;
    --install-packages) install_packages=ON ;;
    --lcg-garfield) lcg_garfield=ON ;;
    --daq) daq=ON ;;
    --no-g4) g4=OFF ;;
    --no-garfield) garfield=OFF ;;
    -h | --help) usage; exit 0 ;;
    *) usage; fail "unknown option $1" ;;
    esac
    shift
done
[[ "$mode" =~ ^(system|cvmfs|source)$ ]] || fail "--mode must be system, cvmfs or source"
[[ -z "$env_script" || "$mode" == system ]] || fail "--env is only for --mode system"
[[ "$install_packages" == OFF || "$mode" == source ]] || fail "--install-packages is only for --mode source"
[[ "$lcg_garfield" == OFF || "$mode" == cvmfs ]] || fail "--lcg-garfield is only for --mode cvmfs"

mkdir -p "$prefix"
prefix=$(cd "$prefix" && pwd)
deps="$prefix/deps"
mkdir -p "$deps"
exec > >(tee -a "$prefix/install.log") 2>&1
info "$(date): mode $mode, prefix $prefix"

if [[ -z "$jobs" ]]; then
    cpus=$(nproc 2>/dev/null || sysctl -n hw.ncpu)
    memory=$(awk '/MemTotal/ { print int($2 / 2097152) }' /proc/meminfo 2>/dev/null \
        || echo $(($(sysctl -n hw.memsize) / 2147483648)))
    jobs=$((cpus < memory ? cpus : memory))
    jobs=$((jobs < 1 ? 1 : jobs))
fi
info "$jobs parallel jobs"
generator=("Unix Makefiles")
if command -v ninja >/dev/null; then generator=(Ninja); fi

# Lines of setup.sh before thisREST.sh: the environment of the dependencies
setup_lines=()

# A script of the environment of a dependency: sourced now and in setup.sh (they use unset variables)
load() {
    set +u
    # shellcheck disable=SC1090
    source "$1"
    set -u
    setup_lines+=("source \"$1\"")
}

# Clones a git repository at a tag, retried (gitlab.cern.ch and github.com sometimes reset connections)
clone() {
    local url=$1 tag=$2 dir=$3
    for attempt in 1 2 3 4 5; do
        git clone -q --depth 1 --branch "$tag" "$url" "$dir" && return 0
        rm -rf "$dir"
        info "clone of $url failed (attempt $attempt), retrying"
        sleep 30
    done
    fail "cannot clone $url"
}

# --- Dependencies -------------------------------------------------------------------------------------

# Packages of the system for ROOT, Geant4, Garfield++ and REST (those of the rest-env image), with the
# package manager found: apt (Debian, Ubuntu), dnf (RHEL, AlmaLinux, Rocky, Fedora; some packages in EPEL)
# or brew (macOS). The packages that the distribution does not have are skipped (the builds use their
# builtin copies or go without the feature)
APT_PACKAGES="build-essential gfortran cmake ninja-build git ca-certificates curl wget file
    python3 python3-dev python3-numpy libssl-dev
    libx11-dev libxext-dev libxft-dev libxpm-dev libxmu-dev libxi-dev
    libgl-dev libglu1-mesa-dev libglew-dev libftgl-dev libgl2ps-dev
    libtbb-dev libvdt-dev libgif-dev libjpeg-dev libpng-dev libtiff-dev
    libfftw3-dev libgsl-dev libxml2-dev libpcre3-dev libpcre2-dev libcfitsio-dev
    libxxhash-dev libzstd-dev liblz4-dev liblzma-dev zlib1g-dev nlohmann-json3-dev
    libxerces-c-dev libexpat1-dev libcurl4-openssl-dev
    libyaml-cpp-dev libgtest-dev libgmock-dev libreadline-dev"
DNF_PACKAGES="gcc gcc-c++ gcc-gfortran make cmake ninja-build git ca-certificates curl wget file
    python3 python3-devel python3-numpy openssl-devel
    libX11-devel libXext-devel libXft-devel libXpm-devel libXmu-devel libXi-devel
    mesa-libGL-devel mesa-libGLU-devel glew-devel ftgl-devel gl2ps-devel
    tbb-devel vdt-devel giflib-devel libjpeg-turbo-devel libpng-devel libtiff-devel
    fftw-devel gsl-devel libxml2-devel pcre2-devel cfitsio-devel
    xxhash-devel libzstd-devel lz4-devel xz-devel zlib-devel json-devel
    xerces-c-devel expat-devel libcurl-devel
    yaml-cpp-devel gtest-devel gmock-devel readline-devel"
BREW_PACKAGES="cmake ninja git gcc python numpy openssl@3 tbb gsl fftw xxhash zstd lz4 xz
    nlohmann-json xerces-c expat yaml-cpp googletest readline pcre2 libpng jpeg-turbo giflib"

# Checks the packages of the system (installed, or not available in the distribution); the missing ones
# are installed only with --install-packages (sudo), otherwise the command is printed (for the
# administrator) and the installation goes on: the dependencies built in deps/ need no privileges, and a
# missing package needed makes their build fail
install_system_packages() {
    local missing=() package manager command sudo=sudo
    [[ "$(id -u)" -ne 0 ]] || sudo=""
    if command -v apt-get >/dev/null; then
        manager=apt
        # The package lists first: without them (e.g. a new container) no package would be available
        [[ "$install_packages" == OFF ]] || $sudo apt-get update
        # Available: a version to install (apt-cache show exits with 0 also for a package only referred to
        # by others, e.g. libpcre3-dev in Debian 13)
        for package in $APT_PACKAGES; do
            dpkg-query -W -f='${Status}' "$package" 2>/dev/null | grep -q "install ok installed" \
                || { apt-cache show --no-all-versions "$package" 2>/dev/null | grep -q '^Package:' \
                    && missing+=("$package"); }
        done
        command="sudo apt-get install -y --no-install-recommends ${missing[*]:-}"
    elif command -v dnf >/dev/null; then
        manager=dnf
        for package in $DNF_PACKAGES; do
            rpm -q "$package" >/dev/null 2>&1 \
                || { dnf -q info "$package" >/dev/null 2>&1 && missing+=("$package"); }
        done
        command="sudo dnf install -y ${missing[*]:-}"
    elif command -v brew >/dev/null; then
        # Homebrew installs in the user's prefix, without privileges
        manager=brew
        for package in $BREW_PACKAGES; do
            brew list --versions "$package" >/dev/null 2>&1 \
                || { brew info --formula "$package" >/dev/null 2>&1 && missing+=("$package"); }
        done
        [[ ${#missing[@]} -eq 0 ]] || brew install "${missing[@]}"
        export CMAKE_PREFIX_PATH="$(brew --prefix)${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
        missing=()
    else
        info "no apt, dnf or brew: the packages of APT_PACKAGES in this script (or their equivalents) must" \
            "be installed"
    fi
    if [[ ${#missing[@]} -gt 0 ]]; then
        if [[ "$install_packages" == ON ]]; then
            info "system packages ($manager): ${missing[*]}"
            if [[ "$manager" == apt ]]; then
                $sudo env DEBIAN_FRONTEND=noninteractive \
                    apt-get install -y --no-install-recommends "${missing[@]}"
            else
                $sudo dnf install -y "${missing[@]}"
            fi
        else
            info "packages of the system missing (ask the administrator or use --install-packages):"
            echo "    $command"
        fi
    fi
    if command -v ninja >/dev/null; then generator=(Ninja); fi
}

# Adds #include <cstdint> to the installed headers using fixed width integers without it: Geant4 11.0 and
# Garfield++ 4.0 use them through headers that GCC >= 13 no longer includes (as the rest-env image)
include_cstdint() {
    local header
    grep -rlE '\b(std::)?u?int(8|16|32|64)_t\b' "$1" | while read -r header; do
        grep -q '<cstdint>' "$header" && continue
        { echo '#include <cstdint>'; cat "$header"; } >"$header.tmp" && mv "$header.tmp" "$header"
    done
}

# ROOT: on Ubuntu 24.04 x86_64 the binary of root.cern used by the CI (rest-env image); elsewhere the same
# version built from its sources (about an hour) with the components used by REST
install_root() {
    if [[ ! -x "$deps/root/bin/root-config" ]]; then
        if [[ "$(uname -m)" == x86_64 ]] && grep -q '^ID=ubuntu' /etc/os-release 2>/dev/null \
            && grep -q 'VERSION_ID="24.04"' /etc/os-release; then
            info "ROOT $ROOT_VERSION (binary of root.cern, as the CI)"
            local tarball="root_v${ROOT_VERSION}.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz"
            curl -fL --retry 5 -o "$deps/$tarball" "https://root.cern/download/$tarball"
            tar -xzf "$deps/$tarball" -C "$deps" && rm "$deps/$tarball"
        else
            info "ROOT $ROOT_VERSION (built from its sources, it takes a while)"
            rm -rf "$deps/root-src" "$deps/root-build"
            clone https://github.com/root-project/root.git "v${ROOT_VERSION//./-}" "$deps/root-src"
            cmake -S "$deps/root-src" -B "$deps/root-build" -G "${generator[@]}" \
                -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$deps/root" -DCMAKE_CXX_STANDARD=17 \
                -Dgdml=ON -Dopengl=ON -Droofit=ON -Dtmva=ON -Dimt=ON -Dmathmore=ON -Dxrootd=OFF \
                -Dbuiltin_xrootd=OFF
            cmake --build "$deps/root-build" -j "$jobs"
            cmake --install "$deps/root-build"
            rm -rf "$deps/root-src" "$deps/root-build"
        fi
    fi
    load "$deps/root/bin/thisroot.sh"
}

install_geant4() {
    if [[ ! -f "$deps/geant4/bin/geant4.sh" ]]; then
        info "Geant4 $GEANT4_VERSION (with its data, it takes a while)"
        rm -rf "$deps/geant4-src" "$deps/geant4-build"
        clone https://github.com/Geant4/geant4.git "v$GEANT4_VERSION" "$deps/geant4-src"
        # <cstdint> in every file (GCC >= 13); the Expat of Geant4 (with the one of the system, Geant4 11.0
        # writes an empty EXPAT_LIBRARY with CMake >= 3.27, which breaks find_package(Geant4))
        local x11=ON
        [[ "$(uname)" != Darwin ]] || x11=OFF
        cmake -S "$deps/geant4-src" -B "$deps/geant4-build" -G "${generator[@]}" \
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$deps/geant4" \
            -DGEANT4_USE_GDML=ON -DGEANT4_INSTALL_DATA=ON -DGEANT4_BUILD_MULTITHREADED=ON \
            -DGEANT4_USE_SYSTEM_EXPAT=OFF -DGEANT4_USE_OPENGL_X11="$x11" -DCMAKE_CXX_FLAGS="-include cstdint"
        cmake --build "$deps/geant4-build" -j "$jobs"
        cmake --install "$deps/geant4-build"
        include_cstdint "$deps/geant4/include"
        rm -rf "$deps/geant4-src" "$deps/geant4-build"
    fi
    load "$deps/geant4/bin/geant4.sh"
}

install_garfield() {
    if [[ ! -f "$deps/garfield/.installed" ]]; then
        info "Garfield++ $GARFIELD_VERSION"
        rm -rf "$deps/garfield" "$deps/garfield-src" "$deps/garfield-build"
        clone https://gitlab.cern.ch/garfield/garfieldpp.git "$GARFIELD_VERSION" "$deps/garfield-src"
        # Makefiles, not Ninja (the Fortran of Magboltz through the C preprocessor, as the rest-env image).
        # <cstdint> and <algorithm> in every file: its sources use them without including them, reached
        # through other headers only on some systems (e.g. not with the GCC 13 of the LCG views: std::find)
        cmake -S "$deps/garfield-src" -B "$deps/garfield-build" -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$deps/garfield" -DWITH_EXAMPLES=OFF \
            -DCMAKE_CXX_FLAGS="-include cstdint -include algorithm"
        cmake --build "$deps/garfield-build" -j "$jobs"
        cmake --install "$deps/garfield-build"
        include_cstdint "$deps/garfield/include"
        mkdir -p "$deps/garfield/share/Heed"
        [[ -d "$deps/garfield/share/Heed/database" ]] \
            || cp -r "$deps/garfield-src/Heed/heed++/database" "$deps/garfield/share/Heed/database"
        rm -rf "$deps/garfield-src" "$deps/garfield-build"
        touch "$deps/garfield/.installed"
    fi
    local lib="$deps/garfield/lib"
    [[ -d "$deps/garfield/lib64" ]] && lib="$deps/garfield/lib64"
    export GARFIELD_INSTALL="$deps/garfield" HEED_DATABASE="$deps/garfield/share/Heed/database"
    export LD_LIBRARY_PATH="$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    setup_lines+=("export GARFIELD_INSTALL=\"$GARFIELD_INSTALL\" HEED_DATABASE=\"$HEED_DATABASE\"")
    setup_lines+=("export LD_LIBRARY_PATH=\"$lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}\"")
}

install_cry() {
    if [[ ! -f "$deps/cry/lib/libCRY.a" ]]; then
        info "CRY $CRY_VERSION"
        rm -rf "$deps/cry"
        curl -fL --retry 5 "https://nuclear.llnl.gov/simulation/cry_v${CRY_VERSION}.tar.gz" \
            | tar -xz -C "$deps"
        mv "$deps/cry_v${CRY_VERSION}" "$deps/cry"
        make -C "$deps/cry"
    fi
    # thisREST.sh sets CRYHOME and CRYDATAPATH of the CRY found at configure time
    export CRYHOME="$deps/cry"
}

# yaml-cpp, if CMake does not find it (servers and LCG views may not have it)
install_yaml_cpp() {
    local probe="$deps/yaml-cpp-probe"
    mkdir -p "$probe"
    printf 'cmake_minimum_required(VERSION 3.18)\nproject(probe CXX)\nfind_package(yaml-cpp REQUIRED)\n' \
        >"$probe/CMakeLists.txt"
    if cmake -S "$probe" -B "$probe/build" >/dev/null 2>&1; then
        rm -rf "$probe"
        return
    fi
    rm -rf "$probe"
    if [[ ! -d "$deps/yaml-cpp" ]]; then
        info "yaml-cpp $YAML_CPP_VERSION"
        rm -rf "$deps/yaml-cpp-src"
        clone https://github.com/jbeder/yaml-cpp.git "$YAML_CPP_VERSION" "$deps/yaml-cpp-src"
        cmake -S "$deps/yaml-cpp-src" -B "$deps/yaml-cpp-src/build" -G "${generator[@]}" \
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$deps/yaml-cpp" -DYAML_BUILD_SHARED_LIBS=ON \
            -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
        cmake --build "$deps/yaml-cpp-src/build" -j "$jobs"
        cmake --install "$deps/yaml-cpp-src/build"
        rm -rf "$deps/yaml-cpp-src"
    fi
    export CMAKE_PREFIX_PATH="$deps/yaml-cpp${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
    local lib="$deps/yaml-cpp/lib"
    [[ -d "$deps/yaml-cpp/lib64" ]] && lib="$deps/yaml-cpp/lib64"
    export LD_LIBRARY_PATH="$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    setup_lines+=("export LD_LIBRARY_PATH=\"$lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}\"")
}

case "$mode" in
system)
    [[ -z "$env_script" ]] || load "$(cd "$(dirname "$env_script")" && pwd)/$(basename "$env_script")"
    ;;
cvmfs)
    [[ -f "$lcg/setup.sh" ]] || fail "LCG view $lcg not found (is /cvmfs mounted?)"
    load "$lcg/setup.sh"
    ;;
source)
    install_system_packages
    install_root
    [[ "$g4" == OFF ]] || install_geant4
    ;;
esac

command -v root-config >/dev/null || fail "ROOT not found (root-config): load it, or --env, or another mode"
info "ROOT $(root-config --version)"
if [[ "$g4" == ON ]]; then
    if command -v geant4-config >/dev/null; then
        info "Geant4 $(geant4-config --version)"
        [[ -n "${CRYHOME:-}" && -d "${CRYHOME:-}" ]] || install_cry
    else
        info "Geant4 not found: without restG4"
        g4=OFF
    fi
fi
if [[ "$garfield" == ON ]]; then
    if [[ "$lcg_garfield" == ON ]]; then
        # Its setup.sh sets neither GARFIELD_INSTALL nor HEED_DATABASE (needed by Heed)
        [[ -d "$lcg/share/Heed/database" ]] || fail "no Garfield++ in the LCG view $lcg (share/Heed/database)"
        export GARFIELD_INSTALL="$lcg" HEED_DATABASE="$lcg/share/Heed/database"
        setup_lines+=("export GARFIELD_INSTALL=\"$GARFIELD_INSTALL\" HEED_DATABASE=\"$HEED_DATABASE\"")
    elif [[ "$mode" != system ]]; then
        install_garfield
    elif [[ -z "${GARFIELD_INSTALL:-}${GARFIELD_HOME:-}" ]]; then
        info "Garfield++ not found (GARFIELD_INSTALL): without restGarfield"
        garfield=OFF
    fi
fi
install_yaml_cpp
# The Garfield++ of GARFIELD_INSTALL given to CMake: CMAKE_PREFIX_PATH (e.g. of an LCG view, with its own
# Garfield++) comes before the HINTS of find_package
garfield_dir=()
if [[ "$garfield" == ON ]]; then
    for dir in "${GARFIELD_INSTALL:-${GARFIELD_HOME:-}}"/lib{64,}/cmake/Garfield; do
        [[ -f "$dir/GarfieldConfig.cmake" ]] && { garfield_dir=(-DGarfield_DIR="$dir"); break; }
    done
    info "Garfield++ ${garfield_dir[*]:-of the environment}"
fi

# --- rest-framework -----------------------------------------------------------------------------------

if [[ -z "$version" ]]; then
    version=$(git ls-remote --tags --refs "$REPOSITORY" 'v*' | sed 's#.*refs/tags/##' \
        | grep -E '^v[0-9]+\.[0-9]+\.[0-9]+$' | sort -V | tail -1)
    version=${version:-master}
fi
src="$prefix/src/rest-framework"
if [[ -d "$src/.git" ]]; then
    [[ -z "$(git -C "$src" status --porcelain --untracked-files=no)" ]] \
        || fail "$src has local changes: commit or remove them, or use another --prefix"
    git -C "$src" fetch -q --tags origin
    if [[ "$version" == master ]]; then
        git -C "$src" checkout -q master && git -C "$src" merge -q --ff-only origin/master
    else
        git -C "$src" checkout -q "$version"
    fi
else
    mkdir -p "$prefix/src"
    git clone -q "$REPOSITORY" "$src"
    git -C "$src" checkout -q "$version"
fi
info "rest-framework $version ($(git -C "$src" rev-parse --short HEAD))"

# Warnings are errors in the CI (REST_WARNINGS_AS_ERRORS); here another compiler may warn of something new
cmake -S "$src" -B "$prefix/build" -G "${generator[@]}" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix/install" -DREST_WARNINGS_AS_ERRORS=OFF \
    -DREST_G4="$g4" -DREST_GARFIELD="$garfield" -DREST_DAQ="$daq" "${garfield_dir[@]}"
cmake --build "$prefix/build" --target install -j "$jobs"

{
    echo "# Environment of rest-for-physics-toolkit (install.sh, $(date +%F), $version, mode $mode)"
    printf '%s\n' "${setup_lines[@]}"
    echo "source \"$prefix/install/thisREST.sh\""
    echo "export REST_SETUP=\"$prefix/setup.sh\""
} >"$prefix/setup.sh"

# Check: the libraries and macros load, the executables run
(
    set +u
    # shellcheck disable=SC1091
    source "$prefix/setup.sh" >/dev/null
    restManager --version
    restRoot -b -q >/dev/null
)
info "installed. Load it in each shell (or in ~/.bashrc):"
echo "    source $prefix/setup.sh"
exit
}
