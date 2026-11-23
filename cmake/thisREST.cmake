# ---------------------------------------------------------------------------
# Generates thisREST.sh (and thisREST.csh) in the install prefix. Uses
# configure_file() so no escaping gymnastics are needed.
# ---------------------------------------------------------------------------

# Locate thisroot.sh for the ROOT version used at compile time. ROOT_BINDIR
# comes from find_package(ROOT), so root-config does not need to be in PATH.
# Distribution packages (e.g. Fedora) may not ship thisroot.sh: the generated
# script only sources it if it exists.
set(thisROOT "${ROOT_BINDIR}/thisroot.sh")

# CRY of restG4 (found at configure time, REST_CRY): the scripts set CRYHOME and
# CRYDATAPATH as the CRY setup script does, unless already set, so the CRY setup
# does not need to be sourced (a new configuration finds CRY from CRYHOME)
set(thisCRYHOME "")
set(thisCRYDATA "")
if(REST_G4
   AND REST_CRY
   AND CRY_LIBRARY
   AND CRY_DATA_DIR)
  get_filename_component(_cryLibDir "${CRY_LIBRARY}" DIRECTORY)
  get_filename_component(thisCRYHOME "${_cryLibDir}" DIRECTORY)
  set(thisCRYDATA "${CRY_DATA_DIR}")
endif()

# Generate the scripts into the build tree; install() copies them
configure_file("${CMAKE_CURRENT_LIST_DIR}/thisREST.sh.in"
               "${CMAKE_BINARY_DIR}/thisREST.sh" @ONLY)

configure_file("${CMAKE_CURRENT_LIST_DIR}/thisREST.csh.in"
               "${CMAKE_BINARY_DIR}/thisREST.csh" @ONLY)

install(
  FILES "${CMAKE_BINARY_DIR}/thisREST.sh" "${CMAKE_BINARY_DIR}/thisREST.csh"
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
    WORLD_READ
    WORLD_EXECUTE)
