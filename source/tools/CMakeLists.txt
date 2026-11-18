# ---------------------------------------------------------------------------
# source/tools  -  RestTools (shared library)
# ---------------------------------------------------------------------------

file(GLOB_RECURSE TOOLS_SOURCES src/*.cxx)
file(GLOB_RECURSE TOOLS_HEADERS inc/*.h)

add_library(RestTools SHARED ${TOOLS_SOURCES})

# Version of the code (REST_Version, TRestVersion.h): RestVersionInfo.h written
# at each build by cmake/RestVersion.cmake (only rewritten if it changes)
set(REST_GIT_COMMIT
    ""
    CACHE STRING "Commit of a build without git (e.g. a tarball or an image)")
set(_version_header ${CMAKE_CURRENT_BINARY_DIR}/generated/RestVersionInfo.h)
add_custom_target(
  rest_version
  COMMAND
    ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
    -DOUTPUT=${_version_header} -DVERSION=${PROJECT_VERSION}
    -DGIT_COMMIT=${REST_GIT_COMMIT} -P
    ${CMAKE_SOURCE_DIR}/cmake/RestVersion.cmake
  BYPRODUCTS ${_version_header}
  COMMENT "Version of the code")
add_dependencies(RestTools rest_version)
target_include_directories(RestTools
                           PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)

target_include_directories(
  RestTools PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
                   $<INSTALL_INTERFACE:include>)

target_link_libraries(RestTools PUBLIC yaml-cpp::yaml-cpp
                                       ${ROOT_TARGET_LIBRARIES})
# <filesystem> needs libstdc++fs only before GCC 9 (e.g. RHEL 8)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION
                                            VERSION_LESS 9)
  target_link_libraries(RestTools PUBLIC stdc++fs)
endif()

rest_attach_dictionary(RestTools HEADERS ${TOOLS_HEADERS} INCLUDE_DIRS
                       ${CMAKE_CURRENT_SOURCE_DIR}/inc)

install(
  TARGETS RestTools
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/libRestTools_rdict.pcm
        ${CMAKE_CURRENT_BINARY_DIR}/libRestTools.rootmap
  DESTINATION lib
  OPTIONAL)
install(FILES ${TOOLS_HEADERS} DESTINATION include)
