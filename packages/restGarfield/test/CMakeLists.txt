# Unit tests of packages/restGarfield (cmake/RestTests.cmake): executable
# restGarfieldTests, label restGarfield, with the sources of the package (its
# functions, not the executable) and Garfield++
#
# SCAN_SOURCES: Garfield++ prints the Heed database path when it is loaded,
# before main, which the discovery of the tests (--gtest_list_tests) would take
# as a test
rest_add_tests(restGarfield SCAN_SOURCES ${DETECTOR_LIB})
target_sources(restGarfieldTests PRIVATE ${RESTGARFIELD_SOURCES})
target_include_directories(restGarfieldTests
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(restGarfieldTests PRIVATE Garfield::Garfield)
