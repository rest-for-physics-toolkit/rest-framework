# Unit tests of packages/restGarfield (cmake/RestTests.cmake): executable
# restGarfieldTests, label restGarfield, with the sources of the package (its
# functions, not the executable) and Garfield++
rest_add_tests(restGarfield ${DETECTOR_LIB})
target_sources(restGarfieldTests PRIVATE ${RESTGARFIELD_SOURCES})
target_include_directories(restGarfieldTests
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(restGarfieldTests PRIVATE Garfield::Garfield)
