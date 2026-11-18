# Unit tests of packages/restDAQ (cmake/RestTests.cmake): executable
# restDAQTests, label restDAQ, with the sources of the package (no electronics
# needed: frames built in memory, configurations written in the work directory)
rest_add_tests(restDAQ ${RAW_LIB})
target_sources(restDAQTests PRIVATE ${RESTDAQ_SOURCES})
target_include_directories(
  restDAQTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include
                       ${READLINE_INCLUDE_DIR})
target_link_libraries(
  restDAQTests PRIVATE prometheus-cpp::core prometheus-cpp::pull
                       ${READLINE_LIBRARY} Threads::Threads)
