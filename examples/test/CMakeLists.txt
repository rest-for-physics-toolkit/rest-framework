# Pipeline tests (cmake/RestTests.cmake): executable pipelinesTests, label
# pipelines. Each test runs restManager with the YAML of the examples on real
# data of rest-common-data (downloaded by URL, cached in the work directory) and
# compares the analysis tree of the outputs with a reference of reference/
# (REST_UPDATE_REFERENCE=1 writes the references instead)
rest_add_tests(pipelines RestTools)
target_compile_definitions(
  pipelinesTests PRIVATE REPOSITORY_DIR="${CMAKE_SOURCE_DIR}"
                         RESTMANAGER_EXECUTABLE="$<TARGET_FILE:restManager>")
add_dependencies(pipelinesTests restManager)
