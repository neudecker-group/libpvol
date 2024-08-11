# Specify the minimum version of CMake
cmake_minimum_required(VERSION 3.14)


set(
  tests
  "xhcff-C-API"
)
set(
  test-srcs-c
  "main.c"
)

add_executable(
  "${PROJECT_NAME}-c-tester"
  "${test-srcs-c}"
)
target_include_directories(
  "${PROJECT_NAME}-c-tester"
  PRIVATE
  "${CMAKE_BINARY_DIR}/include"  # Path where the .mod files are located
)
target_link_libraries(
  "${PROJECT_NAME}-c-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
)


foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-c-tester" "${t}")
endforeach()

