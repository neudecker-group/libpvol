#
#set(dir "${CMAKE_CURRENT_SOURCE_DIR}")
#
#list(APPEND xhcfftest
#  "${dir}/coffeine.f90"
#  "${dir}/timer.f90"
#  "${dir}/main.F90"
#)
#
#set(xhcfftest ${xhcfftest} PARENT_SCOPE)
#
#list(APPEND pvtest
#  "${dir}/coffeine.f90"
#  "${dir}/lih.f90"
#  "${dir}/supermol.f90"
#  "${dir}/timer.f90"
#  "${dir}/test_pv.F90"
#)
#set(pvtest ${pvtest} PARENT_SCOPE)

# Unit testing
set(
  tests
  "xhcff"
)
set(
  test-srcs
  "coffeine.f90"
  "lih.f90"
  "supermol.f90"
  "timer.f90"
  "main.f90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t)
  list(APPEND test-srcs "test_${t}.F90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)
target_include_directories(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "${CMAKE_BINARY_DIR}/include"  # Path where the .mod files are located
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
  "test-drive::test-drive"
)

# Ensure tester is built after the library
add_dependencies("${PROJECT_NAME}-tester" "${PROJECT_NAME}-lib")

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
endforeach()

