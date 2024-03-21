
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND xhcfftest
  "${dir}/coffeine.f90"
  "${dir}/timer.f90"
  "${dir}/main.F90"
)

set(xhcfftest ${xhcfftest} PARENT_SCOPE)

list(APPEND pvtest
  "${dir}/coffeine.f90"
  "${dir}/lih.f90"
  "${dir}/supermol.f90"
  "${dir}/timer.f90"
  "${dir}/test_pv.F90"
)
set(pvtest ${pvtest} PARENT_SCOPE)
