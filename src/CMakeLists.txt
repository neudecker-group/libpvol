
add_subdirectory("surface")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
   "${dir}/xhcff_engrad.f90"
   "${dir}/pv_engrad.f90"
   "${dir}/calculator.f90"
   "${dir}/interface.f90"
   "${dir}/C_interface.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
