
#add_subdirectory("surface")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
   "${dir}/xhcff_engrad.f90"
   "${dir}/xhcff_interface.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
