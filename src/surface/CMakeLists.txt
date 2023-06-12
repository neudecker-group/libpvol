set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
 "${dir}/vdwrad.f90"
 "${dir}/lebedev.f90"
 "${dir}/bisectsearch.f90"
 "${dir}/sasa.f90"
 "${dir}/math_wrapper.f90"
 "${dir}/surface.f90"
 "${dir}/tesspoints.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
