project(
  'xhcff',
  'fortran', 'c',
  version: '0.0.1',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.63',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=static',
    'c_link_args=-static',
    'fortran_link_args=-static',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
exe_deps = []
subdir('config')

# Documentation
#subdir('docs')

# Collect source of the project
prog = []
srcs = []
subdir('src')


# Create library target
xhcfflib_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: exe_deps,
  include_directories: include_directories('include'),
)

# Export as dependency
xhcfflib_inc = [include_directories('include'), xhcfflib_lib.private_dir_include()]
xhcfflib_dep = declare_dependency(
  link_with: xhcfflib_lib,
  include_directories: xhcfflib_inc,
  dependencies: exe_deps,
  variables: {'includedir': meson.current_source_dir() / 'include'},
)


# Create executable target
if get_option('build_exe')
subdir('testprog')
xhcfflib_exe = executable(
  meson.project_name(),
  sources: prog,
  dependencies: xhcfflib_dep,
  install: install,
  link_language: 'fortran',
)
endif


# add the testsuite separate meson.build 
#subdir('testsuite') #has to be filled with tests, not availabel yet
