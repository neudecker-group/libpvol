project(
  'pvol',
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

# Collect source of the project
prog = []
srcs = []
subdir('src')

# Create library target
libpvol_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: exe_deps,
  include_directories: include_directories('include'),
)

# Export as dependency
libpvol_inc = [include_directories('include'), libpvol_lib.private_dir_include()]
libpvol_dep = declare_dependency(
  link_with: libpvol_lib,
  include_directories: libpvol_inc,
  dependencies: exe_deps,
  variables: {'includedir': meson.current_source_dir() / 'include'},
)


# Create executable target
if get_option('build_exe')
app = []
subdir('app')
libpvol_exe = executable(
  'pvol',
  sources: app,
  dependencies: libpvol_dep,
  install: install,
  link_language: 'fortran',
)
endif


# add the testsuite separate meson.build 
if get_option('WITH_TESTS')
subdir('test')
endif
