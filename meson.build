project('probable', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++14'])

includes = include_directories('src')

shlib = shared_library('probable', 'src/probable.cc',
  include_directories : includes,
  dependencies: dependency('openmp'),
  install : true,
)

# Examples
gallium_oxide = executable('gallium_oxide',
                           'examples/gallium_oxide.cc',
                           include_directories : includes,
                           link_with : shlib)

# Tests
test_rng = executable('rng', 'tests/rng.cc',
                      include_directories : includes,
                      link_with : shlib)
test('rng', test_rng)

# Make this library usable as a Meson subproject.
probable_dep = declare_dependency(
  include_directories: includes,
  link_with : shlib)

# Make this library usable from the system's
# package manager.
install_headers(['src/probable.hh', 'src/vec3.hh'], subdir : 'probable')

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name : 'probable',
  filebase : 'probable',
  description : 'Library for Monte Carlo simulations of electronic transport in solids',
  subdirs : 'probable',
  libraries : shlib,
  version : '0.1',
)
