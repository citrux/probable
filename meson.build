project('probable', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++14'])

shlib = shared_library('probable', 'probable.cc',
  dependencies: dependency('openmp'),
  install : true,
)

# Examples
gallium_oxide = executable('gallium_oxide',
                           'examples/gallium_oxide.cc',
                           link_with : shlib)

# Tests
test_rng = executable('rng', 'tests/rng.cc', link_with : shlib)
test('rng', test_rng)

# Make this library usable as a Meson subproject.
probable_dep = declare_dependency(
  include_directories: include_directories('.'),
  link_with : shlib)

# Make this library usable from the system's
# package manager.
install_headers(['probable.hh', 'vec3.hh'], subdir : 'probable')

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name : 'probable',
  filebase : 'probable',
  description : 'Library for Monte Carlo simulations of electronic transport in solids',
  subdirs : 'probable',
  libraries : shlib,
  version : '0.1',
)
