project('probable', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++14'])

shlib = shared_library('probable', 'probable.cc',
  dependencies: dependency('openmp'),
  install : true,
)

test_exe = executable('probable', 'probable_test.cc',
  link_with : shlib)
test('probable', test_exe)

# Make this library usable as a Meson subproject.
probable_dep = declare_dependency(
  include_directories: include_directories('.'),
  link_with : shlib)

# Make this library usable from the system's
# package manager.
install_headers('probable.hh', subdir : 'probable')

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name : 'probable',
  filebase : 'probable',
  description : 'Meson sample project.',
  subdirs : 'probable',
  libraries : shlib,
  version : '0.1',
)
