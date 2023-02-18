project('solcore',
    'c',
    version: '5.9.1',
    license: 'GNU LGPL',
    meson_version: '>= 0.64.0',
    default_options : [
        'buildtype=release',
        'c_std=c99',
        'cpp_std=c++14',
        'fortran_std=legacy',
    ]
)

cc = meson.get_compiler('c')
# TODO: the below -Wno flags are all needed to silence warnings in
# f2py-generated code. This should be fixed in f2py itself.
_global_c_args = cc.get_supported_arguments(
  '-Wno-unused-but-set-variable',
  '-Wno-unused-function',
  '-Wno-unused-variable',
  '-Wno-surprising',
  '-Wno-conversion',
  '-Wno-misleading-indentation',
  '-Wno-incompatible-pointer-types',
)
add_project_arguments(_global_c_args, language : 'c')

add_languages('cython')
add_languages('fortran', native: false)

cross = meson.is_cross_build()
message('CROSS:', cross)

ff = meson.get_compiler('fortran')
_global_ff_args = ff.get_supported_arguments(
  '-Wno-unused-variable',
  '-Wno-surprising',
  '-Wno-conversion',
  '-Wno-maybe-uninitialized',
  '-Wno-unused-dummy-argument'
)
add_project_arguments(_global_ff_args, language: 'fortran')

is_windows = host_machine.system() == 'windows'
message('WINDOWS:', is_windows)

if is_windows and meson.get_compiler('fortran').get_id() == 'gcc'
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  add_project_link_arguments(gcc_link_args, language: ['fortran'])
endif

py_mod = import('python')
py = py_mod.find_installation('python3', pure: false)
py_dep = py.dependency()
message('py DIR', py.full_path())
message('install DIR', py.get_install_dir())

subdir('solcore')

if get_option('install_test')
    install_subdir('tests', install_dir: py.get_install_dir() / 'solcore')
    install_subdir('examples', install_dir: py.get_install_dir() / 'solcore')
endif

