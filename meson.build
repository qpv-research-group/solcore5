project('solcore',
    'c',
    version: '5.8.1',
    license: 'GNU LGPL',
    meson_version: '>= 0.64.0',
    default_options : [
        'buildtype=debugoptimized',
        'c_std=c99',
        'cpp_std=c++14',
    ]
)

cc = meson.get_compiler('c')
# TODO: the below -Wno flags are all needed to silence warnings in
# f2py-generated code. This should be fixed in f2py itself.
_global_c_args = cc.get_supported_arguments(
  '-Wno-unused-but-set-variable',
  '-Wno-unused-function',
  '-Wno-conversion',
  '-Wno-misleading-indentation',
  '-Wno-incompatible-pointer-types',
)
add_project_arguments(_global_c_args, language : 'c')

add_languages('fortran', native: false)
add_languages('cython')

is_windows = host_machine.system() == 'windows'

py_mod = import('python')
py = py_mod.find_installation('python3', pure: false)
py_dep = py.dependency()
message(py.full_path())
message(py.get_install_dir())

subdir('solcore')

if get_option('install_test')
    install_subdir('tests', install_dir: py.get_install_dir() / 'solcore')
    install_subdir('examples', install_dir: py.get_install_dir() / 'solcore')
endif

