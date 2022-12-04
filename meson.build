project('solcore',
    'c',
    version: '0.1',
    meson_version: '>= 0.64.0',
    default_options : [
        'buildtype=debugoptimized',
        'c_std=c99',
        'cpp_std=c++14',
    ]
)

add_languages('fortran', 'cython')

py_mod = import('python')
py = py_mod.find_installation('python3', pure: false)
py_dep = py.dependency()
message(py.full_path())
message(py.get_install_dir())

subdir('solcore')
install_subdir('tests', install_dir: py.get_install_dir() / 'solcore')

