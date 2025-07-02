project('solcore',
    version: '5.10.1',
    license: 'GNU LGPL',
    meson_version: '>= 0.64.0',
    default_options : [
        'buildtype=release',
        'c_std=c99',
        'cpp_std=c++14',
        'fortran_std=legacy',
    ]
)

is_windows = host_machine.system() == 'windows'
message('WINDOWS:', is_windows)

py_mod = import('python')
py = py_mod.find_installation('python3', pure: false)
message('py DIR', py.full_path())
message('install DIR', py.get_install_dir())

subdir('solcore')

if get_option('install_test')
    install_subdir('tests', install_dir: py.get_install_dir() / 'solcore')
    install_subdir('examples', install_dir: py.get_install_dir() / 'solcore')
endif

