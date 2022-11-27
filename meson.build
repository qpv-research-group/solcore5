project('solcore',
    'c',
    version: '0.1',
    default_options : [
        'warning_level=0',
    ]
)

add_languages('fortran', 'cython')

py_mod = import('python')
py3 = py_mod.find_installation('python3')
py3_dep = py3.dependency()
message(py3.path())
message(py3.get_install_dir())

install_subdir(
  'solcore',
  install_dir: py3.get_install_dir(pure: false)
)

subdir('solcore/poisson_drift_diffusion')
