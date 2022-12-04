if get_option('with_pdd')
    incdir_numpy = run_command(py,
    ['-c', 'import os; import numpy; cwd = os.getcwd(); os.chdir(".."); numpy_path = numpy.get_include(); print(os.path.relpath(numpy_path, cwd));'],
    check : true
    ).stdout().strip()

    incdir_f2py = run_command(py,
        ['-c', 'import os; import numpy.f2py; cwd = os.getcwd(); os.chdir(".."); f2py_path = numpy.f2py.get_include(); print(os.path.relpath(f2py_path, cwd))'],
        check : true
    ).stdout().strip()

    message(incdir_numpy)
    message(incdir_f2py)

    ddModule_source = custom_target('ddModelmodule.c',
                                input : ['DDmodel-current.f95'],
                                output : ['ddModelmodule.c', 'ddModel-f2pywrappers2.f90'],
                                command : [py, '-m', 'numpy.f2py', '-m', 'ddModel', '@INPUT@', '--lower',
                                           '--build-dir', 'solcore/poisson_drift_diffusion']
                                )

    inc_np = include_directories(incdir_numpy, incdir_f2py)

    py.extension_module('ddModel',
                         'DDmodel-current.f95',
                         ddModule_source,
                         incdir_f2py / 'fortranobject.c',
                         include_directories: inc_np,
                         dependencies : py_dep,
                         install_dir: py.get_install_dir() / 'solcore' / 'poisson_drift_diffusion',
                         install : true)
endif


python_sources = [
    '__init__.py',
    'DeviceStructure.py',
    'DriftDiffusionUtilities.py',
    'QWunit.py',
]

py.install_sources(
  python_sources,
  subdir: 'solcore' / 'poisson_drift_diffusion'
)