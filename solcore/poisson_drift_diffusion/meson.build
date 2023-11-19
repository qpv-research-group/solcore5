if get_option('with_pdd')
    # NumPy include directory
    
    add_languages('c')
    add_languages('cython')
    add_languages('fortran', native: false)

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

    ff = meson.get_compiler('fortran')
    _global_ff_args = ff.get_supported_arguments(
        '-Wno-unused-variable',
        '-Wno-surprising',
        '-Wno-conversion',
        '-Wno-maybe-uninitialized',
        '-Wno-unused-dummy-argument'
    )
    add_project_arguments(_global_ff_args, language: 'fortran')

    if is_windows and meson.get_compiler('fortran').get_id() == 'gcc'
        # For mingw-w64, link statically against the UCRT.
        gcc_link_args = ['-lucrt', '-static']
        add_project_link_arguments(gcc_link_args, language: ['fortran'])
    endif

    py_dep = py.dependency()

    if is_windows
        incdir_numpy = run_command(py,
        [
            '-c',
            'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
          ],
          check: true
        ).stdout().strip()

    else
        incdir_numpy = run_command(py,
        [
            '-c',
            'import os; import numpy; cwd = os.getcwd(); os.chdir(".."); numpy_path = numpy.get_include(); print(os.path.relpath(numpy_path, cwd));'
        ],
        check : true
        ).stdout().strip()
    endif

    message('numpy DIR', incdir_numpy)

    inc_np = include_directories(incdir_numpy)
    incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
    inc_f2py = include_directories(incdir_f2py)

    message('f2py DIR', incdir_f2py)

    ddModule_source = custom_target('ddModelmodule.c',
                                input : ['DDmodel-current.f95'],
                                output : ['ddModelmodule.c', 'ddModel-f2pywrappers2.f90'],
                                command : [py, '-m', 'numpy.f2py', '-m', 'ddModel', '@INPUT@', '--lower',
                                           '--build-dir', 'solcore/poisson_drift_diffusion/']
                                )

    py.extension_module('ddModel',
                         'DDmodel-current.f95',
                         ddModule_source,
                         incdir_f2py / 'fortranobject.c',
                         include_directories: [incdir_numpy, incdir_f2py],
                         dependencies : py_dep,
                         link_language: 'fortran',
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
  subdir: 'solcore' / 'poisson_drift_diffusion',
)
