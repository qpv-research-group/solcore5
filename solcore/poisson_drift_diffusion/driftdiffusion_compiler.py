import numpy.f2py
import os
import tempfile
import glob
import platform
import shutil

this_dir = os.path.split(__file__)[0]
module_name = 'ddModel'
fotran_source = 'DDmodel-current.f95'

system = platform.system()
if system == 'Windows':
    lib_file = os.path.join(this_dir, 'ddModel.pyd'.format(module_name))
else:
    lib_file = os.path.join(this_dir, 'ddModel.so'.format(module_name))

code_file = os.path.join(this_dir, fotran_source)


def check_ddModel_library_ok(force=False, test=False):
    """ Check if the fortran PDD library exists and if it must be recompiled. If the library file does not exists or it is older than the source file, we re-compile the code. The compilation can be forced, too.

    :param force: If the compilation and building of the library must be done
    :param test: Use for testing purposes
    :return:
    """

    do_it = (not (os.path.exists(lib_file) and os.path.getmtime(lib_file) > os.path.getmtime(code_file))) or force
    if do_it:

        # Optimization flags for the fortran compiler
        # args = '--opt="-O3"'
        code = open(code_file, encoding='utf-8').read()

        # This is the currect directory. We need it in order to go back to it after compiling the code
        pushdir = os.getcwd()

        # We do the compilation in a temporary directory and then we copy the library back to our poisson_drift_diffusion directory
        with tempfile.TemporaryDirectory(prefix="tmp", suffix="_sc3ddModel") as working_directory:

            os.chdir(working_directory)
            result = numpy.f2py.compile(code, modulename=module_name, verbose=False, source_fn='ddModel.f95')

            if result == 1:
                os.chdir(pushdir)
                if test:
                    return False
                raise SyntaxError("Fortran program did not compile")
            else:
                # And now, we copy the library back to the poisson_drift_diffusion folder
                if system == 'Windows':
                    new_lib = glob.glob('*.pyd')[0]
                else:
                    new_lib = glob.glob('*.so')[0]

                shutil.copy2(new_lib, lib_file)
                os.chdir(pushdir)

    return True
