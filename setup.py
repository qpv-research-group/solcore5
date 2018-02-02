from setuptools import setup, find_packages
import os
import sys
from configparser import ConfigParser

here = os.path.abspath(os.path.dirname(__file__))

default_config = os.path.join(here, 'solcore', 'solcore_config.txt')
config = ConfigParser()
config.read([default_config])

if 'update_manifest' in sys.argv:
    # Update the MANIFEST.in file with all the data from folder and subfolders within solcore
    include = 'solcore'
    exclude = ['__pycache__', 'egg']
    with open('MANIFEST.in', 'w', encoding='utf-8') as f:
        for root, dir, files in os.walk("."):
            if not any(sub in root for sub in exclude) and root[2:9] == include:
                include_line = 'include ' + os.path.join(root[2:], '*') + '\n'
                f.write(include_line)

    sys.exit()

if 'build_pdd' in sys.argv:
    # Compiles the fortran-based poisson-drift-diffusion solver, wrapping the resulting library using F2Py to be accessible from Python.
    from solcore.poisson_drift_diffusion.driftdiffusion_compiler import check_ddModel_library_ok

    check_ddModel_library_ok(force=True)
    sys.exit()

if 'install' in sys.argv:
    with open(os.path.join(here, 'LICENSE'), encoding='utf-8') as f:
        print(f.read().encode("utf-8"))

    answer = input('Do you agree with these conditions (y/n)? ')
    answer = answer.lower()

    if answer != 'y':
        print('Sorry, but you need to accept the license in order to use Solcore.\n')
        sys.exit()

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='solcore',
    version=config.get('Configuration', 'version'),
    description='Python-based solar cell simulator',
    long_description=long_description,
    url='https://www.imperial.ac.uk/quantum-photovoltaics',
    author='The Quantum Photovoltaics Group',
    author_email='d.alonso-alvarez@imperial.ac.uk',
    license='GNU LGPL',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved',
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords='photovoltaics modelling physics',
    packages=find_packages(exclude=[]),
    install_requires=['numpy', 'matplotlib', 'scipy', 'Sphinx', 'tmm', 'natsort', 'regex', 'cycler'],
    include_package_data=True,
    test_suite='nose.collector',
    tests_require=['nose', 'numpy', 'matplotlib', 'scipy', 'Sphinx', 'tmm', 'natsort', 'regex', 'cycler'],

)
