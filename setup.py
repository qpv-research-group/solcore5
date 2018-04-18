from setuptools import find_packages
from numpy.distutils.core import setup, Extension

import os
import sys
from configparser import ConfigParser


def gen_data_files(*dirs):
    """ Creates the list of files (not necessarily python files) that need to be installed together with the rest
    of stuff """
    results = []
    exclude = ['.DS_Store', '__pycache__', 'egg', '.git']
    for src_dir in dirs:
        for root, dirs, files in os.walk(src_dir):
            if not any(sub in root for sub in exclude):
                results.append((root, [os.path.join(root, f) for f in files if f not in exclude]))
    return results


here = os.path.abspath(os.path.dirname(__file__))
default_config = os.path.join(here, 'solcore', 'solcore_config.txt')
config = ConfigParser()
config.read([default_config])

# We give the option of not compiling - and installing - the extension modules
if '--no_pdd' in sys.argv:
    ext = []
    sys.argv.remove('--no_pdd')
else:
    sources = os.path.join('solcore', 'poisson_drift_diffusion', 'DDmodel-current.f95')
    ext = [Extension(name='solcore.poisson_drift_diffusion.ddModel',
                     sources=[sources],
                     f2py_options=['--quiet'])]

# Option for updating the manifest
if 'update_manifest' in sys.argv:
    # Update the MANIFEST.in file with all the data from folder and subfolders within solcore
    include = 'solcore'
    exclude = ['__pycache__', 'egg']
    with open('MANIFEST.in', 'w', encoding='utf-8') as f:
        for root, dir, files in os.walk("."):
            if not any(sub in root for sub in exclude) and root[2:9] == include:
                try:
                    files.remove('.DS_Store')
                except ValueError:
                    pass
                for file in files:
                    include_line = 'include ' + os.path.join(root[2:], file) + '\n'
                    f.write(include_line)

    sys.exit()

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='solcore',
    version=config.get('Configuration', 'version'),
    description='Python-based solar cell simulator',
    long_description=long_description,
    url='https://github.com/dalonsoa/solcore5',
    download_url='https://github.com/dalonsoa/solcore5/archive/v{}.tar.gz'.format(config.get('Configuration', 'version')),
    project_urls={
        'Documentation': 'https://dalonsoa.github.io/solcore5/html/index.html',
        'Solcore research paper': 'https://doi.org/10.1007/s10825-018-1171-3'},
    author='The Quantum Photovoltaics Group',
    author_email='d.alonso-alvarez@imperial.ac.uk',
    license='GNU LGPL',
    python_requires='>=3.4',
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
    ext_modules=ext,
    keywords='photovoltaics modelling physics',
    packages=find_packages(exclude=[]),
    package_data={'': ['*.*']},
    data_files=gen_data_files('solcore'),
    install_requires=['numpy', 'matplotlib', 'scipy', 'Sphinx', 'tmm', 'natsort', 'regex', 'cycler'],
    include_package_data=True,
    test_suite='nose.collector',
    tests_require=['nose', 'numpy', 'matplotlib', 'scipy', 'Sphinx', 'tmm', 'natsort', 'regex', 'cycler'],
)
