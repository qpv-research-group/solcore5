Solcore
=======

**Solcore** was born as a modular set of tools, written (almost) entirely in Python 3, to address some of the task we had to solve more often, such as fitting dark IV curves or luminescence decays. With time, however,  it has evolved as a complete semiconductor solver able of modelling the optical and electrical properties of a wide range of solar cells, from quantum well devices to multi-junction solar cells. Some of the features of Solcore are:

    - k•p band structure solver including strain
    - 1D arbitrary potential Schrödinger equation solver
    - Bulk and QW absorption profile calculator
    - Spectral irradiance model and database
    - Multi-junction quantum effciency and IV calculators
    - Coupled Poisson - Drift-Diffusion solver (PDD)

Installation
------------

After downloading Solcore, either using 'git' or as a zip file using one of the links on the right, installing it and become a user should be as easy as writing in the terminal::

    >>> python3 setup.py install

You will be asked to accept the license and then it will install all Solcore dependencies (except a Fortran compiler) and Solcore itself in the Python3 package tree. Using it will be as easy as making::

    >>> import solcore

        Welcome to Solcore - version 4.1.0
        Copyright (c) 2017, Imperial College, London All rights reserved.
        Software released under the GNU Lesser General Public License.


If you want to test first if Solcore will work in your computer, without actually installing it, or if you want to become a developer and therefore you need to have it in a more accessible place, you can test if solcore works with::

    >>> python3 setup.py test

This will install the Solcore dependencies and run a few tests that probe several of the Solcore tools. If it fails, it will indicate which parts failed to work and why, and you could try to solve them.

Things that can go wrong
^^^^^^^^^^^^^^^^^^^^^^^^

There are several things that can go wrong in the above description, specially in Windows.

1. **The tests associated with the Poisson-Drift-Diffusion solver fail**: This is usually the result of not having a Fortran compiler installed in your system, not being correctly configured or having a temperamental F2PY version, the tool - included in numpy - that makes Fotran code accesible from Python. For the first two problems, make sure you actually have a Fortran compiler installed and in the system path. For the latter, it appears in Windows and we have not been able to solve it, yet. Please, let us know if you have a solution.

2. **Some of the dependencies fail to install**: That is rarely the case, as all dependencies are in the main Python repositories. However, there might be issues (again in Windows) with Numpy, Matplotlib and Scipy. These packages need to be compiled and it is often easy to get them as a scientific bundle. You can check Anaconda <https://www.continuum.io/downloads> which provides all these packages together already configured for the correct OS.

Getting started
---------------

After installing Solcore (or even without installing it), there are a few things you might want to do in order to personalise it and start using it.

1. **Create a user configuration file:** This can be done automatically by importing the config_tools. If you do not already have a (hidden) solcore_config.txt file in your home directory, you will be asked if you want to create it::

    >>> import solcore.config_tools

2. **Create a folder with Solcore examples:** This is the fastest way of getting started. The examples will be created in a subfoler called 'solcore/examples'. Simply use the config_tools, again::

    >>> import solcore.config_tools as config

    >>> config.get_solcore_examples('/location/where/you/want/the/examples')

3. **Set the location of a SPICE executable and the SMARTS folder:** You will need to do this eventually in order to use those tools::

    >>> import solcore.config_tools as config

    >>> config.set_location_of_spice('/path/to/the/SPICE/executable')
    >>> config.set_location_of_smarts('/path/to/the/SMARTS/folder')

4. **Open Solcore documentation:** It should contain a description (even minimal) of all Solcore functions, modules and packages. The idea is for it to be a useful tool although it is quite empty, for now. The documentation will open in a webbrowser and you might want to add it to your Bookmarks::

    >>> import solcore.config_tools as config

    >>> config.open_documentation()

5. **Getting specific information about Solcore:** Even though the documentation "should" be more complete, you can get information about any object in python (including any Solcore function, module and package) using the '__doc__' attribute, for example::

    >>> import solcore.config_tools as config

    >>> print(config.get_current_config.__doc__)

    Prints the current Solcore configuration

        :return: None

6. **Python editor:** Learning Python is easy, but some tools make it even easier. That is the case of PyCharm <https://www.jetbrains.com/pycharm/> (the community eddition is free and the other it is too if you are in academia). Selecting an editor is very personal choice, but PyCharm turns out to be quite useful to teach you good coding practices, reviewing your code for errors and, in general, checking that things will work. It will make your life easier. Give it a try. Solcore in its current form is, in part, the result of using PyCharm.

Known issues
------------

We have developed Solcore as part of our ongoing research activities to solve specific challenges, it has (almost) never been a goal in itself. These means that there are parts of Solcore that might not be as polished as they should, that have been just partly implemented or that are only valid under some assumptions (good for us, but maybe not that good for others).

Some of the Solcore issues we are aware off are:

- The poisson-drift-diffusion solver, written in Fortran, only works under Linux and Mac (which are the OS we mostly use). We have never been able to make F2Py and the Fortran compiler work together under Windows. Any help with this is more than welcome!!
- The equations included in the quantum efficiency calculator need a full revision as it is not fuly clear that they include properly the intrinsic region in PIN and NIP solar cells, nor the dependence of the quantum efficiency with voltage.
- Using an external generation profile in the poisson-drift-diffusion solver only works well for dense, homogeneous meshes. The interpolation scheme of the solver does not considers correctly the situation of fast-varying generation.
- The poisson-drift-diffusion solver, the quantum efficiency calculator and the current-voltage calculator all use slightly different input structures, meaning that you can not use the same solar cell object in the three solvers. That's definetly something we need to improve shortly.
- Documentation is incomplete or obscure, in many cases. Again, something to be solved soon.
- The calculator of the generation profile using the TMM module is really, really slow as soon as the structure is slightly complicated or the mesh density is high. We'll need to do something about it sooner than later.


Changelog
---------

**New in Solcore 4.2.0.dev.1**
- Improved the QE module, now being able of providing the IV curve based on the depletion approximation
- Corrected the instalation of examples, making it compatible with Windows

**New in Solcore 4.1.5**
- The way paths in the config file has been changed, so now they are really compatible with Windows
- Modified spice package to work in Linux, Mac OS and Windows systems (tested with ngspice only)

**New in Solcore 4.1.4**
- Increase compatibility with Windows

**New in Solcore 4.1.3**
- Solved a bug in the Spice module

**New in Solcore 4.1.2**
- Solved a problem with the distributable package, that was created incomplete

**New in Solcore 4.1.1**
- Solved bug in the documentation that prevented it to show properly.

**10th March 2017 - A new stable version has been released! - New in Solcore 4.1.0**
- Documentation of greatly improved, though still incomplete.
- License added: Now Solcore is distributed under the GNU Lesser General Public License.
- The old "database" tools have been removed, as they were not adapted to Solcore v4 working principles.


**New in Solcore 4.1.0.dev.7**
- Solve a bug in the way front surface reflexion was calculated in the PDD model: it was using just the real part of the refractive index and not complex refractive index.
- Added the possibility of calculating the optical properties of a solar cell in the PDD model using the TMM formalism or providing an external generation profile as input for the solver under illumination. Note: It only works properly when using very dense, non-dynamic meshes.
- A structure for the PDD solver now can include layers defined as "optics" and "metal" that can be used to calculate the opticla properties but will be ignored y the electrical solver.
- Added the description of most functions of the PDD solver.
- Created new version of the MJ current voltage calculator, as the previous one didn't work in the dark.
- Added tests for the new MJ current voltage calculator
- Added test for the TMM optics calculator of the PDD solver
- Included an option in the TMM calculator (no_back_reflexion) to prevent reflexion from the back surface.
- Added a filter option for the PDD.Illumination object that accept an arbitrary function as filter.

**New in Solcore 4.1.0.dev.6**
- Removed all mention to the Mod CPPB calculator (test, examples, module and data in “absorption_calculator”)
- Moved all material information to the “materials_data” folder.
- SOPRA database moved to the “materials_data” folder. Module modified accordingly. Test added.
- Added the ellipsometry fitting routine to the "data_analysis_tools" package.
- Added Science_Tracker for the reference_spectra. Reference_spectra changed to function.
- Added documentation concerning the SOPRA database and the transfer matrix calculator

**New in Solcore 4.1.0.dev.5**
- Added new dielectric models to the absorption_calculator package
- Minor bug fixing

**New in Solcore 4.1.0.dev.4**
- Added the Modified_CPPB module to the absorption_calculator package.
- Two example files added to AC_examples.
- test_Modified_CPPB added to tests folder.

**New in Solcore 4.1.0.dev.3**
- Added Custom_Colours module to graphing package for adding a splash of colour to individual and multiple plots.
- Added a Custom_Colours example file to the examples folder.

**New in Solcore 4.1.0.dev.2**
- Added SOPRA database of optical constant data to the absorption calculator package. Package now contains a folder
of SOPRA raw data and documentation and a module file containing the sopra_database.

**New in Solcore 4.1.0.dev.1**
- Added a transfer matrix formalism to the absorption calculator package. It can calculate reflexion, transmisison,
absorption, absorption profile and ellipsometry parameters of a optical stack made of solcore materials, theoretical
dielectric models or available n and k data.

**New in Solcore 4.0.3**
- Version number removed from the solcore_config user file.

**New in Solcore 4.0.2**
- Fix a bug related with import statements that where not using the new package names.
- Implemented an automatic MANIFEST.in creator to avoid missing non-python files when creating a distributable package.

**New in Solcore 4.0.1**
- Adds many files (txt, html, etc) that where missing in the previous version for some estrange reason.

**New in Solcore 4.0.0**
- General code re-writing aiming to comply with the PEP 8 -- Style Guide for Python Code (https://www.python.org/dev/peps/pep-0008/).
- The package has been renamed just 'solcore', droping the version number. Importing Solcore now is done by 'import solcore' rather than 'import solcore3', as it is the usual case for python packages which do not include verison number in the name.
- The pluging system has been removed and the code arranged in the more standard package/module configuration. Most old plugings are now directly imported as subpackages of solcore (eg. import solcore.spice rather than import solcore3.plugings.spice).
- The structure of some modules and packages has changed, grouping them in a (hopefully) more consistent and logical order. The following list summarises some of the changes (likely to be incomplete):
    - solcore3.plugings.qm              -> solcore.quantum_mechanics
    - solcore3.plugings.kp              -> solcore.quantum_mechanics
    - solcore3.plugings.IV              -> solcore.analytic_solar_cells
    - solcore3.plugings.analytic_pin_QE -> solcore.analytic_solar_cells
    - solcore3.plugings.spice           -> solcore.spice
    - solcore3.plugings.spectral2       -> solcore.solar_spectrum
    - solcore3.plugings.smarts          -> solcore.solar_spectrum
    - solcore3.plugings.adachi_alpha    -> solcore.absorption_calculator
    - solcore3.PDD                      -> solcore.poisson_drift_diffusion
- Some functions within the above modules and packages have been renamed to have a more descriptive name.
- All the old plugins related with parameters of materials are now included inside the 'parameter system' and not as packages.
- All the old plugins related with description of materials n and k data are now included inside the 'material system' and not as packages.
- A configuration file and a set configuration tools have been created to simplify the addition of new data or personalised behavior.
- SMARTS has been removed from the Solcore directory tree, as it is a third party software and platform dependent. The path to the SMARTS distribution can be set with the configuration tools. SMARTS can be found in http://www.nrel.gov/rredc/smarts/
- Solcore has been arranged to be used with 'setuptools', which simplifies its distribution, testing an maintenance. If everything goes alright, the command 'python3 setup.py install' will install Solcore and all its dependencies (except the Fortran compiler). Tested in MacOS X.
- A set of test have been written to assert the correct behaviour of Solcore, either before performing the installation or if one of the existing packages is modified. They can be run with 'python3 setup.py test'
- The 'poisson_drift_diffusion' solver now can print the output of the calculation to a log file, rather than the terminal.
- The correct temperature dependence has been incorporated to the analytic IV calculator.
- A set of examples have been created to illustrate the use of Solcore. Such examples can be copied to a user-speficied folder, where they can be easily edited.

Contributors
------------

Solcore is the result of many years of developement by many people, trying to put together in a consistent way the tools that the Quantum Photovoltaics Group at Imperial College London <https://www.imperial.ac.uk/quantum-photovoltaics/> needs for its research in solar cells and photovoltaics. The list of contributors (hopefully I am not forgetting anyone!) is:

- Markus Führer
- Daniel Farrel
- Diego Alonso-Álvarez
- Ned Ekins-Daukes
- Tomos Thomas
- Alvin Chan
- Thomas Wilson
