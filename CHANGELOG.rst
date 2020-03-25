**New in Solcore 5.6.0**
- Add automatic deployment to PyPi
- Add adding custom materials and access to the refractive index database
- Add an alternative formulation for the detail balance calculation

**New in Solcore 5.6.0.dev1**
- Parameters "relative_dielectric_constant" and "dielectric_constant" have been renamed "relative_permittivity" and "permittivity" in the parameter system and solar cell solvers. Changed also in the examples.
- Light holes density of states are now considered in the Depletion Approximation solver by using the ni calculated in the materials system
- The Depletion Approximation solver now checks that the junction is an homojunction
- MORE TO BE ADDED (RAY TRACING, REFRACTIVEINDEX, ETC)
- Added Sunglass, Solcore's graphic user interface

**New in Solcore 5.5.2**
- Updated installation instructions to avoid an error that ocurs when PiP passes the --with_pdd option to all Solcore dependencies.

**New in Solcore 5.5.1**
- Changed the way QWs are solved before using them in a solar cell to avoid problems with not finding the n and k data

**New in Solcore 5.5.0**
- Updated documentation for solar cell solvers
- The *position* option is now set consistently for all optical solvers
- License badge added to the README file
- The Fortran compiler methods have been eliminated, since this is done by distutils, now.
- Eliminated the auto-meshing during the IV calculation of the PDD solver. It was causing inconsistencies.
- Changed name of examples files to be more meaningful
- Corrected an error in light_source.smarts
- Eliminated legacy functions from config_tools
- Added default "pn=True" value to the tunnel junctions
- Added more examples to the Examples page

**New in Solcore 5.4.0**
- Corrected several file names to avoid spaces: all have been replaced by underscores
- Added exe=1 flag in setup.cfg to ensure tests are run, even if they are labelled as executable files.
- Several files adapted to be able to create the documentation in Read the Docs
- Now the default installation does NOT install the PDD solver and, in order to install it, the option --with_pdd must be given.
- Added Examples section in the documentation

**New in Solcore 5.3.0**
- Solcore has been added to Pypi and can now be installed with ```pip install solcore```
- New installer including the compilation of extension modules
- New documentation about installing Solcore, specially under Windows

**New in Solcore 5.2.1**
- Corrected errors in tests and examples.

**New in Solcore 5.2.0**
- Improvements in the quantum solver
- Added plot generation for LDOS
- Corrected an error in the parameters for InP
- Implemented parametric tunnel junction model
- Added support for external tunnel junction IV characteristics
- Improved external optics option
- Added optical properties for GaAsP, GaInP and AlInP, although with just 2 or three-point interpolation
- Started a simple tutorial
- Polished the PDD solver and separated the preparation of the solar cell from the optical solvers
- Material system updated with new properties, and other minor changes
- Radiative recombination coefficient for PDD is calculated from the absorption coefficient
- Improvements in the tunnel junction and the plotting of wavefunctions
- The TMM solver has been vectorized and built directly into Solcore.

**New in Solcore 5.0.0**
- A new SolarCell class has been implemented to create solar cell structures, including not only the junctions but also optical layers.
- All solar cell solvers have been re-written in order to use a common interface - the SolarCell class.
- The calculation of the optical properties of a solar cell has been externalized. Added a transfer matrix algorithm and an interface for the S4 RCWA package.
- Materials in the SOPRA database can now ve imported as normal Solcore materials.
- Added the Modified_CPPB module to the absorption_calculator package.
- The Solcore materials now use, automatically, the mobility calculator if available for that material.

**New in Solcore 4.2.0.dev.3**
- Mobility calculator and data moved to the material_data package

**New in Solcore 4.2.0.dev.2**
- Created the class LightSource as a common interface to all types of light sources.
- solar_spectrum package renamed light_source
- Other Solcore packages adapted to use LightSource, while keeping as much backwards compatibility as possible.

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
