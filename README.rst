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

    >>> python setup.py install

You will be asked to accept the license and then it will install all Solcore dependencies (except a Fortran compiler) and Solcore itself in the Python3 package tree. Using it will be as easy as making::

    >>> import solcore

        Welcome to Solcore - version 5.0.0
        Copyright (c) 2017, Imperial College, London All rights reserved.
        Software released under the GNU Lesser General Public License.


If you want to test first if Solcore will work in your computer, without actually installing it, or if you want to become a developer and therefore you need to have it in a more accessible place, you can test if Solcore works with::

    >>> python setup.py test

This will install the Solcore dependencies and run a few tests that probe several of the Solcore tools. If it fails, it will indicate which parts failed to work and why, and you could try to solve them. At the moment, this only cover some of Solcore's functionality, but it will be expanded with time.

Another thing that you might want to do before installing Solcore 5 is compiling the Poisson-drfit-diffusion solver. Assuming there is a Fortran compiler correctly configured to work with F2Py, compiling the library should be as easy as::

    >>> python setup.py build_pdd

This can also be done afterwards using the *config_tools* but you might need admin privileges, depending on where is the Python packages tree.

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

2. **Create a folder with Solcore examples:** This is the fastest way of getting started. The examples will be created in a subfolder called 'solcore/examples'. These examples include all the scripts used in the main Solcore paper (submitted to Computer Physics Communications, preprint in https://arxiv.org/abs/1709.06741 ) . Simply use the config_tools, again::

    >>> import solcore.config_tools as config

    >>> config.get_solcore_examples('/location/where/you/want/the/examples')

3. **Set the location of a SPICE executable and the SMARTS folder:** You will need to do this eventually in order to use those tools::

    >>> import solcore.config_tools as config

    >>> config.set_location_of_spice('/path/to/the/SPICE/executable')
    >>> config.set_location_of_smarts('/path/to/the/SMARTS/folder')

4. **Open Solcore documentation:** It should contain a description (even minimal) of all Solcore functions, modules and packages. The idea is for it to be a useful tool although it is quite empty, for now. The documentation will open in a web browser and you might want to add it to your Bookmarks::

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

- The poisson-drift-diffusion solver, written in Fortran, has been tested only under Linux and Mac. We have never been successful in making F2Py and the Fortran compiler work together under Windows, although they are supposed to work well. Any help with this is more than welcome!!
- Documentation is incomplete or obscure, in many cases. Again, something to be solved soon.
- The calculator of the generation profile using the TMM module is really, really slow as soon as the structure is slightly complicated or the mesh density is high. We'll need to do something about it sooner than later.

Contributors
------------

Solcore is the result of many years of developement by many people, trying to put together in a consistent way the tools that the Quantum Photovoltaics Group at Imperial College London <https://www.imperial.ac.uk/quantum-photovoltaics/> needs for its research in solar cells and photovoltaics. The list of contributors (hopefully I am not forgetting anyone!) is:

- Diego Alonso-Álvarez
- Thomas Wilson
- Phoebe Pearce
- Markus Führer
- Daniel Farrel
- Tomos Thomas
- Alvin Chan
- Ned Ekins-Daukes

File/directory description
--------------------------

**solcore/**

This is the base directory of Solcore, containing numerous modules that the other, more complex packages depend upon.

- *config_tools.py* and *solcore_config.txt* contain the solcore configuration and tools to modify it. These include the route to the material databases (internal to Solcore or user defined), the SPICE and SMARTS binaries, accessing the documentation or the examples.
- *source_managed_class.py* and *singleton.py* contain classes to manage the sources of Solcore, from the materials to the unit conversion system.
- *smooth.py* and *interpolate.py* contain tools to smooth and interpolate data. Most of their uses in Solcore could be replaced by Numpy built-in functions, but they are maintained for backwards compatibility.
- *crystals.py* calculates the k-points of the Brillouin zone in a given direction.
- *strain_balancing.py* contains functions to calculate the critical thickness of a given strain material and to optimise the thickness and/or composition of a quantum well strcuture in order to satisfy the strain balance condition.
- *constants.py* just contains a list of universal constants.
- *science_tracker.py* contains the functions used along Solcore to track the references to the papers and scientific works presenting the theory Solcore is based on.
- *state.py*, *structure.py* and *solar_cell.py* define data containers: *Layer*, *Structure*, *Junction*, *TunnelJunction* and *SolarCell*. Together with *State*, which is a custom version of a dictionary, they form the building blocks of Solcore.
- *solar_cell_solver.py* is a common interface to all the solar cell solvers, including electrical and electrical solvers.

**solcore/material_data/**

Contains material properties, including refractive indeces, bandstructure and transport properties, among others. The data is organised in subfolders, depending on where the data is coming from. The *mobility.py* module contains the functions needed to calculate the mobility of several semiconductor materials.

**solcore/material_system/**

Contains several modules to deal with materials. The most important one is *material_system.py* which contains Solcore's *material* class and subclasses. The *critical_point_interpolate.py* and *critical_point_picker.py* are used to interpolates the n and k data of an alloy based on the known values at certain specific compositions. This interpolation is done in a smart way based on certain critical points (usually the Adachi critical points) and then filling the gaps in between.

**solcore/parameter_system/**

Contains just one module, *parameter_system.py*, with the tools necessary to access the materials data and, in case of derived parameters, calculating them. This module is used extensively by the *material_system*, but it can be used on its own to access individual properties, without creating a material.

**solcore/units_system/**

Contains the module *units_system.py* with functions to convert magnitudes between different units in a high level fashion.

**solcore/light_source/**

This folder contains the tools for creating and managing an illumination spectra, typically, a solar spectrum.

- The core functions are contained in the module *light_source.py*.
- *spectral2.py* contains a python implementation of the irradiance model developed by the NREL ( http://rredc.nrel.gov/solar/models/spectral/SPCTRAL2/ )
- *smarts.py* is an interface to the SMARTS irradiance model, which needs to be installed separately from http://www.nrel.gov/rredc/smarts/ .

**solcore/absorption_calculator/**

Contains several tools to calculate the optical properties (dielectric constants, n, k absorption coeficients) of bulk materials and quantum wells as well as tools to calculate the propagation of light in stacks of structures.

- *absorption_QW.py* calculates the absorption coefficient of quantum wells from the optput of the Schrödinger calculator.
- *adachi_alpha.py* uses Adachi parametrization to calculate the absorption coefficient of some bulk materials. The CPPB subfolder contains an advanced version of this method.
- *dielectric_constant_models.py* contains a large variety of oscillator models, including Lorentz, Cauchy or Drude, which can be combined to create the total dielectric function of a material.
- *sopra.py* contains the functions used to access the SOPRA database, in the **solcore/material_data/** folder.
- *kramers_kronig.py* is used to calculate Kramers-Kronig compatible dielectric constants.
- *transfer_matrix.py* contains the class OptiStack, used across Solcore to represent optical structures, and functions to interface with the TMM calculator, the *tmm* Python package created by Steven Byrnes and included in the PyPi repository.
- Finally, *rigorous_coupled_wave.py* is used to interface with the RCWA solver S4, from the Stanford University, which needs to be installed separately.

**solcore/optics/**

Includes several modules to calculate the reflected, transmitted and absorbed light in a SolarCell object, using the Beer-Lamber law (*beer_lambert.py*), the transfer matrix method (*tmm.py*) and the rigorous coupled-wave analysis (*rcwa.py*). The later two cases consists mostly on an interface to the OptiStack and the solvers described in **solcore/absorption_calculator/**.

**solcore/analytic_solar_cells/**

This folder contains the tools needed to calculate the IV and QE characteristics of solar cells using analytical or semianalitical models

- *depletion_aproximation.py* for the depletion approximation
- *detailed_balanced.py* for the detailed balance model
- *diode_equation.py* for the 2-diode model
- The module *IV.py* contains the multi-junction solar cell calculator, combining the IV characteristics of the individual junctions with or without radiative coupling, as well as the legacy IV calculator, kept for backwards compatibility purposes.
- The module *QE.py* as also been maintained just for backwards compatibility purposes, as now the QE is calculated in the above modules, depending on the chosen solar cell model.

**solcore/quantum_mechanics/**

This folder includes all the tools related to the quantum properties of materials and structures.

- The *kp_bulk.py* module solves the 8-band kp Hamiltonian for bulk materials under strain.
- *kp_QW.py* solves the 4-band and 6-band kp Hamiltonian for quantum wells (not fully implemented, yet).
- *heterostructure_alignment.py* uses the band offsets to align the conduction and valence bands of a heterostructure before staring any quantum calculation.
- *strain.py* calculates the strain in a heterostructure and shifts the band edges accordingly.
- *potential_utilities.py* contains the time independent, 1D Schrödinger solver.
- *structure_utilities.py* uses the kp and strain calculators mentioned above to calculate the bands and efective mass profile of a heterostructure.
- Finally, *high_level_kp_qw.py* provides a common interface for the above solvers.

**solcore/poisson_drift_diffusion/**

This folders includes the Poisson-Drift diffusion (PDD) solver.

- *DDmodel-current.f95* is the Fortran code, which needs ot be compiled in a library accessible by Python. This compilation is done by *driftdiffusion_compiler.py* using F2Py.
- *DeviceStructure.py* extracts from the materials making the solar cells all the properties needed by the solver, and also include tools for saving and loading structures from external files.
- *QWunit.py* solves the properties of a QW structure calling the relevant functions in the **solcore/quantum_mechanics/** module and creates and effective bulk-like stack of materials, usable by the PDD solver.
- Finally, the *DriftDiffusionUtilities.py* contains all the functions to interface with the fortran library and solve the PDD equations under equilibrium, short circuit, calculate the IV curves or the QE.

**solcore/spice/**

The module *spice.py* contains the tools necessary to interface with the electrical solver SPICE, which needs to be installed independently. The other modules in these folder depend on this one. The module *quasi_3D_solver.py* has the tools to solve calculate the IV curve of a solar cell by modelling it as a 3D network of electrical components. The module *pv_module_solver.py* has the tools for calculating the IV curve of a solar module, with many solar cells connected in series and a potential random dispersion of properties.

**solcore/data_analysis_tools/**

Contains modules designed to fit experimental data. For now, it only has *ellipsometry_analysis.py*, to fit ellipsometry data and get the dielectric functions of a stack of materials.

**solcore/graphing/**

Contains several modules related to plotting data, intended to help creating default graphs for complex data, for example a quantum well with energy levels and wavefunctions.

**solcore/examples/**

This folder and subfolders contain example scripts illustrating Solcore's functionality. Most of them reproduce the figures in the main Solcore paper (submitted to Computer Physics Communications, preprint in https://arxiv.org/abs/1709.06741 ), but there are other examples expanding the rest of the capabilities.

**solcore/tests/**

Contains test scripts to be run using "nosetests", allowing to probe that the main capabilities of Solcore work as expected and that there is not
regression after adding new functionality.

**solcore/documentation/**

Contains Solcore's documentation, created with Sphinx. It has several subfolders needed by Sphinx.