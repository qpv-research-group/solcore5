Solcore
=======

**Solcore** was born as a modular set of tools, written (almost)
entirely in Python 3, to address some of the task we had to solve more
often, such as fitting dark IV curves or luminescence decays. With time,
however, it has evolved as a complete semiconductor solver able of
modelling the optical and electrical properties of a wide range of solar
cells, from quantum well devices to multi-junction solar cells. Some of its features are summarised in the following figure: 

![infographics](docs/SphinxDocs/source/_static/infographics.png)

Please, visit [Solcore\'s
Documentation](http://dalonsoa.github.io/solcore5) and check the
Examples of the examples folder.

Contributors
------------

Solcore is the result of many years of developement by many people,
trying to put together in a consistent way the tools that the [Quantum
Photovoltaics Group at Imperial College London](https://www.imperial.ac.uk/quantum-photovoltaics/) needs for its
research in solar cells and photovoltaics. The list of contributors
(hopefully I am not forgetting anyone!) is:

-   Diego Alonso-Álvarez
-   Thomas Wilson
-   Phoebe Pearce
-   Markus Führer
-   Daniel Farrel
-   Tomos Thomas
-   Alvin Chan
-   Ned Ekins-Daukes

File/directory description
--------------------------

**solcore/**

This is the base directory of Solcore, containing numerous modules that
the other, more complex packages depend upon.

-   *config\_tools.py* and *solcore\_config.txt* contain the solcore
    configuration and tools to modify it. These include the route to the
    material databases (internal to Solcore or user defined), the SPICE
    and SMARTS binaries, accessing the documentation or the examples.
-   *source\_managed\_class.py* and *singleton.py* contain classes to
    manage the sources of Solcore, from the materials to the unit
    conversion system.
-   *smooth.py* and *interpolate.py* contain tools to smooth and
    interpolate data. Most of their uses in Solcore could be replaced by
    Numpy built-in functions, but they are maintained for backwards
    compatibility.
-   *crystals.py* calculates the k-points of the Brillouin zone in a
    given direction.
-   *strain\_balancing.py* contains functions to calculate the critical
    thickness of a given strain material and to optimise the thickness
    and/or composition of a quantum well strcuture in order to satisfy
    the strain balance condition.
-   *constants.py* just contains a list of universal constants.
-   *science\_tracker.py* contains the functions used along Solcore to
    track the references to the papers and scientific works presenting
    the theory Solcore is based on.
-   *state.py*, *structure.py* and *solar\_cell.py* define data
    containers: *Layer*, *Structure*, *Junction*, *TunnelJunction* and
    *SolarCell*. Together with *State*, which is a custom version of a
    dictionary, they form the building blocks of Solcore.
-   *solar\_cell\_solver.py* is a common interface to all the solar cell
    solvers, including electrical and electrical solvers.

**solcore/material\_data/**

Contains material properties, including refractive indeces,
bandstructure and transport properties, among others. The data is
organised in subfolders, depending on where the data is coming from. The
*mobility.py* module contains the functions needed to calculate the
mobility of several semiconductor materials.

**solcore/material\_system/**

Contains several modules to deal with materials. The most important one
is *material\_system.py* which contains Solcore\'s *material* class and
subclasses. The *critical\_point\_interpolate.py* and
*critical\_point\_picker.py* are used to interpolates the n and k data
of an alloy based on the known values at certain specific compositions.
This interpolation is done in a smart way based on certain critical
points (usually the Adachi critical points) and then filling the gaps in
between.

**solcore/parameter\_system/**

Contains just one module, *parameter\_system.py*, with the tools
necessary to access the materials data and, in case of derived
parameters, calculating them. This module is used extensively by the
*material\_system*, but it can be used on its own to access individual
properties, without creating a material.

**solcore/units\_system/**

Contains the module *units\_system.py* with functions to convert
magnitudes between different units in a high level fashion.

**solcore/light\_source/**

This folder contains the tools for creating and managing an illumination
spectra, typically, a solar spectrum.

-   The core functions are contained in the module *light\_source.py*.
-   *spectral2.py* contains a python implementation of the irradiance
    model developed by the NREL (
    <http://rredc.nrel.gov/solar/models/spectral/SPCTRAL2/> )
-   *smarts.py* is an interface to the SMARTS irradiance model, which
    needs to be installed separately from
    <http://www.nrel.gov/rredc/smarts/> .

**solcore/absorption\_calculator/**

Contains several tools to calculate the optical properties (dielectric
constants, n, k absorption coeficients) of bulk materials and quantum
wells as well as tools to calculate the propagation of light in stacks
of structures.

-   *absorption\_QW.py* calculates the absorption coefficient of quantum
    wells from the optput of the Schrödinger calculator.
-   *adachi\_alpha.py* uses Adachi parametrization to calculate the
    absorption coefficient of some bulk materials. The CPPB subfolder
    contains an advanced version of this method.
-   *dielectric\_constant\_models.py* contains a large variety of
    oscillator models, including Lorentz, Cauchy or Drude, which can be
    combined to create the total dielectric function of a material.
-   *sopra.py* contains the functions used to access the SOPRA database,
    in the **solcore/material\_data/** folder.
-   *kramers\_kronig.py* is used to calculate Kramers-Kronig compatible
    dielectric constants.
-   *transfer\_matrix.py* contains the class OptiStack, used across
    Solcore to represent optical structures, and functions to interface
    with the TMM calculator, the *tmm* Python package created by Steven
    Byrnes and included in the PyPi repository.
-   Finally, *rigorous\_coupled\_wave.py* is used to interface with the
    RCWA solver S4, from the Stanford University, which needs to be
    installed separately.

**solcore/optics/**

Includes several modules to calculate the reflected, transmitted and
absorbed light in a SolarCell object, using the Beer-Lamber law
(*beer\_lambert.py*), the transfer matrix method (*tmm.py*) and the
rigorous coupled-wave analysis (*rcwa.py*). The later two cases consists
mostly on an interface to the OptiStack and the solvers described in
**solcore/absorption\_calculator/**.

**solcore/analytic\_solar\_cells/**

This folder contains the tools needed to calculate the IV and QE
characteristics of solar cells using analytical or semianalitical models

-   *depletion\_aproximation.py* for the depletion approximation
-   *detailed\_balanced.py* for the detailed balance model
-   *diode\_equation.py* for the 2-diode model
-   The module *IV.py* contains the multi-junction solar cell
    calculator, combining the IV characteristics of the individual
    junctions with or without radiative coupling, as well as the legacy
    IV calculator, kept for backwards compatibility purposes.
-   The module *QE.py* as also been maintained just for backwards
    compatibility purposes, as now the QE is calculated in the above
    modules, depending on the chosen solar cell model.

**solcore/quantum\_mechanics/**

This folder includes all the tools related to the quantum properties of
materials and structures.

-   The *kp\_bulk.py* module solves the 8-band kp Hamiltonian for bulk
    materials under strain.
-   *kp\_QW.py* solves the 4-band and 6-band kp Hamiltonian for quantum
    wells (not fully implemented, yet).
-   *heterostructure\_alignment.py* uses the band offsets to align the
    conduction and valence bands of a heterostructure before staring any
    quantum calculation.
-   *strain.py* calculates the strain in a heterostructure and shifts
    the band edges accordingly.
-   *potential\_utilities.py* contains the time independent, 1D
    Schrödinger solver.
-   *structure\_utilities.py* uses the kp and strain calculators
    mentioned above to calculate the bands and efective mass profile of
    a heterostructure.
-   Finally, *high\_level\_kp\_qw.py* provides a common interface for
    the above solvers.

**solcore/poisson\_drift\_diffusion/**

This folders includes the Poisson-Drift diffusion (PDD) solver.

-   *DDmodel-current.f95* is the Fortran code, which needs ot be
    compiled in a library accessible by Python. This compilation is done
    by *driftdiffusion\_compiler.py* using F2Py.
-   *DeviceStructure.py* extracts from the materials making the solar
    cells all the properties needed by the solver, and also include
    tools for saving and loading structures from external files.
-   *QWunit.py* solves the properties of a QW structure calling the
    relevant functions in the **solcore/quantum\_mechanics/** module and
    creates and effective bulk-like stack of materials, usable by the
    PDD solver.
-   Finally, the *DriftDiffusionUtilities.py* contains all the functions
    to interface with the fortran library and solve the PDD equations
    under equilibrium, short circuit, calculate the IV curves or the QE.

**solcore/spice/**

The module *spice.py* contains the tools necessary to interface with the
electrical solver SPICE, which needs to be installed independently. The
other modules in these folder depend on this one. The module
*quasi\_3D\_solver.py* has the tools to solve calculate the IV curve of
a solar cell by modelling it as a 3D network of electrical components.
The module *pv\_module\_solver.py* has the tools for calculating the IV
curve of a solar module, with many solar cells connected in series and a
potential random dispersion of properties.

**solcore/data\_analysis\_tools/**

Contains modules designed to fit experimental data. For now, it only has
*ellipsometry\_analysis.py*, to fit ellipsometry data and get the
dielectric functions of a stack of materials.

**solcore/graphing/**

Contains several modules related to plotting data, intended to help
creating default graphs for complex data, for example a quantum well
with energy levels and wavefunctions.

**solcore/tests/**

Contains test scripts to be run using \"nosetests\", allowing to probe
that the main capabilities of Solcore work as expected and that there is
not regression after adding new functionality.

**docs/**

Contains Solcore\'s documentation, created with Sphinx. It has several
subfolders needed by Sphinx.

**examples/**

This folder and subfolders contain example scripts illustrating
Solcore\'s functionality. Most of them reproduce the figures in the main
Solcore paper (submitted to Computer Physics Communications, preprint in
<https://arxiv.org/abs/1709.06741> ), but there are other examples
expanding the rest of the capabilities.
