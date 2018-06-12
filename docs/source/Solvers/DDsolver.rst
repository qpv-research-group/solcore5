Poisson - Drift-Diffusion solver (PDD)
======================================

- Example 1: :doc:`Example of a simple 2J solar cell calculated with the PDD solver <../Examples/example_PDD_solver>`
- Example 2: :doc:`Tutorial: 2J solar cell with QWs in the bottom cell <../Examples/tutorial>`

The PDD package provide all tools necesary to build a solar cell structure and calculate its properties by solving simultaneously the Poisson equation and the drfit diffusion equations. Normally, these functions will not need to be accessed directly, but are called internally by :literal:`Solcore` when using the higher level methods in the :doc:`solar cell solver <solving_solar_cells>`.

For using the PDD package, it is enough to include the following line in your code:

.. code-block:: python

	import solcore.poisson_drift_diffusion as PDD
	
With this, all the functionality of the package will be available to the user. The actual functions and calculations are spread in several modules:
	
1.	:doc:`Drift Diffusion Utilities <DriftDiffusionUtilities>`
    **File:** solcore/PDD/DriftDiffusionUtilities.py
	Contains the python interface that dumps all the information of the device structure into the fortran variables and that executes the chosen calculation, getting the data from the fortran variables at the end. This module contains the following methods called internally by the solar_cell_solver:

	- equilibrium_pdd
	- short_circuit_pdd
	- iv_pdd
	- qe_pdd

2.	:doc:`Device Structure <DeviceStructure>`
    **File:** solcore/PDD/DeviceStructure.py
	Contains several functions necessary to build a structure that can be read by the PDD solver. The most important of this functions is **CreateDeviceStructure** that scans the junction layers and extracts from the materials database all the parameters required by the PDD solver. Finally, it contains the default properties that are used if they are not found for a particular material. In general, these default properties correspond to those of GaAs at 293 K.

3.	:doc:`QW unit creator <QWunit>`
    **File:** solcore/PDD/QWunit.py
	Contains utilities that transform the sequence of layers into a *structure* that, in turn, can be used to solve the Schrodinger equation and the kp model. It also prepares the properties of the structure (bandedges, efective density of states, etc) in order to have a meaningful set of properties for the DD solver, using the **GetEffectiveQW**.

4.	**Drift Diffusion Fortran solver**

	**File:** solcore/PDD/DDmodel-current.f95

      This is the current version of the Fortran code that, once compiled, performs all the heavy numerical calculations. Normally you would not need to care about this file except if you intend to modify the numerical solver itself.
	
\---------------------

.. toctree::
   :maxdepth: 2

   DeviceStructure
   DriftDiffusionUtilities
   QWunit