Poisson - Drift-Diffusion solver (PDD)
======================================

The PDD package provide all tools necesary to build a solar cell structure and calculate its properties by solving simultaneously the Poisson equation and the drfit diffusion equations. Normally, these functions will not need to be accessed directly, but are called internally by :literal:`Solcore` when using the higher level methods in the :doc:`solar cell solver <solving_solar_cells>`.

For using the PDD package, it is enough to include the following line in your code:
::

	import solcore.PDD as PDD
	
With this, all the functionality of the package will be available to the user. The actual functions and calculations are spread in several modules, but they don't need to be imported individually (although it can be done that way, too!).
	
#.	:doc:`Device Structure <DeviceStructure>`

	File: solcore/PDD/DeviceStructure.py

	Contains the functions necesary to build a sample structure that can be read by the DD solver.
	
#.	:doc:`Drift Diffusion Utilities <DriftDiffusionUtilities>`

	File: solcore/PDD/DriftDiffusionUtilities.py

	Contains the python interface that dumps all the information of the device structure into the fortran variables and that executes the chosen "virtual expereriment". 

#.	:doc:`QW unit creator <QWunit>`

	File: solcore/PDD/QWunit.py
	
	Contains utilities that transform the device structure into a *structure* that, in turn, can be used to solve the Schrodinger equation and the kp model. It also prepares the properties of the structure (bandedges, efective density of states, etc) in order to have a meaningful set of properties for the DD solver.

#.	**Drift Diffusion Fortran solver**

	File: solcore/PDD/DDmodel-current.f95
	
	This is the current version of the Fortran code that (once compiled) performs all the heavy numerical calculations. Normally you would not need to care about this file except if you intend to modify the numerical solver itself. The compilation should be done automatically by :literal:`Solcore`, which  will also use *F2Py* to produce a python-readable module with all its functions and variables, any time you modify this file.
	
\---------------------

.. toctree::
   :maxdepth: 2

   DeviceStructure
   DriftDiffusionUtilities
   QWunit