Poisson - Drift-Diffusion solver (PDD)
======================================

The PDD package provide all tools necesary to build a solar cell structure and calculate its properties by solving simultaneously the Poisson equation and the drfit diffusion equations. 

**TO BE ADDED**: Some theory of drift diffusion and basics of the numerical solver.

For using the PDD package, it is enough to include the following line in your code:
::

	import solcore3.PDD as PDD
	
With this, all the functionality of the package will be available to the user. The actual functions and calculations are spread in several modules, but they don't need to be imported individually (although it can be done that way, too!).
	
#.	:doc:`Device Structure <DeviceStructure>`

	File: solcore3/PDD/DeviceStructure.py

	Contains the functions necesary to build a sample structure that can be read by the DD solver, as well functions for saving and loading that structure to/from a file. 
	
#.	:doc:`Drift Diffusion Utilities <DriftDiffusionUtilities>`

	File: solcore3/PDD/DriftDiffusionUtilities.py

	Contains the python interface that dumps all the information of the device structure into the fortran variables and that executes the chosen "virtual expereriment". 

#.	:doc:`QW unit creator <QWunit>`

	File: solcore3/PDD/QWunit.py
	
	Contains utilities that transform the device structure into a *solcore3.Structure* that, in turn, can be used to solve the Schrodinger equation and the kp model. It also prepares the properties of the structure (bandedges, efective density of states, etc) in order to have a meaningful set of properties for the DD solver. 
	
#.	:doc:`Light source <Illumination>`

	File: solcore3/PDD/Illumination.py
	
	Contains the information to create an illumination source to be used by th PDD solver. It can use either a standard spectrum or get a custom one as input. It includes utilities for changing the total irradiance or appliting a filter to the source.  
	
#.	:doc:`The mobility module <Mobility>`

	File: solcore3/PDD/Mobility.py
	
	Contains a mobility calculator that is used to build the device structure. It can also be used standalone.  

#.	**Drift Diffusion Fortran solver**

	File: solcore3/plugins/DDmodel/DDmodel-current.f95
	
	This is the current version of the Fortran code that (once compiled) performs all the heavy numerical calculations. Normally you would not need to care about this file except if you intend to modify the numerical solver itself. The compilation should be done automatically by :literal:`solcore3`, which  will also use *F2Py* to produce a python-readable module with all its functions and variables, any time you modify this file. 
	
	This python module is not imported with the PDD package but it is controlled separatelly by the plugins system of :literal:`solcore3`. Normally, the user will not need to access directly the functions it contains - doing that is precisly the purpose of the above modules - but if necesary, they can be used by doing:
	::

		import solcore3.plugins.drift_diffusion_fortran as dd	
	 
	
\---------------------

.. toctree::
   :maxdepth: 2

   DeviceStructure
   DriftDiffusionUtilities
   QWunit
   Illumination
   Mobility