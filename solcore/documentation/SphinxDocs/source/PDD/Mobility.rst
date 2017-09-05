===================
The mobility module
===================

This module allows to calculate the carrier mobilities based on the material composition, temperature (T>100K) and impurity concentration. It is an implementation of the mobility model by:

M. Sotoodeh, A. H. Khalid, and A. A. Rezazadeh, “Empirical low-field mobility model for III–V compounds applicable in device simulation codes,” *J. Appl. Phys.*, **87**, 2890, (2000).

Supported materials
-------------------

The material parameters used in the model are included in the file *mobility_parameters.json*. At the moment, the supported materials are:

	  - **Binaries**: AlAs, GaAs, InAs, InP, GaP
	  - **Ternaries**: InAlAs, AlGaAs, InGaAs, GaInP, (GaAsP), (InAsP)
	  - **Quaternaries**: InGaAsP

The last two ternaries are simply calculated as a linear interpolation of the parameters of the corresponding binaries, so the transition from direct to indirect bandgap and other dependencies might not be very accurate. *InGaAsP*, in turn, is calculated only based on the *InGaAs* and *GaInP* data and it is mostly valid for compositions lattice matched to *InP*. 

Functions
---------

.. py:function:: def calculate_mobility(material, holes, N, [x=0.0, y=0.0, T=300])

	Calculates the mobility using the model by Sotoodeh *et al.* ([#Ref6]_). If the material is not in the database, then the function returns the mobility for *GaAs* at that temperature, *T*, and impurity concentration, *N*. *holes* = 0 calculates the mobility for electrons and *holes* = 1 calculate the mobility for holes.  *x* is the fractional composition in the case of ternaries and *y* in the case of quaternaries.  

	**Output**: Mobility

All functions description
-------------------------

.. automodule:: solcore.poisson_drift_diffusion.mobility
	:members:
	
References
----------

.. [#Ref6] M. Sotoodeh, A. H. Khalid, and A. A. Rezazadeh, “Empirical low-field mobility model for III–V compounds applicable in device simulation codes,” *J. Appl. Phys.*, **87**, 2890, (2000).