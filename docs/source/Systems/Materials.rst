The Material System
===================

The material class
------------------

.. automodule:: solcore.material_system.material_system
    :members:
    :undoc-members:

The mobility module
-------------------

This module allows to calculate the carrier mobilities based on the material composition, temperature (T>100K) and impurity concentration. It is an implementation of the mobility model by Sotoodeh *et al.* ([#Ref6]_). The material class described above uses this module internally to get the mobililty of the materials it is implemented for.

The material parameters used in the model are included in the file *mobility_parameters.json*. At the moment, the supported materials are:

	  - **Binaries**: AlAs, GaAs, InAs, InP, GaP
	  - **Ternaries**: InAlAs, AlGaAs, InGaAs, GaInP, (GaAsP), (InAsP)
	  - **Quaternaries**: InGaAsP

The last two ternaries are simply calculated as a linear interpolation of the parameters of the corresponding binaries, so the transition from direct to indirect bandgap and other dependencies might not be very accurate. *InGaAsP*, in turn, is calculated only based on the *InGaAs* and *GaInP* data and it is mostly valid for compositions lattice matched to *InP*.

.. automodule:: solcore.material_data.mobility
    :members:
    :undoc-members:

The critical point picker
-------------------------

.. automodule:: solcore.material_system.critical_point_picker
    :members:
    :undoc-members:

References
----------

.. [#Ref6] M. Sotoodeh, A. H. Khalid, and A. A. Rezazadeh, “Empirical low-field mobility model for III–V compounds applicable in device simulation codes,” *J. Appl. Phys.*, **87**, 2890, (2000).