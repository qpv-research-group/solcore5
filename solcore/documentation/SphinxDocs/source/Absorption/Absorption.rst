Absorption calculator
=====================

This package contains the tools necessary to claculate the absorption coefficient of either bulk materials or quantum wells. This package also containes the tools needed to import optical data from the SOPRA database and to calculate the optical properties of a stack of materials using a transfer matrix formalism.

Some of these functions can be accessed directly when importing this package:

- create_adachi_alpha (for bulk materials)
- SOPRA_DB (for bulk materials)
- calc_alpha (for QWs)
- calc_emission (for QWs)
- Calculate_absorption_profile, calculate_rat, calculate_ellipsometry, OptiStack, DielectricConstantModel (transfer matrix tools)


Absorption of bulk materials using Adachi's formulation
-------------------------------------------------------
.. automodule:: solcore.absorption_calculator.adachi_alpha
    :members:
    :undoc-members:

Absorption of bulk materials from the SOPRA database
----------------------------------------------------
.. automodule:: solcore.absorption_calculator.sopra_db
    :members:
    :undoc-members:

Absorption of quantum wells
---------------------------
.. automodule:: solcore.absorption_calculator.absorption_QW
    :members:
    :undoc-members:

Optical properties of a stack of materials
------------------------------------------
.. automodule:: solcore.absorption_calculator.transfer_matrix
    :members:
    :undoc-members:

.. automodule:: solcore.absorption_calculator.dielectric_constant_models
    :members:
    :undoc-members:


