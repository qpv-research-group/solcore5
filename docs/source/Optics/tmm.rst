Transfer matrix method
======================

Example 1: :doc:`Using the TMM solver to calculate the reflextion of a multilayered ARC <../Examples/example_RAT_of_ARC>`

Example 2: :doc:`Looking at the effect of substrate and the no_back_reflection option in the TMM solver <../Examples/substrate_example>`

Example 3: :doc:`Coherent and incoherent layers in the TMM solver<../Examples/example_coherency>`


.. toctree::
    :maxdepth: 2

    tmm
	

The TMM interface for the solar cell solver
-------------------------------------------

This is the method actually called from the solar cell solver, serving as interface between the solar cell and the lower level TMM formalism. 
The :doc:`Beer-Lambert calculator <other_methods>`, the :doc:`RCWA calculator <S4doc>` and the :doc:`external optics calculator <other_methods>` 
(where the user simply adds the reflection and the absorption profile manually) have similar interfaces.

The TMM solver can handle different input angles and polarizations, specified in the options (specifically, options.theta and options.pol).
The input angle is in degrees, while pol is 's', 'p'  or 'u' (computationally, 'u' is simply the average of 's' and 'p' and thus requires
two calculations.)

Using the default options, all layers (except very thick layers, with the layer thickness more than 10 times the maximum incident wavelength)
will be considered as coherent, i.e. phase information in the wave-optical calculation is retained. However, layers can also be specified to be
incoherent using a coherency_list option with one entry (either 'c' or 'i') per layer. 

.. automodule:: solcore.optics.tmm
    :members:
    :undoc-members:

The OptiStack and the high level TMM optical methods
-----------------------------------------------------

.. automodule:: solcore.absorption_calculator.transfer_matrix
    :members:
    :undoc-members:

Transfer matrix engine
----------------------

.. automodule:: solcore.absorption_calculator.tmm_core_vec
    :members:
    :undoc-members:
