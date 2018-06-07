Transfer matrix method
======================

Example 1: :doc:`Using the TMM solver to calculate the reflextion of a multilayered ARC <../Examples/example_RAT_of_ARC>`

.. toctree::
    :maxdepth: 2

    tmm

The OpticStack and the high level TMM optical methods
-----------------------------------------------------

.. automodule:: solcore.absorption_calculator.transfer_matrix
    :members:
    :undoc-members:

The TMM interface for the solar cell solver
-------------------------------------------

This is the method actually called from the solar cell solver, serving as interface between the solar cell and the lower level TMM formalism. The :doc:`Beer-Lambert calculator <other_methods>`, the :doc:`RCWA calculator <S4doc>` and the :doc:`external optics calculator <other_methods>` (where the user simply adds the reflection and the absorption profile manually) have similar interfaces.

.. automodule:: solcore.optics.tmm
    :members:
    :undoc-members:


Transfer matrix engine
----------------------

.. automodule:: solcore.absorption_calculator.tmm_core_vec
    :members:
    :undoc-members:
