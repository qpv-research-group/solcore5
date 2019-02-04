SPICE-based solar cell solver
=============================

- Example 1: :doc:`PV module calculator <../Examples/example_pv_module>`
- Example 2: :doc:`Quasi-3D 3J solar cell <../Examples/example_quasi3D_cell>`

﻿When the two diode model is used to define the junctions in a MJ solar cell, then larger scale circuits can be constructed. Solcore includes two levels of large scale equivalent circuits: quasi-3D solar cell modelling and solar array modelling. Both solvers are based on the interface between Solcore and SPICE, allowing for a fast calculation of complex structures with many elements.

This solver has been tested with NGSPICE, only. Check the installation :doc:`instructions for Windows <../Installation/Solcore_on_Windows>` and :doc:`for MacOS <../Installation/Solcore_on_MacOSX>`.

.. toctree::
    :maxdepth: 0

    pv_panel
    quasi3D

Spice solver files
------------------

.. automodule:: solcore.spice.spice
    :members:
    :undoc-members:
