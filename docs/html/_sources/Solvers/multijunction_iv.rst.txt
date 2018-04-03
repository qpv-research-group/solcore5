Multi-junction electrical solver
================================

﻿A complete photovoltaic solar cell can include one or more junctions, metal contacts, optical layers (including anti-reflective coatings and nano-photonic structures) and tunnel junctions. The junctions, in turn, might range from simple PN homojunctions to complex heterojunctions, including multi-quantum well structures. The solvers described so far only calculate the properties of single junction devices. To combine them into a multi-junction device, it is necessary to consider that the individual junctions are electrically connected in series and the potential coupling of light emitted by the wider bandgap junctions into those with smaller bandgap.

Tunnel junctions
----------------

.. automodule:: solcore.analytic_solar_cells.tunnel_junctions
    :members:
    :undoc-members:

Multi-junction IV calculator
----------------------------

.. automodule:: solcore.analytic_solar_cells.IV
    :members:
    :undoc-members:
