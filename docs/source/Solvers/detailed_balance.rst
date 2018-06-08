Detailed balance approximation
==============================

- Example: :doc:`MJ solar cell efficiency map <../Examples/example_MJ_efficiency_map>`

This solver calculates the electrical properties of the junction by
balancing the elementary processes taking place in the solar cell,
carrier generation and radiative recombination, using the formalism
described by Araújo and Martí (1996) ([#Ref20]_). The
method is widely used by the photovoltaic community to calculate the
limiting conversion efficiencies of the different solar cell
architectures or materials. The simplest DB formulation only needs an
absorption edge energy and an absorptivity value above that edge. Out of
this, the carrier generation and radiative recombination are calculated
for different internal chemical potentials, equal to the external
electrical bias, in the ideal case. Solcore includes this basic model,
but also allows the user to provide a more complex absorption profile.

For a more detailed description of the implementation of the DB solver in
Solcore, `refer to the main Solcore paper (open access) <https://doi.org/10.1007/s10825-018-1171-3>`_ and references therein.

Detailed balance functions
--------------------------

.. automodule:: solcore.analytic_solar_cells.detailed_balance
    :members:
    :undoc-members:

References
----------

.. [#Ref20] Martí, A., Araújo, G.L.: Limiting efficiencies for photovoltaic energy conversion in multigap systems. Sol. Energy Mater. Sol. Cells 43(2), 203–222 (1996)