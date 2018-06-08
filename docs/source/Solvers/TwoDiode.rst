2-diode equation
================

- Example 1: :doc:`MJ solar cell with radiative coupling <../Examples/example_radiative_coupling>`
- Example 2: :doc:`Quasi-3D 3J solar cell <../Examples/example_quasi3D_cell>`
- Example 3: :doc:`PV module calculator <../Examples/example_pv_module>`

This is the simplest method for simulating the behaviour of a solar
cell, using electrical components to model the different transport and
recombination mechanisms of the device. The 2D model is widely applied
when modelling solar cells at the most engineering end of the topic,
when a detailed knowledge of the solar cell structure (layers,
absorption coefficient, etc.) are not known or sought. It is often used
to fit experimental IV curves and find approximate, general information
on the solar cell quality without entering on the fundamental processes.
It can provide valuable information to engineers, when designing solar
modules for example, or for diagnostic purposes The complete form of the
equation is:

.. math::

   \label{eq:2diode}
   \begin{split}
   J = J_{sc} & - J_{01} \left( e^ \frac{q(V-R_sJ)}{n_1 k_b T} - 1 \right) \\
   & - J_{02} \left( e^ \frac{q(V-R_sJ)}{n_2 k_b T} - 1 \right) \\
   & - \frac{V-R_sJ}{R_{sh}}
   \end{split}

Generally, the photocurrent is modelled as a current source
(:math:`J_{sc}`), with radiative and non-radiative recombination
modelled as two diodes with reverse saturation currents :math:`J_{01}`
and :math:`J_{02}`, and ideality factors :math:`n_1\approx 1` and
:math:`n_2\approx 2`, respectively. The shunt resistance :math:`R_{sh}`
accounts for alternative current paths between the contacts of the solar
cell, being infinite in the ideal case, and the series resistance
:math:`R_s` accounts for the other transport losses. The values of the
saturation currents and ideality factors can, ultimately, be calculated
from the material properties and device structure, as is done in the
:doc:`depletion approximation model <depletion>`, but the 2D model
allows them to be provided directly as input, obtained from a fit to
experimental data, for example. They can also be calculated internally,
:doc:`using the DB solver <detailed_balance>` to obtain :math:`J_{01}` and :math:`J_{sc}`, and
then using a radiative efficiency coefficient to obtain :math:`J_{02}`.
The radiative efficiency :math:`\eta` is defined as the fraction of
radiative current :math:`J_{rad}` at a given reference total current
:math:`J_{ref}`:

.. math::

   \label{eq:eta_rad}
   \eta = \frac{J_{rad}}{J_{ref}} = \frac{J_{01}}{J_{ref}} \left( e ^{\frac{qV_{ref}}{n_1k_bT}} - 1 \right)

The reference voltage :math:`V_{ref}` can be written as a function of
:math:`\eta` and :math:`J_{ref}` as:

.. math::

   \label{eq:vref}
   V_{ref} = \frac{n_1k_bT}{q} \log \left( \frac{\eta J_{ref}}{J_{01}} + 1 \right)

On the other hand, the radiative coefficient can also be written as:

.. math::

   \label{eq:eta_nrad}
   \eta = \frac{J_{ref} - J_{nrad} - V_{ref}/R_{sh}}{J_{ref}}

Combining the last two equations and using the expression for
the diode with ideality factor :math:`n_2`, :math:`J_{02}` can be
written as:

.. math::

   \label{eq:J02}
   J_{02} = \frac{(1-\eta) J_{ref} - V_{ref} / R_{sh}}{e^ {\frac{qV_{ref}}{n_2k_bT} } - 1 }

In the common situation of very large shunt resistance and
:math:`V_{ref} >> k_bT/q`, this equation further simplifies to:

.. math::

   \label{eq:J02_simple}
   J_{02} = (1-\eta) J_{ref} \left( \frac{J_{01} }{ J_{ref} \eta } \right)^{n_1/n_2}

This process can, of course, be reversed to use knowledge of
:math:`J_{01}` and :math:`J_{02}` at a given reference current to
calculate the radiative efficiency of a solar cell, which is useful to
compare different materials, technologies or processing methods. This
was done by Chan et al. using :math:`J_{ref} = 30` mA/cm\ :math:`^2`,
obtaining :math:`\eta` values of 20% for InGaP, 22% for GaAs, and 27%
for InGaAs devices ([#Ref21]_). It should be pointed
out that this method is only valid under the assumption that
:math:`J_{01}` corresponds only to radiative recombination and
:math:`J_{02}` only to non-radiative recombination, which is generally
true for QW solar cells and some III-V solar cells, like those made of
GaAs or InGaP, but not for Si or Ge, for example. Other definitions of
the radiative efficiency are based on the external quantum efficiency,
the I\ :math:`_{sc}` and V\ :math:`_{oc}` of the cell, as described by
Green (2011) ([#Ref22]_).

Despite the simplicity of the 2-diode model, it is very useful to guide
the design of new solar cells and explore the performance of new
materials, such as dilute bismuth
alloys ([#Ref23]_), or to asses the performance of
large arrays of solar cells ([#Ref24]_).


2-diode equation functions
--------------------------

.. automodule:: solcore.analytic_solar_cells.diode_equation
    :members:
    :undoc-members:

References
----------

.. [#Ref21] Chan, N.L.A., Ekins-Daukes, N.J., Adams, J.G.J., Lumb, M.P., Gonzalez, M., Jenkins, P.P., Vurgaftman, I., Meyer, J.R., Walters, R.J.: Optimal bandgap combinations—does material quality mat- ter? IEEE J. Photovolt. 2(2), 202–208 (2012)
.. [#Ref22] Green, M.A.: Radiative efficiency of state-of-the-art photovoltaic cells. Prog. Photovolt. Res. Appl. 20(4), 472–476 (2011)
.. [#Ref23] Thomas, T., Mellor, A., Hylton, N.P., Führer, M., Alonso-Álvarez, D., Braun, A., Ekins-Daukes, N.J., David, J.P.R., Sweeney, S.J.: Requirements for a GaAsBi 1 eV sub-cell in a GaAs-based multi- junction solar cell. Semicond. Sci. Technol. 30(9), 094010-6 (2015)
.. [#Ref24] Ekins-Daukes, N.J., Kemmoku, Y., Araki, K., Betts, T.R., Gottschalg, R., Infield, D.G., Yamaguchi, M.: The design specifica- tion for Syracuse; a multi-junction concentrator system computer model. In: Proceedings of the 19th European Photovoltaic Solar Energy Conference, pp. 1–4 (2004)