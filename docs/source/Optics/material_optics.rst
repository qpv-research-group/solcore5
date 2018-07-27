Optical properties of materials
===============================

Solcore has several ways of accessing the optical properties of materials: databases and parametric dielectric functions.

.. toctree::
    :maxdepth: 2

    material_optics

Dielectric constants models and Adachi parametrization
------------------------------------------------------

Understanding the optical response of both established and novel
materials is crucial to effective solar cell design. To efficiently
model the complex dielectric function of a material Solcore incorporates
an optical constant calculator based on the well-known Critical-Point
Parabolic-Band (CPPB) formalism popularised by Adachi ([#Ref11]_, [#Ref12]_, [#Ref13]_). In this
model, contributions to *:math:`\epsilon_2(\omega)`* from critical
points in the Brillouin Zone at which the probability for optical
transitions is large (van Hove singularities) are considered. The
transition probability for such transitions is proportional to the joint
density of states (JDOS) :math:`\textbf{J}_{cv}(\omega)`, which
describes the number of available electronic states between the valence
and conduction bands at given photon energy. The imaginary part of the
complex dielectric function is related to the JDOS by:

.. math::

   \label{eqn:JDOS}
   \epsilon_2(\omega) = \frac{4 \hbar^2 e^2}{\pi \mu_{0}^{2} \omega^2} \left| \left \langle c | p | v \right \rangle \right|^2 \textbf{J}_{cv}(\omega)

Where :math:`\left| \left \langle c | p | v \right \rangle \right|` is
the momentum matrix element for transitions from the valence band
(:math:`v`) to the conduction band (:math:`c`). Critical point
transitions are considered at the following points of symmetry in the
band structure: *:math:`E_0`* corresponds to the optical transition at
the *:math:`\Gamma`* point and *:math:`E_0 + \Delta_0`* to the
transition from the spin-orbit split off band to the conduction band at
the *:math:`\Gamma`* point. *:math:`E_1`* and *:math:`E_1 + \Delta_1`*
denote the transitions from the valence heavy-hole (HH) band and the
valence light-hole (LH) band respectively to the conduction band at the
*L* point. The *:math:`E'_0`* triplet and *:math:`E_2`* transitions
occur at higher energies, between the HH band and the split conduction
bands at the :math:`\Gamma` point as well as across the wide gap *X*
valley. The model also includes contributions from the lowest energy
indirect band-gap transition and the exciton absorption at the
*:math:`E_0`* critical point. The contributions listed above are summed
to compute the overall value of :math:`\epsilon_2(\omega)`. The real and
imaginary components of the overall complex dielectric function
:math:`\epsilon(\omega) = \epsilon_1(\omega) - i\epsilon_2(\omega)` are
then related via the Kramers-Kronig relations;

.. math::

   \label{eqn:KKR_1}
       \epsilon_1(\omega) = 1 + \frac{2}{\pi} \int_{0}^{\infty} \frac{\omega' \epsilon_2(\omega')}{(\omega')^2 - \omega^2} d\omega'

.. math::

   \label{eqn:KKR_2}
       \epsilon_2(\omega) = - \frac{2}{\pi} \int_{0}^{\infty} \frac{\epsilon_1(\omega')}{(\omega')^2 - \omega^2} d\omega'

The CPPB model included with Solcore also incorporates a modification to
the critical point broadening present in Adachi’s description, which is
shown to produce a poor fit to experimental data in the vicinity of the
:math:`E_0` and :math:`E_1` critical points ([#Ref14]_). To give a more accurate description of
the broadening of the optical dielectric function, Kim et al. proposed
that a frequency-dependent damping parameter be used to replace the
damping constant given by Adachi at each critical point
([#Ref15]_, [#Ref16]_);

.. math::

   \label{eqn:Kim_damping}
       \Gamma'(\omega) = \Gamma exp \left[ -\alpha \left( \frac{\hbar \omega - E_0}{\Gamma}\right) ^2 \right]

Where :math:`\Gamma` is the damping constant used by Adachi and
:math:`\alpha` describes the shape of the lineshape broadening with
:math:`\alpha = 0` producing purely Lorentzian character and
:math:`\alpha = 0.3` producing a good approximation to Gaussian
broadening.

The Solcore module ``absorption_calculator`` contains the CPPB model
within the ``Custom_CPPB`` class. The class offers a flexible way to
build up the optical constant model by adding individual critical point
contributions through the *Oscillator* structure type within Solcore. In
addition to the oscillator functions described by Adachi the
``Custom_CPPB`` class also provides additional oscilator models and the
Sellmeier equation for describing the real part of the dielectric
function for non-absorbing materials ([#Ref17]_).

**Description of functions in this module**

.. automodule:: solcore.absorption_calculator.dielectric_constant_models
    :members:
    :undoc-members:

.. automodule:: solcore.absorption_calculator.adachi_alpha
    :members:
    :undoc-members:

SOPRA database
--------------

n order to calculate and model the optical response of potential solar
cell architectures and material systems, access to a library of accurate
optical constant data is essential. Therefore, Solcore incorporates a
resource of freely available optical constant data measured by Sopra S.
A. and provided by `Software Spectra Inc <http://www.sspectra.com/sopra.html>`_.
The refractive index :math:`n` and
extinction coefficient :math:`k` are provided for over 200 materials,
including many III-V, II-VI and group IV compounds in addition to a
range of common metals, glasses and dielectrics.

Any material within the Sopra S.A. optical constant database can be
used with the :doc:`“material” function <../Systems/Materials>`, but they will have
only the optical parameters :math:`n` and :math:`k`. In the case of
materials that are in both databases, the keyword “sopra” will need to
be set to “True” when creating the material. Once a material is loaded
its :math:`n`, :math:`k` and absorption coefficient data is returned by
calling the appropriate method, for example ``SiO2.n(wavelength)`` and
``SiO2.k(wavelength)``. For certain materials in the database, the
optical constants are provided for a range of alloy compositions. In
these cases, any desired composition within the range can be specified
and the interpolated :math:`n` and :math:`k` data is returned.

.. automodule:: solcore.absorption_calculator.sopra_db
    :members:
	
Manually changing optical constants of a material
------------------------------------------------------

If you would like to define a material with optical constant data from a file, 
you can do this by telling Solcore the path to the optical constant data, e.g.:

.. code-block:: python

    this_dir = os.path.split(__file__)[0]

    SiGeSn = material('Ge')(T=T, electron_mobility=0.05, hole_mobility=3.4e-3)

    SiGeSn.n_path = this_dir + '/SiGeSn_n.txt'
    SiGeSn.k_path = this_dir + '/SiGeSn_k.txt'
	
In this case, we have defined a material which is like the built-in Solcore germanium
material, but with new data for the refractive index and extinction coefficient from 
the files SiGeSn_n.txt and SiGeSn_k.txt, respectively, which are in the same folder 
as the Python script. The format of these files is tab-separated, with the first column
being wavelength (in nm) and the second column n or k. 

References
----------

.. [#Ref11] Adachi, S.: Model dielectric constants of GaP, GaAs, GaSb, InP, InAs, and InSb. Phys. Rev. B 35(14), 7454–7463 (1987)
.. [#Ref12] Adachi, S.: Optical dispersion relations for GaP, GaAs, GaSb, InP, InAs, InSb, Alx Ga1−x As, and In1−x Gax Asy P1−y . J. Appl. Phys. 66(12), 6030–6040 (1989)
.. [#Ref13] Adachi, S.: Optical dispersion relations for Si and Ge. J. Appl. Phys. 66(7), 3224–3231 (1989)
.. [#Ref14] Rakic ́, A.D., Majewski, M.L.: Modeling the optical dielectric func- tion of GaAs and AlAs: extension of Adachi’s model. J. Appl. Phys. 80(10), 5909–5914 (1996)
.. [#Ref15] Kim, C.C., Garland, J.W., Abad, H., Raccah, P.M.: Modeling the optical dielectric function of semiconductors: extension of the critical-point parabolic-band approximation. Phys. Rev. B 45(20), 11 749–11 767 (1992)
.. [#Ref16] Kim, C.C., Garland, J.W., Raccah, P.M.: Modeling the optical dielectric function of the alloy system AlxGa1-xAs. Phys. Rev. B 47(4), 1876–1888 (1993)
.. [#Ref17] Woollam, J.A.: Guide to using WVASE 32: spectroscopic ellipsometry data acquisition and analysis software. J. A. Woollam Company (2008). `https://books.google.co.uk/books? id=xOupYgEACAAJ`_