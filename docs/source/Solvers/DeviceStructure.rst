Device Structure
================

This module contains the functions necessary to build a sample structure that can be read by the PDD solver. Typically, you will not need to use any of these except in the case of including quantum wells in the solar cell. In this case, you will need to use this methods to solve the quantum properties and create an effective medium before executing the PDD solver.

.. _main-functions:

Main functions
--------------

.. py:function:: CreateDeviceStructure(name [, role='device', T=293, layers=[], comments='', repeat=1, substrate=DefaultMaterial, reflection=None])

	Creates a dictionary with subdctionaries and lists storing all the parameters and material information necessary to solve the Poisson - Drift-Diffusion equations. It might be useful, also, to get all the properties of a given material that will be used in the PDD solver.
	
	*Layers* is a list of objects of class *solcore.Layer*. A *solcore.Structure* or a *solcore.Junction* are also valid inputs since all of them are made of a list of layers. Any property not defined explicitly in the definition of the materials for the layers and that can not be calculated by :literal:`Solcore` is taken from the default material (see section :ref:`default-material`). *substrate* must be a *solcore.material*, otherwise the default material is used.

   	**Output**: dictionary with all the information of the structure.

.. code-block:: python

    # This example illustrates the creation of a p-i-n Junction including AlGaAs window and back surface field layers.
    # A custom absorption coefficient for the GaAs of the intrinsic layer is used.

    import solcore.poisson_drift_diffusion as PDD
    from solcore.structure import Layer, Junction
    from solcore.import material

    # First, we create the materials, overriding any default property we want, such as the doping or the absorption coefficient
    window      = material('AlGaAs')    (T=300, Na = 1e24, Al=0.8)
    p_GaAs      = material('GaAs')      (T=300, Na = 1e24)
    i_GaAs      = material('GaAs')      (T=300)
    n_GaAs      = material('GaAs')      (T=300, Nd = 1e23)
    bsf         = material('AlGaAs')    (T=300, Nd = 1e24, Al=0.4)

    # We put everything together in a Junction. We include the surface recombination velocities,
    # sn and sp, although they are not necessary in this case.
    MyJunction = Junction([ Layer(width = 30e-9,         material = window,           role="Window"),
                            Layer(width = 400e-9,        material = p_GaAs,           role="Emitter"),
                            Layer(width = 400e-9,        material = i_GaAs,           role="Intrinsic"),
                            Layer(width = 2000e-9,       material = n_GaAs,           role="Base"),
                            Layer(width = 200e-9,        material = bsf,              role="BSF")],
                            sn=1e3, sp=1e3, T=300, kind='PDD')

    # Then, we create the structure. What actually this "creation" does is getting all the information of the materials
    # from the materials database relevant for the PDD solver and storing them in a dictionary.
    MyDevice = PDD.CreateDeviceStructure( 'TestDevice', T = MyJunction.T, layers = MyJunction)

Now "MyDevice" is a dictionary with several entries with the specific information that will be used by the PDD solver (if required). For example :code:`MyDevice['layers'][1]['properties']` will contain the following (the order might be different):

.. code-block:: python

    {'composition': {'material': 'GaAs'},
    'width': 4e-07,
    'band_gap': 2.279067404056071e-19,
    'electron_affinity': 6.62903371354393e-19,
    'eff_mass_electron_Gamma': 0.067,
    'eff_mass_hh_z': 0.34129421032745344,
    'eff_mass_lh_z': 0.0879370668050916,
    'electron_mobility': 0.9399999990563772,
    'hole_mobility': 0.017374281562790292,
    'ni': 1767480124457.733,
    'Nc': 4.3519622483962564e+23,
    'Nv': 5.657793654818883e+24,
    'electron_minority_lifetime': 3e-06,
    'hole_minority_lifetime': 2.5e-07,
    'permittivity': 12.9,
    'electron_auger_recombination': 1e-42,
    'hole_auger_recombination': 1e-42,
    'radiative_recombination': 2.8744203894058427e-17,
    'Nd': 1, 'Na': 1e+24,
    'sn': 1000000.0,
    'sp': 1000000.0}

The information for layers 2 and 3 wil be similar since we have GaAs in all cases, except for the mobilities as these depend on the doping, which is different.


.. _default-material:

Default material
----------------

The default material is **GaAs at 293K**. In general, all the intrinsic properties of the materials are taken from the literature and calculated by the :literal:`Solcore` material system. If not found there (for example if the material or alloy is not in the database) or if they are extrinsic, the default values below are used. The extrinsic properties, that depend on the quality of the material or its doping, are just asigned a 'reasonable' value. The user must make sure that this 'reasonable' value is indeed resonable for the intended application and override it with his/her own, otherwise.

The total list of parameters and the default values are:

============================ ============================ ====================== =============================================
Parameter                    Default value                Units                  Notes
============================ ============================ ====================== =============================================
band_gap                     Calculated for GaAs at 293K  J                      \-
electron_affinity            Calculated for GaAs at 293K  J                      \-
eff_mass_electron_Gamma	     Calculated for GaAs at 293K  relative to m0         \-
eff_mass_hh_z                Calculated for GaAs at 293K  relative to m0         \-
eff_mass_lh_z                Calculated for GaAs at 293K  relative to m0         \-
permittivity                 12.9                         relative to epsilon0   \-
electron_mobility            Calculated for GaAs at 293K  m2 V-1 s-1             \-
hole_mobility                Calculated for GaAs at 293K  m2 V-1 s-1             \-
electron_minority_lifetime   3e-6                         s                      SRH [#f1]_ recombination time for electrons
hole_minority_lifetime       2.5e-7                       s                      SRH recombination time for holes
electron_auger_recombination 1e-42                        m6 s-1                 \-
hole_auger_recombination     1e-42                        m6 s-1                 \-
radiative_recombination	     Calculated for GaAs          m3 s-1                 Derived from the absorption coefficient
Nd                           1                            m-3                    Density of donors
Na                           1                            m-3                    Density of acceptors
============================ ============================ ====================== =============================================

.. [#f1] Shockley-Read-Hall

All functions description
-------------------------

.. automodule:: solcore.poisson_drift_diffusion.DeviceStructure
    :members:
    :undoc-members:
