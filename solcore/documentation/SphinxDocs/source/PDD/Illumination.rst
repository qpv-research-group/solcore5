============
Light source
============

Contains the information to create an illumination source to be used by th PDD solver. It can use either a standard spectrum (default) or get a custom one as input. It includes utilities for changing the total irradiance or appliting a filter to the source.  

The Illumination class
----------------------

.. automodule:: solcore.poisson_drift_diffusion.Illumination
    :members:
    :undoc-members:

Example 1:
::

    # We use the above methods to create several light sources often found in PV

    import numpy as np
    import solcore3.PDD as PDD
    # If we only want access to this module, we could do:
    # 	import solcore3.PDD.Illumination as Illu

    # Spectrum for space PV: AM0, irradiance of 1 sun, unfiltered, default wavelength range given by the input file 'astm173.csv'.
    MySpaceSun = PDD.Illumination( spectrum='AM0')

    # Spectrum for flat solar panels: AM1.5g, irradiance of 1 sun, unfiltered, custom wavelength range
    wl = np.linspace(400, 1300, 4)
    MyFlatPanelSun = PDD.Illumination( spectrum='AM1.5g', wavelengths=wl)

    # Spectrum for CPV: AM1.5d (default), irradiance of 500 sun, filtered by a layer with OD=3 and absorption edge at 800 nm,
    # default wavelength range given by the input file 'astm173.csv'.
    MyCPVPanelSun = PDD.Illumination( irradiance=500 )
    MyCPVPanelSun.filter( edge=800e-9, OD=3 )


References
----------

.. [#Ref4] Reference AM1.5 Spectra. NREL. http://rredc.nrel.gov/solar/spectra/am1.5/