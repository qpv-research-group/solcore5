Light Sources
=============

- Example: :doc:`Example of the use of the light_sources module <../Examples/example_light_sources>`

The Solcore module **light_source** is designed to deal easily with
different light sources. It has direct support for:

-  Gaussian emission, typical of lasers and light emitting diodes.

-  Black-body radiation, characteristic of halogen lamps defined by a
   temperature, but also used very often to simulate the spectrum of the
   Sun, very close to a black body source at 5800 K.

-  Standard solar spectra: the extraterrestial spectrum AM0 and the two
   terrestial ones, AM1.5D and AM1.5G as defined by the ASTM G173 -
   03(2008) standard.

-  Irradiance models, using location, time and atmospheric parameters to
   calculate a synthetic solar spectrum. Solcore includes two models:
   SPECTRAL2, fully implemented in Python, and an interface to SMARTS
   binaries (which need to be installed separately), which greatly
   simplifies its use in batch mode.

-  User-defined irradiances, provided externally from a database or any
   other source, allowing for maximum flexibility.

The syntax in all cases is simple and intuitive considering the type of
source that needs to be created. In the case of the irradiance models,
which often have a large number of inputs, Solcore defines a set of
default values, so only those that are different need to be provided.
These default values can be obtained by:

.. code-block:: Python

    import solcore.light_source as ls

    ls.get_default_spectral2_object()
    ls.get_default_smarts_object()

Once created, specific parameters of the light sources can be easily
modified without the need for creating the source from scratch. That is
particularly useful for the irradiance models, where we might be
interested in getting the spectrum as a function of a certain parameter
(e.g. the hour of the day, or the humidity) without changing the others.
For example, ``smarts.spectrum(HOUR=11)`` and
``smarts.spectrum(HOUR=17)`` will provide the spectrum of the SMART
light source (assuming it has already been created) calculated at 11h and at 17h, respectively;
all additional parameters have the default values. This method has been
used to model experimental solar irradiances measured by different
spectroradiometers based on the local atmospheric
conditions.

A final, very convenient feature of the LightSource class is the ability
to request the spectrum in a range of different units. The default is
power density per nanometer, but other common units are power density
per eV or photon flux per nanometer, among others. While these unit
conversions are straightforward, it is often an initial source of errors
due to missing constants or incompatible magnitudes.

The **light_source** module has been described in the context of the
solar spectrum, but it can be applied broadly where there is spectral
data involved, such as the fitting of photoluminescence,
electroluminescence or Raman spectra.


Description of functions and modules
------------------------------------

.. automodule:: solcore.light_source.light_source
    :members:
    :undoc-members:

.. automodule:: solcore.light_source.smarts
    :members:
    :undoc-members:

.. automodule:: solcore.light_source.spectral2
    :members:
    :undoc-members:
