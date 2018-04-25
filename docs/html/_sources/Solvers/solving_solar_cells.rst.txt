Solar cell solvers
==================

Solcore includes several methods to solve the properties of a solar cell, including optics and electrical properties. To solve the optics, :literal:`Solcore` have built in a transfer matrix solver and can be linked to S4, a rigorous couple wave analysis solver. For the electrical proerties, :literal:`Solcore` includes from the fundamental detailed balance solver to the more rigorous Poisson-drift-diffusion equation solver. The electrical solvers apply to the *individual junctions* separately, and then their ouput are combined to get the properties of a multi-junction device.

This section contains the description of the **solar cell solver module** where all calculations start, and links to the specific solvers. Visit the  :doc:`optical methods section <../Optics/optics>`for information specific to the optical solvers.

.. toctree::
    :maxdepth: 1

    detailed_balance
    TwoDiode
    depletion
    DDsolver
    multijunction_iv

The solar cell solver module
----------------------------

.. automodule:: solcore.solar_cell_solver
    :members:
    :undoc-members:

Solver options
--------------

All options that need to be passed to the solvers (either electrical or optical) are passed as a dictionary via de *user_options* keyword to the *solar_cell_solver* method described above. If no options are provided, the solver will try to use the default ones, which might or might not be adequate to your specific problem. These options are common for all calculations.

The options available as well as the default values are:

- General
    - T_ambient = 298
    - T = 298

- :doc:`Illumination spectrum <../spectral/spectral>`
    - wavelength = np.linspace(300, 1800, 251) * 1e-9
    - light_source = LightSource(source_type='standard', version='AM1.5g', x=default_options.wavelength,
                                           output_units='photon_flux_per_m')

- IV control
    - voltages = np.linspace(0, 1.2, 100)
    - mpp = False   # Calculates the maximum power point parameters
    - light_iv = False
    - internal_voltages = np.linspace(-6, 4, 1000)
    - position = np.linspace(0, 4660, 4661)
    - radiative_coupling = False

- Optics control
    - optics_method = 'BL'

- :doc:`Rigorous couple wave analysis options <../Optics/S4doc>`
    - size = [500, 500]
    - orders = 4
    - theta = 0
    - phi = 0
    - pol = 'u'

- :doc:`Detailed balance solver options <detailed_balance>`
    - db_mode = 'boltzmann'

- :doc:`Poisson-drift diffusion solver options <DriftDiffusionUtilities>`
    - Mesh control
        - meshpoints = -400
        - growth_rate = 0.7
        - coarse = 20e-9
        - fine = 1e-9
        - ultrafine = 0.2e-9

    - Convergence control
        - clamp = 20
        - nitermax = 100
        - ATol = 1e-14
        - RTol = 1e-6

    - Recombination control
        - srh = 1
        - rad = 1
        - aug = 0
        - sur = 1
        - gen = 0

