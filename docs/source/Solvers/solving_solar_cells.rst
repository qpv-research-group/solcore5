Solar cell solvers
==================

Solcore includes several methods to solve the properties of a solar cell, including optics and electrical properties. To solve the optics, :literal:`Solcore` have built in a transfer matrix solver and can be linked to `S4 <http://github.com/phoebe-p/S4>`_, a rigorous couple wave analysis solver. For the electrical properties, :literal:`Solcore` includes from the fundamental detailed balance (DB) solver to the more rigorous Poisson-drift-diffusion (PDD) equation solver. The electrical solvers apply to the *individual junctions* separately, and then their ouput are combined to get the properties of a multi-junction device.

The two most important elements of the solar cell solver module are the **solar_cell_solver** function and the **default_options** variable (see :ref:`solver-options`). The former is the function to be called to calculate any property of any solar cell, regardless of how the junctions have been defined or the specific property of interest. It provides a common interface for any calculation. To use it, simply do:
    ::

        from solcore.solar_cell_solver import solar_cell_solver

        solar_cell_solver(my_solar_cell_object, task, user_options)

The *task* has to be "optics", "iv", "qe", "equilibrium" or "short_circuit", the last two only working for PDD junctions. The sequence followed by the solver is:

#.  The default options are updated with the user defined options
#.  The solar cell structure is scanned, calculating derived information like the width of the junctions or the offset of each layer or junction with respect the front of the cell.
#.  The solar cell object and the updated options are sent to the corresponding solver ("optics", "iv", "qe", "equilibrium" or "short_circuit"), depending on the chosen task.
#.  Except in the case of "Equilibrium", all other tasks will require calculating the optics of the cell. This is done according to the chosen *optics_method* selected in the options. The default is "BL" (Beer-Lambert law). Alternative *optics_method* values are "TMM", "RCWA", "external" or None, if optics do not need to be calculated. Visit the  :doc:`optical methods section <../Optics/optics>` for information specific to the optical solvers.
#.  For the "iv" and "qe" tasks, each of the junctions present in the solar cell will be solved, according to the method used in their definition ("PDD", "DA", "2D" or "DB"). Details of each of this methods are described in the corresponding section.
#.  Finally, for the "iv" task, the individual IV curves of the junctions and tunnel junctions, if present, are combined to calculate the IV curve of the multi-junction solar cell, taking into account radiative coupling, if required.

All of the above calculations modify the original solar cell object, adding attributes or methods to its structure. For example, after calculating the IV curve of my_solar_cell_object, this will have a new attribute called "iv" that is a dictionary with the IV of the total solar cell, the IV curves of each of the junctions, information of the Voc and Isc, if relevant, etc.

More details of the specific electrical solvers included in Solcore can be found in:

.. toctree::
    :maxdepth: 0

    detailed_balance
    TwoDiode
    depletion
    DDsolver
    multijunction_iv


.. _solver-options:

Solver Options
--------------

All options that need to be passed to the solvers (either electrical or optical) are passed as a dictionary via de *user_options* keyword to the *solar_cell_solver* method described above. If no options are provided, the solver will try to use the default ones, which might or might not be adequate to your specific problem. These options are common for all calculations.

The options available as well as the default values are:

- General
    - **T_ambient** = 298
        Ambient temperature (K)

    - **T** = 298
        Cell temperature (K). It is actually made equal to the temperature in the solar solar cell definition: my_solar_cell.T .

- :doc:`Illumination spectrum <../spectral/spectral>`
    - **wavelength** = np.linspace(300, 1800, 251) * 1e-9
        Wavelengths of the illumination spectrum (m)

    - **light_source** = LightSource(source_type='standard', version='AM1.5g', x=default_options.wavelength, output_units='photon_flux_per_m')
        The illumination spectrum Air Mass 1.5 global spectrum, provided at the above wavelengths in units of photons•s\ :sup:`-1` m\ :sup:`-2`.

- IV control
    - **voltages** = np.linspace(0, 1.2, 100)
        Voltages at which to calculate the IV curve of the *complete* solar cell.

    - **mpp** = False
        If the parameters of the solar cell under illumination (Isc, Voc, Pmpp, Impp, Vmpp, fill factor and efficiency) have to be calculated. If False, all of them will have the value None.

    - **light_iv** = False
        If the light IV curve is to be simulated.

    - **internal_voltages** = np.linspace(-6, 4, 1000)
        The voltages at which the IV curve of *each of the junctions* inside the cell have to be calculated. This range has to be wider than the **voltages** above, in general. The same voltage range will be used in all junctions.

    - **position** = np.arange(0, solar_cell.width, 1e-10)
        Positions inside the solar cell structure in which to calculate the absorption. By default, it is calculated each angstrom for the whole width of the cell. To control the depth spacing, the user can pass:

        #. a vector which specifies each position (in m) at which the depth should be calculated
        #. a single number which specifies the spacing (in m) to generate the position vector, e.g. 1e-9 for 1 nm spacing
        #. a list of numbers which specify the spacing (in m) to be used in each layer. This list can have EITHER the length of the number of individual layers + the number of junctions in the cell object, OR the length of the total number of individual layers including layers inside junctions.

    - **radiative_coupling** = False
        If radiative coupling has to be included in the calculation.

- Optics control
    - **optics_method** = 'BL'
        Default method to calculate the optical properties of the cell. Other possible values are "TMM", "RCWA", "external" or None.

- :doc:`Rigorous couple wave analysis options <../Optics/S4doc>`
    Check the RCWA section for details on this parameters.

    - **size** = [500, 500]
    - **orders** = 4
    - **theta** = 0
    - **phi** = 0
    - **pol** = 'u'

- :doc:`Detailed balance solver options <detailed_balance>`
    - **db_mode** = 'boltzmann'
        If the Boltzmann approximation should be used in the detailed balance solver. Any other choice will result in using the full Plank equation, which will be slower, in general.

- :doc:`Depletion approximation solver options <depletion>`
    - **da_mode** = 'bvp'
        Selects the numerical approximation method for the drift-diffusion equation in the depletion approximation solver. Possible values are “bvp” for numerical solution using the `solve_bvp` method of the `scipy.integrate` module or 'green' for a semi-analytic solution using Green's functions. The latter is expected to be faster. 

- :doc:`Poisson-drift diffusion solver options <DriftDiffusionUtilities>`
    Check the PDD section for details on this parameters.

    - Mesh control
        - **meshpoints** = -400
        - **growth_rate** = 0.7
        - **coarse** = 20e-9
        - **fine** = 1e-9
        - **ultrafine** = 0.2e-9

    - Convergence control
        - **clamp** = 20
        - **nitermax** = 100
        - **ATol** = 1e-14
        - **RTol** = 1e-6

    - Recombination control
        - **srh** = 1
        - **rad** = 1
        - **aug** = 0
        - **sur** = 1
        - **gen** = 0


Solar cell solver module functions
----------------------------------

.. automodule:: solcore.solar_cell_solver
    :members:
    :undoc-members:

