Sesame Poisson Drift-Diffusion solver
======================================

This solver provides an interface to `Sesame <https://sesame.readthedocs.io/en/latest/>`_ to solve the
Drift-Diffusion equations. While it serves the same purpose as the legacy Fortran PDD solver, this solver has
the advantage of being written entirely in Python, thus being more transparent to most users.
While the two solvers should in principle have the same functionality, the Fortran solver was developed
specifically to deal with solar cells containing quantum wells (QWs), where the abrupt and frequent change in material parameters
can cause convergence issues. The Sesame solver has not been tested for the QW use case. The Fortran solver is less
suitable for cells containing thick (> 10 microns or so) layers, as it was designed for III-V cells. The Sesame solver
was designed for silicon-based cells and thus can handle thicker layers.

Unlike the Fortran-based solver, this solver can handle both constant values of doping per layer,
or depth-dependent doping profiles defined by the user.

Material constants
-------------------

The following material parameters will be extracted from each layer in a junction and used by
Sesame:

- ``Nc``: Conduction band effective density of states
- ``Nv``: Valence band effective density of states
- ``band_gap``
- ``electron_affinity``
- ``relative_permittivity``
- ``electron_mobility``
- ``hole_mobility``
- ``electron_minority_lifetime``
- ``hole_minority_lifetime``
- ``bulk_recombination_energy``: this is an optional parameter (will be set to 0 by default).
  It should be set (in joules!), in the ``material`` definition, if the user wants to use a different value
- ``radiative_recombination``: radiative recombination rate (m\ :sup:`3` s\ :sup:`-1`)
- ``electron_auger_recombination``: Auger recombination rate for electrons (m\ :sup:`6` s\ :sup:`-1`)
- ``hole_auger_recombination``: Auger recombination rate for holes (m\ :sup:`6` s\ :sup:`-1`)

Note that this list uses the names of the parameters used by Solcore's material system, not the
names used internally by Sesame (which can be found in the Sesame documentation). While Sesame internally uses different units, values should be passed to
materials/layers in base SI units, as for the rest of Solcore (m, s, kg, etc).

Mesh
-----

Solcore will try to construct a reasonable mesh to solve the PDD equations, putting
more mesh points close to the front surface of the cell and in regions where the doping
is changing. However, the user can also provide a custom mesh, defined in terms of
distance *from the front surface of the junction* in m.
This is is passed as a 1-dimensional array through the ``mesh`` argument to ``Junction``).

Doping profile
---------------

As for the Fortran PDD solver and the depletion approximation (DA) solver, the user can set fixed
doping levels for each layer in a junction using the ``Nd`` argument for n-type doping and the ``Na``
argument for p-type. However, it is also possible to define depth-dependent doping profiles. These should be a function
which accepts an array of positions (in m) and returns the doping at that position in units of
m\ :sup:`-3`. Sesame will interpret positive values as n-type doping and negative values as p-type.
The doping profile can be specified for the whole junction by setting the ``doping_profile`` argument
of ``Junction``, or for individual layers by setting the ``doping_profile`` argument for ``Layer``. It is
possible to mix constant doping in some layers with depth-dependent doping in others. The position argument
of the function should always be in terms of distance from the front surface of the junction or layer it is for.

Solver options
---------------

Sesame has a number of options which can be set by the user to control the convergence
conditions of the Newton-Raphson solver, the maximum number of iterations, verbosity,
and the use of periodic boundary conditions. These options can be set as other Solcore
options passed to the ``solar_cell_solver`` function and are listed below, with the
name of the option in Solcore and, in brackets, the corresponding keyword argument in
Sesame:

- ``sesame_max_iterations`` (maxiter): integer. Maximum number of steps taken by the
 Newton-Raphson scheme. Default: 300
- ``sesame_tol`` (tol): float. Accepted error made by the Newton-Raphson scheme.
  Default: 1e-6
- ``sesame_verbose`` (verbose): boolean. The solver returns the step number and the associated
  error at every step, and this function prints the current applied
  voltage if set to True (default).
- ``sesame_htp`` (htp):  Number of homotopic Newton loops to perform. Default: 1.
- ``sesame_periodic`` (periodic_bcs): boolean. Defines the choice of boundary conditions in the
  y-direction. True (False) corresponds to periodic (abrupt) boundary conditions.
  Default: True

In addition, the following options control how Solcore interacts with Sesame during
quantum efficiency (QE) calculations:

- ``sesame_qe_parallel``: boolean. If True, the quantum efficiency calculation will be performed in parallel
  using the number of cores specified by ``sesame_qe_n_jobs``. If False, the QE calculation will be done
  sequentially. The downside of parallel execution is that some wavelengths may not converge without a good guess,
  usually provided from the solution of the closest wavelength. Default: False.
- ``sesame_qe_n_jobs``: integer. Number of cores to use for parallel quantum efficiency calculations.
  Can be set to -1 to use all available cores. If set to a positive integer, it will use that many cores.
- ``sesame_use_previous_wl``: boolean. If True, the quantum efficiency calculation will use the solution of the closest wavelength
  as a guess for the next wavelength (assuming the calculation is not running in parallel). If set to False,
  will use the guess provided by the user in the ``guess`` attribute of the junction.
  Default: True.
- ``sesame_qe_flux``: float, array-like, or LightSource. The flux of photons (in m\ :sup:`-2` s\ :sup:`-1`) to use
  in the quantum efficiency calculation. This is used to calculate the generation rate in the junction.
  This can be a float (constant flux for all wavelengths), an array-like object with the same length as the number of wavelengths,
  or a `LightSource` object which will be used to calculate the flux at each wavelength.
  If not provided, it will be set to a constant value of 1e4 m\ :sup:`-2` s\ :sup:`-1` for all wavelengths.

You can provide an initial guess to the IV, QE or equilibrium solver by setting the ``guess``
attribute of a ``Junction`` object. This should be a dictionary with keys 'potential', 'Efe' and 'Efh',
each of which is an array of length sys.nx (for IV calculations) or dimensions (n_wavelengths, sys.nx)
(for EQE calculations). These correspond to the ``potential``, ``Efe`` and ``Efh`` outputs of the solver
(see below).
The units for the potential are volts (V), while the quasi-Fermi levels are in electronvolts (eV),
which will be converted by Solcore to the unit system used by Sesame.


Outputs
-------
In addition to updating the individual junctions and the solar cell as a whole with the current-voltage output
(if ``solar_cell_solver`` is called with task ``iv``) and quantum efficiency
(if ``solar_cell_solver`` is called with task ``qe``), the Sesame solver will also update each junction it solves the
IV for with an attribute called ``pdd_output``, which contains the following:

- ``pdd_output.potential``: the potential (V)
- ``pdd_output.n``: the electron density (m\ :sup:`-3`)
- ``pdd_output.p``: the hole density (m\ :sup:`-3`)
- ``pdd_output.Ec``: the level of the conduction band (eV)
- ``pdd_output.Ev``: the level of the valence band (eV)
- ``pdd_output.Efe``: the electron quasi-Fermi level (eV)
- ``pdd_output.Efh``: the hole quasi-Fermi level (eV)
- ``pdd_output.G``: the generation rate as a function of position at the locations in the mesh (``junction.mesh``) (m\ :sup:`-3` s\ :sup:`-1`)
- ``pdd_output.Rrad``: the radiative recombination rate (m\ :sup:`-3` s\ :sup:`-1`)
- ``pdd_output.Raug``: the Auger recombination rate (m\ :sup:`-3` s\ :sup:`-1`)
- ``pdd_output.Rsrh``: the bulk recombination due to Shockley-Read-Hall processes (m\ :sup:`-3` s\ :sup:`-1`)

Each of these is a 2-dimensional array, with dimensions ``(len(options.internal_voltages), len(mesh))``,
with the exception of ``G`` which is 1D (only position-dependent, not a function of the voltage).

Sub-function documentation
---------------------------

Note that the functions below should generally not be accessed directly, but will be called by
``solar_cell_solver`` for junctions with ``kind="sesame_PDD"``.

.. automodule:: solcore.sesame_drift_diffusion.solve_pdd
    :members:
    :undoc-members: