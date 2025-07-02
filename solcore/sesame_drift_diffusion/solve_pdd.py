from __future__ import annotations

from solsesame.builder import Builder
from solsesame.solvers import Solver
from solsesame.analyzer import Analyzer
import numpy as np
from solcore.constants import q
from scipy.interpolate import interp1d
from solcore.state import State
from solcore.structure import Junction
from solcore.sesame_drift_diffusion.process_structure import process_structure
from solcore.light_source import LightSource

import warnings
from joblib import Parallel, delayed

from solcore.registries import register_iv_solver, register_equilibrium_solver


def process_sesame_options(options):

    sesame_kwargs = {}

    if "sesame_max_iterations" in options:
        sesame_kwargs["maxiter"] = options.sesame_max_iterations

    if "sesame_tol" in options:
        sesame_kwargs["tol"] = options.sesame_tol

    if "sesame_verbose" in options:
        sesame_kwargs["verbose"] = options.sesame_verbose

    if "sesame_htp" in options:
        sesame_kwargs["htp"] = options.sesame_htp

    if "sesame_periodic" in options:
        sesame_kwargs["periodic_bcs"] = options.sesame_periodic

    return sesame_kwargs


@register_equilibrium_solver("sesame_PDD")
def equilibrium(junction: Junction, **kwargs):
    """Solve at equilibrium (no illumination, no applied voltage) using the Sesame
    solver.

    :param junction: a Junction object
    :param **kwargs: options, passed as keyword arguments
    """

    options = State(
        kwargs
    )  # needs to be passed as kwargs to be compatible with Fortran equilibrium solver
    # in registries

    solver_class = Solver()
    IVcurve = solver_class.IVcurve

    if not hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    if hasattr(junction, "guess"):
        # guess for the non-equilibrium solution at V = 0
        guess_sesame = process_guess(junction.guess, junction.sesame_sys)

    else:
        guess_sesame = None

    sesame_kwargs = process_sesame_options(options)

    j, result = IVcurve(junction.sesame_sys, [0], guess=guess_sesame, **sesame_kwargs)

    if np.any(np.isnan(j)):
        warnings.warn(
            "Current calculation did not converge at all voltages", UserWarning
        )

    # j = j * junction.sesame_sys.scaling.current * 1e4  # cm-2 -> m-2

    junction.sesame_output = result


def process_guess(guess, sys):
    """Process the guess for the Sesame solver. This is necessary because the Sesame
    solver requires the guess to be in the form of a dictionary with keys 'v', 'efn'
    and 'efp', each of which is an array of length sys.nx (for IV calculations) or
    (n_wavelengths, sys.nx) (for EQE calculations). Units also need to be
    converted to those used by Sesame internally

    :param guess: a dictionary with keys 'v', 'efn'
      and 'efp', each of which is an array of length sys.nx (for IV calculations) or
      (n_wavelengths, sys.nx) (for EQE calculations)
    :param sys: a Sesame Builder object
    """

    if guess is not None:

        v = guess["potential"] / sys.scaling.energy
        efn = guess["Efe"] / sys.scaling.energy
        efp = guess["Efh"] / sys.scaling.energy

        guess_sesame = {"v": v, "efn": efn, "efp": efp}

    else:
        guess_sesame = None

    return guess_sesame


@register_iv_solver("sesame_PDD")
def iv_sesame(junction, options):
    """Solve the dark or light IV using the Sesame solver. This will scan through the
    voltages  in options.internal_voltages and call sesame.IVcurve. If your calculation
    is failing to converge, make sure your calculation includes 0 V and try scanning
    more voltage points or setting a denser mesh. Note that the solver will not
    calculate anything at voltages > the bandgap + 10*kb*T, since it will fail to
    converge for the high injection regime.

    :param junction: a Junction object
    :param options: a State object containing options for the solver
    """

    solver_class = Solver()
    IVcurve = solver_class.IVcurve

    if not hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    if hasattr(junction, "guess"):
        # guess for the non-equilibrium solution at V = 0
        guess_sesame = process_guess(junction.guess, junction.sesame_sys)

    else:
        guess_sesame = None

    sesame_kwargs = process_sesame_options(options)

    user_defined_generation = options.user_defined_generation if (
            "user_defined_generation" in options) else False

    if options.light_iv and not user_defined_generation:
        gen_wl = junction.absorbed(junction.mesh) / 100  # m-1 -> cm-1

        if hasattr(junction, "layer_absorption"):
            A = np.trapezoid(
                np.nan_to_num(junction.absorbed(junction.mesh), nan=0.0), junction.mesh, axis=0
            )  # total absorption per wavelength using Sesame mesh
            profile_scale = junction.layer_absorption / A
            profile_scale[junction.layer_absorption < 1e-6] = 0

            gen_wl = gen_wl * profile_scale[None, :]

        else:
            warnings.warn(
                "layer_absorption of junction not provided, no generation"
                "profile scaling correction.",
                UserWarning,
            )

        wls = options.wavelength

        gg = (
            options.light_source.spectrum(wls, output_units="photon_flux_per_m")[1][
                :, None
            ]
            * gen_wl.T
        )

        g_vs_z = np.trapezoid(gg, wls, axis=0) / 1e4  # m^2 -> cm^2
        g_vs_z[np.isnan(g_vs_z)] = 0

        # can also pass a function to generation - more flexible?
        junction.sesame_sys.generation(g_vs_z)

    voltages = options.internal_voltages

    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, "R_shunt") else 1e14

    max_Eg = np.max(junction.sesame_sys.Eg * junction.sesame_sys.scaling.energy)
    max_V = (
        max_Eg + 10 * 8.617e-05 * options.T
    )  # do not go into the high injection regime, will not get convergence

    # voltages need to go from 0 (or close) to highest applied +ve or -ve voltage,
    # otherwise do not get convergence; need to go from V = 0 to high applied voltage
    # so that Sesame can use the previous solution as a guess for the next voltage.

    if junction.sesame_sys.rho[junction.sesame_sys.nx - 1] < 0:
        # this is necessary because Sesame will internally flip the sign for an n-p
        # junction
        voltages_for_solve = -voltages

        if np.all(options.voltages >= 0):
            warnings.warn(
                "All voltages are positive, but junction has been identified "
                "as n-p, so the  open-circuit voltage (Voc) of the junction will be "
                "negative.",
                UserWarning,
            )

    else:
        voltages_for_solve = voltages

    voltages_for_solve = voltages_for_solve[voltages_for_solve <= max_V]

    warnings.filterwarnings("ignore")
    # split +ve and -ve voltages if necessary:
    if np.any(voltages_for_solve < 0):
        if np.any(voltages_for_solve > 0):
            # positive and negative voltages

            negative_voltages = voltages_for_solve[voltages_for_solve <= 0]
            positive_voltages = voltages_for_solve[voltages_for_solve >= 0]

            negative_voltages_order = np.argsort(negative_voltages)[::-1]
            positive_voltages_order = np.argsort(positive_voltages)

            negative_voltages = negative_voltages[negative_voltages_order]
            positive_voltages = positive_voltages[positive_voltages_order]

            j_positive, result_positive = IVcurve(
                junction.sesame_sys,
                positive_voltages,
                guess=guess_sesame,
                **sesame_kwargs,
            )
            j_negative, result_negative = IVcurve(
                junction.sesame_sys,
                negative_voltages,
                guess=guess_sesame,
                **sesame_kwargs,
            )

            j_negative = j_negative[::-1]

            result_negative = {
                key: result_negative[key][::-1] for key in result_negative.keys()
            }

            negative_voltages = negative_voltages[::-1]

            if np.any(voltages_for_solve == 0):
                # V = 0 would have been included in both the +ve and -ve voltages, so
                # exclude it from the negative voltage results when concatenating

                j = np.concatenate((j_negative[:-1], j_positive))
                result = {
                    key: np.concatenate(
                        (result_negative[key][:-1], result_positive[key])
                    )
                    for key in result_positive.keys()
                }
                final_voltages = np.concatenate(
                    (negative_voltages[:-1], positive_voltages)
                )

            else:
                j = np.concatenate((j_negative, j_positive))
                result = {
                    key: np.concatenate((result_negative[key], result_positive[key]))
                    for key in result_positive.keys()
                }
                final_voltages = np.concatenate((negative_voltages, positive_voltages))

            # this results in j and result in order of increasing values for
            # voltages_for_solve
        else:
            # negative voltages only
            voltage_order = np.argsort(voltages_for_solve)[::-1]
            final_voltages = voltages_for_solve[voltage_order]
            j, result = IVcurve(
                junction.sesame_sys,
                final_voltages,
                guess=guess_sesame,
                **sesame_kwargs,
            )

    else:
        # positive voltages only
        voltage_order = np.argsort(voltages_for_solve)

        final_voltages = voltages_for_solve[voltage_order]
        j, result = IVcurve(
            junction.sesame_sys,
            final_voltages,
            guess=guess_sesame,
            **sesame_kwargs,
        )
    warnings.resetwarnings()

    if np.any(np.isnan(j)):
        warnings.warn(
            "Current calculation did not converge at all voltages", UserWarning
        )

    # final_voltages are the voltages corresponding to the entries in j and result,
    # using Sesame's sign convention. So if the voltage sign was flipped above, need
    # to flip it back for Solcore

    if junction.sesame_sys.rho[junction.sesame_sys.nx - 1] < 0:
        result_voltage = -final_voltages
        sort_result = np.argsort(result_voltage)
        j = j[sort_result]
        result = {key: result[key][sort_result, :] for key in result.keys()}
        result_voltage = result_voltage[sort_result]

    else:
        result_voltage = final_voltages

    j = j * junction.sesame_sys.scaling.current * 1e4  # cm-2 -> m-2

    shunted_current = j + result_voltage / R_shunt

    non_nans = np.where(~np.isnan(shunted_current))[0]

    if len(non_nans) > 0:
        # may be NaNs in j - find values closest to edge which are not NaN:
        first_non_nan = shunted_current[non_nans[0]]
        last_non_nan = shunted_current[non_nans[-1]]

    else:
        raise Exception(
            "No solutions found for IV curve. "
            "Try increasing the number of voltage points scanned."
        )

    junction.sesame_output = result

    junction.iv = interp1d(
        result_voltage[non_nans],
        shunted_current[non_nans],
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(first_non_nan, last_non_nan),
    )

    junction.voltage = options.internal_voltages
    junction.current = junction.iv(options.internal_voltages)
    junction.pdd_output = process_sesame_results(junction.sesame_sys, result)


def j_per_wl(
    system,
    solve,
    sesame_kwargs,
    guess=None,
):
    """Solve the Drift Diffusion Poisson equations at V=0.
    Parameters
    ----------
    system: Builder
        The discretized system.
    solve: solsesame.solvers.Solver.solve function
    sesame_kwargs: dictionary of keyword arguments to pass to Sesame's solve function,
        which can contain:
            - tol: float. Accepted error made by the Newton-Raphson scheme.
            - periodic_bcs: boolean. Defines the choice of boundary conditions in the
                y-direction. True
                (False) corresponds to periodic (abrupt) boundary conditions.
            - maxiter: integer. Maximum number of steps taken by the Newton-Raphson scheme.
            - verbose: boolean. The solver returns the step number and the associated
                error at every step, and this function prints the current applied
                voltage if set to True (default).
            - htp: integer
                Number of homotopic Newton loops to perform.

    guess: dictionary of numpy arrays of floats (optional)
        Starting point of the solver. Keys of the dictionary must be 'efn',
        'efp', 'v' for the electron and quasi-Fermi levels, and the
        electrostatic potential respectively.

    Returns
    -------
    J: numpy array of floats
        Steady state current computed for each voltage value.
    result: dictionary of numpy arrays of floats with the electron and hole Fermi levels
       (efn, efp) and the electrostatic potential (v) at each point in the mesh. Note
       this will be in Sesame's internal units.

    """

    # create a dictionary 'result' with efn and efp

    # Call the Drift Diffusion Poisson solver
    result = solve(
        system,
        guess=guess,
        **sesame_kwargs,
    )

    if result is not None:

        try:
            az = Analyzer(system, result)
            J = az.full_current()

        except Exception:
            J = np.nan
            result = {
                "efn": np.full(system.nx, np.nan),
                "efp": np.full(system.nx, np.nan),
                "v": np.full(system.nx, np.nan),
            }

    else:
        warnings.warn("The solver failed to converge.", UserWarning)
        J = np.nan
        result = {
            "efn": np.full(system.nx, np.nan),
            "efp": np.full(system.nx, np.nan),
            "v": np.full(system.nx, np.nan),
        }

    return J, result


def qe_sesame(junction: Junction, options: State):
    """Calculate the quantum efficiency of a junction using Sesame. This will scan
    through the wavelengths set in options.wavelength. It will scan from long
    wavelengths to short wavelengths, to improve the chance of convergence, since
    carrier generation profiles will be less steep at longer wavelengths.

    :param junction: a Junction object
    :param options: a State object containing options for the solver
    """

    def process_qe_sesame_options(options):

        """Process options for how Solcore interacts with the Sesame solver for quantum
        efficiency calculations."""

        if "sesame_use_previous_wl" in options:
            use_previous_wl = options.sesame_use_previous_wl

        else:
            use_previous_wl = True

        if "sesame_qe_flux" in options:

            if isinstance(options.sesame_qe_flux, (int, float)):
                flux = (
                    float(options.sesame_qe_flux)
                    * np.ones_like(options.wavelength)
                    / 1e4
                )

            elif isinstance(options.sesame_qe_flux, LightSource):
                flux = (
                    options.sesame_qe_flux.spectrum(
                        x=options.wavelength, output_units="photon_flux_per_m"
                    )[1]
                    / 1e4
                )
                # convert from m-2 -> cm-2

            elif isinstance(options.sesame_qe_flux, np.ndarray):
                if options.sesame_qe_flux.ndim == 1 and len(
                    options.sesame_qe_flux
                ) == len(options.wavelength):
                    # assume this is a 1D array of flux values
                    flux = options.sesame_qe_flux / 1e4

            else:
                raise ValueError(
                    "sesame_qe_flux must be a Solcore LightSource object,"
                    "a 1D array with the same length as options.wavelength,"
                    " or a float/int."
                )

        else:
            flux = 1e4 * np.ones_like(options.wavelength)

        if "sesame_qe_parallel" in options:
            parallel = options.sesame_qe_parallel

        else:
            parallel = False

        if "sesame_qe_n_jobs" in options:
            if isinstance(options.sesame_qe_n_jobs, int):
                n_jobs = options.sesame_qe_n_jobs

                if n_jobs > 1 or n_jobs == -1:
                    parallel = True
            else:
                raise ValueError("sesame_qe_n_jobs must be an integer.")

        else:
            n_jobs = -1

        if "sesame_qe_voltage" in options:
            if isinstance(options.sesame_qe_voltage, (int, float)):
                voltage = float(options.sesame_qe_voltage)
            else:
                raise ValueError("sesame_qe_voltage must be a float or int.")

        else:
            voltage = 0.0

        return flux, use_previous_wl, parallel, n_jobs, voltage

    flux, use_previous_wl, parallel, n_jobs, voltage = process_qe_sesame_options(
        options
    )

    solver_class = Solver()
    solve = solver_class.solve

    if not hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    if hasattr(junction, "guess"):
        # guess for the non-equilibrium solution at V = 0
        guess_sesame = process_guess(junction.guess, junction.sesame_sys)

    else:
        guess_sesame = None

    wls = options.wavelength

    profile_func = junction.absorbed

    sesame_kwargs = process_sesame_options(options)

    # the mesh used by Sesame and the mesh used in the optical solver are generally
    # not the same. This means the total generation (integrated over depth) may not
    # be the same, especially at short wavelengths where the absorption profile is very
    # steep. Calculate the mismatch:

    A = np.trapezoid(
        np.nan_to_num(junction.absorbed(junction.mesh), nan=0.0), junction.mesh, axis=0
    )  # total absorption per wavelength using Sesame mesh

    if hasattr(junction, "layer_absorption"):
        profile_scale = junction.layer_absorption / A

    else:
        # if no layer absorption is defined, do not scale
        profile_scale = np.ones_like(wls)
        warnings.warn(
            "layer_absorption of junction not provided, no generation"
            "profile scaling correction.",
            UserWarning,
        )

    # do not solve EQE if absorption is ~ 0

    EQE_threshold = 1e-5

    wl_solve = np.where(A >= EQE_threshold)[0][::-1]

    def make_gfcn_fun(wl_index, flux):
        def gcfn_fun(x, y):
            return (
                profile_scale[wl_index]
                * flux
                * np.nan_to_num(
                    profile_func(np.array([x / 100]))[0, wl_index] / 100, nan=0.0
                )
            )  # convert to cm-1 from m-1

        return gcfn_fun

    # more code for potential parallel implementation
    if parallel:

        def qe_i(system, i1, gen_profile):

            if guess_sesame is not None:
                # if there is a guess, use it as a starting point for the next
                # wavelength
                guess = {
                    "v": guess_sesame["v"][i1],
                    "efn": guess_sesame["efn"][i1],
                    "efp": guess_sesame["efp"][i1],
                }

            else:
                guess = None

            system.generation(gen_profile)

            j, result = j_per_wl(
                system,
                solve,
                sesame_kwargs,
                guess=guess,
            )

            eqe = np.abs(j) / (q * flux[i1])

            return eqe, result

        # make array for generation profile for each wavelength

        gen_wl_x = (
            profile_scale[:, None]
            * flux[:, None]
            * np.nan_to_num(profile_func(junction.mesh), nan=0.0).T
        ) / 100

        allres = Parallel(n_jobs=n_jobs)(
            delayed(qe_i)(junction.sesame_sys, i1, gen_profile=gen_wl_x[i1])
            for i1 in wl_solve[::-1]
        )

        eqe = np.array([item[0] for item in allres])

        efn_result = np.stack([item[1]["efn"] for item in allres])
        efp_result = np.stack([item[1]["efp"] for item in allres])
        v_result = np.stack([item[1]["v"] for item in allres])

    else:
        eqe = np.zeros_like(wls)

        # go in backwards order through wavelengths - since generation profile tends to
        # be flatter at longer wavelengths, this increases the change of convergence,
        # since the solution for the previous wavelength is always used as a guess for
        # the next wavelength. Having a good guess can help the short wavelength
        # solutions converge

        warnings.filterwarnings("ignore")
        # this is to prevent warnings from Sesame flooding the output. Not ideal but
        # unsure on best way to solve.

        efn_result = np.full((len(wl_solve), junction.sesame_sys.nx), np.nan)
        efp_result = np.full((len(wl_solve), junction.sesame_sys.nx), np.nan)
        v_result = np.full((len(wl_solve), junction.sesame_sys.nx), np.nan)

        result = None

        for i1 in wl_solve:
            junction.sesame_sys.generation(make_gfcn_fun(i1, flux[i1]))

            # if i1 == wl_solve[0] and guess_sesame is not None:
            if use_previous_wl and i1 > 0 and result is not None:
                guess = result

            elif guess_sesame is not None:
                # if there is a guess, use it as a starting point for the next
                # wavelength
                guess = {
                    "v": guess_sesame["v"][i1],
                    "efn": guess_sesame["efn"][i1],
                    "efp": guess_sesame["efp"][i1],
                }

            else:
                guess = None

            j, result = j_per_wl(
                junction.sesame_sys,
                solve,
                sesame_kwargs,
                guess=guess,
            )

            eqe[i1] = np.abs(j) / (q * flux[i1])

            # if not np.isnan(j):
            efn_result[i1] = result["efn"]
            efp_result[i1] = result["efp"]
            v_result[i1] = result["v"]

    if np.any(np.isnan(eqe)):
        warnings.warn(
            "EQE calculation did not converge at all wavelengths", UserWarning
        )

    warnings.resetwarnings()

    eqe = eqe * junction.sesame_sys.scaling.current
    iqe = np.divide(eqe, A, out=np.zeros_like(eqe), where=A > 0)

    # line = ((0, 0), (np.max(junction.mesh_cm), 0))
    # scale_sesame_result(junction.sesame_sys, result, line)

    # convert dimensionless current to dimension-ful current

    junction.iqe = interp1d(wls, iqe)

    junction.eqe = interp1d(
        wls,
        eqe,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(eqe[0], eqe[-1]),
    )

    # PDD output for process_sesame_results should be dictionary with keys 'efn', 'efp',
    # 'v', each entry has shape (n_wavelengths, mesh_points)
    result_dict = {
        "efn": efn_result,
        "efp": efp_result,
        "v": v_result,
    }

    pdd_output = process_sesame_results(junction.sesame_sys, result_dict)

    junction.qe = State(
        {
            "WL": wls,
            "IQE": junction.iqe(wls),
            "EQE": junction.eqe(wls),
            "pdd_output": pdd_output,
        }
    )


def process_sesame_results(sys: Builder, result: dict):
    r"""
        Scale the result of a Sesame calculation to SI units, and calculate other
        quantities like the positions of the conduction and valence band and the
        recombination rates. Produces a State object with entries:

            - potential: the potential (V)
            - n: the electron density (m\ :sup:`-3`)
            - p: the hole density (m\ :sup:`-3`)
            - Ec: the level of the conduction band (eV)
            - Ev the level of the valence band (eV)
            - Efe: the electron quasi-Fermi level (eV)
            - Efh the hole quasi-Fermi level (eV)
            - Rrad: the radiative recombination rate (m\ :sup:`-3` s\ :sup:`-1`)
            - Raug: the Auger recombination rate (m\ :sup:`-3` s\ :sup:`-1`)
            - Rsrh: the bulk recombination due to
              Shockley-Read-Hall processes (m\ :sup:`-3` s\ :sup:`-1`)

    Each of these is a 2-dimensional array, with dimensions
    ``(len(options.internal_voltages), len(mesh))``.

        :param sys: a Sesame Builder object
        :param result: a dictionary containing the results from a Sesame calculation
    """

    line = ((0, 0), (np.max(sys.xpts), 0))
    n_voltages = len(result["v"])

    generation = sys.g * sys.scaling.generation * 1e6  # multiply by internal
    # sesame scaling factor,
    # convert from cm-3 to m-3 (area and depth are in m, in sesame they are in cm)

    potential = result["v"] * sys.scaling.energy
    Efe = result["efn"] * sys.scaling.energy
    Efh = result["efp"] * sys.scaling.energy
    Ec = -(result["v"] + sys.bl) * sys.scaling.energy
    Ev = -(result["v"] + sys.bl + sys.Eg) * sys.scaling.energy

    n = np.zeros((n_voltages, sys.nx))
    p = np.zeros((n_voltages, sys.nx))

    Rrad = np.zeros((n_voltages, sys.nx))
    Raug = np.zeros((n_voltages, sys.nx))
    Rsrh = np.zeros((n_voltages, sys.nx))

    for i1 in range(n_voltages):
        result_loop = {key: result[key][i1, :] for key in result.keys()}

        analyzer = Analyzer(sys, result_loop)

        n[i1] = (
            analyzer.electron_density(location=line) * sys.scaling.density * 1e6
        )  # m-3
        p[i1] = analyzer.hole_density(location=line) * sys.scaling.density * 1e6  # m-3

        Rsrh[i1] = (
            analyzer.bulk_srh_rr(location=line) * sys.scaling.generation * 1e6
        )  # m-3
        Raug[i1] = (
            analyzer.auger_rr(location=line) * sys.scaling.generation * 1e6
        )  # m-3
        Rrad[i1] = (
            analyzer.radiative_rr(location=line) * sys.scaling.generation * 1e6
        )  # m-3

    output = State(
        G=generation,
        potential=potential,
        Efe=Efe,
        Efh=Efh,
        Ec=Ec,
        Ev=Ev,
        n=n,
        p=p,
        Rrad=Rrad,
        Raug=Raug,
        Rsrh=Rsrh,
    )

    return output
