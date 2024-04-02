from __future__ import annotations

from sesame import Builder, IVcurve, Analyzer
import numpy as np
from scipy.optimize import root
from solcore.constants import q, kb
from scipy.interpolate import interp1d
from solcore.state import State
from solcore.structure import Junction, Layer
from solcore.sesame_drift_diffusion.process_structure import process_structure
import warnings

from solcore.registries import register_iv_solver, register_equilibrium_solver


@register_equilibrium_solver("sesame_PDD")
def equilibrium(junction: Junction, **kwargs):
    """
    Solve at equilibrium (no illumination, no applied voltage) using the Sesame solver.

    :param junction: a Junction object
    :param options: a State object containing options for the solver
    """

    # TODO: pass output from 'result' parameter to the user, similar to the Fortran PDD solver, so that
    # band structures etc. can be plotted.

    options = State(
        kwargs
    )  # needs to be passed as kwargs to be compatible with Fortran equilibrium solver
    # in registries

    if not hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    j, result = IVcurve(junction.sesame_sys, [0])

    if np.any(np.isnan(j)):
        warnings.warn(
            "Current calculation did not converge at all voltages", UserWarning
        )

    j = j * junction.sesame_sys.scaling.current * 1e4  # cm-2 -> m-2

    junction.sesame_output = result


@register_iv_solver("sesame_PDD")
def iv_sesame(junction, options):
    """
    Solve the dark or light IV using the Sesame solver. This will scan through the voltages
    in options.internal_voltages and call sesame.IVcurve. If your calculation is failing to converge,
    make sure your calculation includes 0 V and try scanning more voltage points or setting a denser
    mesh. Note that the solver will not calculate anything at voltages > the bandgap + 10*kb*T, since
    it will fail to converge for the high injection regime.

    :param junction: a Junction object
    :param options: a State object containing options for the solver
    """

    # TODO: pass output from 'result' parameter to the user, similar to the Fortran PDD solver, so that
    # band structures etc. can be plotted.

    if not hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    if options.light_iv:  # and np.all(junction.sesame_sys.g == 0):
        gen_wl = junction.absorbed(junction.mesh) / 100  # m-1 -> cm-1
        wls = options.wavelength

        gg = (
            options.light_source.spectrum(wls, output_units="photon_flux_per_m")[1][
                :, None
            ]
            * gen_wl.T
        )

        g_vs_z = np.trapz(gg, wls, axis=0) / 1e4  # m^2 -> cm^2
        g_vs_z[np.isnan(g_vs_z)] = 0

        # can also pass a function to generation - more flexible?
        junction.sesame_sys.generation(g_vs_z)

    voltages = options.internal_voltages

    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, "R_shunt") else 1e14

    max_Eg = np.max(junction.sesame_sys.Eg * junction.sesame_sys.scaling.energy)
    max_V = (
        max_Eg + 10 * 8.617e-05 * options.T
    )  # do not go into the high injection regime, will not get convergence

    # voltages need to go from 0 (or close) to highest applied +ve or -ve voltage, otherwise
    # do not get convergence; need to go from V = 0 to high applied voltage so that Sesame can
    # use the previous solution as a guess for the next voltage.

    if junction.sesame_sys.rho[junction.sesame_sys.nx - 1] < 0:
        # this is necessary because Sesame will internally flip the sign for an n-p junction
        voltages_for_solve = -voltages

        if np.all(options.voltages >= 0):
            warnings.warn(
                "All voltages are positive, but junction has been identified as n-p, so the "
                "open-circuit voltage (Voc) of the junction will be negative.",
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
                junction.sesame_sys, positive_voltages
            )
            j_negative, result_negative = IVcurve(
                junction.sesame_sys, negative_voltages
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

            # this results in j and result in order of increasing values for voltages_for_solve
        else:
            # negative voltages only
            voltage_order = np.argsort(voltages_for_solve)[::-1]
            final_voltages = voltages_for_solve[voltage_order]
            j, result = IVcurve(junction.sesame_sys, final_voltages)

    else:
        # positive voltages only
        voltage_order = np.argsort(voltages_for_solve)

        final_voltages = voltages_for_solve[voltage_order]
        j, result = IVcurve(junction.sesame_sys, final_voltages)  # , verbose=False)

    warnings.resetwarnings()

    if np.any(np.isnan(j)):
        warnings.warn(
            "Current calculation did not converge at all voltages", UserWarning
        )

    # final_voltages are the voltages corresponding to the entries in j and result, USING
    # Sesame's sign convention. So if the voltage sign was flipped above, need to flip it back
    # for Solcore

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
        Exception(
            "No solutions found for IV curve. Try increasing the number of voltage points scanned."
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


def qe_sesame(junction: Junction, options: State):
    """
    Calculate the quantum efficiency of a junction using Sesame. This will scan through the wavelengths
    set in options.wavelength. It will scan from long wavelengths to short wavelengths, to improve
    the chance of convergence, since carrier generation profiles will be less steep at longer wavelengths.

    :param junction: a Junction object
    :param options: a State object containing options for the solver
    """
    if not hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    # if not options.parallel:
    #     n_jobs = 1
    #
    # else:
    #     n_jobs = options.n_jobs if hasattr(options, "n_jobs") else -1
    # n_jobs = 1
    # parallel implementation does not work due to pickling error. Unsure of cause.

    wls = options.wavelength

    voltages = [0]

    # profile_func = interp1d(bulk_positions_cm, 1e7 * Si_profile, kind='linear', bounds_error=False, fill_value=0)
    profile_func = (
        junction.absorbed
    )  # this returns an array of shape (mesh_points, wavelengths)

    A = np.trapz(junction.absorbed(junction.mesh), junction.mesh, axis=0)

    def make_gfcn_fun(wl_index, flux):
        def gcfn_fun(x, y):
            return (
                flux * profile_func(np.array([x / 100]))[0, wl_index] / 100
            )  # convert to cm-1 from m-1

        return gcfn_fun

    # def qe_i(i1):
    #
    #     junction.sesame_sys.generation(make_gfcn_fun(i1, flux))
    #
    #     j, result = IVcurve(junction.sesame_sys, voltages)
    #
    #     eqe = np.abs(j) / (q * flux)
    #
    #     return eqe, result
    #
    # allres = Parallel(n_jobs=n_jobs)(
    #     delayed(qe_i)(
    #         i1,
    #     )
    #     for i1 in range(len(wls))
    # )
    #
    # eqe = np.concatenate([item[0] for item in allres])
    #
    # efn = np.concatenate([item[1]['efn'] for item in allres])
    # efp = np.concatenate([item[1]['efp'] for item in allres])
    # vres = np.concatenate([item[1]['v'] for item in allres])

    flux = 1e20

    eqe = np.zeros_like(wls)

    # do not solve EQE if absorption is ~ 0

    EQE_threshold = 1e-5

    wl_solve = np.where(A >= EQE_threshold)[0][::-1]
    # go in backwards order through wavelengths - since generation profile tends to be flatter at longer wavelengths,
    # this increases the change of convergence, since the solution for the previous wavelength is always used as a guess
    # for the next wavelength. Gaving a good guess helps the short wavelength solutions converge

    warnings.filterwarnings("ignore")
    # this is to prevent warning from Sesame flooding the output. Not ideal but unsure on best way to solve.

    for i1 in wl_solve:
        junction.sesame_sys.generation(make_gfcn_fun(i1, flux))

        if i1 == wl_solve[0]:
            guess = None

        else:
            guess = {key: result[key][0, :] for key in result.keys()}

        j, result = IVcurve(junction.sesame_sys, voltages, guess=guess)

        eqe[i1] = np.abs(j) / (q * flux)

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

    junction.qe = State({"WL": wls, "IQE": junction.iqe(wls), "EQE": junction.eqe(wls)})


def process_sesame_results(sys: Builder, result: dict):
    """
        Scale the result of a Sesame calculation to SI units, and calculate other quantities like
        the positions of the conduction and valence band and the recombination rates. Produces a
        State object with entries:

            - potential: the potential (V)
            - n: the electron density (m\ :sup:`-3`)
            - p: the hole density (m\ :sup:`-3`)
            - Ec: the level of the conduction band (eV)
            - Ev the level of the valence band (eV)
            - Efe: the electron quasi-Fermi level (eV)
            - Efh the hole quasi-Fermi level (eV)
            - Rrad: the radiative recombination rate (m\ :sup:`-3` s\ :sup:`-1`)
            - Raug: the Auger recombination rate (m\ :sup:`-3` s\ :sup:`-1`)
            - Rsrh: the bulk recombination due to Shockley-Read-Hall processes (m\ :sup:`-3` s\ :sup:`-1`)

    Each of these is a 2-dimensional array, with dimensions ``(len(options.internal_voltages), len(mesh))``.


        :param sys: a Sesame Builder object
        :param result: a dictionary containing the results from a Sesame calculation
    """

    line = ((0, 0), (np.max(sys.xpts), 0))
    n_voltages = len(result["v"])

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
