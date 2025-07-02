import numpy as np
from pytest import approx, raises


def test_constant_doping():
    # test sesame PDD process_structure when constant doping levels are passed
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.sesame_drift_diffusion.process_structure import process_structure
    from solcore.state import State

    Si_n = material("Si")(
        Nd=1e24, electron_minority_lifetime=1e-6, hole_minority_lifetime=1e-7
    )
    Si_i = material("Si")(
        Na=0, Nd=0, electron_minority_lifetime=2e-6, hole_minority_lifetime=2e-7
    )
    Si_p = material("Si")(
        Na=2e24, electron_minority_lifetime=1.5e-6, hole_minority_lifetime=1.5e-7
    )

    junction = Junction(
        [Layer(si("200nm"), Si_n), Layer(si("2000nm"), Si_i), Layer(si("2000nm"), Si_p)]
    )

    options = State()
    options.T = 300

    process_structure(junction, options)

    sesame_obj = junction.sesame_sys

    assert len(np.unique(sesame_obj.rho)) == 3
    assert sesame_obj.rho[0] == 1e24 * 1e-6 / sesame_obj.scaling.density
    assert sesame_obj.rho[-1] == -2e24 * 1e-6 / sesame_obj.scaling.density


def test_doping_profile():
    # test process_structure when a doping profile is passed for one of the layers,
    # or to the whole junction
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.sesame_drift_diffusion.process_structure import process_structure
    from solcore.state import State
    from scipy.interpolate import interp1d

    Si_n = material("Si")(
        Nd=1e25, electron_minority_lifetime=1e-6, hole_minority_lifetime=1e-7
    )
    Si_i = material("Si")(electron_minority_lifetime=2e-6, hole_minority_lifetime=2e-7)
    Si_p = material("Si")(
        Na=2e24, electron_minority_lifetime=1.5e-6, hole_minority_lifetime=1.5e-7
    )

    doping_profile = np.linspace(Si_n.Nd, -Si_p.Na, 2000)
    x_pos = np.linspace(0, 2000, 2000) * 1e-9
    doping_profile_func = interp1d(
        x_pos, doping_profile, fill_value=(Si_n.Nd, -Si_p.Na), bounds_error=False
    )

    junction = Junction(
        [
            Layer(si("200nm"), Si_n),
            Layer(si("2000nm"), Si_i, doping_profile=doping_profile_func),
            Layer(si("2000nm"), Si_p),
        ]
    )

    options = State()
    options.T = 300

    process_structure(junction, options)

    rho_1 = interp1d(junction.mesh, junction.sesame_sys.rho)

    x_pos_2 = np.linspace(200, 2200, 2000) * 1e-9
    doping_profile_func = interp1d(
        x_pos_2, doping_profile, fill_value=(Si_n.Nd, -Si_p.Na), bounds_error=False
    )

    junction = Junction(
        [
            Layer(si("200nm"), Si_n),
            Layer(si("2000nm"), Si_i),
            Layer(si("2000nm"), Si_p),
        ],
        doping_profile=doping_profile_func,
    )

    process_structure(junction, options)

    rho_2 = interp1d(junction.mesh, junction.sesame_sys.rho)

    assert np.allclose(rho_1(junction.mesh), rho_2(junction.mesh))


def test_parameter_extraction():
    # test all material parameters are correctly extracted by process_structure
    from solcore.sesame_drift_diffusion.process_structure import get_material_parameters
    from solcore import material
    from solcore.constants import q, kb

    test_mat = material("Si")(Nd=1e24, Na=2e24)

    with raises(ValueError) as excinfo:
        get_material_parameters(test_mat)
    assert excinfo.match("nor in calculable parameters")

    test_mat.hole_diffusion_length, test_mat.electron_diffusion_length = 100e-9, 200e-9

    expected_tau_e = (
        q
        * test_mat.electron_diffusion_length**2
        / (kb * test_mat.T * test_mat.electron_mobility)
    )
    expected_tau_h = (
        q
        * test_mat.hole_diffusion_length**2
        / (kb * test_mat.T * test_mat.hole_mobility)
    )

    result = get_material_parameters(test_mat)

    expected_dict = {
        "Nc": test_mat.Nc * 1e-6,
        "Nv": test_mat.Nv * 1e-6,
        "Eg": test_mat.band_gap / q,
        "affinity": test_mat.electron_affinity / q,
        "epsilon": test_mat.relative_permittivity,
        "mu_e": test_mat.electron_mobility * 1e4,
        "mu_h": test_mat.hole_mobility * 1e4,
        "tau_e": expected_tau_e,
        "tau_h": expected_tau_h,  # hole bulk lifetime (s)
        "Et": 0,  # energy level of bulk recombination centres (eV)
        "B": test_mat.radiative_recombination * 1e6,
        # radiative recombination constant (cm3/s)
        "Cn": test_mat.electron_auger_recombination * 1e12,
        # Auger recombination constant for electrons (cm6/s)
        "Cp": test_mat.hole_auger_recombination * 1e12,
        # Auger recombination constant for holes (cm6/s)
    }

    assert result == expected_dict

    test_mat.electron_minority_lifetime = 1e-6
    test_mat.hole_minority_lifetime = 1e-5

    expected_dict["tau_e"] = 1e-6
    expected_dict["tau_h"] = 1e-5

    result = get_material_parameters(test_mat)

    assert result == expected_dict

    test_mat.bulk_recombination_energy = 0.5

    expected_dict["Et"] = test_mat.bulk_recombination_energy / q

    result = get_material_parameters(test_mat)

    assert result == expected_dict


def test_compare_DA_np():

    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.analytic_solar_cells import iv_depletion
    from solcore.sesame_drift_diffusion.solve_pdd import iv_sesame
    from solcore.sesame_drift_diffusion.process_structure import carrier_constants
    from solcore.state import State
    from solcore.light_source import LightSource

    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    GaAs_p.electron_diffusion_length = carrier_constants(
        "electron_diffusion_length", GaAs_p
    )
    GaAs_n.hole_diffusion_length = carrier_constants("hole_diffusion_length", GaAs_n)

    options = State()
    options.wavelength = np.linspace(300, 950, 100) * 1e-9
    options.voltages = np.linspace(-1.3, 0.5, 30)
    options.internal_voltages = np.linspace(-1.3, 0.5, 30)
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(
        source_type="standard",
        version="AM1.5g",
        x=options.wavelength,
        output_units="photon_flux_per_m",
    )
    options.da_mode = "green"
    options.optics_method = "TMM"

    mesh = np.linspace(0, 2200, 500) * 1e-9

    np_junction = Junction(
        [
            Layer(si("200nm"), GaAs_n, role="emitter"),
            Layer(si("2000nm"), GaAs_p, role="base"),
        ],
        kind="DA",
        R_shunt=0.1,
        mesh=mesh,
        sn=10,
        sp=10,
    )

    solar_cell_solver(SolarCell([np_junction]), "optics", options)

    iv_depletion(np_junction, options)
    depletion_current = np_junction.current
    depletion_current_interp = np_junction.iv(options.voltages)

    iv_sesame(np_junction, options)
    sesame_current = np_junction.current
    sesame_current_interp = np_junction.iv(options.voltages)

    assert depletion_current[-1] == approx(sesame_current[-1], rel=0.05)
    assert np.sign(depletion_current[0]) == np.sign(sesame_current[0])

    assert depletion_current_interp[-1] == approx(sesame_current_interp[-1], rel=0.05)
    assert np.sign(depletion_current_interp[0]) == np.sign(sesame_current_interp[0])

    assert np.all(sesame_current == sesame_current_interp)


def test_compare_DA_pn():

    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.analytic_solar_cells import iv_depletion
    from solcore.sesame_drift_diffusion.solve_pdd import iv_sesame
    from solcore.sesame_drift_diffusion.process_structure import carrier_constants
    from solcore.state import State
    from solcore.light_source import LightSource

    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    GaAs_p.electron_diffusion_length = carrier_constants(
        "electron_diffusion_length", GaAs_p
    )
    GaAs_n.hole_diffusion_length = carrier_constants("hole_diffusion_length", GaAs_n)

    options = State()
    options.wavelength = np.linspace(300, 950, 100) * 1e-9
    options.voltages = np.linspace(-0.5, 1.3, 30)
    options.internal_voltages = np.linspace(-0.5, 1.3, 30)
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(
        source_type="standard",
        version="AM1.5g",
        x=options.wavelength,
        output_units="photon_flux_per_m",
    )
    options.da_mode = "green"
    options.optics_method = "TMM"

    mesh = np.linspace(0, 2200, 500) * 1e-9

    pn_junction = Junction(
        [
            Layer(si("200nm"), GaAs_p, role="emitter"),
            Layer(si("2000nm"), GaAs_n, role="base"),
        ],
        kind="DA",
        R_shunt=0.1,
        mesh=mesh,
    )

    solar_cell_solver(SolarCell([pn_junction]), "optics", options)

    iv_depletion(pn_junction, options)
    depletion_current = pn_junction.current
    depletion_current_interp = pn_junction.iv(options.voltages)

    iv_sesame(pn_junction, options)
    sesame_current = pn_junction.current
    sesame_current_interp = pn_junction.iv(options.voltages)

    assert depletion_current[0] == approx(sesame_current[0], rel=0.05)
    assert np.sign(depletion_current[-1]) == np.sign(sesame_current[-1])

    assert depletion_current_interp[0] == approx(sesame_current_interp[0], rel=0.05)
    assert np.sign(depletion_current_interp[-1]) == np.sign(sesame_current_interp[-1])

    assert np.all(sesame_current == sesame_current_interp)


def test_qe():

    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.analytic_solar_cells import qe_depletion
    from solcore.sesame_drift_diffusion.solve_pdd import qe_sesame
    from solcore.sesame_drift_diffusion.process_structure import carrier_constants
    from solcore.state import State

    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    GaAs_p.electron_diffusion_length = carrier_constants(
        "electron_diffusion_length", GaAs_p
    )
    GaAs_n.hole_diffusion_length = carrier_constants("hole_diffusion_length", GaAs_n)

    wls = np.linspace(300, 950, 100) * 1e-9
    options = State()
    options.wavelength = wls
    options.T = 300
    options.light_iv = True
    options.da_mode = "green"
    options.optics_method = "BL"

    mesh = np.linspace(0, 2500, 1000) * 1e-9

    pn_junction = Junction(
        [
            Layer(si("500nm"), GaAs_p, role="emitter"),
            Layer(si("2000nm"), GaAs_n, role="base"),
        ],
        kind="DA",
        R_shunt=0.1,
        mesh=mesh,
        sn=1e5,
        sp=1e5,
    )

    solar_cell_solver(SolarCell([pn_junction]), "optics", options)

    qe_sesame(pn_junction, options)

    sesame_EQE = pn_junction.eqe(wls)
    sesame_IQE = pn_junction.iqe(wls)

    qe_depletion(pn_junction, options)
    depletion_EQE = pn_junction.eqe(wls)
    depletion_IQE = pn_junction.iqe(wls)

    assert sesame_EQE == approx(depletion_EQE, rel=0.1)
    assert sesame_IQE == approx(depletion_IQE, rel=0.1)


def test_mesh_generation():
    from solcore.sesame_drift_diffusion.process_structure import make_mesh
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.state import State

    GaAs_n = material("GaAs")()
    GaAs_p = material("GaAs")()

    options = State(minimum_spacing=1e-9, maximum_spacing=1e-9)

    pn_junction = Junction([Layer(si("500nm"), GaAs_n), Layer(si("2000nm"), GaAs_p)])

    layer_width = [layer.width * 1e2 for layer in pn_junction]

    make_mesh(pn_junction, layer_width, options, [500e-7])

    assert pn_junction.mesh == approx(np.arange(0, 2500.01, 1) * 1e-9)


def test_get_srv_np(np_junction):
    from solcore.sesame_drift_diffusion.process_structure import (
        get_srv,
        process_structure,
    )
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.state import State

    options = State(T=300)

    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    junction = Junction(
        [
            Layer(si("200nm"), GaAs_n, role="emitter"),
            Layer(si("2000nm"), GaAs_p, role="base"),
        ]
    )

    junction.sn = 5
    junction.sp = 8

    process_structure(junction, options)

    Sfront_e, Sfront_h, Sback_e, Sback_h = get_srv(junction)

    assert Sfront_e == junction.sn * 100
    assert Sfront_h == junction.sn * 100

    assert Sback_e == junction.sp * 100
    assert Sback_e == junction.sp * 100

    junction.sn_e = 2
    junction.sn_h = 3
    junction.sp_e = 4
    junction.sp_h = 5

    Sfront_e, Sfront_h, Sback_e, Sback_h = get_srv(junction)

    assert Sfront_e == junction.sn_e * 100
    assert Sfront_h == junction.sn_h * 100

    assert Sback_e == junction.sp_e * 100
    assert Sback_h == junction.sp_h * 100


def test_get_srv_pn(np_junction):
    from solcore.sesame_drift_diffusion.process_structure import (
        get_srv,
        process_structure,
    )
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.state import State

    options = State(T=300)

    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    junction = Junction(
        [
            Layer(si("200nm"), GaAs_p, role="emitter"),
            Layer(si("2000nm"), GaAs_n, role="base"),
        ]
    )

    junction.sn = 5
    junction.sp = 8

    process_structure(junction, options)

    Sfront_e, Sfront_h, Sback_e, Sback_h = get_srv(junction)

    assert Sfront_e == junction.sp * 100
    assert Sfront_h == junction.sp * 100

    assert Sback_e == junction.sn * 100
    assert Sback_e == junction.sn * 100

    junction.sn_e = 2
    junction.sn_h = 3
    junction.sp_e = 4
    junction.sp_h = 5

    Sfront_e, Sfront_h, Sback_e, Sback_h = get_srv(junction)

    assert Sfront_e == junction.sp_e * 100
    assert Sfront_h == junction.sp_h * 100

    assert Sback_e == junction.sn_e * 100
    assert Sback_h == junction.sn_h * 100


def test_voltages_np():
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.sesame_drift_diffusion import iv_sesame
    from solcore.state import State
    from solcore.light_source import LightSource

    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    options = State()
    options.wavelength = np.linspace(300, 950, 100) * 1e-9
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(
        source_type="standard",
        version="AM1.5g",
        x=options.wavelength,
        output_units="photon_flux_per_m",
    )
    options.optics_method = "TMM"

    np_junction = Junction(
        [
            Layer(si("200nm"), GaAs_n, role="emitter"),
            Layer(si("2000nm"), GaAs_p, role="base"),
        ],
        kind="sesame_PDD",
    )

    solar_cell_solver(SolarCell([np_junction]), "optics", options)

    voltage_end = [0, 1.3, 1.3]
    voltage_points = [20, 41, 40]

    interp_voltages = np.linspace(-1, 0, 10)

    interp_results = np.zeros((len(voltage_end), len(interp_voltages)))

    for i, V_end in enumerate(voltage_end):
        options.voltages = np.linspace(-1.3, V_end, voltage_points[i])
        options.internal_voltages = np.linspace(-1.3, V_end, voltage_points[i])

        iv_sesame(np_junction, options)

        interp_results[i] = np_junction.iv(interp_voltages)

    assert interp_results[0] == approx(interp_results[1], rel=0.02)
    assert interp_results[1] == approx(interp_results[2], rel=0.02)


def test_voltages_pn():
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.sesame_drift_diffusion import iv_sesame
    from solcore.state import State
    from solcore.light_source import LightSource

    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    options = State()
    options.wavelength = np.linspace(300, 950, 100) * 1e-9
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(
        source_type="standard",
        version="AM1.5g",
        x=options.wavelength,
        output_units="photon_flux_per_m",
    )
    options.optics_method = "TMM"

    pn_junction = Junction(
        [
            Layer(si("200nm"), GaAs_p, role="emitter"),
            Layer(si("2000nm"), GaAs_n, role="base"),
        ],
        kind="sesame_PDD",
    )

    solar_cell_solver(SolarCell([pn_junction]), "optics", options)

    voltage_start = [0, -1.3, -1.3]
    voltage_points = [20, 41, 40]

    interp_voltages = np.linspace(0, 1, 10)

    interp_results = np.zeros((len(voltage_start), len(interp_voltages)))

    for i, V_start in enumerate(voltage_start):
        options.voltages = np.linspace(V_start, 1.3, voltage_points[i])
        options.internal_voltages = np.linspace(V_start, 1.3, voltage_points[i])

        iv_sesame(pn_junction, options)

        interp_results[i] = pn_junction.iv(interp_voltages)

    assert interp_results[0] == approx(interp_results[1], rel=0.02)
    assert interp_results[1] == approx(interp_results[2], rel=0.02)


def test_reverse_bias():
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.sesame_drift_diffusion import iv_sesame
    from solcore.state import State
    from solcore.light_source import LightSource

    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    options = State()
    options.wavelength = np.linspace(300, 950, 100) * 1e-9
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(
        source_type="standard",
        version="AM1.5g",
        x=options.wavelength,
        output_units="photon_flux_per_m",
    )
    options.optics_method = "TMM"

    np_junction = Junction(
        [
            Layer(si("200nm"), GaAs_n, role="emitter"),
            Layer(si("2000nm"), GaAs_p, role="base"),
        ],
        kind="sesame_PDD",
    )

    solar_cell_solver(SolarCell([np_junction]), "optics", options)
    options.voltages = np.linspace(0, 1.3, 10)
    options.internal_voltages = options.voltages

    iv_sesame(np_junction, options)

    # this is reverse bias for an n-p junction, so check that all currents are > 0:

    assert np.all(np_junction.current > 0)
    assert np.all(np_junction.iv(options.voltages) > 0)


def test_generation_rate():

    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.sesame_drift_diffusion.process_structure import carrier_constants
    from solcore.state import State
    from solcore.light_source import LightSource

    GaAs_p = material("GaAs")(
        T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )
    GaAs_n = material("GaAs")(
        T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6
    )

    GaAs_p.electron_diffusion_length = carrier_constants(
        "electron_diffusion_length", GaAs_p
    )
    GaAs_n.hole_diffusion_length = carrier_constants("hole_diffusion_length", GaAs_n)

    options = State()
    options.wavelength = np.linspace(300, 950, 100) * 1e-9
    options.voltages = np.linspace(-0.5, 1.3, 30)
    options.internal_voltages = np.linspace(-0.5, 1.3, 30)
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(
        source_type="standard",
        version="AM1.5g",
        x=options.wavelength,
        output_units="photon_flux_per_m",
    )
    options.da_mode = "green"
    options.optics_method = "TMM"

    mesh = np.linspace(0, 2200, 500) * 1e-9

    pn_junction = Junction(
        [
            Layer(si("200nm"), GaAs_p, role="emitter"),
            Layer(si("2000nm"), GaAs_n, role="base"),
        ],
        kind="sesame_PDD",
        R_shunt=0.1,
        mesh=mesh,
    )
    solar_cell = SolarCell([pn_junction])
    solar_cell_solver(solar_cell, "iv", options)

    z_pos = np.linspace(0, 2200, 500) * 1e-9

    gen_per_wl = solar_cell(0).absorbed(z_pos)

    # correction for layer absorption:

    A = np.trapezoid(
        np.nan_to_num(pn_junction.absorbed(pn_junction.mesh), nan=0.0),
        pn_junction.mesh, axis=0
    )  # total absorption per wavelength using Sesame mesh
    profile_scale = pn_junction.layer_absorption / A
    profile_scale[pn_junction.layer_absorption < 1e-6] = 0

    gen_total = np.trapezoid(
        profile_scale[None, :] * gen_per_wl
        * options.light_source.spectrum(
            options.wavelength, output_units="photon_flux_per_m"
        )[1][None, :],
        options.wavelength,
    )

    gen_G = solar_cell(0).pdd_output.G
    assert gen_total == approx(gen_G)


def test_process_sesame_options():
    from solcore.sesame_drift_diffusion.solve_pdd import process_sesame_options
    from solcore.state import State

    options = State(
        sesame_max_iterations=1000,
        sesame_tol=1e-4,
        sesame_verbose=True,
        sesame_htp=2,
        sesame_periodic=False,
    )

    dict_out = process_sesame_options(options)

    assert dict_out["maxiter"] == 1000
    assert dict_out["tol"] == 1e-4
    assert dict_out["verbose"] is True
    assert dict_out["htp"] == 2
    assert dict_out["periodic_bcs"] is False


def test_process_guess():
    from solcore.sesame_drift_diffusion.solve_pdd import process_guess
    from solcore.state import State

    sys_mock = State(scaling=State(energy=2))

    guess = process_guess(None, sys_mock)

    assert guess is None

    v_guess = np.random.rand(3, 5)
    efn_guess = np.random.rand(3, 5)
    efp_guess = np.random.rand(3, 5)

    guess_in = {"potential": v_guess, "Efe": efn_guess, "Efh": efp_guess}

    guess_out = process_guess(guess_in, sys_mock)

    assert np.all(guess_out["v"] == v_guess / sys_mock.scaling.energy)
    assert np.all(guess_out["efn"] == efn_guess / sys_mock.scaling.energy)
    assert np.all(guess_out["efp"] == efp_guess / sys_mock.scaling.energy)
