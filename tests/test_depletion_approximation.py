from pytest import approx, mark, raises
import numpy as np
from solcore.constants import kb, q, vacuum_permittivity



def test_get_j_dark():
    from solcore.analytic_solar_cells.depletion_approximation import get_j_dark

    x = np.power(10, np.random.uniform(-8, -5))
    xi = np.power(10, np.random.uniform(-8, -5))

    l = np.power(10, np.random.uniform(-9, -6))
    s = np.power(10, np.random.uniform(1, 3))
    d = np.power(10, np.random.uniform(-5, 0))
    Vbi = np.random.uniform(0.1,5)
    minor = np.power(10, np.random.uniform(-7, -4))
    T = np.random.uniform(0.1,400)
    es = np.random.uniform(1,20)*vacuum_permittivity
    Na = np.power(10, np.random.uniform(22, 25))
    Nd = np.power(10, np.random.uniform(22, 25))

    V = np.linspace(-6, 4, 20)
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)

    wn = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
    wp = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    w = wn + wp + xi

    expected = (q*d*minor/l)*(np.exp(q*V/(kb*T))-1)*((((s*l)/d)*np.cosh((x-w)/l)+np.sinh((x-w)/l))/
                                                  (((s*l)/d)*np.sinh((x-w)/l)+np.cosh((x-w)/l)))

    result = get_j_dark(x, w, l, s, d, V, minor, T)

    assert result == approx(expected, nan_ok=True)


def test_factor():
    from solcore.analytic_solar_cells.depletion_approximation import factor
    T = np.random.uniform(0.1,400)
    Vbi = np.random.uniform(0.1,5)
    tp = np.power(10, np.random.uniform(-10, -5))
    tn = np.power(10, np.random.uniform(-10, -5))
    dEt = 0

    kT = kb*T
    V = np.linspace(-6, 4, 20)
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)

    m = V >= -1120 * kT / q
    V = V[m]

    b= np.exp(-V*q/(2*kT))*np.cosh((q*dEt/kT)+(1/2)*np.log(tp/tn))

    z1 = np.sqrt(tp/tn)*np.exp(-q*(Vbi - V)/(kT * 2))
    z2 = np.sqrt(tp/tn)*np.exp(q*(Vbi - V)/(kT * 2))

    expected = np.zeros(b.shape)

    # For b values < 1
    inds = b < 1  # use tan-1 formulation without issue

    upper = np.arctan((b[inds]+z2[inds])/np.sqrt(1-b[inds]**2))
    lower = np.arctan((b[inds]+z1[inds])/np.sqrt(1-b[inds]**2))

    expected[inds] = (upper - lower)/np.sqrt(1-b[inds]**2)

    # for b values >=1, log formulation
    inds = (b >= 1) & (b <= 1e7)

    upper = np.log(abs(z2[inds] + b[inds] - np.sqrt(b[inds]**2-1))) - \
            np.log(abs(z2[inds] + b[inds] + np.sqrt(b[inds]**2-1)))

    lower = np.log(abs(z1[inds] + b[inds] - np.sqrt(b[inds]**2-1))) - \
            np.log(abs(z1[inds] + b[inds] + np.sqrt(b[inds]**2-1)))
    expected[inds] = (upper - lower)/(2*np.sqrt(b[inds]**2 - 1))

    # limit as b gets very large: sqrt(b**2-1) = b

    inds = b > 1e7

    upper = np.log(z2[inds]) - np.log(z2[inds] + 2*b[inds])
    lower = np.log(z1[inds]) - np.log(z1[inds] + 2*b[inds])

    expected[inds] = (upper - lower) / (2*b[inds])

    result = factor(V, Vbi, tp, tn, kT, dEt)

    assert result == approx(expected, nan_ok=True)


def test_forward():
    from solcore.analytic_solar_cells.depletion_approximation import factor, forward

    T = np.random.uniform(0.1,400)
    Vbi = np.random.uniform(0.1,5)
    tp = np.power(10, np.random.uniform(-10, -5))
    tn = np.power(10, np.random.uniform(-10, -5))
    dEt = 0

    kT = kb*T
    V = np.linspace(-6, 4, 20)
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)

    ni = np.power(10, np.random.uniform(2, 9))
    es = np.random.uniform(1,20)*vacuum_permittivity
    Na = np.power(10, np.random.uniform(22, 25))
    Nd = np.power(10, np.random.uniform(22, 25))
    xi = np.power(10, np.random.uniform(-8, -5))

    m = V >= -1120 * kT / q
    V = V[m]

    wn = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
    wp = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    w = wn + wp + xi

    f_b = factor(V, Vbi, tp, tn, kT, dEt)
    expected = 2 * q * ni * w / np.sqrt(tn * tp) * \
          np.sinh(q*V / (2 * kT)) / (q * (Vbi - V) / kT) * f_b

    result = forward(ni, V, Vbi, tp, tn, w, kT, dEt)

    assert result == approx(expected, nan_ok=True)


def test_get_J_srh():
    from solcore.analytic_solar_cells.depletion_approximation import forward, get_Jsrh

    T = np.random.uniform(0.1,400)
    Vbi = np.random.uniform(0.1,5)
    tp = np.power(10, np.random.uniform(-10, -5))
    tn = np.power(10, np.random.uniform(-10, -5))
    dEt = 0

    kT = kb*T
    V = np.linspace(-6, 4, 20)
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)

    ni = np.power(10, np.random.uniform(2, 9))
    es = np.random.uniform(1,20)*vacuum_permittivity
    Na = np.power(10, np.random.uniform(22, 25))
    Nd = np.power(10, np.random.uniform(22, 25))
    xi = np.power(10, np.random.uniform(-8, -5))

    wn = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
    wp = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    w = wn + wp + xi

    expected = np.zeros(V.shape)

    inds = V >= -1120 * kT / q
    expected[inds] = forward(ni, V[inds], Vbi, tp, tn, w[inds], kT)

    result = get_Jsrh(ni, V, Vbi, tp, tn, w, kT, dEt)

    assert result == approx(expected, nan_ok=True)


def test_get_J_sc_diffusion_top():
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_diffusion
    from scipy.integrate import solve_bvp
    from solcore.interpolate import interp1d
    from solcore.light_source import LightSource

    D = np.power(10, np.random.uniform(-5, 0)) # Diffusion coefficient
    L = np.power(10, np.random.uniform(-9, -6)) # Diffusion length
    minority = np.power(10, np.random.uniform(-7, -4))# minority carrier density
    s = np.power(10, np.random.uniform(0, 3)) # surface recombination velocity

    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = np.random.uniform(1, 1000)
    xa = xa_nm*1e-9
    xb = np.random.uniform(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e8, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)


    zz = np.linspace(xa, xb, 1002)[:-1]
    gg = gen_prof(zz) * phg

    g_vs_z = np.trapz(gg, wl, axis=1)

    A = lambda x: np.interp(x, zz, g_vs_z) / D + minority / L ** 2

    def fun(x, y):
        out1 = y[1]
        out2 = y[0] / L ** 2 - A(x)
        return np.vstack((out1, out2))

    def bc(ya, yb):
        left = ya[1] - s / D * (ya[0] - minority)
        right = yb[0]
        return np.array([left, right])


    guess = minority * np.ones((2, zz.size))
    guess[1] = np.zeros_like(guess[0])

    solution = solve_bvp(fun, bc, zz, guess)

    expected = solution.y[1][-1]

    result = get_J_sc_diffusion(xa, xb, gen_prof, D, L, minority, s, wl, phg, side='top')

    assert result == approx(expected)


def test_get_J_sc_diffusion_bottom():
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_diffusion
    from scipy.integrate import solve_bvp
    from solcore.interpolate import interp1d
    from solcore.light_source import LightSource

    D = np.power(10, np.random.uniform(-5, 0)) # Diffusion coefficient
    L = np.power(10, np.random.uniform(-9, -6)) # Diffusion length
    minority = np.power(10, np.random.uniform(-7, -4))# minority carrier density
    s = np.power(10, np.random.uniform(0, 3)) # surface recombination velocity

    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = np.random.uniform(1, 1000)
    xa = xa_nm*1e-9
    xb = np.random.uniform(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e8, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)


    zz = np.linspace(xa, xb, 1002)[:-1]
    gg = gen_prof(zz) * phg

    g_vs_z = np.trapz(gg, wl, axis=1)

    A = lambda x: np.interp(x, zz, g_vs_z) / D + minority / L ** 2

    def fun(x, y):
        out1 = y[1]
        out2 = y[0] / L ** 2 - A(x)
        return np.vstack((out1, out2))

    def bc(ya, yb):
        left = ya[0]
        right = yb[1] - s / D * (yb[0] - minority)
        return np.array([left, right])


    guess = minority * np.ones((2, zz.size))
    guess[1] = np.zeros_like(guess[0])

    solution = solve_bvp(fun, bc, zz, guess)

    expected = solution.y[1][0]

    result = get_J_sc_diffusion(xa, xb, gen_prof, D, L, minority, s, wl, phg, side='bottom')

    assert result == approx(expected)


def test_get_J_sc_SCR():
    from solcore.light_source import LightSource
    from solcore.interpolate import interp1d
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_SCR

    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = np.random.uniform(1, 1000)
    xa = xa_nm*1e-9
    xb = np.random.uniform(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e5, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)

    zz = np.linspace(xa, xb, 1002)[:-1]
    gg = gen_prof(zz) * phg
    expected = np.trapz(np.trapz(gg, wl, axis=1), zz)

    result = get_J_sc_SCR(xa, xb, gen_prof, wl, phg)

    assert  expected == approx(result)


def test_get_J_sc_SCR_vs_WL():
    from solcore.light_source import LightSource
    from solcore.interpolate import interp1d
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_SCR_vs_WL

    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = np.random.uniform(1, 1000)
    xa = xa_nm*1e-9
    xb = np.random.uniform(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e5, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)

    zz = np.linspace(xa, xb, 1002)[:-1]
    gg = gen_prof(zz) * phg
    expected = np.trapz(gg, zz, axis=0)

    result = get_J_sc_SCR_vs_WL(xa, xb, gen_prof, wl, phg)

    assert  expected == approx(result)


def test_get_J_sc_diffusion_vs_WL_top():
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_diffusion_vs_WL
    from scipy.integrate import solve_bvp
    from solcore.interpolate import interp1d
    from solcore.light_source import LightSource

    D = np.power(10, np.random.uniform(-5, 0)) # Diffusion coefficient
    L = np.power(10, np.random.uniform(-9, -6)) # Diffusion length
    minority = np.power(10, np.random.uniform(-7, -4))# minority carrier density
    s = np.power(10, np.random.uniform(0, 3)) # surface recombination velocity

    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = np.random.uniform(1, 1000)
    xa = xa_nm*1e-9
    xb = np.random.uniform(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e8, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)


    zz =  np.linspace(xa, xb, 1002)[:-1]
    gg = gen_prof(zz) * phg

    expected = np.zeros_like(wl)

    for i in range(len(wl)):

        if np.all(gg[:,i] == 0): # no reason to solve anything if no generation at this wavelength
            expected[i] = 0

        A = lambda x: np.interp(x, zz, gg[:, i]) / D + minority / L ** 2

        def fun(x, y):
            out1 = y[1]
            out2 = y[0] / L ** 2 - A(x)
            return np.vstack((out1, out2))

        def bc(ya, yb):
            left = ya[1] - s / D * (ya[0] - minority)
            right = yb[0]
            return np.array([left, right])

        guess = minority * np.ones((2, zz.size))
        guess[1] = np.zeros_like(guess[0])
        solution = solve_bvp(fun, bc, zz, guess)

        expected[i] = solution.y[1][-1]

    result = get_J_sc_diffusion_vs_WL(xa, xb, gen_prof, D, L, minority, s, wl, phg, side='top')


    assert result == approx(expected)


def test_get_J_sc_diffusion_vs_WL_bottom():
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_diffusion_vs_WL
    from scipy.integrate import solve_bvp
    from solcore.interpolate import interp1d
    from solcore.light_source import LightSource

    D = np.power(10, np.random.uniform(-5, 0)) # Diffusion coefficient
    L = np.power(10, np.random.uniform(-9, -6)) # Diffusion length
    minority = np.power(10, np.random.uniform(-7, -4))# minority carrier density
    s = np.power(10, np.random.uniform(0, 3)) # surface recombination velocity

    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50) * 1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = np.random.uniform(1, 1000)
    xa = xa_nm * 1e-9
    xb = np.random.uniform(xa_nm + 1, 1100) * 1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e8, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None, :])

    output = alphas[:, None] * expn
    output = output.T
    gen_prof = interp1d(dist, output, axis=0)

    zz = np.linspace(xa, xb, 1002)[:-1]
    gg = gen_prof(zz) * phg

    expected = np.zeros_like(wl)

    for i in range(len(wl)):
        if np.all(gg[:, i] == 0):  # no reason to solve anything if no generation at this wavelength
            expected[i] = 0

        A = lambda x: np.interp(x, zz, gg[:, i]) / D + minority / L ** 2

        def fun(x, y):
            out1 = y[1]
            out2 = y[0] / L ** 2 - A(x)
            return np.vstack((out1, out2))

        def bc(ya, yb):
            left = ya[0]
            right = yb[1] - s / D * (yb[0] - minority)
            return np.array([left, right])

        guess = minority * np.ones((2, zz.size))
        guess[1] = np.zeros_like(guess[0])
        solution = solve_bvp(fun, bc, zz, guess)

        expected[i] = solution.y[1][0]

    result = get_J_sc_diffusion_vs_WL(xa, xb, gen_prof, D, L, minority, s, wl, phg, side='bottom')

    assert result == approx(expected)


def test_identify_layers_exceptions():
    from solcore.analytic_solar_cells.depletion_approximation import identify_layers
    from solcore import material
    from solcore.structure import Layer, Junction

    Na = np.power(10, np.random.uniform(22, 25))
    Nd = np.power(10, np.random.uniform(22, 25))

    Lp = np.power(10, np.random.uniform(-9, -6))# Diffusion length
    Ln = np.power(10, np.random.uniform(-9, -6)) # Diffusion length

    GaAs_n = material("GaAs")(Nd = Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na = Na, electron_diffusion_length=Lp)
    GaAs_i = material("GaAs")()
    Ge_n = material("Ge")(Nd = Nd, hole_diffusion_length=Ln)

    n_width = np.random.uniform(500, 1000)*1e-9
    p_width = np.random.uniform(3000, 5000)*1e-9
    i_width = np.random.uniform(300, 500) * 1e-9

    test_junc  = Junction([Layer(n_width, GaAs_n,role="emitter"),
                           Layer(p_width, GaAs_p, role="neither")])

    with raises(RuntimeError):
        identify_layers(test_junc)

    test_junc =  Junction([Layer(n_width, GaAs_n,role="emitter"),
                           Layer(i_width, GaAs_i, role="intrinsic"),
                           Layer(p_width, GaAs_p, role="nothing")])

    with raises(RuntimeError):
        identify_layers(test_junc)

    test_junc  = Junction([Layer(n_width, Ge_n,role="emitter"),
                           Layer(p_width, GaAs_p, role="base")])

    with raises(AssertionError):
        identify_layers(test_junc)


def test_process_junction_np():
    from solcore.analytic_solar_cells.depletion_approximation import identify_layers, identify_parameters
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Na = np.power(10, np.random.uniform(23, 26))
    Nd = np.power(10, np.random.uniform(22, 25))

    options = State()
    options.T = np.random.uniform(250, 350)

    Lp = np.power(10, np.random.uniform(-9, -6))# Diffusion length
    Ln = np.power(10, np.random.uniform(-9, -6)) # Diffusion length

    GaAs_window = material("GaAs")()
    GaAs_n = material("GaAs")(Nd = Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na = Na, electron_diffusion_length=Lp)

    n_width = np.random.uniform(500, 1000)*1e-9
    p_width = np.random.uniform(3000, 5000)*1e-9
    window_width = np.random.uniform(25, 200)*1e-9

    test_junc  = Junction([Layer(window_width, GaAs_window, role="window"),
                           Layer(n_width, GaAs_n,role="emitter"),
                           Layer(p_width, GaAs_p, role="base")])

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(test_junc)
    xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd_c, Na_c, ni, es = identify_parameters(test_junc, options.T, pRegion, nRegion, iRegion)

    ni_expect = GaAs_n.ni

    assert [id_top, id_bottom] == approx([1, 2])
    assert pn_or_np == 'np'
    assert [xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd_c, Na_c, ni, es] == approx([n_width, p_width, 0, 0, 0, Lp, Ln, GaAs_n.hole_mobility * kb * options.T / q,
                                    GaAs_p.electron_mobility * kb * options.T / q, Nd, Na, ni_expect, GaAs_n.permittivity])




def test_process_junction_pn():
    from solcore.analytic_solar_cells.depletion_approximation import identify_layers, identify_parameters
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Na = np.power(10, np.random.uniform(22, 25))
    Nd = np.power(10, np.random.uniform(23, 26))

    options = State()
    options.T = np.random.uniform(250, 350)

    Lp = np.power(10, np.random.uniform(-9, -6))  # Diffusion length
    Ln = np.power(10, np.random.uniform(-9, -6))  # Diffusion length

    GaAs_n = material("GaAs")(Nd=Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na=Na, electron_diffusion_length=Lp)

    p_width = np.random.uniform(500, 1000)*1e-9
    n_width = np.random.uniform(3000, 5000)*1e-9

    test_junc = Junction([Layer(p_width, GaAs_p, role="emitter"), Layer(n_width, GaAs_n, role="base")])

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(test_junc)
    xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd_c, Na_c, ni, es = identify_parameters(test_junc, options.T, pRegion, nRegion, iRegion)

    ni_expect = GaAs_n.ni

    assert [id_top, id_bottom] == approx([0, 1])
    assert pn_or_np == 'pn'
    assert [xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd_c, Na_c, ni, es] == approx([n_width, p_width, 0, 0, 0, Lp, Ln, GaAs_n.hole_mobility * kb * options.T / q,
                                    GaAs_p.electron_mobility * kb * options.T / q, Nd, Na, ni_expect, GaAs_n.permittivity])


def test_process_junction_nip():
    from solcore.analytic_solar_cells.depletion_approximation import identify_layers, identify_parameters
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Na = np.power(10, np.random.uniform(23, 26))
    Nd = np.power(10, np.random.uniform(22, 25))

    options = State()
    options.T = np.random.uniform(250, 350)

    Lp = np.power(10, np.random.uniform(-9, -6)) # Diffusion length
    Ln = np.power(10, np.random.uniform(-9, -6)) # Diffusion length

    GaAs_n = material("GaAs")(Nd = Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na = Na, electron_diffusion_length=Lp)
    GaAs_i = material("GaAs")()

    n_width = np.random.uniform(500, 1000)*1e-9
    p_width = np.random.uniform(3000, 5000)*1e-9
    i_width = np.random.uniform(100, 300)*1e-9

    test_junc  = Junction([Layer(n_width, GaAs_n,role="emitter"),
                           Layer(i_width, GaAs_i, role="intrinsic"),
                           Layer(p_width, GaAs_p, role="base")])

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(test_junc)
    xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd_c, Na_c, ni, es = identify_parameters(test_junc, options.T, pRegion, nRegion, iRegion)

    ni_expect = GaAs_n.ni

    assert [id_top, id_bottom] == approx([0, 2])
    assert pn_or_np == 'np'
    assert [xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd_c, Na_c, ni, es] == approx([n_width, p_width, i_width, 0, 0, Lp, Ln, GaAs_n.hole_mobility * kb * options.T / q,
                                    GaAs_p.electron_mobility * kb * options.T / q, Nd, Na, ni_expect, GaAs_n.permittivity])



def test_process_junction_pin():
    from solcore.analytic_solar_cells.depletion_approximation import identify_layers, identify_parameters
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Na = np.power(10, np.random.uniform(22, 25))
    Nd = np.power(10, np.random.uniform(23, 26))

    options = State()
    options.T = np.random.uniform(250, 350)

    Lp = np.power(10, np.random.uniform(-9, -6))  # Diffusion length
    Ln = np.power(10, np.random.uniform(-9, -6)) # Diffusion length

    GaAs_n = material("GaAs")(Nd=Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na=Na, electron_diffusion_length=Lp)
    GaAs_i = material("GaAs")()

    p_width = np.random.uniform(500, 1000)
    n_width = np.random.uniform(3000, 5000)
    i_width = np.random.uniform(100, 300)*1e-9

    test_junc = Junction([Layer(p_width, GaAs_p, role="emitter"),
                          Layer(i_width, GaAs_i, role="intrinsic"),
                          Layer(n_width, GaAs_n, role="base")])

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(test_junc)
    xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd_c, Na_c, ni, es = identify_parameters(test_junc, options.T, pRegion, nRegion, iRegion)

    ni_expect = GaAs_n.ni

    assert [id_top, id_bottom] == approx([0, 2])
    assert pn_or_np == 'pn'
    assert [xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd_c, Na_c, ni, es] == approx([n_width, p_width, i_width, 0, 0, Lp, Ln, GaAs_n.hole_mobility * kb * options.T / q,
                                    GaAs_p.electron_mobility * kb * options.T / q, Nd, Na, ni_expect, GaAs_n.permittivity])


def test_process_junction_set_in_junction():
    from solcore.analytic_solar_cells.depletion_approximation import identify_layers, identify_parameters
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    options = State()
    options.T = np.random.uniform(250, 350)

    Lp = np.power(10, np.random.uniform(-9, -6)) # Diffusion length
    Ln = np.power(10, np.random.uniform(-9, -6)) # Diffusion length

    sn = np.power(10, np.random.uniform(0, 3))
    sp = np.power(10, np.random.uniform(0, 3))

    se = np.random.uniform(1, 20)*vacuum_permittivity

    GaAs_n = material("GaAs")()
    GaAs_p = material("GaAs")()
    GaAs_i = material("GaAs")()

    p_width = np.random.uniform(500, 1000)
    n_width = np.random.uniform(3000, 5000)
    i_width = np.random.uniform(100, 300)*1e-9

    mun = np.power(10, np.random.uniform(-5, 0))
    mup = np.power(10, np.random.uniform(-5, 0))

    Vbi = np.random.uniform(0, 3)

    test_junc = Junction([Layer(p_width, GaAs_p, role="emitter"),
                          Layer(i_width, GaAs_i, role="intrinsic"),
                          Layer(n_width, GaAs_n, role="base")],
                         sn = sn, sp = sp, permittivity = se,
                         ln = Ln, lp= Lp,
                         mup = mup, mun = mun, Vbi = Vbi)

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(test_junc)
    xn, xp, xi, sn_c, sp_c, ln, lp, dn, dp, Nd_c, Na_c, ni, es = identify_parameters(test_junc, options.T, pRegion, nRegion, iRegion)

    ni_expect = GaAs_n.ni

    assert [id_top, id_bottom] == approx([0, 2])
    assert pn_or_np == 'pn'
    assert [xn, xp, xi, sn_c, sp_c, ln, lp, dn, dp, Nd_c, Na_c, ni, es] == approx([n_width, p_width, i_width, sn, sp, Ln, Lp, mun * kb * options.T / q,
                                    mup * kb * options.T / q, 1, 1, ni_expect, se])


def test_get_depletion_widths():
    from solcore.analytic_solar_cells.depletion_approximation import get_depletion_widths
    from solcore.structure import Junction

    xi  = np.power(10, np.random.uniform(-10, -6))

    Vbi = np.random.uniform(0, 3)
    es = np.random.uniform(1, 20)*vacuum_permittivity
    Na = np.power(10, np.random.uniform(22, 25))
    Nd = np.power(10, np.random.uniform(22, 25))

    V = np.linspace(-6, 4, 20)
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)

    test_junc = Junction()

    wn_e = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
    wp_e = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    wn_r, wp_r = get_depletion_widths(test_junc, es, Vbi, V, Na, Nd, xi)

    assert wn_r == approx(wn_e)
    assert wp_r == approx(wp_e)


def test_get_depletion_widths_onesided():
    from solcore.analytic_solar_cells.depletion_approximation import get_depletion_widths
    from solcore.structure import Junction

    xi  = np.power(10, np.random.uniform(-10, -6))

    Vbi = np.random.uniform(0, 3)
    es = np.random.uniform(1, 20)*vacuum_permittivity
    Na = np.power(10, np.random.uniform(22, 25))
    Nd = np.power(10, np.random.uniform(22, 25))

    V = np.linspace(-6, 4, 20)
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)

    test_junc = Junction(depletion_approximation="one-sided abrupt")

    wn_e = np.sqrt(2 * es * (Vbi - V) / (q * Nd))
    wp_e = np.sqrt(2 * es * (Vbi - V) / (q * Na))

    wn_r, wp_r = get_depletion_widths(test_junc, es, Vbi, V, Na, Nd, xi)

    assert wn_r == approx(wn_e)
    assert wp_r == approx(wp_e)


def test_get_depletion_widths_set_in_junction():
    from solcore.analytic_solar_cells.depletion_approximation import get_depletion_widths
    from solcore.structure import Junction

    wn = np.random.uniform(1,100)
    wp = np.random.uniform(1,100)
    test_junc = Junction(wn=wn, wp=wp)

    wn_r, wp_r = get_depletion_widths(test_junc, 0, 0, 0, 0, 0, 0)

    assert wn_r == approx(wn)
    assert wp_r == approx(wp)


def test_dark_iv_depletion_pn(pn_junction):
    from solcore.analytic_solar_cells.depletion_approximation import iv_depletion, get_depletion_widths, get_j_dark, get_Jsrh, identify_layers, identify_parameters
    from scipy.interpolate import interp1d

    test_junc, options = pn_junction
    options.light_iv = False
    T = options.T

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(test_junc[0])
    xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd, Na, ni, es = identify_parameters(test_junc[0], T, pRegion, nRegion, iRegion)

    niSquared = ni**2

    kbT = kb * T

    Vbi = (kbT / q) * np.log(Nd * Na / niSquared)

    test_junc[0].voltage = options.internal_voltages

    V = np.where(test_junc[0].voltage < Vbi - 0.001, test_junc[0].voltage, Vbi - 0.001)

    wn, wp = get_depletion_widths(test_junc[0], es, Vbi, V, Na, Nd, xi)

    w = wn + wp + xi

    l_top, l_bottom = ln, lp
    x_top, x_bottom = xp, xn
    w_top, w_bottom = wp, wn
    s_top, s_bottom = sp, sn
    d_top, d_bottom = dp, dn
    min_top, min_bot = niSquared / Na, niSquared / Nd

    JtopDark = get_j_dark(x_top, w_top, l_top, s_top, d_top, V, min_top, T)
    JbotDark = get_j_dark(x_bottom, w_bottom, l_bottom, s_bottom, d_bottom, V, min_bot, T)

    JnDark, JpDark = JbotDark, JtopDark

    lifetime_n = ln ** 2 / dn
    lifetime_p = lp ** 2 / dp  # Jenny p163

    Jrec = get_Jsrh(ni, V, Vbi, lifetime_p, lifetime_n, w, kbT)

    J_sc_top = 0
    J_sc_bot = 0
    J_sc_scr = 0

    current = Jrec + JnDark + JpDark + V / 1e14- J_sc_top - J_sc_bot - J_sc_scr
    iv = interp1d(test_junc[0].voltage, current, kind='linear', bounds_error=False, assume_sorted=True,
                           fill_value=(current[0], current[-1]), copy=True)

    iv_depletion(test_junc[0], options)

    assert test_junc[0].iv(options.internal_voltages) == approx(iv(options.internal_voltages), nan_ok=True)


def test_dark_iv_depletion_np(np_junction):
    from solcore.analytic_solar_cells.depletion_approximation import iv_depletion, get_depletion_widths, get_j_dark, get_Jsrh, identify_layers, identify_parameters
    from scipy.interpolate import interp1d

    test_junc, options = np_junction
    options.light_iv = False
    T = options.T

    test_junc[0].voltage = options.internal_voltages

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(test_junc[0])
    xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd, Na, ni, es = identify_parameters(test_junc[0], T, pRegion, nRegion, iRegion)

    niSquared = ni**2

    kbT = kb * T

    Vbi = (kbT / q) * np.log(Nd * Na / niSquared)

    V = np.where(test_junc[0].voltage < Vbi - 0.001, test_junc[0].voltage, Vbi - 0.001)

    wn, wp = get_depletion_widths(test_junc[0], es, Vbi, V, Na, Nd, xi)

    w = wn + wp + xi

    l_bottom, l_top = ln, lp
    x_bottom, x_top = xp, xn
    w_bottom, w_top = wp, wn
    s_bottom, s_top = sp, sn
    d_bottom, d_top = dp, dn
    min_bot, min_top = niSquared / Na, niSquared / Nd

    JtopDark = get_j_dark(x_top, w_top, l_top, s_top, d_top, V, min_top, T)
    JbotDark = get_j_dark(x_bottom, w_bottom, l_bottom, s_bottom, d_bottom, V, min_bot, T)

    JpDark, JnDark = JbotDark, JtopDark

    lifetime_n = ln ** 2 / dn
    lifetime_p = lp ** 2 / dp  # Jenny p163

    Jrec = get_Jsrh(ni, V, Vbi, lifetime_p, lifetime_n, w, kbT)

    J_sc_top = 0
    J_sc_bot = 0
    J_sc_scr = 0

    current = Jrec + JnDark + JpDark + V / 1e14- J_sc_top - J_sc_bot - J_sc_scr
    iv = interp1d(test_junc[0].voltage, current, kind='linear', bounds_error=False, assume_sorted=True,
                           fill_value=(current[0], current[-1]), copy=True)

    iv_depletion(test_junc[0], options)

    assert test_junc[0].iv(options.internal_voltages) == approx(iv(options.internal_voltages), nan_ok=True)


def test_qe_depletion_np(np_junction):

    from solcore.analytic_solar_cells import qe_depletion

    test_junc, options = np_junction

    wl = options.wavelength

    qe_depletion(test_junc[0], options)

    assert np.all(test_junc[0].eqe(wl) < 1)
    assert np.all(test_junc[0].eqe_emitter(wl) < 1)
    assert np.all(test_junc[0].eqe_base(wl) < 1)
    assert np.all(test_junc[0].eqe_scr(wl) < 1)
    assert np.all(test_junc[0].eqe(wl)[test_junc[0].eqe(wl) > 1e-3]
                  <= test_junc.absorbed[test_junc[0].eqe(wl) > 1e-3])
    assert np.all(test_junc[0].eqe_emitter(wl) + test_junc[0].eqe_base(wl) +
                  test_junc[0].eqe_scr(wl) == approx(test_junc[0].eqe(wl)))
    assert np.all(test_junc[0].iqe(wl) >= test_junc[0].eqe(wl))



def test_qe_depletion_pn(pn_junction):
    from solcore.analytic_solar_cells import qe_depletion

    test_junc, options = pn_junction

    wl = options.wavelength

    qe_depletion(test_junc[0], options)

    assert np.all(test_junc[0].eqe(wl) < 1)
    assert np.all(test_junc[0].eqe_emitter(wl) < 1)
    assert np.all(test_junc[0].eqe_base(wl) < 1)
    assert np.all(test_junc[0].eqe_scr(wl) < 1)
    assert np.all(test_junc[0].eqe(wl)[test_junc[0].eqe(wl) > 1e-3]
                  <= test_junc.absorbed[test_junc[0].eqe(wl) > 1e-3])
    assert np.all(test_junc[0].eqe_emitter(wl) + test_junc[0].eqe_base(wl) +
                  test_junc[0].eqe_scr(wl) == approx(test_junc[0].eqe(wl)))
    assert np.all(test_junc[0].iqe(wl) >= test_junc[0].eqe(wl))


def test_iv_depletion_np(np_junction):

    from solcore.analytic_solar_cells import iv_depletion

    test_junc, options = np_junction
    options.light_iv = True
    V = options.internal_voltages
    wl = options.wavelength

    iv_depletion(test_junc[0], options)

    wl_sp, ph = options.light_source.spectrum(output_units='photon_flux_per_m', x=wl)
    Jph = q*np.trapz(test_junc.absorbed*ph, wl)

    approx_Voc = V[np.argmin(abs(test_junc[0].iv(V)))]

    quadrant = (V > 0) * (V < approx_Voc)

    power = abs(test_junc[0].iv(V[quadrant])*V[quadrant])[:-1]


    assert abs(test_junc[0].iv(0)) <= Jph
    assert approx_Voc < test_junc[0][1].material.band_gap/q
    assert np.all(power < options.light_source.power_density)


def test_iv_depletion_pn(pn_junction):
    from solcore.analytic_solar_cells import iv_depletion

    test_junc, options = pn_junction
    options.light_iv = True
    V = options.internal_voltages
    wl = options.wavelength

    iv_depletion(test_junc[0], options)

    wl_sp, ph = options.light_source.spectrum(output_units='photon_flux_per_m', x=wl)
    Jph = q*np.trapz(test_junc.absorbed*ph, wl)

    approx_Voc = V[np.argmin(abs(test_junc[0].iv(V)))]

    quadrant = (V > 0) * (V < approx_Voc)

    power = abs(test_junc[0].iv(V[quadrant])*V[quadrant])[:-1]

    assert abs(test_junc[0].iv(0)) <= Jph
    assert approx_Voc < test_junc[0][1].material.band_gap/q
    assert np.all(power < options.light_source.power_density)

