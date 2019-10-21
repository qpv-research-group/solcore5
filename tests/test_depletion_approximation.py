from pytest import approx, mark, raises
from random import randrange
import numpy as np

from solcore.constants import kb, q, vacuum_permittivity


def test_get_j_dark():
    from solcore.analytic_solar_cells.depletion_approximation import get_j_dark

    xnm = randrange(2, 1000)
    x = xnm*1e-9
    xi  = randrange(0, 1000)*1e-9

    l = randrange(5, 10000)*1e-9
    s = randrange(1, 1000)
    d = randrange(1, 1e5)*1e-5
    Vbi = randrange(1, 50)*1e-1
    minor = randrange(1, 1e6)*1e-8
    T = randrange(1, 300)
    es = randrange(1, 20)*vacuum_permittivity
    Na = randrange(1, 1e5)*1e19
    Nd = randrange(1, 1e5)*1e19

    V = np.linspace(-6, 4, 20)
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)

    wn = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
    wp = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    w = wn + wp + xi

    expected = (q*d*minor/l)*(np.exp(q*V/(kb*T))-1)*((((s*l)/d)*np.cosh((x-w)/l)+np.sinh((x-w)/l))/
                                                  (((s*l)/d)*np.sinh((x-w)/l)+np.cosh((x-w)/l)))

    result = get_j_dark(x, w, l, s, d, V, minor, T)

    assert result == approx(expected)

def test_factor():
    from solcore.analytic_solar_cells.depletion_approximation import factor
    T = randrange(1, 300)
    kT = kb*T
    V = np.linspace(-6, 4, 20)
    Vbi = randrange(1, 50)*1e-1
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)
    tp = randrange(1, 1e5)*1e-10
    tn = randrange(1, 1e5)*1e-10
    dEt = 0

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

    T = randrange(1, 300)
    kT = kb*T
    V = np.linspace(-6, 4, 20)
    Vbi = randrange(1, 50)*1e-1
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)
    tp = randrange(1, 1e5)*1e-10
    tn = randrange(1, 1e5)*1e-10
    dEt = 0
    ni = randrange(1, 1e7)*1e2
    es = randrange(1, 20)*vacuum_permittivity
    Na = randrange(1, 1e5)*1e19
    Nd = randrange(1, 1e5)*1e19
    xi  = randrange(0, 1000)*1e-9

    m = V >= -1120 * kT / q
    V = V[m]

    wn = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
    wp = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    w = wn + wp + xi
    #w = w[m]

    f_b = factor(V, Vbi, tp, tn, kT, dEt)
    expected = 2 * q * ni * w / np.sqrt(tn * tp) * \
          np.sinh(q*V / (2 * kT)) / (q * (Vbi - V) / kT) * f_b

    result = forward(ni, V, Vbi, tp, tn, w, kT, dEt)

    assert result == approx(expected, nan_ok=True)


# this test failed one, not sure on conditions (rng)
def test_get_J_srh():
    from solcore.analytic_solar_cells.depletion_approximation import forward, get_Jsrh
    T = randrange(1, 300)
    kT = kb*T
    V = np.linspace(-6, 4, 20)
    Vbi = randrange(1, 50)*1e-1
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)
    tp = randrange(1, 1e5)*1e-10
    tn = randrange(1, 1e5)*1e-10
    dEt = 0
    ni = randrange(1, 1e7)*1e2
    es = randrange(1, 20)*vacuum_permittivity
    Na = randrange(1, 1e5)*1e19
    Nd = randrange(1, 1e5)*1e19
    xi  = randrange(0, 1000)*1e-9

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

    D = randrange(1, 1e5)*1e-5 # Diffusion coefficient
    L = randrange(5, 10000)*1e-9 # Diffusion length
    minority = randrange(1, 70)# minority carrier density
    s = randrange(1, 1000) # surface recombination velocity


    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = randrange(1, 1000)
    xa = xa_nm*1e-9
    xb = randrange(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e8, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)


    zz = np.linspace(xa, xb, 1001)
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

    print(result)

    assert result == approx(expected)


def test_get_J_sc_diffusion_bottom():
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_diffusion
    from scipy.integrate import solve_bvp
    from solcore.interpolate import interp1d
    from solcore.light_source import LightSource

    D = randrange(1, 1e5)*1e-5 # Diffusion coefficient
    L = randrange(5, 10000)*1e-9 # Diffusion length
    minority = randrange(1, 70)# minority carrier density
    s = randrange(1, 1000) # surface recombination velocity


    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = randrange(1, 1000)
    xa = xa_nm*1e-9
    xb = randrange(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e8, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)


    zz = np.linspace(xa, xb, 1001)
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

    xa_nm = randrange(1, 1000)
    xa = xa_nm*1e-9
    xb = randrange(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e5, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)

    #plt.figure()
    #plt.imshow(output)
    #plt.show()

    zz = np.linspace(xa, xb, 1001)
    gg = gen_prof(zz) * phg
    expected = np.trapz(np.trapz(gg, wl, axis=1), zz)

    result = get_J_sc_SCR(xa, xb, gen_prof, wl, phg)
    # think the units might be wrong (factor of 1e9?) but it doesn't really matter

    assert  expected == approx(result)

# get_Jsrh: what is the 1120?

## QE

# no need to pass wl to this function
def test_get_J_sc_SCR_vs_WL():
    from solcore.light_source import LightSource
    from solcore.interpolate import interp1d
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_SCR_vs_WL

    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = randrange(1, 1000)
    xa = xa_nm*1e-9
    xb = randrange(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e5, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)

    #plt.figure()
    #plt.imshow(output)
    #plt.show()

    zz = np.linspace(xa, xb, 1001)
    gg = gen_prof(zz) * phg
    expected = np.trapz(gg, zz, axis=0)

    result = get_J_sc_SCR_vs_WL(xa, xb, gen_prof, wl, phg)
    # think the units might be wrong (factor of 1e9?) but it doesn't really matter

    assert  expected == approx(result)


def test_get_J_sc_diffusion_vs_WL_top():
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_diffusion_vs_WL
    from scipy.integrate import solve_bvp
    from solcore.interpolate import interp1d
    from solcore.light_source import LightSource

    D = randrange(1, 1e5)*1e-5 # Diffusion coefficient
    L = randrange(5, 10000)*1e-9 # Diffusion length
    minority = randrange(1, 70)# minority carrier density
    s = randrange(1, 1000) # surface recombination velocity


    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50)*1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = randrange(1, 1000)
    xa = xa_nm*1e-9
    xb = randrange(xa_nm+1, 1100)*1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e8, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = alphas[:, None]*expn
    output = output.T
    gen_prof = interp1d(dist, output, axis = 0)


    zz = np.linspace(xa, xb, 1001)
    gg = gen_prof(zz) * phg

    expected = np.zeros_like(wl)

    for i in range(len(wl)):
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
        #print(zz)
        solution = solve_bvp(fun, bc, zz, guess)

        expected[i] = solution.y[1][-1]

    result = get_J_sc_diffusion_vs_WL(xa, xb, gen_prof, D, L, minority, s, wl, phg, side='top')


    assert result == approx(expected)


# very similar to just get_J_sc_diffusion.

def test_get_J_sc_diffusion_vs_WL_bottom():
    from solcore.analytic_solar_cells.depletion_approximation import get_J_sc_diffusion_vs_WL
    from scipy.integrate import solve_bvp
    from solcore.interpolate import interp1d
    from solcore.light_source import LightSource

    D = randrange(1, 1e5) * 1e-5  # Diffusion coefficient
    L = randrange(5, 10000) * 1e-9  # Diffusion length
    minority = randrange(1, 70)  # minority carrier density
    s = randrange(1, 1000)  # surface recombination velocity

    light_source = LightSource(source_type="standard", version="AM1.5g")
    wl = np.linspace(300, 1800, 50) * 1e-9
    wl_ls, phg = light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    xa_nm = randrange(1, 1000)
    xa = xa_nm * 1e-9
    xb = randrange(xa_nm + 1, 1100) * 1e-9

    ## make a simple Beer-Lambert profile
    dist = np.linspace(0, xb, 1000)
    alphas = np.linspace(1e8, 10, len(wl))

    expn = np.exp(- alphas[:, None] * dist[None, :])

    output = alphas[:, None] * expn
    output = output.T
    gen_prof = interp1d(dist, output, axis=0)

    zz = np.linspace(xa, xb, 1001)
    gg = gen_prof(zz) * phg

    expected = np.zeros_like(wl)

    for i in range(len(wl)):
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
        # print(zz)
        solution = solve_bvp(fun, bc, zz, guess)

        expected[i] = solution.y[1][0]

    result = get_J_sc_diffusion_vs_WL(xa, xb, gen_prof, D, L, minority, s, wl, phg, side='bottom')

    assert result == approx(expected)


def test_process_junction_exceptions():
    from solcore.analytic_solar_cells.depletion_approximation import process_junction
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Nd = randrange(1, 10)*1e18
    Na = randrange(1, 9)*1e17

    options = State()
    options.T = randrange(1, 350)

    Lp = randrange(5, 10000)*1e-9 # Diffusion length
    Ln = randrange(5, 10000)*1e-9 # Diffusion length

    GaAs_n = material("GaAs")(Nd = Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na = Na, electron_diffusion_length=Lp)
    GaAs_i = material("GaAs")()

    n_width = randrange(500, 1000)*1e-9
    p_width = randrange(3000, 5000)*1e-9
    i_width = randrange(300, 500) * 1e-9

    test_junc  = Junction([Layer(n_width, GaAs_n,role="emitter"),
                           Layer(p_width, GaAs_p, role="neither")])

    with raises(RuntimeError):
        results = process_junction(test_junc, options)

    test_junc =  Junction([Layer(n_width, GaAs_n,role="emitter"),
                           Layer(i_width, GaAs_i, role="intrinsic"),
                           Layer(p_width, GaAs_p, role="nothing")])

    with raises(RuntimeError):
        results = process_junction(test_junc, options)


def test_process_junction_np():
    from solcore.analytic_solar_cells.depletion_approximation import process_junction
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Nd = randrange(1, 10)*1e18
    Na = randrange(1, 9)*1e17

    options = State()
    options.T = randrange(1, 350)

    Lp = randrange(5, 10000)*1e-9 # Diffusion length
    Ln = randrange(5, 10000)*1e-9 # Diffusion length

    GaAs_n = material("GaAs")(Nd = Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na = Na, electron_diffusion_length=Lp)

    n_width = randrange(500, 1000)*1e-9
    p_width = randrange(3000, 5000)*1e-9

    test_junc  = Junction([Layer(n_width, GaAs_n,role="emitter"), Layer(p_width, GaAs_p, role="base")])

    results = process_junction(test_junc, options)

    ni_expect = GaAs_n.ni
    niSquared_expect = ni_expect**2

    Vbi_expect = (kb*options.T / q) * np.log(Nd * Na / niSquared_expect)

    assert results[0:17] == approx((Na, Nd, ni_expect, niSquared_expect, 0,
                             Lp, Ln, n_width, p_width,
                             0, 0, GaAs_n.hole_mobility * kb*options.T / q,
                                    GaAs_p.electron_mobility * kb * options.T / q,
                              GaAs_n.permittivity, 0, 1, Vbi_expect))

    assert results[17] == 'np'


def test_process_junction_pn():
    from solcore.analytic_solar_cells.depletion_approximation import process_junction
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Nd = randrange(1, 9) * 1e17
    Na = randrange(1, 10) * 1e18

    options = State()
    options.T = randrange(1, 350)

    Lp = randrange(5, 10000) * 1e-9  # Diffusion length
    Ln = randrange(5, 10000) * 1e-9  # Diffusion length

    GaAs_n = material("GaAs")(Nd=Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na=Na, electron_diffusion_length=Lp)

    p_width = randrange(500, 1000)*1e-9
    n_width = randrange(3000, 5000)*1e-9

    test_junc = Junction([Layer(p_width, GaAs_p, role="emitter"), Layer(n_width, GaAs_n, role="base")])

    results = process_junction(test_junc, options)

    ni_expect = GaAs_n.ni
    niSquared_expect = ni_expect ** 2

    Vbi_expect = (kb * options.T / q) * np.log(Nd * Na / niSquared_expect)

    assert results[0:17] == approx((Na, Nd, ni_expect, niSquared_expect, 0,
                                    Lp, Ln, n_width, p_width,
                                    0, 0, GaAs_n.hole_mobility * kb * options.T / q,
                                    GaAs_p.electron_mobility * kb * options.T / q,
                                    GaAs_n.permittivity, 0, 1, Vbi_expect))

    assert results[17] == 'pn'


def test_process_junction_nip():
    from solcore.analytic_solar_cells.depletion_approximation import process_junction
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Nd = randrange(1, 10)*1e18
    Na = randrange(1, 9)*1e17

    options = State()
    options.T = randrange(1, 350)

    Lp = randrange(5, 10000)*1e-9 # Diffusion length
    Ln = randrange(5, 10000)*1e-9 # Diffusion length

    GaAs_n = material("GaAs")(Nd = Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na = Na, electron_diffusion_length=Lp)
    GaAs_i = material("GaAs")()

    n_width = randrange(500, 1000)*1e-9
    p_width = randrange(3000, 5000)*1e-9
    i_width = randrange(100, 300)*1e-9

    test_junc  = Junction([Layer(n_width, GaAs_n,role="emitter"),
                           Layer(i_width, GaAs_i, role="intrinsic"),
                           Layer(p_width, GaAs_p, role="base")])

    results = process_junction(test_junc, options)

    ni_expect = GaAs_n.ni
    niSquared_expect = ni_expect**2

    Vbi_expect = (kb*options.T / q) * np.log(Nd * Na / niSquared_expect)

    assert results[0:17] == approx((Na, Nd, ni_expect, niSquared_expect, i_width,
                             Lp, Ln, n_width, p_width,
                             0, 0, GaAs_n.hole_mobility * kb*options.T / q,
                                    GaAs_p.electron_mobility * kb * options.T / q,
                              GaAs_n.permittivity, 0, 2, Vbi_expect))

    assert results[17] == 'np'


def test_process_junction_pin():
    from solcore.analytic_solar_cells.depletion_approximation import process_junction
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    Nd = randrange(1, 9) * 1e17
    Na = randrange(1, 10) * 1e18

    options = State()
    options.T = randrange(1, 350)

    Lp = randrange(5, 10000) * 1e-9  # Diffusion length
    Ln = randrange(5, 10000) * 1e-9  # Diffusion length

    GaAs_n = material("GaAs")(Nd=Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na=Na, electron_diffusion_length=Lp)
    GaAs_i = material("GaAs")()

    p_width = randrange(500, 1000)
    n_width = randrange(3000, 5000)
    i_width = randrange(100, 300)*1e-9

    test_junc = Junction([Layer(p_width, GaAs_p, role="emitter"),
                          Layer(i_width, GaAs_i, role="intrinsic"),
                          Layer(n_width, GaAs_n, role="base")])

    results = process_junction(test_junc, options)

    ni_expect = GaAs_n.ni
    niSquared_expect = ni_expect ** 2

    Vbi_expect = (kb * options.T / q) * np.log(Nd * Na / niSquared_expect)

    assert results[0:17] == approx((Na, Nd, ni_expect, niSquared_expect, i_width,
                             Lp, Ln, n_width, p_width,
                             0, 0, GaAs_n.hole_mobility * kb*options.T / q,
                                    GaAs_p.electron_mobility * kb * options.T / q,
                              GaAs_n.permittivity, 0, 2, Vbi_expect))

    assert results[17] == 'pn'



def test_process_junction_set_in_junction():
    from solcore.analytic_solar_cells.depletion_approximation import process_junction
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State

    options = State()
    options.T = randrange(1, 350)

    Lp = randrange(5, 10000) * 1e-9  # Diffusion length
    Ln = randrange(5, 10000) * 1e-9  # Diffusion length

    sn = randrange(1, 1000)
    sp = randrange(1, 1000)

    se = randrange(1, 20)*vacuum_permittivity

    GaAs_n = material("GaAs")()
    GaAs_p = material("GaAs")()
    GaAs_i = material("GaAs")()

    p_width = randrange(500, 1000)
    n_width = randrange(3000, 5000)
    i_width = randrange(100, 300)*1e-9

    mun = randrange(1, 1e5)*1e-5
    mup = randrange(1, 1e5)*1e-5

    Vbi = randrange(0,3)

    test_junc = Junction([Layer(p_width, GaAs_p, role="emitter"),
                          Layer(i_width, GaAs_i, role="intrinsic"),
                          Layer(n_width, GaAs_n, role="base")],
                         sn = sn, sp = sp, permittivity = se,
                         ln = Ln, lp= Lp,
                         mup = mup, mun = mun, Vbi = Vbi)

    results = process_junction(test_junc, options)

    ni_expect = GaAs_n.ni
    niSquared_expect = ni_expect ** 2

    assert results[0:17] == approx((GaAs_n.Na, GaAs_p.Nd, ni_expect, niSquared_expect, i_width,
                             Ln, Lp, n_width, p_width,
                             sn, sp, mun* kb*options.T / q, mup* kb*options.T / q,
                              se, 0, 2, Vbi))

    assert results[17] == 'pn'


def test_get_depletion_widths():
    from solcore.analytic_solar_cells.depletion_approximation import get_depletion_widths
    from solcore.structure import Junction

    xnm = randrange(2, 1000)
    xi  = randrange(0, 1000)*1e-9

    Vbi = randrange(1, 50)*1e-1
    es = randrange(1, 20)*vacuum_permittivity
    Na = randrange(1, 1e5)*1e19
    Nd = randrange(1, 1e5)*1e19

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

    xi  = randrange(0, 1000)*1e-9

    Vbi = randrange(1, 50)*1e-1
    es = randrange(1, 20)*vacuum_permittivity
    Na = randrange(1, 1e5)*1e19
    Nd = randrange(1, 1e5)*1e19

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

    wn = randrange(1,100)
    wp = randrange(1,100)
    test_junc = Junction(wn=wn, wp=wp)

    wn_r, wp_r = get_depletion_widths(test_junc, 0, 0, 0, 0, 0, 0)

    assert wn_r == approx(wn)
    assert wp_r == approx(wp)


def test_dark_iv_depletion_pn():
    from solcore.analytic_solar_cells.depletion_approximation import iv_depletion, get_depletion_widths, get_j_dark, get_Jsrh, process_junction
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State
    from scipy.interpolate import interp1d

    Nd = randrange(1, 9) * 1e17
    Na = randrange(1, 10) * 1e18

    options = State()

    options.T = 270
    options.wavelength = np.linspace(290, 700, 150)*1e-9
    options.internal_voltages =  np.linspace(-6, 4, 200)
    options.light_iv = False

    Lp = randrange(5, 500) * 1e-9  # Diffusion length
    Ln = randrange(5, 500) * 1e-9  # Diffusion length

    GaAs_n = material("GaAs")(Nd=Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na=Na, electron_diffusion_length=Lp)

    p_width = randrange(500, 1000)*1e-9
    n_width = randrange(3000, 5000)*1e-9

    test_junc = Junction([Layer(p_width, GaAs_p, role="emitter"),
                          Layer(n_width, GaAs_n, role="base")])

    test_junc.voltage = options.internal_voltages
    T = options.T

    Na, Nd, ni, niSquared, xi, ln, lp, xn, xp, sn, sp, dn, dp, es, id_top, id_bottom, \
    Vbi, pn_or_np = process_junction(test_junc, options)

    kbT = kb * T

    # And now we account for the possible applied voltage, which can be, at most, equal to Vbi
    V = np.where(test_junc.voltage < Vbi - 0.001, test_junc.voltage, Vbi - 0.001)

    wn, wp = get_depletion_widths(test_junc, es, Vbi, V, Na, Nd, xi)

    w = wn + wp + xi

    # Now it is time to calculate currents
    l_top, l_bottom = ln, lp
    x_top, x_bottom = xp, xn
    w_top, w_bottom = wp, wn
    s_top, s_bottom = sp, sn
    d_top, d_bottom = dp, dn
    min_top, min_bot = niSquared / Na, niSquared / Nd

    #print(min_bot, min_top)
    JtopDark = get_j_dark(x_top, w_top, l_top, s_top, d_top, V, min_top, T)
    JbotDark = get_j_dark(x_bottom, w_bottom, l_bottom, s_bottom, d_bottom, V, min_bot, T)

    # hereby we define the subscripts to refer to the layer in which the current is generated:
    JnDark, JpDark = JbotDark, JtopDark

    # These might not be the right lifetimes. Actually, they are not as they include all recombination processes, not
    # just SRH recombination, which is what the equation in Jenny, p159 refers to. Let´ leave them, for now.
    lifetime_n = ln ** 2 / dn
    lifetime_p = lp ** 2 / dp  # Jenny p163

    # Here we use the full version of the SRH recombination term as calculated by Sah et al. Works for positive bias
    # and moderately negative ones.

    Jrec = get_Jsrh(ni, V, Vbi, lifetime_p, lifetime_n, w, kbT)

    J_sc_top = 0
    J_sc_bot = 0
    J_sc_scr = 0

    test_junc.current = Jrec + JnDark + JpDark + V / 1e14- J_sc_top - J_sc_bot - J_sc_scr
    iv = interp1d(test_junc.voltage, test_junc.current, kind='linear', bounds_error=False, assume_sorted=True,
                           fill_value=(test_junc.current[0], test_junc.current[-1]), copy=True)

    iv_depletion(test_junc, options)

    #print(test_junc.iv(options.internal_voltages))
    assert test_junc.iv(options.internal_voltages) == approx(iv(options.internal_voltages), nan_ok=True)



def test_dark_iv_depletion_np():
    from solcore.analytic_solar_cells.depletion_approximation import iv_depletion, get_depletion_widths, get_j_dark, get_Jsrh, process_junction
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore.state import State
    from scipy.interpolate import interp1d

    Nd = randrange(1, 10) * 1e18
    Na = randrange(1, 9) * 1e17

    options = State()

    options.T = 270
    options.wavelength = np.linspace(290, 700, 150)*1e-9
    options.internal_voltages =  np.linspace(-6, 4, 200)
    options.light_iv = False

    Lp = randrange(5, 500) * 1e-9  # Diffusion length
    Ln = randrange(5, 500) * 1e-9  # Diffusion length

    GaAs_n = material("GaAs")(Nd=Nd, hole_diffusion_length=Ln)
    GaAs_p = material("GaAs")(Na=Na, electron_diffusion_length=Lp)

    n_width = randrange(500, 1000)*1e-9
    p_width = randrange(3000, 5000)*1e-9

    test_junc = Junction([Layer(n_width, GaAs_n, role="emitter"),
                          Layer(p_width, GaAs_p, role="base")])

    test_junc.voltage = options.internal_voltages
    T = options.T

    Na, Nd, ni, niSquared, xi, ln, lp, xn, xp, sn, sp, dn, dp, es, id_top, id_bottom, \
    Vbi, pn_or_np = process_junction(test_junc, options)

    kbT = kb * T

    # And now we account for the possible applied voltage, which can be, at most, equal to Vbi
    V = np.where(test_junc.voltage < Vbi - 0.001, test_junc.voltage, Vbi - 0.001)

    wn, wp = get_depletion_widths(test_junc, es, Vbi, V, Na, Nd, xi)

    w = wn + wp + xi

    # Now it is time to calculate currents
    l_bottom, l_top = ln, lp
    x_bottom, x_top = xp, xn
    w_bottom, w_top = wp, wn
    s_bottom, s_top = sp, sn
    d_bottom, d_top = dp, dn
    min_bot, min_top = niSquared / Na, niSquared / Nd

    #print(min_bot, min_top)
    JtopDark = get_j_dark(x_top, w_top, l_top, s_top, d_top, V, min_top, T)
    JbotDark = get_j_dark(x_bottom, w_bottom, l_bottom, s_bottom, d_bottom, V, min_bot, T)

    # hereby we define the subscripts to refer to the layer in which the current is generated:
    JpDark, JnDark = JbotDark, JtopDark

    # These might not be the right lifetimes. Actually, they are not as they include all recombination processes, not
    # just SRH recombination, which is what the equation in Jenny, p159 refers to. Let´ leave them, for now.
    lifetime_n = ln ** 2 / dn
    lifetime_p = lp ** 2 / dp  # Jenny p163

    # Here we use the full version of the SRH recombination term as calculated by Sah et al. Works for positive bias
    # and moderately negative ones.

    Jrec = get_Jsrh(ni, V, Vbi, lifetime_p, lifetime_n, w, kbT)

    J_sc_top = 0
    J_sc_bot = 0
    J_sc_scr = 0

    test_junc.current = Jrec + JnDark + JpDark + V / 1e14- J_sc_top - J_sc_bot - J_sc_scr
    iv = interp1d(test_junc.voltage, test_junc.current, kind='linear', bounds_error=False, assume_sorted=True,
                           fill_value=(test_junc.current[0], test_junc.current[-1]), copy=True)

    iv_depletion(test_junc, options)

    #print(test_junc.iv(options.internal_voltages))
    assert test_junc.iv(options.internal_voltages) == approx(iv(options.internal_voltages), nan_ok=True)


