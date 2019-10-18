from pytest import approx, mark
from random import randrange
import numpy as np

from solcore.constants import kb, q, vacuum_permittivity

# get_j_top and get_j_bot can be one function
# separate function which calculates the depletion widths

# w is not a constant but has the same length as V! because it's a function of V

# pre-processing for both qe and iv done in separate function


## IV

@mark.skip
def test_get_j_top():
    from solcore.analytic_solar_cells.depletion_approximation import get_j_top, get_depletion_widths

    xnm = randrange(2, 1000)
    x = xnm*1e-9

    xi  = randrange(0, 1000)*1e-9

    es = randrange(1, 20)*vacuum_permittivity
    l = randrange(5, 10000)*1e-9
    s = randrange(1, 1000)
    d = randrange(1, 1e5)*1e-5
    Vbi = randrange(1, 50)*1e-1
    minor = randrange(1, 1e6)*1e8
    T = randrange(1, 300)
    Na = randrange(1, 1e5)*1e19
    Nd = randrange(1, 1e5)*1e19

    V = np.linspace(-6, 4, 20)
    V = np.where(V < Vbi - 0.001, V, Vbi - 0.001)

    wn, wp = get_depletion_widths(es, Vbi, V, Na, Nd, xi, one_sided=False)
    w = wn + wp + xi

    expected = (q*d*minor/l)*(np.exp(q*V/(kb*T))-1)*((((s*l)/d)*np.cosh((x-w)/l)+np.sinh((x-w)/l))/
                                                  (((s*l)/d)*np.sinh((x-w)/l)+np.cosh((x-w)/l)))

    result = get_j_top(x, w, l, s, d, V, minor, T)

    assert result == approx(expected)

@mark.skip
def test_get_j_bot():
    from solcore.analytic_solar_cells.depletion_approximation import get_j_bot, get_depletion_widths

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

    wn, wp = get_depletion_widths(es, Vbi, V, Na, Nd, xi, one_sided=False)
    w = wn + wp + xi

    expected = (q*d*minor/l)*(np.exp(q*V/(kb*T))-1)*((((s*l)/d)*np.cosh((x-w)/l)+np.sinh((x-w)/l))/
                                                  (((s*l)/d)*np.sinh((x-w)/l)+np.cosh((x-w)/l)))

    result = get_j_bot(x, w, l, s, d, V, minor, T)

    assert result == approx(expected)

def test_get_j_dark():
    from solcore.analytic_solar_cells.depletion_approximation import get_j_dark, get_depletion_widths

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

    wn, wp = get_depletion_widths(es, Vbi, V, Na, Nd, xi, one_sided=False)
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
    from solcore.analytic_solar_cells.depletion_approximation import factor, forward, get_depletion_widths

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

    wn, wp = get_depletion_widths(es, Vbi, V, Na, Nd, xi, one_sided=False)
    w = wn + wp + xi
    w = w[m]

    f_b = factor(V, Vbi, tp, tn, kT, dEt)
    expected = 2 * q * ni * w / np.sqrt(tn * tp) * \
          np.sinh(q*V / (2 * kT)) / (q * (Vbi - V) / kT) * f_b

    result = forward(ni, V, Vbi, tp, tn, w, kT, dEt)

    assert result == approx(expected, nan_ok=True)


# this test failed one, not sure on conditions (rng)
def test_get_J_srh():
    from solcore.analytic_solar_cells.depletion_approximation import forward, get_Jsrh, get_depletion_widths
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

    wn, wp = get_depletion_widths(es, Vbi, V, Na, Nd, xi, one_sided=False)
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


def test_qe_depletion():
    from solcore.structure import Junction, Layer
    from solcore import si, material
    from solcore.solar_cell import SolarCell
    from solcore.solar_cell_solver import solar_cell_solver
    from solcore.state import State
    from solcore.light_source import LightSource
    from solcore.analytic_solar_cells import qe_depletion


    AlInP = material("AlInP")
    InGaP = material("GaInP")
    window_material = AlInP(Al=0.52)
    top_cell_n_material = InGaP(In=0.48, Nd=si(2e18, "cm-3"), hole_diffusion_length=si("300nm"))
    top_cell_p_material = InGaP(In=0.48, Na=si(1.5e17, "cm-3"), electron_diffusion_length=si("2um"))

    for mat in [top_cell_n_material, top_cell_p_material]:
        mat.permittivity = 11.75 * vacuum_permittivity

    test_junc = SolarCell([Junction([Layer(si("25nm"), material=window_material, role='window'),
                  Layer(si("80nm"), material=top_cell_n_material, role='emitter'),
                  Layer(si("700nm"), material=top_cell_p_material, role='base'),
                 ], sn=1, sp=1, kind='DA')])


    light_source = LightSource(source_type="standard", version="AM1.5g")


    options = State()

    options.T = 270
    options.wavelength = np.linspace(290, 700, 150)*1e-9
    options.light_source = light_source
    options.optics_method = 'TMM'

    solar_cell_solver(test_junc, 'optics', options)

    qe_depletion(test_junc[0], options)

    print('CHECK', test_junc[0].eqe(options.wavelength)[0], 6.28783762e-02)

    assert test_junc[0].eqe(options.wavelength) == approx(np.array([6.28783762e-02, 6.56312772e-02, 6.77987661e-02, 6.92704234e-02,
       7.02852937e-02, 7.07563896e-02, 7.09991841e-02, 7.12213622e-02,
       7.15543586e-02, 7.22625861e-02, 7.35141632e-02, 7.56851451e-02,
       7.88330859e-02, 8.31660933e-02, 8.87859279e-02, 9.56381165e-02,
       1.03869319e-01, 1.13323257e-01, 1.23917013e-01, 1.35609558e-01,
       1.48368991e-01, 1.61906966e-01, 1.76290253e-01, 1.91191711e-01,
       2.06569878e-01, 2.22295620e-01, 2.38440155e-01, 2.55491133e-01,
       2.72770422e-01, 2.89933735e-01, 3.06859046e-01, 3.23813615e-01,
       3.41419300e-01, 3.59100080e-01, 3.77057108e-01, 3.94879261e-01,
       4.12908679e-01, 4.28412383e-01, 4.43286335e-01, 4.56456921e-01,
       4.69114578e-01, 4.80615083e-01, 4.91420477e-01, 5.01658044e-01,
       5.11259849e-01, 5.19927699e-01, 5.28129905e-01, 5.35536613e-01,
       5.42547029e-01, 5.49438168e-01, 5.56056603e-01, 5.62501174e-01,
       5.68419898e-01, 5.74110234e-01, 5.79389925e-01, 5.84708510e-01,
       5.89888234e-01, 5.94772225e-01, 5.99431212e-01, 6.03708440e-01,
       6.08135459e-01, 6.12718701e-01, 6.16815591e-01, 6.21206834e-01,
       6.25496073e-01, 6.29321181e-01, 6.33532612e-01, 6.37802695e-01,
       6.41566535e-01, 6.45919749e-01, 6.50107008e-01, 6.54232088e-01,
       6.58823060e-01, 6.63431720e-01, 6.67830990e-01, 6.72797625e-01,
       6.77010858e-01, 6.81324008e-01, 6.84306321e-01, 6.86362457e-01,
       6.89366529e-01, 6.91952134e-01, 6.93877708e-01, 6.96142534e-01,
       6.97826330e-01, 6.99477241e-01, 7.00426721e-01, 7.01736362e-01,
       7.02399588e-01, 7.03029068e-01, 7.03628581e-01, 7.04239271e-01,
       7.03844750e-01, 7.03039162e-01, 7.02070760e-01, 7.01006923e-01,
       6.99600578e-01, 6.98003177e-01, 6.95714479e-01, 6.93114506e-01,
       6.90294959e-01, 6.87285580e-01, 6.85180955e-01, 6.85184524e-01,
       6.85178042e-01, 6.85582386e-01, 6.85633019e-01, 6.85618741e-01,
       6.85492757e-01, 6.85422334e-01, 6.85753348e-01, 6.85649853e-01,
       6.85460449e-01, 6.85349526e-01, 6.84279221e-01, 6.81840211e-01,
       6.80227683e-01, 6.77785518e-01, 6.76303705e-01, 6.74463213e-01,
       6.71749698e-01, 6.67701114e-01, 6.64125641e-01, 6.59566124e-01,
       6.54211689e-01, 6.48485996e-01, 6.40837924e-01, 6.28725582e-01,
       6.16941627e-01, 5.99958985e-01, 5.74994985e-01, 5.40382562e-01,
       4.96420617e-01, 4.28038672e-01, 3.37684277e-01, 2.11523519e-01,
       1.01066369e-01, 7.86512303e-02, 5.59183882e-02, 3.32479392e-04,
       8.37778441e-05, 8.16799923e-06, 3.92426495e-06, 4.13542516e-32,
       4.53130213e-32, 5.01468542e-32, 4.47113089e-32, 4.45697014e-32,
       4.18585110e-32, 4.34094702e-32]))

    assert test_junc[0].iqe(options.wavelength) == approx(np.array([0.97367538, 0.97375406, 0.97383056, 0.97390309, 0.97397857,
       0.97403364, 0.97407053, 0.97410751, 0.97413875, 0.97416899,
       0.97419971, 0.97423027, 0.97426069, 0.9742916 , 0.9743225 ,
       0.97434943, 0.97437674, 0.97440381, 0.97443064, 0.97445785,
       0.97449835, 0.9745412 , 0.97458427, 0.97462849, 0.97467567,
       0.97472903, 0.97478323, 0.97489724, 0.97503793, 0.9751862 ,
       0.97534233, 0.97550724, 0.97574234, 0.97605329, 0.9764076 ,
       0.97692643, 0.97781994, 0.97845083, 0.97911714, 0.97975648,
       0.98045403, 0.98106622, 0.98172135, 0.98241978, 0.98289109,
       0.98332447, 0.98376793, 0.98422025, 0.98467967, 0.98509539,
       0.98541639, 0.98575053, 0.98607405, 0.98639245, 0.98664093,
       0.98687221, 0.98709617, 0.98731145, 0.98751654, 0.98771363,
       0.98784935, 0.9879432 , 0.98802966, 0.98810915, 0.98818312,
       0.98825134, 0.98830842, 0.98837472, 0.98844157, 0.98849472,
       0.98853316, 0.98855585, 0.98856173, 0.98854969, 0.98851863,
       0.98846743, 0.98839548, 0.98830227, 0.98818634, 0.98804673,
       0.98800044, 0.98796452, 0.98792435, 0.98788994, 0.98785128,
       0.98781176, 0.98776822, 0.98773013, 0.9876849 , 0.98764206,
       0.98760152, 0.98755707, 0.98746511, 0.98735737, 0.98724168,
       0.98712689, 0.98700426, 0.98687672, 0.98674174, 0.98660731,
       0.98646298, 0.98631657, 0.9861978 , 0.98618513, 0.98616995,
       0.98615485, 0.98613742, 0.98612484, 0.98610744, 0.98609254,
       0.98607767, 0.98606282, 0.98604806, 0.98603327, 0.98598835,
       0.98589035, 0.9858101 , 0.98571717, 0.98566245, 0.9855933 ,
       0.98549286, 0.98535319, 0.98523557, 0.98508699, 0.98492475,
       0.98475758, 0.98454257, 0.98422515, 0.98393297, 0.98354185,
       0.98302706, 0.98239727, 0.98171046, 0.98081462, 0.97985356,
       0.97878668, 0.97802882, 0.97789095, 0.97775466, 0.97743861,
       0.97743726, 0.97743687, 0.97743599,        np.inf,        np.inf,
              np.inf,        np.inf,        np.inf,        np.inf,        np.inf]))

    assert test_junc[0].eqe_base(options.wavelength) == approx(np.array([2.15645289e-09, 3.74405444e-09, 6.09307808e-09, 9.28198091e-09,
       1.38941539e-08, 1.82545978e-08, 2.17099703e-08, 2.57260906e-08,
       2.96794808e-08, 3.41083720e-08, 3.93934418e-08, 4.59400226e-08,
       5.40848898e-08, 6.43563971e-08, 7.72866748e-08, 9.20518745e-08,
       1.10362754e-07, 1.32710448e-07, 1.59700782e-07, 1.92051289e-07,
       2.40788643e-07, 3.02680542e-07, 3.78845530e-07, 4.71345674e-07,
       5.88444279e-07, 7.41390420e-07, 9.28839820e-07, 1.35848469e-06,
       2.07902321e-06, 3.14835818e-06, 4.71158366e-06, 6.99200940e-06,
       1.15302532e-05, 2.06976235e-05, 3.66976919e-05, 7.71613934e-05,
       2.17998331e-04, 4.07071098e-04, 7.26455695e-04, 1.19148412e-03,
       1.93313118e-03, 2.85230426e-03, 4.17904139e-03, 6.08199838e-03,
       7.76391938e-03, 9.62942482e-03, 1.18968491e-02, 1.46343320e-02,
       1.79389719e-02, 2.15037768e-02, 2.48166786e-02, 2.85775600e-02,
       3.28182776e-02, 3.76063832e-02, 4.19170396e-02, 4.64584918e-02,
       5.14135652e-02, 5.67965039e-02, 6.26409759e-02, 6.90846191e-02,
       7.43319958e-02, 7.85564402e-02, 8.28113874e-02, 8.72058157e-02,
       9.17618047e-02, 9.64262688e-02, 1.01332480e-01, 1.06880513e-01,
       1.14034376e-01, 1.21677641e-01, 1.29690315e-01, 1.38099279e-01,
       1.47031023e-01, 1.56410309e-01, 1.66187516e-01, 1.76565867e-01,
       1.87153120e-01, 1.98146504e-01, 2.09183151e-01, 2.20331409e-01,
       2.24503639e-01, 2.27761991e-01, 2.30814526e-01, 2.33989215e-01,
       2.36978452e-01, 2.39962259e-01, 2.42710091e-01, 2.45584270e-01,
       2.48234069e-01, 2.50868975e-01, 2.53492222e-01, 2.56117263e-01,
       2.60800990e-01, 2.65930230e-01, 2.71012651e-01, 2.76064974e-01,
       2.80986203e-01, 2.85828776e-01, 2.90374832e-01, 2.94764534e-01,
       2.99033119e-01, 3.03182025e-01, 3.06551016e-01, 3.07093330e-01,
       3.07631142e-01, 3.08349118e-01, 3.08904574e-01, 3.09427311e-01,
       3.09899376e-01, 3.10391229e-01, 3.11060716e-01, 3.11533175e-01,
       3.11963848e-01, 3.12425791e-01, 3.13572414e-01, 3.15706350e-01,
       3.17641732e-01, 3.19551222e-01, 3.20699408e-01, 3.22037652e-01,
       3.23830300e-01, 3.26173106e-01, 3.27958219e-01, 3.30057409e-01,
       3.31999121e-01, 3.33715332e-01, 3.35523638e-01, 3.37303845e-01,
       3.38038980e-01, 3.37574968e-01, 3.34305238e-01, 3.25909305e-01,
       3.10600033e-01, 2.79682354e-01, 2.30077483e-01, 1.50267042e-01,
       7.37753996e-02, 5.76874883e-02, 4.12040010e-02, 2.47577771e-04,
       6.23872022e-05, 6.08258194e-06, 2.92234450e-06, 4.07658114e-32,
       4.46682507e-32, 4.94333018e-32, 4.40751003e-32, 4.39355078e-32,
       4.12628955e-32, 4.27917858e-32]))

    assert test_junc[0].eqe_emitter(options.wavelength) == approx(np.array([6.27553737e-02, 6.54767793e-02, 6.76104699e-02, 6.90479502e-02,
       7.00252606e-02, 7.04679885e-02, 7.06913399e-02, 7.08932597e-02,
       7.12077158e-02, 7.18956615e-02, 7.31232507e-02, 7.52640225e-02,
       7.83744653e-02, 8.26606090e-02, 8.82226811e-02, 9.50088134e-02,
       1.03160830e-01, 1.12522440e-01, 1.23010294e-01, 1.34582659e-01,
       1.47188395e-01, 1.60550633e-01, 1.74736651e-01, 1.89420495e-01,
       2.04552910e-01, 2.19996989e-01, 2.35831313e-01, 2.52360832e-01,
       2.68960805e-01, 2.85328124e-01, 3.01330179e-01, 3.17209131e-01,
       3.33225952e-01, 3.48633708e-01, 3.63761733e-01, 3.76940606e-01,
       3.85990089e-01, 3.93835746e-01, 3.99676032e-01, 4.03280377e-01,
       4.04609926e-01, 4.05173947e-01, 4.03527875e-01, 3.99616662e-01,
       3.98415972e-01, 3.96605988e-01, 3.93648778e-01, 3.89280893e-01,
       3.83782928e-01, 3.78530372e-01, 3.74664903e-01, 3.70193207e-01,
       3.64878763e-01, 3.58920327e-01, 3.54268657e-01, 3.49681843e-01,
       3.44705738e-01, 3.39252539e-01, 3.33370196e-01, 3.26796051e-01,
       3.22457908e-01, 3.19867208e-01, 3.17045413e-01, 3.14385484e-01,
       3.11588217e-01, 3.08478812e-01, 3.05475012e-01, 3.01910997e-01,
       2.96398942e-01, 2.91006896e-01, 2.85385525e-01, 2.79588416e-01,
       2.73840040e-01, 2.67943612e-01, 2.61816638e-01, 2.55766762e-01,
       2.49340553e-01, 2.42889906e-01, 2.35866532e-01, 2.28460438e-01,
       2.27033338e-01, 2.26070747e-01, 2.24894702e-01, 2.23835255e-01,
       2.22592043e-01, 2.21345267e-01, 2.19883083e-01, 2.18545515e-01,
       2.17013859e-01, 2.15485208e-01, 2.13959149e-01, 2.12446998e-01,
       2.08949463e-01, 2.04956987e-01, 2.00952882e-01, 1.96964558e-01,
       1.92921172e-01, 1.88869931e-01, 1.84682001e-01, 1.80471559e-01,
       1.76260138e-01, 1.72063792e-01, 1.68861978e-01, 1.68528759e-01,
       1.68192605e-01, 1.67960311e-01, 1.67643033e-01, 1.67313983e-01,
       1.66956573e-01, 1.66617120e-01, 1.66378123e-01, 1.66034013e-01,
       1.65671134e-01, 1.65330238e-01, 1.64073166e-01, 1.61508461e-01,
       1.59503602e-01, 1.57095772e-01, 1.55647149e-01, 1.53904743e-01,
       1.51452127e-01, 1.48009021e-01, 1.45151171e-01, 1.41629585e-01,
       1.37820811e-01, 1.33978522e-01, 1.29154605e-01, 1.22188001e-01,
       1.16022040e-01, 1.08046292e-01, 9.78349346e-02, 8.58661520e-02,
       7.32157386e-02, 5.72944508e-02, 4.07012349e-02, 2.26532425e-02,
       9.93406613e-03, 7.60897211e-03, 5.32560941e-03, 3.05263356e-05,
       7.69070684e-06, 7.49769997e-07, 3.60217832e-07, 5.88440218e-34,
       6.44770563e-34, 7.13552408e-34, 6.36208644e-34, 6.34193675e-34,
       5.95615452e-34, 6.17684448e-34]))

    assert test_junc[0].eqe_scr(options.wavelength) == approx(np.array([1.23000436e-04, 1.54494130e-04, 1.88290129e-04, 2.22463916e-04,
       2.60019248e-04, 2.88382916e-04, 3.07822492e-04, 3.28076797e-04,
       3.46613144e-04, 3.66890460e-04, 3.90873062e-04, 4.21076704e-04,
       4.58566598e-04, 5.05420002e-04, 5.63169459e-04, 6.29211135e-04,
       7.08378462e-04, 8.00684097e-04, 9.06560110e-04, 1.02670710e-03,
       1.18035487e-03, 1.35602946e-03, 1.55322287e-03, 1.77074422e-03,
       2.01638010e-03, 2.29788887e-03, 2.60791356e-03, 3.12894155e-03,
       3.80753839e-03, 4.60246208e-03, 5.52415508e-03, 6.59749231e-03,
       8.18181722e-03, 1.04456743e-02, 1.32586774e-02, 1.78614935e-02,
       2.67005917e-02, 3.41695657e-02, 4.28838477e-02, 5.19850597e-02,
       6.25715213e-02, 7.25888321e-02, 8.37135603e-02, 9.59593836e-02,
       1.05079958e-01, 1.13692286e-01, 1.22584278e-01, 1.31621388e-01,
       1.40825130e-01, 1.49404020e-01, 1.56575021e-01, 1.63730406e-01,
       1.70722857e-01, 1.77583524e-01, 1.83204229e-01, 1.88568176e-01,
       1.93768931e-01, 1.98723182e-01, 2.03420040e-01, 2.07827770e-01,
       2.11345555e-01, 2.14295053e-01, 2.16958791e-01, 2.19615534e-01,
       2.22146051e-01, 2.24416100e-01, 2.26725120e-01, 2.29011185e-01,
       2.31133217e-01, 2.33235213e-01, 2.35031168e-01, 2.36544393e-01,
       2.37951997e-01, 2.39077799e-01, 2.39826836e-01, 2.40464996e-01,
       2.40517185e-01, 2.40287599e-01, 2.39256638e-01, 2.37570610e-01,
       2.37829552e-01, 2.38119396e-01, 2.38168480e-01, 2.38318065e-01,
       2.38255834e-01, 2.38169715e-01, 2.37833547e-01, 2.37606576e-01,
       2.37151660e-01, 2.36674885e-01, 2.36177211e-01, 2.35675010e-01,
       2.34094297e-01, 2.32151945e-01, 2.30105226e-01, 2.27977391e-01,
       2.25693203e-01, 2.23304470e-01, 2.20657645e-01, 2.17878413e-01,
       2.15001702e-01, 2.12039764e-01, 2.09767962e-01, 2.09562435e-01,
       2.09354295e-01, 2.09272958e-01, 2.09085412e-01, 2.08877447e-01,
       2.08636808e-01, 2.08413985e-01, 2.08314509e-01, 2.08082666e-01,
       2.07825466e-01, 2.07593498e-01, 2.06633642e-01, 2.04625400e-01,
       2.03082349e-01, 2.01138524e-01, 1.99957149e-01, 1.98520817e-01,
       1.96467272e-01, 1.93518987e-01, 1.91016250e-01, 1.87879130e-01,
       1.84391757e-01, 1.80792141e-01, 1.76159681e-01, 1.69233736e-01,
       1.62880607e-01, 1.54337725e-01, 1.42854812e-01, 1.28607106e-01,
       1.12604846e-01, 9.10618672e-02, 6.69055594e-02, 3.86032338e-02,
       1.73569035e-02, 1.33547699e-02, 9.38877780e-03, 5.43752850e-05,
       1.36999351e-05, 1.33564729e-06, 6.41702618e-07, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00]))



# fix infs when there is no light source at that wavelength
# infs at long wavelength in IQE


def test_iv_depletion():
    from solcore.structure import Junction, Layer
    from solcore import si, material
    from solcore.solar_cell import SolarCell
    from solcore.solar_cell_solver import solar_cell_solver
    from solcore.state import State
    from solcore.light_source import LightSource
    from solcore.analytic_solar_cells import iv_depletion


    AlInP = material("AlInP")
    InGaP = material("GaInP")
    window_material = AlInP(Al=0.52)
    top_cell_n_material = InGaP(In=0.48, Nd=si(2e18, "cm-3"), hole_diffusion_length=si("300nm"))
    top_cell_p_material = InGaP(In=0.48, Na=si(1.5e17, "cm-3"), electron_diffusion_length=si("2um"))

    for mat in [top_cell_n_material, top_cell_p_material]:
        mat.permittivity = 11.75 * vacuum_permittivity

    test_junc = SolarCell([Junction([Layer(si("25nm"), material=window_material, role='window'),
                  Layer(si("80nm"), material=top_cell_n_material, role='emitter'),
                  Layer(si("700nm"), material=top_cell_p_material, role='base'),
                 ], sn=1, sp=1, kind='DA')])


    light_source = LightSource(source_type="standard", version="AM1.5g")


    options = State()

    options.T = 270
    options.wavelength = np.linspace(290, 700, 150)*1e-9
    options.light_source = light_source
    options.optics_method = 'TMM'
    options.internal_voltages =  np.linspace(-6, 4, 200)
    options.light_iv = True

    solar_cell_solver(test_junc, 'optics', options)

    iv_depletion(test_junc[0], options)
    Vs = np.linspace(0, 1.3, 100)


    #plt.figure()
    #plt.plot(Vs, test_junc[0].iv(Vs))
    #plt.xlim(0, 1.5)
    #plt.ylim(-125, 0)
    #plt.show()
    print('CHECK', test_junc[0].iv(Vs), -108.58213025)

    assert test_junc[0].iv(Vs) == approx(np.array([-108.58213025, -108.58213025, -108.58213025, -108.58213025,
       -108.58213025, -108.58213025, -108.58213025, -108.58213025,
       -108.58213025, -108.58213025, -108.58213025, -108.58213025,
       -108.58213025, -108.58213024, -108.58213024, -108.58213024,
       -108.58213023, -108.58213023, -108.58213022, -108.5821302 ,
       -108.58213019, -108.58213017, -108.58213014, -108.58213009,
       -108.58213004, -108.58213   , -108.58212988, -108.58212974,
       -108.5821296 , -108.58212945, -108.58212905, -108.58212863,
       -108.5821282 , -108.5821277 , -108.58212642, -108.58212514,
       -108.58212387, -108.58212191, -108.58211806, -108.58211422,
       -108.58211037, -108.58210314, -108.58209155, -108.58207995,
       -108.58206836, -108.58204252, -108.58200755, -108.58197257,
       -108.58193759, -108.58184735, -108.58174171, -108.58163606,
       -108.58153042, -108.58122051, -108.58090101, -108.58058152,
       -108.58017883, -108.5792111 , -108.57824338, -108.57727565,
       -108.57571444, -108.57277832, -108.56984219, -108.56690607,
       -108.56112694, -108.5522015 , -108.54327606, -108.53435062,
       -108.51359116, -108.48639924, -108.45920732, -108.4320154 ,
       -108.3589581 , -108.27590284, -108.19284758, -108.10979232,
       -107.85630752, -107.60182878, -107.34735003, -107.00442706,
       -106.22152977, -105.43863249, -104.65573521, -103.31430026,
       -100.89129296,  -98.46828566,  -96.04527836,  -90.97466173,
        -83.39788068,  -75.82109963,  -68.24431857,  -49.25528419,
        -25.06426601,   -0.87324784,   23.31777033,   96.25210425,
        177.11681133,  257.98151841,  346.07933073,  644.525936  ]))

    assert test_junc[0].current == approx(np.array([-1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582130e+02,
       -1.08582130e+02, -1.08582130e+02, -1.08582130e+02, -1.08582129e+02,
       -1.08582128e+02, -1.08582123e+02, -1.08582108e+02, -1.08582064e+02,
       -1.08581930e+02, -1.08581526e+02, -1.08580303e+02, -1.08576600e+02,
       -1.08565364e+02, -1.08531208e+02, -1.08427149e+02, -1.08109311e+02,
       -1.07135465e+02, -1.04139453e+02, -9.48670243e+01, -6.58719983e+01,
        2.67028297e+01,  3.36158029e+02,  1.47826139e+03,  6.52808449e+03,
        3.48823191e+04,  2.32541710e+05,  1.83901183e+06,  1.57703678e+07,
        1.38102361e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08,
        5.19416018e+08,  5.19416018e+08,  5.19416018e+08,  5.19416018e+08]))