from pytest import approx
from random import randrange
import numpy as np

from solcore.constants import kb, q, vacuum_permittivity

# get_j_top and get_j_bot can be one function
# separate function which calculates the depletion widths

# w is not a constant but has the same length as V! because it's a function of V


## IV

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

    assert result == approx(expected)


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

    assert result == approx(expected)


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

    assert result == approx(expected)


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

