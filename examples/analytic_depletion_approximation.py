import numpy as np
from solcore import material
from solcore.constants import q, kb
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp, quad_vec
from functools import partial
from solcore.interpolate import interp1d
import warnings
from pytest import approx

GaAs = material("GaAs")()

# make n and p the same
# width = 850e-9
width_top = 300e-9
width_bottom = 1000000e-9


def get_J_sc_diffusion_vs_WL(xa, xb, g, D, L, y0, S, wl, ph, side="top"):
    zz = np.linspace(xa, xb, 1001, endpoint=False)
    # excluding the last point - depending on the mesh/floating point errors,
    # sometimes this is actually in the next layer

    gg = g(zz) * ph  # generation rate * photon flux

    out = np.zeros_like(wl)
    sol_success = np.zeros_like(wl)  # keep track of whether solve_bvp is converging
    solutions = []

    for i in range(len(wl)):
        if np.all(
            gg[:, i] == 0
        ):  # no reason to solve anything if no generation at this wavelength
            out[i] = 0
            sol_success[i] = 0

        else:

            def A(x):
                return np.interp(x, zz, gg[:, i]) / D + y0 / L**2

            # generation and n0/p0 term in differential equation
            # (eq. 6.15 & 6.20 in Jenny Nelson, Physics of Solar Cells)

            def fun(x, y):
                # differential equation (eq. 6.15 & 6.20 in Jenny Nelson,
                # Physics of Solar Cells)
                # solve_bvp solves equation of form:
                # dy / dx = f(x, y, p), a <= x <= b
                # in this case y = [n or p, dn/dx or dp/dx]
                # y[0] = carrier concentration (n or p)
                # y[1] = carrier concentration gradient (dn/dx or dp/dx)
                out1 = y[1]  # by definition! dy/dx = dy/dx

                out2 = y[0] / L**2 - A(x)
                # actually solving the differential equation (6.15 & 6.20)

                return np.vstack((out1, out2))

            # boundary conditions for solve_bvp:
            if side == "top":

                def bc(ya, yb):
                    left = ya[1] - S / D * (ya[0] - y0)
                    # eq. 6.18 - b.c. at front of junction (surface recombination)

                    right = yb[0] - y0
                    # eq. 6.17 - b.c. edge of depletion region, top half of junction
                    # added - y0 (generally very small so makes almost no difference)
                    return np.array([left, right])

            else:

                def bc(ya, yb):
                    left = ya[0] - y0
                    print("left", left, y0)
                    # eq. 6.21 - b.c. edge of depletion region, bottom half of junction
                    # added - y0 (generally very small so makes almost no difference)

                    right = yb[1] + S / D * (yb[0] - y0)
                    # eq. 6.22 - b.c. at back of junction (surface recombination)
                    # changed sign! Current is going the other way
                    return np.array([left, right])

            guess = y0 * np.ones((2, zz.size))
            guess[1] = np.zeros_like(guess[0])

            solution = solve_bvp(fun, bc, zz, guess, max_nodes=2 * zz.shape[0])
            # increase max_nodes to avoid "too many mesh points" message

            sol_success[i] = solution.status

            if side == "top":
                out[i] = solution.y[1][-1]
                # current at edge of depletion region (top half of junction), eq. 6.33
            else:
                out[i] = solution.y[1][0]
                # current at edge of depletion region (bottom half of junction), eq 6.38

            solutions.append(solution)

        # give a warning f any of the solution statuses are not 0 using warnings.warn:
        if np.any(sol_success != 0):
            warnings.warn(
                "Depletion approximation (DA) EQE calculation: "
                "solve_bvp did not converge as expected for some wavelengths",
                RuntimeWarning,
            )

    return out, solutions


def get_J_sc_diffusion_green(xa, xb, g, D, L, y0, S, ph, side="top"):
    """Computes the derivative of the minority carrier concentration at the edge of the
    junction by approximating the convolution integral resulting from applying the
    Green's function method to the drift-diffusion equation.

    :param xa: Coordinate at the start the junction.
    :param xb: Coordinate at the end the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param D: Diffusion constant.
    :param L: Diffusion length.
    :param y0: Carrier equilibrium density.
    :param S: Surface recombination velocity.
    :param ph: Light spectrum.
    :param side: String to indicate the edge of interest. Either 'top' or 'bottom'.

    :return: The derivative of the minority carrier concentration at the edge of the
        junction.
    """

    xbL = (xb - xa) / L
    crvel = S / D * L
    ph_over_D = ph / D

    print("xbL", xbL)

    # if L too low in comparison to junction width, avoid nan's
    if xbL > 1.0e2:
        print("low L")
        if side == "top":
            fun = partial(_conv_exp_top, xa=xa, xb=xb, g=g, L=L, phoD=ph_over_D)
        else:
            fun = partial(_conv_exp_bottom, xa=xa, g=g, L=L, phoD=ph_over_D)
        cp = 1.0

    else:
        if side == "top":
            cp = -np.cosh(xbL) - crvel * np.sinh(xbL)

            fun = partial(
                _conv_green_top, xa=xa, xb=xb, g=g, L=L, phoD=ph_over_D, crvel=crvel
            )
        else:
            cp = np.cosh(xbL) + crvel * np.sinh(xbL)

            fun = partial(
                _conv_green_bottom, xb=xb, g=g, L=L, phoD=ph_over_D, crvel=-crvel
            )

    out, err, info = quad_vec(fun, xa, xb, epsrel=1.0e-5, full_output=True)
    # print(info)
    return out.squeeze() / cp


def _conv_exp_top(x, xa, xb, g, L, phoD):
    """Convolution of the carrier generation rate with the approximate Green's function
    kernel at point x. To be used with the numerical integration routine to compute the
    minority carrier derivative on the top edge. This kernel approximates the original
    one when the diffusion length is 2 orders of magnitude higher than the junction
    width by assuming that sinh(x) = cosh(x) = .5 * exp(x).

    :param x: Coordinate in the junction (variable to be integrated).
    :param xa: Coordinate at the start the junction.
    :param xb: Coordinate at the end the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param L: Diffusion length.
    :param phoD: Light spectrum divided by the diffusion constant D.
    """
    xc = (xa - x) / L
    xv = np.array(
        [
            xa + xb - x,
        ]
    )
    Pkern = -np.exp(xc)
    Gx = g(xv) * phoD
    return Pkern * Gx


def _conv_exp_bottom(x, xa, g, L, phoD):
    """Convolution of the carrier generation rate with the approximate Green's function
    kernel at point x. To be used with the numerical integration routine to compute the
    minority carrier derivative on the bottom edge. This kernel approximates the
    original one when the diffusion length is 2 orders of magnitude higher than
    the junction width by assuming that sinh(x) = cosh(x) = .5 * exp(x).

    :param x: Coordinate in the junction (variable to be integrated).
    :param xa: Coordinate at the start the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param L: Diffusion length.
    :param phoD: Light spectrum divided by the diffusion constant D.
    """
    xc = (xa - x) / L
    xv = np.array(
        [
            x,
        ]
    )
    Pkern = np.exp(xc)
    Gx = g(xv) * phoD
    return Pkern * Gx


def _conv_green_top(x, xa, xb, g, L, phoD, crvel):
    """Convolution of the carrier generation rate with the Green's function kernel at
    point x. To be used with the numerical integration routine to compute the minority
    carrier derivative on the top edge.

    :param x: Coordinate in the junction (variable to be integrated).
    :param xa: Coordinate at the start the junction.
    :param xb: Coordinate at the end the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param L: Diffusion length.
    :param phoD: Light spectrum divided by the diffusion constant D.
    :param crvel: Coefficient computed as S / D * L, with S the surface recombination
         velocity.
    """
    xc = (xb - x) / L
    xv = np.array(
        [
            xa + xb - x,
        ]
    )
    Pkern = np.cosh(xc) + crvel * np.sinh(xc)
    Gx = g(xv) * phoD
    return Pkern * Gx


def _conv_green_bottom(x, xb, g, L, phoD, crvel):
    """Convolution of the carrier generation rate with the Green's function kernel at
    point x. To be used with the numerical integration routine to compute the minority
    carrier derivative on the bottom edge.

    :param x: Coordinate in the junction (variable to be integrated).
    :param xb: Coordinate at the end the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param L: Diffusion length.
    :param phoD: Light spectrum divided by the diffusion constant D.
    :param crvel: Coefficient computed as S / D * L, with S the surface recombination
        velocity.
    """
    xc = (xb - x) / L
    xv = np.array(
        [
            x,
        ]
    )
    Pkern = np.cosh(xc) - crvel * np.sinh(xc)
    Gx = g(xv) * phoD
    return Pkern * Gx


bs = 1e20

Dn, Ln, Sn = 0.001263, 7e-06, 1000

alpha = GaAs.alpha(800e-9)
alpha = 1000


V = 0

T = 273

kT = kb * T

wp = 100e-9
R = 0

xp = 0

Na = 6e16 * 1e6
Nd = 3e17 * 1e6

ni = GaAs.ni
xi = 0

Dp, Lp, Sp = 0.001263, 7e-06, 400
p0 = ni**2 / Nd

es = GaAs.permittivity

Vbi = (kT / q) * np.log(Nd * Na / ni**2)

wn = (-xi + np.sqrt(xi**2 + 2.0 * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (
    1 + Nd / Na
)
wp = (-xi + np.sqrt(xi**2 + 2.0 * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (
    1 + Na / Nd
)

n0 = ni**2 / Na

wl = np.array([400, 500]) * 1e-9
alpha_vec = np.array([alpha, alpha * 1.1])
zz = np.linspace(0, width_top + width_bottom, 2001)
generation = interp1d(
    zz, (1 - R) * alpha_vec * np.exp(-alpha_vec * zz[:, None]), axis=0
)

all_sols_n = get_J_sc_diffusion_vs_WL(
    0, width_top - wp, generation, Dn, Ln, n0, Sn, wl, bs, "top"
)[1]
all_sols_p = get_J_sc_diffusion_vs_WL(
    width_top + wn,
    width_top + width_bottom,
    generation,
    Dp,
    Lp,
    p0,
    Sp,
    wl,
    bs,
    "bottom",
)[1]
n_sol_solcore = all_sols_n[0].y[0]
p_sol_solcore = all_sols_p[0].y[0]

xbb = width_top - wp - (width_top - wp) / 1001.0

all_sols_n_green = get_J_sc_diffusion_green(
    0, xbb, generation, Dn, Ln, n0, Sn, bs, "top"
)[0]

xbb = width_top + width_bottom - (width_top + width_bottom - width_top + wn) / 1001.0
all_sols_p_green = get_J_sc_diffusion_green(
    width_top + wn, xbb, generation, Dp, Lp, p0, Sp, bs, "bottom"
)[0]


def n_mathematica(x, L, alpha, xa, D, S, n0, A):
    aL = alpha * L
    a1 = -aL * np.exp(2 * alpha * (x + xa) + (2 * x / L)) + aL * np.exp(
        2 * alpha * (x + xa) + (2 * xa / L)
    )

    a2 = np.exp(2 * alpha * x + alpha * xa + (xa / L)) - np.exp(
        alpha * x + 2 * alpha * xa + (x / L)
    )

    a3 = np.exp(((aL + 1) * (2 * x + xa)) / L) - np.exp(((aL + 1) * (x + 2 * xa)) / L)

    b1 = np.exp(((aL + 1) * (2 * x + xa)) / L) - np.exp(((aL + 1) * (x + 2 * xa)) / L)

    b2 = np.exp(alpha * x + 2 * alpha * xa + (x / L)) - np.exp(
        2 * alpha * x + alpha * xa + (xa / L)
    )

    b3 = -np.exp(2 * alpha * (x + xa) + (2 * x / L)) + np.exp(
        2 * alpha * (x + xa) + (2 * xa / L)
    )

    numerator = (
        A
        * L**2
        * np.exp(-2 * alpha * (x + xa) - (x / L))
        * (D * (a1 + a2 + a3) + L * S * (b1 + b2 + b3))
    )
    denominator = (
        D
        * (alpha**2 * L**2 - 1)
        * (D * (np.exp(2 * xa / L) + 1) + L * S * (np.exp(2 * xa / L) - 1))
    )
    result = numerator / denominator + n0

    return result


def p_mathematica(x, L, alpha, xb, D, S, p0, A):
    aL = alpha * L
    denominator = (
        D
        * (-1 + alpha**2 * L**2)
        * (D * (1 + np.exp((2 * xb) / L)) + (-1 + np.exp((2 * xb) / L)) * L * S)
    )

    a1 = np.exp((2 * xb) / L + 2 * alpha * (x + xb)) - np.exp(
        (((1 + aL) * (x + 2 * xb)) / L)
    )

    a2 = np.exp((2 * x) / L + 2 * alpha * (x + xb)) - np.exp(
        alpha * x + x / L + 2 * alpha * xb
    )

    a3 = (
        alpha * np.exp(2 * alpha * x + alpha * xb + xb / L) * L
        - alpha * np.exp(((1 + aL) * (2 * x + xb)) / L) * L
    )

    b1 = np.exp(((1 + aL) * (2 * x + xb)) / L) - np.exp(((1 + aL) * (x + 2 * xb)) / L)

    b2 = +np.exp(alpha * x + x / L + 2 * alpha * xb) - np.exp(
        2 * alpha * x + alpha * xb + xb / L
    )

    b3 = np.exp((2 * xb) / L + 2 * alpha * (x + xb)) - np.exp(
        (2 * x) / L + 2 * alpha * (x + xb)
    )

    numerator = (
        A
        * np.exp(-(x / L) - 2 * alpha * (x + xb))
        * L**2
        * (D * (a1 + a2 + a3) + (b1 + b2 + b3) * L * S)
    )

    p = p0 + numerator / denominator
    return p


x_0 = np.linspace(0, width_top - wp, 1001, endpoint=False)

pref = (1 - R) * alpha * bs

n_ma = n_mathematica(x_0, Ln, alpha, (width_top - wp), Dn, Sn, n0, pref)

pref = (1 - R) * alpha * bs * np.exp(-alpha * (width_top + wn))
x_1 = np.linspace(0, width_bottom - wn, 1001, endpoint=False)

p_ma = p_mathematica(x_1, Lp, alpha, width_bottom - wn, Dp, Sp, p0, pref)

plt.figure()
plt.semilogy(x_0 * 1e9, n_ma)
plt.semilogy(x_0 * 1e9, n_sol_solcore, "--")
plt.tight_layout()
plt.show()

plt.figure()
plt.semilogy(x_1 * 1e9, p_ma)
plt.semilogy(x_1 * 1e9, p_sol_solcore, "--")
plt.tight_layout()
plt.show()

# plt.figure()
# plt.plot(x_1*1e9, p_ma/p_sol_solcore)
# plt.show()

print(n_ma == approx(n_sol_solcore, rel=0.02))
print(p_ma == approx(p_sol_solcore, rel=0.02))

plt.figure()
plt.semilogy(x_0 * 1e9, all_sols_n[0].y[1])
plt.show()

plt.figure()
plt.semilogy(x_1 * 1e9, all_sols_p[0].y[1])
plt.show()

print(all_sols_n[0].y[1][-1], all_sols_n_green)
print(all_sols_p[0].y[1][0], all_sols_p_green)
