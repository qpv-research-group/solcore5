from scipy.optimize import bisect
from numpy import sqrt, log, pi, linspace, array, ones
from solcore import convert, asUnit
from solcore import ParameterSystem
from solcore.science_tracker import science_reference


# Zero stress criterion from:
# Ekins-Daukes et al. Strain-balanced Criteria for Multiple Quantum Well Structures and Its
# Signature in X-ray Rocking Curves. Crystal growth & design (2002) vol. 2 (4) pp. 287-292

# Critical Thickness from:
# Matthews, J., & Blakeslee, A. (1974). Defects in epitaxial multilayers: I. Misfit dislocations.
# Journal of Crystal Growth, 27, 118-125.


def optimise(well_material, well_thickness, well_fraction, barrier_material, barrier_thickness, barrier_fraction,
        lattice_material, lattice_material_fraction=0, parameter_to_optimise="well_thickness", T=300):
    """ This can't handle sol-core materials yet. Needs layers."""

    options = ["barrier_thickness", "well_thickness", "barrier_fraction", "well_fraction"]
    if parameter_to_optimise not in options:
        print('ERROR: Wring option in strain_balancing.optimise. '
              'The only valid options are: {}, {}, {}, {}'.format(*options))

    science_reference("strain-balancing",
                      "Ekins-Daukes et al. Strain-balanced Criteria for Multiple Quantum Well Structures and Its "
                      "Signature in X-ray Rocking Curves. Crystal growth & design (2002) vol. 2 (4) pp. 287-292")

    get_vurgaftman_parameter = ParameterSystem().get_parameter
    t1 = well_thickness
    t2 = barrier_thickness

    x1 = well_fraction
    x2 = barrier_fraction

    C11_well = get_vurgaftman_parameter(well_material, "c11", x=x1, T=T)
    C12_well = get_vurgaftman_parameter(well_material, "c12", x=x1, T=T)
    A1 = C11_well + C12_well - 2 * C12_well ** 2 / C11_well
    a1 = get_vurgaftman_parameter(well_material, "lattice_constant", x=x1, T=T)

    C11_barrier = get_vurgaftman_parameter(barrier_material, "c11", x=x2, T=T)
    C12_barrier = get_vurgaftman_parameter(barrier_material, "c12", x=x2, T=T)
    A2 = C11_barrier + C12_barrier - 2 * C12_barrier ** 2 / C11_barrier
    a2 = get_vurgaftman_parameter(barrier_material, "lattice_constant", x=x2, T=T)

    a0 = get_vurgaftman_parameter(lattice_material, "lattice_constant", x=lattice_material_fraction, T=T)

    a0 = convert(a0, "angstrom", "m")
    a1 = convert(a1, "angstrom", "m")
    a2 = convert(a2, "angstrom", "m")
    if parameter_to_optimise == "barrier_thickness":
        t2 = - (A1 * t1 * a2 ** 2 * (a0 - a1)) / (A2 * a1 ** 2 * (a0 - a2))
    if parameter_to_optimise == "well_thickness":
        t1 = - (A2 * t2 * a1 ** 2 * (a0 - a2)) / (A1 * a2 ** 2 * (a0 - a1))

    if parameter_to_optimise == "well_fraction":
        def func(x):
            C11_well = get_vurgaftman_parameter(well_material, "c11", x=x, T=T)
            C12_well = get_vurgaftman_parameter(well_material, "c12", x=x, T=T)
            A1 = C11_well + C12_well - 2 * C12_well ** 2 / C11_well
            a1 = get_vurgaftman_parameter(well_material, "lattice_constant", x=x, T=T)
            a1 = convert(a1, "angstrom", "m")

            misfit = - (A1 * t1 * a2 ** 2 * (a0 - a1)) / (A2 * a1 ** 2 * (a0 - a2)) - t2
            return misfit

        print(func(linspace(0, 1, 10)))

        x1 = bisect(func, 0, 1)

    if parameter_to_optimise == "barrier_fraction":
        def func(x):
            C11_barrier = get_vurgaftman_parameter(barrier_material, "c11", x=x, T=T)
            C12_barrier = get_vurgaftman_parameter(barrier_material, "c12", x=x, T=T)
            A2 = C11_barrier + C12_barrier - 2 * C12_barrier ** 2 / C11_barrier
            a2 = get_vurgaftman_parameter(barrier_material, "lattice_constant", x=x, T=T)
            a2 = convert(a2, "angstrom", "m")

            if a0 == a2 and x == 0:
                a2 = get_vurgaftman_parameter(barrier_material, "lattice_constant", x=0.001, T=T)
                a2 = convert(a2, "angstrom", "m")

            misfit = - (A1 * t1 * a2 ** 2 * (a0 - a1)) / (A2 * a1 ** 2 * (a0 - a2)) - t2
            return misfit

        x2 = bisect(func, 0, 1)

    return {
        "well_material": well_material,
        "well_thickness": t1,
        "well_fraction": x1,
        "barrier_material": barrier_material,
        "barrier_thickness": t2,
        "barrier_fraction": x2,
        "lattice_material": lattice_material,
        "lattice_fraction": lattice_material_fraction
    }


def critical_thickness(layer_material="GaInAs", lattice_material="GaAs", layer_fraction=0, lattice_fraction=None, T=300,
                       during_growth=True, final_unit="m", bowed_material="In"):
    science_reference("critical thickness",
                      "Matthews, J., & Blakeslee, A. (1974). Defects in epitaxial multilayers: I. Misfit dislocations. "
                      "Journal of Crystal Growth, 27, 118-125.")
    get_vurgaftman_parameter = ParameterSystem().get_parameter

    other_params = {
        "T": T,
        bowed_material: layer_fraction,
    }
    c11 = get_vurgaftman_parameter(layer_material, "c11", **other_params)
    c12 = get_vurgaftman_parameter(layer_material, "c12", **other_params)
    a = get_vurgaftman_parameter(layer_material, "lattice_constant", **other_params)
    a_lattice = get_vurgaftman_parameter(lattice_material, "lattice_constant", **other_params)

    b = array(a / sqrt(2))
    v = array(abs(c12 / (c11 + c12)))
    f = array(abs(a / a_lattice - 1))

    longest = max([each.size for each in [v, b, f]])

    if longest != 1:
        if b.size != longest:
            assert b.size == 1, "array arguments do not match length"
            b = ones(longest) * b
        if v.size != longest:
            assert v.size == 1, "array arguments do not match length"
            v = ones(longest) * v
        if f.size != longest:
            assert f.size == 1, "array arguments do not match length"
            f = ones(longest) * f
    else:
        v, b, f = [v], [b], [f]

    result = []
    for vi, fi, bi in zip(v, f, b):
        roots_at_soln = lambda hc_over_b: (hc_over_b) / (log(hc_over_b) + 1) - (1 - vi / 4) / (
            fi * pi * (1 + vi))  # matthews & blakeslee
        try:
            hc = bisect(roots_at_soln, 1e-10, 1e10) * bi
        except ValueError:
            raise ValueError("Critical thickness infinite?")

        result.append(convert(hc / 4, "Ang", final_unit)) if during_growth else result.append(
            convert(hc, "Ang", final_unit))
    return array(result)


if __name__ == "__main__":
    T = linspace(300, 1000, 10)
    from solcore.graphing import *

    result = critical_thickness(layer_material="GaInAs", lattice_material="GaAs", layer_fraction=0.24, T=T,
                                final_unit="nm")
    g = Graph(
        GraphData(T, asUnit(result, 'nm'), label="Matthews-Blakeslee $\\frac{h_c}{4}$"),  # , "-", 1,"red"),
        # yscale="log", 
        # xlim=(0.03,0.53), 
        # ylim=(1,1000), 
        xlabel="Temperature (K)",
        ylabel="Critical Thickness (nm)",
        # grid=True,
        title="InGaAs critical thickness"
        # labels = ([(x_at_024,y_at_024,"x=0.24, $\\frac{h_c}{4}$=%.2fnm"%y_at_024)])
    )
    g.draw()
    T = 300
    x = linspace(0, 1, 100)[1:]
    result = critical_thickness(layer_material="GaInAs", lattice_material="GaAs", layer_fraction=x, T=T,
                                final_unit="nm")
    print(result)
    g = Graph(
        GraphData(x, asUnit(result, 'nm'), label="Matthews-Blakeslee $\\frac{h_c}{4}$"),
        # (f[0], f[1], "Previous calculation", "-g")
        # yscale="log",
        # size="4 square",
        # figsize=(4, 8),
        xlim=(0.03, 0.53),
        ylim=(0, 10),
        xlabel="Indium Fraction",
        ylabel="Critical Thickness (nm)",
        # grid=True,
        title="InGaAs critical thickness"
        # labels = ([(x_at_024,y_at_024,"x=0.24, $\\frac{h_c}{4}$=%.2fnm"%y_at_024)])
    )
    g.draw()
