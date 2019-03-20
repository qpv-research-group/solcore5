import sys
from random import random

import numpy as np

from solcore.constants import *
from solcore.smooth import smooth
from solcore.structure import *
from solcore.quantum_mechanics.heterostructure_alignment import VBO_align
from solcore.quantum_mechanics.strain import strain_calculation_parameters

from solcore.quantum_mechanics.kp_QW import kp6x6, kp4x4
from solcore.quantum_mechanics.kp_bulk import kp8x8_bulk

m0 = electron_mass


def assemble_qw_structure(repeats, well, bulk_l_top, bulk_l_bottom, barrier, well_interlayer=None,
                          structure_label="QW Structure", shift_wells=0):
    half_barrier = Layer(barrier.width / 2, barrier.material)
    qw_structure = Structure()
    # qw_structure.append(bulk, layer_label="bulk")
    qw_structure.append_multiple([bulk_l_top], layer_labels=["bulk"])
    qw_structure.append_multiple([half_barrier], layer_labels=["half barrier"])
    if well_interlayer:
        qw_structure.append_multiple([well_interlayer, well, well_interlayer, barrier], layer_labels=["interlayer",
                                                                                                      "well",
                                                                                                      "interlayer",
                                                                                                      "barrier"],
                                     repeats=repeats - 1)
        qw_structure.append_multiple([well_interlayer, well, well_interlayer, half_barrier, bulk_l_bottom],
                                     layer_labels=[
                                         "interlayer", "well", "interlayer", "half barrier", "bulk"])
    else:
        qw_structure.append_multiple([well, barrier], layer_labels=["well", "barrier"], repeats=repeats - 1)
        qw_structure.append_multiple([well, half_barrier, bulk_l_bottom], layer_labels=["well", "half barrier", "bulk"])

    return qw_structure


def vary_well_widths(structure, fraction=0.5, region_sought="well"):
    """Returns an array of indices which correspond to the region that are quantum wells in the structure."""
    depth = 0

    for i, (layer, label) in enumerate(zip(structure, structure.labels)):
        if label != region_sought:
            continue

        layer.width *= (1 + (2 * random() - 1) * fraction)
        print(layer.width)
    return structure


def locate_regions(x, structure, region_sought="well"):
    """Returns an array of indices which correspond to the region that are quantum wells in the structure."""
    depth = 0

    well_structure_indices = set([item for item in range(
        len(structure.labels)) if structure.labels[
                                      item] == region_sought])  # <---- Uses new labels in the structure object
    well_z_indices = np.array([], dtype=np.dtype(int))
    for i, layer in enumerate(structure):
        if well_structure_indices.issuperset([i]):
            well_z_indices = np.hstack((well_z_indices, np.where((depth <= x) * (x <= depth + layer.width) == True)[0]))
        depth += layer.width
    return well_z_indices


def well_regions(x, structure):
    return locate_regions(x, structure, "well")


def text_render(structure, resolution=100):
    x = np.linspace(0, structure.width(), resolution)
    bulk = locate_regions(x, structure, "bulk")
    barrier = set(locate_regions(x, structure, "barrier")) | set(locate_regions(x, structure, "half barrier"))
    interlayer = locate_regions(x, structure, "interlayer")
    well = locate_regions(x, structure, "well")

    chars = ["―" if i in bulk else "_" if i in well else "-" if i in interlayer else "‾" if i in barrier else "?" for i
             in range(resolution)]
    return ("".join(chars))


def structure_to_potentials(structure, step_size=None, minimum_step_size=0, smallest_feature_steps=20,
                            blur=False, blurmode="even", return_qw_boolean_for_layer=False, Efield=0,
                            mode='kp4x4'):
    """ Discretizes the structure as a function of the position, providing the potential for the electrons, HH, LH and
    SO bands as well as the Luttinger parameters and the effective masses.

    The result depend on the chosen mode:
        - kp4x4, calculates the bands and effective masses at k=0 coupling the HH and LH bands
        - kp6x6, calculates the bands and effective masses at k=0 coupling the HH, LH and SO bands
        - kp8x8_bulk, calculates the four bands and fit the effective masses with a parabola around k=0
        - strain, just shifts the bands according to the strain
        - relaxed, do nothing and things are calculated as if we had the bulk, unstrained materials

    Notice that although the kp8x8_bulk couples the 4 (8) bands, the effective mass is the result of a fitting around
    k=0. Therefore, it is not valid to solve 1D problems, but just a bulk-like problem under a parabolic aproximation

    Additionlly, the band structure can be blured - to simulate the intermixing of materials - as well as being under
    an electric field.
    """

    # We calculate the mesh based on the minimum feature of the structure and the number of points desired for that
    # feature

    available_modes = ['kp6x6', 'kp4x4', 'kp8x8_bulk', 'strain', 'relaxed']
    assert mode in available_modes, 'ERROR: Calculation mode must be {}'.format(available_modes)


    widths = [layer.width for layer in structure]
    total_width = sum(widths)
    if step_size is None:
        smallest_feature_width = min(widths)
        step_size = int(smallest_feature_steps / smallest_feature_width * total_width)
        step_size = max(step_size, minimum_step_size)

    x = np.linspace(0, total_width, step_size)
    Ve = np.zeros(len(x))  # Electron potential
    Vhh = np.zeros(len(x))  # Heavy hole    "
    Vlh = np.zeros(len(x))  # Light hole    "
    Vso = np.zeros(len(x))  # SO band    "
    me = np.zeros(len(x))  # Electron effective mass
    mhh_p = np.zeros(len(x))  # HH paralell effective mass
    mhh_t = np.zeros(len(x))  # HH transverse effective mass
    mlh_p = np.zeros(len(x))  # LH paralell effective mass
    mlh_t = np.zeros(len(x))  # LH transverse effective mass
    g1 = np.zeros(len(x))  # Gamma 1
    g2 = np.zeros(len(x))  # Gamma 2
    g3 = np.zeros(len(x))  # Gamma 3
    isWell = np.zeros(len(x), dtype=bool)  # indicate where we have the QW    "

    # Spatially varying effective mass
    depth = 0

    # At each position of the structure, we have to calculate the effective mass of electrons and holes, and the band
    # edges.
    for i, layer in enumerate(structure):

        positions = (depth <= x) * (x <= depth + layer.width)

        me[positions] = layer.material.get("eff_mass_electron") * electron_mass
        g1[positions] = layer.material.get("gamma1")
        g2[positions] = layer.material.get("gamma2")
        g3[positions] = layer.material.get("gamma3")
        mhh_p[positions] = 1. / (g1[positions] - 2 * g2[positions]) * m0
        mhh_t[positions] = 1. / (g1[positions] + 2 * g2[positions]) * m0
        mlh_p[positions] = 1. / (g1[positions] + g2[positions]) * m0
        mlh_t[positions] = 1. / (g1[positions] - g2[positions]) * m0

        if return_qw_boolean_for_layer:
            isWell[positions] = (layer == return_qw_boolean_for_layer)

        try:
            Ve[positions] = layer.Ec
            Vhh[positions] = layer.Ev
            Vlh[positions] = layer.Ev
            Vso[positions] = layer.Ev - layer.material.spin_orbit_splitting

        except KeyError:
            print("There was a problem converting the structure to a potential because some information was missing "
                  "from the structure state.")
            print("Try calling the function 'align_heterostructure_using_Vurgaftman' before calling "
                  "'structure_to_potentials'.")
            sys.exit()

        except Exception as inst:
            print(inst)
            raise

        substrate = structure.substrate

        # Now we modify the band profile and effective masses according to the chosen mode:
        # - kp4x4, calculates the bands and effective masses at k=0 coupling the HH and LH bands
        # - kp6x6, calculates the bands and effective masses at k=0 coupling the HH, LH and SO bands
        # - kp8x8_bulk, calculates the four bands and fit the effective masses with a parabola around k=0
        # - strain, just shifts the bands according to the strain
        # - relaxed, do nothing and things are calculated as if we had the bulk, unstrained materials
        if mode == 'kp6x6':
            # The band edges are shifted and the effective mass of the LHs also changes.
            c, hh, lh, so, mc_kp, mhh_p_kp, mlh_p_kp, mso_p_kp, mhh_t_kp, mlh_t_kp, mso_t_kp = kp6x6(layer.material,
                                                                                                         substrate)

            mlh_p[positions] = mlh_p_kp
            mlh_t[positions] = mlh_t_kp

            Ve[positions] = c
            Vhh[positions] = hh
            Vlh[positions] = lh
            Vso[positions] = so

        elif mode == 'kp4x4':
            # The band edges are shifted but the effective masses around k=0 are not affected.
            c, hh, lh, mc_kp, mhh_p_kp, mlh_p_kp, mhh_t_kp, mlh_t_kp = kp4x4(layer.material, substrate)

            Ve[positions] = c
            Vhh[positions] = hh
            Vlh[positions] = lh

        elif mode == 'kp8x8_bulk':
            # The band edges and the effective mass of Electrons, HH, LH and the SO holes changes.
            c, hh, lh, so, mc, mhh, mlh, mso = kp8x8_bulk(layer.material, substrate)

            me[positions] = mc
            mhh_p[positions] = mhh
            mlh_p[positions] = mlh
            mhh_t[positions] = mhh
            mlh_t[positions] = mlh

            Ve[positions] = c
            Vhh[positions] = hh
            Vlh[positions] = lh
            Vso[positions] = so

        elif mode == 'strain':
            # The bands are just shifted according to strain
            strain_parameters = strain_calculation_parameters(substrate, layer.material, SO=True)
            Ve[positions] = Ve[positions] + strain_parameters.delta_Ec
            Vhh[positions] = Vhh[positions] + strain_parameters.delta_Ehh
            Vlh[positions] = Vlh[positions] + strain_parameters.delta_Elh

        else:
            pass

        depth += layer.width

    # We blur the structure if required
    if blur:
        Ve = smooth(Ve, blur / total_width * step_size, blurmode=blurmode)
        Vhh = smooth(Vhh, blur / total_width * step_size, blurmode=blurmode)
        Vlh = smooth(Vlh, blur / total_width * step_size, blurmode=blurmode)
        Vso = smooth(Vso, blur / total_width * step_size, blurmode=blurmode)

        me = smooth(me, blur / total_width * step_size, blurmode=blurmode)
        g1 = smooth(g1, blur / total_width * step_size, blurmode=blurmode)
        g2 = smooth(g2, blur / total_width * step_size, blurmode=blurmode)
        g3 = smooth(g3, blur / total_width * step_size, blurmode=blurmode)

        mhh_p = smooth(mhh_p, blur / total_width * step_size, blurmode=blurmode)
        mhh_t = smooth(mhh_t, blur / total_width * step_size, blurmode=blurmode)
        mlh_p = smooth(mlh_p, blur / total_width * step_size, blurmode=blurmode)
        mlh_t = smooth(mlh_t, blur / total_width * step_size, blurmode=blurmode)

    # We tilt the structure if there is an electric field
    if Efield != 0:
        Ve = Ve + q * x * Efield
        Vhh = Vhh + q * x * Efield
        Vlh = Vlh + q * x * Efield
        Vso = Vso + q * x * Efield

    if mode in ['kp4x4', 'strain', 'relaxed']:
        Vso = Vso * 0


    if mode in ['kp4x4', 'kp6x6']:
        return {"x": x,
                "Ve": Ve,
                "me": me,
                "Vhh": Vhh,
                "g1": g1,
                "Vlh": Vlh,
                "g2": g2,
                "Vso": Vso,
                "g3": g3,
                "isWell": isWell,
                "mhh_p": mhh_p,
                "mhh_t": mhh_t,
                "mlh_p": mlh_p,
                "mlh_t": mlh_t}
    else:
        return {"x": x,
                "Ve": Ve,
                "me": me,
                "Vhh": Vhh,
                "Vlh": Vlh,
                "Vso": Vso,
                "isWell": isWell,
                "mhh": mhh_p,
                "mlh": mlh_p}


if __name__ == "__main__":
    from solcore import si, material
    from solcore.structure import Layer
    import matplotlib.pyplot as plt

    bulk = material("GaAs")(T=293)
    QW = material("InGaAs")(T=293, In=0.147)
    barrier = material("GaAsP")(T=293, P=0.1)

    bulk.strained = False
    QW.strained = True
    barrier.strained = True

    top_layer = Layer(width=si("50nm"), material=bulk)
    well_layer = Layer(width=si("7.2nm"), material=QW)
    barrier_layer = Layer(width=si("29nm"), material=barrier)
    bottom_layer = top_layer

    print("Here is a ASC_examples structure with no materials:")

    test_structure = assemble_qw_structure(
        repeats=3,
        well=well_layer,
        bulk_l_top=top_layer,
        bulk_l_bottom=bottom_layer,
        barrier=barrier_layer,
    )

    test_structure.substrate = bulk

    VBO_align(test_structure)
    result = structure_to_potentials(test_structure, mode='kp8x8_bulk')
    print(test_structure)

    # plt.plot(result['x']*1e9, result['Ve']/q, result['x']*1e9, result['Vhh']/q, result['x']*1e9, result['Vlh']/q, result['x']*1e9, result['Vso']/q)
    #    plt.ylim(-1.5, 1.0)
    # plt.plot(result['x']*1e9, result['me']/m0, result['x']*1e9, result['mhh_p']/m0, result['x']*1e9, result['mlh_p']/m0)
    # plt.plot(result['x'] * 1e9, result['me'] / m0, result['x'] * 1e9, result['mhh_t'] / m0, result['x'] * 1e9,
    #          result['mlh_t'] / m0)

    plt.show()
