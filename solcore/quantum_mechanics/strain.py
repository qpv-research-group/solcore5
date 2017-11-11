import numpy as np
import solcore
from solcore.state import State

"""
A comment on strain calculations. -- DJF July 2012

It is really easy to make mistakes when applying the strain code because the 
sign of the parameters is crucial and there appears to be no standard sign 
convention in the literature.

The procedure to align the heterointerfaces:

0) Align unstrained heterosructure using valence band offsets.

1) Shift conduction band using 'a_c' deformation potential (compressive/tensile = up/down)

2) Shift HH and LH bands using 'a_v' deformation potential (compressive/tensile = down/up)

3) Split the valence bands using the 'b' deformation potential (compressive/tensile = (HH)up, (LH)down/(HH)down, (LH)up)

"""

def strain_calculation_asserts(k, should_print=False):
    """Will check that the sign of parameters for the strain calcuations is correct."""
    if should_print: print(k)
    assert_sign = dict()

    # Sign for compressive strain
    assert_sign['epsilon'] = -1
    assert_sign['e_xx'] = -1
    assert_sign['e_yy'] = -1
    assert_sign['e_zz'] = 1
    assert_sign['Pe'] = 1
    assert_sign['Qe'] = -1
    assert_sign['delta_Ec'] = 1
    # assert_sign['deltaEhh'] = -1 # Sign of delta_Ehh not guaranteed
    # assert_sign['deltaElh'] = -1 # Sign of delta_Elh not guaranteed

    # Filp sign of these parameters for tensile strain
    if k.epsilon > 0:
        for key in assert_sign.keys():
            assert_sign[key] = assert_sign[key] * -1

    # Sign constant regardless of compressive or tensile strain
    assert_sign['ac'] = -1
    assert_sign['av'] = 1
    assert_sign['b'] = -1

    if should_print:
        print()
        print("Checking sign of strain parameters:")
    for key in assert_sign.keys():

        # If the material of the substrate form part of the layers, the sign makes no sense. We skip that.
        if k[key] == 0:
            continue

        if np.sign(k[key]) == assert_sign[key]:
            if should_print:
                print(key, "\tOK")
        else:
            if should_print:
                print(key, "\tPROBLEM")
            raise ValueError("The parameter `%s' used by the qwell module appears to have the wrong sign." % key)


def strain_calculation_parameters(substrate_material, layer_material, should_print=False, SO=False):
    """Will extract materials parameters and perform a bit of conditioning to the values.
    
    Returns a solcore State object with the following keys:
    
        -- ac, the conduction band deformation potential
        -- av, the valence band deformation potential
        -- b,  the valence band splitting deformation potential
        -- C11, element of the elastic stiffness tensor
        -- C12, element of the elastic stiffness tensor
        -- a0, the unstrained substrate lattice constant
        -- a, the unstrained layer lattice constant
        -- epsilon, in-plain strain
        -- epsilon_perp, out-plain strain
        -- e_xx, in-plain strain (e_xx = epsilon)
        -- e_yy, in-plain strain (e_yy = epsilon)
        -- e_zz, out-plain strain (e_zz = epsilon_perp)
        -- Tre, element of come matrix (Tre = e_xx + e_yy + e_zz)
        -- Pe, parameter use by Chuang
        -- Qe, parameter use by Chuang
        -- cb_hydrostatic_shift, CB moved by this amount
        -- vb_hydrostatic_shift, VB moved by this amount
        -- vb_shear_splitting, VB split by this amount (i.e. HH/LH separation)
        -- delta_Ec, final conduction band shift
        -- delta_Elh, final light hole band shift
        -- delta_Ehh, final heavy hole band shift
    
    Care has to be taken when calculating the energy shifts because different 
    sign conversion are used by different authors. Here we use the approach of
    S. L. Chuang, 'Physics of optoelectronic devices'. 
    
    Note that this requires that the 'a_v' deformation potential to be 
    positive, where as Vurgaftman defines this a negative!
    """

    sub = substrate_material
    mat = layer_material
    k = State()

    # deformation potentials
    k.av = abs(mat.a_v)  # make sure that av is positive for this calculation
    k.ac = mat.a_c  # Vurgaftman uses the convention that this is negative
    k.b = mat.b

    # Matrix elements from elastic stiffness tensor
    k.C11 = mat.c11
    k.C12 = mat.c12
    if should_print: print(sub, mat)
    # Strain fractions
    k.a0 = sub.lattice_constant
    k.a = mat.lattice_constant
    k.epsilon = (k.a0 - k.a) / k.a  # in-plain strain
    k.epsilon_perp = -2 * k.C12 / k.C11 * k.epsilon  # out-plain
    k.e_xx = k.epsilon
    k.e_yy = k.epsilon
    k.e_zz = k.epsilon_perp
    k.Tre = (k.e_xx + k.e_yy + k.e_zz)
    k.Pe = -k.av * k.Tre
    k.Qe = -k.b / 2 * (k.e_xx + k.e_yy - 2 * k.e_zz)

    k.cb_hydrostatic_shift = k.ac * k.Tre
    k.vb_hydrostatic_shift = k.av * k.Tre
    k.vb_shear_splitting = 2 * k.b * (1 + 2 * k.C12 / k.C11) * k.epsilon

    # Shifts and splittings
    k.delta_Ec = k.ac * k.Tre

    if should_print: print(k.ac, k.Tre)

    k.delta_Ehh = -k.Pe - k.Qe
    k.delta_Elh = -k.Pe + k.Qe
    k.delta_Eso = 0.0

    if SO:
        k.delta = mat.spin_orbit_splitting
        shift = k.delta ** 2 + 2 * k.delta * k.Qe + 9 * k.Qe ** 2
        k.delta_Elh = -k.Pe + 0.5 * (k.Qe - k.delta + np.sqrt(shift))
        k.delta_Eso = -k.Pe + 0.5 * (k.Qe - k.delta - np.sqrt(shift))

    strain_calculation_asserts(k, should_print=should_print)

    if should_print:
        print()
        print("Lattice:")
        print("a0", k.a0)
        print("a", k.a)
        print()
        print("Deformation potentials:")
        print("ac = ", solcore.asUnit(k.ac, 'eV'))
        print("av = ", solcore.asUnit(k.av, 'eV'))
        print("ac - av = ", solcore.asUnit(k.ac - k.av, 'eV'))
        print("b = ", solcore.asUnit(k.b, 'eV'))
        print()
        print("Matrix elements from elastic stiffness tensor:")
        print("C_11 = ", solcore.asUnit(k.C11, "GPa"))
        print("C_12 = ", solcore.asUnit(k.C12, "GPa"))
        print()
        print("Strain fractions:")
        print("e_xx = e_yy = epsilon = ", k.epsilon)
        print("e_zz = epsilon_perp = ", k.epsilon_perp)
        print("e_xx + e_yy + e_zz = Tre = ", k.Tre)
        print()
        print("Shifts and splittings:")
        print("Pe = -av * Tre = ", solcore.asUnit(k.Pe, 'eV'))
        print("Qe = -b/2*(e_xx + e_yy - 2*e_zz) = ", solcore.asUnit(k.Qe, 'eV'))
        print("dEc = ac * Tre = ", solcore.asUnit(k.delta_Ec, 'eV'))
        print("dEhh = av * Tre + b[1 + 2*C_11/C_12]*epsilon = -Pe - Qe = ", solcore.asUnit(k.delta_Ehh, 'eV'))
        print("dElh = av * Tre - b[1 + 2*C_11/C_12]*epsilon = -Pe + Qe = ", solcore.asUnit(k.delta_Elh, 'eV'))
        print()

    return k
