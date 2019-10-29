from pytest import approx, fixture
import random
from solcore import material
from solcore.poisson_drift_diffusion.DeviceStructure import CreateDeviceStructure
from solcore.structure import Junction, Layer, SolcoreMaterialToStr, Structure, ToSolcoreMaterial, TunnelJunction
from solcore.structure import InLineComposition, ToLayer, ToStructure

# Materials
t = 300
qw_mat_material_name = 'InGaAs'
qw_mat_in = 0.2
qw_mat_structure = {'material': qw_mat_material_name, 'element': 'In', 'fraction': qw_mat_in}
b_mat_material_name = 'GaAsP'
b_mat_p = 0.1
b_mat_structure = {'material': b_mat_material_name, 'element': 'P', 'fraction': b_mat_p}
i_gaas_material_name = 'GaAs'
i_gaas_structure = {'material': i_gaas_material_name}

# Layers
available_roles = ['Barrier', 'Base', 'BSF', 'Emitter', 'Intrinsic', 'Window']
wkt_box = 'POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))'
width1 = random.uniform(1e-9, 1e-8)
role1 = random.choice(available_roles)
width2 = random.uniform(1e-9, 1e-8)
role2 = random.choice(available_roles)
width3 = random.uniform(1e-9, 1e-8)
role3 = random.choice(available_roles)

@fixture
def qw_mat():
    return material(qw_mat_material_name)(T=t, In=qw_mat_in)


@fixture
def b_mat():
    return material(b_mat_material_name)(T=t, P=b_mat_p)


@fixture
def i_gaas():
    return material(i_gaas_material_name)(T=t)


@fixture
def layer1(qw_mat):
    return Layer(width1, qw_mat, role1, wkt_box, new_property='new_property')


@fixture
def layer2(b_mat):
    return Layer(width2, b_mat, role2, wkt_box)


@fixture
def layer3(i_gaas):
    return Layer(width3, i_gaas, role3, wkt_box)


@fixture
def device(layer1, layer2, layer3):
    return CreateDeviceStructure(layers=[layer1, layer2, layer3])


def test_layer_and_junction(qw_mat, b_mat, i_gaas, layer1, layer2, layer3):
    assert layer1.width == width1
    assert layer1.role == role1
    assert layer1.material == qw_mat
    assert layer1.geometry == wkt_box
    assert layer1.__dict__['new_property'] == 'new_property'

    assert layer2.width == width2
    assert layer2.role == role2
    assert layer2.material == b_mat
    assert layer2.geometry == wkt_box

    assert layer3.width == width3
    assert layer3.role == role3
    assert layer3.material == i_gaas
    assert layer3.geometry == wkt_box

    sn = random.uniform(1e6, 1e7)
    sp = random.uniform(1e6, 1e7)
    junction1 = Junction([layer1, layer2, layer3], sn=sn, sp=sp, T=t, kind='PDD')

    assert junction1.__len__() == 3
    assert junction1.__dict__ == {'sn': sn, 'sp': sp, 'T': t, 'kind': 'PDD'}
    assert junction1[0] == layer1
    assert junction1[1] == layer2
    assert junction1[2] == layer3

    tunnel1 = TunnelJunction([layer1, layer2, layer3], sn=sn, sp=sp, T=t, kind='PDD')

    assert tunnel1.__len__() == 3
    assert tunnel1.__dict__ == {'sn': sn, 'sp': sp, 'T': t, 'kind': 'PDD', 'R': 1e-16, 'pn': True}
    assert tunnel1[0] == layer1
    assert tunnel1[1] == layer2
    assert tunnel1[2] == layer3


def test_structure(b_mat, layer1, layer2, layer3):
    structure1 = Structure([layer1, layer2], substrate=b_mat, T=t)

    assert structure1.__len__() == 2
    assert structure1.__dict__ == {'substrate': b_mat, 'T': t, 'labels': [None, None]}
    assert structure1.width() == approx(width1 + width2)
    assert structure1[0] == layer1
    assert structure1[1] == layer2

    structure2 = Structure([], substrate=b_mat, T=t)
    structure2.append(layer1, layer_label='layer1')

    assert structure2.__len__() == 1
    assert structure2.__dict__ == {'substrate': b_mat, 'T': t, 'labels': ['layer1']}
    assert structure2[0] == layer1

    structure2.append_multiple([layer2, layer3], layer_labels=['layer2', 'layer3'])

    assert structure2.__len__() == 3
    assert structure2.__dict__ == {'substrate': b_mat, 'T': t, 'labels': ['layer1', 'layer2', 'layer3']}
    assert structure2[0] == layer1
    assert structure2[1] == layer2
    assert structure2[2] == layer3
    assert structure2.width() == approx(width1 + width2 + width3)
    assert structure2.relative_widths()['layer1'] == approx(width1 / structure2.width())
    assert structure2.relative_widths()['layer2'] == approx(width2 / structure2.width())
    assert structure2.relative_widths()['layer3'] == approx(width3 / structure2.width())

    structure3 = Structure([])
    structure3.append(layer1, layer_label='layer1', repeats=2)

    assert structure3.__len__() == 2
    assert structure3.__dict__ == {'labels': ['layer1', 'layer1']}
    assert structure3[0] == layer1
    assert structure3[1] == layer1
    assert structure3.width() == approx(width1 * 2)
    # Below currently fails due to a bug when specifying labels and repeating more than once
    # assert structure3.relative_widths()['layer1'] == width1 / structure3.width()

    structure4 = Structure([])
    structure4.append_multiple([layer1, layer2], layer_labels=['layer1', 'layer2'], repeats=2)

    assert structure4.__len__() == 4
    assert structure4.__dict__ == {'labels': ['layer1', 'layer2', 'layer1', 'layer2']}
    assert structure4[0] == layer1
    assert structure4[1] == layer2
    assert structure4[2] == layer1
    assert structure4[3] == layer2
    assert structure4.width() == approx(width1 * 2 + width2 * 2)
    # Below currently fails due to a bug when specifying labels and repeating more than once
    # assert structure4.relative_widths()['layer1'] == width1 / structure4.width()
    # assert structure4.relative_widths()['layer2'] == width2 / structure4.width()


def test_material_to_str(qw_mat, b_mat, i_gaas):
    assert SolcoreMaterialToStr(qw_mat) == qw_mat_structure
    assert SolcoreMaterialToStr(b_mat) == b_mat_structure
    assert SolcoreMaterialToStr(i_gaas) == i_gaas_structure


def test_to_material(qw_mat, b_mat, i_gaas):
    qw_material = ToSolcoreMaterial(qw_mat_structure, t, True)
    assert qw_material.__class__.__name__ == qw_mat_material_name
    assert qw_material.__dict__ == qw_mat.__dict__

    b_material = ToSolcoreMaterial(b_mat_structure, t, True)
    assert b_material.__class__.__name__ == b_mat_material_name
    assert b_material.__dict__ == b_mat.__dict__

    i_gaas_material = ToSolcoreMaterial(i_gaas_structure, t, True)
    assert i_gaas_material.__class__.__name__ == i_gaas_material_name
    assert i_gaas_material.__dict__ == i_gaas.__dict__


def test_inline_composition(device):
    assert InLineComposition(device['layers'][0]) == 'In0.2GaAs'
    assert InLineComposition(device['layers'][1]) == 'GaAsP0.1'
    assert InLineComposition(device['layers'][2]) == 'GaAs'


def test_to_layer(qw_mat, device):
    device_layer = device['layers'][0]
    composition = device_layer['properties']['composition']
    layer = ToLayer(width1, ToSolcoreMaterial(composition, t), role1)

    assert layer.width == width1
    assert layer.role == role1
    assert layer.material.__str__() == qw_mat.__str__()


def test_to_structure(device):
    structure = ToStructure(device)

    assert structure.__len__() == 3
    assert structure.width() == approx(width1 + width2 + width3)
