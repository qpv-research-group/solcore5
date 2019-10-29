from pytest import fixture
import random
from solcore import material
from solcore.poisson_drift_diffusion.DeviceStructure import CreateDeviceStructure
from solcore.structure import Junction, Layer, SolcoreMaterialToStr, Structure, ToSolcoreMaterial, TunnelJunction
from solcore.structure import InLineComposition, ToLayer, ToStructure

# Materials
T = 300
QWmat_material_name = 'InGaAs'
QWmat_in = 0.2
QWmat_structure = {'material': QWmat_material_name, 'element': 'In', 'fraction': QWmat_in}
Bmat_material_name = 'GaAsP'
Bmat_p = 0.1
Bmat_structure = {'material': Bmat_material_name, 'element': 'P', 'fraction': Bmat_p}
i_GaAs_material_name = 'GaAs'
i_GaAs_structure = {'material': i_GaAs_material_name}

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
def QWmat():
    return material(QWmat_material_name)(T=T, In=QWmat_in)


@fixture
def Bmat():
    return material(Bmat_material_name)(T=T, P=Bmat_p)


@fixture
def i_GaAs():
    return material(i_GaAs_material_name)(T=T)


@fixture
def layer1(QWmat):
    return Layer(width1, QWmat, role1, wkt_box, new_property='new_property')


@fixture
def layer2(Bmat):
    return Layer(width2, Bmat, role2, wkt_box)


@fixture
def layer3(i_GaAs):
    return Layer(width3, i_GaAs, role3, wkt_box)


@fixture
def device(layer1, layer2, layer3):
    return CreateDeviceStructure(layers=[layer1, layer2, layer3])


def test_layer_and_junction(QWmat, layer1, Bmat, layer2, i_GaAs, layer3):
    assert layer1.width == width1
    assert layer1.role == role1
    assert layer1.material == QWmat
    assert layer1.geometry == wkt_box
    assert layer1.__dict__['new_property'] == 'new_property'

    assert layer2.width == width2
    assert layer2.role == role2
    assert layer2.material == Bmat
    assert layer2.geometry == wkt_box

    assert layer3.width == width3
    assert layer3.role == role3
    assert layer3.material == i_GaAs
    assert layer3.geometry == wkt_box

    sn = random.uniform(1e6, 1e7)
    sp = random.uniform(1e6, 1e7)
    junction1 = Junction([layer1, layer2, layer3], sn=sn, sp=sp, T=T, kind='PDD')

    assert junction1.__len__() == 3
    assert junction1.__dict__ == {'sn': sn, 'sp': sp, 'T': T, 'kind': 'PDD'}
    assert junction1[0] == layer1
    assert junction1[1] == layer2
    assert junction1[2] == layer3

    tunnel1 = TunnelJunction([layer1, layer2, layer3], sn=sn, sp=sp, T=T, kind='PDD')

    assert tunnel1.__len__() == 3
    assert tunnel1.__dict__ == {'sn': sn, 'sp': sp, 'T': T, 'kind': 'PDD', 'R': 1e-16, 'pn': True}
    assert tunnel1[0] == layer1
    assert tunnel1[1] == layer2
    assert tunnel1[2] == layer3


def test_structure(Bmat, layer1, layer2, layer3):
    structure1 = Structure([layer1, layer2], substrate=Bmat, T=T)

    assert structure1.__len__() == 2
    assert structure1.__dict__ == {'substrate': Bmat, 'T': T, 'labels': [None, None]}
    assert structure1.width() == width1 + width2
    assert structure1[0] == layer1
    assert structure1[1] == layer2

    structure2 = Structure([], substrate=Bmat, T=T)
    structure2.append(layer1, layer_label='layer1')

    assert structure2.__len__() == 1
    assert structure2.__dict__ == {'substrate': Bmat, 'T': T, 'labels': ['layer1']}
    assert structure2[0] == layer1

    structure2.append_multiple([layer2, layer3], layer_labels=['layer2', 'layer3'])

    assert structure2.__len__() == 3
    assert structure2.__dict__ == {'substrate': Bmat, 'T': T, 'labels': ['layer1', 'layer2', 'layer3']}
    assert structure2[0] == layer1
    assert structure2[1] == layer2
    assert structure2[2] == layer3
    assert structure2.width() == width1 + width2 + width3
    assert structure2.relative_widths()['layer1'] == width1 / structure2.width()
    assert structure2.relative_widths()['layer2'] == width2 / structure2.width()
    assert structure2.relative_widths()['layer3'] == width3 / structure2.width()

    structure3 = Structure([])
    structure3.append(layer1, layer_label='layer1', repeats=2)

    assert structure3.__len__() == 2
    assert structure3.__dict__ == {'labels': ['layer1', 'layer1']}
    assert structure3[0] == layer1
    assert structure3[1] == layer1
    assert structure3.width() == width1 * 2
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
    assert structure4.width() == width1 * 2 + width2 * 2
    # Below currently fails due to a bug when specifying labels and repeating more than once
    # assert structure4.relative_widths()['layer1'] == width1 / structure4.width()
    # assert structure4.relative_widths()['layer2'] == width2 / structure4.width()


def test_material_to_str(QWmat, Bmat, i_GaAs):
    assert SolcoreMaterialToStr(QWmat) == QWmat_structure
    assert SolcoreMaterialToStr(Bmat) == Bmat_structure
    assert SolcoreMaterialToStr(i_GaAs) == i_GaAs_structure


def test_to_material(QWmat, Bmat, i_GaAs):
    QWmat_material = ToSolcoreMaterial(QWmat_structure, T, True)
    assert QWmat_material.__class__.__name__ == QWmat_material_name
    assert QWmat_material.__dict__ == QWmat.__dict__

    Bmat_material = ToSolcoreMaterial(Bmat_structure, T, True)
    assert Bmat_material.__class__.__name__ == Bmat_material_name
    assert Bmat_material.__dict__ == Bmat.__dict__

    i_GaAs_material = ToSolcoreMaterial(i_GaAs_structure, T, True)
    assert i_GaAs_material.__class__.__name__ == i_GaAs_material_name
    assert i_GaAs_material.__dict__ == i_GaAs.__dict__


def test_inline_composition(device):
    assert InLineComposition(device['layers'][0]) == 'In0.2GaAs'
    assert InLineComposition(device['layers'][1]) == 'GaAsP0.1'
    assert InLineComposition(device['layers'][2]) == 'GaAs'


def test_to_layer(QWmat, device):
    device_layer = device['layers'][0]
    composition = device_layer['properties']['composition']
    layer = ToLayer(width1, ToSolcoreMaterial(composition, T), role1)

    assert layer.width == width1
    assert layer.role == role1
    assert layer.material.__str__() == QWmat.__str__()


def test_to_structure(device):
    structure = ToStructure(device)

    assert structure.__len__() == 3
    assert structure.width() == width1 + width2 + width3
