import random
from solcore import material
from solcore.structure import Junction, Layer, SolcoreMaterialToStr, Structure, ToLayer, ToSolcoreMaterial, TunnelJunction

T = 300
QWmat_material_name = 'InGaAs'
QWmat_in = 0.2
QWmat = material(QWmat_material_name)(T=T, In=QWmat_in)
QWmat_structure = {'material': QWmat_material_name, 'element': 'In', 'fraction': QWmat_in}
Bmat_material_name = 'GaAsP'
Bmat_p = 0.1
Bmat = material(Bmat_material_name)(T=T, P=Bmat_p)
Bmat_structure = {'material': Bmat_material_name, 'element': 'P', 'fraction': Bmat_p}
i_GaAs_material_name = 'GaAs'
i_GaAs = material(i_GaAs_material_name)(T=T)
i_GaAs_structure = {'material': i_GaAs_material_name}

available_roles = ['Barrier', 'Base', 'BSF', 'Emitter', 'Intrinsic', 'Window']

wkt_box = 'POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))'

def test_layer_and_junction():
    width1 = random.uniform(1e-9, 1e-8)
    role1 = random.choice(available_roles)
    layer1 = Layer(width1, QWmat, role1, wkt_box, new_property='new_property')

    assert layer1.width == width1
    assert layer1.role == role1
    assert layer1.material == QWmat
    assert layer1.geometry == wkt_box
    assert layer1.__dict__['new_property'] == 'new_property'
    # TODO: see if the following can work:
    #  assert layer1 == ToLayer(width1, {'material': 'InGaAs', 'element': 'In', 'fraction': 0.2}, role1)

    width2 = random.uniform(1e-9, 1e-8)
    role2 = random.choice(available_roles)
    layer2 = Layer(width2, Bmat, role2, wkt_box)

    assert layer2.width == width2
    assert layer2.role == role2
    assert layer2.material == Bmat
    assert layer2.geometry == wkt_box

    width3 = random.uniform(1e-9, 1e-8)
    role3 = random.choice(available_roles)
    layer3 = Layer(width3, i_GaAs, role3, wkt_box)

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

def test_material_to_str():
    assert SolcoreMaterialToStr(QWmat) == QWmat_structure
    assert SolcoreMaterialToStr(Bmat) == Bmat_structure
    assert SolcoreMaterialToStr(i_GaAs) == i_GaAs_structure

def test_to_material():
    QWmat_material = ToSolcoreMaterial(QWmat_structure, T, True)
    assert QWmat_material.__class__.__name__ == QWmat_material_name
    assert QWmat_material.__dict__ == QWmat.__dict__

    Bmat_material = ToSolcoreMaterial(Bmat_structure, T, True)
    assert Bmat_material.__class__.__name__ == Bmat_material_name
    assert Bmat_material.__dict__ == Bmat.__dict__

    i_GaAs_material = ToSolcoreMaterial(i_GaAs_structure, T, True)
    assert i_GaAs_material.__class__.__name__ == i_GaAs_material_name
    assert i_GaAs_material.__dict__ == i_GaAs.__dict__
