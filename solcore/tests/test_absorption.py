""" Absorption calculator related tests
"""
from pytest import approx, mark

from solcore import material, si
from solcore.structure import Structure, Layer
from solcore.absorption_calculator import create_adachi_alpha
from solcore.absorption_calculator import (
    calculate_rat,
    calculate_ellipsometry,
    calculate_absorption_profile,
)
from solcore.state import State
from solcore.solar_cell_solver import prepare_solar_cell
from solcore.absorption_calculator.dielectric_constant_models import (
    DielectricConstantModel,
    Drude,
)

from solcore.solar_cell_solver import solar_cell_solver
from solcore.solar_cell import SolarCell
from solcore.material_system import create_new_material
from solcore.absorption_calculator import download_db, search_db
from solcore.absorption_calculator.nk_db import nkdb_load_n
from solcore.config_tools import add_source
from solcore.optics import solve_tmm

import os

import numpy as np


def test_adachi_absorption():
    material_name = "InGaAs"
    solcore_material = material(material_name)(T=300, In=0.1)

    E, nn, kk, adachi_alpha_data = create_adachi_alpha(solcore_material, T=300)

    E /= 1.6e-19

    idx = np.argmin(abs(E - 2))

    out = [E[idx], nn[idx], kk[idx]]

    data = [2.0003267755918639, 3.6288603412514253, 0.30754355994545018]

    assert all([d == approx(o) for d, o in zip(data, out)])


def test_TMM_rat():
    GaAs = material("GaAs")(T=300)

    my_structure = Structure([Layer(si(3000, "nm"), material=GaAs)])

    wavelength = np.linspace(450, 1100, 300)
    idx = np.argmin(abs(wavelength - 800))

    out = calculate_rat(
        my_structure, wavelength, coherent=True, no_back_reflexion=False
    )
    out = (out["R"][idx], out["A"][idx], out["T"][idx])

    data = (0.33328918841743332, 0.65996607786373396, 0.0067447337188326914)

    assert all([d == approx(o) for d, o in zip(data, out)])


def test_TMM_absorption_profile():
    GaAs = material("GaAs")(T=300)

    my_structure = Structure(
        [Layer(si(3000, "nm"), material=GaAs), Layer(si(300, "um"), material=GaAs)]
    )

    out = calculate_absorption_profile(
        my_structure, np.array([800]), z_limit=3000, steps_size=20
    )

    data = (
        0.00093920198054733134,
        0.00091329190431268755,
        0.00088809661793623094,
        0.00086359640227330484,
        0.00083977208217782804,
        0.00081660501149482764,
        0.0007940770584669893,
        0.00077217059154379222,
        0.00075086846558214263,
        0.00073015400842769127,
        0.00071001100786633451,
        0.00069042369893569059,
        0.00067137675158661923,
        0.00065285525868512457,
        0.00063484472434525627,
        0.00061733105258387328,
        0.00060030053628839001,
        0.00058373984648887986,
        0.00056763602192612503,
        0.00055197645890746013,
        0.00053674890144246557,
        0.0005219414316507909,
        0.00050754246043459913,
        0.00049354071840833585,
        0.00047992524707871981,
        0.00046668539026805409,
        0.0004538107857741483,
        0.00044129135726031623,
        0.00042911730636911071,
        0.00041727910505361793,
        0.00040576748812031177,
        0.00039457344597762961,
        0.00038368821758459675,
        0.00037310328359397839,
        0.00036281035968459497,
        0.00035280139007758007,
        0.00034306854123150837,
        0.00033360419571145685,
        0.00032440094622720458,
        0.00031545158983589954,
        0.00030674912230466134,
        0.00029828673262870248,
        0.00029005779770068137,
        0.00028205587712711315,
        0.00027427470818778237,
        0.00026670820093421067,
        0.00025935043342334711,
        0.00025219564708274435,
        0.00024523824220360293,
        0.00023847277355814448,
        0.00023189394613789383,
        0.00022549661100952902,
        0.0002192757612850577,
        0.00021322652820316562,
        0.00020734417731867015,
        0.00020162410479709653,
        0.00019606183381147725,
        0.00019065301103855347,
        0.0001853934032516379,
        0.00018027889400747007,
        0.00017530548042447405,
        0.00017046927004989423,
        0.00016576647781335848,
        0.00016119342306448588,
        0.00015674652669221659,
        0.00015242230832361372,
        0.00014821738359994165,
        0.00014412846152789071,
        0.00014015234190387394,
        0.00013628591280938143,
        0.00013252614817542982,
        0.0001288701054142039,
        0.00012531492311603331,
        0.00012185781880990432,
        0.00011849608678575335,
        0.00011522709597683633,
        0.00011204828790051968,
        0.00010895717465587773,
        0.00010595133697653208,
        0.00010302842233720751,
        0.00010018614311252307,
        9.7422274786577231e-05,
        9.4734654211925496e-05,
        9.2121177916588886e-05,
        8.9579800457766819e-05,
        8.7108532820967001e-05,
        8.470544086329888e-05,
        8.2368643799712795e-05,
        8.009631273099949e-05,
        7.7886669212397987e-05,
        7.573798386169243e-05,
        7.3648575005706959e-05,
        7.1616807364140996e-05,
        6.9641090769713272e-05,
        6.7719878923614451e-05,
        6.5851668185292527e-05,
        6.4034996395625842e-05,
        6.2268441732560816e-05,
        6.0550621598319883e-05,
        5.8880191537308663e-05,
        5.7255844183874585e-05,
        5.5676308239094561e-05,
        5.41403474757902e-05,
        5.2646759770991723e-05,
        5.1194376165094236e-05,
        4.9782059946968693e-05,
        4.8408705764312587e-05,
        4.7073238758543685e-05,
        4.5774613723559514e-05,
        4.4511814287704784e-05,
        4.3283852118305772e-05,
        4.2089766148149713e-05,
        4.0928621823303459e-05,
        3.9799510371682887e-05,
        3.8701548091800482e-05,
        3.7633875661134488e-05,
        3.6595657463578247e-05,
        3.5586080935443529e-05,
        3.4604355929505788e-05,
        3.3649714096593824e-05,
        3.2721408284239568e-05,
        3.1818711951917646e-05,
        3.0940918602416916e-05,
        3.0087341228898868e-05,
        2.9257311777210295e-05,
        2.8450180623029266e-05,
        2.7665316063435288e-05,
        2.6902103822505617e-05,
        2.6159946570550893e-05,
        2.5438263456613905e-05,
        2.4736489653865175e-05,
        2.4054075917540213e-05,
        2.3390488155071715e-05,
        2.2745207008081107e-05,
        2.211772744590139e-05,
        2.1507558370314028e-05,
        2.0914222231189692e-05,
        2.0337254652732693e-05,
        1.9776204070036316e-05,
        1.9230631375664536e-05,
        1.8700109575983667e-05,
        1.8184223456974879e-05,
        1.7682569259266036e-05,
        1.7194754362128619e-05,
        1.6720396976192213e-05,
        1.6259125844636267e-05,
        1.5810579952625165e-05,
        1.5374408244759177e-05,
        1.4950269350320186e-05,
        1.4537831316097222e-05,
    )

    assert all([d == approx(o) for d, o in zip(data, out["absorption"][0])])


def test_TMM_ellipsometry():
    GaAs = material("GaAs")(T=300)

    my_structure = Structure(
        [Layer(si(3000, "nm"), material=GaAs), Layer(si(300, "um"), material=GaAs)]
    )

    wavelength = np.linspace(450, 1100, 300)
    idx = np.argmin(abs(wavelength - 800))

    angles = [60, 65, 70]
    out = calculate_ellipsometry(my_structure, wavelength, angle=angles)

    data = (
        22.2849089096,
        181.488417672,
        16.4604621886,
        182.277656469,
        9.10132195668,
        184.509752582,
    )

    for i in range(len(angles)):
        assert data[2 * i] == approx(out["psi"][idx, i])
        assert data[2 * i + 1] == approx(out["Delta"][idx, i])


def test_TMM_dielectric_model():
    drud = Drude(An=24.317000, Brn=0.125740)

    model = DielectricConstantModel(e_inf=3.4837, oscillators=[drud])

    wavelength = 2 * np.logspace(3, 4, 10)
    n = model.n_and_k(wavelength)

    data = (
        0.37377710 + 2.0726883j,
        0.53555835 + 3.03640188j,
        0.81383835 + 4.12781828j,
        1.25563998 + 5.39395751j,
        1.92396319 + 6.8504966j,
        2.88306008 + 8.48061105j,
        4.17692618 + 10.23820188j,
        5.81068443 + 12.06785126j,
        7.75184851 + 13.93546615j,
        9.95389660 + 15.84816722j,
    )

    assert all([d == approx(o) for d, o in zip(data, n)])


def test_sopra_absorption():
    from solcore.absorption_calculator import sopra_database

    # Import material constant data for Gallium Arsenide :: Do this by placing the
    # material name as the sole argument...
    SOPRA_Material = sopra_database("GaAs")

    # Can also load alpha data...
    GaAs_alpha = SOPRA_Material.load_alpha()

    out = GaAs_alpha[1][10]

    data = 163666134.03339368

    assert data == approx(out)

def test_substrate_presence_A():
    wavelength = np.linspace(300, 800, 3) * 1e-9

    GaAs = material("GaAs")(T=300)

    my_structure = SolarCell([Layer(si(700, "nm"), material=GaAs)], substrate=GaAs)

    solar_cell_solver(my_structure, 'optics',
                      user_options={'wavelength': wavelength, 'optics_method': 'TMM',
                                    'no_back_reflexion': False})

    z_pos = np.linspace(0, my_structure.width, 10)

    A_subs = my_structure[0].layer_absorption

    my_structure = SolarCell([Layer(si(700, "nm"), material=GaAs)])

    solar_cell_solver(my_structure, 'optics',
                      user_options={'wavelength': wavelength, 'optics_method': 'TMM',
                                    'no_back_reflexion': False})

    A_nosubs = my_structure[0].layer_absorption


    A = np.vstack((A_subs, A_nosubs))

    A_data = np.array([[0.56610281, 0.62692985, 0.41923175],
                       [0.56610281, 0.62711355, 0.37837737]])

    assert all([d == approx(o) for d, o in zip(A, A_data)])

def test_BL_correction():

    wl = np.linspace(800, 950, 4) * 1e-9

    GaAs = material('GaAs')()

    thick_cell = SolarCell([Layer(material=GaAs, width=si('20um'))])

    opts = State()
    opts.position = None
    prepare_solar_cell(thick_cell, opts)
    position = np.arange(0, thick_cell.width, 1e-9)
    opts.position = position
    opts.recalculate_absorption = True
    opts.no_back_reflexion = False

    opts.BL_correction = False
    opts.wavelength = wl
    solve_tmm(thick_cell, opts)

    no_corr = thick_cell.absorbed

    opts.BL_correction = True

    solve_tmm(thick_cell, opts)

    with_corr = thick_cell.absorbed

    assert with_corr == approx(np.array([ 6.71457872e-01,  6.75496354e-01,  2.09738887e-01, 0]))
    assert no_corr == approx(np.array([ 6.71457872e-01,  6.75496071e-01,  2.82306407e-01, 0]))


def test_substrate_presence_profile():
    wavelength = np.linspace(300, 800, 3) * 1e-9

    GaAs = material("GaAs")(T=300)

    my_structure = SolarCell([Layer(si(700, "nm"), material=GaAs)], substrate=GaAs)

    solar_cell_solver(my_structure, 'optics',
                      user_options={'wavelength': wavelength, 'optics_method': 'TMM',
                                    'no_back_reflexion': False})

    z_pos = np.linspace(0, my_structure.width, 10)

    profile_subs = my_structure[0].absorbed(z_pos)

    my_structure = SolarCell([Layer(si(700, "nm"), material=GaAs)])

    solar_cell_solver(my_structure, 'optics',
                      user_options={'wavelength': wavelength, 'optics_method': 'TMM',
                                    'no_back_reflexion': False})

    profile_nosubs = my_structure[0].absorbed(z_pos)

    profile = np.vstack((profile_subs, profile_nosubs))

    profile_data = np.array([[4.69390061e+07, 4.31690443e+06, 9.39070619e+05],
       [7.00847453e+04, 2.53652234e+06, 8.42216253e+05],
       [1.04643340e+02, 1.49040720e+06, 7.55351304e+05],
       [1.56242133e-01, 8.75731920e+05, 6.77445478e+05],
       [2.33283058e-04, 5.14561641e+05, 6.07574744e+05],
       [3.54211961e-07, 3.02759196e+05, 5.45062846e+05],
       [5.28873825e-10, 1.77894940e+05, 4.88845864e+05],
       [7.89658664e-13, 1.04527326e+05, 4.38427018e+05],
       [1.17903096e-15, 6.14180570e+04, 3.93208298e+05],
       [0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                             [4.69390061e+07, 4.28007419e+06, 4.79821200e+05],
                             [7.00847453e+04, 2.54503334e+06, 6.73708148e+05],
                             [1.04643340e+02, 1.52813445e+06, 9.69600309e+05],
                             [1.56242133e-01, 9.06948301e+05, 4.94806809e+05],
                             [2.33283058e-04, 5.11344674e+05, 2.57215183e+05],
                             [3.54211961e-07, 2.67269785e+05, 7.11553179e+05],
                             [5.28873827e-10, 1.38657803e+05, 6.48896832e+05],
                             [7.89660581e-13, 9.55930119e+04, 1.30180135e+05],
                             [1.18135121e-15, 9.46668507e+04, 3.29375923e+05],
                             [0.00000000e+00, 0.00000000e+00, 0.00000000e+00]]
                            )



    assert all([d == approx(o) for d, o in zip(profile, profile_data)])

# TODO: the following tests for custom materials do not work as they require changes to the user config file.
# It is possible the downloading of the database for test_database_materials is also an issue.

@mark.skip
def test_define_material():
    home_folder = os.path.expanduser('~')
    custom_nk_path = os.path.join(home_folder, 'Solcore/custommats')
    param_path = os.path.join(home_folder, 'Solcore/custom_params.txt')

    add_source('Others', 'custom_mats', custom_nk_path)
    add_source('Parameters', 'custom', param_path)
    this_dir = os.path.split(__file__)[0]
    create_new_material('SiGeSn', os.path.join(this_dir, 'SiGeSn_n.txt'), os.path.join(this_dir, 'SiGeSn_k.txt'), os.path.join(this_dir, 'SiGeSn_params.txt'))

@mark.skip
def test_use_material():
    SiGeSn = material('SiGeSn')()
    assert SiGeSn.n(400e-9) == approx(4.175308391752484)
    assert SiGeSn.k(400e-9) == approx(2.3037424963866306)

@mark.skip
def test_database_materials():
    home_folder = os.path.expanduser('~')
    nk_db_path = os.path.join(home_folder, 'Solcore/NK.db')

    add_source('Others', 'nk', nk_db_path)
    download_db(confirm=True)
    wl, n = nkdb_load_n(2683) # Should be carbon, from Phillip
    n_data = np.array([0.58321493, 0.57867586, 0.57339939, 0.56591405, 0.56065836,
                       0.5563835 , 0.55229037, 0.54790994, 0.54244807, 0.53583666,
                       0.52780304, 0.52220281, 0.51743344, 0.5120067 , 0.50729163,
                       0.5046621 , 0.50166726, 0.49931595, 0.49501672, 0.49060884,
                       0.48855459, 0.48678095, 0.48522576, 0.48404383, 0.48359425,
                       0.48357126, 0.48406997, 0.48538076, 0.48643123, 0.48607976,
                       0.48596997, 0.48630729, 0.48839097, 0.49108562, 0.49244456,
                       0.49254289, 0.49350922, 0.4974711 , 0.50123933, 0.50465428,
                       0.50913919, 0.51323753, 0.51583738, 0.51988073, 0.52478609,
                       0.52834711, 0.53210744, 0.53641754, 0.541486  , 0.54850356,
                       0.55434881, 0.55925238, 0.56516219, 0.5689907 , 0.57441913,
                       0.57926509, 0.58287789, 0.58646844, 0.58975747, 0.59105466,
                       0.58655149, 0.5786466 , 0.56872678, 0.55518287, 0.53712377,
                       0.51098487, 0.49047403, 0.47911691, 0.45812142, 0.46345176,
                       0.44772251, 0.45534214, 0.4714259 , 0.48511583, 0.49380479,
                       0.4982652 , 0.49772907, 0.50612376, 0.51912528, 0.53467711,
                       0.55029993, 0.56862873, 0.58502876, 0.60229328, 0.61909644,
                       0.63547993, 0.65545356, 0.67805267, 0.69756522, 0.71659019,
                       0.73623614, 0.75563388, 0.76949687, 0.78444575, 0.81244567,
                       0.83941617, 0.86525298, 0.88673291, 0.90652684, 0.92235589,
                       0.94089495, 0.9708571 , 0.98590931, 1.00269955, 1.02097684,
                       1.03623883, 1.05626117, 1.08581851, 1.10194641, 1.1354526 ,
                       1.21260589, 1.34344321, 1.68147323, 2.08187171, 2.55437527,
                       2.84776234, 3.05673423, 3.27833604, 3.43571649, 3.50768429,
                       3.49100363, 3.5007246 , 3.49998604, 3.47187414, 3.45372454,
                       3.42135457, 3.37725267, 3.32462109, 3.29181609, 3.262685  ,
                       3.24502588, 3.23124178, 3.22611077, 3.22317221, 3.23477168,
                       3.25217791, 3.28216011, 3.32475433, 3.37108112, 3.42448917,
                       3.28385656, 3.10966919, 2.98894561, 2.9223206 , 2.86291625,
                       2.81756759, 2.77148488, 2.727706  , 2.69074153, 2.66695519,
                       2.63873644, 2.60595664, 2.5817649 , 2.56088461, 2.54132446,
                       2.5236323 , 2.50581324, 2.48786656, 2.47332368, 2.46498682,
                       2.45505193, 2.44125173, 2.43147897, 2.42707643, 2.42085729,
                       2.41142904, 2.40654316, 2.4062024 , 2.40003542, 2.39052295,
                       2.38812479, 2.38778349, 2.38324149, 2.37852055, 2.37817577,
                       2.37783305])

    assert all([d == approx(o) for d, o in zip(n, n_data)])
