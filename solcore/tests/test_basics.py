from unittest import TestCase

import solcore

class TestConstants(TestCase):
    def test_01_constants_correctly_imported(self):
        q = 1.60217646e-19
        self.assertTrue(q == solcore.constants.q)

class TestMaterial(TestCase):
    def test_02_material_correctly_imported(self):
        GaAs = solcore.material('GaAs')(T=300)
        beta_Gamma_GaAs = 204
        self.assertTrue(GaAs.beta_Gamma == beta_Gamma_GaAs)

class TestUnits(TestCase):
    def test_03_units_correctly_calculated(self):
        a_nm = 1239.8417166827828
        self.assertAlmostEqual(a_nm, solcore.eVnm(1))


### Quantum mechanics tests
import solcore
from solcore import material, si
from solcore.constants import vacuum_permittivity, q
from solcore.structure import Layer, Structure
import solcore.quantum_mechanics as QM
import numpy as np

# Energies in meV
my_energies = {'Elh': np.array([-770.]), 'Ee': np.array([577.]),
               'Ehh': np.array([-677., -731.])}

my_absorption = [1.267056856187291, 988133.91705224814]  # Energy in meV and absorption coefficent in m-1


class TestQM(TestCase):
    def test_61_kp_bulk_solvers(self):
        GaAs = solcore.material("GaAs")(T=300)
        GaAsP = solcore.material("GaAsP")(P=0.3, T=300)
        InGaAs = solcore.material("InGaAs")(In=0.2, T=300)

        # Testing QM.kp_bands
        bands_GaAsP = QM.kp_bands(GaAsP, GaAs, graph=False, fit_effective_mass=True, effective_mass_direction="X",
                                  return_so=True)
        bands_InGaAs = QM.kp_bands(InGaAs, GaAs, graph=False, fit_effective_mass=True, effective_mass_direction="X",
                                   return_so=True)

        expected_bands_GaAsP = (1.2168382480631407e-19, -1.5452519153253004e-19, -1.4042149045828435e-19,
                                -1.9192138182935611e-19, 8.0093555784846597e-32, 1.2472835955929216e-31,
                                2.6423749777877535e-31, 1.2393634061184521e-31)
        expected_bands_InGaAs = (8.6764014773233634e-20, -1.0573103504669798e-19, -1.1984351916698905e-19,
                                 -1.6993543036257329e-19, 6.6193386922591731e-32, 1.3576713980579555e-31,
                                 8.0904387208083259e-32, 1.2268187169919973e-31)

        for i in range(len(bands_GaAsP)):
            self.assertAlmostEqual(bands_GaAsP[i], expected_bands_GaAsP[i])
            self.assertAlmostEqual(bands_InGaAs[i], expected_bands_InGaAs[i])

        # Testing QM.KPbands and QM.fit_effective_masses
        edges_GaAsP = QM.KPbands(GaAsP, GaAs, fraction=0.2, return_edges_only=True)
        bands_GaAsP = QM.KPbands(GaAsP, GaAs, fraction=0.2)
        masses_GaAsP = QM.fit_effective_masses(bands_GaAsP, GaAsP, GaAs, plot_result=False)

        edges_InGaAs = QM.KPbands(InGaAs, GaAs, fraction=0.2, return_edges_only=True)
        bands_InGaAs = QM.KPbands(InGaAs, GaAs, fraction=0.2)
        masses_InGaAs = QM.fit_effective_masses(bands_InGaAs, InGaAs, GaAs, plot_result=False)

        expected_edges_GaAsP = (1.2168382480631407e-19, -1.5452519153253004e-19, -1.4042149045828435e-19,
                                -1.9192138182935611e-19)

        expected_edges_InGaAs = (8.6764014773233634e-20, -1.0573103504669798e-19, -1.1984351916698905e-19,
                                 -1.6993543036257329e-19)

        expected_masses_GaAsP = (8.049577422084102e-32, 1.2627430248682043e-31, 2.6577242586172804e-31,
                                 1.2305748108835472e-31)

        expected_masses_InGaAs = (6.6895885457875e-32, 1.3994390560400583e-31, 8.142667105522975e-32,
                                  1.2060355194525871e-31)

        for i in range(len(edges_GaAsP)):
            self.assertAlmostEqual(edges_GaAsP[i], expected_edges_GaAsP[i])
            self.assertAlmostEqual(edges_InGaAs[i], expected_edges_InGaAs[i])
            self.assertAlmostEqual(masses_GaAsP[i], expected_masses_GaAsP[i])
            self.assertAlmostEqual(masses_InGaAs[i], expected_masses_InGaAs[i])

        # Testing QM.kp8x8_bulk
        bands_GaAsP = QM.kp8x8_bulk(GaAsP, GaAs)
        bands_InGaAs = QM.kp8x8_bulk(InGaAs, GaAs)

        expected_bands_GaAsP = (1.2168382480631407e-19, -1.5452519153253004e-19, -1.4042149045828435e-19,
                                -1.9192138182935611e-19, 8.138445281947437e-32, 1.3543674202428507e-31,
                                2.848952594319033e-31, 1.145125159048442e-31)

        expected_bands_InGaAs = (8.6764014773233634e-20, -1.0573103504669798e-19, -1.1984351916698905e-19,
                                 -1.6993543036257329e-19, 6.77957483053393e-32, 1.6988821114765817e-31,
                                 7.820551493038613e-32, 1.1461138067300424e-31)

        for i in range(len(bands_GaAsP)):
            self.assertAlmostEqual(bands_GaAsP[i], expected_bands_GaAsP[i])
            self.assertAlmostEqual(bands_InGaAs[i], expected_bands_InGaAs[i])

    def test_64_quantum_mechanics_schrodinger(self):
        bulk = material("GaAs")(T=293)
        barrier = material("GaAsP")(T=293, P=0.1)

        bulk.strained = False
        barrier.strained = True

        top_layer = Layer(width=si("30nm"), material=bulk)
        inter = Layer(width=si("3nm"), material=bulk)
        barrier_layer = Layer(width=si("15nm"), material=barrier)
        bottom_layer = top_layer

        E = np.linspace(1.15, 1.5, 300) * q
        alfas = np.zeros((len(E), 6))

        alfas[:, 0] = E / q

        alpha_params = {
            "well_width": si("7.2nm"),
            "theta": 0,
            "eps": 12.9 * vacuum_permittivity,
            "espace": E,
            "hwhm": si("6meV"),
            "dimensionality": 0.16,
            "line_shape": "Gauss"
        }

        QW = material("InGaAs")(T=293, In=0.2)
        QW.strained = True
        well_layer = Layer(width=si("7.2nm"), material=QW)

        my_structure = Structure([top_layer, barrier_layer, inter] + 1 * [well_layer, inter, barrier_layer, inter] +
                                   [bottom_layer], substrate = bulk)

        band_edge, bands = QM.schrodinger(my_structure, quasiconfined=0,
                             num_eigenvalues=20, alpha_params=alpha_params, calculate_absorption=True)

        for key in band_edge['E']:
            for i in range(len(band_edge['E'][key])):
                band_edge['E'][key][i] = solcore.asUnit(band_edge['E'][key][i], 'eV') * 1000
                band_edge['E'][key][i] = round(band_edge['E'][key][i])

        Ehh = np.all(np.equal(band_edge['E']['Ehh'], my_energies['Ehh']))
        Elh = np.all(np.equal(band_edge['E']['Elh'], my_energies['Elh']))
        Ee = np.all(np.equal(band_edge['E']['Ee'], my_energies['Ee']))

        idx = 100
        out = [band_edge['alpha'][0][idx] / q, band_edge['alpha'][1][idx]]

        # Test over the energies
        self.assertTrue(Ehh and Elh and Ee)
        # Test over the absorption coefficent at a given energy
        for i, data in enumerate(out):
            self.assertAlmostEqual(out[i], my_absorption[i])


### Absorption calculator related tests

from solcore import material, si
from solcore.structure import Structure, Layer
from solcore.absorption_calculator import create_adachi_alpha
from solcore.absorption_calculator import calculate_rat, calculate_ellipsometry, calculate_absorption_profile
from solcore.absorption_calculator.dielectric_constant_models import DielectricConstantModel, Poles, Drude
import numpy as np


class TestAbsorption(TestCase):
    def test_41_adachi_absorption(self):
        material_name = "InGaAs"
        solcore_material = material(material_name)(T=300, In=0.1)

        E, nn, kk, adachi_alpha_data = create_adachi_alpha(solcore_material, T=300)

        E /= 1.6e-19

        idx = np.argmin(abs(E - 2))

        out = [E[idx], nn[idx], kk[idx]]

        data = [2.0003267755918639, 3.6288603412514253, 0.30754355994545018]

        self.assertAlmostEqual(data, out)

    def test_42_TMM_rat(self):
        GaAs = material('GaAs')(T=300)

        my_structure = Structure([
            Layer(si(3000, 'nm'), material=GaAs),
        ])

        wavelength = np.linspace(450, 1100, 300)
        idx = np.argmin(abs(wavelength - 800))

        out = calculate_rat(my_structure, wavelength, coherent=True, no_back_reflexion=False)
        out = (out['R'][idx], out['A'][idx], out['T'][idx])

        data = (0.33328918841743332, 0.65996607786373396, 0.0067447337188326914)

        for i in range(3):
            self.assertAlmostEqual(data[i], out[i])

    def test_43_TMM_absorption_profile(self):
        GaAs = material('GaAs')(T=300)

        my_structure = Structure([
            Layer(si(3000, 'nm'), material=GaAs),
            Layer(si(300, 'um'), material=GaAs),
        ])

        out = calculate_absorption_profile(my_structure, [800], z_limit=3000, steps_size=20)

        data = (0.00093920198054733134, 0.00091329190431268755, 0.00088809661793623094, 0.00086359640227330484,
                0.00083977208217782804, 0.00081660501149482764, 0.0007940770584669893, 0.00077217059154379222,
                0.00075086846558214263, 0.00073015400842769127, 0.00071001100786633451, 0.00069042369893569059,
                0.00067137675158661923, 0.00065285525868512457, 0.00063484472434525627, 0.00061733105258387328,
                0.00060030053628839001, 0.00058373984648887986, 0.00056763602192612503, 0.00055197645890746013,
                0.00053674890144246557, 0.0005219414316507909, 0.00050754246043459913, 0.00049354071840833585,
                0.00047992524707871981, 0.00046668539026805409, 0.0004538107857741483, 0.00044129135726031623,
                0.00042911730636911071, 0.00041727910505361793, 0.00040576748812031177, 0.00039457344597762961,
                0.00038368821758459675, 0.00037310328359397839, 0.00036281035968459497, 0.00035280139007758007,
                0.00034306854123150837, 0.00033360419571145685, 0.00032440094622720458, 0.00031545158983589954,
                0.00030674912230466134, 0.00029828673262870248, 0.00029005779770068137, 0.00028205587712711315,
                0.00027427470818778237, 0.00026670820093421067, 0.00025935043342334711, 0.00025219564708274435,
                0.00024523824220360293, 0.00023847277355814448, 0.00023189394613789383, 0.00022549661100952902,
                0.0002192757612850577, 0.00021322652820316562, 0.00020734417731867015, 0.00020162410479709653,
                0.00019606183381147725, 0.00019065301103855347, 0.0001853934032516379, 0.00018027889400747007,
                0.00017530548042447405, 0.00017046927004989423, 0.00016576647781335848, 0.00016119342306448588,
                0.00015674652669221659, 0.00015242230832361372, 0.00014821738359994165, 0.00014412846152789071,
                0.00014015234190387394, 0.00013628591280938143, 0.00013252614817542982, 0.0001288701054142039,
                0.00012531492311603331, 0.00012185781880990432, 0.00011849608678575335, 0.00011522709597683633,
                0.00011204828790051968, 0.00010895717465587773, 0.00010595133697653208, 0.00010302842233720751,
                0.00010018614311252307, 9.7422274786577231e-05, 9.4734654211925496e-05, 9.2121177916588886e-05,
                8.9579800457766819e-05, 8.7108532820967001e-05, 8.470544086329888e-05, 8.2368643799712795e-05,
                8.009631273099949e-05, 7.7886669212397987e-05, 7.573798386169243e-05, 7.3648575005706959e-05,
                7.1616807364140996e-05, 6.9641090769713272e-05, 6.7719878923614451e-05, 6.5851668185292527e-05,
                6.4034996395625842e-05, 6.2268441732560816e-05, 6.0550621598319883e-05, 5.8880191537308663e-05,
                5.7255844183874585e-05, 5.5676308239094561e-05, 5.41403474757902e-05, 5.2646759770991723e-05,
                5.1194376165094236e-05, 4.9782059946968693e-05, 4.8408705764312587e-05, 4.7073238758543685e-05,
                4.5774613723559514e-05, 4.4511814287704784e-05, 4.3283852118305772e-05, 4.2089766148149713e-05,
                4.0928621823303459e-05, 3.9799510371682887e-05, 3.8701548091800482e-05, 3.7633875661134488e-05,
                3.6595657463578247e-05, 3.5586080935443529e-05, 3.4604355929505788e-05, 3.3649714096593824e-05,
                3.2721408284239568e-05, 3.1818711951917646e-05, 3.0940918602416916e-05, 3.0087341228898868e-05,
                2.9257311777210295e-05, 2.8450180623029266e-05, 2.7665316063435288e-05, 2.6902103822505617e-05,
                2.6159946570550893e-05, 2.5438263456613905e-05, 2.4736489653865175e-05, 2.4054075917540213e-05,
                2.3390488155071715e-05, 2.2745207008081107e-05, 2.211772744590139e-05, 2.1507558370314028e-05,
                2.0914222231189692e-05, 2.0337254652732693e-05, 1.9776204070036316e-05, 1.9230631375664536e-05,
                1.8700109575983667e-05, 1.8184223456974879e-05, 1.7682569259266036e-05, 1.7194754362128619e-05,
                1.6720396976192213e-05, 1.6259125844636267e-05, 1.5810579952625165e-05, 1.5374408244759177e-05,
                1.4950269350320186e-05, 1.4537831316097222e-05)

        for i in range(len(out['absorption'][0])):
            self.assertAlmostEqual(data[i], out['absorption'][0][i])

    def test_44_TMM_ellipsometry(self):
        GaAs = material('GaAs')(T=300)

        my_structure = Structure([
            Layer(si(3000, 'nm'), material=GaAs),
            Layer(si(300, 'um'), material=GaAs),
        ])

        wavelength = np.linspace(450, 1100, 300)
        idx = np.argmin(abs(wavelength - 800))

        angles = [60, 65, 70]
        out = calculate_ellipsometry(my_structure, wavelength, angle=angles)

        data = (22.2849089096, 181.488417672, 16.4604621886, 182.277656469, 9.10132195668, 184.509752582)

        for i in range(len(angles)):
            self.assertAlmostEqual(data[2 * i], out['psi'][idx, i])
            self.assertAlmostEqual(data[2 * i + 1], out['Delta'][idx, i])

    def test_45_TMM_dielectric_model(self):
        drud = Drude(An=24.317000, Brn=0.125740)

        model = DielectricConstantModel(e_inf=3.4837, oscillators=[drud])

        wavelength = 2 * np.logspace(3, 4, 10)
        n = model.n_and_k(wavelength)

        data = (0.37377710 + 2.0726883j, 0.53555835 + 3.03640188j, 0.81383835 + 4.12781828j, 1.25563998 + 5.39395751j,
                1.92396319 + 6.8504966j, 2.88306008 + 8.48061105j, 4.17692618 + 10.23820188j, 5.81068443 + 12.06785126j,
                7.75184851 + 13.93546615j, 9.95389660 + 15.84816722j)

        for i in range(len(wavelength)):
            self.assertAlmostEqual(data[i], n[i])

    def test_46_sopra_absorption(self):
        from solcore.absorption_calculator import sopra_database

        # Import material constant data for Gallium Arsenide :: Do this by placing the material name as the sole argument...
        SOPRA_Material = sopra_database("GaAs")

        # Can also load alpha data...
        GaAs_alpha = SOPRA_Material.load_alpha()

        out = GaAs_alpha[1][10]

        data = 163666134.03339368

        self.assertAlmostEqual(data, out)


#### Poisson_drift_diffusion related tests
import numpy as np
import os
import tempfile

import solcore.poisson_drift_diffusion as PDD
from solcore.solar_cell import SolarCell, default_GaAs
from solcore.structure import Layer, Junction
from solcore import si
from solcore import material
from solcore.light_source import LightSource
from solcore.solar_cell_solver import solar_cell_solver

T = 298

substrate = material('GaAs')(T=T)


def AlGaAs(T):
    # We create the other materials we need for the device
    window = material('AlGaAs')(T=T, Na=5e24, Al=0.8)
    p_AlGaAs = material('AlGaAs')(T=T, Na=1e24, Al=0.4)
    n_AlGaAs = material('AlGaAs')(T=T, Nd=8e22, Al=0.4)
    bsf = material('AlGaAs')(T=T, Nd=2e24, Al=0.6)

    output = Junction([Layer(width=si('30nm'), material=window, role="Window"),
                       Layer(width=si('150nm'), material=p_AlGaAs, role="Emitter"),
                       Layer(width=si('1000nm'), material=n_AlGaAs, role="Base"),
                       Layer(width=si('200nm'), material=bsf, role="BSF")], sn=1e6, sp=1e6, T=T, kind='PDD')

    return output


Vin = np.linspace(-2, 2.61, 201)
V = np.linspace(0, 2.6, 300)
wl = np.linspace(350, 1000, 301) * 1e-9
light_source = LightSource(source_type='standard', version='AM1.5g', x=wl,
                           output_units='photon_flux_per_m', concentration=1)


class TestPDD(TestCase):
    def test_92_light_iv(self):
        answer = [141.434980729, 2.46952886616, 0.91434377329, 319.359951949, 2.29565217391, 139.115130584,
                  0.319241623262]
        with tempfile.TemporaryDirectory(prefix="tmp", suffix="_sc3TESTS") as working_directory:
            filename = os.path.join(working_directory, 'solcore_log.txt')
            PDD.log(filename)

            my_solar_cell = SolarCell([AlGaAs(T), default_GaAs(T)], T=T, R_series=0, substrate=substrate)
            solar_cell_solver(my_solar_cell, 'iv',
                              user_options={'T_ambient': T, 'db_mode': 'boltzmann', 'voltages': V, 'light_iv': True,
                                            'wavelength': wl, 'optics_method': 'BL', 'mpp': True,
                                            'internal_voltages': Vin,
                                            'light_source': light_source})

            output = [my_solar_cell.iv.Isc, my_solar_cell.iv.Voc, my_solar_cell.iv.FF, my_solar_cell.iv.Pmpp,
                      my_solar_cell.iv.Vmpp, my_solar_cell.iv.Impp, my_solar_cell.iv.Eta]

        for i in range(len(output)):
            self.assertAlmostEqual(output[i], answer[i])

    def test_93_qe(self):
        answer = [0.9831923128532823, 1.315965183418519e-13, 0.9672990699170962, 0.032767290395462376]
        with tempfile.TemporaryDirectory(prefix="tmp", suffix="_sc3TESTS") as working_directory:
            filename = os.path.join(working_directory, 'solcore_log.txt')
            PDD.log(filename)

            my_solar_cell = SolarCell([AlGaAs(T), default_GaAs(T)], T=T, R_series=0, substrate=substrate)

            solar_cell_solver(my_solar_cell, 'qe',
                              user_options={'T_ambient': T, 'db_mode': 'boltzmann', 'voltages': V, 'light_iv': True,
                                            'wavelength': wl, 'optics_method': 'BL', 'mpp': True,
                                            'internal_voltages': Vin,
                                            'light_source': light_source})

            output = [my_solar_cell[0].eqe(500e-9), my_solar_cell[0].eqe(800e-9), my_solar_cell[1].eqe(700e-9),
                   my_solar_cell[1].eqe(900e-9)]

        for i in range(len(output)):
            self.assertAlmostEqual(output[i], answer[i])
