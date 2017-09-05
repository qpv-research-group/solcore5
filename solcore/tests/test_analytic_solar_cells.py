from unittest import TestCase

import numpy as np
import os
from solcore import siUnits, material, si
from solcore.constants import h, c, vacuum_permittivity, kb, q
from solcore.structure import Structure, Junction, Layer
from solcore.state import State
from solcore.interpolate import interp1d
import solcore.analytic_solar_cells as ASC
from solcore.light_source import calculate_spectrum_spectral2

# ********************************
# Triple junction cell QE calculation

all_materials = []


# We need to build the solar cell layer by layer.
# We start from the AR coating. In this case, we load it from an an external file
def this_dir_file(f):
    return os.path.join(os.path.split(__file__)[0], f)


refl_nm = np.loadtxt(this_dir_file("MgF-ZnS_AR.csv"), unpack=True, delimiter=",")
refl_J = np.array((h * c / siUnits(refl_nm[0], "nm")[::-1], refl_nm[1][::-1]))

# Next is the window layer, made of AnInP. We load the absorption coefficent of AlInP from an external file
AlInP = material("AlInP")
window_material = AlInP(Al=0.52)
window_alpha = np.loadtxt(this_dir_file("alinp.csv"), unpack=True, delimiter=",")
window_material.alphaE = interp1d(x=siUnits(window_alpha[0], "eV"), y=window_alpha[1], bounds_error=False, fill_value=0)

all_materials.append(window_material)

# TOP CELL - GaInP
# Now we build the top cell, which requires the n and p sides of GaInP.
# We also load the absorption coefficient from an external file. We also add some extra parameters neded for the
# calculation such as the minority carriers diffusion lengths
InGaP = material("GaInP")
top_cell_n_material = InGaP(In=0.49, Nd=siUnits(2e18, "cm-3"), role="n")
top_cell_p_material = InGaP(In=0.49, Na=siUnits(1e17, "cm-3"), role="p")
top_cell_alpha = np.loadtxt(this_dir_file("in048ga052p.csv"), unpack=True, delimiter=",")
top_cell_n_material.alphaE = interp1d(x=siUnits(top_cell_alpha[0], "eV"), y=top_cell_alpha[1], bounds_error=False,
                                      fill_value=0)
top_cell_p_material.alphaE = top_cell_n_material.alphaE
top_cell_n_material.hole_minority_carrier_diffusion_length = si("200nm")
top_cell_p_material.electron_minority_carrier_diffusion_length = si("1um")

all_materials.append(top_cell_n_material)
all_materials.append(top_cell_p_material)

# MID CELL  - InGaAs
# The same thing.
InGaAs = material("InGaAs")
mid_cell_n_material = InGaAs(In=0.01, Nd=siUnits(3e18, "cm-3"), role="n")
mid_cell_p_material = InGaAs(In=0.01, Na=siUnits(1e17, "cm-3"), role="p")
mid_cell_alpha = np.loadtxt(this_dir_file("in01gaas.csv"), unpack=True, delimiter=",")
mid_cell_n_material.alphaE = interp1d(x=siUnits(mid_cell_alpha[0], "eV"), y=mid_cell_alpha[1], bounds_error=False,
                                      fill_value=0)
mid_cell_p_material.alphaE = mid_cell_n_material.alphaE
mid_cell_n_material.hole_minority_carrier_diffusion_length = si("500nm")
mid_cell_p_material.electron_minority_carrier_diffusion_length = si("5um")

all_materials.append(mid_cell_n_material)
all_materials.append(mid_cell_p_material)

# BOTTOM CELL - Ge
# Idem
Ge = material("Ge")
bot_cell_n_material = Ge(Nd=siUnits(2e18, "cm-3"), role="n")
bot_cell_p_material = Ge(Na=siUnits(1e17, "cm-3"), role="p")
Ge_alpha = np.loadtxt(this_dir_file("Ge-Palik.csv"), unpack=True, delimiter=",")
Ge.alphaE = interp1d(x=siUnits(Ge_alpha[0], 'eV'), y=Ge_alpha[1])
bot_cell_n_material.hole_minority_carrier_diffusion_length = si("800nm")
bot_cell_p_material.electron_minority_carrier_diffusion_length = si("50um")

all_materials.append(bot_cell_n_material)
all_materials.append(bot_cell_p_material)

# We add some other properties to the materials, assumed the same in all cases.
# If different, we should have added them above.
for mat in all_materials:
    mat.hole_mobility = 5e-2
    mat.electron_mobility = 3.4e-3
    mat.hole_mobility = 3.4e-3
    mat.electron_mobility = 5e-2
    mat.dielectric_constant = 9 * vacuum_permittivity

# And, finally, we put everithing together, adding also the surface recombination velocities.
triplejunction = Structure(
    [
        Layer(material=window_material, width=si("25nm")),
        Junction(
            (
                Layer(si("100nm"), material=top_cell_n_material),
                Layer(si("600nm"), material=top_cell_p_material),
            ),
            sn=1,
            sp=1,
        ),
        Junction(
            (
                Layer(si("100nm"), material=mid_cell_n_material),
                Layer(si("3.5um"), material=mid_cell_p_material),
            ),
            sn=1,
            sp=1,
        ),
        Junction(
            (
                Layer(si("400nm"), material=bot_cell_n_material),
                Layer(si("100um"), material=bot_cell_p_material),
            ),
            sn=1,
            sp=1,
        ),
    ]
)

# Now the 'solar_cell' object below contains all the information regarding the solar cell. So, we use it
solar_cell = State()
solar_cell.structure = triplejunction


class TestASC(TestCase):
    def test_81_quantum_efficiency(self):
        # Now we define the energies in which we want to calculate the QE.
        energy_bounds = (0.6, 3.4, 1000)
        E = si(np.linspace(*energy_bounds), "eV")

        # We add this energy to the solar cell structure.
        solar_cell.energies = E

        # We run the simulation
        qe_result = ASC.spectral_response_all_junctions(solar_cell, verbose=False)

        qe_top = qe_result["junctions"][0]["qe_tot"]
        qe_mid = qe_result["junctions"][1]["qe_tot"]
        qe_bot = qe_result["junctions"][2]["qe_tot"]

        top = np.max(qe_top)
        mid = np.max(qe_mid)
        bot = np.max(qe_bot)

        expected_top = 0.951927820731
        expected_mid = 0.976962091259
        expected_bot = 0.964448331241

        self.assertAlmostEqual(top, expected_top)
        self.assertAlmostEqual(mid, expected_mid)
        self.assertAlmostEqual(bot, expected_bot)

    def test_82_short_circuit_currents(self):
        # We create a solar spectrum using SPECTRAL2 and the default configuration
        # - See documentation and example of that package for more info
        spectrum = calculate_spectrum_spectral2()

        # We use the spectrum is in SI units: Watts m-2 joule-1.
        cell_area = 0.7 * 0.7 / 1e4  # expressed in m-2
        incident_x_J_y_per_J = spectrum["incident spectrum energy si"]
        incident_function = interp1d(y=incident_x_J_y_per_J[1] * cell_area, x=incident_x_J_y_per_J[0])

        # Now we define the energies in which we want to calculate the QE.
        energy_bounds = (0.6, 3.4, 1000)
        E = si(np.linspace(*energy_bounds), "eV")  # In SI units
        solar_cell.incident_light = E, incident_function(E) / E

        # We run the simulation
        qe_result = ASC.spectral_response_all_junctions(solar_cell, verbose=False)

        top = qe_result["junctions"][0]["J"]
        mid = qe_result["junctions"][1]["J"]
        bot = qe_result["junctions"][2]["J"]

        expected_top = 0.00496429316706
        expected_mid = 0.00597510188549
        expected_bot = 0.0136391295003

        self.assertAlmostEqual(top, expected_top)
        self.assertAlmostEqual(mid, expected_mid)
        self.assertAlmostEqual(bot, expected_bot)

    def test_83_solar_cell_parameters(self):
        # We create a solar spectrum using SPECTRAL2 and the default configuration
        # - See documentation and example of that package for more info
        spectrum = calculate_spectrum_spectral2()

        # We use the spectrum is in SI units: Watts m-2 joule-1.
        cell_area = 0.7 * 0.7 / 1e4  # expressed in m-2
        incident_x_J_y_per_J = spectrum["incident spectrum energy si"]
        incident_function = interp1d(y=incident_x_J_y_per_J[1] * cell_area, x=incident_x_J_y_per_J[0])
        power_density = spectrum["incident power density"]

        # Now we define the energies in which we want to calculate the QE.
        energy_bounds = (0.6, 3.4, 1000)
        E = si(np.linspace(*energy_bounds), "eV")  # In SI units
        solar_cell.incident_light = E, incident_function(E) / E

        # We run the simulation
        qe_result = ASC.spectral_response_all_junctions(solar_cell, verbose=False)

        top = qe_result["junctions"][0]["J"]
        mid = qe_result["junctions"][1]["J"]
        bot = qe_result["junctions"][2]["J"]

        cell_temp = 60  # ºC
        ref_temp = 25  # ºC

        Tcell = 273 + cell_temp
        Tref = 273 + ref_temp

        # The IV data will be stored in a State object. We create it, including the cell and reference temperatures.
        IV_calculation_state = State(T=Tcell, Tref=Tref)

        # From the QE object we get the short circuit currents
        Isc_array = (top, mid, bot)
        I01_array = [4.93e-24, 1.0e-21, 4.93e-6]
        I02_array = [3.28e-15, 2.7e-10, 1.0e-5]

        # This is the structure to calculate.
        IV_calculation_state.structure = Structure(
            [
                Junction(Eg=1.9, j01=I01_array[0], j02=I02_array[0], R_shunt=3e6,
                         R_series=0.0236, n1=1.00, n2=2.0, photocurrent=Isc_array[0]),
                Junction(Eg=1.4, j01=I01_array[1], j02=I02_array[1], R_shunt=1.5e6,
                         R_series=0.0012, n1=1.00, n2=2.0, photocurrent=Isc_array[1]),
                Junction(Eg=0.66, j01=I01_array[2], j02=I02_array[2], R_shunt=115,
                         R_series=8e-4, n1=1.00, n2=2.0, photocurrent=Isc_array[2]),
            ]
        )

        # We solve it, including explicitely the range of voltages we are interested
        IV_result = ASC.multijunctionIV(IV_calculation_state, V=np.linspace(0, 5, 10000))

        # Finally, we calculate the interesting solar cell parameters.
        I = IV_result["I"]
        V = IV_result["V"]
        P = np.max(I * V)

        Isc = IV_result["Isc"]
        Voc = IV_result["Voc"]
        FF = IV_result["FF"]
        eta = P / (power_density * cell_area)

        expected_Isc = 0.00496465142334
        expected_Voc = 2.57966686261
        expected_FF = 0.89393626837
        expected_eta = 0.261968183842

        self.assertAlmostEqual(Isc, expected_Isc)
        self.assertAlmostEqual(Voc, expected_Voc)
        self.assertAlmostEqual(FF, expected_FF)
        self.assertAlmostEqual(eta, expected_eta)

