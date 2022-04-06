import numpy as np
import matplotlib.pyplot as plt

from solcore import material, si

from solcore.solar_cell import SolarCell, Layer, Junction
from solcore.solar_cell_solver import solar_cell_solver
from solcore.state import State

GaInP = material("GaInP")(In=0.5)
GaAs = material("GaAs")()
Ge = material("Ge")()

optical_struct = SolarCell(
    [
        Layer(material=GaInP, width=si("5000nm")),
        Junction(
            [
                Layer(material=GaAs, width=si("200nm")),
                Layer(material=GaAs, width=si("5um")),
            ],
            kind="DA",
        ),
        Layer(material=Ge, width=si("50um")),
    ]
)

wl = np.linspace(250, 1700, 300) * 1e-9

options = State()
options.wavelength = wl
options.optics_method = "TMM"
options.no_back_reflection = False
options.BL_correction = True
options.recalculate_absorption = True
options.positions = [1e-8, 1e-9, 1e-8, 1e-7]
options.theta = 0


c_list = [
    ["c", "c", "c", "c"],
    ["c", "c", "c", "i"],
    ["c", "i", "i", "c"],
    ["i", "i", "i", "i"],
]

titles = [
    "All coherent",
    "Bottom Ge layer explicity incoherent",
    "Both layers of GaAs junction incoherent",
    "All layers incoherent",
]

for i1, cl in enumerate(c_list):
    plt.figure(i1)
    options.coherency_list = cl
    solar_cell_solver(optical_struct, "optics", options)
    plt.plot(wl * 1e9, optical_struct[0].layer_absorption)
    plt.plot(wl * 1e9, optical_struct[1].layer_absorption)
    plt.plot(wl * 1e9, optical_struct[2].layer_absorption)
    plt.plot(wl * 1e9, optical_struct.reflected, "--")
    plt.plot(wl * 1e9, optical_struct.transmitted, "--")
    plt.plot(
        wl * 1e9,
        optical_struct[0].layer_absorption
        + optical_struct[1].layer_absorption
        + optical_struct[2].layer_absorption
        + optical_struct.reflected
        + optical_struct.transmitted,
    )
    plt.legend(["GaInP", "GaAs", "Ge", "R", "T", "R+A+T"], loc="upper left")
    plt.title(titles[i1])
    plt.xlabel("Wavelength (nm)")
    plt.show()

# plt.figure(1)
# plt.plot(np.linspace(0, 3000, 100),
#         optical_struct[1].absorbed(np.linspace(0, 3000, 100)*1e-9)[:,0])
# plt.show()
