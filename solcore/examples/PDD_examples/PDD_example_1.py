# Example 2 - Device creation
# ---------------------------
# This example illustrates the creation of a QW structure with 40 wells, solve its quantum properties and
# add it to the intrinsic region of a p-i-n solar cell

import solcore.poisson_drift_diffusion as PDD
from solcore.structure import Layer
from solcore import material

# First, we create the materials of the QW
QWmat = material('InGaAs')(T=300, In=0.2)
Bmat = material('GaAsP')(T=300, P=0.1)
i_GaAs = material('GaAs')(T=300)

# The QW is 7 nm wide, with GaAs interlayers 2 nm thick at each side and GaAsP barriers 10 nm thick.
# The final device will have 40 of these QWs.
QW = PDD.CreateDeviceStructure('QW', T=300, repeat=40, layers=[
    Layer(width=10e-9, material=Bmat, role="barrier"),
    Layer(width=2e-9, material=i_GaAs, role="interlayer"),
    Layer(width=7e-9, material=QWmat, role="well"),
    Layer(width=2e-9, material=i_GaAs, role="interlayer"),
    Layer(width=10e-9, material=Bmat, role="barrier"),
])

# We solve the quantum properties of the QW, leaving the default values of all parameters
PDD.SolveQWproperties(QW)

# We create the other materials we need for the device
window = material('AlGaAs')(T=300, Na=1e24, Al=0.8)
p_GaAs = material('GaAs')(T=300, Na=1e24)
n_GaAs = material('GaAs')(T=300, Nd=1e23)
bsf = material('AlGaAs')(T=300, Nd=1e24, Al=0.4)

# And finally we create another p-i-n structure incorporating the QWs in the intrinsic region.
MyDevice = PDD.CreateDeviceStructure('TestDevice', T=300, layers=[
    Layer(width=30e-9, material=window, role="Window"),
    Layer(width=400e-9, material=p_GaAs, role="Emitter"),
    Layer(width=10e-9, material=i_GaAs, role="Intrinsic"),
    QW,
    Layer(width=10e-9, material=i_GaAs, role="Intrinsic"),
    Layer(width=2000e-9, material=n_GaAs, role="Base"),
    Layer(width=200e-9, material=bsf, role="BSF")
])

# Example 3 - Device creation
# ---------------------------
# Following the previous example, we can save the properties of the QW and the final device.
# We choose to keep the absorption coeficients in separate files but leave the default directory.
# The directory containing the absorption coefficients will be 'myQW_inputs' and the name of the files will have the following pattern:
#
#       <deviceName>_<layerNumber>_<layerComposition>_<layerRole>.dat
#

import os

# First we load the device structure stored in MyDevice.json
this_dir = os.path.split(__file__)[0]
file = 'MyDevice'
QW_file = 'myQW'
full_path = os.path.join(this_dir, file)
full_path_QW = os.path.join(this_dir, QW_file)

# The absorption of the well itself will be: 'QW_2_In0.2GaAs_well.dat'
PDD.Save(QW, full_path_QW, save_absorptions_individually=True, remove_absorption_from_json=True)

# And the same for the complete device:
PDD.Save(MyDevice, full_path, save_absorptions_individually=True, remove_absorption_from_json=True)
