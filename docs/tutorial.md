Solcore tutorial
================

This tutorial will guide you, step-by-step, in the creation of a solar cell and the calculation of their properties using *Solcore*.

**What do we have?**

- Dual junction InGaP/GaAs solar cell, lattice matched to GaAs. 
- The bottom cell has 30 strained-balanced quantum wells (QW), made of GaAsP/InGaAs.
- There is a distributed Brag reflector on the back of the cell, designed to reflect the radiation in the wavelength range of the QWs. 
- There is a tunnel junction in between the subcells.
- There is a dual layer anti-reflecting coating on the front made of MgF-ZnS

**What do we want to find out?**

- The absorption properties of our QWs
- The QE of the solar cell
- Its efficiency (Eff), short circuit current (Isc), open circuit voltage (Voc) and fill factor (FF) as a function of the concentration. 

Let's get started!

Defining the structure
----------------------

First we need to create the solar cell structure. It is made of several bits and pieces: the QWs, 2xJunctions, 1x tunnel junction and the ARC. We start with the first.

**Defining the QWs:**

```python
from solcore import material
from solcore.structure import Layer
import solcore.poisson_drift_diffusion as PDD

T = 300 

# First, we create the materials of the QW
QWmat = material('InGaAs')(T=T, In=0.2, strained=True)
Bmat = material('GaAsP')(T=T, P=0.1, strained=True)
i_GaAs = material('GaAs')(T=T)

# The QW is 7 nm wide, with GaAs interlayers 2 nm thick at each side and GaAsP barriers 10 nm thick.
# The final device will have 30 of these QWs.
QW = PDD.CreateDeviceStructure('QW', T=T, repeat=30, substrate=i_GaAs, layers=[
    Layer(width=10e-9, material=Bmat, role="barrier"),
    Layer(width=2e-9, material=i_GaAs, role="interlayer"),
    Layer(width=7e-9, material=QWmat, role="well"),
    Layer(width=2e-9, material=i_GaAs, role="interlayer"),
    Layer(width=10e-9, material=Bmat, role="barrier") ])

# We solve the quantum properties of the QW, leaving the default values of all parameters
QW_list = PDD.SolveQWproperties(QW)
```

The first few lines import the *Solcore* utilities needed to define the structure, including the *material* function, the *Layer* class and the Poisson-Drift-Diffusion solver. After that, we create the materials that will made the QWs and create a "Device structure". We will use these structure just to be able to solve the QWs afterwards and transform it into a sequence of layers with effective properties that the PDD solver understands. The call to "CreateDeviceStructure" has several inputs, including the temperature, the substrate, the number of repetition of the QWs and the structure of the layers. The call to "SolveQWproperties" will, indeed, use the utilities within the *quantum_mechanics* module to calculate the band structure of the QWs, their absorption coefficient and, finally, will calculate effective bandgap, density of states, etc. that the PDD solver will use. Although the device will have 30 quantum wells, only one unit (the one indicated in *Layers*) will be modelled as an isolated QW.  

If we only want to solve the properties of the QWs, without creating an effective structure for the PDD solver, we will use, instead:

```python
from solcore import material
from solcore.structure import Layer, Structure
import solcore.quantum_mechanics as QM

T = 300 

# First, we create the materials of the QW
QWmat = material('InGaAs')(T=T, In=0.2, strained=True)
Bmat = material('GaAsP')(T=T, P=0.1, strained=True)
i_GaAs = material('GaAs')(T=T)

# The QW is 7 nm wide, with GaAs interlayers 2 nm thick at each side and GaAsP barriers 10 nm thick.
QW = Structure([Layer(width=10e-9, material=Bmat, role="barrier"),
                Layer(width=2e-9, material=i_GaAs, role="interlayer"),
                Layer(width=7e-9, material=QWmat, role="well"),
                Layer(width=2e-9, material=i_GaAs, role="interlayer"),
                Layer(width=10e-9, material=Bmat, role="barrier") ], 
                substrate = i_GaAs)

# Finally, the quantum properties are calculated here
output = QM.schrodinger(QW)                     
```

While *Solcore* can solve the Schrödinger equation in a structure with any number of layers, the absorption calculator for QWs can only deal properly with single QWs. That is the reason of modelling only 1 QW despite having 30 in the structure. This will clearly represent a limitation when modelling the absorption of superlattices, where there is a strong coupling between neighbouring QWs. 

In the code above, we have used the "PDD.SolveQWproperties" and "QM.schrodinger" functions with the default values, but they both can have a number of optional input parameters to define the number of confine states to calculate, the energy of quasiconfine states, electric field, boundary conditions, etc. Please, visit the documentation of those functions to find out all the available options. 
