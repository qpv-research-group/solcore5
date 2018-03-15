Solcore tutorial
================

This tutorial will guide you, step-by-step, in the creation of a solar cell and the calculation of their properties using *Solcore*.

**What do we have?**

- Dual junction GaInP/GaAs solar cell, lattice matched to GaAs. 
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

In the code above, we have used the "PDD.SolveQWproperties" and "QM.schrodinger" functions with the default values, but they both can have a number of optional input parameters to define the number of confined states to calculate, the energy of quasiconfined states, electric field, boundary conditions, etc. Please, visit the documentation of those functions to find out all the available options.

**Defining the junctions:**

In order to calculate the properties of a solar junction using the PDD solver, we need to give all the layers and materials the junciton is made of, in a similar way we have done for the QWs. One thing to note is that if *Solcore* cannot find a property it needs to solve the PDD equations, *it will take the corresponding property for GaAs as a default value*. So, be sure you provide all the required values or that you are happy with the defaults. 

***NOTE***: The different code snippets are additive in order to get a final, complete script. Normally, all the "import" statements would be packed together at the beginning. 

```python
from solcore.structure import Junction

T = 300 

## Materials for the BOTTOM junction
window_bottom = material('GaInP')(T=T, Nd=5e24, In=0.49)
n_GaAs = material('GaAs')(T=T, Nd=1e24)
p_GaAs = material('GaAs')(T=T, Na=8e22)
bsf_bottom = material('GaInP')(T=T, Na=5e24, In=0.49)

GaAs_junction = Junction([Layer(width=10e-9, material=window_bottom, role="Window"),
                   Layer(width=150e-9, material=n_GaAs, role="Emitter")] + 
                   QW_list + 
                   [Layer(width=3000e-9, material=p_GaAs, role="Base"),
                   Layer(width=200e-9, material=bsf_bottom, role="BSF")], sn=1e6, sp=1e6, T=T, kind='PDD')
                   
## Materials for the TOP junction
window_top = material('AlGaAs')(T=T, Nd=5e24, Al=0.8)
n_GaInP = material('GaInP')(T=T, Nd=5e24, In=0.49)
p_GaInP = material('GaInP')(T=T, Na=8e22, In=0.49)
bsf_top = material('AlGaAs')(T=T, Na=5e24, Al=0.8)

GaInP_junction = Junction([Layer(width=40e-9, material=window_top, role="Window"),
                   Layer(width=120e-9, material=p_GaAs, role="Emitter"),
                   Layer(width=400e-9, material=n_GaAs, role="Base"),
                   Layer(width=30e-9, material=bsf_top, role="BSF")], sn=1e6, sp=1e6, T=T, kind='PDD')
```

