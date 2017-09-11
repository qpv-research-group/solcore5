"""
In this example we calculate the bandstructure of a QW by solving the Schrodinger equation. We used two methods,
altough the interface is the same changing only the value of one input argument:

- Bulk-like solution (mode = 'kp8x8_bulk'): We calculate the band edges and effective masses as if they were bulk
materials and then the Schrodinger equation is solved using those as input. This method uses the last of the options
described in QM_example_1.

- 1D kp solver (mode = 'kp4x4'): In this case we solve the 4 bands kp hamiltonian in 1 dimension. The split-off band is
ignored.
"""
import solcore
from solcore import material, si
from solcore.structure import Layer
import solcore.quantum_mechanics as QM

GaAs = material("GaAs")
InGaAs = material("InGaAs")
GaAsP = material("GaAsP")
GaP = material("GaP")

InGaAs.strained = True
GaAsP.strained = True
GaP.strained = True
GaAs.strained = False

bulkMaterial = GaAs(T=300)
wellMaterial = InGaAs(In=0.245, T=300)
barrierMaterial = GaAsP(P=0.1, T=300)
interlayer_material = InGaAs(In=0.14, T=300)

my_structure = QM.assemble_qw_structure(
    repeats=1,
    well=Layer(si("7.2nm"), wellMaterial),
    bulk_l_top=Layer(si("1nm"), barrierMaterial),
    bulk_l_bottom=Layer(si("1nm"), barrierMaterial),
    barrier=Layer(si("28nm"), barrierMaterial),
    well_interlayer=Layer(si("3nm"), interlayer_material)
)
my_structure.substrate = bulkMaterial

# The argument graphtype='potentials' produces a plot with the bandstructure and the wavefunctions
band_edge, bands = QM.schrodinger(my_structure, mode='kp8x8_bulk', graphtype='potentials')

# We using the kp4x4 mode, the 'potentials' plot produces the band-edge potentials. If we want the actual 1D-energy
# bands, we can set bands=True. In the later plot, it can be seen that sometimes the band solver fails in finding the
# correct solution or it works but the bands are sorted incorrectly. In particular, band crossing is not always treated
# properly.
band_edge, bands = QM.schrodinger(my_structure, mode='kp4x4', graphtype='potentials', plot_bands=True)

# The band_edge dictionary contains the bands, effective masses and wavefunctions as a function of the position.
# The bands dictionary (in the case of kp4x4) contains everithing, including the bands versus k and the wavefunctions
# at each k.




