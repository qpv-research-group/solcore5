"""
This example shows several options for obtaining the band structure of an strained bulk material using the solution of
the 8 bands Luttinger-Kohn hamiltonian. The Hamiltonian is the one shown in the following paper, but removing the two
bands corresponding to nitrogen-related levels.

Stanko Tomic et al., Electronic structure of InyGa1yAs1xNxGaAs(N) quantum dots by ten-band k.p theory.
Phys. Rev. B 73, 125348 (2006)

The bands are solve in the same way in all cases, being the difference the way the effective masses are calculated and
how the output is presented.
"""
import solcore
import solcore.quantum_mechanics as QM
from solcore import asUnit
from solcore.constants import electron_mass as m0

# Material parameters
GaAs = solcore.material("GaAs")(T=300)
GaAsP = solcore.material("GaAsP")(P=0.3, T=300)
InGaAs = solcore.material("InGaAs")(In=0.2, T=300)

""" QM.kp_bands solve the 8x8 kp hamiltonian for all k-values in a specific direction of the reciprocal space,
given by the effective_mass_direction. Optionally, the bands can be fitted around k=0 (the gamma point) to have an
effective mass, affected by the strain.

While the bands are calculated at all k-points in the specified direction, the output is just the band-edges
(the energies at k=0) and effective masses, if requested.
"""
bands_GaAsP = QM.kp_bands(GaAsP, GaAs, graph=False, fit_effective_mass=True, effective_mass_direction="X", return_so=True)
bands_InGaAs = QM.kp_bands(InGaAs, GaAs, graph=False, fit_effective_mass=True, effective_mass_direction="X", return_so=True)

print('Using QM.kp_bands:\n')
print('\t\t\tGaAsP\tInGaAs')
print('Ec (eV)\t\t = {:.3f}\t{:.3f}'.format(asUnit(bands_GaAsP[0], 'eV'), asUnit(bands_InGaAs[0], 'eV')))
print('Ehh (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(bands_GaAsP[1], 'eV'), asUnit(bands_InGaAs[1], 'eV')))
print('Elh (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(bands_GaAsP[2], 'eV'), asUnit(bands_InGaAs[2], 'eV')))
print('Eso (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(bands_GaAsP[3], 'eV'), asUnit(bands_InGaAs[3], 'eV')))
print('me (m0)\t\t = {:.3f}\t{:.3f}'.format(bands_GaAsP[4]/m0, bands_InGaAs[4]/m0))
print('mhh (m0)\t = {:.3f}\t{:.3f}'.format(bands_GaAsP[5]/m0, bands_InGaAs[5]/m0))
print('mlh (m0)\t = {:.3f}\t{:.3f}'.format(bands_GaAsP[6]/m0, bands_InGaAs[6]/m0))
print('mso (m0)\t = {:.3f}\t{:.3f}'.format(bands_GaAsP[7]/m0, bands_InGaAs[7]/m0))
print()
print()

""" QM.KPbands is similar to the above but it has the option of returning the whole band (for all k-points) or only
the band edges at k=0. The direction of the k-points is defined by two angles in the reciprocal space (t and p) and you
(it defaults to the X direction) and you can choose to calculate just a certain fraction of that direction and not all
the k-points untile the boundary of the Brillouin zone.

The effective masses are not calculated and therefore a separate function has to be used for that. The routine for that
is slightly different, so the effective mases might also differ a bit from those of the QM.kp_bands.
"""

edges_GaAsP = QM.KPbands(GaAsP, GaAs, fraction=0.2, return_edges_only=True)
bands_GaAsP = QM.KPbands(GaAsP, GaAs, fraction=0.2)
masses_GaAsP = QM.fit_effective_masses(bands_GaAsP, GaAsP, GaAs, plot_result=False)

edges_InGaAs = QM.KPbands(InGaAs, GaAs, fraction=0.2, return_edges_only=True)
bands_InGaAs = QM.KPbands(InGaAs, GaAs, fraction=0.2)
masses_InGaAs = QM.fit_effective_masses(bands_InGaAs, InGaAs, GaAs, plot_result=False)

print('Using QM.KPbands and QM.fit_effective_masses:\n')
print('\t\t\tGaAsP\tInGaAs')
print('Ec (eV)\t\t = {:.3f}\t{:.3f}'.format(asUnit(edges_GaAsP[0], 'eV'), asUnit(edges_InGaAs[0], 'eV')))
print('Ehh (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(edges_GaAsP[1], 'eV'), asUnit(edges_InGaAs[1], 'eV')))
print('Elh (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(edges_GaAsP[2], 'eV'), asUnit(edges_InGaAs[2], 'eV')))
print('Eso (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(edges_GaAsP[3], 'eV'), asUnit(edges_InGaAs[3], 'eV')))
print('me (m0)\t\t = {:.3f}\t{:.3f}'.format(masses_GaAsP[0]/m0, masses_InGaAs[0]/m0))
print('mhh (m0)\t = {:.3f}\t{:.3f}'.format(masses_GaAsP[1]/m0, masses_InGaAs[1]/m0))
print('mlh (m0)\t = {:.3f}\t{:.3f}'.format(masses_GaAsP[2]/m0, masses_InGaAs[2]/m0))
print('mso (m0)\t = {:.3f}\t{:.3f}'.format(masses_GaAsP[3]/m0, masses_InGaAs[3]/m0))
print()
print()

"""
The last option is QM.kp8x8_bulk. This one uses KPbands to calculate the band edges and then calculates the effective
masses by averaging the effective mass in all possible k-directions in the gamma-X-L plane. Yep, let's call this
an "effective" effective mass.

This is the method used by default when solving QWs in a bulk-like way
"""
bands_GaAsP = QM.kp8x8_bulk(GaAsP, GaAs)
bands_InGaAs = QM.kp8x8_bulk(InGaAs, GaAs)

print('Using QM.kp8x8_bulk:\n')
print('\t\t\tGaAsP\tInGaAs')
print('Ec (eV)\t\t = {:.3f}\t{:.3f}'.format(asUnit(bands_GaAsP[0], 'eV'), asUnit(bands_InGaAs[0], 'eV')))
print('Ehh (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(bands_GaAsP[1], 'eV'), asUnit(bands_InGaAs[1], 'eV')))
print('Elh (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(bands_GaAsP[2], 'eV'), asUnit(bands_InGaAs[2], 'eV')))
print('Eso (eV)\t = {:.3f}\t{:.3f}'.format(asUnit(bands_GaAsP[3], 'eV'), asUnit(bands_InGaAs[3], 'eV')))
print('me (m0)\t\t = {:.3f}\t{:.3f}'.format(bands_GaAsP[4]/m0, bands_InGaAs[4]/m0))
print('mhh (m0)\t = {:.3f}\t{:.3f}'.format(bands_GaAsP[5]/m0, bands_InGaAs[5]/m0))
print('mlh (m0)\t = {:.3f}\t{:.3f}'.format(bands_GaAsP[6]/m0, bands_InGaAs[6]/m0))
print('mso (m0)\t = {:.3f}\t{:.3f}'.format(bands_GaAsP[7]/m0, bands_InGaAs[7]/m0))
print()
print()