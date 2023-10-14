# from sesame import Builder, IVcurve, Analyzer
import numpy as np
from solcore.constants import q

from solcore.registries import (
    register_iv_solver,
    register_equilibrium_solver
)

try:
    import sesame

    reason_to_exclude = None
except ImportError:
    reason_to_exclude = (
        "Sesame was not installed."
    )

def process_structure(junction):

    # get material parameters from the junction and convert them to format required by Sesame (dictionary)
    # convert from Solcore units (i.e. base SI units like m) to Sesame's units (cm, eV).
    # Note that internally in sesame, many quantities are scaled/dimensionless.

    material_list = []
    for layer in junction:

        bulk_recombination_energy = layer.mat.bulk_recombination_energy if \
            hasattr(layer.mat, "bulk_recombination_energy") else 0

        new_mat = {
            'Nc': layer.mat.Nc * 1e-6, # effective density of states at CB edge (cm-3)
            'Nv': layer.mat.Nv * 1e-6, # effective density of states at VB edge (cm-3)
            'Eg': layer.mat.band_gap / q, # material bandgap (eV)
            'affinity': layer.mat.electron_affinity / q, # electron affinity (eV)
            'epsilon': layer.mat.relative_permittivity, # relative permittivity (dimensionless)
            'mu_e': layer.mat.electron_mobility*1e4, # electron mobility (cm2/(V s))
            'mu_h': layer.mat.hole_mobility*1e4, # hole mobility (cm2/(V s))
            'tau_e': layer.mat.electron_minority_lifetime, # electron bulk lifetime (s)
            'tau_h': layer.mat.hole_minority_lifetime, # hole bulk lifetime (s)
            'Et': bulk_recombination_energy, # energy level of bulk recombination centres (eV)
            'B': layer.mat.radiative_recombination*1e6, # radiative recombination constant (cm3/s)
            'Cn': layer.mat.electron_auger_recombination*1e12, # Auger recombination constant for electrons (cm6/s)
            'Cp': layer.mat.hole_auger_recombination*1e12, # Auger recombination constant for holes (cm6/s)
        }
        # Note: all of these can be functions of position; should eventually give users a way to use this
        # functionality

    return material_list

@register_equilibrium_solver("sesame_PDD", reason_to_exclude=reason_to_exclude)
def equilibrium():
    pass

@register_iv_solver("sesame_PDD", reason_to_exclude=reason_to_exclude)
def iv():
    pass

def eqe():
    pass