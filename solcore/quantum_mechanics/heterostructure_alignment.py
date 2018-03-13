def VBO_align(structure):
    """Takes in the the structure layers and sets Ec and Ev energies in the state which aligns the heterostructure."""

    for layer in structure:
        layer.Ev = layer.material.valence_band_offset
        layer.Ec = layer.material.valence_band_offset + layer.material.band_gap
    return structure
