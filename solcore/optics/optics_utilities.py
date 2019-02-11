import numpy as np


def calculate_absorption(tmm_out):
    """It returns the diff_absorption (function), and the integration of it
    all_absorbed (float)
    ::tmm_out (dict)
    The tmm_out dictionary must have the position values (np.array, [nm]) and the 
    absorption matrix [wl,z].
    """
    all_z = tmm_out['position'] * 1e-9
    all_abs = tmm_out['absorption'] / 1e-9

    def diff_absorption(z):
        idx = all_z.searchsorted(z)
        idx = np.where(idx <= len(all_z) - 2, idx, len(all_z) - 2)
        try:
            z1 = all_z[idx]
            z2 = all_z[idx + 1]

            f = (z - z1) / (z2 - z1)

            out = f * all_abs[:, idx] + (1 - f) * all_abs[:, idx + 1]

        except IndexError:
            out = all_abs[:, idx]

        return out

    all_absorbed = np.trapz(diff_absorption(all_z), all_z)

    return diff_absorption, all_absorbed

def absorbed(self, z):
    out = self.diff_absorption(self.offset + z) * (z < self.width)
    return out.T

