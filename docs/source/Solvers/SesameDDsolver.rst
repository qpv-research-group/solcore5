Sesame Poisson Drift-Diffusion solver
======================================

This solver provides an interface to `Sesame <https://sesame.readthedocs.io/en/latest/>`_ to solve the
Drift-Diffusion equations. While it serves the same purpose as the legacy Fortran PDD solver, this solver has
the advantage of being written entirely in Python, thus being more transparent to most users.

Unlike the Fortran-based solver, this solver can handle both constant values of doping per layer,
or depth-dependent doping profiles defined by the user.

.. automodule:: solcore.sesame_drift_diffusion.solve_pdd
    :members:
    :undoc-members: