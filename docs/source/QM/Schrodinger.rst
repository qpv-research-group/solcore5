Quantum Solvers
===============

- Example 1: :doc:`KP calculator and effective mass fitting <../Examples/example_K_and_effective_mass>`
- Example 2: :doc:`Calculating the bandstructure and local density of states in QWs and MQWs <../Examples/example_QWs>`

﻿The electronic band structure of semiconductor materials is responsible for their light absorption and emission properties as well as for many of their transport properties, ultimately depending on the carriers' effective masses. These properties are not intrinsic to the material, but depend on external factors, too, most notably the strain and the quantum confinement.

Given the crystalline nature of most semiconductor materials, there will be strain whenever two materials with different crystal lattice constants are grown on top of each other pseudomorphically. Even those with the same lattice constant might be under stress due to other effects such as the presence of impurities or if used at different temperatures having dissimilar thermal expansion coefficients. Quantum confinement, in turn, takes place when the size of the semiconductor material in one or more dimensions is reduced to a few nanometres. In that situation, the energy levels available to the carriers become quantized in the direction of confinement, also changing the density of states. Both conditions take place simultaneously when dealing with strain-balanced quantum wells (QW).

Quantum wells - and more recently quantum wires - have been employed to tune the absorption properties of high efficiency solar cells for the past two decades. The need for appropriate tools to study them in the context of photovoltaics led to the development of the simulation models that were the seed of Solcore. As strained materials with quantum confinement, special care must be taken to obtain a sensible set of parameters for the QW structures, including the band edges with confined energy levels, the effective masses and the absorption coefficient.

Solcore's approach for evaluating the properties of QWs involves calculating first the effect of strain using a 8-band Pikus-Bir Hamiltonian , treating each material in the structure as bulk, and then using the shifted bands and effective masses to solve the 1D Schödinger equation, after a proper alignment between layers. Finally, the absorption coefficient is calculated based on the 2D density of states, including the effect of excitons.

This section describes the tools used for the calculation of the band structure. The QW absorption profile is discussed in :doc:`the optics section section <../Optics/material_optics>`.

Bulk 8-band kp solver
---------------------

.. automodule:: solcore.quantum_mechanics.kp_bulk
    :members:
    :undoc-members:

Quantum well 1D-Schrodinger calculator
--------------------------------------

.. automodule:: solcore.quantum_mechanics.high_level_kp_QW
    :members:
    :undoc-members:

.. automodule:: solcore.quantum_mechanics.structure_utilities
    :members:
    :undoc-members:

.. automodule:: solcore.quantum_mechanics.potential_utilities
    :members:
    :undoc-members:

.. automodule:: solcore.quantum_mechanics.graphics
    :members:
    :undoc-members:

4-bands kp QW solver (experimental)
-----------------------------------

.. automodule:: solcore.quantum_mechanics.kp_QW
    :members:
    :undoc-members: