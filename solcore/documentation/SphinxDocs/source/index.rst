Solcore
=======

:literal:`Solcore` was born as a modular set of tools, written (almost) entirely in Python 3, to address some of the task we had to solve more often, such as fitting dark IV curves or luminescence decays. With time, however,  it has evolved as a complete semiconductor solver able of modelling the optical and electrical properties of a wide range of solar cells, from quantum well devices to multi-junction solar cells. Some of the features of Solcore are:

- :doc:`k•p band structure solver including strain <QM/Schrodinger>`
- :doc:`1D arbitrary potential Schrödinger equation solver <QM/Schrodinger>`
- :doc:`Bulk and QW absorption profile calculator <Absorption/Absorption>`
- :doc:`Spectral irradiance model and database <spectral/spectral>`
- :doc:`Multi-junction quantum effciency and IV calculators <ASC/ASC>`
- :doc:`Coupled Poisson - Drift-Diffusion solver (PDD) <PDD/DDsolver>`


Contents:
---------

.. toctree::
    :maxdepth: 2

    installation
    structure
    Systems/Units
    Systems/Parameters
    Systems/Materials
    QM/Schrodinger
    Absorption/Absorption
    spectral/spectral
    ASC/ASC
    PDD/DDsolver


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

