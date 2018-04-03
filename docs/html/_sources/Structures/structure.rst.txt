Structures and support classes
==============================

While Solcore is mostly about physics, it also needs a lot of tools that keep everiting working together properly or that support the creation and management of solar cell structure. Here you can find (some) information about all those extra bits and pieces Solcore is made of.

.. toctree::
    :maxdepth: 2

    structure

Structure
---------

Solcore calculates optical and electrical properties of solar cells or, in other words, of a certain combination of layers made of different materials and serving to specific purposes. The *structure* module contains tha basic building blocks that allow you to make a solar cell structure and calculate their properties. Their use will be explained with examples in other packages - such as in the quantum efficiency calculator or the Poisson-drift-diffusion solver.

.. automodule:: solcore.structure
    :members:
    :undoc-members:

Another, very useful class is the *state_object* which makes the use of dictionaries very easy.

.. automodule:: solcore.state
    :members:
    :undoc-members:

Solar Cells
-----------

Finally, the higher level building block of solar cells is contained in the *solar_cell* module.

.. automodule:: solcore.solar_cell
    :members:
    :undoc-members:

Science tracker
---------------

Solcore is an original work, but the equations it implements and the data it uses have been often published elsewere. The Science tracker allows you to keep track of those references and check yourself their origin and assumpions.

.. automodule:: solcore.science_tracker
    :members:
    :undoc-members:

Other classes
-------------

Solcore depends on many other classes and functions that work in the background and keep everything together. Some of them are used directly, others are completely transparent but make Solcore work they way it does. Unfortunately, their documentation is not very extensive or explicit, so feel free to explore the source code to understand what they are doing and how.

.. automodule:: solcore.interpolate
    :members:
    :undoc-members:

.. automodule:: solcore.crystals
    :members:
    :undoc-members:

.. automodule:: solcore.smooth
    :members:
    :undoc-members:

.. automodule:: solcore.singleton
    :members:
    :undoc-members:

.. automodule:: solcore.source_managed_class
    :members:
    :undoc-members:

