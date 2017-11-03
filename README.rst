Solcore
=======

**Solcore** was born as a modular set of tools, written (almost) entirely in Python 3, to address some of the task we had to solve more often, such as fitting dark IV curves or luminescence decays. With time, however,  it has evolved as a complete semiconductor solver able of modelling the optical and electrical properties of a wide range of solar cells, from quantum well devices to multi-junction solar cells. Some of the features of Solcore are:

    - k•p band structure solver including strain
    - 1D arbitrary potential Schrödinger equation solver
    - Bulk and QW absorption profile calculator
    - Spectral irradiance model and database
    - Multi-junction quantum effciency and IV calculators
    - Coupled Poisson - Drift-Diffusion solver (PDD)

Installation
------------

After downloading Solcore, either using 'git' or as a zip file using one of the links on the right, installing it and become a user should be as easy as writing in the terminal::

    >>> python setup.py install

You will be asked to accept the license and then it will install all Solcore dependencies (except a Fortran compiler) and Solcore itself in the Python3 package tree. Using it will be as easy as making::

    >>> import solcore

        Welcome to Solcore - version 5.0.0
        Copyright (c) 2017, Imperial College, London All rights reserved.
        Software released under the GNU Lesser General Public License.


If you want to test first if Solcore will work in your computer, without actually installing it, or if you want to become a developer and therefore you need to have it in a more accessible place, you can test if Solcore works with::

    >>> python setup.py test

This will install the Solcore dependencies and run a few tests that probe several of the Solcore tools. If it fails, it will indicate which parts failed to work and why, and you could try to solve them. At the moment, this only cover some of Solcore's functionality, but it will be expanded with time.

Another thing that you might want to do before installing Solcore 5 is compiling the Poisson-drfit-diffusion solver. Assuming there is a Fortran compiler correctly configured to work with F2Py, compiling the library should be as easy as::

    >>> python setup.py build_pdd

This can also be done afterwards using the *config_tools* but you might need admin privileges, depending on where is the Python packages tree.

Things that can go wrong
^^^^^^^^^^^^^^^^^^^^^^^^

There are several things that can go wrong in the above description, specially in Windows.

1. **The tests associated with the Poisson-Drift-Diffusion solver fail**: This is usually the result of not having a Fortran compiler installed in your system, not being correctly configured or having a temperamental F2PY version, the tool - included in numpy - that makes Fotran code accesible from Python. For the first two problems, make sure you actually have a Fortran compiler installed and in the system path. For the latter, it appears in Windows and we have not been able to solve it, yet. Please, let us know if you have a solution.

2. **Some of the dependencies fail to install**: That is rarely the case, as all dependencies are in the main Python repositories. However, there might be issues (again in Windows) with Numpy, Matplotlib and Scipy. These packages need to be compiled and it is often easy to get them as a scientific bundle. You can check Anaconda <https://www.continuum.io/downloads> which provides all these packages together already configured for the correct OS.

Getting started
---------------

After installing Solcore (or even without installing it), there are a few things you might want to do in order to personalise it and start using it.

1. **Create a user configuration file:** This can be done automatically by importing the config_tools. If you do not already have a (hidden) solcore_config.txt file in your home directory, you will be asked if you want to create it::

    >>> import solcore.config_tools

2. **Create a folder with Solcore examples:** This is the fastest way of getting started. The examples will be created in a subfolder called 'solcore/examples'. These examples include all the scripts used in the main Solcore paper (https://arxiv.org/abs/1709.06741 ). Simply use the config_tools, again::

    >>> import solcore.config_tools as config

    >>> config.get_solcore_examples('/location/where/you/want/the/examples')

3. **Set the location of a SPICE executable and the SMARTS folder:** You will need to do this eventually in order to use those tools::

    >>> import solcore.config_tools as config

    >>> config.set_location_of_spice('/path/to/the/SPICE/executable')
    >>> config.set_location_of_smarts('/path/to/the/SMARTS/folder')

4. **Open Solcore documentation:** It should contain a description (even minimal) of all Solcore functions, modules and packages. The idea is for it to be a useful tool although it is quite empty, for now. The documentation will open in a web browser and you might want to add it to your Bookmarks::

    >>> import solcore.config_tools as config

    >>> config.open_documentation()

5. **Getting specific information about Solcore:** Even though the documentation "should" be more complete, you can get information about any object in python (including any Solcore function, module and package) using the '__doc__' attribute, for example::

    >>> import solcore.config_tools as config

    >>> print(config.get_current_config.__doc__)

    Prints the current Solcore configuration

        :return: None

6. **Python editor:** Learning Python is easy, but some tools make it even easier. That is the case of PyCharm <https://www.jetbrains.com/pycharm/> (the community eddition is free and the other it is too if you are in academia). Selecting an editor is very personal choice, but PyCharm turns out to be quite useful to teach you good coding practices, reviewing your code for errors and, in general, checking that things will work. It will make your life easier. Give it a try. Solcore in its current form is, in part, the result of using PyCharm.

Known issues
------------

We have developed Solcore as part of our ongoing research activities to solve specific challenges, it has (almost) never been a goal in itself. These means that there are parts of Solcore that might not be as polished as they should, that have been just partly implemented or that are only valid under some assumptions (good for us, but maybe not that good for others).

Some of the Solcore issues we are aware off are:

- The poisson-drift-diffusion solver, written in Fortran, has been tested only under Linux and Mac. We have never been successful in making F2Py and the Fortran compiler work together under Windows, although they are supposed to work well. Any help with this is more than welcome!!
- Documentation is incomplete or obscure, in many cases. Again, something to be solved soon.
- The calculator of the generation profile using the TMM module is really, really slow as soon as the structure is slightly complicated or the mesh density is high. We'll need to do something about it sooner than later.

Contributors
------------

Solcore is the result of many years of developement by many people, trying to put together in a consistent way the tools that the Quantum Photovoltaics Group at Imperial College London <https://www.imperial.ac.uk/quantum-photovoltaics/> needs for its research in solar cells and photovoltaics. The list of contributors (hopefully I am not forgetting anyone!) is:

- Diego Alonso-Álvarez
- Thomas Wilson
- Phoebe Pearce
- Markus Führer
- Daniel Farrel
- Tomos Thomas
- Alvin Chan
- Ned Ekins-Daukes
