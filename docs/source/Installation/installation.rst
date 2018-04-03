Installation and configuration
==============================

Installing Solcore
^^^^^^^^^^^^^^^^^^

Written mostly in Python, Solcore should work in all OS. The following instructions are expected to work in more systems, but check the sections below for OS specific instructions and caveats.

In order to install Solcore in your computer, you will need the following:

- Python >3.4
- pip
- setuptools

*Pip* is normally part of the standard Python 3 installation, but you might need to install *setuptools* manually with (pip3 only indicates that is the *pip* that comes with Python 3)::

    pip3 install setuptools

After downloading Solcore, either using 'git' or as a zip file using one of the links on the right, installing it and become a user should be as easy as writing in the terminal (within the Solcore5 folder, where there is the *setup.py* file)::

    python3 setup.py install

You will be asked to accept the license and then it will install all Solcore dependencies (except a Fortran compiler) and Solcore itself in the Python3 package tree. Using it will be as easy as making::

    >>> import solcore

        Welcome to Solcore - version 5.2.0
        Copyright (c) 2017, Imperial College, London All rights reserved.
        Software released under the GNU Lesser General Public License.


If you want to test first if Solcore will work in your computer, without actually installing it, or if you want to become a developer and therefore you need to have it in a more accessible place, you can test if Solcore works with::

    python3 setup.py test

This will also install the Solcore dependencies and run a few tests that probe several of the Solcore tools. If it fails, it will indicate which parts failed to work and why, and you could try to solve them. At the moment, this only cover some of Solcore's functionality, but it will be expanded with time.

Another thing that you might want to do before installing Solcore 5 is compiling the Poisson-drfit-diffusion solver. Assuming there is a Fortran compiler correctly configured to work with F2Py, compiling the library should be as easy as::

    >>> python3 setup.py build_pdd

This can also be done afterwards using the *config_tools* (see below) but you might need admin privileges, depending on where the Python packages tree is located. Solcore will try to compile the fortran code automatically the first time that any of the modules using it are imported. During the compilation, you will see many warnings printed. Don't worry, it is normal and, hopefully, the compilation should complete without any critical error.


Getting started
^^^^^^^^^^^^^^^

After installing Solcore (or even without installing it), there are a few things you might want to do in order to personalise it and start using it. These are general instructions that, hopefully, should work in most *NIX systems. Check specific instructions for:

- :doc:`Windows 10 <Solcore_on_Windows>`
- :doc:`Mac OS X <Solcore_on_MacOSX>`

1. **Create a user configuration file:** This can be done automatically by importing the config_tools. If you do not already have a (hidden) solcore_config.txt file in your home directory, you will be asked if you want to create it::

    >>> import solcore.config_tools

2. **Check Solcore examples:** This is the fastest way of getting started. These examples include all the scripts used in the main Solcore paper (submitted to Journal of Computational Electronics, preprint in https://arxiv.org/abs/1709.06741 ) and a few others to explore other functionality. We hope to increase the number and usefulness of these examples over time.

3. **Check the Tutorial:** The :doc:`Tutorial <../tutorial>` (with the code available in the Examples folder) contains a step-by-step example on how to create a solar cell with multiple junctions and quantum wells, and calculate its quantum efficiency and current-voltage characteristics. While it does not cover every option - not even close! - it should give you a good idea on how things are done within Solcore and its possibilities.

3. **Set the location of a SPICE executable and the SMARTS folder:** You will need to do this eventually in order to use those tools::

    >>> import solcore.config_tools as config

    >>> config.set_location_of_spice('/path/to/the/SPICE/executable')
    >>> config.set_location_of_smarts('/path/to/the/SMARTS/folder')

4. **Installing S4:** The rigorous-coupled wave analysis (RCWA) solver S4, created by the Stanford University, needs to be installed separately using `the fork by Phoebe Pearce <http://github.com/phoebe-p/S4>`_. The original version **do not work** with Python 3.x, only with Python 2. You can find more information about the solver in the `project webpage <http://web.stanford.edu/group/fan/S4/>`_. An example of its use with Solcore is included in the examples folder, *Figure9.py*.

5. **Getting specific information about Solcore:** Even though the documentation "should" be more complete, you can get information about any object in python (including any Solcore function, module and package) using the '__doc__' attribute, for example::

    >>> import solcore.config_tools as config

    >>> print(config.get_current_config.__doc__)

    Prints the current Solcore configuration

        :return: None

6. **Python editor:** Learning Python is easy, but some tools make it even easier. That is the case of PyCharm <https://www.jetbrains.com/pycharm/> (the community eddition is free and the other it is too if you are in academia). Selecting an editor is very personal choice, but PyCharm turns out to be quite useful to teach you good coding practices, reviewing your code for errors and, in general, checking that things will work. It will make your life easier. Give it a try. Solcore in its current form is, in part, the result of using PyCharm.

Problems with the installation
------------------------------

There are several things that can go wrong in the above description, specially in Windows.

1. **The tests associated with the Poisson-Drift-Diffusion solver fail**: This is usually the result of not having a Fortran compiler installed in your system, not being correctly configured or having a temperamental F2PY version, the tool - included in numpy - that makes Fotran code accesible from Python. For the first two problems, make sure you actually have a Fortran compiler installed and in the system path. For the latter, it appears in Windows and we have not been able to solve it, yet. Please, let us know if you have a solution.

2. **Some of the dependencies fail to install**: That is rarely the case, as all dependencies are in the main Python repositories. However, there might be issues (again in Windows) with Numpy, Matplotlib and Scipy. These packages need to be compiled and it is often easy to get them as a scientific bundle. You can check Anaconda <https://www.continuum.io/downloads> which provides all these packages together already configured for the correct OS.

Other Known issues
------------------

We have developed Solcore as part of our ongoing research activities to solve specific challenges, it has (almost) never been a goal in itself. These means that there are parts of Solcore that might not be as polished as they should, that have been just partly implemented or that are only valid under some assumptions (good for us, but maybe not that good for others).

Some of the Solcore issues we are aware off are:

- The poisson-drift-diffusion solver, written in Fortran, has been tested only under Linux and Mac. We have never been successful in making F2Py and the Fortran compiler work together under Windows, although they are supposed to work well. Any help with this is more than welcome!!
- Documentation is incomplete or obscure, in many cases. Again, something to be solved soon.


The config_tools
^^^^^^^^^^^^^^^^

This module contains all the functions that will help you to setup and configure your Solcore installation, as it has been highlighted above. The full description of the funcitons included in this module are:

.. automodule:: solcore.config_tools
    :members:
    :undoc-members:
