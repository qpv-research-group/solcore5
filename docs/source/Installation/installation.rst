Installation and configuration
==============================

Installing Solcore
^^^^^^^^^^^^^^^^^^

Solcore is written mostly in Python, but the Poisson-Drift-diffusion (PDD) solver is written in Fortran to make it more efficient. The following instructions are expected to work in most systems, but check the sections below for OS specific instructions and caveats. In order to install Solcore in your computer, you will need the following:

- Python >3.7
- pip
- setuptools
- numpy
- :doc:`a suitable Fortran compiler (only needed for the PDD solver) <compilation>`

*Pip* is normally part of the standard Python installation, but you might need to install *setuptools* and *numpy* manually with::

    pip install setuptools numpy

Now, installing Solcore should be as easy as writing in the terminal::

    pip install solcore

This will download Solcore form the Pypi repository and install the package within the Python packages tree, but it will *NOT* install the PDD solver, for which you need a suitable Fortran compiler :doc:`(read this to install your Fortran compiler) <compilation>`. Assuming you have a Fortran compiler correctly installed and configured, you can install Solcore with the PDD solver by doing::

    pip install solcore
	pip install --no-deps --force-reinstall --upgrade --install-option="--with_pdd" solcore

**NOTE**: Pip passes the --install-option to all the dependencies of the package and therefore the installation will fail since those dependencies do not have the option "--with_pdd". That is the reason why, for now, Solcore needs to be installed twice: the first one install Solcore normally with all the dependencies and the second one re-installs Solcore with PDD support. 

And that's all!! Solcore should be available to be used as with any Python package::

    >>> import solcore

        Welcome to Solcore - version 5.6.0
        Copyright (c) 2018, Imperial College, London All rights reserved.
        Software released under the GNU Lesser General Public License.

Alternative installation method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, you can `download the source from the Solcore GitHub repository <https://github.com/dalonsoa/solcore5>`_, either using 'git' or as a zip file using one of the links on the right. If you want to install it, unpack it and run from the directory where *setup.py* is located (you still need *pip*, *setuptools* and *numpy* installed; see above)::

    pip install .

or::

    pip install . --install-option="--with_pdd"

If you want to test first if Solcore will work in your computer, without actually installing it, or if you want to become a developer and therefore you need to have it in a more accessible place, you can test if Solcore works with::

    python setup.py test

This will also install the Solcore dependencies and run a few tests that probe several of the Solcore tools. If it fails, it will indicate which parts failed to work and why, and you could try to solve them. At the moment, this only cover some of Solcore's functionality, but it will be expanded with time. The tests related to the PDD solver will fail, of course, since the PDD solver will not be compiled. If you want to test everything, compiling the the PDD solver, just run::

    python setup.py test --with_pdd

Install in development mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are planning to develop Solcore further, you would want to have all the files in an accessible place but still being able to use the package from other places and examples, behaving as if it were truly installed. For this, install Solcore in development mode. `Download the source from the Solcore GitHub repository <https://github.com/dalonsoa/solcore5>`_ as above and then::

    pip install -e .[dev]

or::

    pip install -e .[dev] --install-option="--with_pdd"
    
Solcore uses *pre-commit* to do a few things before commiting the changes (for example, clearing the output of Jupyter Notebooks). The *pre-commit* package is installed automatically with the above commands, but each user needs to be initialise it before it can work. This can be done with::

    pre-commit install
   
Check the `pre-commit webpage <https://pre-commit.com/#3-install-the-git-hook-scripts>`_ for more information on how it works. 

Getting started
^^^^^^^^^^^^^^^

After installing Solcore (or even without installing it), there are a few things you might want to do in order to personalise it and start using it. These are general instructions that, hopefully, should work in most *NIX systems. Check specific instructions for:

- :doc:`Windows 10 <Solcore_on_Windows>`
- :doc:`Mac OS X <Solcore_on_MacOSX>`

1. **Create a user configuration file:** This can be done automatically by importing the config_tools. If you do not already have a (hidden) solcore_config.txt file in your home directory, you will be asked if you want to create it::

    >>> import solcore.config_tools

2. **Check Solcore examples:** This is the fastest way of getting started. These examples include all the scripts used in the main Solcore paper in the `Journal of Computational Electronics (2018) <https://doi.org/10.1007/s10825-018-1171-3>`_ and a few others to explore other functionality. We hope to increase the number and usefulness of these examples over time.

3. **Check the Tutorial:** The :doc:`Tutorial <../Examples/tutorial>` (with the code available in the Examples folder) contains a step-by-step example on how to create a solar cell with multiple junctions and quantum wells, and calculate its quantum efficiency and current-voltage characteristics. While it does not cover every option - not even close! - it should give you a good idea on how things are done within Solcore and its possibilities.

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

1. **The tests associated with the Poisson-Drift-Diffusion solver fail**: This is usually the result of not having a Fortran compiler installed in your system, not being correctly configured or having a temperamental F2PY version, the tool - included in numpy - that makes Fotran code accesible from Python. Please, make sure you follow all the steps indicated in the :doc:`Fortran compiler section <compilation>` and above to have the PDD solver up and running.

2. **Some of the dependencies fail to install**: That is rarely the case, as all dependencies are in the main Python repositories. However, there might be issues with Numpy, Matplotlib and Scipy. These packages need to be compiled and it is often easy to get them as a scientific bundle. You can check `Anaconda <https://www.continuum.io/downloads>`_ which provides all these packages together already configured for the correct OS.

The config_tools
^^^^^^^^^^^^^^^^

This module contains all the functions that will help you to setup and configure your Solcore installation, as it has been highlighted above. The full description of the funcitons included in this module are:

.. automodule:: solcore.config_tools
    :members:
    :undoc-members:
