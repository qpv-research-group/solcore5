Solcore on Windows 10
=====================
The Python part of Solcore, which is the majority, should work under Windows with no problems. However, we have had some
trouble trying to make the parts and complements that require compilation to work properly (the PDD solver, S4, SMARTS 
and SPICE). For those who need to use those tools and encounter issues on Windows, please follow the instructions below. 

Note that in principle, out of the packages listed above (the PDD solver, S4, SMARTS and SPICE), all can be
installed on Windows **except S4**, which to date we have not been able to produce Windows installation instructions for. 
To compile the PDD, you need to do [a few extra steps](compilation.md) to get a suitable compiler working.
S4 is only required to use the rigorous coupled-wave analysis (RCWA) functionality in Solcore, used to calculate diffraction
from periodic structures (e.g. diffraction gratings, photonic crystals). Thus, if you do not need this functionality, you 
can stick with Windows.

Installing Solcore
------------------

It is possible to run all the parts of Solcore from a Windows 10 environment without using a virtual machine or dual-booting a UNIX operating
system (e.g. Ubuntu). This can be done by using the Ubuntu shell that comes with Windows 10. **It should be noted that this Windows
Subsystem on Linux (WSL) does not support graphical applications, so it is purely a command line environment.** This can be 
inconvenient for actually writing code and viewing results, so please bear this in mind before continuing with the instructions below.
If you are having problems on Windows (or want to use the RCWA functionality/S4) and do not want to be limited to a command-line only
environment, we recommend installing Ubuntu either as a virtual machine, or as a dual boot alongside Windows 10. If you do want to 
install Solcore on WSL, follow these steps:

All steps on fresh install of Ubuntu (using the Ubuntu terminal on
Windows 10 distributed by Canonical Group, Ubuntu 16.04.3 LTS, Codename:
xenial)

- Install git if not done already
- Install python 3.x if not already done
- Install pip3 (for installing Python3 packages; you may need to update the package list first: `sudo apt-get update`):

    ```bash
    sudo apt install python3-pip
    ```

- Install matplotlib (and tk, otherwise get an error later):
    
    ```bash
    pip3 install matplotlib
    sudo apt-get install python3-tk 
    ```

- Install Solcore5, which can be done according to the [standard instructions](installation.rst).


Installing S4
-------------

- The “make” command must be available:

    ```bash
    sudo apt install make
    ```
   
- You need LAPACK, BOOST, FFT and BLAS libraries to compile S4. The developers of S4 recommend OpenBLAS (you can find installation 
   instructions by Googling), but this also works and is simpler:
    
    ```bash
    sudo apt install libopenblas-dev libfftw3-dev libsuitesparse-dev libboost-all-dev
    ```

  (if you have any issues, more detailed installation instructions can be found [here](https://github.com/phoebe-p/S4).
- You must use the fork of S4 at [https://github.com/phoebe-p/S4](https://github.com/phoebe-p/S4), as the
      main branch is not compatible with Python 3:

    ```bash
    git clone https://github.com/phoebe-p/S4.git
    cd S4
    make S4_pyext
    ```

Checking if everything works 
--------------------------

- To see if the PDD is working:

    ```bash
    python3
    >>> import solcore.solar_cell_solver
    ```

-  Run tests:

    ```bash
    sudo python3 setup.py test
    ```

This might result in an error saying that quantum mechanics failed because 5\ :sup:`th` decimal place of result doesn’t match. 
You can simply ignore this error.

- Other information:
    - gcc/g++/gfortran versions used here: gcc version 5.4.0 20160609 (Ubuntu 5.4.0-6ubuntu1~16.04.6)

    -  If everything is working correctly, should be able to import
      solcore.solar_cell_solver without getting warnings about the RCWA
      solver or PDD.

    -  Sometimes, an issue arises with the SOPRA database (and occasionally
   importing solar_cell_solver?) where the permissions are set such that
   you cannot access the database – running sudo python3, or sudo
   [whatever python editor you’re using] fixes this.

Installing SMARTS Windows10/Ubuntu Shell 
----------------------------------------
(thanks Andrew!)

* Register and download SMARTS FOR MAC from http://www.nrel.gov/rredc/smarts/    Yes, we are going to use a Mac software on Windows. Hold on tight.
* Move the SMARTS_xxx_Mac.tar.gz file from Windows to the dedicated Ubuntu folder (e.g. ~$ mv /mnt/c/Users/Andrew/Downloads ~/)
* Unpack the file (~$ tar zxf SMARTS_xxx_Mac.tar.gz)
* The executable script "smarts295.command" in the folder won't work (try running ~$ ./smarts295.command). For example on 64bit Windows the following error is thrown, requiring the source to be recompiled:

```
Source_code/smarts295_cmd: 1: Source_code/smarts295_cmd: dJ: not found
Source_code/smarts295_cmd: 1: Source_code/smarts295_cmd: ȅ8__PAGEZERO__TEXT: not found
Source_code/smarts295_cmd: 1: Source_code/smarts295_cmd: Syntax error: word unexpected (expecting ")")
```

#### Recompiling SMARTS   

* Move to the directory /SMARTS_295_Mac/Source_code/
* Delete the file "smarts295_cmd" (type ~$ rm smarts_cmd). This will need to be recompiled on for your CPU architecture
* Assuming you have gfortran installed (either as part of Anaconda or another third party distribution) you should compile smarts (type ~$ gfortran smarts295.f) which results in an executable called "a.out"
	* N.B. Although the SMARTS documentation suggests changing some flags in the source code you _should not modify_ the SMARTS source code.  Solcore expects the default behaviour with respect to batch mode and file overwrites.  
* Rename the file "a.out" to "smarts295_cmd". Now SMARTS has been recompiled for your CPU architecture.
* IMPORTANT FOR WINDOWS: Solcore on Windows10/Ubuntu looks for the script called "smarts295", instead of "smarts295.command", due to different naming systems in Linux and Mac. Rename "smarts295.command" to just "smarts295".
* You can test the installation of SMARTS using the example file Figure8.py
* Please note that running SMARTS from the command line (typing ~$ ./smarts295) will produce two text files, "smarts295.ext.txt" and "smarts295.out.txt".
   If these are in the SMARTS root folder, they will mess up running SMARTS on Solcore. If you do run SMARTS in the command line, remove these files when you are done.

