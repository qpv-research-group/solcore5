Solcore on Windows 10
=====================

These instructions should be relevant only if you are a Windows user who would like to use Solcore's rigorous coupled-wave
analysis (RCWA) functionality, requiring the S4 package, for which we unfortunately do not have installation instructions 
for Windows. All other parts of Solcore, and the optional SMARTS and SPICE dependencies, can be installed on Windows. S4 
is only required to use the rigorous coupled-wave analysis (RCWA) functionality in Solcore, used to calculate diffraction
from periodic structures (e.g. diffraction gratings, photonic crystals). Thus, if you do not need this functionality, you 
can stick with Windows.

Installing Solcore
------------------

It is possible to run all the parts of Solcore from a Windows 10 environment without using a virtual machine or dual-booting a UNIX operating
system (e.g. Ubuntu). This can be done by using the Windows Subsystem for Linux (WSL) which can be installed on Windows 10. 
By itself,
the WSL does not support graphical applications, but it is possible to use Python installed in the WSL from a development
environment running in Windows (such as VSCode), including plotting capabilities. Step-by-step instructions for
setting up the WSL and accessing it from VSCode are given [here](https://code.visualstudio.com/docs/remote/wsl-tutorial).
Note that to be able to see plots generated during code execution you can right-click on the Python file you want to execute
in VSCode and select 'Run Current File in Python Interactive window'. This has been tested on 64-bit Windows 10 using
the Ubuntu20.04 version of the
WSL (using [WSL 2](https://docs.microsoft.com/en-us/windows/wsl/install#upgrade-version-from-wsl-1-to-wsl-2)); both S4
and the PDD install and run without issues.

If you would like to use a different development environment, you will have to check if it supports accessing Python from
the WSL from a Windows application; it
is possible in [PyCharm](https://www.jetbrains.com/help/pycharm/using-wsl-as-a-remote-interpreter.html).

However, depending on your preferences you may still want to use a virtual machine or dual boot your computer to run a 
complete Linux operating system. Whichever option you choose, you should be able to install Solcore in your preferred way
using the normal installation
instructions; instructions for installing S4 on Ubuntu, if you want to use it, are given below.

Installing S4
-------------

These are the steps for installing S4 in an Ubuntu WSL environment (through the WSL command line environment). These steps
should be the same for any Ubuntu environment.

- The “make” and “git“ commands must be available:

    ```bash
    sudo apt install make git
    ```
   
- You need LAPACK, BOOST, FFT and BLAS libraries to compile S4. The developers of S4 recommend OpenBLAS (you can find installation 
   instructions by Googling), but this also works and is simpler:
    
    ```bash
    sudo apt install libopenblas-dev libfftw3-dev libsuitesparse-dev libboost-all-dev
    ```

  (if you have any issues, more detailed installation instructions can be found [here](https://github.com/phoebe-p/S4).)
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

This might result in an error saying that quantum mechanics failed because the 5th decimal place of result doesn’t match. 
You can simply ignore this error. If everything is working correctly, should be able to `import solcore.solar_cell_solver` 
without getting warnings about the RCWA solver or PDD.
  

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

