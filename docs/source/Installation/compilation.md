# Fortran Compilers

In order to use the Poisson-Drift-diffusion solver, it will be necessary to have a suitable fortran compiler installed and correctly configured in your system. After that, you do not need to worry about the compilation process: it is done automatically by the installation script. 

## Linux

Most linux systems have a Fortran compiler already installed, typically part of GCC, the GNU compiler collection. If not, you will need to check how to install it for your particular linux distribution. 

## Mac OS X

For Mac OS X we have used the *gfortran* compiler installed together with GCC using [Homebrew](https://brew.sh), a MacOS package manager.  

```brew install gcc```

Other package managers like [MacPorts](https://www.macports.org) or [Fink](http://www.finkproject.org) might also work, but we have not tried. 

## Windows

To get a Fortran compiler properly working under Windows with F2Py, we have followed the [detailed instructions written by Michael Hirsch](https://www.scivision.co/f2py-running-fortran-code-in-python-on-windows/). This instructions recommend to use [Anaconda Python](https://www.anaconda.com/download/) (or [Miniconda](https://conda.io/miniconda.html)) and so do we: this is the only way we have found to make the PDD solver to work under Windows, appart from using the Ubuntu Shell that comes with Windows 10 and that [we describe here](Solcore_on_Windows.md). Read Michael Hirsch's instructions in full for a complete picture, but the important bits are:

1. Install [MinGW-W64](https://sourceforge.net/projects/mingw-w64/) to ```c:\mingw``` with:
    
    - Architecture: x86_64
    - Threads: win32        <-- posix, as indicated in the instructions, did not work for us
    - Exception: seh
    
2. Insert the *mingw* bin path to your system path permanently: Control Panel → Advanced System Settings → Environment Variables → System variables → Path (you might need to scroll down the System variables window to find *Path*). To the end of the list of values for the variable *Path* add the path to *gfortran.exe*, such as:

    ```c:\mingw\mingw64\bin```
    
3. Tell Python to use MinGW: create file ```c:\path_to_Anaconda\Lib\distutils\distutils.cfg``` containing:

    ```
    [build]
    compiler=mingw32
    ```

And that should be it. With this, your will be able to call the gfortran compiler from the Windows terminal and also the Solcore installer should be able to find it. You can try installing Solcore itself or with the *lowtran* package described in the instructions by [Michael Hirsch](https://www.scivision.co/f2py-running-fortran-code-in-python-on-windows/). 