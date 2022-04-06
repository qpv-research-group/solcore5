# Fortran Compilers

In order to use the Poisson-Drift-diffusion solver, it will be necessary to have a suitable fortran compiler installed and correctly configured in your system. After that, you do not need to worry about the compilation process: it is done automatically by the installation script. 

## Linux

Most Linux systems have a Fortran compiler already installed, typically part of GCC, the GNU compiler collection. If not, you will need to check how to install it for your particular linux distribution. 

## Mac OS X

For Mac OS X we have used the *gfortran* compiler installed together with GCC using [Homebrew](https://brew.sh), a MacOS package manager.  

```brew install gcc```

Other package managers like [MacPorts](https://www.macports.org) or [Fink](http://www.finkproject.org) might also work, but we have not tried. 

There is currently a known issue involving the compilation of the PDD solver on new Macs with M1/ARM/Apple silicon chips. 
The homebrew version of gfortran (installed above alongside gcc as indicated above) does not work. For a workaround please
see [here](https://github.com/qpv-research-group/solcore5/issues/209).

## Windows

To get a Fortran compiler properly working under Windows with F2Py, we have followed the [detailed instructions written by Michael Hirsch](https://www.scivision.co/f2py-running-fortran-code-in-python-on-windows/). Read Michael Hirsch's instructions in full for a complete picture, but the important bits are:

1. Install [MSYS2 and MinGW](https://www.scivision.co/install-msys2-windows)

	Follow the instructions there to update all packages and to setup your environment.

2. Install the fortran compiler executing from the PowerShell:

	```bash
	pacman -S mingw-w64-x86_64-gcc-fortran
	````
    
3. Tell Python to use MinGW: create file ```~/pydistutils.cfg``` containing:

    ```
    [build]
    compiler=mingw32
    ```
	
	You can do this executing in the PowerShell
	
	```bash
	echo "[build]`ncompiler=mingw32" | Out-File -Encoding ASCII ~/pydistutils.cfg
	````

And that should be it. With this, your will be able to call the gfortran compiler from the Windows terminal and also the Solcore installer should be able to find it. You can try installing Solcore itself or with the *lowtran* package described in the instructions by [Michael Hirsch](https://www.scivision.co/f2py-running-fortran-code-in-python-on-windows/). 