Installing S4
--------------

Installing S4 is only necessary if you would like to use Solcore's rigorous coupled-wave analysis (RCWA) functionality. 
This is a wave optical method for solving Maxwell's equations in a periodic structure, and can be used for structures which
contain e.g. a diffraction grating or other photonic structures. S4 is written in C++ and must be compiled on your system
before you are able to import it as a Python package. Up-to-date instructions on how to do this in MacOS
and Ubuntu are maintained [here](https://github.com/phoebe-p/S4) (scroll down to the README section). We do not currently have 
installation instructions for other Linux distributions, but it should  be possible to compile S4 on any UNIX operating system 
as long as the relevant packages are installed first (see the link above for requirements). Unfortunately,
there are currently no instructions available for installing S4 on Windows; if you want to use the RCWA functionality from
a Windows machine there are [a few options](Solcore_on_Windows.md).
