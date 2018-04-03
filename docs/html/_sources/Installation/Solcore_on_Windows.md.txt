Solcore on Windows 10
=====================
The Python part of Solcore, which is the majority, should work under Windows with no problems. However, we have found a lot trouble trying to make the parts and complements that require compilation to work properly (the PDD solver, S4, SMARTS and SPICE). For those who need to use those tools, please follow the instructions below. 

Installing Solcore
------------------
(thanks Phoebe!)

After a lot of effort, we have Solcore fully up and running in Windows 10... more or less, using the Ubuntu shell that comes with Windows 10. To install Solcore there, follow this steps:

All steps on fresh install of Ubuntu (using the Ubuntu terminal on
Windows 10 distributed by Canonical Group, Ubuntu 16.04.3 LTS, Codename:
xenial)

-  Install git if not done already
-  Install python 3.x if not already done
-  Install pip3 (for installing Python3 packages; you may need to update the package list first: sudo apt-get update)::

```bash
sudo apt install python3-pip
```
- You need LAPACK and BLAS libraries linked to the -llapack and -lblas library flags – these are used for scipy and S4. The developers of S4 recommend OpenBLAS (you can find installation instructions by Googling), but this also works and is simpler:

```bash
sudo apt-get install libblas-dev liblapack-dev
```

-  Install NumPy:
    
```bash
pip3 install numpy
```
-  Install matplotlib (and tk, otherwise get an error later):
    
```bash
pip3 install matplotlib
sudo apt-get install python3-tk 
```

-  Other dependencies install automatically & successfully when
      installing Solcore5... hopefully.

-  Now, we actually Install Solcore::

```bash
git clone https://github.com/dalonsoa/solcore5.git
cd solcore5
sudo python3 setup.py install
```

Installing S4
-------------

-  The “make” command must be available:

```bash
sudo apt install make
```

-  You must use the fork of S4 at https://github.com/phoebe-p/S4; the
      main branch is not compatible with Python 3.x:

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

This might result in an error saying that quantum mechanics failed because 5\ :sup:`th` decimal place of result doesn’t match. Simply, ignore it.

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

