Solcore on Mac OS X
=====================

Solcore has been heavily tested under Mac OS X, so there should be no problems. However, using NGSPICE requires some details, as shown below.

Installing NGSpice on MacOSX
----------------------------

* Download NGSpice from 
[http://ngspice.sourceforge.net/download.html](http://ngspice.sourceforge.net/download.html)
* Installing the binary package on mac os X will result in the executable residing in ```/Applications/ngspice/bin/ngspice```
* You an set the path using 
``` `\>\>\>import solcore.config\_tools as config
> > > config.set_location_of_spice('/Applications/ngspice/bin/ngspice')


* If NG spice throws an error
``` `dyld: Library not loaded: /opt/X11/lib/libXaw.7.dylib```
`* then you will likely need to install xquartz.  It can be downloaded from from [https://www.xquartz.org](https://www.xquartz.org)
	* run the installer and restart the computer.
* You can test the NGSpice installation using Example\_PV\_module.py



