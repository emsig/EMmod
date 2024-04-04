# EMmod

This repository contains the code that comes with the publication

> Hunziker, J., J. Thorbecke, and E. Slob, 2015,  
> The electromagnetic response in a layered vertical transverse isotropic medium:  
> A new look at an old problem:  
> Geophysics, 80(1), F1-F18; 
> DOI: [10.1190/geo2013-0411.1](https://doi.org/10.1190/geo2013-0411.1);  

The original code was published with the article in a static archive by the SEG, see [wiki.seg.org/wiki/Software:emmod](https://wiki.seg.org/wiki/Software:emmod). The version in the SEG archive is (c) 2015 by the Society of Exploration Geophysicists, for more info consult the file [wiki.seg.org/wiki/Software:Disclaimer](https://wiki.seg.org/wiki/Software:Disclaimer).

We[^1] release our source code here under the CPL-1.0 license, see the file `Common_Public_License.txt`. Please cite the above article if you use the code in your research.


## License Info

THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THIS COMMON PUBLIC LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT.

A copy of this license can be found in the file `Common_Public_License.txt` in the directory where you have found this README.

http://www.opensource.org/licenses/cpl1.0.php

Some routines are from Seismic Unix. For those the `SU LEGAL_STATEMENT`, which is included in the  source code, applies and not the common public license. 

Some routines are from slatec: http://www.netlib.org/slatec/. For those, the slatec guidelines apply (http://www.netlib.org/slatec/guide) and not the common public license. Routines from slatec include one of the following lines in the corresponding source file: 

```
    LIBRARY   SLATEC
    LIBRARY   SLATEC (PCHIP)
    LIBRARY   SLATEC (QUADPACK)
    LIBRARY   SLATEC (XERROR)
```

Note: Some of these routines are slightly modified. Therefore, using the original library from slatec will not work with EMmod.

## Compile the code using windows

1. Go to www.cygwin.org and download the newest version of cygwin

2. Run the setup-file, which will install cygwin. During the installation, you can select packages to be installed. Next to the basic install, also select the following packages from the devel (development) section:

   - make: The GNU version of the 'make' utility
   - gcc-core: GNU compiler collection (C, OpenMP)
   - gcc-fortran: GNU compiler collection (Fortran)

   If you have cygwin installed already, but not these packages, just run the setup-file again and select the missing packages. 

3. In the directoy, in which you have installed cygwin, you will find a home-directory. This home directory contains another directory with your user name. Copy all files related to the EMmod code into a folder in the directory with your user name. 

4. Open the cygwin terminal and follow the steps described below on how to "compile the code using a linux or a linux-like system".

## Compile the code using a linux or a linux-like system

1. To compile and link the code you first have to set the ROOT variable in the `Make_include` file which can be found in the directory where you have found this README.

2. Check the compiler, `CFLAGS` and `FFLAGS` options in the file `Make_include` and adapt to the system you are using. The default options are set for the GNU C and Fortran compiler on a Linux system. The Makefile has been tested with GNU make.

3. If the compiler options are set in the `Make_include` file you can type 

```
make 
```

If the compilation has finished without errors and produced an executable called emmod you can run one of the following demo programs found in the subdirectory emmod by running 

```
./halfspacemod.scr
./simplemod.scr
./gprloop_twointerface.scr
```

These three scripts and variations thereof were used to compute the three examples presented in the accompanying GEOPHYSICS paper "The electromagnetic response in a layered VTI medium: A new look at an old problem". 

Information about the parameters are found in the file `manual.pdf` in the subdirectory doc. 

If you get the error `emmod: command not found`, write a dot and a forward slash `./` in front of emmod in the scr-file.

In the subfolder matlab, Matlab scripts can be found to load and plot EMmod data as well as to perform the Fourier transformation over frequency to create time-domain GPR datasets. 

If you use this code for a scientific publication, please cite the accompanying paper mentioned above. 

Depending on the compiler an error about missing object files for the Bessel functions J0, J1 can be encountered during linking on the emmod executable. These missing symbols are defined in `libm.a`; the mathematical C library. Some C-compilers add this library automatically during linking, others don't.  When you encounter this error you can solve it by adding `-lm` to line 25 in the `Make_include` file: 

```
  LIBS = -lgfortran -lm
```

[^1]: 2024-04-04: The code was uploaded to
  [github.com/emsig/EMmod](https://github.com/emsig/EMmod) by
  [@prisae](https://github.com/prisae),
  upon the request of the main author Jürg Hunziker.
