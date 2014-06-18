# PIAACMC coronagraph design 


## Overview

Diffraction-based PIAACMC simulation / optimization
- Uses Fresnel propagations engine ( OptSystProp.c OptSystProp.h ) between optical elements
- Computes linear perturbations around current design for optimization
- Automatically computes multiple Lyot stop to follow variable conjugation in different parts of the beam
- Polychromatic propagations
- Fits aspheric PIAA shapes on basis of radial cosines and 2-D Fourier modes


## Design Steps

### Overview

Directories "./piaacmcconf<nnn>/" hold the default/current configuration settings and files for the PIAACMC. The directory will be automatically created if it does not exist. By copying from/to this directory, you can save/load PIAACMC designs.\n

By default, if no configuration file exists, a monochromatic PIAACMC for a centrally obscured pupil will be created. This is meant as a starting point for the PIAACMC, which will then be optimized further\n

The main command to run the PIAACMC simulation is
\verbatim
piaacmcsimrun <nnn> <mode>
\endverbatim

where <nnn> [long] is the configuration index and <mode> [long] defines the operation to be performed.\n
See function dsPIAACMCsimul_run(long confindex, long mode) for list of modes\n

### Initialization rules (function  PIAAsimul_initpiaacmc() )

-# if configuration directory exists, use it and load configuration file ( function  PIAAsimul_loadpiaacmcconf ), otherwise, create it
-# load/create Cmodes 
-# load/create Fmodes
-# load mode coefficients for piaa shapes if they exist. If not:
	-# create radial apodization for centrally obscured idealized monochromatic PIAACMC
	-# fit / extrapolate radial apodization profile with cosines
	-# using above fit, create 2D radial sag for both PIAA optics ( -> PIAA_Mshapes.txt)
	-# make 2D sag maps for both optics ( -> piaa0z.fits, piaa1z.fits)
	-# fit 2D sag maps with Cmodes and Fmodes coefficients ( -> piaa0Cmodes, piaa0Fmodes, piaa1Cmodes, piaa1Fmodes )
-# load/create focal plane mask zone map. This is the map that defines the geometry (which ring is where)
-# load/create focal plane mask thickness array
-# load/create focal plane mask transmission array
-# load/create Lyot stops


### STEP 1: Create an idealized centrally obscured apodized PIAACMC monochromatic design

This is meant as a starting point for the PIAACMC, which will then be optimized further\n

\verbatim
piaacmcsimrun 0 10
\endverbatim




-# specify input pupil geometry\n
Input to design is a pupil map, FITS file, size x size pixel, valued 1.0 inside the pupil and 0.0 outside. The file name should be pupa_<size>.fits (for example pupa_1024.fits).



### Create piaacmcparams.conf file



