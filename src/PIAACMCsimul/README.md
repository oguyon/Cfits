# PIAACMC coronagraph design 


## Overview

Diffraction-based PIAACMC simulation / optimization
- Uses Fresnel propagation engine ( OptSystProp.c OptSystProp.h ) between optical elements
- Computes linear perturbations around current design for optimization
- Automatically computes multiple Lyot stop to follow variable conjugation in different parts of the beam
- Polychromatic propagations
- Fits aspheric PIAA shapes on basis of radial cosines and 2-D Fourier modes


## Design Steps for PIAACMC 

### Overview

Directories "./piaacmcconf<nnn>/" hold the default/current configuration settings and files for the PIAACMC. The directory will be automatically created if it does not exist. By copying from/to this directory, you can save/load PIAACMC designs.\n

By default, if no configuration file exists, a monochromatic PIAACMC for a centrally obscured pupil will be created. This is meant as a starting point for the PIAACMC, which will then be optimized further\n

The main command to run the PIAACMC simulation is
\verbatim
<executable>
piaacmcsimrun <nnn> <mode>
exit
\endverbatim

where <nnn> [long] is the configuration index and <mode> [long] defines the operation to be performed.\n
See function PIAACMCsimul_run(long confindex, long mode) for list of modes\n

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
Optional parameters are:
- PIAACMC_centobs0 : central obstruction in input pupil (default = 0.3)
- PIAACMC_centobs1 : central obstruction after remapping (default = 0.2)
- PIAACMC_fpmradld : focal plane mask radius for monochromatic idealized design, in l/D system unit (default = 0.9)

\verbatim
<executable>
PIAACMC_centobs0=0.3
PIAACMC_centobs1=0.15
PIAACMC_fpmradld=1.0
piaacmcsimrun 0 0
exit
\endverbatim

This will create the prolate function for a centrally obscured pupil, the idealized focal plane mask, and run the diffraction propagation.\n

This step takes a few minutes - most of the time is spent on iterations to compute the 2D apodization prolate function.\n




### STEP 2: Specify input pupil geometry

Input to design is a pupil map, FITS file, size x size pixel, valued 1.0 inside the pupil and 0.0 outside. The file name should be pupa_<size>.fits (for example pupa_1024.fits).\n
Replace the file ./piaacmcconf000/pupa0_<size>.fits by the desired pupil.\n

\verbatim
<executable>
PIAACMC_dftgrid=2
piaacmcsimrun 0 0
exit
\endverbatim



### STEP 3: Compute Lyot Stops

For circular centrally obscured pupils, the default Lyot stops configuration consists of two stops: one that masks the outer part of the beam, and one that masks the central obstruction.\n

The following command computes the Lyot stops shapes and approximate locations. It is required for non-circular pupils, and can be re-run at any time.\n

\verbatim
<executable>
PIAACMC_nblstop=3
PIAACMC_lstransm=0.85
piaacmcsimrun 0 5
exit
\endverbatim

Fine optimization of the stops locations is done with a separate command. The optional lsoptrange value is the range (unit: m) for the mask position search (default = 2 m).\n

\verbatim
<executable>
PIAACMC_dftgrid=2
PIAACMC_lsoptrange=0.2
piaacmcsimrun 0 1
exit
\endverbatim


### STEP 4: Optimize PIAA shapes for high contrast


\verbatim
<executable>
PIAACMC_dftgrid=2
piaacmcsimrun 0 4
exit
\endverbatim


This command will run for a large number of iterations. Converging is monitored by reading ASCII file val.opt, and the script can be stopped (ctrl-C) anytime.



