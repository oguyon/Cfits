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

The PIAACMC design process is as follows:
-# design an idealized monochromatic PIAACMC for a centrally obscured aperture (steps 1-4)
-# modify the design for the pupil aperture (steps 5-)


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
- PIAACMC_nblambda : number of wavelengths (1 for monochromatic)
- PIAACMC_dftgrid : dft grid sampling gap [pix] in pupil plane

\verbatim
<executable>
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
PPIAACMC_centobs0=0.3
PIAACMC_centobs1=0.15
PIAACMC_fpmradld=1.0
PIAACMC_dftgrid=2
piaacmcsimrun 0 0
exit
\endverbatim

This will create the prolate function for a centrally obscured pupil, the idealized focal plane mask, and run the diffraction propagation.\n
This step takes a few minutes - most of the time is spent on iterations to compute the 2D apodization prolate function.\n

Run this twice: once to set up the configuration, and once to run the on-axis PSF.\n

After this step, the contrast will likely be around 1e-7. The next step to improve this nominal design is to find the optimal locations for the Lyot stops.\n

Example result (for 1 l/D mask, 2048 array size):
\verbatim
Peak constrast (rough estimate)= 2.02762e-06
Total light in scoring field = 9.49968e-06  -> Average contrast = 6.88685e-08
\endverbatim


### STEP 2: Find optimal location of the two Lyot stops


For circular centrally obscured pupils, the default Lyot stops configuration consists of two stops: one that masks the outer part of the beam, and one that masks the central obstruction.\n

Fine optimization of the stops locations is done with a separate command. The optional lsoptrange value is the range (unit: m) for the mask position search (default = 2 m).\n

\verbatim
<executable>
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
PIAACMC_lsoptrange=2.0
piaacmcsimrun 0 1
exit
\endverbatim

Example result (for 1 l/D mask, 2048 array size):
\verbatim
Peak constrast (rough estimate)= 1.78874e-06
Total light in scoring field = 7.75534e-06  -> Average contrast = 5.62228e-08
\endverbatim



Optional: you can try to tune the geometry of the masks with the following routine:
\verbatim
<executable>
PIAACMC_nblstop=2
PIAACMC_lstransm=0.9
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
piaacmcsimrun 0 5
exit
\endverbatim



### STEP 3 (optional): Optimize focal plane mask transmission


\verbatim
<executable>
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
piaacmcsimrun 0 2
exit
\endverbatim

See file result.log to track progress:
-# col 1: mask transm
-# col 2: average contrast
-# col 3: scan iteration number
-# col 4: test range for scan
-# col 5: step size for scan


### STEP 4 (optional): Optimize PIAA shapes for high contrast


\verbatim
<executable>
PIAACMC_nbiter=5
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
piaacmcsimrun 0 4
exit
\endverbatim


This command will run for the specified number of iterations. Converging is monitored by reading ASCII file val.opt, and the script can be stopped (ctrl-C) anytime. You should see the contrast (col 2) improve. 



### Repeat steps 2-4 as needed to optimize the idealized PIAACMC system







### STEP 5: Specify input pupil geometry

Input to design is a pupil map, FITS file, size x size pixel, valued 1.0 inside the pupil and 0.0 outside. The file name should be pupa_<size>.fits (for example pupa_2048.fits).\n
Replace the file ./piaacmcconf000/pupa0_<size>.fits by the desired pupil.\n

Then, run the propagation for an on-axis point source:
\verbatim
<executable>
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
piaacmcsimrun 0 0
exit
\endverbatim

The contrat is now much poorer. For example:
\verbatim
Peak constrast (rough estimate)= 0.00209981
Total light in scoring field = 0.0107966  -> Average contrast = 7.82706e-05
\endverbatim

### STEP 6: Compute Lyot Stops

For circular centrally obscured pupils, the default Lyot stops configuration consists of two stops: one that masks the outer part of the beam, and one that masks the central obstruction.\n

The following command computes the Lyot stops shapes and approximate locations. It is required for non-circular pupils, and can be re-run at any time.\n

\verbatim
<executable>
PIAACMC_nblstop=4
PIAACMC_lstransm=0.8
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
piaacmcsimrun 0 5
exit
\endverbatim

Fine optimization of the stops locations is done with a separate command. The optional lsoptrange value is the range (unit: m) for the mask position search (default = 2 m).\n

\verbatim
<executable>
PIAACMC_dftgrid=2
PIAACMC_nblambda=1
PIAACMC_lsoptrange=1.0
piaacmcsimrun 0 1
exit
\endverbatim





### STEP 4: Optimize PIAA shapes for high contrast


\verbatim
<executable>
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
piaacmcsimrun 0 4
exit
\endverbatim


This command will run for a large number of iterations. Converging is monitored by reading ASCII file val.opt, and the script can be stopped (ctrl-C) anytime.



