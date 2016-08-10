% PIAACMC design
% Olivier Guyon
% 2015-present







# Overview


## Scope

Diffraction-based PIAACMC simulation / optimization

- Uses Fresnel propagation engine ( OptSystProp.c OptSystProp.h ) between optical elements
- Computes linear perturbations around current design for optimization
- Automatically computes multiple Lyot stop to follow variable conjugation in different parts of the beam
- Polychromatic propagations
- Fits aspheric PIAA shapes on basis of radial cosines and 2-D Fourier modes


## Usage

Scripts to run the software are located within the source code directory:

	./src/PIAACMCsimul/scripts/

The scripts can be linked to your working directory by executing the following command:

	ln -s $PWD/syncscripts /myworkdirectory/syncscripts

Then, execute in your work directory:

	./syncscripts
	
This will install all required scripts in workdirectory and install any packages required.


Code is composed of a several layers (from high to low) :

-------------------- -----------------------------------------------------------
Script                   Description
-------------------- -----------------------------------------------------------
**run**                  Top level script, calls **runopt** script

**runopt**               Optimize a PIAACMC design or run an existing design
						calls **sim** script

**sim**                   calls **runPIAACMC** script

**runPIAACMC**            lower-level scropt calling C-written executable
-------------------- -----------------------------------------------------------




# Scripts

## Script run

~~~~~~~~~~~ {.bash .numberLines}
!INCLUDE "../scripts/run"
~~~~~~~~~~~

## Script runopt

~~~~~~~~~~~ {.bash .numberLines}
!INCLUDE "../scripts/runopt"
~~~~~~~~~~~

## Script sim

~~~~~~~~~~~ {.bash .numberLines}
!INCLUDE "../scripts/sim"
~~~~~~~~~~~

## Script runPIAACMC

~~~~~~~~~~~ {.bash .numberLines}
!INCLUDE "../scripts/runPIAACMC"
~~~~~~~~~~~



# C code description

## PIAACMC_designcodes


The main function in the source code is PIAACMCsimul_exec(), which takes two arguments: the configuration index (usually a 3 digit integer) and the mode (integer) which describes the operation to be performed to the PIAACMC design.

Mode		| Description
----------------|-------------------------------------
0	*	| Compute on-axis propagation for specified configuration. If configuration does not exist, create idealized monochromatic PIAACMC (intended to create a new index) and compute on-axis propagation
1	*	| Optimize Lyot stop(s) locations (scan) 
2	*	| Optimize focal plane mask transmission for idealized monochromatic PIAACMC (scan)
3		| Run without focal plane mask (for testing and calibration)
4		| Linear optimization around current design, free parameters = PIAA optics cosines shapes (track progress by looking at val.opt file)
5	*	| Optimize Lyot stops shapes and positions
10	*	| Setup polychromatic optimization
11	*	| Compute polychromatic response to zones, store result in FPMresp
12		| Search for best mask solution using FPMresp, random search
13	*	| Optimize PIAA optics shapes and focal plane mask transmission (idealized PIAACMC)
40	*	| Optimize PIAA optics shapes and focal plane mask transmission (idealized PIAACMC)
41		| Optimize PIAA optics shapes and focal plane mask zones (polychromatic)
100 	*	| Evaluate current design: polychromatic contrast, pointing sensitivity
101 	*	| Transmission as a function of angular separation
200 	*	| Make focal plane mask OPD map



The following variables can be set 

Variable		| Description
------------------------|-------------------------------------
PIAACMC_size		| Array size (1024, 2048, etc...)
PIAACMC_pixscale	| Pixel scale [m]
PIAACMC_dftgrid		| Sampling interval in DFTs
PIAACMC_centobs0	| Input central obstruction
PIAACMC_centobs1	| Output central obstruction
PIAACMC_nblambda	| Number of wavelength points
PIAACMC_resolved	| 1 if resolved source (3 points at r = 0.01 l/D, 120 deg apart)
PIAACMC_fpmtype		| 1 if physical mask, 0 if idealized mask
PIAACMC_FPMsectors	| Number of sectors in focal plane mask
PIAACMC_NBrings		| Number of rings in focal plane mask
PIAACMC_fpmradld	| Focal plane mask outer radius




# Design Steps for PIAACMC 



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








# Initialization rules (function  PIAAsimul_initpiaacmc() )


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








## MODE 000: Create an idealized centrally obscured apodized PIAACMC monochromatic design

This is meant as a starting point for the PIAACMC, which will then be optimized further\n
Optional parameters are:
- PIAACMC_size : array size (default = 1024)
- PIAACMC_pixscale : pixel scale (m/pix)
- PIAACMC_nblambda : number of wavelengths (1 for monochromatic)
- PIAACMC_dftgrid : DFT grid sampling gap [pix] in pupil plane
- PIAACMC_centobs0 : central obstruction in input pupil (default = 0.3)
- PIAACMC_centobs1 : central obstruction after remapping (default = 0.2)
- PIAACMC_fpmradld : focal plane mask radius for monochromatic idealized design, in l/D system unit (default = 0.9)


\verbatim
<executable>
PIAACMC_size=1024
PIAACMC_pixscale=0.00011
PIAACMC_nblambda=1
PIAACMC_dftgrid=2
PIAACMC_centobs0=0.3
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


### Code breakdown

- PIAACMCsimul_exec() :
	- PIAAsimul_initpiaacmcconf(): Load/Creates/initializes piaacmcconf structure and directory
		- perform default initialization
		- PIAAsimul_loadpiaacmcconf(): Loading PIAACMC configuration from "piaacmcconfxxx/piaacmcparams.conf" if it exists
		- Creating/loading Cmodes and Fmodes
		- IMPORT / CREATE PIAA SHAPES
			- create 2D prolate iteratively
			- PIAACMCsimul_load2DRadialApodization():fit PIAA shapes with Cosine modes
			- PIAACMCsimul_init_geomPIAA_rad(): compute radial PIAA sag from cosine apodization fit
			- PIAACMCsimul_mkPIAAMshapes_from_RadSag(): Make 2D sag shapes from radial PIAA sag
		- MAKE FOCAL PLANE MASK
		- MAKE LYOT STOPS
		- PIAAsimul_savepiaacmcconf(): save configuration
	- PIAACMCsimul_makePIAAshapes(): construct PIAA shapes from fitting coefficients
		- construct 2D PIAA mirror shapes from piaa0Cmodescoeff, piaa0Fmodescoeff, piaa1Cmodescoeff, piaa1Fmodescoeff
	- PIAACMCsimul_computePSF(): Compute PSF

### Output Files


APLCmaskCtransm.txt ??
fpm_ampl.fits
fpm_pha.fits
FPmask.tmp.fits 


#### Configuration : 

Output file	| Notes
----------------|-------------------------------------
./piaaconfxxx/piaacmcparams.conf | Configuration parameters
./piaaconfxxx/conjugations.txt	| Conjugations
./piaaconfxxx/lambdalist.txt  | list of wavelength values
./piaaconfxxx/pupa0_[size].fits	| input pupil (created by default if does not exist)

#### Wavefront Modes :

Output file	| Notes
----------------|-------------------------------------
Cmodes.fits	| circular radial cosine modes (40 modes, hard coded)
Fmodes.fits	| Fourier modes (625 modes = 10 CPA, hard coded)
./piaaconfxxx/ModesExpr_CPA.txt | modes definition
./piaaconfxxx/APOmodesCos.fits	| Cosine modes for fitting 2D apodization profile


#### PIAA mirrors, apodization, fits:

Output file	| Notes
----------------|-------------------------------------
./piaaconfxxx/APLCapo.1.400.0.300.info | file written by prolate generation function coronagraph_make_2Dprolate in coronagraphs.c
./piaaconfxxx/apo2Drad.fits	| idealized PIAACMC 2D apodization
./piaaconfxxx/piaam0z.fits	| PIAA M0 shape (2D sag)
./piaaconfxxx/piaam1z.fits	| PIAA M1 shape (2D sag)
./piaaconfxxx/PIAA_Mshapes.txt	| PIAA shapes (radial txt file, cols: r0, z0, r1, z1)
./piaaconfxxx/piaa0Fz.fits	| PIAA M0 shape, Fourier components (2D file)	
./piaaconfxxx/piaa1Fz.fits	| PIAA M1 shape, Fourier components (2D file)
./piaaconfxxx/piaa0Cmodes.fits  | idealized PIAACMC mirror 0 cosine modes (copied from ./piaaref/)
./piaaconfxxx/piaa0Fmodes.fits  | idealized PIAACMC mirror 0 Fourier modes (copied from ./piaaref/)
./piaaconfxxx/piaa1Cmodes.fits  | idealized PIAACMC mirror 1 cosine modes (copied from ./piaaref/)
./piaaconfxxx/piaa1Fmodes.fits  | idealized PIAACMC mirror 1 Fourier modes (copied from ./piaaref/)
./piaaconfxxx/piaa0Cres.fits	| idealized PIAA M0 cosine fit residual
./piaaconfxxx/piaa1Cres.fits	| idealized PIAA M1 cosine fit residual
./piaaconfxxx/piaa0Cz.fits	| idealized PIAA M0 cosine fit sag
./piaaconfxxx/piaa1Cz.fits	| idealized PIAA M1 cosine fit sag
./piaaconfxxx/piaa0Fz.fits	| idealized PIAA M0 Fourier fit sag
./piaaconfxxx/piaa1Fz.fits	| idealized PIAA M1 Fourier fit sag


#### Idealized PIAACMC reference point:

Output file	| Notes
----------------|-------------------------------------
./piaaconfxxx/piaaref/APLCmaskCtransm.txt | idealized PIAACMC focal plane mask transmission
./piaaconfxxx/piaaref/apo2Drad.fits	| idealized PIAACMC output apodization
./piaaconfxxx/piaaref/piaa0Cmodes.fits  | idealized PIAACMC mirror 0 cosine modes
./piaaconfxxx/piaaref/piaa0Fmodes.fits  | idealized PIAACMC mirror 0 Fourier modes
./piaaconfxxx/piaaref/piaa1Cmodes.fits  | idealized PIAACMC mirror 1 cosine modes
./piaaconfxxx/piaaref/piaa1Fmodes.fits  | idealized PIAACMC mirror 1 Fourier modes


#### Focal plane mask:

Focal plane mask design defined by :
- [s] Sectors flag (0: no sectors, 1: sectors), variable PIAACMC_FPMsectors
- [r] Resolved target flag (0: point source, 1: resolved source)
- [mr] Mask radius in units of 0.1 l/D
- [rrr] number of rings, variable piaacmc[0].NBrings
- [zzz] number of zones

Output file	| Notes
----------------|-------------------------------------
fpmzmap[s]_[rrr]_[zzz].fits | Zones map
fpm_zonea[r][s]_[mr]_[rrr]_[zzz].fits | amplitude for each zone
fpm_zonea[r][s]_[mr]_[rrr]_[zzz].fits | thickness for each zone


IDEALIZED OR PHYSICAL MASK\n

Idealized mask is a single zone mask with thickness adjusted for lambda/2 phase shift and a (non-physical) partial transmission.
Physical mask consist of 1 or more zones with full transmission. Each zone can have a different thickness.

By default, a non-physical mask is first created with transmission piaacmc[0].fpmaskamptransm read from piaacmcparams.conf.
Computations indices using a physical mask:
- set transmission to 1.0:  piaacmc[0].fpmaskamptransm = 1.0.
- set focal plane mask radius to larger value: piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio * PIAACMC_MASKRADLD


#### Lyot stops:

Output file	| Notes
----------------|-------------------------------------
./piaaconfxxx/LyotStop0.fits	| Lyot Stop 0
./piaaconfxxx/LyotStop1.fits	| Lyot Stop 1


#### Amplitude & Phase at planes:

Files are /piaaconfxxx/WFamp_nnn.fits and WFpha_nnn.fits, where nnn is the plane index.\n
Complex amplitude is shown AFTER the element has been applied, in the plane of the element.\n

Plane index	| description
----------------|-------------------------------------
000	|	Input pupil
001	|	Fold mirror used to induce pointing offsets
002	|	PIAA M0
003	|	PIAA M1
004	|	PIAAM1 edge opaque mask
005	|	post-focal plane mask pupil
006	|	Lyot Stop 0
007	|	Lyot Stop 1
008	|	invPIAA1
009	|	invPIAA0
010	|	back end mask


#### Performance Evaluation:

Plane index	| description
----------------|-------------------------------------
./piaaconfxxx/scoringmask0.fits | Evaluation points in focal plane, hardcoded in PIAACMCsimul_computePSF()
./piaaconfxxx/CnormFactor.txt | PSF normalization factor used to compute contrast
./piaaconfxxx/flux.txt	| total intensity at each plane






## STEP 002: Specify input pupil geometry

The pupil geometry is copied to file ./piaacmcconf[nnn]/pupa0_[size].fits




## STEP 003 (mode = 0): compute on-axis PSF for new pupil geometry




## STEP 004 (mode = 5): Compute Lyot stops shapes and locations, 1st pass

For circular centrally obscured pupils, the default Lyot stops configuration consists of two stops: one that masks the outer part of the beam, and one that masks the central obstruction.\n

Fine optimization of the stops locations is done with a separate command. The optional lsoptrange value is the range (unit: m) for the mask position search (default = 2 m).\n


Example result (for 1 l/D mask, 2048 array size):
\verbatim
Peak constrast (rough estimate)= 1.78874e-06
Total light in scoring field = 7.75534e-06  -> Average contrast = 5.62228e-08
\endverbatim



## STEP 005 (mode = 2): Optimize focal plane mask transmission, 1st pass



## STEP 006 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput



## STEP 007 (mode = 40): Tune PIAA shapes and focal plane mask transm, 10 cosine modes, 5 Fourier modes

This takes approximately 20mn for size = 1024.\n
Progress can be tracked by watching file :
\verbatim
tail -f linoptval.txt
\endverbatim


## STEP 008 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes


## STEP 009 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput



## STEP 010 (mode = 1): Tune Lyot stops conjugations



## STEP 011 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes



## STEP 012 (mode = 40): Tune PIAA shapes and focal plane mask transm,  40 cosine modes, 150 Fourier modes

The total number of free parameters is 380 = (40+150)*2, so this routine takes a long time to complete (hours).


## STEP 013 (mode = 5): Compute Lyot stops shapes and locations, 3nd pass, 70% throughput






