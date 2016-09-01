# PIAACMC design

Tools for design of Phase-Induced Amplitude Apodization Complex Mask Coronagraph (PIAACMC).

Design scripts are in `./script/`  

Documentation files are in `./doc/`

For detailed decumentations and step-by-step PIAACMC design instructions, see: `./doc/PIAACMCsimul.html`

# Description

- Uses Fresnel propagation engine ( OptSystProp.c OptSystProp.h ) between optical elements
- Computes linear perturbations around current design for optimization
- Automatically computes multiple Lyot stop to follow variable conjugation in different parts of the beam
- Polychromatic propagations
- Fits aspheric PIAA shapes on basis of radial cosines and 2-D Fourier modes



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






\subsection step002 3.002. STEP 002: Specify input pupil geometry

The pupil geometry is copied to file ./piaacmcconf[nnn]/pupa0_[size].fits




\subsection step003 3.003. STEP 003 (mode = 0): compute on-axis PSF for new pupil geometry




\subsection step004 3.004. STEP 004 (mode = 5): Compute Lyot stops shapes and locations, 1st pass

For circular centrally obscured pupils, the default Lyot stops configuration consists of two stops: one that masks the outer part of the beam, and one that masks the central obstruction.\n

Fine optimization of the stops locations is done with a separate command. The optional lsoptrange value is the range (unit: m) for the mask position search (default = 2 m).\n


Example result (for 1 l/D mask, 2048 array size):
\verbatim
Peak constrast (rough estimate)= 1.78874e-06
Total light in scoring field = 7.75534e-06  -> Average contrast = 5.62228e-08
\endverbatim



\subsection step005 3.005. STEP 005 (mode = 2): Optimize focal plane mask transmission, 1st pass



\subsection step006 3.006. STEP 006 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput



\subsection step007: 3.007. STEP 007 (mode = 40): Tune PIAA shapes and focal plane mask transm, 10 cosine modes, 5 Fourier modes

This takes approximately 20mn for size = 1024.\n
Progress can be tracked by watching file :
\verbatim
tail -f linoptval.txt
\endverbatim

\subsection step008: 3.008. STEP 008 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes


\subsection step009: 3.009. STEP 009 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput



\subsection step010: 3.010. STEP 010 (mode = 1): Tune Lyot stops conjugations



\subsection step011: 3.011. STEP 011 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes



\subsection step012: 3.012. STEP 012 (mode = 40): Tune PIAA shapes and focal plane mask transm,  40 cosine modes, 150 Fourier modes

The total number of free parameters is 380 = (40+150)*2, so this routine takes a long time to complete (hours).


\subsection step013: 3.013. STEP 013 (mode = 5): Compute Lyot stops shapes and locations, 3nd pass, 70% throughput






