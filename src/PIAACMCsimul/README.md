# PIAACMC coronagraph design 


## Overview

Diffraction-based PIAACMC simulation / optimization
- Uses Fresnel propagations engine ( OptSystProp.c OptSystProp.h ) between optical elements
- Computes linear perturbations around current design for optimization
- Automatically computes multiple Lyot stop to follow variable conjugation in different parts of the beam
- Polychromatic propagations
- Fits aspheric PIAA shapes on basis of radial cosines and 2-D Fourier modes


## Design Steps

### STEP 1: Create an idealized centrally obscured apodized PIAACMC monochromatic design

This is meant as a starting point for the PIAACMC, which will then be optimized further\n

-# Create the 2D apodization function for the centrally obscured aperture\n
The following set of commands can be used, and builds up the apodization function iteratively\n
\verbatim
DFTZFACTOR=2
PNBITER=20
cormk2Dprolateld 1.2 100.0 0.2 apo 512
resizeim apo apostart 1024 1024
rm apo
DFTZFACTOR=4
PNBITER=10
cormk2Dprolateld 1.2 200.0 0.2 apo 1024
mv apo apostart
DFTZFACTOR=8
PNBITER=5
cormk2Dprolateld 1.2 200.0 0.2 apo 1024
savefits apo "!apo2D.fits"
\endverbatim

-# specify input pupil geometry\n
Input to design is a pupil map, FITS file, size x size pixel, valued 1.0 inside the pupil and 0.0 outside. The file name should be pupa_<size>.fits (for example pupa_1024.fits).


### Starting point for design: centrally obscured circular aperture


### Create piaacmcparams.conf file



