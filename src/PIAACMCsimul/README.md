# Overview

Tools for design of Phase-Induced Amplitude Apodization Complex Mask Coronagraph (PIAACMC) and Apodized Pupil Lyot Complex Mask Coronagraphs (APLCMC).

Design scripts are in `./script/`  

Documentation files are in `./doc/`

Example files are in `./examples/`

For detailed decumentations and step-by-step PIAACMC design instructions, see: `./doc/PIAACMCsimul.html`


# Description

- Uses Fresnel propagation engine ( OptSystProp.c OptSystProp.h ) between optical elements
- Computes linear perturbations around current design for optimization
- Automatically computes multiple Lyot stop to follow variable conjugation in different parts of the beam
- Polychromatic propagations
- Fits aspheric PIAA shapes on basis of radial cosines and 2-D Fourier modes


