#ifndef _OPTSYSTPROP_H
#define _OPTSYSTPROP_H


// ************************************************************************
// ------------------- DEFINITION OF OPTICAL ELEMENTS ---------------------
// ************************************************************************

//
// -------- DM -----------
// square grid geom only
//
typedef struct {
  long NBact1D;
  double pitch; // [m]
  double maxstroke; // max deviation; actuator moves from -maxstroke to +maxstroke [m]
  long dispID; // points to displacement matrix

  long IF_ID; // points to influence function map
  double IFpixscale; // influence function pixel scale [m]
  double IFsize; // influence function map size (linear)
  double pixscale; // map pixel scale [m/pix]
  long dispmapID; // points to displacement map    

} DM_SIM;



//
// --------- aspheric surface mirror -------------- 
// using same sampling as nominal beam
//
typedef struct {
  long surfID; // surface Z sag
} ASPHSURFM;


//
// ----------- aspheric surface, refractive ------------ 
// using same sampling as nominal beam
//
typedef struct {
  long mat0; // material before surface
  long mat1; // material after surface
  long surfID; // surface Z sag
} ASPHSURFR;




// ------- Focal plane mask ------------
typedef struct {
  long fpmID; // 1-focal plane mask complex amplitude cube
  double zfactor; // oversampling factor  
  int mode; // 1 if 1-CA, 0 if CA
} FOCMASK;


//
// *****************************************************************************************************
// ------------------------------ structure defining an optical system ---------------------------------
// *****************************************************************************************************

// Fresnel propagation used for simulations
// All elements except focal plane mask are placed in equivalent collimated beam space
// Input pupil is adopted as reference
//
// optical elements are applied sequentially, and consist of amplitude and phase cubes (1 slize per lambda)
// types of elements :
// 1: static opaque mask, defined by image identifyer 
// 2: DM, defined by index in array of DM structures
// 3: reflective aspheric surface. Defined by OPD map (used for mirrors)
// 4: refractive aspheric surface. Defined by OPD map, material type (index of refraction) before and after surface
// 5: focal plane mask 
//
typedef struct {

  int nblambda;
  double lambdaarray[2000];

  double beamrad; // beam radius at input in collimated space [m]
  double pixscale; // pixel scale in collimated beam [m/pix]
  long size; // array size
  long DFTgridpad; // 0 for full res DFT, >0 for faster coarser DFTs
  
  // =============== OPTICAL ELEMENTS ===================
  
  // predefined arrays that can be used for optical elements
  long NB_DM; // max number of DMs
  DM_SIM DMarray[4];

  long NB_asphsurfm; // max number of aspheric mirrors
  ASPHSURFM ASPHSURFMarray[4];

  long NB_asphsurfr; // max number of aspheric refractive surfaces
  ASPHSURFM ASPHSURFRarray[4];

  long NB_focmask; // max number of focal plane masks
  FOCMASK FOCMASKarray[4];
  
 
  
  long NBelem; // number of optical elements
  int elemtype[100]; // element type
  int elemarrayindex[100]; // if element is DM or aspheric surface, this is the index in the corresponding array of elements, otherwise, this is the image index
  double flux[100]; // total flux AFTER element
  double elemZpos[100]; // position along beam

  // this is what is used for propagations, created from info above
  long elem_amp_ID_array[100]; // amplitude map identifyer, multiplicative
  long elem_pha_ID_array[100]; // phase map identifyer, additive


} OPTSYST;



int OptSystProp_propagateCube(OPTSYST *optsyst, long index, char *IDin_amp_name, char *IDin_pha_name, char *IDout_amp_name, char *IDout_pha_name, double zprop);

int OptSystProp_run(OPTSYST *optsyst, long index, long elemstart, long elemend, char *savedir);


#endif
