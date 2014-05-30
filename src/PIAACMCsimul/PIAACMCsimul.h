#ifndef _PIAACMCSIMUL_H
#define _PIAACMCSIMUL_H



#define ApoFitCosFact 1.0


//
// *****************************************************************************************************
// -------------------------- structure defining a reflective PIAACMC system ---------------------------
// *****************************************************************************************************



//
// this structure holds parameters to be optimized in the PIAACMC diffractive design
//
typedef struct {


  // ======= SEED RADIAL PIAACMC PARAMETERS ======

  double centObs0; // input central obstruction
  double centObs1; // output central obstruction
  double r0lim; // outer radius after extrapolation, piaa mirror 0
  double r1lim; // outer radius after extrapolation, piaa mirror 1
  long NBradpts; // number of points for common r0, r1, piaa sags 1D table
 

  // Wavelength
  int nblambda;



  // ====== Overall OPTICAL Geometry ===============

  float beamrad; // [m]
  long size;
  float pixscale; // [m/pix]
  float piaasep;// separation between PIAA surfaces [m]

  // ========= LYOT STOPS ============
  long NBLyotStop;
  long IDLyotStop[10];
  double LyotStop_zpos[10];

  // ======= Optics shapes modes ============

  long CmodesID; // Cosine radial mode
  long Cmsize; // cosine modes size 
  long NBCmodes;

  long FmodesID; // Fourier 2D modes
  long Fmsize; 
  long NBFmodes;

  long piaa0CmodesID;
  long piaa0FmodesID;
  long piaa1CmodesID;
  long piaa1FmodesID;
  


  // ========= Focal Plane Mask =============

  long focmNBzone; // number of zones
  double Fratio; // beam Fratio at focal plane
  double maskscale; // m/pixel
  long zonezID;  // focm zone material thickness, double precision image
  long zoneaID;  // focm zone amplitude transmission, double precision image
  double fpzfactor; // focal plane mask DFT zoom factor

  double fpmRad; // outer radius
  long NBrings; // number of rings
  long fpmarraysize;
  int fpmmaterial; // materials:  1: SiO2  2: Si  3: PMGI  4: PMMA

} MIRRORPIAACMCDESIGN;






// module initialization
int init_PIAACMCsimul();
void  PIAACMCsimul_free( void );

// Focal plane mask
long PIAACMCsimul_mkFPM_zonemap(char *IDname);
long PIAACMCsimul_mkFocalPlaneMask(char *IDzonemap_name, char *ID_name,  int mode);

// initializes the optsyst structure to simulate reflective PIAACMC system
void PIAACMCsimul_init( MIRRORPIAACMCDESIGN *design, long index, double TTxld, double TTyld);

// PIAA optics (geometrical optics) tools
int PIAACMCsimul_loadRadialApodization(char *IDapo_name, long beamradpix, char *IDapofit_name);
int PIAACMCsimul_mkPIAAMshapes_from_RadSag(char *fname, char *ID_PIAAM0_name, char *ID_PIAAM1_name);


// misc optimization tools
int PIAACMCsimul_unwrapPhase();


int PIAACMCsimul_run();


#endif
