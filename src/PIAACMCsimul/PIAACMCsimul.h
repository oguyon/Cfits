#ifndef _PIAACMCSIMUL_H
#define _PIAACMCSIMUL_H




//
// *****************************************************************************************************
// -------------------------- structure defining a reflective PIAACMC system ---------------------------
// *****************************************************************************************************

//
// this structure holds parameters to be optimized in the design
//
typedef struct {
  float beamrad; // [m]
  long size;
  float pixscale; // [m/pix]

  // PIAA mirror shapes
  
  long CmodesID; // Cosine radial mode
  long NBCmodes;
  long FmodesID; // Fourier 2D modes
  long NBFmodes;

  long piaa0CmodesID;
  long piaa0FmodesID;
  long piaa1CmodesID;
  long piaa1FmodesID;
  
  
  
} MIRRORPIAACMCDESIGN;




int init_PIAACMCsimul();


void PIAACMCsimul_init( MIRRORPIAACMCDESIGN *design, long index );
void  PIAACMCsimul_free( void );

int PIAACMCsimul_mkPIAAMshapes_from_RadSag(char *fname, double radius_edge, char *ID_PIAAM1_name, char *ID_PIAAM2_name);
int PIAACMCsimul_run();

#endif
