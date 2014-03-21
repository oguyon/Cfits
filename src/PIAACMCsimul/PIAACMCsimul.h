#ifndef _PIAACMCSIMUL_H
#define _PIAACMCSIMUL_H


int init_PIAACMCsimul();


void PIAACMCsimul_init( void );
void  PIAACMCsimul_free( void );

int PIAACMCsimul_mkPIAAMshapes_from_RadSag(char *fname, double radius_edge, char *ID_PIAAM1_name, char *ID_PIAAM2_name);

int PIAACMCsimul_run();

#endif
