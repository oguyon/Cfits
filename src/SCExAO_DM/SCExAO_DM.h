#ifndef _AOSYSTSIM_H
#define _AOSYSTSIM_H



#define DISPCOMB_FILENAME_CONF "/tmp/dmdispcombconf.conf.shm"


int init_AOsystSim();

int AOsystSim_run();



typedef struct
{
  int ON;

  long loopcnt;
  long updatecnt;
  
  int busy; // if set to 1, hold off and wait
  float MAXVOLT; // maximum voltage on DM

  struct timespec tstart;
  struct timespec tend;

  long moninterval; // [us]  
} SCEXAO_DISPCOMB_CONF;


#endif
