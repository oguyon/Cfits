#ifndef _AOSYSTSIM_H
#define _AOSYSTSIM_H



#define DISPCOMB_FILENAME_CONF "/tmp/dmdispcombconf.conf.shm"
#define DMTURBCONF_FILENAME "/tmp/dmturb.conf.shm"

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



typedef struct
{
  int on;
  long cnt;


  double wspeed; // wind speed [m/s]
  double ampl; // [um RMS]
  double LOcoeff; // 0 for full correction of low orders, 1 for no correction

  long tint; // interval between consecutive DM updates [us]


  double simtime;

  struct timespec tstart;
  struct timespec tend;

} SCEXAO_DMTURBCONF;


#endif
