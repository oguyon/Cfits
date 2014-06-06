#ifndef _AOSYSTSIM_H
#define _AOSYSTSIM_H



#define DISPCOMB_FILENAME_CONF "/tmp/dmdispcombconf.conf.shm"
#define DMTURBCONF_FILENAME "/tmp/dmturb.conf.shm"



typedef struct
{
  int ON;

  long loopcnt;
  long updatecnt;
  
  int busy; // if set to 1, hold off and wait
  float MAXVOLT; // maximum voltage on DM

  int status;

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




int init_AOsystSim();
int SCExAO_DM_disp2V(long IDdisp, long IDvolt);
int SCEXAO_DM_createconf();
int SCEXAO_DM_loadconf();
int SCEXAO_DM_unloadconf();
int SCExAO_DM_CombineChannels(int mode);
int SCExAO_DM_dmdispcombstatus();
int SCExAO_DM_dmdispcomboff();
int SCExAO_DM_dmtrigoff();

int SCEXAO_DMturb_createconf();
int SCEXAO_DMturb_loadconf();
int SCExAO_DM_dmturboff();
int SCExAO_DM_dmturb_wspeed(double wspeed);
int SCExAO_DM_dmturb_ampl(double ampl);
int SCExAO_DM_dmturb_LOcoeff(double LOcoeff);
int SCExAO_DM_dmturb_tint(long tint);
int SCExAO_DM_dmturb_printstatus();
int SCExAO_DM_turb();


#endif
