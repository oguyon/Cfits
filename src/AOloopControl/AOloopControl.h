#ifndef _AOLOOPCONTROL_H
#define _AOLOOPCONTROL_H



#define AOLOOPCONTROL_CONF_FNAME "/tmp/aoloopcontrol.conf.shm"

#define maxNBMB 100


typedef struct
{
  struct timespec tnow;  // computed at time of sending DM commands
  double time_sec; // converted in second
  
  // SETUP
  int init; // has been initialized
  unsigned long long cnt;
  unsigned long long cntmax;
  unsigned long long DMupdatecnt;
  int kill; // set to 1 to kill computation loop

  char name[80];
  
  char WFSname[80];   
  float DarkLevel;
  long sizexWFS;
  long sizeyWFS;
  long sizeWFS;
  long long WFScnt;  

  char DMname[80];
  char DMnameRM[80];
  long sizexDM;
  long sizeyDM;
  long sizeDM;

  int init_refWFS;    // WFS reference image loaded 
  int init_RM;        // Response Matrix loaded
  int init_CM;        // Control Matrix loaded

  
  long NBDMmodes;  
  float maxlimit; // maximum absolute value for mode values
  

  char respMname[80];
  char contrMname[80];
  

  

  
  // LOOP CONTROL
  int on;  // goes to 1 when loop starts, put to 0 to turn loop off
  float gain; // overall loop gain
  long framesAve; // number of frames to average

  int status;
  // -1: PROGRAM NOT RUNNING
  // 0: LOOP OFF
  // 1: WAITING FOR WFS IMAGE
  // 2: REMOVING DARK
  // 3: NORMALIZING WFS IMAGE
  // 4: REMOVING REF
  // 5: MULTIPLYING BY CONTROL MATRIX -> MODE VALUES
  // 6: MULTIPLYING BY GAINS
  // 7: MULTIPLYING BY MODE MATRIX -> COMMANDS SENT TO DM
  
  
  // LOOP TUNING
  // BLOCKS OF MODES
  long NBMblocks; // number of mode blocks
  long indexmaxMB[maxNBMB];
  float gainMB[maxNBMB];
  float limitMB[maxNBMB];
  float multfMB[maxNBMB];



  int GPU; // 1 if computation done by GPU

  // LOOP TELEMETRY
  double RMSmodes;
  double RMSmodesCumul;
  long long RMSmodesCumulcnt;

  // logs
  char logdir[80];
  int logon; // 1 if log is on, 0 if off
  long logsize;  // # of entries per log
  long IDlog0;  // image identifyer for log file #1
  long IDlog1;  // image identifyer for log file #2
  int logcnt; // current position in log 
  int logfnb; // current log file number (0 or 1)
  char userLOGstring[80];
  long timeorigin_sec;

} AOLOOPCONTROL_CONF;

long AOloopControl_makeTemplateAOloopconf(long loopnb);
int init_AOloopControl();
long AOloopControl_mkModes(char *ID_name, long msize, float CPAmax, float deltaCPA, double xc, double yx, double r0, double r1);
int AOloopControl_camimage_extract2D_sharedmem_loop(char *in_name, char *out_name, long size_x, long size_y, long xstart, long ystart);
int compute_ControlMatrix(long loop, long NB_MODE_REMOVED, char *ID_Rmatrix_name, char *ID_Cmatrix_name, char *ID_VTmatrix_name, double Beta);
int AOloopControl_InitializeMemory();
int Average_cam_frames(long loop, long NbAve);
long AOloopControl_MakeDMModes(long loop, long NBmodes, char *IDname);
long AOloopControl_loadCM(long loop, char *CMfname);
int AOloopControl_loadconfigure(long loop, char *config_fname, int mode);
int set_DM_modes(long loop);
int Measure_Resp_Matrix(long loop, long NbAve, float amp, long nbloop, long fDelay, long NBiter);
int ControlMatrixMultiply( float *cm_array, float *imarray, long m, long n, float *outvect);
int AOcompute(long loop);
int AOloopControl_run();

int AOloopControl_printloopstatus(long loop, long nbcol);
int AOloopControl_loopMonitor(long loop, double frequ, long nbcol);
int AOloopControl_statusStats();
int AOloopControl_showparams(long loop);

int AOloopControl_setLoopNumber(long loop);
int AOloopControl_loopkill();
int AOloopControl_loopon();
int AOloopControl_loopoff();
int AOloopControl_logon();
int AOloopControl_loopstep(long loop, long NBstep);
int AOloopControl_logoff();
int AOloopControl_loopreset();
int AOloopControl_setgain(float gain);
int AOloopControl_setmaxlimit(float maxlimit);
int AOloopControl_setframesAve(long nbframes);
int AOloopControl_setgainrange(long m0, long m1, float gainval);
int AOloopControl_setlimitrange(long m0, long m1, float limval);
int AOloopControl_setmultfrange(long m0, long m1, float multfval);
int AOloopControl_setgainblock(long mb, float gainval);
int AOloopControl_setlimitblock(long mb, float limitval);
int AOloopControl_setmultfblock(long mb, float multfval);
int AOloopControl_resetRMSperf();
int AOloopControl_scanGainBlock(long NBblock, long NBstep, float gainStart, float gainEnd, long NBgain);
int AOloopControl_InjectMode( long index, float ampl );
int AOloopControl_AutoTune();

#endif
