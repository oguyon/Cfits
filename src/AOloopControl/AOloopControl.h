#ifndef _AOLOOPCONTROL_H
#define _AOLOOPCONTROL_H



#define AOLOOPCONTROL_CONF_FNAME "/tmp/aoloopcontrol.conf.shm"



typedef struct
{
  struct timespec tnow;  // computed at time of sending DM commands
  double time_sec; // converted in second
  
  // SETUP
  int init; // has been initialized
  unsigned long long cnt;
  unsigned long long DMupdatecnt;
  int kill; // set to 1 to kill computation loop

  char name[80];
  
  long ID_WFS;    // Camera input
  char WFSname[80];   
  long sizexWFS;
  long sizeyWFS;
  long sizeWFS;
  long long WFScnt;  

  long ID_WFS1; // averaged, dark-subtracted frame (float), name = "imWFS1_%ld", loop;
  long ID_WFS2; // WFS1 image minus WFS reference (float), name = "imWFS2_%ld", loop;
  long ID_refWFS; // WFS reference

  long ID_DM;     // DM
  char DMname[80];
  long ID_DMRM;     // DM channel for RM
  char DMnameRM[80];
   long sizexDM;
  long sizeyDM;
  long sizeDM;

  int init_refWFS;    // WFS reference image loaded 
  int init_RM;        // Response Matrix loaded
  int init_CM;        // Control Matrix loaded

  
  long ID_DMmodes;  // DM modes
  long NBDMmodes;  
  float maxlimit; // maximum absolute value for mode values
  

  long ID_cmd_modes; // modal commands to be applied to DM
  long ID_cmd1_modes; // modal values at output of control matrix multiplication

  char respMname[80];
  long ID_respM;  // Response matrix
  char contrMname[80];
  long ID_contrM; // Control matrix
  

  
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

  int GPU; // 1 if computation done by GPU

  // LOOP TELEMETRY
  double RMSmodes;

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


int init_AOloopControl();
int AOloopControl_camimage_extract2D_sharedmem_loop(char *in_name, char *out_name, long size_x, long size_y, long xstart, long ystart);
int compute_ControlMatrix(long loop, long NB_MODE_REMOVED, char *ID_Rmatrix_name, char *ID_Cmatrix_name, char *ID_VTmatrix_name);
int AOloopControl_InitializeMemory();
int Average_cam_frames(long loop, long NbAve);
long AOloopControl_MakeDMModes(long loop, long NBmodes, char *IDname);
long AOloopControl_loadCM(long loop, char *CMfname);
int AOloopControl_loadconfigure(long loop, char *config_fname);
int set_DM_modes(long loop);
int Measure_Resp_Matrix(long loop, long NbAve, float amp, long nbloop, long fDelay, long NBiter);
int ControlMatrixMultiply( float *cm_array, float *imarray, long m, long n, float *outvect);
int AOcompute(long loop);
int AOloopControl_run();

int AOloopControl_showparams(long loop);
int AOloopControl_setLoopNumber(long loop);
int AOloopControl_loopkill();
int AOloopControl_loopon();
int AOloopControl_loopoff();
int AOloopControl_logon();
int AOloopControl_logoff();
int AOloopControl_setgain(float gain);
int AOloopControl_setmaxlimit(float maxlimit);
int AOloopControl_setframesAve(long nbframes);


#endif
