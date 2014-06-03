#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <err.h>
#include <fcntl.h>
#include <sched.h>
#include <ncurses.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>




#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "linopt_imtools/linopt_imtools.h"
#include "AOloopControl/AOloopControl.h"
#include "image_filter/image_filter.h"

#include "ZernikePolyn/ZernikePolyn.h"



int wcol, wrow; // window size


long aoconfID_WFS = -1;
long aoconfID_WFS1 = -1;
long aoconfID_WFS2 = -1;
long aoconfID_refWFS = -1;
long aoconfID_DM = -1;
long aoconfID_DMRM = -1;
long aoconfID_DMmodes = -1;


// Fourier Modes
long aoconfID_cmd_modes = -1;
long aoconfID_cmd1_modes = -1;
long aoconfID_RMS_modes = -1;
long aoconfID_AVE_modes = -1;
long aoconfID_GAIN_modes = -1;
long aoconfID_LIMIT_modes = -1;
long aoconfID_MULTF_modes = -1;




long aoconfID_cmd_modesRM = -1;

long aoconfID_respM = -1;
long aoconfID_contrM = -1;

long aoconfIDlog0 = -1;
long aoconfIDlog1 = -1;



extern DATA data;


#define NB_AOloopcontrol 10 // max number of loops
long LOOPNUMBER = 0; // current loop index
int AOloopcontrol_meminit = 0;
int AOlooploadconf_init = 0;

#define AOconfname "/tmp/AOconf.shm"
AOLOOPCONTROL_CONF *AOconf; // configuration - this can be an array
int *AOconf_loaded = 0;
int *AOconf_fd; 

float *arrayftmp;
unsigned short *arrayutmp;
int avcamarraysInit = 0;


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//




int AOloopControl_mkModes_cli()
{ 
  if(CLI_checkarg(1,3)+CLI_checkarg(2,2)+CLI_checkarg(3,1)+CLI_checkarg(4,1)+CLI_checkarg(5,1)+CLI_checkarg(6,1)+CLI_checkarg(7,1)+CLI_checkarg(8,1)==0)
    AOloopControl_mkModes(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numf, data.cmdargtoken[5].val.numf, data.cmdargtoken[6].val.numf, data.cmdargtoken[7].val.numf, data.cmdargtoken[8].val.numf);
  else
    return 1;      
}

int AOloopControl_camimage_extract2D_sharedmem_loop_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,2)+CLI_checkarg(4,2)+CLI_checkarg(5,2)+CLI_checkarg(6,2)==0)
    {
      AOloopControl_camimage_extract2D_sharedmem_loop(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string , data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl, data.cmdargtoken[5].val.numl, data.cmdargtoken[6].val.numl);
      return 0;
    }
  else
    return 1;
}


int AOloopControl_loadconfigure_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,3)==0)
    {
      AOloopControl_loadconfigure(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.string, 0);
      return 0;
    }
  else
    return 1;
}

int AOloopControl_setLoopNumber_cli()
{
  if(CLI_checkarg(1,2)==0)
    {
      AOloopControl_setLoopNumber(data.cmdargtoken[1].val.numl);
      return 0;
    }
  else
    return 1;
}


int AOloopControl_setgain_cli()
{
  if(CLI_checkarg(1,1)==0)
    {
      AOloopControl_setgain(data.cmdargtoken[1].val.numf);
      return 0;
    }
  else
    return 1;
}

int AOloopControl_setmaxlimit_cli()
{
  if(CLI_checkarg(1,1)==0)
    {
      AOloopControl_setmaxlimit(data.cmdargtoken[1].val.numf);
      return 0;
    }
  else
    return 1;
}

int AOloopControl_Measure_Resp_Matrix_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,1)+CLI_checkarg(3,2)+CLI_checkarg(4,2)+CLI_checkarg(5,2)==0)
    {
      Measure_Resp_Matrix(LOOPNUMBER, data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl, data.cmdargtoken[5].val.numl);
      return 0;
    }
  else
    return 1;
}




int AOloopControl_setframesAve_cli()
{
  if(CLI_checkarg(1,2)==0)
    {
      AOloopControl_setframesAve(data.cmdargtoken[1].val.numl);
      return 0;
    }
  else
    return 1;
}


int AOloopControl_computeCM_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,4)+CLI_checkarg(3,3)==0)
    {
      compute_ControlMatrix(LOOPNUMBER, data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string, "evecM");
      save_fits("evecM","!evecM.fits");
      delete_image_ID("evecM");
    }
  else
    return 1;
}


int AOloopControl_loadCM_cli()
{
  if(CLI_checkarg(1,3)==0)
    {
      AOloopControl_loadCM(LOOPNUMBER, data.cmdargtoken[1].val.string);
      return 0;
    }
  else
    return 1;
}




int AOloopControl_loopstep_cli()
{
  if(CLI_checkarg(1,2)==0)
    {
      AOloopControl_loopstep(LOOPNUMBER, data.cmdargtoken[1].val.numl);
      return 0;
    }
  else
    return 1;
}




int AOloopControl_loopMonitor_cli()
{
 if(CLI_checkarg(1,1)+CLI_checkarg(2,2)==0)
   {
     AOloopControl_loopMonitor(LOOPNUMBER, data.cmdargtoken[1].val.numf, data.cmdargtoken[2].val.numl);
     return 0;
   }
 else
   {
     AOloopControl_loopMonitor(LOOPNUMBER, 10.0, 3);
     return 0;
   }
}


int AOloopControl_setgainrange_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
    {
      AOloopControl_setgainrange(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numf);
      return 0;
    }
  else
    return 1; 
}


int AOloopControl_setlimitrange_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
    {
      AOloopControl_setlimitrange(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numf);
      return 0;
    }
  else
    return 1; 
}


int AOloopControl_setmultfrange_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
    {
      AOloopControl_setmultfrange(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numf);
      return 0;
    }
  else
    return 1; 
}



int AOloopControl_setgainblock_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
    {
      AOloopControl_setgainblock(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf);
      return 0;
    }
  else
    return 1; 
}


int AOloopControl_setlimitblock_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
    {
      AOloopControl_setlimitblock(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf);
      return 0;
    }
  else
    return 1; 
}


int AOloopControl_setmultfblock_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
    {
      AOloopControl_setmultfblock(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf);
      return 0;
    }
  else
    return 1; 
}





int init_AOloopControl()
{
  FILE *fp;
  int r;

  if((fp=fopen("loopnb.txt","r"))!=NULL)
    {
      r = fscanf(fp,"%ld", &LOOPNUMBER);
      printf("LOOP NUMBER = %ld\n", LOOPNUMBER);
      fclose(fp);
    }
  else
    LOOPNUMBER = 0;
  

  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "AO loop control");
  data.NBmodule++;


  strcpy(data.cmd[data.NBcmd].key,"aolmkmodes");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_mkModes_cli;
  strcpy(data.cmd[data.NBcmd].info,"make control modes");
  strcpy(data.cmd[data.NBcmd].syntax,"<output modes> <size> <max CPA> <delta CPA> <cx> <cy> <r0> <r1>");
  strcpy(data.cmd[data.NBcmd].example,"aolmkmodes modes 50 5.0 0.8");
  strcpy(data.cmd[data.NBcmd].Ccall,"long AOloopControl_mkModes(char *ID_name, long msize, float CPAmax, float deltaCPA, double xc, double yx, double r0, double r1)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"cropshim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_camimage_extract2D_sharedmem_loop_cli;
  strcpy(data.cmd[data.NBcmd].info,"crop shared mem image");
  strcpy(data.cmd[data.NBcmd].syntax,"<input image> <output image>");
  strcpy(data.cmd[data.NBcmd].example,"cropshim imin imout");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_camimage_extract2D_sharedmem_loop(char *in_name, char *out_name, long size_x, long size_y, long xstart, long ystart)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolloadconf");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loadconfigure_cli;
  strcpy(data.cmd[data.NBcmd].info,"load AO loop configuration from file");
  strcpy(data.cmd[data.NBcmd].syntax,"<loop #> <conf file>");
  strcpy(data.cmd[data.NBcmd].example,"AOlooploadconf 1 aoloop1.conf");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loadconfigure(long loopnb, char *fname, 0)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolacqresp");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_Measure_Resp_Matrix_cli;
  strcpy(data.cmd[data.NBcmd].info,"acquire AO response matrix and WFS reference");
  strcpy(data.cmd[data.NBcmd].syntax,"<ave# [long]> <ampl [float]> <nbloop [long]> <frameDelay [long]> <NBiter [long]>");
  strcpy(data.cmd[data.NBcmd].example,"aolacqresp 50 0.1 5 2");
  strcpy(data.cmd[data.NBcmd].Ccall,"int Measure_Resp_Matrix(long loop, long NbAve, float amp, long nbloop, long fDelay, long NBiter)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolrun");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_run;
  strcpy(data.cmd[data.NBcmd].info,"run AO loop");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"AOlooprun");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_run()");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolkill");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopkill;
  strcpy(data.cmd[data.NBcmd].info,"kill AO loop");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aolkill");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setLoopNumber()");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolnb");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setLoopNumber_cli;
  strcpy(data.cmd[data.NBcmd].info,"set AO loop #");
  strcpy(data.cmd[data.NBcmd].syntax,"<loop nb>");
  strcpy(data.cmd[data.NBcmd].example,"AOloopnb 0");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setLoopNumber(long loop)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolon");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopon;
  strcpy(data.cmd[data.NBcmd].info,"turn loop on");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aolon");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loopon()");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aolstep");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopstep_cli;
  strcpy(data.cmd[data.NBcmd].info,"turn loop on for N steps");
  strcpy(data.cmd[data.NBcmd].syntax,"<nbstep>");
  strcpy(data.cmd[data.NBcmd].example,"aolstep");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loopstep(long loop, long NBstep)");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aoloff");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopoff;
  strcpy(data.cmd[data.NBcmd].info,"turn loop off");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aoloff");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loopoff()");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aolreset");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopreset;
  strcpy(data.cmd[data.NBcmd].info,"reset loop, and turn it off");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aolreset");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loopreset()");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aollogon");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_logon;
  strcpy(data.cmd[data.NBcmd].info,"turn log on");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aollogon");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_logon()");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aollogoff");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_logoff;
  strcpy(data.cmd[data.NBcmd].info,"turn log off");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aollogoff");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_logoff()");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aolsetgain");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setgain_cli;
  strcpy(data.cmd[data.NBcmd].info,"set gain");
  strcpy(data.cmd[data.NBcmd].syntax,"<gain value>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetgain 0.1");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setgain(float gain)");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aolsetmaxlim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setmaxlimit_cli;
  strcpy(data.cmd[data.NBcmd].info,"set max limit for AO mode correction");
  strcpy(data.cmd[data.NBcmd].syntax,"<limit value>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetmaxlim 0.01");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setmaxlimit(float maxlimit)");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aolsetnbfr");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setframesAve_cli;
  strcpy(data.cmd[data.NBcmd].info,"set number of frames to be averaged");
  strcpy(data.cmd[data.NBcmd].syntax,"<nb frames>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetnbfr 10");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setframesAve(long nbframes)");
  data.NBcmd++;
  

  strcpy(data.cmd[data.NBcmd].key,"aolcmmake");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_computeCM_cli;
  strcpy(data.cmd[data.NBcmd].info,"make control matrix");
  strcpy(data.cmd[data.NBcmd].syntax,"<NBmodes removed> <RespMatrix> <ContrMatrix>");
  strcpy(data.cmd[data.NBcmd].example,"aolcmmake 8 respm cmat");
  strcpy(data.cmd[data.NBcmd].Ccall,"int compute_ControlMatrix(long loop, long NB_MODE_REMOVED, char *ID_Rmatrix_name, char *ID_Cmatrix_name, char *ID_VTmatrix_name)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolloadcm");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loadCM_cli;
  strcpy(data.cmd[data.NBcmd].info,"load new control matrix from file");
  strcpy(data.cmd[data.NBcmd].syntax,"<fname>");
  strcpy(data.cmd[data.NBcmd].example,"aolloadcm cm32.fits");
  strcpy(data.cmd[data.NBcmd].Ccall,"long AOloopControl_loadCM(long loop, char *CMfname)");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aolmon");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopMonitor_cli;
  strcpy(data.cmd[data.NBcmd].info,"monitor loop");
  strcpy(data.cmd[data.NBcmd].syntax,"<frequ> <Nbcols>");
  strcpy(data.cmd[data.NBcmd].example,"aolmon 10.0 3");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loopMonitor(long loop, double frequ)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"aolsetgainr");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setgainrange_cli;
  strcpy(data.cmd[data.NBcmd].info,"set modal gains");
  strcpy(data.cmd[data.NBcmd].syntax,"<modemin [long]> <modemax [long]> <gainval>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetgainr 20 30 0.2");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setgainrange(long m0, long m1, float gainval)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolsetlimitr");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setlimitrange_cli;
  strcpy(data.cmd[data.NBcmd].info,"set modal limits");
  strcpy(data.cmd[data.NBcmd].syntax,"<modemin [long]> <modemax [long]> <limval>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetlimitr 20 30 0.02");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setlimitrange(long m0, long m1, float gainval)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolsetmultfr");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setmultfrange_cli;
  strcpy(data.cmd[data.NBcmd].info,"set modal multf");
  strcpy(data.cmd[data.NBcmd].syntax,"<modemin [long]> <modemax [long]> <multfval>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetmultfr 10 30 0.98");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setmultfrange(long m0, long m1, float multfval)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolsetgainb");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setgainblock_cli;
  strcpy(data.cmd[data.NBcmd].info,"set modal gains by block");
  strcpy(data.cmd[data.NBcmd].syntax,"<modemin [long]> <modemax [long]> <gainval>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetgainb 2 0.2");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setgainblock(long m0, long m1, float gainval)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolsetlimitb");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setlimitblock_cli;
  strcpy(data.cmd[data.NBcmd].info,"set modal limits by block");
  strcpy(data.cmd[data.NBcmd].syntax,"<modemin [long]> <modemax [long]> <limval>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetlimitb 2 0.02");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setlimitblock(long m0, long m1, float gainval)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolsetmultfb");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setmultfblock_cli;
  strcpy(data.cmd[data.NBcmd].info,"set modal multf by block");
  strcpy(data.cmd[data.NBcmd].syntax,"<modemin [long]> <modemax [long]> <multfval>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetmultfb 2 0.98");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setmultfblock(long m0, long m1, float multfval)");
  data.NBcmd++;



  // add atexit functions here
  // atexit((void*) SCEXAO_DM_unloadconf);
  
  return 0;
}




long AOloopControl_mkModes(char *ID_name, long msize, float CPAmax, float deltaCPA, double xc, double yc, double r0, double r1)
{
  long ID0;
  long ID;
  long k, ii, jj;

  long IDmask;

  double ave;
  double offset;
  double totm;
  double totvm;
  
  double a0=0.88;
  double b0=40.0;

  double a1=1.2;
  double b1=12.0;

  double x, y, r;
  double val0, val1, rms;
  
  long IDtm, IDem;
  long IDeModes;

  long kelim = 20;
  double coeff;
  long citer;
  long NBciter = 200;
  long IDg;

  // if Mmask exists, use it, otherwise create it
  //
  IDmask = image_ID("Mmask");
  if(IDmask==-1)
    {
      IDmask = create_2Dimage_ID("Mmask", msize, msize);
      for(ii=0;ii<msize;ii++)
	for(jj=0;jj<msize;jj++)
	  {
	    x = 1.0*ii-xc;
	    y = 1.0*jj-yc;
	    r = sqrt(x*x+y*y)/r1;
	    val1 = 1.0-exp(-pow(a1*r,b1));
	    r = sqrt(x*x+y*y)/r0;
	    val0 = exp(-pow(a0*r,b0));
	    data.image[IDmask].array.F[jj*msize+ii] = val0*val1;
	  }
      //      save_fits("Mmask", "!Mmask.fits");    
    }
      
  totm = arith_image_total("Mmask");
  msize = data.image[IDmask].md[0].size[0];



  




  linopt_imtools_makeCPAmodes("CPAmodes", msize, CPAmax, deltaCPA, 0.5*msize, 1.2, 0);
  ID0 = image_ID("CPAmodes");
  list_image_ID();
  printf("  %ld %ld %ld\n", msize, msize, data.image[ID0].md[0].size[2]-1 );
  ID = create_3Dimage_ID(ID_name, msize, msize, data.image[ID0].md[0].size[2]-1);
  
  for(k=0;k<data.image[ID0].md[0].size[2]-1;k++)
    {

      // Remove excluded modes
      IDeModes = image_ID("emodes");
      if(IDeModes!=-1)
	{
	  IDtm = create_2Dimage_ID("tmpmode", msize, msize);
	  
	  

	  for(ii=0;ii<msize*msize;ii++)
	    data.image[IDtm].array.F[ii] = data.image[ID0].array.F[(k+1)*msize*msize+ii];
	  linopt_imtools_image_fitModes("tmpmode", "emodes", "Mmask", 1.0e-5, "lcoeff", 0);
	  linopt_imtools_image_construct("emodes", "lcoeff", "em00");
	  delete_image_ID("lcoeff");
	  IDem = image_ID("em00");
	  
	  coeff = 1.0-exp(-pow(1.0*k/kelim,6.0));
	  if(k>2.0*kelim)
	    coeff = 1.0;
	  for(ii=0;ii<msize*msize;ii++)
	    {
	      data.image[ID].array.F[k*msize*msize+ii] = data.image[IDtm].array.F[ii] - coeff*data.image[IDem].array.F[ii];	  	  
	    }
	  
	  
	  delete_image_ID("em00");
	  delete_image_ID("tmpmode");
	} 
      
      ave = 0.0;
      totvm = 0.0;
      for(ii=0;ii<msize*msize;ii++)
	{	  
	  //	  data.image[ID].array.F[k*msize*msize+ii] = data.image[ID0].array.F[(k+1)*msize*msize+ii];	  
	  totvm += data.image[ID].array.F[k*msize*msize+ii]*data.image[IDmask].array.F[ii];
	}
      offset = totvm/totm;
      
      for(ii=0;ii<msize*msize;ii++)
	{
	  data.image[ID].array.F[k*msize*msize+ii] -= offset;
	  data.image[ID].array.F[k*msize*msize+ii] *= data.image[IDmask].array.F[ii];	
	}

      offset = 0.0;
      for(ii=0;ii<msize*msize;ii++)
	offset += data.image[ID].array.F[k*msize*msize+ii];
      
      rms = 0.0;
      for(ii=0;ii<msize*msize;ii++)
	{
	  data.image[ID].array.F[k*msize*msize+ii] -= offset/msize/msize;
	  rms += data.image[ID].array.F[k*msize*msize+ii]*data.image[ID].array.F[k*msize*msize+ii];
	}
      rms = sqrt(rms/totm);

      for(ii=0;ii<msize*msize;ii++)
	data.image[ID].array.F[k*msize*msize+ii] /= sqrt(rms);	  
    }

  for(citer=0;citer<NBciter;citer++)
    {
      printf("Convolution [%3ld/%3ld]\n", citer, NBciter);
      gauss_filter(ID_name, "modeg", 4.0*(NBciter-citer)/NBciter, 5);
      IDg = image_ID("modeg");
      for(k=0;k<data.image[ID].md[0].size[2];k++)
	{
	  for(ii=0;ii<msize*msize;ii++)
	    if(data.image[IDmask].array.F[ii]<0.98)
	      data.image[ID].array.F[k*msize*msize+ii] = data.image[IDg].array.F[k*msize*msize+ii];
	}
      delete_image_ID("modeg");
    }
    

  return(ID);
}



//
// every time im_name changes (counter increments), crop it to out_name in shared memory
//
int AOloopControl_camimage_extract2D_sharedmem_loop(char *in_name, char *out_name, long size_x, long size_y, long xstart, long ystart)
{
  long iiin,jjin, iiout, jjout;
  long IDin, IDout;
  int atype;
  long *sizeout;
  long long cnt0;
  
  sizeout = (long*) malloc(sizeof(long)*2);
  sizeout[0] = size_x;
  sizeout[1] = size_y;

  IDin = image_ID(in_name);
  atype = data.image[IDin].md[0].atype;

  // Create shared memory output image 
  IDout = create_image_ID(out_name, 2, sizeout, atype, 1, 0);
  
  cnt0 = -1;

  switch (atype) {
  case USHORT :
    while(1)
      {
	usleep(100);
	if(data.image[IDin].md[0].cnt0!=cnt0)
	  {
	    cnt0 = data.image[IDin].md[0].cnt0;
	    for(iiout=0; iiout<size_x; iiout++)
	      for(jjout=0; jjout<size_y; jjout++)
		{
		  iiin = xstart + iiout;
		  jjin = ystart + jjout;
		  data.image[IDout].array.U[jjout*size_x+iiout] = data.image[IDin].array.U[jjin*data.image[IDin].md[0].size[0]+iiin];
		}
	    data.image[IDout].md[0].cnt0 = cnt0;
	  }
      }
    break;
  case FLOAT :
    while(1)
      {
	usleep(100);
	if(data.image[IDin].md[0].cnt0!=cnt0)
	  {
	    cnt0 = data.image[IDin].md[0].cnt0;
	    for(iiout=0; iiout<size_x; iiout++)
	      for(jjout=0; jjout<size_y; jjout++)
		{
		  iiin = xstart + iiout;
		  jjin = ystart + jjout;
		  data.image[IDout].array.F[jjout*size_x+iiout] = data.image[IDin].array.F[jjin*data.image[IDin].md[0].size[0]+iiin];
		}
	    data.image[IDout].md[0].cnt0 = cnt0;
	  }
      }
    break;
  default :
     printf("ERROR: DATA TYPE NOT SUPPORTED\n");
      exit(0);
    break;
  }
  free(sizeout);

  return(0);
}


//
// Computes control matrix
// Conventions:
//   m: number of actuators (= NB_MODES)
//   n: number of sensors  (= # of pixels)
//
int compute_ControlMatrix(long loop, long NB_MODE_REMOVED, char *ID_Rmatrix_name, char *ID_Cmatrix_name, char *ID_VTmatrix_name) /* works even for m != n */
{
  FILE *fp;
  long ii1, jj1, k, ii;
  gsl_matrix *matrix_D; /* this is the response matrix */
  gsl_matrix *matrix_Ds; /* this is the pseudo inverse of D */ 
  gsl_matrix *matrix_Dtra;
  gsl_matrix *matrix_DtraD;
  gsl_matrix *matrix_DtraDinv;
  gsl_matrix *matrix_DtraD_evec;
  gsl_matrix *matrix1;
  gsl_matrix *matrix2;
  gsl_vector *matrix_DtraD_eval;
  gsl_eigen_symmv_workspace *w;
 
  gsl_matrix *matrix_save;

  long m;
  long n;
  long ID_Rmatrix, ID_Cmatrix, ID_VTmatrix;
  long *arraysizetmp;


  long IDeigenmodes;


  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();


  arraysizetmp = (long*) malloc(sizeof(long)*3);


  ID_Rmatrix = image_ID(ID_Rmatrix_name);

  n = data.image[ID_Rmatrix].md[0].size[0]*data.image[ID_Rmatrix].md[0].size[1]; //AOconf[loop].NBDMmodes;
  m = data.image[ID_Rmatrix].md[0].size[2]; //AOconf[loop].sizeWFS;

  /* in this procedure, m=number of actuators/modes, n=number of WFS elements */
  //  long m = smao[0].NBmode;
  // long n = smao[0].NBwfselem;

  printf("m = %ld actuators (modes), n = %ld sensors\n", m, n);
  fflush(stdout);

  matrix_DtraD_eval = gsl_vector_alloc (m); 
  matrix_D = gsl_matrix_alloc (n,m);
  matrix_Ds = gsl_matrix_alloc (m,n);
  matrix_Dtra = gsl_matrix_alloc (m,n);
  matrix_DtraD = gsl_matrix_alloc (m,m); 
  matrix_DtraDinv = gsl_matrix_alloc (m,m); 
  matrix_DtraD_evec = gsl_matrix_alloc (m,m);
  


  /* write matrix_D */
  for(k=0;k<m;k++)
    for(ii=0;ii<n;ii++)
      gsl_matrix_set (matrix_D, ii, k, data.image[ID_Rmatrix].array.F[k*n+ii]);

  /* compute DtraD */
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, matrix_D, matrix_D, 0.0, matrix_DtraD);


  /* compute the inverse of DtraD */

  /* first, compute the eigenvalues and eigenvectors */
  w =   gsl_eigen_symmv_alloc (m);
  matrix_save = gsl_matrix_alloc (m,m);
  gsl_matrix_memcpy(matrix_save, matrix_DtraD);
  gsl_eigen_symmv (matrix_save, matrix_DtraD_eval, matrix_DtraD_evec, w);
  gsl_matrix_free(matrix_save);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort (matrix_DtraD_eval, matrix_DtraD_evec, GSL_EIGEN_SORT_ABS_DESC);

  printf("Eigenvalues\n");
  fflush(stdout);

  // Write eigenvalues
  if((fp=fopen("eigenv.dat","w"))==NULL)
    {
      printf("ERROR: cannot create file \"eigenv.dat\"\n");
      exit(0);
    }
  for(k=0; k<m; k++)
    fprintf(fp,"%ld %g\n", k, gsl_vector_get(matrix_DtraD_eval,k));
  fclose(fp);
  
  for(k=0; k<m; k++)
    printf("Mode %ld eigenvalue = %g\n", k, gsl_vector_get(matrix_DtraD_eval,k));
  
  
  // Write rotation matrix to go from DM modes to eigenmodes
  arraysizetmp[0] = m;
  arraysizetmp[1] = m;
  ID_VTmatrix = create_image_ID(ID_VTmatrix_name, 2, arraysizetmp, FLOAT, 0, 0);
  for(ii=0;ii<m;ii++) // modes
    for(k=0;k<m;k++) // modes
      data.image[ID_VTmatrix].array.F[k*m+ii] = (float) gsl_matrix_get( matrix_DtraD_evec, k, ii);
  

  
  /* second, build the "inverse" of the diagonal matrix of eigenvalues (matrix1) */
  matrix1 = gsl_matrix_alloc (m, m);
  for(ii1=0; ii1<m; ii1++)
    for(jj1=0; jj1<m; jj1++)
      {
	if(ii1==jj1)
	  {
	    if((m-ii1-1)<NB_MODE_REMOVED)
	      gsl_matrix_set(matrix1, ii1, jj1, 0.0);
	    else
	      gsl_matrix_set(matrix1, ii1, jj1, 1.0/gsl_vector_get(matrix_DtraD_eval,ii1));
	  }
	else
	  gsl_matrix_set(matrix1, ii1, jj1, 0.0);
      }
 
  printf("Compute inverse\n");
  fflush(stdout);

  /* third, compute the "inverse" of DtraD */
  matrix2 = gsl_matrix_alloc (m, m);
  //  printf("step 0\n");
  // fflush(stdout);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, matrix_DtraD_evec, matrix1, 0.0, matrix2);
  // printf("step 1\n");
  //fflush(stdout);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, matrix2, matrix_DtraD_evec, 0.0, matrix_DtraDinv);
  //printf("step 2\n");
  //fflush(stdout);
  gsl_matrix_free(matrix1);
  gsl_matrix_free(matrix2);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, matrix_DtraDinv, matrix_D, 0.0, matrix_Ds);

  printf("Write result\n");
  fflush(stdout);
  arraysizetmp[0] = AOconf[loop].sizexWFS;
  arraysizetmp[1] = AOconf[loop].sizexWFS;
  arraysizetmp[2] = m;

  ID_Cmatrix = create_image_ID(ID_Cmatrix_name, 3, arraysizetmp, FLOAT, 0, 0);


  /* write result */
  for(ii=0;ii<n;ii++) // sensors
    for(k=0;k<m;k++) // actuator modes
      data.image[ID_Cmatrix].array.F[k*n+ii] = (float) gsl_matrix_get(matrix_Ds, k, ii);
      

  gsl_vector_free(matrix_DtraD_eval);
  gsl_matrix_free(matrix_D);
  gsl_matrix_free(matrix_Ds);
  gsl_matrix_free(matrix_Dtra);
  gsl_matrix_free(matrix_DtraD);
  gsl_matrix_free(matrix_DtraDinv);
  gsl_matrix_free(matrix_DtraD_evec);

  free(arraysizetmp);

  return(ID_Cmatrix);
}










int AOloopControl_InitializeMemory()
{
  int SM_fd;
  struct stat file_stat;
  int create = 0;
  int result;
  long loop;

  SM_fd = open(AOconfname, O_RDWR);
  if(SM_fd==-1)
    {
      printf("Cannot import file -> creating file\n");
      create = 1;
    }
  else
    {
      fstat(SM_fd, &file_stat);
      printf("File %s size: %zd\n", AOconfname, file_stat.st_size);
      if(file_stat.st_size!=sizeof(AOLOOPCONTROL_CONF)*NB_AOloopcontrol)
	{
	  printf("File size wrong -> recreating file\n");
	  create = 1;
	  close(SM_fd);
	}
    }
  
  if(create==1)
    {
      SM_fd = open(AOconfname, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
	 
      if (SM_fd == -1) {
	perror("Error opening file for writing");
	exit(0);
      }

      result = lseek(SM_fd, sizeof(AOLOOPCONTROL_CONF)*NB_AOloopcontrol-1, SEEK_SET);
      if (result == -1) {
	close(SM_fd);
	perror("Error calling lseek() to 'stretch' the file");
	exit(0);
      }
      
      result = write(SM_fd, "", 1);
      if (result != 1) {
	close(SM_fd);
	perror("Error writing last byte of the file");
	exit(0);
      }
    }

  AOconf = (AOLOOPCONTROL_CONF*) mmap(0, sizeof(AOLOOPCONTROL_CONF)*NB_AOloopcontrol, PROT_READ | PROT_WRITE, MAP_SHARED, SM_fd, 0);
  if (AOconf == MAP_FAILED) {
    close(SM_fd);
    perror("Error mmapping the file");
    exit(0);
  }
  
  AOconf[loop].on = 0;
  AOconf[loop].cnt = 0;	  
  AOconf[loop].cntmax = 0;	  

  if(create==1)
    {
      for(loop=0; loop<NB_AOloopcontrol; loop++)
	{
	  AOconf[loop].init = 0;
	  AOconf[loop].on = 0;
	  AOconf[loop].cnt = 0;	  
	  AOconf[loop].cntmax = 0;	  
	  AOconf[loop].maxlimit = 0.3;
	  AOconf[loop].gain = 0.0;
	  AOconf[loop].framesAve = 1;
	}
    }
  else
    {
      for(loop=0; loop<NB_AOloopcontrol; loop++)
	if(AOconf[loop].init == 1)
	  {
	    printf("LIST OF ACTIVE LOOPS:\n");
	    printf("----- Loop %ld   (%s) ----------\n", loop, AOconf[loop].name);
	    printf("  WFS:  %s  [%ld]  %ld x %ld\n", AOconf[loop].WFSname, aoconfID_WFS, AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS);
	    printf("   DM:  %s  [%ld]  %ld x %ld\n", AOconf[loop].DMname, aoconfID_DM, AOconf[loop].sizexDM, AOconf[loop].sizeyDM);
	    printf("DM RM:  %s  [%ld]  %ld x %ld\n", AOconf[loop].DMnameRM, aoconfID_DM, AOconf[loop].sizexDM, AOconf[loop].sizeyDM);
	  }
    }
  
  AOloopcontrol_meminit = 1;


  return 0;
}




//
// supports ring buffer
//
int Average_cam_frames(long loop, long NbAve)
{
  long imcnt;
  long ii;
  double total;
  char name[200];
  int atype;
  long slice;
  char *ptrv;
 

  atype = data.image[aoconfID_WFS].md[0].atype;

  if(avcamarraysInit==0)
    {
      arrayftmp = (float*) malloc(sizeof(float)*AOconf[loop].sizeWFS);
      arrayutmp = (unsigned short*) malloc(sizeof(unsigned short)*AOconf[loop].sizeWFS);
      avcamarraysInit = 1;
    }
  
  if(NbAve>1)
    for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
      data.image[aoconfID_WFS1].array.F[ii] = 0.0;

  if(data.image[aoconfID_WFS].md[0].naxis==2) // single buffer
    {
      switch (atype) {
      case FLOAT :
	imcnt = 0;
	while(imcnt<NbAve)
	  {
	    usleep(50);	  
	    if(data.image[aoconfID_WFS].md[0].write == 0)
	      {
		if(AOconf[loop].WFScnt!=data.image[aoconfID_WFS].md[0].cnt0)
		  {	      
		    memcpy (arrayftmp, data.image[aoconfID_WFS].array.F, sizeof(float)*AOconf[loop].sizeWFS);
		    AOconf[loop].WFScnt = data.image[aoconfID_WFS].md[0].cnt0;
		    if(NbAve>1)
		      {
			
			for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
			  data.image[aoconfID_WFS1].array.F[ii] += arrayftmp[ii];
		      }
		    else
		      memcpy(data.image[aoconfID_WFS1].array.F, arrayftmp,  sizeof(float)*AOconf[loop].sizeWFS);
		    imcnt++;
		  }      
	      }
	  }
	break;
      case USHORT :
	imcnt = 0;
	while(imcnt<NbAve)
	  {
	    usleep(50);
	    if(data.image[aoconfID_WFS].md[0].write == 0)
	      {
		if(AOconf[loop].WFScnt!=data.image[aoconfID_WFS].md[0].cnt0)
		  {
		    memcpy (arrayutmp, data.image[aoconfID_WFS].array.U, sizeof(unsigned short)*AOconf[loop].sizeWFS);
		    AOconf[loop].WFScnt = data.image[aoconfID_WFS].md[0].cnt0;
		    if(NbAve>1)
		      {
			for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
			  data.image[aoconfID_WFS1].array.F[ii] += arrayutmp[ii]; 
		      }
		    else
		      {
			for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
			  data.image[aoconfID_WFS1].array.F[ii] = arrayutmp[ii]; 
		      }
		    imcnt++;
		  }
	      }
	  }         
	break;
      default :
	printf("ERROR: DATA TYPE NOT SUPPORTED\n");
	exit(0);
	break;
      }
    }
  else // ring buffer mode, only works with NbAve = 1
    {
      while(AOconf[loop].WFScnt==data.image[aoconfID_WFS].md[0].cnt0) // new frame exists
	{
	  usleep(50);
	  // do nothing, wait
	}
      slice = data.image[aoconfID_WFS].md[0].cnt1-1;
      if(slice==-1)
	slice = data.image[aoconfID_WFS].md[0].size[2]-1;
      
      //      printf("READING SLICE %ld\n", slice);
      switch (atype) {
      case FLOAT :
	ptrv = (char*) data.image[aoconfID_WFS].array.F;
	ptrv += sizeof(float)*slice* AOconf[loop].sizeWFS;
	memcpy(data.image[aoconfID_WFS1].array.F, ptrv,  sizeof(float)*AOconf[loop].sizeWFS);	
	break;
      case USHORT :
	ptrv = (char*) data.image[aoconfID_WFS].array.U;
	ptrv += sizeof(unsigned short)*slice* AOconf[loop].sizeWFS;
       	memcpy (arrayutmp, ptrv, sizeof(unsigned short)*AOconf[loop].sizeWFS);
	for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
	  data.image[aoconfID_WFS1].array.F[ii] = (float) arrayutmp[ii]; 
	break;
      default :
	printf("ERROR: DATA TYPE NOT SUPPORTED\n");
	exit(0);
	break;
      }

      AOconf[loop].WFScnt = data.image[aoconfID_WFS].md[0].cnt0;
	
    }
       
     
  
  // Normalize  
  sprintf(name, "imWFS1_%ld", loop);
  total = arith_image_total(data.image[aoconfID_WFS1].md[0].name);
     
     
  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
    data.image[aoconfID_WFS1].array.F[ii] /= total;
  
  total = arith_image_total(data.image[aoconfID_WFS1].md[0].name);

  /*
  printf("SLEEPING ... ");
  fflush(stdout);
  sleep(4);
  printf("done\n");
  fflush(stdout);
  */


  // save_fits(data.image[aoconfID_WFS].md[0].name, "!testim.fits");
  //sleep(5);


  return(0);
}









long AOloopControl_MakeDMModes(long loop, long NBmodes, char *IDname)
{
  long ID;
  long IDtmp;
  long ii, jj;
  double x, y;
  long size = AOconf[loop].sizexDM;
  long m;
  float rpix;
  long size2;

  size2 = size*size;
  rpix = 0.5*size;

  ID = create_3Dimage_ID_float(IDname, size, size, NBmodes);
  
  for(m=0;m<NBmodes;m++)
    {
      IDtmp = mk_zer("zertmp", size, m+1, rpix);
      for(ii=0;ii<size2;ii++)
	data.image[ID].array.F[size2*m+ii] = data.image[IDtmp].array.F[ii];
      delete_image_ID("zertmp");
    }

  return(ID);
}




long AOloopControl_loadCM(long loop, char *CMfname)
{
  long ID = -1;
  int vOK;
  char name[200];
  long ID0;
  long ii;

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  if( (ID=load_fits(CMfname, "tmpcontrM")) != -1 )
    {
      // check size is OK
      vOK = 1;
      if(data.image[ID].md[0].naxis!=3)
	{
	  printf("Control matrix has wrong dimension\n");
	  vOK = 0;
	}
      if(data.image[ID].md[0].atype!=FLOAT)
	{
	  printf("Control matrix has wrong type\n");
	  vOK = 0;
	}
      if(vOK==1)
	{
	  if(data.image[ID].md[0].size[0]!=AOconf[loop].sizexWFS)
	    {
	      printf("Control matrix has wrong x size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexWFS);
	      vOK = 0;
	    }
	  if(data.image[ID].md[0].size[1]!=AOconf[loop].sizeyWFS)
	    {
	      printf("Control matrix has wrong y size\n");
	      vOK = 0;
	    }
	  if(data.image[ID].md[0].size[2]!=AOconf[loop].NBDMmodes)
	    {
	      printf("Control matrix has wrong z size\n");
	      vOK = 0;
	    }	  
	}
      if(vOK==1)
	{
	  AOconf[loop].init_CM = 1;
	  sprintf(name, "ContrM_%ld", loop);
	  ID = image_ID(name);
	  if(ID==-1)
	    ID = read_sharedmem_image(name);
	  ID0 = image_ID("tmpcontrM");
	  for(ii=0;ii<AOconf[loop].sizexWFS*AOconf[loop].sizeyWFS*AOconf[loop].NBDMmodes;ii++)
	    data.image[ID].array.F[ii] = data.image[ID0].array.F[ii];
	}
      delete_image_ID("tmpcontrM");
    }
 
  return(ID);
}


//
// read loop configuration file
// mode = 0: read
// mode = 1: create
int AOloopControl_loadconfigure(long loop, char *config_fname, int mode)
{
  FILE *fp;
  char content[200];
  char name[200];
  char fname[200];
  long ID;
  long *sizearray;
  int vOK;
  int kw;
  long k;


  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();



 
  
  


  if(mode==1)
    {
 


      sizearray = (long*) malloc(sizeof(long)*3);
      
      
      if(read_config_parameter(config_fname, "loopname", content)==0)
	exit(0);
      printf("loop name : %s\n", content);
      fflush(stdout);
      strcpy(AOconf[loop].name, content);
      
     

      if(read_config_parameter(config_fname, "GPU", content)==0)
	exit(0);
      printf("GPU : %d\n", atoi(content));
      fflush(stdout);
      AOconf[loop].GPU = atoi(content);
      
      
  
      
      // Connect to WFS camera
      if(read_config_parameter(config_fname, "WFSname", content)==0)
	exit(0);
      printf("WFS file name : %s\n", content);
      fflush(stdout);
      strcpy(AOconf[loop].WFSname, content);
      aoconfID_WFS = read_sharedmem_image(content);
      AOconf[loop].sizexWFS = data.image[aoconfID_WFS].md[0].size[0];
      AOconf[loop].sizeyWFS = data.image[aoconfID_WFS].md[0].size[1];
      //  AOconf[loop].sizezWFS = data.image[aoconfID_WFS].md[0].size[2];
      AOconf[loop].sizeWFS = AOconf[loop].sizexWFS*AOconf[loop].sizeyWFS;
      	
      
      
      // Create WFS1 image memory (averaged, dark-subtracted)
      sprintf(name, "imWFS1_%ld", loop);
      sizearray[0] =  AOconf[loop].sizexWFS;
      sizearray[1] =  AOconf[loop].sizeyWFS;
      aoconfID_WFS1 = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
      
      
      
      
      // Create WFS1 image memory (averaged, dark-subtracted)
      sprintf(name, "imWFS2_%ld", loop);
      sizearray[0] =  AOconf[loop].sizexWFS;
      sizearray[1] =  AOconf[loop].sizeyWFS;
      aoconfID_WFS2 = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
            
      
      
      // Create WFSref image memory (averaged, dark-subtracted)
      sprintf(name, "refWFS_%ld", loop);
      sizearray[0] =  AOconf[loop].sizexWFS;
      sizearray[1] =  AOconf[loop].sizeyWFS;
      aoconfID_refWFS = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
            
      
      if(read_config_parameter(config_fname, "WFSrefim", content)==0)
	exit(0);
      printf("WFS ref file name : %s\n", content);
      vOK = 0;
      AOconf[loop].init_refWFS = 0;
      if(is_fits_file(content)==1)
	if((ID=load_fits(content, "tmprefwfs"))!=-1)
	  {
	    printf("Verifying file\n");
	    vOK = 1;
	    if(data.image[ID].md[0].naxis!=2)
	      {
		printf("refWFS has wrong dimension\n");
		vOK = 0;
	      }
	    if(data.image[ID].md[0].atype!=FLOAT)
	      {
		printf("refWFS has wrong type\n");
		vOK = 0;
	      }
	    if(vOK==1)
	      {
		if(data.image[ID].md[0].size[0]!=AOconf[loop].sizexWFS)
		  {
		    printf("refWFS has wrong x size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexWFS);
		    vOK = 0;
		  }
		if(data.image[ID].md[0].size[1]!=AOconf[loop].sizeyWFS)
		  {
		    printf("refWFS has wrong y size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexWFS);
		    vOK = 0;
		  }
	      }
	    if(vOK==1)
	      {
		sprintf(name, "refWFS_%ld", loop);
		sizearray[0] =  AOconf[loop].sizexWFS;
		sizearray[1] =  AOconf[loop].sizeyWFS;
		
		//    aoconfID_refWFS = read_sharedmem_image(name);
		//if(aoconfID_refWFS == -1)
		if(mode==1)
		  aoconfID_refWFS = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
		else
		  aoconfID_refWFS = read_sharedmem_image(name);
		

		copy_image_ID("tmprefwfs", name);
		AOconf[loop].init_refWFS = 1;
	      }
	    delete_image_ID("tmprefwfs");
	  }  
      
      
      // Connect to DM 
      if(read_config_parameter(config_fname, "DMname", content)==0)
	exit(0);
      printf("DM file name : %s\n", content);
      strcpy(AOconf[loop].DMname, content);
      aoconfID_DM = read_sharedmem_image(content);
      AOconf[loop].sizexDM = data.image[aoconfID_DM].md[0].size[0];
      AOconf[loop].sizeyDM = data.image[aoconfID_DM].md[0].size[1];
      AOconf[loop].sizeDM = AOconf[loop].sizexDM*AOconf[loop].sizeyDM;
      
      if(read_config_parameter(config_fname, "DMnameRM", content)==0)
	exit(0);
      printf("DM file name RM: %s\n", content);
      strcpy(AOconf[loop].DMnameRM, content);
      aoconfID_DMRM = read_sharedmem_image(content);
      



      // Load DM modes (will exit if not successful)
      sprintf(fname,"DMmodes_%ld.fits", loop);
      sprintf(name,"DMmodes_%ld", loop);
      if(read_config_parameter(config_fname, "DMmodes", content)==0)
	exit(0);
      printf("DM modes file name : %s\n", content);
      
      vOK = 0;
      if(is_fits_file(content)==1)
	if((ID=load_fits(content, "tmpdmmodes"))!=-1)
	  {
	    printf("Verifying file\n");
	    vOK = 1;
	    if(data.image[ID].md[0].naxis != 3)
	      {
		printf("DM modes has wrong dimension\n");
		vOK = 0;
	      }
	    if(data.image[ID].md[0].atype != FLOAT)
	      {
		printf("DM modes has wrong type\n");
		vOK = 0;
	      }
	    if(vOK==1)
	      {
		if(data.image[ID].md[0].size[0]!=AOconf[loop].sizexDM)
		  {
		    printf("DM modes has wrong x size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexDM);
		    vOK = 0;
		  }
		if(data.image[ID].md[0].size[1]!=AOconf[loop].sizeyDM)
		  {
		    printf("DM modes has wrong y size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexDM);
		    vOK = 0;
		  }
	      }
	    if(vOK==1)
	      {
		AOconf[loop].NBDMmodes = data.image[ID].md[0].size[2];
		sprintf(name, "DMmodes_%ld", loop);
		sizearray[0] =  AOconf[loop].sizexDM;
		sizearray[1] =  AOconf[loop].sizeyDM;
		sizearray[2] =  AOconf[loop].NBDMmodes;
		aoconfID_DMmodes = create_image_ID(name, 3, sizearray, FLOAT, 1, 0);
	
		copy_image_ID("tmpdmmodes", name);
	      }
	    delete_image_ID("tmpdmmodes");
	  }
      if(vOK == 0)
	{
	  printf("\n");
	  printf("========== ERROR: NEED DM MODES TO START AO LOOP ===========\n");	
	  printf("\n");
	  exit(0);
	}
      


      // modes blocks

      if(read_config_parameter(config_fname, "NBMblocks", content)==0)
	AOconf[loop].NBMblocks = 1;
      printf("GPU : %d\n", atoi(content));
      fflush(stdout);
      AOconf[loop].NBMblocks = atoi(content);
      
      if(AOconf[loop].NBMblocks==1)
	AOconf[loop].indexmaxMB[0] = AOconf[loop].NBDMmodes;
      else
	{
	  for(k=0;k<AOconf[loop].NBMblocks;k++)
	    AOconf[loop].indexmaxMB[k] = (long) (pow(1.0*(k+1.0)/AOconf[loop].NBMblocks,2.0)*AOconf[loop].NBDMmodes);
	  AOconf[loop].indexmaxMB[AOconf[loop].NBMblocks-1] = AOconf[loop].NBMblocks;
	}
      

      for(k=0;k<AOconf[loop].NBMblocks;k++)
	{
	  AOconf[loop].gainMB[k] = 1.0;
	  AOconf[loop].limitMB[k] = 1.0;
	  AOconf[loop].multfMB[k] = 1.0;
	}


      
      // Allocate / create logging data files/memory
      if(read_config_parameter(config_fname, "logdir", content)==0)
	{
	  printf("parameter logdir missing\n");
	  exit(0);
	}
      strcpy(AOconf[loop].logdir, content);
      if(read_config_parameter(config_fname, "logsize", content)==0)
	{
	  printf("parameter logsize missing\n");
	  exit(0);
	}
      AOconf[loop].logsize = atol(content);
      // time [s]       (1)
      // gains          ( AOconf[loop].NBDMmodes )
      // ID_cmd_modes   ( AOconf[loop].NBDMmodes )
      // ID_cmd1_modes  ( AOconf[loop].NBDMmodes )
      sizearray[0] = 1+3*AOconf[loop].NBDMmodes;
      sizearray[1] = AOconf[loop].logsize;
      sprintf(name, "loop%ldlog0", loop);
      aoconfIDlog0 = create_image_ID(name, 2, sizearray, FLOAT, 1, 10);
      ID = aoconfIDlog0;
      data.image[ID].md[0].NBkw = 1;
      kw = 0;
      strcpy(data.image[ID].kw[kw].name, "TIMEORIGIN");	
      data.image[ID].kw[kw].type = 'L';
      data.image[ID].kw[kw].value.numl = 0;
      strcpy(data.image[ID].kw[kw].comment, "time offset [sec]");
      
      sprintf(name, "loop%ldlog1", loop);
      aoconfIDlog1 = create_image_ID(name, 2, sizearray, FLOAT, 1, 10);
      ID = aoconfIDlog1;
      data.image[ID].md[0].NBkw = 1;
      kw = 0;
      strcpy(data.image[ID].kw[kw].name, "TIMEORIGIN");	
      data.image[ID].kw[kw].type = 'L';
      data.image[ID].kw[kw].value.numl = 0;
      strcpy(data.image[ID].kw[kw].comment, "time offset [sec]");
      
      
      
      AOconf[loop].logcnt = 0;
      AOconf[loop].logfnb = 0;
      strcpy(AOconf[loop].userLOGstring, "");
      
      // AOconf[loop].ID_DMmodes = AOloopControl_MakeDMModes(loop, 5, name);
      
      printf("%ld modes\n", AOconf[loop].NBDMmodes);
      







      // Create modal command vector memory
      sprintf(name, "DMmode_cmd_%ld", loop);
      sizearray[0] =  AOconf[loop].NBDMmodes;
      sizearray[1] =  1;
      aoconfID_cmd_modes = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
      
      
      sprintf(name, "DMmode_cmd1_%ld", loop);
      sizearray[0] =  AOconf[loop].NBDMmodes;
      sizearray[1] =  1;
      aoconfID_cmd1_modes = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
      
      sprintf(name, "DMmode_RMS_%ld", loop);
      sizearray[0] =  AOconf[loop].NBDMmodes;
      sizearray[1] =  1;
      aoconfID_RMS_modes = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
      for(k=0;k<AOconf[loop].NBDMmodes;k++)
	data.image[aoconfID_RMS_modes].array.F[k] = 0.0;
     
      sprintf(name, "DMmode_AVE_%ld", loop);
      sizearray[0] =  AOconf[loop].NBDMmodes;
      sizearray[1] =  1;
      aoconfID_AVE_modes = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
      for(k=0;k<AOconf[loop].NBDMmodes;k++)
	data.image[aoconfID_AVE_modes].array.F[k] = 0.0;

      sprintf(name, "DMmode_GAIN_%ld", loop);
      sizearray[0] =  AOconf[loop].NBDMmodes;
      sizearray[1] =  2;
      aoconfID_GAIN_modes = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
      for(k=0;k<AOconf[loop].NBDMmodes;k++)
	data.image[aoconfID_GAIN_modes].array.F[k] = 1.0;
      
      sprintf(name, "DMmode_LIMIT_%ld", loop);
      sizearray[0] =  AOconf[loop].NBDMmodes;
      sizearray[1] =  1;
      aoconfID_LIMIT_modes = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
      for(k=0;k<AOconf[loop].NBDMmodes;k++)
	data.image[aoconfID_LIMIT_modes].array.F[k] = 1.0;
 
      sprintf(name, "DMmode_MULTF_%ld", loop);
      sizearray[0] =  AOconf[loop].NBDMmodes;
      sizearray[1] =  1;
      aoconfID_MULTF_modes = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
      for(k=0;k<AOconf[loop].NBDMmodes;k++)
	data.image[aoconfID_MULTF_modes].array.F[k] = 1.0;
    

 





      // Create Response Matrix shared memory
      sprintf(name, "RespM_%ld", loop);
      sizearray[0] =  AOconf[loop].sizexWFS;
      sizearray[1] =  AOconf[loop].sizeyWFS;
      sizearray[2] =  AOconf[loop].NBDMmodes;
      aoconfID_respM = create_image_ID(name, 3, sizearray, FLOAT, 1, 0);
     
      
      // Load response matrix if possible
      AOconf[loop].init_RM = 0;  
      if(read_config_parameter(config_fname, "RespMatrix", content)==0)
	exit(0);
      printf("Response Matrix file name : %s\n", content);
      if((ID=load_fits(content, "tmprespM"))!=-1)
	{
	  // check size is vOK
	  vOK = 1;
	  if(data.image[ID].md[0].naxis!=3)
	    {
	      printf("Response matrix has wrong dimension\n");
	      vOK = 0;
	    }
	  if(data.image[ID].md[0].atype!=FLOAT)
	    {
	      printf("Response matrix has wrong type\n");
	      vOK = 0;
	    }
	  if(vOK==1)
	    {
	      if(data.image[ID].md[0].size[0]!=AOconf[loop].sizexWFS)
		{
		  printf("Response matrix has wrong x size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexWFS);
		  vOK = 0;
		}
	      if(data.image[ID].md[0].size[1]!=AOconf[loop].sizeyWFS)
		{
		  printf("Response matrix has wrong y size\n");
		  vOK = 0;
		}
	      if(data.image[ID].md[0].size[2]!=AOconf[loop].NBDMmodes)
		{
		  printf("Response matrix has wrong z size\n");
		  vOK = 0;
		}	  
	    }
	  if(vOK==1)
	    {
	      AOconf[loop].init_RM = 1;
	      sprintf(name, "RespM_%ld", loop);
	      copy_image_ID("tmprespM", name);
	    }
	  delete_image_ID("tmprespM");
	}
      
      

      
      // Create Control Matrix shared memory
      sprintf(name, "ContrM_%ld", loop);
      sizearray[0] =  AOconf[loop].sizexWFS;
      sizearray[1] =  AOconf[loop].sizeyWFS;
      sizearray[2] =  AOconf[loop].NBDMmodes;
      aoconfID_contrM = create_image_ID(name, 3, sizearray, FLOAT, 1, 0);
      
      
      // Load control matrix if possible
      AOconf[loop].init_CM = 0;
      if(read_config_parameter(config_fname, "ContrMatrix", content)==0)
	exit(0);
      printf("Control Matrix file name : %s\n", content);
      AOloopControl_loadCM(loop, content);

  
      free(sizearray);

      AOconf[loop].init = 1;


      printf("   init_WFSref    %d\n", AOconf[loop].init_refWFS);
      printf("   init_RM        %d\n", AOconf[loop].init_RM);
      printf("   init_CM        %d\n", AOconf[loop].init_CM);
      
      AOconf[loop].init = 1;
    }
  else
    {
      aoconfID_WFS = read_sharedmem_image(AOconf[loop].WFSname);
    

      sprintf(name, "imWFS1_%ld", loop);
      aoconfID_WFS1 = read_sharedmem_image(name);
      
      sprintf(name, "imWFS2_%ld", loop);
      aoconfID_WFS2 = read_sharedmem_image(name);
      
      sprintf(name, "refWFS_%ld", loop);
      aoconfID_refWFS = read_sharedmem_image(name);

   
      aoconfID_DM = read_sharedmem_image(AOconf[loop].DMname);
      aoconfID_DMRM = read_sharedmem_image(AOconf[loop].DMnameRM);

	
      sprintf(name, "DMmodes_%ld", loop);
      aoconfID_DMmodes = read_sharedmem_image(name);
 
      sprintf(name, "loop%ldlog0", loop);
      aoconfIDlog0 = read_sharedmem_image(name);

      sprintf(name, "loop%ldlog1", loop);
      aoconfIDlog1 = read_sharedmem_image(name);


      sprintf(name, "DMmode_cmd_%ld", loop);
      aoconfID_cmd_modes = read_sharedmem_image(name);
  
      sprintf(name, "DMmode_cmd1_%ld", loop);
      aoconfID_cmd1_modes = read_sharedmem_image(name);

      sprintf(name, "DMmode_RMS_%ld", loop);
      aoconfID_RMS_modes = read_sharedmem_image(name);

      sprintf(name, "DMmode_AVE_%ld", loop);
      aoconfID_AVE_modes = read_sharedmem_image(name);

      sprintf(name, "DMmode_GAIN_%ld", loop);
      aoconfID_GAIN_modes = read_sharedmem_image(name);

      sprintf(name, "DMmode_LIMIT_%ld", loop);
      aoconfID_LIMIT_modes = read_sharedmem_image(name);

      sprintf(name, "DMmode_MULTF_%ld", loop);
      aoconfID_MULTF_modes = read_sharedmem_image(name);



      sprintf(name, "RespM_%ld", loop);
      aoconfID_respM = read_sharedmem_image(name); 

      sprintf(name, "ContrM_%ld", loop);
      aoconfID_contrM = read_sharedmem_image(name);
    }
  
  AOconf[loop].init = 1;

  return(0);
}




int set_DM_modes(long loop)
{
  long k;
  long i, j;
  float *arrayf;

  arrayf = (float*) malloc(sizeof(float)*AOconf[loop].sizeDM);

  for(j=0;j<AOconf[loop].sizeDM;j++)
    arrayf[j] = 0.0;

  for(k=0; k < AOconf[loop].NBDMmodes; k++)
    {
      for(i=0;i<AOconf[loop].sizeDM;i++)
	arrayf[i] += data.image[aoconfID_cmd_modes].array.F[k] * data.image[aoconfID_DMmodes].array.F[k*AOconf[loop].sizeDM+i];
    }

  data.image[aoconfID_DM].md[0].write = 1;
  memcpy (data.image[aoconfID_DM].array.F, arrayf, sizeof(float)*AOconf[loop].sizeDM);
  data.image[aoconfID_DM].md[0].cnt0++;
  data.image[aoconfID_DM].md[0].write = 0;

  free(arrayf);
  AOconf[loop].DMupdatecnt ++;

  return(0);
}


int set_DM_modesRM(long loop)
{
  long k;
  long i, j;
  float *arrayf;

  arrayf = (float*) malloc(sizeof(float)*AOconf[loop].sizeDM);



  for(j=0;j<AOconf[loop].sizeDM;j++)
    arrayf[j] = 0.0;

  for(k=0; k < AOconf[loop].NBDMmodes; k++)
    {
      for(i=0;i<AOconf[loop].sizeDM;i++)
	arrayf[i] += data.image[aoconfID_cmd_modesRM].array.F[k] * data.image[aoconfID_DMmodes].array.F[k*AOconf[loop].sizeDM+i];
    }


  data.image[aoconfID_DMRM].md[0].write = 1;
  memcpy (data.image[aoconfID_DMRM].array.F, arrayf, sizeof(float)*AOconf[loop].sizeDM);
  data.image[aoconfID_DMRM].md[0].cnt0++;
  data.image[aoconfID_DMRM].md[0].write = 0;

  free(arrayf);
  AOconf[loop].DMupdatecnt ++;


  return(0);
}




// measures response matrix AND reference 

int Measure_Resp_Matrix(long loop, long NbAve, float amp, long nbloop, long fDelay, long NBiter)
{
  long NBloops;
  long kloop;
  long delayus = 0; // delay in us
  long ii, i, imax;
  int Verbose = 1;
  long k1, k, k2;
  char fname[200];
  char name0[200];
  char name[200];
  
  long kk;
  long RespMatNBframes;
  long IDrmc;
  long kc;

  int recordCube = 1;
  long IDeigenmodes;

  long frameDelay = 0;
  long frameDelayMax = 50;
  long double RMsig;
  long double RMsigold;
  long kc0;
  FILE *fp;
  long NBexcl = 2; // number of frames excluded between DM mode changes
  long kc0min, kc0max;
  long IDrmtest;
  int vOK;


  long iter;
  long IDrmi;
  float beta = 0.0;
  float gain = 0.001;
  long IDrmcumul;
  long IDrefi;
  long IDrefcumul;

  long *sizearray;

  long IDrespM;
  long IDrefWFS;


  sizearray = (long*) malloc(sizeof(long)*3);


  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();


 

  if(fDelay==-1)
    {
      kc0min = 0;
      kc0max = frameDelayMax;
    }
  else
    {
      kc0min = fDelay;
      kc0max = fDelay+1;
    }


  RespMatNBframes = nbloop*2*AOconf[loop].NBDMmodes*NbAve;
  //  printf("%ld frames total\n");
  // fflush(stdout);

  if(recordCube == 1)
    IDrmc = create_3Dimage_ID("RMcube", AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS, RespMatNBframes+kc0max);



  
  printf("SETTING UP...\n");
  fflush(stdout);
  sprintf(fname, "AOloop%ld.conf", LOOPNUMBER);






  //  AOloopControl_loadconfigure(LOOPNUMBER, fname, 0);
  //exit(0);

  aoconfID_DMRM = read_sharedmem_image(AOconf[loop].DMnameRM);
  aoconfID_WFS = read_sharedmem_image(AOconf[loop].WFSname);

  sprintf(name, "imWFS1RM_%ld", loop);
  sizearray[0] = AOconf[loop].sizexWFS;
  sizearray[1] = AOconf[loop].sizeyWFS;
  aoconfID_WFS1 = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);





  
  sprintf(name, "DMmodes_%ld", loop);
  aoconfID_DMmodes = read_sharedmem_image(name);
  

  IDrmi = create_3Dimage_ID("RMiter", AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS, AOconf[loop].NBDMmodes);
  IDrmcumul = create_3Dimage_ID("RMcumul", AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS, AOconf[loop].NBDMmodes);

  IDrefi = create_2Dimage_ID("REFiter", AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS);
  IDrefcumul = create_2Dimage_ID("REFcumul", AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS);
 
  sizearray[0] =  AOconf[loop].NBDMmodes;
  printf("size = %ld\n", sizearray[0]);
  fflush(stdout);
  aoconfID_cmd_modesRM = create_image_ID("RMmodes", 1, sizearray, FLOAT, 1, 0);
  
  
  
  IDrespM = create_3Dimage_ID("respmacq", AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS, AOconf[loop].NBDMmodes);
  IDrefWFS = create_2Dimage_ID("refwfsacq", AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS);


 
  for(iter=0;iter<NBiter;iter++)
    {  
      NBloops = nbloop;
      
      // initialize RMiter to zero
      for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	for(k=0;k<AOconf[loop].NBDMmodes;k++)
	  data.image[IDrmi].array.F[k*AOconf[loop].sizeWFS+ii] = 0.0;
      
      
      // initialize reference to zero
      
      for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	data.image[IDrefi].array.F[ii] = 0.0;
      
      
      printf("\n");
      printf("Testing (in measure_resp_matrix function) :,  NBloops = %ld, NBmode = %ld\n",  NBloops, AOconf[loop].NBDMmodes);
      fflush(stdout);
      sleep(1);
      
      
      
      for(k2 = 0; k2 < AOconf[loop].NBDMmodes; k2++)
	data.image[aoconfID_cmd_modesRM].array.F[k2] = 0.0;
      
      

  
      kc = 0;
      for (kloop = 0; kloop < NBloops; kloop++)
	{
	  if(Verbose)
	    {
	      printf("\n Loop %ld / %ld (%f)\n", kloop, NBloops, amp);
	      fflush(stdout);
	    }
	  
	  for(k1 = 0; k1 < AOconf[loop].NBDMmodes; k1++)
	    {	 
	      for(k2 = 0; k2 < AOconf[loop].NBDMmodes; k2++)
		data.image[aoconfID_cmd_modesRM].array.F[k2] = 0.0;
	      
	      // positive 
	      data.image[aoconfID_cmd_modesRM].array.F[k1] = amp;
	      

	  
	      set_DM_modesRM(loop);
	      usleep(delayus);
	      
	
	      for(kk=0;kk<NbAve;kk++)
		{
		  Average_cam_frames(loop, 1);
		  
		  
		  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
		    {		      
		      data.image[IDrefi].array.F[ii] += data.image[aoconfID_WFS1].array.F[ii];	
		      data.image[IDrmc].array.F[kc*AOconf[loop].sizeWFS+ii] = data.image[aoconfID_WFS1].array.F[ii];
		    }
		  kc++;
		}



	      // negative
	      data.image[aoconfID_cmd_modesRM].array.F[k1] = 0.0-amp;
	      set_DM_modesRM(loop);
	      
	      usleep(delayus);
	      
	      for(kk=0;kk<NbAve;kk++)
		{
		  Average_cam_frames(loop, 1);
		  
		  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
		    {
		      data.image[IDrefi].array.F[ii] += data.image[aoconfID_WFS1].array.F[ii];
		      data.image[IDrmc].array.F[kc*AOconf[loop].sizeWFS+ii] = data.image[aoconfID_WFS1].array.F[ii];
		    }	
		  kc++;
		}
	    }
	}
      for(kk=0;kk<kc0max;kk++)
	{
	  printf("additional frame %ld [%ld/%ld]... ", kk, kc, RespMatNBframes+kc0max);
	  fflush(stdout);
	  Average_cam_frames(loop, 1);
	  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	    {
	      data.image[IDrefi].array.F[ii] += data.image[aoconfID_WFS1].array.F[ii];
	      data.image[IDrmc].array.F[kc*AOconf[loop].sizeWFS+ii] = data.image[aoconfID_WFS1].array.F[ii];	      
	    }
	  kc++;
	  printf("done\n");
	  fflush(stdout);
	}
      
      
      for(k2 = 0; k2 < AOconf[loop].NBDMmodes; k2++)
	data.image[aoconfID_cmd_modesRM].array.F[k2] = 0.0;
      set_DM_modesRM(loop);
      
      
      




      for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	for(k1=0;k1<AOconf[loop].NBDMmodes;k1++)
	  data.image[IDrmi].array.F[k1*AOconf[loop].sizeWFS+ii] /= (NBloops*2.0*amp*NbAve);
      
      for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	data.image[IDrefi].array.F[ii] /= RespMatNBframes+kc0max; //(NBloops*2.0*AOconf[loop].NBDMmodes*NbAve);
      
      printf("Acquisition done, compiling results...");
      fflush(stdout);
      
   
      
      
      // PROCESS RMCUBE  
      fp = fopen("TimeDelayRM.txt", "w");
      RMsig = 0.0;
      vOK = 1;
      IDrmtest = create_3Dimage_ID("rmtest", AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS, AOconf[loop].NBDMmodes);
      for(kc0=kc0min;kc0<kc0max;kc0++)
	{
	  // initialize RM to zero
	  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	    for(k=0;k<AOconf[loop].NBDMmodes;k++)
	      data.image[IDrmtest].array.F[k*AOconf[loop].sizeWFS+ii] = 0.0;
	  
	  
	  // initialize reference to zero  
	  kc = kc0;
	  for (kloop = 0; kloop < NBloops; kloop++)
	    {
	      for(k1 = 0; k1 < AOconf[loop].NBDMmodes; k1++)
		{	 
		  // positive
		  for(kk=0;kk<NbAve-NBexcl;kk++)
		    {		  
		      for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
			data.image[IDrmtest].array.F[k1*AOconf[loop].sizeWFS+ii] += data.image[IDrmc].array.F[kc*AOconf[loop].sizeWFS+ii];
		      kc++;
		    }
		  kc+=NBexcl;
		  
		  // negative	  
		  for(kk=0;kk<NbAve-NBexcl;kk++)
		    {
		      for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
			data.image[IDrmtest].array.F[k1*AOconf[loop].sizeWFS+ii] -= data.image[IDrmc].array.F[kc*AOconf[loop].sizeWFS+ii];
		      kc++;
		    }
		  kc+=NBexcl;
		}
	    }
	  RMsigold = RMsig; 
	  RMsig = 0.0;
	  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	    for(k1=0;k1<AOconf[loop].NBDMmodes;k1++)
	      {
		data.image[IDrmtest].array.F[k1*AOconf[loop].sizeWFS+ii] /= (NBloops*2.0*amp*(NbAve-NBexcl));
		RMsig += data.image[IDrmtest].array.F[k1*AOconf[loop].sizeWFS+ii]*data.image[IDrmtest].array.F[k1*AOconf[loop].sizeWFS+ii];
	      }
	  
	  if(RMsig<RMsigold)
	    vOK = 0;
	  printf("Delay = %ld frame(s)   ->  RM signal = %lf   %d\n", kc0, (double) RMsig, vOK);
	  fprintf(fp, "%ld %lf\n", kc0, (double) RMsig);      
	  if(RMsig<RMsigold)
	    vOK = 0;

	  if(vOK==1) // ADOPT THIS MATRIX
	    {
	      for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
		for(k1=0;k1<AOconf[loop].NBDMmodes;k1++)
		  data.image[IDrmi].array.F[k1*AOconf[loop].sizeWFS+ii] += data.image[IDrmtest].array.F[k1*AOconf[loop].sizeWFS+ii];      
	    }
	}
      fclose(fp);




      printf("\n");
      fflush(stdout);
      delete_image_ID("rmtest");
  
      beta = (1.0-gain)*beta + gain;
      for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	{
	  data.image[IDrefcumul].array.F[ii] = (1.0-gain)*data.image[IDrefcumul].array.F[ii] + gain*data.image[IDrefi].array.F[ii];
	  data.image[IDrefWFS].array.F[ii] = data.image[IDrefcumul].array.F[ii]/beta;
	  
	  for(k1=0;k1<AOconf[loop].NBDMmodes;k1++)
	    {
	      data.image[IDrmcumul].array.F[k1*AOconf[loop].sizeWFS+ii] = (1.0-gain)*data.image[IDrmcumul].array.F[k1*AOconf[loop].sizeWFS+ii] + gain*data.image[IDrmi].array.F[k1*AOconf[loop].sizeWFS+ii];
	      data.image[IDrespM].array.F[k1*AOconf[loop].sizeWFS+ii] = data.image[IDrmcumul].array.F[k1*AOconf[loop].sizeWFS+ii]/beta;
	    }
	}

      save_fits("refwfsacq", "!refwfs.fits");
      save_fits("respmacq", "!respm.fits");
    }

  

  printf("Done\n");
  free(sizearray);

  return(0);
}







int ControlMatrixMultiply( float *cm_array, float *imarray, long m, long n, float *outvect)
{
  long i;

  cblas_sgemv (CblasRowMajor, CblasNoTrans, m, n, 1.0, cm_array, n, imarray, 1, 0.0, outvect, 1);


  return(0);
}




int AOcompute(long loop)
{
  float total = 0.0;
  long k, k1, k2;
  long ii;
  long i;
  long m, n;
  long index;
  //  long long wcnt;
  // long long wcntmax;
  long cnttest;
  double a;
 

  // get dark-subtracted image
  Average_cam_frames(loop, AOconf[loop].framesAve);
  
  AOconf[loop].status = 4;  // 4: REMOVING REF
  


  for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
    data.image[aoconfID_WFS2].array.F[ii] = data.image[aoconfID_WFS1].array.F[ii] - data.image[aoconfID_refWFS].array.F[ii];
  cnttest = data.image[aoconfID_cmd1_modes].md[0].cnt0;
  data.image[aoconfID_WFS2].md[0].cnt0 ++;


  //  save_fits(data.image[aoconfID_WFS].md[0].name, "!testim.fits");
  // sleep(5);

  AOconf[loop].status = 5; // MULTIPLYING BY CONTROL MATRIX -> MODE VALUES


  if(AOconf[loop].GPU == 0)
    {
      ControlMatrixMultiply( data.image[aoconfID_contrM].array.F, data.image[aoconfID_WFS2].array.F, AOconf[loop].NBDMmodes, AOconf[loop].sizeWFS, data.image[aoconfID_cmd1_modes].array.F);
      data.image[aoconfID_cmd1_modes].md[0].cnt0 ++;
    }
  else
    {
      a = 0.1;
      while(cnttest==data.image[aoconfID_cmd1_modes].md[0].cnt0)
	{
	  a = sqrt(a+0.1);
	  //  usleep(10);
	  // printf(".");
	}
      //      printf("\n");
      
    }
 
  if(0)
    {
      printf("GSL MATRIX MULT  :  ");
      for(k=0; k<AOconf[loop].NBDMmodes; k++)
	printf(" %6.3f ",data.image[aoconfID_cmd1_modes].array.F[k]);
      printf("\n");       
    }
  
  if(0)
    {       
      for(m=0; m<AOconf[loop].NBDMmodes; m++)
	{
	  data.image[aoconfID_cmd1_modes].array.F[m] = 0.0;
	  for(n=0; n<AOconf[loop].sizeWFS; n++)
	    {
	      index = m*AOconf[loop].sizeWFS+n;
	      data.image[aoconfID_cmd1_modes].array.F[m] += data.image[aoconfID_contrM].array.F[index]*data.image[aoconfID_WFS2].array.F[n];
	    }
	}	  

      printf("CONV MATRIX MULT :  ");
      for(k=0; k<AOconf[loop].NBDMmodes; k++)
	printf(" %6.3f ",data.image[aoconfID_cmd1_modes].array.F[k]);
      printf("\n");          
   
      printf("cMat  : ");
      for(i=0;i<5;i++)
	printf("%f ", data.image[aoconfID_contrM].array.F[i]);
      printf(" ... ");
      for(i=AOconf[loop].sizeWFS*AOconf[loop].NBDMmodes-5;i<AOconf[loop].sizeWFS*AOconf[loop].NBDMmodes;i++)
	printf("%f ", data.image[aoconfID_contrM].array.F[i]);
      printf("\n");
      
      printf("wfsVec: ");
      for(n=0;n<5;n++)
	printf("%f ", data.image[aoconfID_WFS2].array.F[n]);
      printf(" ... ");
      for(n=AOconf[loop].sizeWFS-5;n<AOconf[loop].sizeWFS;n++)
	printf("%f ", data.image[aoconfID_WFS2].array.F[n]);
      printf("\n");
    }

  AOconf[loop].RMSmodes = 0;
  for(k=0; k<AOconf[loop].NBDMmodes; k++)
    AOconf[loop].RMSmodes += data.image[aoconfID_cmd1_modes].array.F[k]*data.image[aoconfID_cmd1_modes].array.F[k];
    

  AOconf[loop].status = 6; //  MULTIPLYING BY GAINS
   

  for(k=0; k<AOconf[loop].NBDMmodes; k++)
    {
      data.image[aoconfID_RMS_modes].array.F[k] = 0.99*data.image[aoconfID_RMS_modes].array.F[k] + 0.01*data.image[aoconfID_cmd1_modes].array.F[k]*data.image[aoconfID_cmd1_modes].array.F[k];
      data.image[aoconfID_AVE_modes].array.F[k] = 0.99*data.image[aoconfID_AVE_modes].array.F[k] + 0.01*data.image[aoconfID_cmd1_modes].array.F[k];
      
      data.image[aoconfID_cmd_modes].array.F[k] -= AOconf[loop].gain * data.image[aoconfID_GAIN_modes].array.F[k] * data.image[aoconfID_cmd1_modes].array.F[k];
      
      if(data.image[aoconfID_cmd_modes].array.F[k] < -AOconf[loop].maxlimit * data.image[aoconfID_LIMIT_modes].array.F[k])
	data.image[aoconfID_cmd_modes].array.F[k] = -AOconf[loop].maxlimit * data.image[aoconfID_LIMIT_modes].array.F[k];
      
      if(data.image[aoconfID_cmd_modes].array.F[k] > AOconf[loop].maxlimit * data.image[aoconfID_LIMIT_modes].array.F[k])
	data.image[aoconfID_cmd_modes].array.F[k] = AOconf[loop].maxlimit * data.image[aoconfID_LIMIT_modes].array.F[k];

      data.image[aoconfID_cmd_modes].array.F[k] *= data.image[aoconfID_MULTF_modes].array.F[k];
      

      // update total gain
      data.image[aoconfID_GAIN_modes].array.F[k+AOconf[loop].NBDMmodes] = AOconf[loop].gain * data.image[aoconfID_GAIN_modes].array.F[k];
    }




      //data.image[aoconfID_cmd_modes].array.F[k] *= 1.0 - 0.05*(1.0*k/AOconf[loop].NBDMmodes);


  data.image[aoconfID_cmd_modes].md[0].cnt0 ++;

  /*
    for(k=0; k<smao[0].NBmode; k++)
      {
	smao[0].mon_mode_residual[k] += smao[0].cmmode_array1[k];
	smao[0].mon_mode_residual_cnt[k]++;
	
	smao[0].mon_mode_residualRMS[k] += smao[0].cmmode_array1[k]*smao[0].cmmode_array1[k];
	smao[0].mon_mode_residualRMS_cnt[k]++;
	
	smao[0].mon_mode_applied[k] += smao[0].command_modes[k];
	smao[0].mon_mode_applied_cnt[k]++;
      }
  */

    AOconf[loop].status = 7; // fill in covariance matrix
    /*    if(smao[0].covM_on == 1)
      {
	for(k1=0; k1<smao[0].NBmode; k1++)
	  for(k2=0; k2<k1+1; k2++)
	    smao[0].covM[k2*smao[0].NBmode+k1] += smao[0].cmmode_array1[k1] * smao[0].cmmode_array1[k2];
	smao[0].covMcnt++;
	}*/

  return(0);
}


int AOloopControl_run()
{
  char fname[200];
  long loop;
  int vOK;

  long ID;
  long j, m;
  struct tm *uttime;
  time_t t;
  struct timespec *thetime = (struct timespec *)malloc(sizeof(struct timespec));
  char logfname[1000];
  char command[1000];
  int r;
  int RT_priority = 90; //any number from 0-99
  struct sched_param schedpar;
 

  schedpar.sched_priority = RT_priority;
  r = seteuid(euid_called); //This goes up to maximum privileges
  sched_setscheduler(0, SCHED_FIFO, &schedpar); //other option is SCHED_RR, might be faster

  loop = LOOPNUMBER;
  
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  


  printf("SETTING UP...\n");
  sprintf(fname, "AOloop%ld.conf", LOOPNUMBER);
  AOloopControl_loadconfigure(LOOPNUMBER, fname, 1);
 

  vOK = 1;
  if(AOconf[loop].init_refWFS==0)
    {
      printf("ERROR: CANNOT RUN LOOP WITHOUT WFS REFERENCE\n");
      vOK = 0;
    }
  if(AOconf[loop].init_CM==0)
    {
      printf("ERROR: CANNOT RUN LOOP WITHOUT CONTROL MATRIX\n");
      vOK = 0;
    }
  
  if(vOK==1)
    {
      
      AOconf[loop].kill = 0;
      AOconf[loop].on = 0;
      printf("\n");
      while( AOconf[loop].kill == 0)
	{
	  printf(" WAITING                    \r");
	  fflush(stdout);
	  usleep(100);
	  
	  while(AOconf[loop].on == 1)
	    {
	      //	      printf("LOOP IS RUNNING  %llu  %g      Gain = %f \r", AOconf[loop].cnt, AOconf[loop].RMSmodes, AOconf[loop].gain);
	      //	      fflush(stdout);
	      //usleep(10000);
	      
	      AOcompute(loop);
	    
	        
	      if(fabs(AOconf[loop].gain)>1.0e-6)
		set_DM_modes(loop);

	      clock_gettime(CLOCK_REALTIME, &AOconf[loop].tnow);
	      AOconf[loop].time_sec = 1.0*((long) AOconf[loop].tnow.tv_sec) + 1.0e-9*AOconf[loop].tnow.tv_nsec;

	      if(AOconf[loop].logfnb==0)
		ID = aoconfIDlog0;
	      else
		ID = aoconfIDlog1;

	      if(AOconf[loop].logcnt==0)
		{
		  AOconf[loop].timeorigin_sec = (long) AOconf[loop].tnow.tv_sec;
		  data.image[ID].kw[0].value.numl = AOconf[loop].timeorigin_sec;
		}
	 	   	      
	      data.image[ID].array.F[AOconf[loop].logcnt*data.image[ID].md[0].size[0]+0] = AOconf[loop].time_sec - 1.0*AOconf[loop].timeorigin_sec;
	      j = 1;
	     
	      for(m=0;m<AOconf[loop].NBDMmodes;m++)
		{
		  data.image[ID].array.F[AOconf[loop].logcnt*data.image[ID].md[0].size[0]+j] = AOconf[loop].gain;
		  j++;
		  data.image[ID].array.F[AOconf[loop].logcnt*data.image[ID].md[0].size[0]+j] = data.image[aoconfID_cmd1_modes].array.F[m];
		  j++;
		  data.image[ID].array.F[AOconf[loop].logcnt*data.image[ID].md[0].size[0]+j] = data.image[aoconfID_cmd_modes].array.F[m];
		  j++;		  
		}
	  	      

	      AOconf[loop].logcnt++;
	      if(AOconf[loop].logcnt==AOconf[loop].logsize)
		{
		  if(AOconf[loop].logon == 1) 
		    {
		      printf("Saving to disk...\n");
		      fflush(stdout);
		      		    
		      t = time(NULL);
		      uttime = gmtime(&t);
		      clock_gettime(CLOCK_REALTIME, thetime);
		      printf("writing file name\n");
		      fflush(stdout);
		      sprintf(logfname, "%s/LOOP%ld_%04d%02d%02d-%02d:%02d:%02d.%09ld%s.log", AOconf[loop].logdir, LOOPNUMBER, 1900+uttime->tm_year, 1+uttime->tm_mon, uttime->tm_mday, uttime->tm_hour, uttime->tm_min, uttime->tm_sec, thetime->tv_nsec, AOconf[loop].userLOGstring);		      
		      printf("writing file name\n");
		      fflush(stdout);
		      sprintf(command, "cp /tmp/loop%ldlog%d.im.shm %s &", loop, AOconf[loop].logfnb, logfname);
		      printf("Executing command : %s\n", command);
		      fflush(stdout);
		      r = system(command);
		    }
		  AOconf[loop].logfnb++;
		  AOconf[loop].logcnt = 0;
		}
	      if(AOconf[loop].logfnb == 2)
		  AOconf[loop].logfnb = 0;
	      
	      AOconf[loop].cnt++;
	    
	      if(AOconf[loop].cnt == AOconf[loop].cntmax)
		AOconf[loop].on = 0;
	    }
	  
	}
    }

  free(thetime);

  r = seteuid(euid_real);//Go back to normal privileges

  return(0);
}






int AOloopControl_printloopstatus(long loop, long nbcol)
{
  long k, kmax;
  long col;
  float val;
  
  float AVElim = 0.01;
  float RMSlim = 0.01;
  
  printw("loop number %ld    ", loop);

  
  if(AOconf[loop].on == 1)
    printw("loop is ON     ");
  else
    printw("loop is OFF    ");
  if(AOconf[loop].logon == 1)
    printw("log is ON   ");
  else
    printw("log is OFF  ");
  
  printw("STATUS = %d  ", AOconf[loop].status);
  
  kmax = (wrow-3)*(nbcol-1);
  printw("Gain = %f   maxlim = %f     GPU = %d    kmax=%ld\n", AOconf[loop].gain, AOconf[loop].maxlimit, AOconf[loop].GPU, kmax);
  printw("CNT : %lld  / %lld\n", AOconf[loop].cnt, AOconf[loop].cntmax);
  
  
  if(kmax>AOconf[loop].NBDMmodes)
    kmax = AOconf[loop].NBDMmodes;

  col = 0;
  for(k=0;k<kmax;k++)
    {
      attron(A_BOLD);
      printw("%4ld ", k);
      attroff(A_BOLD);
   
      printw("[%4.2f %4.2f %4.2f] ", data.image[aoconfID_GAIN_modes].array.F[k], data.image[aoconfID_LIMIT_modes].array.F[k], data.image[aoconfID_MULTF_modes].array.F[k]);
      
      val = data.image[aoconfID_cmd_modes].array.F[k];
      if(fabs(val)>0.99*AOconf[loop].maxlimit)
	{
	  attron(A_BOLD | COLOR_PAIR(2));
	  printw("%7.4f ", val);
	  attroff(A_BOLD | COLOR_PAIR(2));
	}
      else	
	{
	  if(fabs(val)>0.99*AOconf[loop].maxlimit*data.image[aoconfID_LIMIT_modes].array.F[k])
	    {
	      attron(COLOR_PAIR(1));
	      printw("%7.4f ", val);
	      attroff(COLOR_PAIR(1));
	    }
	  else
	    printw("%7.4f ", val);
	}

      printw("%7.4f ", data.image[aoconfID_cmd1_modes].array.F[k]);
      

      val = data.image[aoconfID_AVE_modes].array.F[k];
      if(fabs(val)>AVElim)
	{
	  attron(A_BOLD | COLOR_PAIR(2));
	  printw("%7.4f ", val);
	  attroff(A_BOLD | COLOR_PAIR(2));
	}
      else
	printw("%7.4f ", val);

        val = data.image[aoconfID_RMS_modes].array.F[k];
      if(fabs(val)>RMSlim)
	{
	  attron(A_BOLD | COLOR_PAIR(2));
	  printw("%7.4f ", val);
	  attroff(A_BOLD | COLOR_PAIR(2));
	}
      else
	printw("%7.4f ", val);

      col++;
      if(col==nbcol-1)
	{
	  col = 0;
	  printw("\n");
	}
      else
	printw(" | ");
    }


  return(0);
}



int AOloopControl_loopMonitor(long loop, double frequ, long nbcol)
{
  char name[200];

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  printf("MEMORY HAS BEEN INITIALIZED\n");
  fflush(stdout);

  // load arrays that are required
  if(aoconfID_cmd_modes==-1)
    {
      sprintf(name, "DMmode_cmd_%ld", loop);
      aoconfID_cmd_modes = read_sharedmem_image(name);
    }
     
  if(aoconfID_cmd1_modes==-1)
    {
      sprintf(name, "DMmode_cmd1_%ld", loop);
      aoconfID_cmd1_modes = read_sharedmem_image(name);
    }


   if(aoconfID_RMS_modes==-1)
     {
       sprintf(name, "DMmode_RMS_%ld", loop);
       aoconfID_RMS_modes = read_sharedmem_image(name);
     }
   
   if(aoconfID_AVE_modes==-1)
    {
      sprintf(name, "DMmode_AVE_%ld", loop);
      aoconfID_AVE_modes = read_sharedmem_image(name);
    }
   
   if(aoconfID_GAIN_modes==-1)
     {
       sprintf(name, "DMmode_GAIN_%ld", loop);
       aoconfID_GAIN_modes = read_sharedmem_image(name);
    }
   
   if(aoconfID_LIMIT_modes==-1)
     {
       sprintf(name, "DMmode_LIMIT_%ld", loop);
       aoconfID_LIMIT_modes = read_sharedmem_image(name);
     }

   if(aoconfID_MULTF_modes==-1)
     {
       sprintf(name, "DMmode_MULTF_%ld", loop);
       aoconfID_MULTF_modes = read_sharedmem_image(name);
     }

  
   initscr();		
   getmaxyx(stdscr, wrow, wcol);
   

   start_color();
   init_pair(1, COLOR_BLUE, COLOR_BLACK); 
   init_pair(2, COLOR_RED, COLOR_BLACK);
   init_pair(3, COLOR_GREEN, COLOR_BLACK);
   init_pair(4, COLOR_RED, COLOR_BLACK);
   
   while( !kbdhit() )
     {
       usleep((long) (1000000.0/frequ));
       clear();
       attron(A_BOLD);
       print_header(" PRESS ANY KEY TO STOP MONITOR ", '-');
       attroff(A_BOLD);
       
       AOloopControl_printloopstatus(loop, nbcol);
       
       refresh();
     }
   endwin();	
   
   return 0;
}







int AOloopControl_showparams(long loop)
{
  printf("loop number %ld\n", loop);
  if(AOconf[loop].on == 1)
    printf("loop is ON\n");
  else
    printf("loop is OFF\n");
  if(AOconf[loop].logon == 1)
    printf("log is ON\n");
  else
    printf("log is OFF\n");
  printf("Gain = %f   maxlim = %f\n  GPU = %d\n", AOconf[loop].gain, AOconf[loop].maxlimit, AOconf[loop].GPU);

  return 0;
}



int AOloopControl_setLoopNumber(long loop)
{
  printf("LOOPNUMBER = %ld\n", loop);
  LOOPNUMBER = loop;
  //  AOloopControl_showparams(loop);

  return 0;
}

int AOloopControl_loopkill()
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].kill = 1;

  return 0;
}

int AOloopControl_loopon()
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].cntmax = AOconf[LOOPNUMBER].cnt-1;

  AOconf[LOOPNUMBER].on = 1;
  AOloopControl_showparams(LOOPNUMBER);

  return 0;
}

int AOloopControl_loopstep(long loop, long NBstep)
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[loop].cntmax = AOconf[loop].cnt + NBstep;

  printf("\nLOOP %ld STEP    %lld %ld %lld\n\n", loop, AOconf[loop].cnt, NBstep, AOconf[loop].cntmax);
  fflush(stdout);
  
  AOconf[loop].on = 1;
  // AOloopControl_showparams(loop);

  return 0;
}



int AOloopControl_loopoff()
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].on = 0;
  AOloopControl_showparams(LOOPNUMBER);

  return 0;
}



int AOloopControl_loopreset()
{
  char name[200];
  long k;

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  if(aoconfID_cmd_modes==-1)
    {
      sprintf(name, "DMmode_cmd_%ld", LOOPNUMBER);
      aoconfID_cmd_modes = read_sharedmem_image(name);
    }


  AOconf[LOOPNUMBER].on = 0;
  for(k=0; k<AOconf[LOOPNUMBER].NBDMmodes; k++)
    data.image[aoconfID_cmd_modes].array.F[k] = 0.0;

  return 0;
}






int AOloopControl_logon()
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].logon = 1;
  AOloopControl_showparams(LOOPNUMBER);

  return 0;
}

int AOloopControl_logoff()
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].logon = 0;
  AOloopControl_showparams(LOOPNUMBER);

  return 0;
}

int AOloopControl_setgain(float gain)
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].gain = gain;
  AOloopControl_showparams(LOOPNUMBER);

  return 0;
}

int AOloopControl_setmaxlimit(float maxlimit)
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].maxlimit = maxlimit;
  AOloopControl_showparams(LOOPNUMBER);

  return 0;
}


int AOloopControl_setframesAve(long nbframes)
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].framesAve = nbframes;
  AOloopControl_showparams(LOOPNUMBER);

  return 0;
}



int AOloopControl_setgainrange(long m0, long m1, float gainval)
{
  long k;
  long kmax;
  char name[200];

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  if(aoconfID_GAIN_modes==-1)
    {
      sprintf(name, "DMmode_GAIN_%ld", LOOPNUMBER);
      aoconfID_GAIN_modes = read_sharedmem_image(name);
    }

  kmax = m1+1;
  if(kmax>AOconf[LOOPNUMBER].NBDMmodes)
    kmax = AOconf[LOOPNUMBER].NBDMmodes-1;

  for(k=m0;k<kmax;k++)
    data.image[aoconfID_GAIN_modes].array.F[k] = gainval;

  return 0;
}



int AOloopControl_setlimitrange(long m0, long m1, float limval)
{
  long k;
  long kmax;
  char name[200];

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  if(aoconfID_LIMIT_modes==-1)
    {
      sprintf(name, "DMmode_LIMIT_%ld", LOOPNUMBER);
      aoconfID_LIMIT_modes = read_sharedmem_image(name);
    }

  kmax = m1+1;
  if(kmax>AOconf[LOOPNUMBER].NBDMmodes)
    kmax = AOconf[LOOPNUMBER].NBDMmodes-1;

  for(k=m0;k<kmax;k++)
    data.image[aoconfID_LIMIT_modes].array.F[k] = limval;

  return 0;
}


int AOloopControl_setmultfrange(long m0, long m1, float multfval)
{
  long k;
  long kmax;
  char name[200];

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  if(aoconfID_MULTF_modes==-1)
    {
      sprintf(name, "DMmode_MULTF_%ld", LOOPNUMBER);
      aoconfID_MULTF_modes = read_sharedmem_image(name);
    }

  kmax = m1+1;
  if(kmax>AOconf[LOOPNUMBER].NBDMmodes)
    kmax = AOconf[LOOPNUMBER].NBDMmodes-1;

  for(k=m0;k<kmax;k++)
    data.image[aoconfID_MULTF_modes].array.F[k] = multfval;

  return 0;
}


int AOloopControl_setgainblock(long mb, float gainval)
{
  long k;
  char name[200];
  long kmin, kmax;

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  if(aoconfID_GAIN_modes==-1)
    {
      sprintf(name, "DMmode_GAIN_%ld", LOOPNUMBER);
      aoconfID_GAIN_modes = read_sharedmem_image(name);
    }

  if(mb<AOconf[LOOPNUMBER].NBMblocks)
    {
      if(mb==0)
	kmin = 0;
      else
	kmin = AOconf[LOOPNUMBER].indexmaxMB[mb-1];
      kmax = AOconf[LOOPNUMBER].indexmaxMB[mb];

      for(k=kmin; k<kmax; k++)
	data.image[aoconfID_GAIN_modes].array.F[k] = gainval;
    }

  return 0;
}


int AOloopControl_setlimitblock(long mb, float limitval)
{
  long k;
  char name[200];
  long kmin, kmax;

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  if(aoconfID_LIMIT_modes==-1)
    {
      sprintf(name, "DMmode_LIMIT_%ld", LOOPNUMBER);
      aoconfID_LIMIT_modes = read_sharedmem_image(name);
    }

  if(mb<AOconf[LOOPNUMBER].NBMblocks)
    {
      if(mb==0)
	kmin = 0;
      else
	kmin = AOconf[LOOPNUMBER].indexmaxMB[mb-1];
      kmax = AOconf[LOOPNUMBER].indexmaxMB[mb];

      for(k=kmin; k<kmax; k++)
	data.image[aoconfID_LIMIT_modes].array.F[k] = limitval;
    }

  return 0;
}


int AOloopControl_setmultfblock(long mb, float multfval)
{
  long k;
  char name[200];
  long kmin, kmax;

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  if(aoconfID_MULTF_modes==-1)
    {
      sprintf(name, "DMmode_MULTF_%ld", LOOPNUMBER);
      aoconfID_MULTF_modes = read_sharedmem_image(name);
    }

  if(mb<AOconf[LOOPNUMBER].NBMblocks)
    {
      if(mb==0)
	kmin = 0;
      else
	kmin = AOconf[LOOPNUMBER].indexmaxMB[mb-1];
      kmax = AOconf[LOOPNUMBER].indexmaxMB[mb];

      for(k=kmin; k<kmax; k++)
	data.image[aoconfID_MULTF_modes].array.F[k] = multfval;
    }

  return 0;
}

