#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "Cfits.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "AOsystSim/AOsystSim.h"

extern DATA data;


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


int AOsystSim_run_cli()
{
  
  AOsystSim_run();

  return 0;
}




int init_AOsystSim()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "conversion between image format, I/O");
  data.NBmodule++;


  strcpy(data.cmd[data.NBcmd].key,"AOsystsim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOsystSim_run_cli;
  strcpy(data.cmd[data.NBcmd].info,"run fake AO system");
  strcpy(data.cmd[data.NBcmd].syntax,"no argument");
  strcpy(data.cmd[data.NBcmd].example,"AOsystsim");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOsystSim_run()");
  data.NBcmd++;


 // add atexit functions here

  return 0;
}




int AOsystSim_run()
{
  long dmxsize = 50;
  long dmysize = 50;
  long *dmsize;
  long twait = 1; // us
  long long cnt = 0;
  struct timespec sleep_time = { 0, 1000000};
  long ID_dm;
  double a;

  printf("Running fake AO system simulation\n");

  dmsize = (long*) malloc(sizeof(long)*2);
  
  // INITIALIZE
  
  // create DM
  dmsize[0] = dmxsize;
  dmsize[1] = dmysize;
  ID_dm = create_image_ID("dm_sim", 2, dmsize, FLOAT, 1, 0);

  a = 0.0;
  while(cnt<1000)
    {
      a += cos(cnt);
      cnt ++;
      clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &sleep_time, NULL);
      //nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, time, NULL);
      //usleep(1);
    }
 
  printf("%ld %lf\n", (long) cnt, a);
  free(dmsize);

  return 0;
}
