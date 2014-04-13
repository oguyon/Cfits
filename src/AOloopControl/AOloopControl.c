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

#include "Cfits.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "AOloopControl/AOloopControl.h"

extern DATA data;


long NB_AOcontrol = 2; // max number of loops
long LOOPNUMBER; // current loop index



AOLOOPCONTROL_CONF *AOconf; // configuration - this can be an array
int *AOconf_loaded = 0;
int *AOconf_fd; 



int AOloopControl_run();

// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//








int init_AOloopControl()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "AO loop control");
  data.NBmodule++;


  strcpy(data.cmd[data.NBcmd].key,"AOlooprun");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_run;
  strcpy(data.cmd[data.NBcmd].info,"run AO loop");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"AOlooprun");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_run()");
  data.NBcmd++;





  // add atexit functions here
  // atexit((void*) SCEXAO_DM_unloadconf);
  
  return 0;
}










int AOloopControl_run()
{
  printf("running loop\n");
  
  return(0);
}
