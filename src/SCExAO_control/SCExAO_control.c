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

#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "linopt_imtools/linopt_imtools.h"
#include "image_filter/image_filter.h"

#include "ZernikePolyn/ZernikePolyn.h"



# ifdef _OPENMP
# include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
# endif


/** SCExAO instrument control
 * 
 * These routines only work on SCExAO computer, and are specific to SCExAO system
 * 
 */



extern DATA data;




long SCExAO_DM_STAGE_Xpos = 0;
long SCExAO_DM_STAGE_Ypos = 0;


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//

int SCExAOcontrol_mv_DMstage_cli()
{
	 if(CLI_checkarg(1,2)+CLI_checkarg(2,2)==0)
    {
      SCExAOcontrol_mv_DMstage(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl);
      return 0;
    }
  else
    return 1;
}


int init_SCExAO_control()
{

  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "SCExAO control");
  data.NBmodule++;

  strcpy(data.cmd[data.NBcmd].key,"scexaottdmpos");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = SCExAOcontrol_mv_DMstage_cli;
  strcpy(data.cmd[data.NBcmd].info,"move DM TT stage to position");
  strcpy(data.cmd[data.NBcmd].syntax,"<x pos> <y pos>");
  strcpy(data.cmd[data.NBcmd].example,"scexaottdmpos");
  strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_mv_DMstage(long stepXpos, long stepYpos)");
  data.NBcmd++;

  

  // add atexit functions here
  
  
  return 0;
}






/** \brief Move DM stage 
 * 
 *  Absolute position
 * 
 */

int SCExAOcontrol_mv_DMstage(long stepXpos, long stepYpos)
{
    char command[200];
    long ABoffset = 200; /// anti-backlash offset - rule: go negative first, and then positive
    int r;
    long delayus = 2000000;
    long stepX, stepY;

    stepX = stepXpos - SCExAO_DM_STAGE_Xpos;
    stepY = stepYpos - SCExAO_DM_STAGE_Ypos;

	printf("Moving by  %ld x %ld\n", stepX, stepY);

    if((fabs(stepX)>500.0)||(fabs(stepY)>500)||(fabs(stepXpos)>1000.0)||(fabs(stepYpos)>1000))
    {
        printf("ERROR: motion is too large or out of range: ignoring command\n");
    }
    else
    {

        if(stepX!=0)
        {
            if(stepX>ABoffset)
            {
                sprintf(command, "dm_stage x push %ld\n", stepX);
                printf("command : %s\n", command);
                r = system(command);
                usleep(delayus);
            }
            else
            {
                sprintf(command, "dm_stage x push %ld\n", stepX-ABoffset);
                printf("command : %s\n", command);
                r = system(command);
                usleep(delayus);

                sprintf(command, "dm_stage x push %ld\n", ABoffset);
                printf("command : %s\n", command);
                r = system(command);
                usleep(delayus);
            }

            SCExAO_DM_STAGE_Xpos += stepX;
        }

        if(stepY!=0)
        {
            if(stepY>ABoffset)
            {
                sprintf(command, "dm_stage y push %ld\n", stepY);
                printf("command : %s\n", command);
                r = system(command);
                usleep(delayus);
            }
            else
            {
                sprintf(command, "dm_stage y push %ld\n", stepY-ABoffset);
                printf("command : %s\n", command);
               r = system(command);
                usleep(delayus);

                sprintf(command, "dm_stage y push %ld\n", ABoffset);
                printf("command : %s\n", command);
                r = system(command);
                usleep(delayus);
            }

            SCExAO_DM_STAGE_Ypos += stepY;
        }
    }

    return(0);
}










int SCExAOcontrol_PyramidWFS_AutoAlign_TT(char *IDref_name, char *IDWFScam_name)
{
	long IDref, IDWFScam;
	
	
	return(0);
}
