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

#include "SCExAO_control/SCExAO_control.h"

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



/// CONFIGURATION
char *WFScam_name = "zyladata";
long long WFScnt = 0;



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

   strcpy(data.cmd[data.NBcmd].key,"scexaopywfsttalign");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = SCExAOcontrol_PyramidWFS_AutoAlign_TT;
  strcpy(data.cmd[data.NBcmd].info,"move DM TT stage to center pyrWFS");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"scexaopywfsttalign");
  strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_PyramidWFS_AutoAlign_TT();");
  data.NBcmd++;


  // add atexit functions here
  
  
  return 0;
}




long SCExAOcontrol_TakePyrWFS_image(char *IDname, long NbAve)
{
    long ID;
    long IDWFScam;
    int slice;
    long k;
    long xsize, ysize, xysize;
    char *ptrv;
    unsigned short *arrayutmp;
    long ii;


    IDWFScam = image_ID(WFScam_name);
    xsize = data.image[IDWFScam].md[0].size[0];
    ysize = data.image[IDWFScam].md[0].size[1];
    xysize = xsize*ysize;

    ID = create_2Dimage_ID(IDname, xsize, ysize);

    arrayutmp = (unsigned short*) malloc(sizeof(unsigned short)*xysize);


    for(k=0; k<NbAve; k++)
    {
        while(WFScnt==data.image[IDWFScam].md[0].cnt0) // test if new frame exists
        {
            usleep(50);
            // do nothing, wait
        }

        slice = data.image[IDWFScam].md[0].cnt1-1;
        if(slice==-1)
            slice = data.image[IDWFScam].md[0].size[2]-1;

        ptrv = (char*) data.image[IDWFScam].array.U;
        ptrv += sizeof(unsigned short)*slice* xysize;
        memcpy (arrayutmp, ptrv, sizeof(unsigned short)*xysize);
        for(ii=0; ii<xysize; ii++)
            data.image[ID].array.F[ii] += (float) arrayutmp[ii];


        WFScnt = data.image[IDWFScam].md[0].cnt0;
    }

    for(ii=0; ii<xysize; ii++)
        data.image[ID].array.F[ii] /= NbAve;

    free(arrayutmp);

    return(ID);
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

   printf("X: %ld -> %ld      Y: %ld -> %ld\n", SCExAO_DM_STAGE_Xpos, stepXpos, SCExAO_DM_STAGE_Ypos, stepYpos);
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









/** auto aligns tip-tilt by equalizing fluxes between quadrants */

int SCExAOcontrol_PyramidWFS_AutoAlign_TT()
{
    long ID;
    long xsize, ysize;
    long ii, jj;
    double tot00, tot01, tot10, tot11, tot;
	double xsig, ysig;


    ID = SCExAOcontrol_TakePyrWFS_image("imwfs", 10);
    xsize = data.image[ID].md[0].size[0];
    ysize = data.image[ID].md[0].size[1];

	printf("%ld x %ld image\n", xsize, ysize);

    tot00 = 0.0;
    tot01 = 0.0;
    tot10 = 0.0;
    tot11 = 0.0;

    for(ii=0; ii<xsize/2; ii++)
        for(jj=0; jj<ysize/2; jj++)
            tot00 += data.image[ID].array.F[jj*xsize+ii];

    for(ii=xsize/2; ii<xsize; ii++)
        for(jj=0; jj<ysize/2; jj++)
            tot10 += data.image[ID].array.F[jj*xsize+ii];

    for(ii=0; ii<xsize/2; ii++)
        for(jj=ysize/2; jj<ysize; jj++)
            tot01 += data.image[ID].array.F[jj*xsize+ii];

    for(ii=xsize/2; ii<xsize; ii++)
        for(jj=ysize/2; jj<ysize; jj++)
            tot11 += data.image[ID].array.F[jj*xsize+ii];

    tot = tot00+tot10+tot01+tot11;
    tot00 /= tot;
    tot10 /= tot;
    tot01 /= tot;
    tot11 /= tot;

    printf("  %6.4f   %6.4f\n", tot01, tot11);
    printf("  %6.4f   %6.4f\n", tot00, tot10);

	xsig = tot01-tot10;
	ysig = tot11-tot00;
	printf(" sig = %6.4f  x %6.4f\n", xsig, ysig);

    save_fits("imwfs", "!imwfs.fits");

    return(0);
}

