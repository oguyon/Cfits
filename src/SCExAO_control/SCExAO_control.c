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

long pXsize = 120;
long pYsize = 120;

long SCExAO_DM_STAGE_Xpos = 0;
long SCExAO_DM_STAGE_Ypos = 0;
long SCExAO_Pcam_Xpos = 59000;
long SCExAO_Pcam_Ypos = 56000;

// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


int SCExAOcontrol_TakePyrWFS_image_cli()
{
	 if(CLI_checkarg(1,3)+CLI_checkarg(2,2)==0)
    {
      SCExAOcontrol_TakePyrWFS_image(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl);
      return 0;
    }
  else
    return 1;
}


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


  strcpy(data.cmd[data.NBcmd].key,"scexaotakepwfsim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = SCExAOcontrol_TakePyrWFS_image_cli;
  strcpy(data.cmd[data.NBcmd].info,"take pyr WFS camera image");
  strcpy(data.cmd[data.NBcmd].syntax,"<outname> <nbcoadd>");
  strcpy(data.cmd[data.NBcmd].example,"scexaotakepwfsim imp 100");
  strcpy(data.cmd[data.NBcmd].Ccall,"long SCExAOcontrol_TakePyrWFS_image(char *IDname, long NbAve)");
  data.NBcmd++;


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

  strcpy(data.cmd[data.NBcmd].key,"scexaopywfscamalign");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = SCExAOcontrol_PyramidWFS_AutoAlign_cam;
  strcpy(data.cmd[data.NBcmd].info,"move Camera to center pyrWFS");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"scexaopywfscamalign");
  strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_PyramidWFS_AutoAlign_cam();");
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
    long NBcoadd;
    long kw;
	double darkv;
	long IDv;
	
	
    IDWFScam = image_ID(WFScam_name);
    if(IDWFScam ==-1)
        IDWFScam = read_sharedmem_image(WFScam_name);

    kw = image_read_keyword_L(WFScam_name, "NBcoadd", &NBcoadd);
	
    if(kw==-1)
        printf("keyword not found\n");
    else
        printf("found [%ld] NBcoadd = %ld\n", kw, NBcoadd);

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
        data.image[ID].array.F[ii] /= NbAve*NBcoadd;

	if((IDv=variable_ID("AOLCAMDARK"))!=-1)
		{
			
			darkv = data.variable[IDv].value.f;
			printf("REMOVING DARK VALUE %f\n", darkv);
			
			for(ii=0; ii<xysize; ii++)
				data.image[ID].array.F[ii] -= darkv;
		}
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
	long ttxpos, ttypos;
	double gain = 1.0;


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

	/// 100 steps -> sig = 0.055 for modulation = 1.0
	ttxpos = (long) (SCExAO_DM_STAGE_Xpos - gain*100.0*(xsig/0.055));
	ttypos = (long) (SCExAO_DM_STAGE_Ypos - gain*100.0*(ysig/0.055));
	
	SCExAOcontrol_mv_DMstage(ttxpos, ttypos);

    save_fits("imwfs", "!imwfs.fits");

    return(0);
}

/** assumes imref has been loaded */
int SCExAOcontrol_PyramidWFS_AutoAlign_cam()
{
    FILE *fp;
    long ID, IDc;
    long IDref;
    long ii, jj;
    long brad = 30; // box radius
    double totx, toty, tot;
    double alpha = 20.0;
    double peak, v;
	double gain = 1.0;
	long stepx, stepy;
	int r;
	char command[200];
	long delayus = 1000000;
	long NBframes = 10000;
	
	/// read position of stages
	if((fp = fopen("pcampos.txt", "r"))!=NULL)
	{
		r = fscanf(fp, "%ld %ld\n", &SCExAO_Pcam_Xpos, &SCExAO_Pcam_Ypos);
		fclose(fp);
	}
	
    IDref = image_ID("imref");
    ID = SCExAOcontrol_TakePyrWFS_image("imwfs", NBframes);
	save_fits("imwfs", "!imwfs.fits");
	
    tot = 0.0;
    for(ii=0; ii<pXsize*pYsize; ii++)
        tot += data.image[ID].array.F[ii];
    for(ii=0; ii<pXsize*pYsize; ii++)
        data.image[ID].array.F[ii] /= tot;

    /** compute offset */
    fft_correlation("imwfs", "imref", "outcorr");
    IDc = image_ID("outcorr");
	peak = 0.0;
	for(ii=0; ii<pXsize*pYsize; ii++)
		if(data.image[IDc].array.F[ii]>peak)
			peak = data.image[IDc].array.F[ii];
    
    for(ii=0; ii<pXsize*pYsize; ii++)
        if(data.image[IDc].array.F[ii]>0.0)
            data.image[IDc].array.F[ii] = pow(data.image[IDc].array.F[ii]/peak,alpha);
        else
            data.image[IDc].array.F[ii] = 0.0;

    totx = 0.0;
    toty = 0.0;
    tot = 0.0;
    for(ii=pXsize/2-brad; ii<pXsize/2+brad; ii++)
        for(jj=pXsize/2-brad; jj<pXsize/2+brad; jj++)
        {
            v = data.image[IDc].array.F[jj*pXsize+ii];
            totx += 1.0*(ii-pXsize/2)*v;
            toty += 1.0*(jj-pXsize/2)*v;
            tot += v;
        }
    totx /= tot;
    toty /= tot;

    save_fits("outcorr", "!outcorr.fits");
    delete_image_ID("outcorr");

    printf("  %6.4f  x  %6.4f\n", totx, toty);

	stepx = (long) (gain*toty/3.0*400.0);
	stepy = (long) (gain*totx/3.0*400.0);

	printf("STEP : %ld %ld\n", stepx, stepy);

	SCExAO_Pcam_Xpos += stepx;
	SCExAO_Pcam_Ypos += stepy;

	/// write stages position
	fp = fopen("pcampos.txt", "w");
	fprintf(fp, "%ld %ld\n", SCExAO_Pcam_Xpos, SCExAO_Pcam_Ypos);
	fclose(fp);

	sprintf(command, "pywfs cam x goto %ld\n", SCExAO_Pcam_Xpos);
	r = system(command);
	usleep(delayus);
	
	sprintf(command, "pywfs cam y goto %ld\n", SCExAO_Pcam_Ypos);
	r = system(command);
	usleep(delayus);


    return(0);
}


