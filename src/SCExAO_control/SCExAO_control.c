#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include <math.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <err.h>
#include <fcntl.h>
#include <sched.h>
#include <ncurses.h>
#include <semaphore.h>

#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "linopt_imtools/linopt_imtools.h"
#include "image_filter/image_filter.h"
#include "info/info.h"
#include "fft/fft.h"
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

float PcamPixScaleAct = 0.7*10000.0; // pyramid re-image actuators: step per pixel

/// CONFIGURATION
/// CONFIGURATION
char WFScam_name[200];
long long WFScnt = 0;

long pXsize = 120;
long pYsize = 120;

long SCExAO_DM_STAGE_Xpos = 0;
long SCExAO_DM_STAGE_Ypos = 0;

long SCExAO_Pcam_Xpos0 = 60000;
long SCExAO_Pcam_Ypos0 = 62000;
long SCExAO_Pcam_Xpos = 60000;
long SCExAO_Pcam_Ypos = 62000;
long SCExAO_Pcam_Range = 50000;

float SCExAO_PZT_STAGE_Xpos = -5.0;
float SCExAO_PZT_STAGE_Xpos_ref = -5.0;
float SCExAO_PZT_STAGE_Xpos_min = -7.0;
float SCExAO_PZT_STAGE_Xpos_max = -3.0;

float SCExAO_PZT_STAGE_Ypos = -5.0;
float SCExAO_PZT_STAGE_Ypos_ref = -5.0;
float SCExAO_PZT_STAGE_Ypos_min = -7.0;
float SCExAO_PZT_STAGE_Ypos_max = -3.0;

// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


int SCExAOcontrol_Average_image_cli()
{
    if(CLI_checkarg(2,2)+CLI_checkarg(3,3)+CLI_checkarg(4,2)==0)
    {

        SCExAOcontrol_Average_image(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.string, data.cmdargtoken[4].val.numl);
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



int SCExAOcontrol_PyramidWFS_AutoAlign_TT_cli()
{
    if(CLI_checkarg(1,4)==0)
    {
        SCExAOcontrol_PyramidWFS_AutoAlign_TT(data.cmdargtoken[1].val.string);
        return 0;
    }
    else
        return 1;
}


int SCExAOcontrol_PyramidWFS_AutoAlign_cam_cli()
{
    if(CLI_checkarg(1,4)==0)
    {
        SCExAOcontrol_PyramidWFS_AutoAlign_cam(data.cmdargtoken[1].val.string);
        return 0;
    }
    else
        return 1;
}


int SCExAOcontrol_PyramidWFS_Pcenter_cli()
{
    if(CLI_checkarg(1,4)+CLI_checkarg(2,1)+CLI_checkarg(3,1)==0)
    {
        SCExAOcontrol_PyramidWFS_Pcenter(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numf);
        return 0;
    }
    else
        return 1;
}


int SCExAOcontrol_Pyramid_flattenRefWF_cli()
{
    if(CLI_checkarg(1,4)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
    {
        SCExAOcontrol_Pyramid_flattenRefWF(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numf);
        return 0;
    }
    else
        return 1;
}


int SCExAOcontrol_optPSF_cli()
{
    if(CLI_checkarg(1,4)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
    {
        SCExAOcontrol_optPSF(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numf);
        return 0;
    }
    else
        return 1;
}





int SCExAOcontrol_SAPHIRA_cam_process_cli()
{
    if(CLI_checkarg(1,4)+CLI_checkarg(2,3)==0)
    {
        SCExAOcontrol_SAPHIRA_cam_process(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);
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


    strcpy(data.cmd[data.NBcmd].key,"scexaoaveim");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = SCExAOcontrol_Average_image_cli;
    strcpy(data.cmd[data.NBcmd].info,"take averaged camera image. Image in shared mem is <imname>.im.shm");
    strcpy(data.cmd[data.NBcmd].syntax,"<imname> <nbcoadd> <output image> <semaphore index>");
    strcpy(data.cmd[data.NBcmd].example,"scexaoaveim cam1 100 outave 3");
    strcpy(data.cmd[data.NBcmd].Ccall,"long SCExAOcontrol_Average_image(char *imname, long NbAve, char *IDnameout, long semindex)");
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
    data.cmd[data.NBcmd].fp = SCExAOcontrol_PyramidWFS_AutoAlign_TT_cli;
    strcpy(data.cmd[data.NBcmd].info,"move TT to center pyrWFS");
    strcpy(data.cmd[data.NBcmd].syntax,"<wfscamname>");
    strcpy(data.cmd[data.NBcmd].example,"scexaopywfsttalign wfscam");
    strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_PyramidWFS_AutoAlign_TT(char *WFScam_name);");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"scexaopywfscamalign");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = SCExAOcontrol_PyramidWFS_AutoAlign_cam_cli;
    strcpy(data.cmd[data.NBcmd].info,"move Camera to center pyrWFS");
    strcpy(data.cmd[data.NBcmd].syntax,"<wfscamname>");
    strcpy(data.cmd[data.NBcmd].example,"scexaopywfscamalign");
    strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_PyramidWFS_AutoAlign_cam();");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"scexaopypcent");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = SCExAOcontrol_PyramidWFS_Pcenter_cli;
    strcpy(data.cmd[data.NBcmd].info,"center pyrWFS pupil");
    strcpy(data.cmd[data.NBcmd].syntax,"<wfsimname> <pup radius [float]>");
    strcpy(data.cmd[data.NBcmd].example,"scexaopypcent wfsim 25.0");
    strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_PyramidWFS_Pcenter(char *IDwfsname, float prad)");
    data.NBcmd++;


    strcpy(data.cmd[data.NBcmd].key,"scexaopyflatten");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = SCExAOcontrol_Pyramid_flattenRefWF_cli;
    strcpy(data.cmd[data.NBcmd].info,"flatten  pyrWFS");
    strcpy(data.cmd[data.NBcmd].syntax,"<wfscamname> <NB zern> <ampl>");
    strcpy(data.cmd[data.NBcmd].example,"scexaopyflatten wfsim 20 0.05");
    strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_Pyramid_flattenRefWF(char *WFScam_name, long zimaxmax, float ampl0);");
    data.NBcmd++;


	strcpy(data.cmd[data.NBcmd].key,"scexaoPSFopt");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = SCExAOcontrol_optPSF_cli;
    strcpy(data.cmd[data.NBcmd].info,"optimize PSF shape");
    strcpy(data.cmd[data.NBcmd].syntax,"<wfscamname> <NB modes> <ampl>");
    strcpy(data.cmd[data.NBcmd].example,"scexaoPSFopt wfsim 20 0.05");
    strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_optPSF(char *WFScam_name, long zimaxmax, float alpha)");
    data.NBcmd++;
     

    strcpy(data.cmd[data.NBcmd].key,"scexaosaphiraproc");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = SCExAOcontrol_SAPHIRA_cam_process_cli;
    strcpy(data.cmd[data.NBcmd].info,"process saphira camera images");
    strcpy(data.cmd[data.NBcmd].syntax,"<input> <output>");
    strcpy(data.cmd[data.NBcmd].example,"scexaosaphiraproc");
    strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAOcontrol_SAPHIRA_cam_process(char *IDinname, char *IDoutname)");
    data.NBcmd++;




    // add atexit functions here


    return 0;
}





long SCExAOcontrol_Average_image(char *imname, long NbAve, char *IDnameout, long semindex)
{
    long ID;
    long IDdark;
    long IDcam;
    int slice;
    long k;
    long xsize, ysize, xysize;
    char *ptrv;
    unsigned short *arrayutmp;
    float *arraytmp;
    long ii;
    long kw;
    double darkv;
    long IDv;
    char imnameave[200];
    long long cntref;
    int semval;
    

    cntref = -1;

    IDcam = image_ID(imname);
    if(IDcam ==-1)
        IDcam = read_sharedmem_image(imname);


	

    xsize = data.image[IDcam].md[0].size[0];
    ysize = data.image[IDcam].md[0].size[1];
    xysize = xsize*ysize;

    ID = create_2Dimage_ID(IDnameout, xsize, ysize);

//    list_image_ID();

	if(data.image[IDcam].md[0].atype == FLOAT)
		arraytmp = (float*) malloc(sizeof(float)*xysize);
	else
		arrayutmp = (unsigned short*) malloc(sizeof(unsigned short)*xysize);


    for(k=0; k<NbAve; k++)
    {
        //printf("k = %ld\n", k);
        //fflush(stdout);
        
        if(data.image[IDcam].sem<semindex)
        {
            while(cntref==data.image[IDcam].md[0].cnt0) // test if new frame exists
                usleep(10);
        }
        else
        {
//            sem_getvalue(data.image[IDcam].semptr[4], &semval);
            sem_wait(data.image[IDcam].semptr[(int) semindex]);
        }

        //slice = data.image[IDcam].md[0].cnt1;
        //if(slice==-1)
         //   slice = data.image[IDcam].md[0].size[2]-1;

        
		if(data.image[IDcam].md[0].atype == FLOAT)
			{
				ptrv = (char*) data.image[IDcam].array.F;
				memcpy (arraytmp, ptrv, sizeof(float)*xysize);
				for(ii=0; ii<xysize; ii++)
					data.image[ID].array.F[ii] += arraytmp[ii];
			}
        else
			{
				ptrv = (char*) data.image[IDcam].array.U;
				memcpy (arrayutmp, ptrv, sizeof(unsigned short)*xysize);
				for(ii=0; ii<xysize; ii++)
					data.image[ID].array.F[ii] += (float) arrayutmp[ii];
			}
        
        
        //ptrv += sizeof(unsigned short)*slice*xysize;
    

        cntref = data.image[IDcam].md[0].cnt0;
    }

    for(ii=0; ii<xysize; ii++)
        data.image[ID].array.F[ii] /= NbAve;


    if((IDdark=image_ID("wfsdark"))!=-1)
    {

        for(ii=0; ii<xysize; ii++)
            data.image[ID].array.F[ii] -= data.image[IDdark].array.F[ii];
    }

	if(data.image[IDcam].md[0].atype == FLOAT)
		free(arraytmp);
    else
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

int SCExAOcontrol_PyramidWFS_AutoAlign_TT_DM(char *WFScam_name)
{
    long ID;
    long xsize, ysize;
    long ii, jj;
    double tot00, tot01, tot10, tot11, tot;
    double xsig, ysig;
    long ttxpos, ttypos;
    double gain = 1.0;


    ID = SCExAOcontrol_Average_image(WFScam_name, 5000, "imwfs", 7);
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
    printf("tot = %f\n", tot);
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



int SCExAOcontrol_PyramidWFS_AutoAlign_TT(char *WFScam_name)
{
    FILE *fp;
    long ID;
    long xsize, ysize;
    long ii, jj;
    double tot00, tot01, tot10, tot11, tot;
    double tot00x, tot01x, tot10x, tot11x;
    double tot00y, tot01y, tot10y, tot11y;
    double xsig, ysig;
    long ttxpos, ttypos;
    double gain = 1.0;
    char command[200];
    int r;
    double x, y;
    double totx, toty;
    char pausefilename[200];
    float v0;
    long IDshm;
    long *sizearray;

    long NBframesAve;
    long NBframesAveMin = 500;
    long NBframesAveMax = 30000;
    long twaitus = 500000; // 0.5 sec

    float gainfactor;

    //        SCExAOcontrol_PyramidWFS_AutoAlign_TT_DM();
    // exit(0);

    SCExAO_PZT_STAGE_Xpos = -5.0;
    SCExAO_PZT_STAGE_Ypos = -5.0;

    IDshm = image_ID("pyrTT");
    if(IDshm == -1)
    {
        sizearray = (long*) malloc(sizeof(long)*2);
        sizearray[0] = 2;
        sizearray[1] = 1;
        IDshm = create_image_ID("pyrTT", 2, sizearray, FLOAT, 1, 0);
    }


    NBframesAve = NBframesAveMin;
    gainfactor = 1.0;
    
    
    
    while(file_exists("stop_PyAlignTT.txt")==0)
    {

        while (file_exists("pause_PyAlignTT.txt"))
            usleep(100000);

        if(file_exists("./status/gain_PyAlignTT.txt"))
        {
            fp = fopen("./status/gain_PyAlignTT.txt", "r");
            r = fscanf(fp, "%f", &v0);
            fclose(fp);
            if((v0>0.0)&&(v0<1.0))
                gain = gainfactor*v0;
        }
        
        gainfactor = 0.99*gainfactor;
               
        gain *= gainfactor;
        if(gain < 0.1)
            gain = 0.1;
        
        printf("================== AVERAGING %6ld FRAMES    gain = %f ================ \n", NBframesAve, gain);
        ID = SCExAOcontrol_Average_image(WFScam_name, NBframesAve, "imwfs", 4);
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];
        
        NBframesAve = (long) (1.1*NBframesAve);
        if (NBframesAve>NBframesAveMax)
            NBframesAve = NBframesAveMax;
        
        
        // save_fits("imwfs", "!imwfs.fits"); // TEST

        printf("%ld x %ld image\n", xsize, ysize);

        tot00 = 0.0;
        tot01 = 0.0;
        tot10 = 0.0;
        tot11 = 0.0;

        tot00x = 0.0;
        tot01x = 0.0;
        tot10x = 0.0;
        tot11x = 0.0;

        tot00y = 0.0;
        tot01y = 0.0;
        tot10y = 0.0;
        tot11y = 0.0;

   


        for(ii=0; ii<xsize/2; ii++)
            for(jj=0; jj<ysize/2; jj++)
            {
                x = 1.0*(0.5+ii)/(xsize/2)-0.5;
                y = 1.0*(0.5+jj)/(ysize/2)-0.5;
                tot00x += x*data.image[ID].array.F[jj*xsize+ii];
                tot00y += y*data.image[ID].array.F[jj*xsize+ii];
                tot00 += data.image[ID].array.F[jj*xsize+ii];
            }

        for(ii=xsize/2; ii<xsize; ii++)
            for(jj=0; jj<ysize/2; jj++)
            {
                x = 1.0*(0.5+ii-xsize/2)/(xsize/2)-0.5;
                y = 1.0*(0.5+jj)/(ysize/2)-0.5;
                tot10x += x*data.image[ID].array.F[jj*xsize+ii];
                tot10y += y*data.image[ID].array.F[jj*xsize+ii];
                tot10 += data.image[ID].array.F[jj*xsize+ii];
            }

        for(ii=0; ii<xsize/2; ii++)
            for(jj=ysize/2; jj<ysize; jj++)
            {
                x = 1.0*(0.5+ii)/(xsize/2)-0.5;
                y = 1.0*(0.5+jj-ysize/2)/(ysize/2)-0.5;
                tot01x += x*data.image[ID].array.F[jj*xsize+ii];
                tot01y += y*data.image[ID].array.F[jj*xsize+ii];
                tot01 += data.image[ID].array.F[jj*xsize+ii];
            }

        for(ii=xsize/2; ii<xsize; ii++)
            for(jj=ysize/2; jj<ysize; jj++)
            {
                x = 1.0*(0.5+ii-xsize/2)/(xsize/2)-0.5;
                y = 1.0*(0.5+jj-ysize/2)/(ysize/2)-0.5;
                tot11x += x*data.image[ID].array.F[jj*xsize+ii];
                tot11y += y*data.image[ID].array.F[jj*xsize+ii];
                tot11 += data.image[ID].array.F[jj*xsize+ii];
            }

        tot = tot00+tot10+tot01+tot11;

        tot00x /= tot00;
        tot10x /= tot10;
        tot01x /= tot01;
        tot11x /= tot11;

        tot00y /= tot00;
        tot10y /= tot10;
        tot01y /= tot01;
        tot11y /= tot11;

        tot00 /= tot;
        tot10 /= tot;
        tot01 /= tot;
        tot11 /= tot;



        printf("  %6.4f   %6.4f\n", tot01, tot11);
        printf("  %6.4f   %6.4f\n", tot00, tot10);
        printf("tot = %f   ave = %f \n", tot, tot/xsize/ysize);

        totx = 0.25*(tot00x+tot10x+tot10x+tot11x);
        toty = 0.25*(tot00y+tot10y+tot10y+tot11y);

        printf(" PUP X   %+6.4f %+6.4f %+6.4f %+6.4f  -> %+6.4f\n", tot00x, tot01x, tot10x, tot11x, totx);
        printf(" PUP Y   %+6.4f %+6.4f %+6.4f %+6.4f  -> %+6.4f\n", tot00y, tot01y, tot10y, tot11y, toty);

        xsig = (tot10+tot11)-(tot00+tot01); // camera coordinates
        ysig = (tot01+tot11)-(tot00+tot10);
        printf(" sig = %6.4f  x %6.4f\n", xsig, ysig);

        //exit(0);

        /// 1 V step -> sig = 0.2 for modulation = 0.3
        //       SCExAO_PZT_STAGE_Xpos += gain*(xsig/0.2);
        //     SCExAO_PZT_STAGE_Ypos -= gain*(ysig/0.2);



        if(tot > 10.0*xsize*ysize)
        {
            SCExAO_PZT_STAGE_Xpos -= gain*((xsig-ysig)/1.0);  // C actuator
            SCExAO_PZT_STAGE_Ypos -= gain*((xsig+ysig)/1.0);  // D actuator



            printf("  --- %f  %f ----\n", SCExAO_PZT_STAGE_Xpos, SCExAO_PZT_STAGE_Ypos);


            if(SCExAO_PZT_STAGE_Xpos<SCExAO_PZT_STAGE_Xpos_min)
                SCExAO_PZT_STAGE_Xpos = SCExAO_PZT_STAGE_Xpos_min;
            if(SCExAO_PZT_STAGE_Xpos>SCExAO_PZT_STAGE_Xpos_max)
                SCExAO_PZT_STAGE_Xpos = SCExAO_PZT_STAGE_Xpos_max;

            if(SCExAO_PZT_STAGE_Ypos<SCExAO_PZT_STAGE_Ypos_min)
                SCExAO_PZT_STAGE_Ypos = SCExAO_PZT_STAGE_Ypos_min;
            if(SCExAO_PZT_STAGE_Ypos>SCExAO_PZT_STAGE_Ypos_max)
                SCExAO_PZT_STAGE_Ypos = SCExAO_PZT_STAGE_Ypos_max;

            // sig X
            //sprintf(command, "analog_output.py voltage C %5.3f\n", SCExAO_PZT_STAGE_Xpos);
            sprintf(command, "./aocustomscripts/SCExAO_analogoutput D %5.3f", SCExAO_PZT_STAGE_Xpos);
            printf("COMMAND: \"%s\"\n", command);
            r = system(command);

            // sig Y
            //sprintf(command, "analog_output.py voltage D %5.3f\n", SCExAO_PZT_STAGE_Ypos);
            sprintf(command, "./aocustomscripts/SCExAO_analogoutput C %5.3f", SCExAO_PZT_STAGE_Ypos);
            printf("COMMAND: \"%s\"\n", command);
            r = system(command);

            data.image[IDshm].md[0].write = 1;
            data.image[IDshm].array.F[0] = SCExAO_PZT_STAGE_Xpos;
            data.image[IDshm].array.F[1] = SCExAO_PZT_STAGE_Ypos;
            data.image[IDshm].md[0].cnt0 ++;
            data.image[IDshm].md[0].write = 0;
        }

        save_fits("imwfs", "!./tmp/imwfs_alignTT.fits");
        usleep(twaitus);
    }

    r = system("rm stop_PyAlignTT.txt");

    return(0);
}











/** assumes imref has been loaded */
int SCExAOcontrol_PyramidWFS_AutoAlign_cam(char *WFScam_name)
{
    FILE *fp;
    long ID, IDc;
    long ii, jj;
    long brad = 30; // box radius
    double totx, toty, tot;
    double alpha = 20.0;
    double peak, v;
    double gain = 0.2;
    long stepx, stepy;
    int r;
    char command[200];
    long delayus = 500000;
    
    long NBframes = 20000;
    
    float v0;
    long maxstep = 3000;
    float ave;
    char pausefilename[200];


    long NBframesAve;
    long NBframesAveMin = 500;
    long NBframesAveMax = 30000;
    float gainfactor;



    /// read position of stages
    if((fp = fopen("./status/pcampos.txt", "r"))!=NULL)
    {
        r = fscanf(fp, "%ld %ld\n", &SCExAO_Pcam_Xpos, &SCExAO_Pcam_Ypos);
        fclose(fp);
    }

    SCExAO_Pcam_Xpos0 = SCExAO_Pcam_Xpos;
    SCExAO_Pcam_Ypos0 = SCExAO_Pcam_Ypos;

    NBframesAve = NBframesAveMin;
    gainfactor = 1.0;
    while(file_exists ("stop_PyAlignCam.txt")==0)
    {
        while (file_exists ("pause_PyAlignCam.txt"))
            usleep(100000);



        if(file_exists("./status/gain_PyAlignCam.txt"))
        {
            fp = fopen("./status/gain_PyAlignCam.txt", "r");
            r = fscanf(fp, "%f", &v0);
            fclose(fp);
            if((v0>0.0)&&(v0<1.0))
                gain = gainfactor*v0;
        }
        gainfactor = 0.95*gainfactor;
        if(gainfactor < 0.1)
            gainfactor = 0.1;


        printf("================== AVERAGING %6ld FRAMES    gain = %f ================ \n", NBframesAve, gain);
        ID = SCExAOcontrol_Average_image(WFScam_name, NBframesAve, "imwfs", 5);
        save_fits("imwfs", "!./tmp/imwfs_aligncam.fits");
  
        NBframesAve = (long) (1.1*NBframesAve);
        if (NBframesAve>NBframesAveMax)
            NBframesAve = NBframesAveMax;
        

        tot = 0.0;
        for(ii=0; ii<pXsize*pYsize; ii++)
            tot += data.image[ID].array.F[ii];
        for(ii=0; ii<pXsize*pYsize; ii++)
            data.image[ID].array.F[ii] /= tot;
        ave =  tot/pXsize/pYsize;
        printf("tot = %f   ave = %f \n", tot, ave);

      

        if(ave > 10.0)
        {
            /** compute offset */
            fft_correlation("imwfs", "imref", "outcorr");
            IDc = image_ID("outcorr");
            list_image_ID();
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

            save_fits("outcorr", "!./tmp/outcorr.fits");
            delete_image_ID("outcorr");

            printf("  %6.4f  x  %6.4f\n", totx, toty);

            stepx = (long) (-gain*totx*PcamPixScaleAct); // 0.7*10000.0);
            stepy = (long) (gain*toty*PcamPixScaleAct); //  0.7*10000.0);

            if(stepx>maxstep)
                stepx = maxstep;
            if(stepx<-maxstep)
                stepx = -maxstep;
            if(stepy>maxstep)
                stepy = maxstep;
            if(stepy<-maxstep)
                stepy = -maxstep;


            printf("STEP     : %ld %ld\n", stepx, stepy);

            SCExAO_Pcam_Xpos += stepx;
            SCExAO_Pcam_Ypos += stepy;

            if (SCExAO_Pcam_Xpos>SCExAO_Pcam_Xpos0+SCExAO_Pcam_Range)
                SCExAO_Pcam_Xpos = SCExAO_Pcam_Xpos0+SCExAO_Pcam_Range;
            if (SCExAO_Pcam_Ypos>SCExAO_Pcam_Ypos0+SCExAO_Pcam_Range)
                SCExAO_Pcam_Ypos = SCExAO_Pcam_Ypos0+SCExAO_Pcam_Range;

            if (SCExAO_Pcam_Xpos<SCExAO_Pcam_Xpos0-SCExAO_Pcam_Range)
                SCExAO_Pcam_Xpos = SCExAO_Pcam_Xpos0-SCExAO_Pcam_Range;
            if (SCExAO_Pcam_Ypos<SCExAO_Pcam_Ypos0-SCExAO_Pcam_Range)
                SCExAO_Pcam_Ypos = SCExAO_Pcam_Ypos0-SCExAO_Pcam_Range;

            /// write stages position
            fp = fopen("./status/pcampos.txt", "w");
            fprintf(fp, "%ld %ld\n", SCExAO_Pcam_Xpos, SCExAO_Pcam_Ypos);
            fclose(fp);

            sprintf(command, "pywfs reimage x goto %ld\n", SCExAO_Pcam_Xpos);
            printf("%s", command);
            r = system(command);
            usleep(delayus);

            sprintf(command, "pywfs reimage y goto %ld\n", SCExAO_Pcam_Ypos);
            printf("%s", command);
            r = system(command);
            usleep(delayus);
        }
        else
        {
            printf("Not enough light on detector... waiting... \n");
        }
    }
    r = system("rm stop_PyAlignCam.txt");

    return(0);
}



/// pupil centering tool
/// watch pcenter stream
int SCExAOcontrol_PyramidWFS_Pcenter(char *IDwfsname, float prad, float poffset)
{
    long IDmask;
    long IDwfs;
    long size;
    float centobs = 0.3;
    long ii, jj;
    float x, y, r;
    long *sizearray;
    long ID;
    long size2;
    long cnt;
    float voltAmpOffset = 2.0;
    char command[200];
    long NBframes = 10000;
    long IDpp, IDpm, IDmp, IDmm;
    long IDpyrTTref;
    float xcmm, ycmm, xcpm, ycpm, xcmp, ycmp, xcpp, ycpp;
    float flim;
    float totmm, totpm, totmp, totpp;
    float p10, p90;
    long ii1, jj1;
    float xave, yave;
    FILE *fp;
    long delayus = 1000000;
  

    IDwfs = image_ID(IDwfsname);
    size = data.image[IDwfs].md[0].size[0];
    size2 = size*size;

    sizearray = (long*) malloc(sizeof(long)*2);
    sizearray[0] = size;
    sizearray[1] = size;


    // Read reference pupil illumination
    IDpyrTTref = read_sharedmem_image("pyrTT");
    SCExAO_PZT_STAGE_Xpos_ref = data.image[IDpyrTTref].array.F[0];
    SCExAO_PZT_STAGE_Ypos_ref = data.image[IDpyrTTref].array.F[1];
    printf("X = %f   Y = %f\n", SCExAO_PZT_STAGE_Xpos_ref, SCExAO_PZT_STAGE_Ypos_ref);

    // + +
    SCExAO_PZT_STAGE_Xpos = SCExAO_PZT_STAGE_Xpos_ref + voltAmpOffset;
    SCExAO_PZT_STAGE_Ypos = SCExAO_PZT_STAGE_Ypos_ref;
    //sprintf(command, "analog_output.py voltage D %5.3f\n", SCExAO_PZT_STAGE_Xpos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput D %5.3f", SCExAO_PZT_STAGE_Xpos);
    printf("COMMAND: \"%s\"\n", command);
    r = system(command);
    //sprintf(command, "analog_output.py voltage C %5.3f\n", SCExAO_PZT_STAGE_Ypos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput C %5.3f", SCExAO_PZT_STAGE_Ypos);
    printf("COMMAND: \"%s\"\n", command);
    r = system(command);
    IDpp = SCExAOcontrol_Average_image(IDwfsname, NBframes, "imwfspp", 4);
    save_fits("imwfspp", "!imwfspp.fits");


    // + -
    SCExAO_PZT_STAGE_Xpos = SCExAO_PZT_STAGE_Xpos_ref + voltAmpOffset;
    SCExAO_PZT_STAGE_Ypos = SCExAO_PZT_STAGE_Ypos_ref + 2.0*voltAmpOffset;
    //sprintf(command, "analog_output.py voltage D %5.3f\n", SCExAO_PZT_STAGE_Xpos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput D %5.3f", SCExAO_PZT_STAGE_Xpos);     
    printf("COMMAND: \"%s\"\n", command);
    r = system(command);
    //sprintf(command, "analog_output.py voltage C %5.3f\n", SCExAO_PZT_STAGE_Ypos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput C %5.3f", SCExAO_PZT_STAGE_Ypos);
    printf("COMMAND: \"%s\"\n", command);
    r = system(command);
    IDpm = SCExAOcontrol_Average_image(IDwfsname, NBframes, "imwfspm", 4);
    save_fits("imwfspm", "!imwfspm.fits");

    // - +
    SCExAO_PZT_STAGE_Xpos = SCExAO_PZT_STAGE_Xpos_ref - voltAmpOffset;
    SCExAO_PZT_STAGE_Ypos = SCExAO_PZT_STAGE_Ypos_ref - 2.0*voltAmpOffset;
    //sprintf(command, "analog_output.py voltage D %5.3f\n", SCExAO_PZT_STAGE_Xpos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput D %5.3f", SCExAO_PZT_STAGE_Xpos);
	printf("COMMAND: \"%s\"\n", command);
    r = system(command);
    //sprintf(command, "analog_output.py voltage C %5.3f\n", SCExAO_PZT_STAGE_Ypos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput C %5.3f", SCExAO_PZT_STAGE_Ypos);
    printf("COMMAND: \"%s\"\n", command);
    r = system(command);
    IDmp = SCExAOcontrol_Average_image(IDwfsname, NBframes, "imwfsmp", 4);
    save_fits("imwfsmp", "!imwfsmp.fits");


    // - -
    SCExAO_PZT_STAGE_Xpos = SCExAO_PZT_STAGE_Xpos_ref - voltAmpOffset;
    SCExAO_PZT_STAGE_Ypos = SCExAO_PZT_STAGE_Ypos_ref;
    //sprintf(command, "analog_output.py voltage D %5.3f\n", SCExAO_PZT_STAGE_Xpos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput D %5.3f", SCExAO_PZT_STAGE_Xpos);
	printf("COMMAND: \"%s\"\n", command);
    r = system(command);
    //sprintf(command, "analog_output.py voltage C %5.3f\n", SCExAO_PZT_STAGE_Ypos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput C %5.3f", SCExAO_PZT_STAGE_Ypos);
	printf("COMMAND: \"%s\"\n", command); 
    r = system(command);
    IDmm = SCExAOcontrol_Average_image(IDwfsname, NBframes, "imwfsmm", 4);
    save_fits("imwfsmm", "!imwfsmm.fits");


    // going back to reference

    SCExAO_PZT_STAGE_Xpos = SCExAO_PZT_STAGE_Xpos_ref;
    SCExAO_PZT_STAGE_Ypos = SCExAO_PZT_STAGE_Ypos_ref;
    //sprintf(command, "analog_output.py voltage D %5.3f\n", SCExAO_PZT_STAGE_Xpos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput D %5.3f", SCExAO_PZT_STAGE_Xpos);
	printf("COMMAND: \"%s\"\n", command);
    r = system(command);
    //sprintf(command, "analog_output.py voltage C %5.3f\n", SCExAO_PZT_STAGE_Ypos);
    sprintf(command, "./aocustomscripts/SCExAO_analogoutput C %5.3f", SCExAO_PZT_STAGE_Ypos);
	printf("COMMAND: \"%s\"\n", command);
    r = system(command);


    // sum the 4 images
    arith_image_add("imwfspp", "imwfspm", "prefsum");
    arith_image_add_inplace("prefsum", "imwfsmp");
    arith_image_add_inplace("prefsum", "imwfsmm");

    delete_image_ID("imwfspp");
    delete_image_ID("imwfspm");
    delete_image_ID("imwfsmp");
    delete_image_ID("imwfsmm");

    save_fits("prefsum", "!prefsum.fits");
    p10 = img_percentile("prefsum", 0.10);
    p90 = img_percentile("prefsum", 0.90);

    xcpp = 0.0;
    ycpp = 0.0;
    totpp = 0.0;

    xcmp = 0.0;
    ycmp = 0.0;
    totmp = 0.0;

    xcpm = 0.0;
    ycpm = 0.0;
    totpm = 0.0;

    xcmm = 0.0;
    ycmm = 0.0;
    totmm = 0.0;

    flim = p10 + 0.3*(p90-p10);
    ID = image_ID("prefsum");
    for(ii=0; ii<size/2; ii++)
        for(jj=0; jj<size/2; jj++)
        {
            if(data.image[ID].array.F[jj*size+ii]>flim)
            {
                totmm += 1.0;
                xcmm += ii;
                ycmm += jj;
            }

            ii1 = ii+size/2;
            jj1 = jj;
            if(data.image[ID].array.F[jj1*size+ii1]>flim)
            {
                totpm += 1.0;
                xcpm += ii;
                ycpm += jj;
            }

            ii1 = ii;
            jj1 = jj+size/2;
            if(data.image[ID].array.F[jj1*size+ii1]>flim)
            {
                totmp += 1.0;
                xcmp += ii;
                ycmp += jj;
            }

            ii1 = ii+size/2;
            jj1 = jj+size/2;
            if(data.image[ID].array.F[jj1*size+ii1]>flim)
            {
                totpp += 1.0;
                xcpp += ii;
                ycpp += jj;
            }

        }

    xcpp /= totpp;
    ycpp /= totpp;

    xcpm /= totpm;
    ycpm /= totpm;

    xcmp /= totmp;
    ycmp /= totmp;

    xcmm /= totmm;
    ycmm /= totmm;



    printf("++ : %f %f\n", xcpp, ycpp);
    printf("+- : %f %f\n", xcpm, ycpm);
    printf("-+ : %f %f\n", xcmp, ycmp);
    printf("-- : %f %f\n", xcmm, ycmm);

    xave = 0.25*(xcpp+xcpm+xcmp+xcmm);
    yave = 0.25*(ycpp+ycpm+ycmp+ycmm);

    xave = xave - 0.25*size + 0.5;
    yave = yave - 0.25*size + 0.5;
    printf("AVERAGE PIXEL OFFSET = %f %f\n", xave, yave);

    delete_image_ID("prefsum");


    /// read position of stages
    if((fp = fopen("./status/pcampos.txt", "r"))!=NULL)
    {
        r = fscanf(fp, "%ld %ld\n", &SCExAO_Pcam_Xpos, &SCExAO_Pcam_Ypos);
        printf("CURRENT POSITION : %ld %ld\n", SCExAO_Pcam_Xpos, SCExAO_Pcam_Ypos);
        fclose(fp);
    }

    SCExAO_Pcam_Xpos += (long) (xave*PcamPixScaleAct);
    SCExAO_Pcam_Ypos -= (long) (yave*PcamPixScaleAct);
    printf("NEW POSITION : %ld %ld\n", SCExAO_Pcam_Xpos, SCExAO_Pcam_Ypos);

    xcpp -= xave;
    xcpm -= xave;
    xcmp -= xave;
    xcmm -= xave;

    ycpp -= yave;
    ycpm -= yave;
    ycmp -= yave;
    ycmm -= yave;

    /// write stages position
    fp = fopen("./status/pcampos.txt", "w");
    fprintf(fp, "%ld %ld\n", SCExAO_Pcam_Xpos, SCExAO_Pcam_Ypos);
    fclose(fp);

    sprintf(command, "pywfs reimage x goto %ld\n", SCExAO_Pcam_Xpos);
    printf("%s", command);
    r = system(command);
    usleep(delayus);

    sprintf(command, "pywfs reimage y goto %ld\n", SCExAO_Pcam_Ypos);
    printf("%s", command);
    r = system(command);
    usleep(delayus);


    ID = create_image_ID("pcenter", 2, sizearray, FLOAT, 1, 0);


    IDmask = create_2Dimage_ID("pmask", size, size);

    for(ii=0; ii<size/2; ii++)
        for(jj=0; jj<size/2; jj++)
        {
            // --
            ii1 = ii;
            jj1 = jj;
            x = xcmm-ii;
            y = ycmm-jj;
            r = sqrt(x*x+y*y);
            r /= prad;
            if((r>centobs)&&(r<1.0))
                data.image[IDmask].array.F[jj1*size+ii1] = 1.0;

            // +-
            ii1 = ii+size/2;
            jj1 = jj;
            x = xcpm-ii;
            y = ycpm-jj;
            r = sqrt(x*x+y*y);
            r /= prad;
            if((r>centobs)&&(r<1.0))
                data.image[IDmask].array.F[jj1*size+ii1] = 1.0;

            // -+
            ii1 = ii;
            jj1 = jj+size/2;
            x = xcmp-ii;
            y = ycmp-jj;
            r = sqrt(x*x+y*y);
            r /= prad;
            if((r>centobs)&&(r<1.0))
                data.image[IDmask].array.F[jj1*size+ii1] = 1.0;

            // ++
            ii1 = ii+size/2;
            jj1 = jj+size/2;
            x = xcpp-ii;
            y = ycpp-jj;
            r = sqrt(x*x+y*y);
            r /= prad;
            if((r>centobs)&&(r<1.0))
                data.image[IDmask].array.F[jj1*size+ii1] = 1.0;
        }



    printf("Applying mask to image ...\n");
    fflush(stdout);
    cnt = data.image[IDwfs].md[0].cnt0;
    while (1)
    {
        usleep(10);
        if(cnt != data.image[IDwfs].md[0].cnt0)
        {
            cnt = data.image[IDwfs].md[0].cnt0;

            data.image[ID].md[0].write = 1;
            for(ii=0; ii<size2; ii++)
                data.image[ID].array.F[ii] = 1.0*data.image[IDwfs].array.U[ii]*data.image[IDmask].array.F[ii];
            data.image[ID].md[0].cnt0++;
            data.image[ID].md[0].write = 1;
        }
    }

    free(sizearray);

    return(0);
}





int SCExAOcontrol_Pyramid_flattenRefWF(char *WFScam_name, long zimaxmax, float ampl0)
{
    long zimax;
    long zi;
    long ID;
    long NBframes = 100;
    double val, valp, valm, val0;
    double ampl;
    double a;
    char command[200];
    int r;
    long IDdm5, IDdm6;
    long ii;
    long dmsize;
    long dmsize2;
    long IDz;
    long IDdisp;
    long sleeptimeus = 1000; // 1ms


    // 60perc of pixels illuminated
    // perc 70 is median over pupil

    double p0, p1, p2;
    float level0, level1, level2;

    level0 = 0.55;
    level1 = 0.70;
    level2 = 0.95;
    

    IDdm5 = read_sharedmem_image("dmdisp5");
    IDdm6 = read_sharedmem_image("dmdisp6");
    IDdisp = read_sharedmem_image("dmdisp");
    dmsize = data.image[IDdm5].md[0].size[0];
    dmsize2 = dmsize*dmsize;

    // prepare modes
    IDz = mk_zer_seriescube("zcube", dmsize, zimaxmax, 0.5*dmsize);
    list_image_ID();
    printf("IDz = %ld\n", IDz);

    zimax = zimaxmax;


    
            ID = SCExAOcontrol_Average_image(WFScam_name, NBframes, "imwfs", 6);
            save_fits("imwfs", "!./tmp/imwfs_pyrflat.fits");
            p0 = img_percentile("imwfs", level0);
            p1 = img_percentile("imwfs", level1);
            p2 = img_percentile("imwfs", level2);
            val = (p2-p0)/p1; //+p90);
            printf("%lf %lf %lf -> %f\n", p0, p1, p2, val);
            val0 = val;

    
    ampl = ampl0;

    while(1)
    {
        ampl *= 0.95;
            if(ampl<0.1*ampl0)
                ampl = 0.1*ampl0;

        
        for(zi=4; zi<zimax; zi++)
        {
//            ampl = ampl0; //*pow((1.0 - 0.9*(zimax/zimaxmax)), 2.0);

            data.image[IDdm5].md[0].write = 1;
            for(ii=0; ii<dmsize2; ii++)
                data.image[IDdm5].array.F[ii] += ampl*data.image[IDz].array.F[zi*dmsize2+ii];
            sem_post(data.image[IDdm5].semptr[0]);
            sem_post(data.image[IDdisp].semptr[1]);
            data.image[IDdm5].md[0].cnt0++;
            data.image[IDdm5].md[0].write = 0;
           
            usleep(sleeptimeus);


            ID = SCExAOcontrol_Average_image(WFScam_name, NBframes, "imwfs", 6);
            save_fits("imwfs", "!./tmp/imwfs_pyrflat.fits");
            p0 = img_percentile("imwfs", level0);
            p1 = img_percentile("imwfs", level1);
            p2 = img_percentile("imwfs", level2);
            val = (p2-p0)/p1; //+p90);
            printf("%lf %lf %lf -> %f\n", p0, p1, p2, val);
            valp = val;


            data.image[IDdm5].md[0].write = 1;
            for(ii=0; ii<dmsize2; ii++)
                data.image[IDdm5].array.F[ii] -= 2.0*ampl*data.image[IDz].array.F[zi*dmsize2+ii];
            sem_post(data.image[IDdm5].semptr[0]);
            sem_post(data.image[IDdisp].semptr[1]);
            data.image[IDdm5].md[0].cnt0++;
            data.image[IDdm5].md[0].write = 0;
           
            usleep(sleeptimeus);


            ID = SCExAOcontrol_Average_image(WFScam_name, NBframes, "imwfs", 6);
            save_fits("imwfs", "!./tmp/imwfs_pyrflat.fits");
            p0 = img_percentile("imwfs", level0);
            p1 = img_percentile("imwfs", level1);
            p2 = img_percentile("imwfs", level2);
            val = (p2-p0)/p1; //+p90);
            printf("%lf %lf %lf -> %f\n", p0, p1, p2, val);
            valm = val;

            /*	if(valm>valp)
            		a = -amp;
            	else
            		a = amp;
            */

            a = (1.0/valp-1.0/valm)/(1.0/valp+1.0/valm)*ampl;
            printf("== ZERNIKE %ld / %ld ========== %f %f -> a = %f  [ampl = %f] ( %f <- %f)\n", zi, zimax, valp, valm, a, ampl, 0.5*(valp+valm), val0);


            data.image[IDdm5].md[0].write = 1;
            for(ii=0; ii<dmsize2; ii++)
                data.image[IDdm5].array.F[ii] += (ampl+a)*data.image[IDz].array.F[zi*dmsize2+ii];
            sem_post(data.image[IDdm5].semptr[0]);
            sem_post(data.image[IDdisp].semptr[1]);
            data.image[IDdm5].md[0].cnt0++;
            data.image[IDdm5].md[0].write = 0;

            usleep(sleeptimeus);

            //			sleep(1);
        }

        printf("%ld -> %ld\n", IDdm5, IDdm6);
        data.image[IDdm5].md[0].write = 1;
        data.image[IDdm6].md[0].write = 1;
        for(ii=0; ii<2500; ii++)
        {
            data.image[IDdm6].array.F[ii] += data.image[IDdm5].array.F[ii];
            data.image[IDdm5].array.F[ii] = 0.0;
        }
        data.image[IDdm5].md[0].cnt0++;
        data.image[IDdm6].md[0].cnt0++;
        data.image[IDdm5].md[0].write = 0;
        data.image[IDdm6].md[0].write = 0;

        zimax ++;
        if(zimax>zimaxmax)
            zimax = zimaxmax;
    }

    return(0);
}







int SCExAOcontrol_optPSF(char *WFScam_name, long NBmodesmax, float alpha)
{
	FILE *fp;
    long NBmodes;
    long mode;
    long modestart = 1;
    long ID;
    long NBframes = 20;
    double val, valp, valm, val0;
    double ampl;
    double a;
    char command[200];
    int r;
    long IDdm5, IDdm6;
    long ii, jj;
    long dmsize;
    long dmsize2;
    long IDm;
    long IDdisp;
    long sleeptimeus = 100000; // 100ms

	long iter = 0;
    double p0, p1, p2;
    float level0, level1, level2;
	double v0;
	double tot, tot1;
	double pv;
	long xsize, ysize;
	long cnt;
	long iter1 = 0;
	
	double ampl0 = 0.005;


	double beta = 0.99999; // regularization

	double ampcoeff = 10.0;


    level0 = 0.1;
    level1 = 0.2;
    level2 = 0.4;
    

    IDdm5 = read_sharedmem_image("dm00disp05");
    IDdm6 = read_sharedmem_image("dm00disp06");
    IDdisp = read_sharedmem_image("dm00disp");
    dmsize = data.image[IDdm5].md[0].size[0];
    dmsize2 = dmsize*dmsize;

    // prepare modes
    
    IDm = image_ID("modes");
    
    if(IDm==-1)
	{    
		IDm = mk_zer_seriescube("modes", dmsize, NBmodesmax, 0.5*dmsize);
		list_image_ID();
		printf("IDm = %ld\n", IDm);
		NBmodes = NBmodesmax;
	}
	else
		NBmodes = data.image[IDm].md[0].size[2];

	if(NBmodes>NBmodesmax)
		NBmodes = NBmodesmax;



	printf("Averaging %ld frames ,,,", NBframes);
	fflush(stdout);
            ID = SCExAOcontrol_Average_image(WFScam_name, NBframes, "impsf", 6);
    printf(" done\n");
    fflush(stdout);
    
            save_fits("impsf", "!./tmp/impsf.fits");
            p0 = img_percentile("impsf", level0);
            p1 = img_percentile("impsf", level1);
            p2 = img_percentile("impsf", level2);
 
	xsize = data.image[ID].md[0].size[0];
	ysize = data.image[ID].md[0].size[1];

    tot1 = 0.0;
    tot = 0.0;
    cnt = 0;
	v0 = 1.0*p2;
    for(ii=0;ii<xsize;ii++)
		for(jj=0;jj<ysize;jj++)
			{
				pv = data.image[ID].array.F[jj*xsize+ii];
				pv -= v0;
				if(pv>0.0)
					{
						tot1 += pow(pv, alpha);
						tot += pv;
						cnt++;
					}
			}
    val = tot1/pow(tot, alpha);
    val0 = val;
    
    
    
    
    ampl = ampl0;

    while(1)
    {
        ampl *= 0.95;
            if(ampl<0.1*ampl0)
                ampl = 0.1*ampl0;

        
        for(mode=modestart; mode<NBmodes; mode++)
        {
            data.image[IDdm5].md[0].write = 1;
            for(ii=0; ii<dmsize2; ii++)
                data.image[IDdm5].array.F[ii] += ampl*data.image[IDm].array.F[mode*dmsize2+ii];
            COREMOD_MEMORY_image_set_sempost_byID(IDdm5, -1);
            sem_post(data.image[IDdisp].semptr[1]);
            data.image[IDdm5].md[0].cnt0++;
            data.image[IDdm5].md[0].write = 0;
           
            usleep(sleeptimeus);


            ID = SCExAOcontrol_Average_image(WFScam_name, NBframes, "impsf", 6);
            
            save_fits("impsf", "!./tmp/impsf.fits");
            p0 = img_percentile("impsf", level0);
            p1 = img_percentile("impsf", level1);
            p2 = img_percentile("impsf", level2);
            tot1 = 0.0;
			tot = 0.0;
			cnt = 0;
			v0 = 1.0*p2;
			for(ii=0;ii<xsize;ii++)
				for(jj=0;jj<ysize;jj++)
				{
				pv = data.image[ID].array.F[jj*xsize+ii];
				pv -= v0;
				if(pv>0.0)
					{
						tot1 += pow(pv, alpha);
						tot += pv;
						cnt++;
					}
				}
			val = tot1/pow(tot, alpha);
               
            valp = val;


            data.image[IDdm5].md[0].write = 1;
            for(ii=0; ii<dmsize2; ii++)
                data.image[IDdm5].array.F[ii] -= 2.0*ampl*data.image[IDm].array.F[mode*dmsize2+ii];            
            COREMOD_MEMORY_image_set_sempost_byID(IDdm5, -1);
            sem_post(data.image[IDdisp].semptr[1]);
            data.image[IDdm5].md[0].cnt0++;
            data.image[IDdm5].md[0].write = 0;
           
            usleep(sleeptimeus);


            ID = SCExAOcontrol_Average_image(WFScam_name, NBframes, "impsf", 6);
            
            save_fits("impsf", "!./tmp/impsf.fits");
            p0 = img_percentile("impsf", level0);
            p1 = img_percentile("impsf", level1);
            p2 = img_percentile("impsf", level2);
            tot1 = 0.0;
			tot = 0.0;
			cnt = 0;
			v0 = 1.0*p2;
			for(ii=0;ii<xsize;ii++)
				for(jj=0;jj<ysize;jj++)
				{
				pv = data.image[ID].array.F[jj*xsize+ii];
				pv -= v0;
				if(pv>0.0)
					{
						tot1 += pow(pv, alpha);
						tot += pv;
						cnt++;
					}
				}
			val = tot1/pow(tot, alpha);
			if(cnt<1)
				{
					val = 0.0;
					printf("ERROR: image is zero \n");
					exit(0);
				}
			printf("=========== %f  %f  %ld  %f\n", tot, tot1, cnt, val);
			
            valm = val;


           
             a = ampcoeff * (valp - valm) / (valp + valm)*ampl;
            if(a<-ampl)
				a = -ampl;
			if(a>ampl)
				a = ampl;
            
            printf("==  MODE %ld / %ld ========== (%ld) %f %f -> a = %f  [ampl = %f] ( %f <- %f)\n", mode, NBmodes, cnt, valp, valm, a, ampl, 0.5*(valp+valm), val0);

			fp = fopen("log.txt", "a");
			fprintf(fp, "%8ld  %8ld  %4ld  %20f  %20f\n", iter1, iter, mode, 0.5*(valp+valm), val0);
			fclose(fp);
			
            data.image[IDdm5].md[0].write = 1;
            for(ii=0; ii<dmsize2; ii++)
                data.image[IDdm5].array.F[ii] += (ampl+a)*data.image[IDm].array.F[mode*dmsize2+ii];
            COREMOD_MEMORY_image_set_sempost_byID(IDdm5, -1);
            sem_post(data.image[IDdisp].semptr[1]);
            data.image[IDdm5].md[0].cnt0++;
            data.image[IDdm5].md[0].write = 0;


            data.image[IDdm5].md[0].write = 1;
            for(ii=0; ii<dmsize2; ii++)
                data.image[IDdm5].array.F[ii] *= beta;
            COREMOD_MEMORY_image_set_sempost_byID(IDdm5, -1);
            sem_post(data.image[IDdisp].semptr[1]);
            data.image[IDdm5].md[0].cnt0++;
            data.image[IDdm5].md[0].write = 0;


            usleep(sleeptimeus);

            iter1++;
        }
        ampcoeff *= 0.9;
        if(ampcoeff<1.0)
			ampcoeff = 1.0;



        printf("%ld -> %ld\n", IDdm5, IDdm6);
		fflush(stdout);

		
		
        data.image[IDdm5].md[0].write = 1;
        data.image[IDdm6].md[0].write = 1;
        for(ii=0; ii<dmsize2; ii++)
        {
            data.image[IDdm6].array.F[ii] += data.image[IDdm5].array.F[ii];
            data.image[IDdm5].array.F[ii] = 0.0;
        }
        COREMOD_MEMORY_image_set_sempost_byID(IDdm5, -1);
        COREMOD_MEMORY_image_set_sempost_byID(IDdm6, -1);
        data.image[IDdm5].md[0].cnt0++;
        data.image[IDdm6].md[0].cnt0++;
        data.image[IDdm5].md[0].write = 0;
        data.image[IDdm6].md[0].write = 0;

		

        NBmodes ++;
        if(NBmodes>NBmodesmax)
            NBmodes = NBmodesmax;
        if(NBmodes>data.image[IDm].md[0].size[2])
			NBmodes = data.image[IDm].md[0].size[2];
    
		
    
    
		iter++;
    }

    return(0);
}





/** SAPHIRA image: process data cube into single frame
 *
 * full linear regression, up to saturation level
 *
 * */
int SCExAOcontrol_SAPHIRA_cam_process(char *IDinname, char *IDoutname)
{
    long IDout;
    long IDin;
    long xsize, ysize, zsize;
    long *sizeoutarray;
    long k;
    long ii, jj;
    float v0;
    long xysize;
    double v1, vk, vt, vv;
    long kk;
    long cnt0, cnt1, cnt2;
    int SATURATION = 65534;
    long iter;
    double eps = 1e-8;
    long k1;
    long IDintmp, IDsatmask, ID2dtmp;
    unsigned short int pvu, vcnt, vaveku;
    long vavevu;
    float vavev, vavek;
    long k1start;

    IDin = image_ID(IDinname);

    xsize = data.image[IDin].md[0].size[0];
    ysize = data.image[IDin].md[0].size[1];
    zsize = data.image[IDin].md[0].size[2];
    xysize = xsize*ysize;


    if(zsize>2)
        k1start = 1;
    else
        k1start = 0;

    sizeoutarray = (long*) malloc(sizeof(long)*3);
    sizeoutarray[0] = xsize;
    sizeoutarray[1] = ysize;
    sizeoutarray[2] = zsize;



    IDintmp = create_image_ID("intmp", 3, sizeoutarray, USHORT, 1, 0); // temporary buffer
    IDsatmask = create_image_ID("satmask", 3, sizeoutarray, USHORT, 1, 0); // saturation mask
    ID2dtmp = create_image_ID("saphira2dtmp", 2, sizeoutarray, FLOAT, 1, 0); // intermediate resutl
    IDout = create_image_ID(IDoutname, 2, sizeoutarray, FLOAT, 1, 0);
    COREMOD_MEMORY_image_set_createsem(IDoutname, 4);

    if(data.image[IDin].sem == 0)
    {
        printf("Error: no semaphore detected\n");
        exit(0);
    }


    // drive semaphore to zero
    while(sem_trywait(data.image[IDin].semptr[0])==0) {}



    iter = 0;


    while(1)
    {
        sem_wait(data.image[IDin].semptr[0]);
        while(sem_trywait(data.image[IDin].semptr[0])==0) {}

        k = data.image[IDin].md[0].cnt1;
        printf("%ld   slice %ld written [%ld] \n      ", iter, k, IDin);
        fflush(stdout);

        if(k == zsize-1)  // process cube
        {
            memcpy(data.image[IDintmp].array.U, data.image[IDin].array.U, sizeof(short)*xysize*zsize);
            for(ii=0; ii<xysize; ii++)
            {
                k1 = 0;
                v0 = 0.0;
                v1 = 0.0;
                vcnt = 0;
                vaveku = 0;
                vavevu = 0;
                for(k1=k1start; k1<zsize; k1++)
                {
                    pvu = data.image[IDintmp].array.U[k1*xysize+ii];
                    //	printf("[%d %u] ", pvu, pvu);
                    if(pvu<SATURATION)
                    {
                        data.image[IDsatmask].array.U[k1*xysize+ii] = 1;
                        vavevu += pvu;
                        vaveku += k1;
                        vcnt++;
                    }
                    else
                    {
                        data.image[IDsatmask].array.U[k1*xysize+ii] = 0;
                    }
                }
                vavev = 1.0*vavevu/vcnt;
                vavek = 1.0*vaveku/vcnt;

                for(k1=k1start; k1<zsize; k1++)
                {
                    pvu = data.image[IDintmp].array.U[k1*xysize+ii];
                    if(data.image[IDsatmask].array.U[k1*xysize+ii] == 1)
                    {
                        vk = 1.0*k1 - vavek;
                        vv = 1.0*pvu - vavev;
                        v0 += vk*vv;
                        v1 += vk*vk;
                    }
                }
                data.image[ID2dtmp].array.F[ii] = v0/(v1+eps);

            }

            iter++;
            printf("\n CUBE COMPLETED -> 2D image ready\n");
            data.image[IDout].md[0].write = 1;
            memcpy(data.image[IDout].array.F, data.image[ID2dtmp].array.F, sizeof(float)*xysize);
            if(data.image[IDout].sem > 0)
                sem_post(data.image[IDout].semptr[0]);
            data.image[IDout].md[0].cnt0 ++;
            data.image[IDout].md[0].write = 0;
        }
    }

    free(sizeoutarray);

    return(IDout);
}





