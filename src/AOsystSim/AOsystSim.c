#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "image_filter/image_filter.h"
#include "fft/fft.h"
#include "info/info.h"
#include "statistic/statistic.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "image_filter/image_filter.h"

#include "OptSystProp/OptSystProp.h"
#include "AOsystSim/AOsystSim.h"

extern DATA data;


// non-null pixels in DM influence function (to speed up computing time)
int DMifpixarray_init = 0;
float *DMifpixarray_val;
long *DMifpixarray_index; // which actuator
long *DMifpixarray_pixindex; // which pixel
long DMifpixarray_NBpix;


OPTSYST *optsystsim;


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


int AOsystSim_simpleAOfilter_cli()
{
  if(CLI_checkarg(1,3)+CLI_checkarg(1,3)==0)
    AOsystSim_simpleAOfilter(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);

  return(0);
}



int AOsystSim_fitTelPup_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)==0)
    AOsystSim_fitTelPup(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);

  return(0);
}



int AOsystSim_run_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)==0)
    {
        AOsystSim_run(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl);
        return 0;
    }
    else
        return 1;
}





int init_AOsystSim()
{
    strcpy(data.module[data.NBmodule].name, __FILE__);
    strcpy(data.module[data.NBmodule].info, "conversion between image format, I/O");
    data.NBmodule++;

    strcpy(data.cmd[data.NBcmd].key,"AOsimfilt");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOsystSim_simpleAOfilter_cli;
    strcpy(data.cmd[data.NBcmd].info,"AO simple filtering");
    strcpy(data.cmd[data.NBcmd].syntax,"<input WF> <output WF>");
    strcpy(data.cmd[data.NBcmd].example,"AOsimfilt wfin wfout");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOsystSim_simpleAOfilter(char *IDin_name, char *IDout_name)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"AOsystsfitpup");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOsystSim_fitTelPup_cli;
    strcpy(data.cmd[data.NBcmd].info,"fit telescope pupil");
    strcpy(data.cmd[data.NBcmd].syntax,"tel pupil file");
    strcpy(data.cmd[data.NBcmd].example,"AOsystfitpup");
    strcpy(data.cmd[data.NBcmd].Ccall,"AOsystSim_fitTelPup(char *ID_name)");
    data.NBcmd++;


    strcpy(data.cmd[data.NBcmd].key,"AOsystsim");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOsystSim_run_cli;
    strcpy(data.cmd[data.NBcmd].info,"run fake AO system");
    strcpy(data.cmd[data.NBcmd].syntax,"no argument");
    strcpy(data.cmd[data.NBcmd].example,"AOsystsim");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOsystSim_run()");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"AOsystexaosim");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOsystSim_extremeAO_contrast_sim;
    strcpy(data.cmd[data.NBcmd].info,"run extremeAO analysis");
    strcpy(data.cmd[data.NBcmd].syntax,"no argument");
    strcpy(data.cmd[data.NBcmd].example,"AOsystexaosim");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOsystSim_extremeAO_contrast_sim()");
    data.NBcmd++;


    // add atexit functions here

    return 0;
}





/** \brief simple AO filtering model using Fourier analysis
 *         simulates WFS integration, delay, noise (as a function of spatial frequency)
 *
 *  Open loop model
 * 
 *  AO simple filter is modeled by several 2D maps in spatial frequency:
 *  aosf_noise : noise per spatial frequency
 *  aosf_mult : signal throughput per spatial frequency
 *  aosf_gain : loop gain to be applied per spatial frequency
 */


int AOsystSim_simpleAOfilter(char *IDin_name, char *IDout_name)
{
    long IDin, IDmask, IDout, IDdm, IDwfe;
    long *sizearray;
    long cnt0;
    long ii, jj;
    double x, y, r;
    long ID, IDin1;

    long IDaosf_noise;
    long IDaosf_mult;
    long IDaosf_gain;

    float loopgain = 0.001;
    float multr0 = 0.2;

    float wfsetime = 0.001; /**< WFS exposure time */
    float timedelay = 0.0005;/**< time delay between wavefront sensor measurement (end of exposure) and DM correction [sec] */
    long size2;
    double tnowdouble = 0.0;
    double time_wfse0 = 0.0;
    long wfsecnt = 0;

    double rmask;
    double r1;
    double WFrms0 = 0.0;
    double WFrms = 0.0;
    long WFrmscnt = 0;

    double noiselevel = 0.3; // noise level per WFS measurement [um] - only white noise is currently supported

    /** wavefront correction buffer to implement time delay */
    long IDwfcbuff;
    long wfcbuff_size = 10000;
    double *wfcbuff_time;
    long *wfcbuff_status; // 1 : waiting to be applied
    long k, k0, k1;
    long wfsecnt1;

    double dmmovetime = 0.001; /**< time it takes for DM to move */
    long dmmoveNBpt = 10;

    sizearray = (long*) malloc(sizeof(long)*2);

    IDin = read_sharedmem_image(IDin_name);  /**< turbulence channel */
    sizearray[0] = data.image[IDin].md[0].size[0];
    sizearray[1] = data.image[IDin].md[0].size[1];
    size2 = sizearray[0]*sizearray[1];

    IDwfcbuff = create_3Dimage_ID("wfcbuff", sizearray[0], sizearray[1], wfcbuff_size);
    wfcbuff_time = (double*) malloc(sizeof(double)*wfcbuff_size);
    wfcbuff_status = (long*) malloc(sizeof(long)*wfcbuff_size);
    for(k=0; k<wfcbuff_size; k++)
        wfcbuff_status[k] = 0;

    rmask = 0.4*sizearray[0];

    IDaosf_noise = create_2Dimage_ID("aosf_noise", sizearray[0], sizearray[1]);
    IDaosf_mult = create_2Dimage_ID("aosf_mult", sizearray[0], sizearray[1]);
    IDaosf_gain = create_2Dimage_ID("aosf_gain", sizearray[0], sizearray[1]);

    IDmask = create_2Dimage_ID("aosf_mask", sizearray[0], sizearray[1]);


    for(ii=0; ii<sizearray[0]; ii++)
        for(jj=0; jj<sizearray[1]; jj++)
        {
            x = 1.0*ii-0.5*sizearray[0];
            y = 1.0*jj-0.5*sizearray[1];
            r = sqrt(x*x+y*y);
            r1 = (r-rmask)/(0.5*sizearray[0]-rmask);

            if(r1<0.0)
                data.image[IDmask].array.F[jj*sizearray[0]+ii] = 1.0;
            else if (r1>1.0)
                data.image[IDmask].array.F[jj*sizearray[0]+ii] = 0.0;
            else
                data.image[IDmask].array.F[jj*sizearray[0]+ii] = 0.5*(cos(r1*M_PI)+1.0);

            data.image[IDaosf_gain].array.F[jj*sizearray[0]+ii] = loopgain;
            data.image[IDaosf_mult].array.F[jj*sizearray[0]+ii] = exp(-pow(r/(multr0*sizearray[0]),8.0));
            if(r>multr0*sizearray[0])
                data.image[IDaosf_mult].array.F[jj*sizearray[0]+ii] = 0.0;
            data.image[IDaosf_noise].array.F[jj*sizearray[0]+ii] = 0.0;
        }
    save_fits("aosf_noise", "!aosf_noise.fits");
    save_fits("aosf_mult", "!aosf_mult.fits");
    save_fits("aosf_gain", "!aosf_gain.fits");

    permut("aosf_mult");
    permut("aosf_noise");
    permut("aosf_gain");




    IDdm = create_image_ID("aofiltdm", 2, sizearray, FLOAT, 1, 0);
    IDwfe = create_image_ID("aofiltwfe", 2, sizearray, FLOAT, 1, 0);

    IDout = create_image_ID(IDout_name, 2, sizearray, FLOAT, 1, 0);
    strcpy(data.image[IDout].kw[0].name, "TIME");
    data.image[IDout].kw[0].type = 'D';
    data.image[IDout].kw[0].value.numf = 0.0;
    strcpy(data.image[IDout].kw[0].comment, "Physical time [sec]");

    IDin1 = create_image_ID("aofiltin", 2, sizearray, FLOAT, 1, 0);

    printf("%s -> %s\n", IDin_name, IDout_name);
    cnt0 = -1;
    time_wfse0 = data.image[IDin].kw[0].value.numf; /** start of WFS exposure */
    wfsecnt = 0;
    wfsecnt1 = 0;
    while(1)
    {
        usleep(10);
        if(data.image[IDin].md[0].cnt0!=cnt0)
        {
            /** input masking */
            for(ii=0; ii<size2; ii++)
                data.image[IDin1].array.F[ii] = data.image[IDin].array.F[ii] * data.image[IDmask].array.F[ii];

            /** Wavefront correction & measurement */
            WFrms0 = 0.0;
            WFrms = 0.0;
            WFrmscnt = 0;
            for(ii=0; ii<sizearray[0]*sizearray[1]; ii++)
            {
                data.image[IDout].array.F[ii] = data.image[IDin1].array.F[ii] - data.image[IDdm].array.F[ii];
                if(data.image[IDmask].array.F[ii]>0.99)
                {
                    WFrms0 += data.image[IDin1].array.F[ii]*data.image[IDin1].array.F[ii];
                    WFrms += data.image[IDout].array.F[ii]*data.image[IDout].array.F[ii];
                    WFrmscnt ++;
                }
            }
            WFrms0 = sqrt(WFrms0/WFrmscnt);
            WFrms = sqrt(WFrms/WFrmscnt);
            data.image[IDout].kw[0].value.numf = tnowdouble;


            tnowdouble = data.image[IDin].kw[0].value.numf;
            printf("\r Time : %.6lf    WFSexposure time:  %.6f   [%5ld]    %10.5lf -> %10.5lf ", tnowdouble, tnowdouble-time_wfse0, wfsecnt, WFrms0, WFrms);
            fflush(stdout);
            cnt0 = data.image[IDin].md[0].cnt0;

            do2drfft(IDout_name, "aosf_tmpfft");
            ID = image_ID("aosf_tmpfft");

            for(ii=0; ii<sizearray[0]*sizearray[1]; ii++)
            {
                data.image[ID].array.CF[ii].re *= data.image[IDaosf_mult].array.F[ii]/size2;
                data.image[ID].array.CF[ii].im *= data.image[IDaosf_mult].array.F[ii]/size2;
            }

            do2dffti("aosf_tmpfft", "testo");
            delete_image_ID("aosf_tmpfft");
            ID = image_ID("testo");

            /** Wavefront estimation */
            for(ii=0; ii<sizearray[0]*sizearray[1]; ii++)
                data.image[IDwfe].array.F[ii] += data.image[ID].array.CF[ii].re;
            delete_image_ID("testo");
            wfsecnt ++;

            if((tnowdouble-time_wfse0) > wfsetime)
            {
                if(wfsecnt>0)
                    for(ii=0; ii<sizearray[0]*sizearray[1]; ii++)
                        data.image[IDwfe].array.F[ii] /= wfsecnt;


                /** write entries in buffer */
                for(k1=0; k1<dmmoveNBpt; k1++)
                {
                    k0 = 0;
                    while((wfcbuff_status[k0]==1)&&(k0<wfcbuff_size))
                        k0++;
                    if(k0>wfcbuff_size-1)
                    {
                        printf("\n WFC buffer full [%ld]\n", k0);
                        exit(0);
                    }

                    for(ii=0; ii<size2; ii++)
                        data.image[IDwfcbuff].array.F[k0*size2+ii] = data.image[IDwfe].array.F[ii]/dmmoveNBpt;
                    wfcbuff_time[k0] = tnowdouble + timedelay + 1.0*k1/(dmmoveNBpt-1)*dmmovetime;
                    wfcbuff_status[k0] = 1;
                }

                for(ii=0; ii<sizearray[0]*sizearray[1]; ii++)
                    data.image[IDwfe].array.F[ii] = 0.0;
                time_wfse0 += wfsetime;
                wfsecnt = 0;
                wfsecnt1++;
            }



            for(k=0; k<wfcbuff_size; k++)
            {
                if(wfcbuff_status[k]==1)
                    if(wfcbuff_time[k]<tnowdouble)
                    {
                        //			printf("Reading WFC buffer slice %ld  [%5ld %5ld]\n", k, wfsecnt1, wfsecnt);
                        /** update DM shape */
                        for(ii=0; ii<sizearray[0]*sizearray[1]; ii++)
                            data.image[IDdm].array.F[ii] += data.image[IDaosf_gain].array.F[ii]*data.image[IDwfcbuff].array.F[k*size2+ii];
                        wfcbuff_status[k] = 0;
                    }
            }
            data.image[IDout].md[0].cnt0 = cnt0;
        }
    }

    delete_image_ID("circbuff");

    free(sizearray);
    free(wfcbuff_time);
    free(wfcbuff_status);

    return(0);
}




















// all sizes in actuators on DM

long AOsystSim_mkTelPupDM(char *ID_name, long msize, double xc, double yc, double rin, double rout, double pupPA, double spiderPA, double spideroffset, double spiderthick, double stretchx)
{
    long ID, IDz;
    long ii,  jj;
    double x, y, x1, y1, r;
    long binfact = 8;
    long size;
    double PA;
    double val;
    long IDindex, IDi;
    long index;
    long ii1, jj1;


    size = msize*binfact;

    ID = create_2Dimage_ID(ID_name, msize, msize);
    IDi = create_3Dimage_ID("TPind", msize, msize, 5);

    IDz = create_2Dimage_ID("telpupDMz", size, size);
    IDindex = create_2Dimage_ID("telpupDMzindex", size, size);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            index = 0;
            val = 1.0;

            x = (1.0*ii/binfact);
            y = (1.0*jj/binfact);
            x -= xc;
            y -= yc;
            r = sqrt(x*x+y*y);
            PA = atan2(y, x);
            PA += pupPA;
            x = r*cos(PA)*stretchx;
            y = r*sin(PA);

            if(r>rout)
                val = 0.0;
            if(r<rin)
                val = 0.0;

            x1 = x-spideroffset-spiderthick/2.0;
            y1 = y;
            PA = atan2(y1,x1);
            if(fabs(PA)<spiderPA)
                index = 1;


            x1 = x+spideroffset+spiderthick/2.0;
            y1 = y;
            PA = atan2(y1,-x1);
            if(fabs(PA)<spiderPA)
                index = 2;




            x1 = x+spideroffset-spiderthick/2.0;
            y1 = y;
            PA = atan2(x1,y1);
            if((fabs(PA)<M_PI/2-spiderPA)&&(x<0))
                index = 3;

            x1 = -x+spideroffset-spiderthick/2.0;
            y1 = y;
            PA = atan2(x1,y1);
            if((fabs(PA)<M_PI/2-spiderPA)&&(x>0))
                index = 3;





            x1 = x+spideroffset-spiderthick/2.0;
            y1 = -y;
            PA = atan2(x1,y1);
            if((fabs(PA)<M_PI/2-spiderPA)&&(x<0))
                index = 4;

            x1 = -x+spideroffset-spiderthick/2.0;
            y1 = -y;
            PA = atan2(x1,y1);
            if((fabs(PA)<M_PI/2-spiderPA)&&(x>0))
                index = 4;

            if(index==0)
                val = 0.0;
            data.image[IDz].array.F[jj*size+ii] = val;
            data.image[IDindex].array.F[jj*size+ii] = index*val;

            ii1 = (long) (ii/binfact);
            jj1 = (long) (jj/binfact);


            data.image[ID].array.F[jj1*msize+ii1] += val/binfact/binfact;

            if(val>0.5)
            {
                data.image[IDi].array.F[jj1*msize+ii1] = 1;
                data.image[IDi].array.F[index*msize*msize+jj1*msize+ii1] = 1;
            }
        }

    delete_image_ID("telpupDMz");
    delete_image_ID("telpupDMzindex");

    save_fits("TPind", "!TPind.fits");

    return(ID);
}




/** fits DM illumination to pupil geometry */

long AOsystSim_fitTelPup(char *ID_name, char *IDtelpup_name)
{
    FILE *fp;
    long ID, ID1, IDt;
    long IDtelpup;
    double xc, yc, pupPA, spiderPA, spideroffset, spiderthick, rin, rout, stretchx;
    long size;
    double vp10, vp90;
    double rms;
    long ii, jj;
    double v1;

    double xc_min, yc_min, pupPA_min, rout_min, stretchx_min, spiderPA_min, spideroffset_min;
    double xc_max, yc_max, pupPA_max, rout_max, stretchx_max, spiderPA_max, spideroffset_max;

    double xc1, yc1, pupPA1, rout1, stretchx1, spiderPA1, spideroffset1;
    long nbstep = 2;
    double rms1;

    double xc_range, yc_range, pupPA_range, rout_range, stretchx_range, spiderPA_range, spideroffset_range;
    long iter;
    long NBiter = 5;

    /** compensate for illumination gradient */
    double coeffx, coeffy, coeffx1, coeffy1;
    double x, y;
    double val, val1;

    double eps = 1.0e-8;


    coeffx = 0.0;
    coeffy = 0.0;

    rout = 22.15;
    rin = 0.315*rout;
    pupPA = -0.105;
    spiderPA = 0.85;
    spiderthick = 0.06*rout;
    spideroffset = 0.15*rout;


    xc = 24.629630;
    yc = 23.518519;
    stretchx = 1.048148;

    rout1 = rout;
    xc1 = xc;
    yc1 = yc;
    pupPA1 = pupPA;
    stretchx1 = stretchx;
    spiderPA1 = spiderPA;
    spideroffset1 = spideroffset;

    /** set percentiles */
    ID = image_ID(ID_name);
    size = data.image[ID].md[0].size[0];


    vp10 = img_percentile_float(ID_name, 0.1);
    vp90 = img_percentile_float(ID_name, 0.9);
    printf("%f %f\n", vp10, vp90);
    for(ii=0; ii<data.image[ID].md[0].size[0]*data.image[ID].md[0].size[1]; ii++)
        data.image[ID].array.F[ii] = (data.image[ID].array.F[ii]-vp10)/(vp90-vp10);



    /** compensate for image gradient */
    ID1 = create_2Dimage_ID("tmpftpim", size, size);
    val1 = 1000000000000000.0;
    for(coeffx=-1.5; coeffx < 1.5; coeffx += 0.05)
        for(coeffy=-1.5; coeffy < 1.5; coeffy += 0.05)
        {
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    x = 1.0*ii/size;
                    y = 1.0*jj/size;
                    data.image[ID1].array.F[jj*size+ii] =  data.image[ID].array.F[jj*size+ii]*(1.0+coeffx*(x-0.5))*(1.0+coeffy*(y-0.5));
                }
            val = img_percentile_float("tmpftpim", 0.9) - img_percentile_float("tmpftpim", 0.7);
            if(val<val1)
            {
                val1 = val;
                coeffx1 = coeffx;
                coeffy1 = coeffy;
            }
            printf("     %f %f   %g       ( %f %f %g )\n", coeffx, coeffy, val, coeffx1, coeffy1, val1);
        }

    printf("COEFF : %f %f\n", coeffx1, coeffy1);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = 1.0*ii/size;
            y = 1.0*jj/size;
            data.image[ID].array.F[jj*size+ii] =  data.image[ID].array.F[jj*size+ii]*(1.0+coeffx1*(x-0.5))*(1.0+coeffy1*(y-0.5));
        }



    vp10 = img_percentile_float(ID_name, 0.1);
    vp90 = img_percentile_float(ID_name, 0.9);
    printf("%f %f\n", vp10, vp90);



    for(ii=0; ii<data.image[ID].md[0].size[0]*data.image[ID].md[0].size[1]; ii++)
        data.image[ID].array.F[ii] = (data.image[ID].array.F[ii]-vp10)/(vp90-vp10);

    for(ii=0; ii<size*size; ii++)
    {
        if(data.image[ID].array.F[ii]>0.05)
            data.image[ID].array.F[ii] = pow(data.image[ID].array.F[ii], 0.5);
        else
            data.image[ID].array.F[ii] = 0.0;
    }

    vp10 = img_percentile_float(ID_name, 0.1);
    vp90 = img_percentile_float(ID_name, 0.9);
    for(ii=0; ii<data.image[ID].md[0].size[0]*data.image[ID].md[0].size[1]; ii++)
        data.image[ID].array.F[ii] = (data.image[ID].array.F[ii]-vp10)/(vp90-vp10);




    rout_min = rout1*0.9;
    rout_max = rout1*1.1;
    xc_min = xc1-2.0;
    xc_max = xc1+2.0;
    yc_min = yc1-2.0;
    yc_max = yc1+2.0;
    pupPA_min = pupPA1-0.1;
    pupPA_max = pupPA1+0.1;
    stretchx_min = stretchx1*0.95;
    stretchx_max = stretchx1*1.05;
    spideroffset_min = spideroffset1*0.9;
    spideroffset_max = spideroffset1*1.1;
    spiderPA_min = spiderPA1-0.1;
    spiderPA_max = spiderPA1+0.1;

    fp = fopen("pupfit.txt", "w");
    fclose(fp);


    rms1 = 1000000000000000000000.0;
    for(iter=0; iter<NBiter; iter++)
    {
        for(spiderPA=spiderPA_min; spiderPA<spiderPA_max+eps; spiderPA += (spiderPA_max-spiderPA_min)/nbstep)
            for(spideroffset=spideroffset_min; spideroffset<spideroffset_max+eps; spideroffset += (spideroffset_max-spideroffset_min)/nbstep)
                for(stretchx=stretchx_min; stretchx<stretchx_max+eps; stretchx += (stretchx_max-stretchx_min)/nbstep)
                    for(rout=rout_min; rout<rout_max+eps; rout += (rout_max-rout_min)/nbstep)
                        for(xc=xc_min; xc<xc_max+eps; xc += (xc_max-xc_min)/nbstep)
                            for(yc=yc_min; yc<yc_max+eps; yc += (yc_max-yc_min)/nbstep)
                                for(pupPA=pupPA_min; pupPA<pupPA_max+eps; pupPA += (pupPA_max-pupPA_min)/nbstep)
                                {
                                    rin = 0.305*rout;
                                    spiderthick = 0.06*rout;
                                    spideroffset = 0.165*rout;
                                    AOsystSim_mkTelPupDM("testpup", size, xc, yc, rin, rout, pupPA, spiderPA, spideroffset, spiderthick, stretchx);
                                    IDt = image_ID("testpup");
                                    // list_image_ID();
                                    // save_fits("testpup", "!testpup.fits");

                                    rms = 0.0;
                                    for(ii=0; ii<size*size; ii++)
                                    {
                                        v1 = data.image[ID].array.F[ii] - data.image[IDt].array.F[ii];
                                        rms += v1*v1;
                                    }

                                    rms = sqrt(rms/size/size);

                                    if(rms<rms1)
                                    {
                                        rms1 = rms;
                                        rout1 = rout;
                                        xc1 = xc;
                                        yc1 = yc;
                                        pupPA1 = pupPA;
                                        stretchx1 = stretchx;
                                        spiderPA1 = spiderPA;
                                        spideroffset1 = spideroffset;
                                    }
                                    delete_image_ID("testpup");
                                    printf("%f %f %f  %f -> rms = %g\n", rout, xc, yc, pupPA, rms);
                                    fp = fopen("pupfit.txt", "a");
                                    fprintf(fp,"%f %f %f %f %g\n", rout, xc, yc, pupPA, rms);
                                    fclose(fp);
                                }

        printf("ITERATION %ld    %f %f  %f  %f  ->  %f\n", iter, rout1, xc1, yc1, pupPA1, rms1);


        spideroffset_range = spideroffset_max-spideroffset_min;
        spideroffset_range /= 3.0;
        spideroffset_min = spideroffset1 - 0.5*spideroffset_range;
        spideroffset_max = spideroffset1 + 0.5*spideroffset_range;

        spiderPA_range = spiderPA_max-spiderPA_min;
        spiderPA_range /= 3.0;
        spiderPA_min = spiderPA1 - 0.5*spiderPA_range;
        spiderPA_max = spiderPA1 + 0.5*spiderPA_range;

        stretchx_range = stretchx_max-stretchx_min;
        stretchx_range /= 3.0;
        stretchx_min = stretchx1 - 0.5*stretchx_range;
        stretchx_max = stretchx1 + 0.5*stretchx_range;

        rout_range = rout_max-rout_min;
        rout_range /= 3.0;
        rout_min = rout1 - 0.5*rout_range;
        rout_max = rout1 + 0.5*rout_range;

        xc_range = xc_max-xc_min;
        xc_range /= 3.0;
        xc_min = xc1 - 0.5*xc_range;
        xc_max = xc1 + 0.5*xc_range;

        yc_range = yc_max-yc_min;
        yc_range /= 3.0;
        yc_min = yc1 - 0.5*yc_range;
        yc_max = yc1 + 0.5*yc_range;

        pupPA_range = xc_max-xc_min;
        pupPA_range /= 3.0;
        pupPA_min = xc1 - 0.5*xc_range;
        pupPA_max = xc1 + 0.5*pupPA_range;
    }



    printf("BEST SOLUTION: \n");
    printf("     xc            %f\n", xc1);
    printf("     yc            %f\n", yc1);
    printf("     rout          %f\n", rout1);
    printf("     pupPA         %f\n", pupPA1);
    printf("     spiderPA      %f\n", spiderPA1);
    printf("     spideroffset  %f\n", spideroffset1);
    printf("     stretchx      %f\n", stretchx1);

    rout = rout1;
    rin = 0.315*rout;
    spiderthick = 0.06*rout;

    AOsystSim_mkTelPupDM(IDtelpup_name, size, xc1, yc1, rin, rout, pupPA1, spiderPA1, spideroffset1, spiderthick, stretchx1);

    return(ID);
}



/** \brief DM control signals to DMshape
 *
 *
 */

int AOsystSim_DMshape(char *IDdmctrl_name, char *IDdmifc_name, char *IDdm_name)
{
    long IDdmctrl, IDdmifc, IDdm;
    long dmsizex, dmsizey, DMnbact;
    long ii, jj;
    long dmact;
    double eps=1.0e-12;
    long k;
    long DMifpixarray_NBpix0;
    

    IDdmctrl = image_ID(IDdmctrl_name);

    IDdmifc = image_ID(IDdmifc_name);
    dmsizex = data.image[IDdmifc].md[0].size[0];
    dmsizey = data.image[IDdmifc].md[0].size[1];
    DMnbact = data.image[IDdmifc].md[0].size[2];

  
    if(DMifpixarray_init==0)
    {
        DMifpixarray_NBpix = 0.0;
        for(dmact=0; dmact<DMnbact; dmact++)
            for(ii=0; ii<dmsizex*dmsizey; ii++)
                if(fabs(data.image[IDdmifc].array.F[dmact*dmsizex*dmsizey+ii])>eps)
                    DMifpixarray_NBpix++;
        
        if(DMifpixarray_val!=NULL)
            {
                printf("WARNING: array DMifpixarray_val is not NULL");
                fflush(stdout);
            }
        if((DMifpixarray_val = (float*) malloc(sizeof(float)*DMifpixarray_NBpix))==NULL)
            {
                printf("ERROR: could not allocate array DMifpixarray_val\n");
                exit(0);
            }
       
        DMifpixarray_index = (long*) malloc(sizeof(long)*DMifpixarray_NBpix);
        DMifpixarray_pixindex = (long*) malloc(sizeof(long)*DMifpixarray_NBpix);
        DMifpixarray_NBpix0 = DMifpixarray_NBpix;
        DMifpixarray_NBpix = 0;
        for(dmact=0; dmact<DMnbact; dmact++)
            for(ii=0; ii<dmsizex*dmsizey; ii++)
                if(fabs(data.image[IDdmifc].array.F[dmact*dmsizex*dmsizey+ii])>eps)
                    {
              //          printf("%ld / %ld \n", DMifpixarray_NBpix, DMifpixarray_NBpix0);
                        DMifpixarray_val[DMifpixarray_NBpix] = data.image[IDdmifc].array.F[dmact*dmsizex*dmsizey+ii];
                        DMifpixarray_index[DMifpixarray_NBpix] = dmact;
                        DMifpixarray_pixindex[DMifpixarray_NBpix] = ii;
                        DMifpixarray_NBpix++;
                    }
    }

   
    IDdm = image_ID(IDdm_name);
    if(IDdm==-1)
        IDdm = create_2Dimage_ID(IDdm_name, dmsizex, dmsizey);
   
 
    for(ii=0; ii<dmsizex*dmsizey; ii++)
        data.image[IDdm].array.F[ii] = 0.0;
/*
    for(dmact=0; dmact<DMnbact; dmact++)
    {
        for(ii=0; ii<dmsizex*dmsizey; ii++)
            data.image[IDdm].array.F[ii] += data.image[IDdmctrl].array.F[dmact]*data.image[IDdmifc].array.F[dmact*dmsizex*dmsizey+jj*dmsizex+ii];
    }
*/

    for(k=0;k<DMifpixarray_NBpix;k++)
        data.image[IDdm].array.F[DMifpixarray_pixindex[k]] += data.image[IDdmctrl].array.F[DMifpixarray_index[k]] * DMifpixarray_val[k];
    
    return (0);
}




int AOsystSim_WFSsim_Pyramid(char *inWFc_name, char *outWFSim_name, double modampl, long modnbpts)
{
    long ID_inWFc, ID_outWFSim;
    long arraysize;
    long *imsize;
    long IDa, IDp;
    long ID_inWFccp;
    long arraysize2;
    long IDpyramp, IDpyrpha;
    double lenssize;
    double pcoeff = 1.4;
    double x, y;
    long ii, jj;
    char pnamea[200];
    char pnamep[200];
    char pfnamea[200];
    char pfnamep[200];

    long PYRMOD_nbpts = 8;
    long pmodpt;
    double PYRMOD_rad = 5.0;
    double xc, yc, PA;
    long ID_outWFSim_tmp;
    
    PYRMOD_nbpts = modnbpts;
    PYRMOD_rad = modampl;
    
    
    
    ID_inWFc = image_ID(inWFc_name);
    arraysize = data.image[ID_inWFc].md[0].size[0];
    arraysize2 = arraysize*arraysize;
    lenssize = 0.4*arraysize;


    for(pmodpt=0; pmodpt<PYRMOD_nbpts; pmodpt++)
    {
        sprintf(pnamea, "pyramp_%03ld", pmodpt);
        sprintf(pnamep, "pyrpha_%03ld", pmodpt);
        IDpyramp = image_ID(pnamea);
        IDpyrpha = image_ID(pnamep);
        if((IDpyramp==-1)||(IDpyrpha==-1))
        {
            imsize = (long*) malloc(sizeof(long)*2);
            imsize[0] = arraysize;
            imsize[1] = arraysize;
            IDpyramp = create_image_ID(pnamea, 2, imsize, FLOAT, 0, 0);
            IDpyrpha = create_image_ID("pyrpha0", 2, imsize, FLOAT, 0, 0);
            free(imsize);

            PA = 2.0*M_PI*pmodpt/PYRMOD_nbpts;
            xc = PYRMOD_rad * cos(PA);
            yc = PYRMOD_rad * sin(PA);

            for(ii=0; ii<arraysize; ii++)
                for(jj=0; jj<arraysize; jj++)
                {
                    x = 1.0*(ii-arraysize/2) - xc;
                    y = 1.0*(jj-arraysize/2) - yc;

                    data.image[IDpyrpha].array.F[jj*arraysize+ii] = pcoeff*(fabs(x)+fabs(y));
                    if((fabs(x)>lenssize)||(fabs(y)>lenssize))
                        data.image[IDpyramp].array.F[jj*arraysize+ii] = 0.0;
                    else
                        data.image[IDpyramp].array.F[jj*arraysize+ii] = 1.0;
                }
            gauss_filter("pyrpha0", pnamep, 1.0, 10);
            delete_image_ID("pyrpha0");

            sprintf(pfnamea, "!pyramp_%03ld.fits", pmodpt);
            sprintf(pfnamep, "!pyrpha_%03ld.fits", pmodpt);

            printf("SAVING: %s -> %s\n", pnamea, pfnamea);
            save_fits(pnamea, pfnamea);
            save_fits(pnamep, pfnamep);
        }
    }


    ID_outWFSim = image_ID(outWFSim_name);
    if(ID_outWFSim==-1)
    {
        imsize = (long*) malloc(sizeof(long)*2);
        imsize[0] = arraysize;
        imsize[1] = arraysize;
        ID_outWFSim = create_image_ID(outWFSim_name, 2, imsize, FLOAT, 1, 0);
        free(imsize);
    }

    ID_outWFSim_tmp = image_ID("outpwfsimtmp");
    if(ID_outWFSim_tmp==-1)
    {
        imsize = (long*) malloc(sizeof(long)*2);
        imsize[0] = arraysize;
        imsize[1] = arraysize;
        ID_outWFSim_tmp = create_image_ID("outpwfsimtmp", 2, imsize, FLOAT, 1, 0);
        free(imsize);
    }

    ID_inWFccp = image_ID("pyrwfcin");
    if(ID_inWFccp==-1)
    {
        imsize = (long*) malloc(sizeof(long)*2);
        imsize[0] = arraysize;
        imsize[1] = arraysize;
        ID_inWFccp = create_image_ID("pyrwfcin", 2, imsize, COMPLEX_FLOAT, 0, 0);
        free(imsize);
    }


    data.image[ID_outWFSim].md[0].write = 1;
    for(ii=0; ii<arraysize2; ii++)
        data.image[ID_outWFSim_tmp].array.F[ii] = 0.0;

    for(pmodpt=0; pmodpt<PYRMOD_nbpts; pmodpt++)
    {
        memcpy(data.image[ID_inWFccp].array.CF, data.image[ID_inWFc].array.CF, sizeof(complex_float)*arraysize*arraysize);



        permut("pyrwfcin");
        do2dfft("pyrwfcin","pyrpsfcin");
        permut("pyrpsfcin");
        mk_amph_from_complex("pyrpsfcin", "pyrpsfa", "pyrpsfp");
        delete_image_ID("pyrpsfcin");

        sprintf(pnamea, "pyramp_%03ld", pmodpt);
        sprintf(pnamep, "pyrpha_%03ld", pmodpt);
        IDpyramp = image_ID(pnamea);
        IDpyrpha = image_ID(pnamep);
        IDa = image_ID("pyrpsfa");
        IDp = image_ID("pyrpsfp");

        for(ii=0; ii<arraysize2; ii++)
        {
            data.image[IDa].array.F[ii] *= data.image[IDpyramp].array.F[ii];
            data.image[IDp].array.F[ii] += data.image[IDpyrpha].array.F[ii];
        }


        mk_complex_from_amph("pyrpsfa", "pyrpsfp", "pyrpsfc");
        delete_image_ID("pyrpsfa");
        delete_image_ID("pyrpsfp");

        permut("pyrpsfc");
        do2dfft("pyrpsfc","pyrwfs_pupc");
        delete_image_ID("pyrpsfc");
        permut("pyrwfs_pupc");
        mk_amph_from_complex("pyrwfs_pupc", "pyrwfs_pupa", "pyrwfs_pupp");

        delete_image_ID("pyrwfs_pupp");
        delete_image_ID("pyrwfs_pupc");

        IDa = image_ID("pyrwfs_pupa");

        for(ii=0; ii<arraysize2; ii++)
            data.image[ID_outWFSim_tmp].array.F[ii] += data.image[IDa].array.F[ii]*data.image[IDa].array.F[ii]/PYRMOD_nbpts;
        delete_image_ID("pyrwfs_pupa");
    }
    memcpy(data.image[ID_outWFSim].array.F, data.image[ID_outWFSim_tmp].array.F, sizeof(float)*arraysize*arraysize);
    data.image[ID_outWFSim].md[0].cnt0++;
    data.image[ID_outWFSim].md[0].write = 0;

    return (0);
}







/** \brief simplified AO system simulator (DM command -> WFS image part)
 *
 * creates a DM map(s) and a WF error input
 * When either DM map or WF error input changes, compute intensity outputs (images)
 *
 * syncmode:
 * 0: sync to turbulence
 * 1: sync to DM
 * 2: sync to both
 * default: use delayus
 * 
 */

int AOsystSim_run(int syncmode, long delayus)
{
    long arraysize = 128;
    long ii, jj;
    double puprad, dmrad;
    double dmifscale = 2.0; // scale magnification between full DM map and DM influence function
    long *imsize;
    long DMsize = 50; // default
    long DMnbact;
    long mx, my;
    double x, y, rx, ry;
    long rxi, ryi;
    double u, t, v00, v01, v10, v11, ii0f, jj0f, rxif, ryif;
    long IDif, IDifc;
    double sig;
    long k;
    long IDdmctrl;
    long IDpupm;
    int elem;
    long IDdm0shape;
    long IDfocmask;
    double r;
    double dftzoomfact = 2.0;
    long IDturb;
    long *dmsizearray;
    long ID;
    long IDout;
    long *IDarray;

    int COROmode = 0; // 1 if coronagraph


    puprad = 0.17*arraysize;
    dmrad = 0.20*arraysize;


    // INITIALIZE DM CONTROL ARRAY IF DOESN'T EXIST
    IDdmctrl = image_ID("aosimdmctrl");
    if(IDdmctrl==-1)
        IDdmctrl = read_sharedmem_image("aosimdmctrl");
    if(IDdmctrl==-1)
    {
        dmsizearray = (long*) malloc(sizeof(long)*2);
        dmsizearray[0] = DMsize;
        dmsizearray[1] = DMsize;
        IDdmctrl = create_image_ID("aosimdmctrl", 2, dmsizearray, FLOAT, 1, 0);
        free(dmsizearray);
        COREMOD_MEMORY_image_set_createsem("aosimdmctrl", 2);
    }
    else
        {
            DMsize = data.image[IDdmctrl].md[0].size[0];
            COREMOD_MEMORY_image_set_createsem("aosimdmctrl", 2);
        }
        
    

    for(k=0; k<DMsize*DMsize; k++)
        data.image[IDdmctrl].array.F[k] = ran1()*2.0e-8;

    dmifscale = 0.5*DMsize;
    

    // MAKE DM INFLUENCE FUNCTIONS

    // DM influence functions stored as a data cube
    DMnbact = DMsize*DMsize;
    imsize = (long*) malloc(sizeof(long)*3);
    imsize[0] = arraysize;
    imsize[1] = arraysize;
    imsize[2] = DMnbact;
    IDif = create_image_ID("dmif0", 2, imsize, FLOAT, 0, 0);
    // construct DM influence function (for 1 actuator)
    // step 1: square
    // actuator size = dmifscale*(arraysize*2.0*dmrad/DMsize) [pix]
    printf("DM size = %ld x %ld actuators\n", DMsize, DMsize);
    printf("actuator pix size = %f pix\n", dmifscale*(2.0*dmrad/DMsize));

    list_image_ID();
    for(ii=0; ii<arraysize; ii++)
        for(jj=0; jj<arraysize; jj++)
        {
            x = (1.0*ii-0.5*arraysize); // [pix]
            y = (1.0*jj-0.5*arraysize); // [pix]
            x /= dmifscale*(2.0*dmrad/DMsize);
            y /= dmifscale*(2.0*dmrad/DMsize);
            if((fabs(x)<0.5)&&(fabs(y)<0.5))
                data.image[IDif].array.F[jj*arraysize+ii] = 1.0e-6;
            else
                data.image[IDif].array.F[jj*arraysize+ii] = 0.0;
        }
    printf("convolve\n");
    fflush(stdout);
    // convolve dmif
    save_fits("dmif0", "!dmif0.fits");
    sig = 0.5*dmifscale*(2.0*dmrad/DMsize);
    printf("gauss filter   %lf %ld\n", sig, (long) (2.0*sig));
    fflush(stdout);
    gauss_filter("dmif0", "dmif", sig, (long) (2.0*sig));
    list_image_ID();
    delete_image_ID("dmif0");
    IDif = image_ID("dmif");


    list_image_ID();
    save_fits("dmif", "!dmif.fits");

    IDifc = create_image_ID("dmifc", 3, imsize, FLOAT, 0, 0);
    printf("\n");
 
 
    list_image_ID();
    printf("dmifc = %ld   %ld %ld %ld    %ld %ld\n", IDifc, data.image[IDifc].md[0].size[0], data.image[IDifc].md[0].size[1], data.image[IDifc].md[0].size[2], arraysize, arraysize);

 
    for(mx=0; mx<DMsize; mx++)
        for(my=0; my<DMsize; my++)
        {
            printf("\r actuator %2ld %2ld    ", mx, my);
            fflush(stdout);
            // actuator center coordinates
            // center: mx = DMsize/2-0.5  4 -> 1.5
            ii0f = 0.5*arraysize + 2.0*(mx-DMsize/2)/DMsize*dmrad;
            jj0f = 0.5*arraysize + 2.0*(my-DMsize/2)/DMsize*dmrad;
            for(ii=0; ii<arraysize; ii++)
                for(jj=0; jj<arraysize; jj++)
                {
                    rx = 1.0*ii-ii0f;
                    ry = 1.0*jj-jj0f;
                    rx *= dmifscale;
                    ry *= dmifscale;
                    rxif = rx+arraysize/2;
                    ryif = ry+arraysize/2;
                    rxi = (long) (rxif);
                    ryi = (long) (ryif);
                    u = rxif - rxi;
                    t = ryif - ryi;

                    if((rxi>0)&&(rxi<arraysize-1)&&(ryi>0)&&(ryi<arraysize-1))
                    {
                        v00 = data.image[IDif].array.F[ryi*arraysize+rxi];
                        v01 = data.image[IDif].array.F[ryi*arraysize+rxi+1];
                        v10 = data.image[IDif].array.F[(ryi+1)*arraysize+rxi];
                        v11 = data.image[IDif].array.F[(ryi+1)*arraysize+rxi+1];
                        data.image[IDifc].array.F[(my*DMsize+mx)*arraysize*arraysize + jj*arraysize + ii] = (1.0-u)*(1.0-t)*v00 + (1.0-u)*t*v10 + u*(1.0-t)*v01 + u*t*v11;
                    }
                }
        }
    free(imsize);
    save_fits("dmifc","!dmifc.fits");
    printf("\n");






    // INITIALIZE TURBULENCE SCREEN
        
    imsize = (long*) malloc(sizeof(long)*2);
    imsize[0] = arraysize;
    imsize[1] = arraysize;
    IDturb = create_image_ID("WFturb", 2, imsize, FLOAT, 1, 0);
    free(imsize);
    COREMOD_MEMORY_image_set_createsem("WFturb", 2);
    list_image_ID();

    AOsystSim_DMshape("aosimdmctrl", "dmifc", "dmdisp");
    IDdm0shape = image_ID("dmdisp");
    save_fits("dmdisp", "!dmdisp.fits");





     // INITIALIZE OPTICAL SYSTEM

    optsystsim = (OPTSYST*) malloc(sizeof(OPTSYST)*1);
    optsystsim[0].nblambda = 1;
    optsystsim[0].lambdaarray[0] = 1.6e-6;
    optsystsim[0].beamrad = 0.008; // 8mm
    optsystsim[0].size = arraysize;
    optsystsim[0].pixscale = optsystsim[0].beamrad/50.0;
    optsystsim[0].DFTgridpad = 0;


    optsystsim[0].NB_asphsurfm = 2;
    optsystsim[0].NB_asphsurfr = 0;
    optsystsim[0].NBelem = 100; // to be updated later

    // 0: INPUT PUPIL
    IDpupm = make_disk("pupmask", arraysize, arraysize, 0.5*arraysize, 0.5*arraysize, puprad);
    elem = 0;
    optsystsim[0].elemtype[elem] = 1; // pupil mask
    optsystsim[0].elemarrayindex[elem] = IDpupm;
    optsystsim[0].elemZpos[elem] = 0.0;
    elem++;



    // 1: Turbulence screen
    optsystsim[0].elemtype[elem] = 3; // reflective surface
    optsystsim[0].elemarrayindex[elem] = 0; // index
    optsystsim[0].ASPHSURFMarray[0].surfID = IDturb;
    optsystsim[0].elemZpos[elem] = 0.0;
    elem++;

    // 2: DM 0
    optsystsim[0].elemtype[elem] = 3; // reflective surface
    optsystsim[0].elemarrayindex[elem] = 1; // index
    optsystsim[0].ASPHSURFMarray[1].surfID = IDdm0shape;
    optsystsim[0].elemZpos[elem] = 0.0;
    optsystsim[0].keepMem[elem] = 1;
    elem++;



    if(COROmode==1)    // FOCAL PLANE MASK
    {
        IDfocmask = create_2DCimage_ID("focpm", arraysize, arraysize);
        for(ii=0; ii<arraysize; ii++)
            for(jj=0; jj<arraysize; jj++)
            {
                x = 1.0*ii-0.5*arraysize;
                y = 1.0*jj-0.5*arraysize;
                r = sqrt(x*x+y*y);
                if(r<20.0*dftzoomfact)
                {
                    data.image[IDfocmask].array.CF[jj*arraysize+ii].re = 1.0;  // 1-(CA) : 1.0=opaque 0=transmissive 2.0=phase shifting
                    data.image[IDfocmask].array.CF[jj*arraysize+ii].im = 0.0;
                }
                else
                {
                    data.image[IDfocmask].array.CF[jj*arraysize+ii].re = 0.0;
                    data.image[IDfocmask].array.CF[jj*arraysize+ii].im = 0.0;
                }
            }

        optsystsim[0].elemtype[elem] = 5; // focal plane mask
        optsystsim[0].FOCMASKarray[0].fpmID = IDfocmask;
        optsystsim[0].FOCMASKarray[0].zfactor = dftzoomfact;
        optsystsim[0].FOCMASKarray[0].mode = 1;
        optsystsim[0].elemZpos[elem] = optsystsim[0].elemZpos[elem-1]; // plane from which FT is done
        elem++;
    }

    optsystsim[0].NBelem = elem;

    optsystsim[0].SAVE = 1;

    // propagate
    OptSystProp_run(optsystsim, 0, 0, optsystsim[0].NBelem, "./testconf/");

    ID = image_ID("psfi0");
    imsize = (long*) malloc(sizeof(long)*2);
    imsize[0] = data.image[ID].md[0].size[0];
    imsize[1] = data.image[ID].md[0].size[1];
    imsize[2] = data.image[ID].md[0].size[2];
    IDout = create_image_ID("aosimpsfout", 3, imsize, FLOAT, 1, 0);
    free(imsize);

    COREMOD_MEMORY_image_set_createsem("aosimpsfout", 2);
    data.image[IDout].md[0].write = 1;
    memcpy(data.image[IDout].array.F, data.image[ID].array.F, sizeof(FLOAT)*data.image[ID].md[0].size[0]*data.image[ID].md[0].size[1]*data.image[ID].md[0].size[2]);
    data.image[IDout].md[0].cnt0++;
    data.image[IDout].md[0].write = 0;
    COREMOD_MEMORY_image_set_sempost("aosimpsfout", -1);

    IDarray = (long*) malloc(sizeof(long)*2);
    IDarray[0] = image_ID("WFturb");
    IDarray[1] = image_ID("aosimdmctrl");



    while(1)
    {
        AOsystSim_DMshape("aosimdmctrl", "dmifc", "dmdisp");

        OptSystProp_run(optsystsim, 0, 0, optsystsim[0].NBelem, "./testconf/");

        mk_complex_from_amph("WFamp0_002", "WFpha0_002", "wfc");

        AOsystSim_WFSsim_Pyramid("wfc", "aosimwfsim", 0.0, 1);
       // COREMOD_MEMORY_image_set_sempost("aosimwfsim", 0);
        delete_image_ID("wfc");

        ID = image_ID("psfi0");
        data.image[IDout].md[0].write = 1;
        memcpy(data.image[IDout].array.F, data.image[ID].array.F, sizeof(FLOAT)*data.image[ID].md[0].size[0]*data.image[ID].md[0].size[1]*data.image[ID].md[0].size[2]);
        data.image[IDout].md[0].cnt0++;
        data.image[IDout].md[0].write = 0;
        COREMOD_MEMORY_image_set_sempost("aosimpsfout", -1);

        switch (syncmode) {
            case 0 : // sync to turbulence
            waitforsemID(IDarray[0]);
            break;
            case 1 : // sync to DM
            waitforsemID(IDarray[1]);
            break;
            case 2 :
            COREMOD_MEMORY_image_set_semwait_OR_IDarray(IDarray, 2);
            break;
            default :
            usleep(delayus);
            break;
        }
        
        COREMOD_MEMORY_image_set_semflush_IDarray(IDarray, 2);
    }


    return(0);
}












int AOsystSim_extremeAO_contrast_sim()
{
    EXAOSIMCONF *exaosimconf;
    long CN2layer;
    double tmpv1, tmpv2, tmpv3, tmpv4;
    double tmpC;

    double lambda_V = 0.545e-6;
    double zeropt_V = 9.9690e10;

    double lambda_R = 0.638e-6;
    double zeropt_R = 7.2384e10;

    double lambda_I = 0.797e-6;
    double zeropt_I = 4.5825e10;

    double lambda_J = 1.22e-6;
    double zeropt_J = 1.9422e10;

    double lambda_H = 1.63e-6;
    double zeropt_H = 9.4440e9;

    double lambda_K = 2.19e-6;
    double zeropt_K = 4.3829e9;

    double lambda_L = 3.45e-6;
    double zeropt_L = 1.2292e9;

    double zeroptWFS;
    double zeroptWFSsci;

    double sourcemag_wfs;
    double sourcemag_sci;
    double lambdaBwfs = 0.4; // dlambda/lambda
    double lambdaBsci = 0.4; // dlambda/lambda
   double systemEfficiency = 0.3;

    double coeff;
    FILE *fp;
    double tmpA, tmpB;
    double dsa;

    double nwfs,nsci; // refractive indices

    double tobs = 3600.0; // observation time
    double crosstime;

    double WFStlim = 0.0000; // WFS min exposure time
    double sciWFStlim = 0.0000; // sci WFS min exposure time
    double IWAld = 1.3;
    int OK;
    double att2;


    exaosimconf = (EXAOSIMCONF*) malloc(sizeof(EXAOSIMCONF));

    // initialization
    exaosimconf[0].lambda0 = 0.55e-6;
    exaosimconf[0].lambdai = 1.6e-6;
    exaosimconf[0].lambdawfs = 0.8e-6;
    exaosimconf[0].D = 30.0;
    exaosimconf[0].r0 = 0.15;
    exaosimconf[0].windspeed = 10.0;
    exaosimconf[0].betapWFS = sqrt(2.0);
    exaosimconf[0].betaaWFS = sqrt(2.0);
    exaosimconf[0].betapWFSsci = 2.0;
    exaosimconf[0].betaaWFSsci = 2.0;
    exaosimconf[0].framedelay = 1.5;
    nwfs = OPTICSMATERIALS_n( OPTICSMATERIALS_code("Air"), exaosimconf[0].lambdawfs);
    nsci = OPTICSMATERIALS_n( OPTICSMATERIALS_code("Air"), exaosimconf[0].lambdai);

    printf("n = %f %f\n", nwfs, nsci);

    for(CN2layer=0; CN2layer<20; CN2layer++)
    {
        exaosimconf[0].CN2layer_h[CN2layer] = 0.0;
        exaosimconf[0].CN2layer_coeff[CN2layer] = 0.0;
    }
    exaosimconf[0].CN2layer_h[0] = 500.0;
    exaosimconf[0].CN2layer_coeff[0] = 0.2283;
    exaosimconf[0].CN2layer_h[1] = 1000.0;
    exaosimconf[0].CN2layer_coeff[1] = 0.0883;
    exaosimconf[0].CN2layer_h[2] = 2000.0;
    exaosimconf[0].CN2layer_coeff[2] = 0.0666;
    exaosimconf[0].CN2layer_h[3] = 4000.0;
    exaosimconf[0].CN2layer_coeff[3] = 0.1458;
    exaosimconf[0].CN2layer_h[4] = 8000.0;
    exaosimconf[0].CN2layer_coeff[4] = 0.3350;
    exaosimconf[0].CN2layer_h[5] = 16000.0;
    exaosimconf[0].CN2layer_coeff[5] = 0.1350;


    zeroptWFS = zeropt_I;
    exaosimconf[0].lambdawfs = lambda_I;
  
    zeroptWFSsci = zeropt_H;
    exaosimconf[0].lambdai = lambda_H;

    sourcemag_wfs = 8.0;
    sourcemag_sci = 6.0;
    exaosimconf[0].Fwfs = zeroptWFS*1.0e6*(exaosimconf[0].lambdawfs*lambdaBwfs)*pow(100.0, -0.2*sourcemag_wfs)*lambdaBwfs*systemEfficiency;
    exaosimconf[0].Fsci = zeroptWFSsci*1.0e6*(exaosimconf[0].lambdai*lambdaBsci)*pow(100.0, -0.2*sourcemag_sci)*lambdaBsci*systemEfficiency;

    exaosimconf[0].alpha_arcsec=0.25;

    fp = fopen("result.out.txt", "w");

    for(exaosimconf[0].alpha_arcsec=0.010; exaosimconf[0].alpha_arcsec<0.25; exaosimconf[0].alpha_arcsec+=0.001)
    {
        OK = 0;
        while(OK==0)
        {
            exaosimconf[0].alpha = exaosimconf[0].alpha_arcsec/3600/180*M_PI;
            exaosimconf[0].alpha_ld = exaosimconf[0].alpha / (exaosimconf[0].lambdai / exaosimconf[0].D);
            exaosimconf[0].f = exaosimconf[0].alpha/exaosimconf[0].lambdai;
            if(exaosimconf[0].alpha_ld>IWAld)
                OK = 1;
            else
                exaosimconf[0].alpha_arcsec+=0.001;
        }

        if(0) // SHWFS
        {
            dsa=0.15;
            exaosimconf[0].betapWFS = 1.48/(exaosimconf[0].f*dsa)*sqrt(1.0+(dsa*dsa/exaosimconf[0].r0/exaosimconf[0].r0));
        }

        exaosimconf[0].f_wfs = exaosimconf[0].alpha/exaosimconf[0].lambdawfs;
        exaosimconf[0].f_0 = exaosimconf[0].alpha/exaosimconf[0].lambda0;

        exaosimconf[0].hf = 0.22 * exaosimconf[0].lambda0 / (pow(exaosimconf[0].f, 11.0/6.0) * exaosimconf[0].D * pow(exaosimconf[0].r0, 5.0/6.0));

        exaosimconf[0].X = 0.0;
        exaosimconf[0].dX = 0.0;
        exaosimconf[0].dY = 0.0;
        for(CN2layer=0; CN2layer<20; CN2layer++)
        {
            tmpv1 = cos( M_PI * exaosimconf[0].CN2layer_h[CN2layer] * exaosimconf[0].f * exaosimconf[0].f * exaosimconf[0].lambdai);
            tmpv2 = cos( M_PI * exaosimconf[0].CN2layer_h[CN2layer] * exaosimconf[0].f * exaosimconf[0].f * exaosimconf[0].lambdawfs);
            tmpv3 = sin( M_PI * exaosimconf[0].CN2layer_h[CN2layer] * exaosimconf[0].f * exaosimconf[0].f * exaosimconf[0].lambdai);
            tmpv4 = sin( M_PI * exaosimconf[0].CN2layer_h[CN2layer] * exaosimconf[0].f * exaosimconf[0].f * exaosimconf[0].lambdawfs);
            exaosimconf[0].X += exaosimconf[0].CN2layer_coeff[CN2layer] * tmpv1 * tmpv1;
            exaosimconf[0].dX += exaosimconf[0].CN2layer_coeff[CN2layer] * (tmpv1-tmpv2) * (tmpv1-tmpv2);
            exaosimconf[0].dY += exaosimconf[0].CN2layer_coeff[CN2layer] * (tmpv3-tmpv4) * (tmpv3-tmpv4);
        }
        exaosimconf[0].Y = sqrt(1.0 - exaosimconf[0].X*exaosimconf[0].X);

        exaosimconf[0].C0 = pow(M_PI*exaosimconf[0].hf/exaosimconf[0].lambdai, 2.0) * exaosimconf[0].X;
        exaosimconf[0].C1 = pow(M_PI*exaosimconf[0].hf/exaosimconf[0].lambdai, 2.0) * exaosimconf[0].Y;

        tmpA = 2.0*M_PI*exaosimconf[0].hf*exaosimconf[0].windspeed*exaosimconf[0].f*exaosimconf[0].framedelay;
        tmpB = exaosimconf[0].lambdawfs/M_PI * exaosimconf[0].betapWFS / sqrt(exaosimconf[0].Fwfs * M_PI) / exaosimconf[0].D;
        exaosimconf[0].twfs_opt = pow(0.5*tmpB/tmpA*tmpB/tmpA, 1.0/3.0);
        exaosimconf[0].twfs = exaosimconf[0].twfs_opt;

        if(exaosimconf[0].twfs<WFStlim)
            exaosimconf[0].twfs = WFStlim;

        exaosimconf[0].hfca = tmpA * exaosimconf[0].twfs;   // lag
        exaosimconf[0].hfcb = tmpB / sqrt(exaosimconf[0].twfs);  // noise

        exaosimconf[0].hfc = sqrt(exaosimconf[0].hfca*exaosimconf[0].hfca + exaosimconf[0].hfcb*exaosimconf[0].hfcb);

        exaosimconf[0].C2 = pow(M_PI*exaosimconf[0].hfc/exaosimconf[0].lambdai, 2.0);
        exaosimconf[0].C2_wfs = exaosimconf[0].C2 * pow( exaosimconf[0].lambdai/exaosimconf[0].lambdawfs,2.0);

        exaosimconf[0].twfs_opt_amp = exaosimconf[0].twfs_opt*pow(exaosimconf[0].X/exaosimconf[0].Y, 1.0/3.0)*pow(exaosimconf[0].betaaWFS/exaosimconf[0].betapWFS,2.0/3.0);
        exaosimconf[0].C3 = exaosimconf[0].C2*pow(exaosimconf[0].Y/exaosimconf[0].X, 1.0/3.0)*pow(exaosimconf[0].betaaWFS/exaosimconf[0].betapWFS,4.0/3.0);

        exaosimconf[0].C4 = exaosimconf[0].C0*exaosimconf[0].dX;
        exaosimconf[0].C5 = exaosimconf[0].C1*exaosimconf[0].dY;

        exaosimconf[0].C6 = exaosimconf[0].C0 * pow((nwfs-nsci)/(1.0-0.5*(nwfs+nsci)), 2.0);



        printf("alpha_ld = %.2lf l/D\n", exaosimconf[0].alpha_ld);
        printf("f = %f  -> p = %g m\n", exaosimconf[0].f, 1.0/exaosimconf[0].f);
        printf("X = %g\n", exaosimconf[0].X);
        printf("Y = %g\n", exaosimconf[0].Y);
        printf("dX = %g\n", exaosimconf[0].dX);
        printf("dY = %g\n", exaosimconf[0].dY);
        printf("OPTIMAL WFS exposure time = %g sec  -> %.3lf kHz\n", exaosimconf[0].twfs_opt, 0.001/exaosimconf[0].twfs_opt);
        printf("single frequ WF error: %g -> %g  (%g + %g)\n", exaosimconf[0].hf, exaosimconf[0].hfc, exaosimconf[0].hfca, exaosimconf[0].hfcb);
        printf("WFS total flux per frame = %lf ph\n", exaosimconf[0].twfs*exaosimconf[0].Fwfs*exaosimconf[0].D*exaosimconf[0].D/4.0);
        crosstime = exaosimconf[0].D/exaosimconf[0].windspeed;
        printf("C0 contrast = %20g\n", exaosimconf[0].C0);
        printf("C1 contrast = %20g\n", exaosimconf[0].C1);
        printf("C2 contrast = %20g \n", exaosimconf[0].C2);
        printf("   WFS  lag = %20g    %20g\n", pow(M_PI*exaosimconf[0].hfca/exaosimconf[0].lambdai, 2.0), pow(M_PI*exaosimconf[0].hfca/exaosimconf[0].lambdai, 2.0)/sqrt(tobs/crosstime)*2.0);
        printf("   WFS noise= %20g    %20g\n", pow(M_PI*exaosimconf[0].hfcb/exaosimconf[0].lambdai, 2.0), pow(M_PI*exaosimconf[0].hfcb/exaosimconf[0].lambdai, 2.0)/sqrt(tobs/exaosimconf[0].twfs)*2.0);
        printf("C3 contrast = %20g    %20g\n", exaosimconf[0].C3, exaosimconf[0].C3/sqrt(tobs/crosstime)*2.0);
        printf("C4 contrast = %20g    %20g\n", exaosimconf[0].C4, exaosimconf[0].C4/sqrt(tobs/crosstime)*2.0);
        printf("C5 contrast = %20g    %20g\n", exaosimconf[0].C5, exaosimconf[0].C5/sqrt(tobs/crosstime)*2.0);
        printf("C6 contrast = %20g    %20g\n", exaosimconf[0].C6, exaosimconf[0].C6/sqrt(tobs/crosstime)*2.0);
        exaosimconf[0].Csum = exaosimconf[0].C2+exaosimconf[0].C3+exaosimconf[0].C4+exaosimconf[0].C5+exaosimconf[0].C6;
        exaosimconf[0].Csum_detection = pow(M_PI*exaosimconf[0].hfca/exaosimconf[0].lambdai, 2.0)/sqrt(tobs/crosstime)*2.0 + pow(M_PI*exaosimconf[0].hfcb/exaosimconf[0].lambdai, 2.0)/sqrt(tobs/exaosimconf[0].twfs)*2.0 + exaosimconf[0].C3/sqrt(tobs/crosstime)*2.0 + exaosimconf[0].C4/sqrt(tobs/crosstime)*2.0 + exaosimconf[0].C5/sqrt(tobs/crosstime)*2.0+exaosimconf[0].C6/sqrt(tobs/crosstime)*2.0;
        printf("TOTAL CONTRAST = %20g   %20g\n", exaosimconf[0].Csum, exaosimconf[0].Csum_detection);
        printf("WFS speckle  contrast = %g    ph per WFS speckle per frame = %g\n", exaosimconf[0].C2_wfs, exaosimconf[0].C2_wfs*exaosimconf[0].twfs_opt*exaosimconf[0].Fwfs*exaosimconf[0].D*exaosimconf[0].D/4.0);
        printf("Time lag speckle lifetime = %.4f sec (intensity), %.4f sec (complex amplitude)\n", exaosimconf[0].D/exaosimconf[0].windspeed, (1.0/exaosimconf[0].f)/exaosimconf[0].windspeed/2.0/M_PI);


        // NEAR-IR LOOP
        

        // time lag attenuation
        tmpA = 2.0*M_PI*exaosimconf[0].hfca*exaosimconf[0].windspeed*exaosimconf[0].f*exaosimconf[0].framedelay;
        tmpB = exaosimconf[0].lambdai/M_PI * exaosimconf[0].betapWFSsci / sqrt(exaosimconf[0].Fsci * M_PI) / exaosimconf[0].D;
        exaosimconf[0].twfssci_opt = pow(0.5*tmpB/tmpA*tmpB/tmpA, 1.0/3.0);
        exaosimconf[0].twfssci = exaosimconf[0].twfssci_opt;
        if(exaosimconf[0].twfssci<sciWFStlim)
            exaosimconf[0].twfssci = sciWFStlim;
        exaosimconf[0].TL_hfca = tmpA * exaosimconf[0].twfssci;   // lag
        exaosimconf[0].TL_hfcb = tmpB / sqrt(exaosimconf[0].twfssci);  // noise
        exaosimconf[0].TL_hfc = sqrt(exaosimconf[0].TL_hfca*exaosimconf[0].TL_hfca + exaosimconf[0].TL_hfcb*exaosimconf[0].TL_hfcb);
        exaosimconf[0].C7 = pow(M_PI*exaosimconf[0].TL_hfc/exaosimconf[0].lambdai, 2.0);
  
        exaosimconf[0].C8 = (exaosimconf[0].C7+pow(M_PI*exaosimconf[0].hfcb/exaosimconf[0].lambdai, 2.0))*pow(exaosimconf[0].Y/exaosimconf[0].X, 1.0/3.0)*pow(exaosimconf[0].betaaWFS/exaosimconf[0].betapWFS,4.0/3.0);
        att2 = (exaosimconf[0].C7+pow(M_PI*exaosimconf[0].hfcb/exaosimconf[0].lambdai, 2.0))/exaosimconf[0].C2;
  
  
        exaosimconf[0].C9 = att2*exaosimconf[0].C4;
        exaosimconf[0].C10 = att2*exaosimconf[0].C5;
 
 
        // refractive index chromaticity attenuation
        tmpA = 2.0*M_PI*exaosimconf[0].hf*(nwfs-nsci)/(1.0-0.5*(nwfs+nsci))*exaosimconf[0].windspeed*exaosimconf[0].f*exaosimconf[0].framedelay;
        tmpB = exaosimconf[0].lambdai/M_PI * exaosimconf[0].betapWFSsci / sqrt(exaosimconf[0].Fsci * M_PI) / exaosimconf[0].D;
        exaosimconf[0].twfssci_opt = pow(0.5*tmpB/tmpA*tmpB/tmpA, 1.0/3.0);
        exaosimconf[0].twfssci = exaosimconf[0].twfssci_opt;
        if(exaosimconf[0].twfssci<sciWFStlim)
            exaosimconf[0].twfssci = sciWFStlim;
        exaosimconf[0].RIC_hfca = tmpA * exaosimconf[0].twfssci;   // lag
        exaosimconf[0].RIC_hfcb = tmpB / sqrt(exaosimconf[0].twfssci);  // noise
        exaosimconf[0].RIC_hfc = sqrt(exaosimconf[0].RIC_hfca*exaosimconf[0].RIC_hfca + exaosimconf[0].RIC_hfcb*exaosimconf[0].RIC_hfcb);
        exaosimconf[0].C11 = pow(M_PI*exaosimconf[0].RIC_hfc/exaosimconf[0].lambdai, 2.0);
        

        exaosimconf[0].Csum2 = 0.0;
        exaosimconf[0].Csum2 += exaosimconf[0].C7;
        exaosimconf[0].Csum2 += pow(M_PI*exaosimconf[0].hfcb/exaosimconf[0].lambdai, 2.0);
        exaosimconf[0].Csum2 += exaosimconf[0].C8;
        exaosimconf[0].Csum2 += exaosimconf[0].C9;
        
        exaosimconf[0].Csum2ave = 0.0;
        exaosimconf[0].Csum2ave += exaosimconf[0].C7/sqrt(tobs/exaosimconf[0].twfssci)*2.0;
        exaosimconf[0].Csum2ave += pow(M_PI*exaosimconf[0].hfcb/exaosimconf[0].lambdai, 2.0)/sqrt(tobs/exaosimconf[0].twfs)*2.0;
        exaosimconf[0].Csum2ave += exaosimconf[0].C8/sqrt(tobs/exaosimconf[0].twfssci)*2.0;
        exaosimconf[0].Csum2ave += exaosimconf[0].C9/sqrt(tobs/exaosimconf[0].twfssci)*2.0;
        // #1 : arcsec
        // #2 : C0
        // #3 : C1
        // #4 : C2
        // #5 : C3
        // #6 : C4
        // #7 : C5
        // #8 : C6
        // #9 : Csum
        // #10 : Csum_detection [5 sig]
        // #11 : wfs etime
        // #12 : ph/frame

        // #13 : C7  time lag correction residual (C2 ->)
        // #14 : C8 (C3->)
        // #15 : C9 (C4->)
        // #16 : C10 (C5->)
        // #17 : C11 (C6->)
        // #18 : RAW CONTRAST
        // #19 : Detection limit [5 sig]
        // #20 : Photon noise limit 1hr [5 sig]
        // #21 : wfs etime
        // #22 : near-IR ph/speckle/frame

        // #23 : photon noise limit loop 1
 
        
        fprintf(fp, "%f %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g %15g\n", exaosimconf[0].alpha_arcsec, log10(exaosimconf[0].C0), log10(exaosimconf[0].C1), log10(exaosimconf[0].C2), log10(exaosimconf[0].C3), log10(exaosimconf[0].C4), log10(exaosimconf[0].C5), log10(exaosimconf[0].C6), log10(exaosimconf[0].Csum), log10(exaosimconf[0].Csum_detection*5), exaosimconf[0].twfs, exaosimconf[0].twfs*exaosimconf[0].Fwfs*exaosimconf[0].D*exaosimconf[0].D/4.0, log10(exaosimconf[0].C7), log10(exaosimconf[0].C8), log10(exaosimconf[0].C9), log10(exaosimconf[0].C10), log10(exaosimconf[0].C11), log10(exaosimconf[0].Csum2), log10(exaosimconf[0].Csum2ave*5), log10(5.0*exaosimconf[0].Csum2/sqrt(tobs*exaosimconf[0].Fsci*exaosimconf[0].D*exaosimconf[0].D/4.0*exaosimconf[0].Csum2)), exaosimconf[0].twfssci, exaosimconf[0].twfssci*exaosimconf[0].Fsci*exaosimconf[0].D*exaosimconf[0].D/4.0*exaosimconf[0].Csum2, log10(5.0*exaosimconf[0].Csum/sqrt(tobs*exaosimconf[0].Fsci*exaosimconf[0].D*exaosimconf[0].D/4.0*exaosimconf[0].Csum)));
        printf("Nphoton Sci = %g\n",tobs*exaosimconf[0].Fsci*exaosimconf[0].D*exaosimconf[0].D/4.0*exaosimconf[0].Csum2);
        //pow(M_PI*exaosimconf[0].hfca/exaosimconf[0].lambdai, 2.0), pow(M_PI*exaosimconf[0].hfcb/exaosimconf[0].lambdai, 2.0));
    }
    fclose(fp);

    fp = fopen("printcmd", "w");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fclose(fp);


    free(exaosimconf);

    return(0);
}





