#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>


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



long NBprobesG = 3;
int CENTERprobe=1; // 1 if center probe included
long NBoptVar;
long double probe_re[100]; // 100 probes maximum
long double probe_im[100];
long double nprobe_re[100]; // 100 probes maximum
long double nprobe_im[100];


long double probe_tflux[100]; // test point computed flux
long double Cflux = 1.0; // coherent flux, ph for CA unity circle
long double probe_nmflux[100]; // measured flux (noisy)
long double probe_nmnoise[100]; // measured flux (noisy)


double tmpvalue1, tmpvalue2;



OPTSYST *optsystsim;


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string (not image)
// 4: existing image
// 5: string


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
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,2)==0)
    {
        AOsystSim_run(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numl);
        return 0;
    }
    else
        return 1;
}

int AOsystSim_FPWFS_mkprobes_CLI()
{
    if(CLI_checkarg(1,3)+CLI_checkarg(2,3)+CLI_checkarg(3,2)+CLI_checkarg(4,2)+CLI_checkarg(5,1)+CLI_checkarg(6,1)+CLI_checkarg(7,1)+CLI_checkarg(8,1)+CLI_checkarg(9,2)==0)
    {
        AOsystSim_FPWFS_mkprobes(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl, data.cmdargtoken[5].val.numf, data.cmdargtoken[6].val.numf, data.cmdargtoken[7].val.numf, data.cmdargtoken[8].val.numf, data.cmdargtoken[9].val.numl);
        return 0;
    }
    else
        return 1;
}


int AOsystSim_FPWFS_sensitivityAnalysis_cli()
{
     if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,2)+CLI_checkarg(4,2)==0)
    {   
        AOsystSim_FPWFS_sensitivityAnalysis(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl);
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
    strcpy(data.cmd[data.NBcmd].syntax,"<syncmode> <DMindex> <delayus>");
    strcpy(data.cmd[data.NBcmd].example,"AOsystsim");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOsystSim_run(int syncmode, long DMindex, long delayus)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"AOsystexaosim");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOsystSim_extremeAO_contrast_sim;
    strcpy(data.cmd[data.NBcmd].info,"run extremeAO analysis");
    strcpy(data.cmd[data.NBcmd].syntax,"no argument");
    strcpy(data.cmd[data.NBcmd].example,"AOsystexaosim");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOsystSim_extremeAO_contrast_sim()");
    data.NBcmd++;
    
    strcpy(data.cmd[data.NBcmd].key,"AOsystmkABprobes");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOsystSim_FPWFS_mkprobes_CLI;
    strcpy(data.cmd[data.NBcmd].info,"make AB probes for focal plane sensing and coherence measurement");
    strcpy(data.cmd[data.NBcmd].syntax, "AOsystmkABprobes prA prB 50 50 0.7 0.1 0.7 0.1 1");
    strcpy(data.cmd[data.NBcmd].example,"AOsystmkABprobes");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOsystSim_FPWFS_mkprobes(char *IDprobeA_name, char *IDprobeB_name, long dmxsize, long dmysize, double CPAmax, double CPArmin, double CPArmax, double RMSampl, long modegeom)");
    data.NBcmd++;
 
    strcpy(data.cmd[data.NBcmd].key,"AOsystFPWFSan");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOsystSim_FPWFS_sensitivityAnalysis_cli;
    strcpy(data.cmd[data.NBcmd].info,"run focal plane WFS sensitivity analysis");
    strcpy(data.cmd[data.NBcmd].syntax,"<mapmode> <mode> <optmode> <NBprobes>");
    strcpy(data.cmd[data.NBcmd].example,"AOsystFPWFSan");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOsystSim_FPWFS_sensitivityAnalysis(int mapmode, int mode, int optmode, int NBprobes)");
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
        COREMOD_MEMORY_image_set_createsem(outWFSim_name, 5);
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
        mk_amph_from_complex("pyrpsfcin", "pyrpsfa", "pyrpsfp", 0);
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


        mk_complex_from_amph("pyrpsfa", "pyrpsfp", "pyrpsfc", 0);
        delete_image_ID("pyrpsfa");
        delete_image_ID("pyrpsfp");

        permut("pyrpsfc");
        do2dfft("pyrpsfc","pyrwfs_pupc");
        delete_image_ID("pyrpsfc");
        permut("pyrwfs_pupc");
        mk_amph_from_complex("pyrwfs_pupc", "pyrwfs_pupa", "pyrwfs_pupp", 0);

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




//AOsystSim_runWFS(2, "aosimwfsim");

int AOsystSim_runWFS(long index, char *IDout_name)
{
    long cnt0;
    long IDinamp;
    char imnameamp[200];
    char imnamepha[200];
    int ret;
    
    ret = sprintf(imnameamp, "WFamp0_%03ld", index);
    ret = sprintf(imnamepha, "WFpha0_%03ld", index);
    IDinamp = image_ID(imnameamp);

    cnt0 = 0;
    
    while(1)
    {
        while(cnt0 == data.image[IDinamp].md[0].cnt0)
            usleep(50);
        cnt0 = data.image[IDinamp].md[0].cnt0;
        
        mk_complex_from_amph(imnameamp, imnamepha, "_tmpwfc", 0);
        AOsystSim_WFSsim_Pyramid("_tmpwfc", IDout_name, 0.0, 1);
        
        delete_image_ID("_tmpwfc");
    }
    
    
    return(0);
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

int AOsystSim_run(int syncmode, long DMindex, long delayus)
{
    long arraysize = 128;
    long ii, jj, ii1, jj1;
    double puprad, dmrad;
    double dmifscale = 2.0; // scale magnification between full DM map and DM influence function
    long *imsize;
    long DMsize = 50; // default
    long DMnbact;
    long mx, my;
    double x, y, rx, ry, ii1ld, rld;
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
    long iter;
    
    long *dhsizearray;
    long IDdh, IDdhmask;
    long dhxsize, dhysize, dhsize;
    long dhxoffset, dhyoffset;
    long IDre, IDim;
    char imdhname[200];
    char name[200];
    int COROmode = 0; // 1 if coronagraph
    
    float pupradcoeff = 0.17;
    float dmradcoeff = 0.20;
    float iwald = 2.0;

    char imnameamp[200];
    char imnamepha[200];
    int ret;
    long index;

    puprad = pupradcoeff*arraysize;
    dmrad = dmradcoeff*arraysize;


    // INITIALIZE DM CONTROL ARRAY IF DOESN'T EXIST
    sprintf(name, "dm%lddisp", DMindex);
    IDdmctrl = image_ID(name);
    if(IDdmctrl==-1)
        IDdmctrl = read_sharedmem_image(name);
    if(IDdmctrl==-1)
    {
        dmsizearray = (long*) malloc(sizeof(long)*2);
        dmsizearray[0] = DMsize;
        dmsizearray[1] = DMsize;
        IDdmctrl = create_image_ID(name, 2, dmsizearray, FLOAT, 1, 0);
        free(dmsizearray);
        COREMOD_MEMORY_image_set_createsem(name, 2);
    }
    else
        {
            DMsize = data.image[IDdmctrl].md[0].size[0];
            COREMOD_MEMORY_image_set_createsem(name, 2);
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

    sprintf(name, "dm%lddisp", DMindex);
    AOsystSim_DMshape(name, "dmifc", "dm2Ddisp");
    IDdm0shape = image_ID("dm2Ddisp");
    save_fits("dm2Ddisp", "!dm2Ddisp.fits");





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
    OptSystProp_run(optsystsim, 0, 0, optsystsim[0].NBelem, "./testconf/", 1);

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
    sprintf(name, "dm%lddisp", DMindex);
    IDarray[1] = image_ID(name);

    sprintf(imdhname, "dhfield");
    dhsize = (long) (0.5/pupradcoeff * DMsize*(pupradcoeff/dmradcoeff)*0.5);
    dhxsize = dhsize;
    dhysize = dhsize*2;
    dhxoffset = arraysize/2; 
    dhyoffset = arraysize/2 - dhsize;
    dhsizearray = (long*) malloc(sizeof(long)*2);
    dhsizearray[0] = dhxsize*2;
    dhsizearray[1] = dhysize;
    
    IDdhmask = create_2Dimage_ID("dhmask", dhxsize, dhysize);
    for(ii=0;ii<dhxsize;ii++)
        for(jj=0;jj<dhysize;jj++)
            {
                data.image[IDdhmask].array.F[jj*dhxsize+ii] = 1.0;
                ii1 = (ii + dhxoffset) - arraysize/2;
                jj1 = (jj + dhyoffset) - arraysize/2;
                ii1ld = 2.0*ii1*pupradcoeff;
                r = sqrt(ii1*ii1+jj1*jj1);
                rld = r*pupradcoeff*2.0;
                if((ii1ld<0.5)||(rld<iwald))
                    data.image[IDdhmask].array.F[jj*dhxsize+ii] = 0.0;
            }

    save_fits("dhmask", "!dhmask.fits");

    iter = 0;
    while(1)
    {
//        printf("ITERATION %6ld   \n", iter);
  //      fflush(stdout);
        sprintf(name, "dm%lddisp", DMindex);
        AOsystSim_DMshape(name, "dmifc", "dm2Ddisp");
        OptSystProp_run(optsystsim, 0, 0, optsystsim[0].NBelem, "./testconf/", 1);
    
        // PYWFS code
        index = 2;
        ret = sprintf(imnameamp, "WFamp0_%03ld", index);
        ret = sprintf(imnamepha, "WFpha0_%03ld", index);
        mk_complex_from_amph(imnameamp, imnamepha, "_tmpwfc", 0);
        AOsystSim_WFSsim_Pyramid("_tmpwfc", "aosimwfsim", 0.0, 1);
        delete_image_ID("_tmpwfc");
        
        COREMOD_MEMORY_image_set_sempost("aosimwfsim", 0);

     
        ID = image_ID("psfi0");
        data.image[IDout].md[0].write = 1;
        memcpy(data.image[IDout].array.F, data.image[ID].array.F, sizeof(FLOAT)*data.image[ID].md[0].size[0]*data.image[ID].md[0].size[1]*data.image[ID].md[0].size[2]);
        data.image[IDout].md[0].cnt0++;
        data.image[IDout].md[0].write = 0;
        COREMOD_MEMORY_image_set_sempost("aosimpsfout", -1);

        

        // CREATE DARK HOLE FIELD
        IDre = image_ID("psfre0");
        IDim = image_ID("psfim0");
        IDdh = create_image_ID(imdhname, 2, dhsizearray, FLOAT, 1, 0);
        data.image[IDdh].md[0].write = 1;
        for(ii=0;ii<dhxsize;ii++)
            for(jj=0;jj<dhysize;jj++)
                {
                    data.image[IDdh].array.F[jj*(2*dhxsize)+ii] = data.image[IDre].array.F[(jj+dhyoffset)*arraysize + (ii+dhxoffset)] * data.image[IDdhmask].array.F[jj*dhxsize+ii];
                    data.image[IDdh].array.F[jj*(2*dhxsize)+(ii+dhxsize)] = data.image[IDim].array.F[(jj+dhyoffset)*arraysize + (ii+dhxoffset)] * data.image[IDdhmask].array.F[jj*dhxsize+ii];
                }
        data.image[IDdh].md[0].cnt0++;
        data.image[IDdh].md[0].write = 0;
        

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
            printf("WAITING %ld us\n", delayus);
            usleep(delayus);
            break;
        }
        
        COREMOD_MEMORY_image_set_semflush_IDarray(IDarray, 2);
        iter++;
    }

    free(dhsizearray);

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











double f_eval (const gsl_vector *v, void *params)
{
    double *p = (double *)params;
    long double value;
    long k;
    long pr;

    long double ptre_test;
    long double ptim_test;
    long double Iflux_test;
    long double are_test = 1.0;
    long double aim_test = 0.0;
    long double e_test = 1.0;
    long double ai;
    long double re, im, x, y, tmp1;


    ptre_test = gsl_vector_get(v, 0);
    ptim_test = gsl_vector_get(v, 1);
    Iflux_test = gsl_vector_get(v, 2);

    are_test = 1.0;
    aim_test = 0.0;
    e_test = 1.0;
    if(NBoptVar>3)
        are_test = gsl_vector_get(v, 3);
    if(NBoptVar>4)
    {
        aim_test = gsl_vector_get(v, 4);
        e_test = gsl_vector_get(v, 5);
    }

    value = 0.0;
    for(pr=0; pr<NBprobesG; pr++)
    {
        re = probe_re[pr] - ptre_test;
        im = probe_im[pr] - ptim_test;
        if(NBoptVar>3)
        {
            x = re*are_test + im*aim_test;
            y = -re*aim_test + im*are_test;
            y = y*e_test;
            probe_tflux[pr] = Iflux_test + Cflux*(x*x+y*y);
        }
        else
            probe_tflux[pr] = Iflux_test + Cflux*(re*re+im*im);

        tmp1 = (probe_tflux[pr]-probe_nmflux[pr])/probe_nmnoise[pr];
        value += tmp1*tmp1;
    }

    tmpvalue1 = value;

    // penalize large (pte)
    value += pow(0.01*ptre_test*ptre_test, 4.0);
    // penalize large (ptim)
    value += pow(0.01*ptim_test*ptim_test, 4.0);


    if(NBoptVar>4)
    {
        ai = are_test*are_test + aim_test*aim_test;
        value += 0.0001*pow((ai-1.0)*(ai-1.0), 4.0);

        tmpvalue2 = value;

        // penalize large (e-1)
        value += 0.01*pow((e_test-1.0)*(e_test-1.0)*1.0, 8.0);
    }
    // printf("%lf %lf %lf    %lf %lf %lf   -> %g\n", (double) ptre_test, (double) ptim_test, (double) Iflux_test, (double) are_test, (double) aim_test, (double) e_test, value);
    // fflush(stdout);
    // usleep(100000);
    //exit(0);

    // value = tmpvalue1;// no regularization

    return ((double) value);

}




long AOsystSim_FPWFS_imsimul(double probeamp, double sepx, double sepy, double contrast, double wferramp, double totFlux, double DMgainErr, double RON, double CnoiseFloor)
{
    long size = 256;
    double puprad = 50.0;
    long IDpupa;
    long IDwf0; // initial WF
    long IDwfA; // probe A
    long IDwfB; // probe B
    long IDwf;
    long IDpsfC; // PSF cube
    long pr;
    long IDa;
    long ii, jj, ii1, jj1;
    double coeffA, coeffB;
    double x, y, x1, y1;
    double sincx, sincy;
    double CPAx, CPAy;
    double CPAstep;
    long ID, ID1;
    long xsize, ysize, xmin, xmax, ymin, ymax;

    double tot1 = 0.0;
    double peak = 0.0;

    // define WFS probe area
    double CPAxmin, CPAxmax, CPAymin, CPAymax;
    double pixscaleld;
    long IDtmp;
    char imname[200];
    char fname[200];
    double probeAmultcoeff = 1.0;
    double probeBmultcoeff = 1.0;
    double probeBphaseoffset = 1.0*M_PI/2.0;

    double dmactgain;
    long IDnoise;
    double val;
    


    printf("Creating images .... Contrast = %g  \n", contrast);
    fflush(stdout);


    IDpupa = make_subpixdisk("pupa", size, size, 0.5*size, 0.5*size, puprad);


    IDwf0 = image_ID("wf0");
    if(IDwf0==-1)
    {
        IDwf0 = create_2Dimage_ID("wf0", size, size);
        for(ii=0; ii<size*size; ii++)
            data.image[IDwf0].array.F[ii] = wferramp*(1.0-2.0*ran1());
        save_fl_fits("wf0", "!wf0.fits");
    }


    IDwfA = create_2Dimage_ID("wfA", size, size);
    IDwfB = create_2Dimage_ID("wfB", size, size);
    IDwf = create_2Dimage_ID("wf", size, size);

    IDpsfC = create_3Dimage_ID("psfC", size, size, NBprobesG);

    // initialize wf0

    // initialize wfA and wfB
    CPAxmin = 4.0;
    CPAxmax = 20.0;
    CPAymin = 4.0;
    CPAymax = 20.0;
    CPAstep = 0.2;

    probeAmultcoeff = 1.0; // 1.1
    probeBmultcoeff = 1.0; // 0.9
    probeBphaseoffset = 1.0*M_PI/2.0; // 0.8


    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = (1.0*ii-0.5*size)/puprad;
            y = (1.0*jj-0.5*size)/puprad;
            for(CPAx=CPAxmin; CPAx<CPAxmax; CPAx += CPAstep)
                for(CPAy=CPAymin; CPAy<CPAymax; CPAy += CPAstep)
                {
                    data.image[IDwfA].array.F[jj*size+ii] += probeAmultcoeff*probeamp*CPAstep*CPAstep*cos(M_PI*(x*CPAx+y*CPAy));
                    data.image[IDwfB].array.F[jj*size+ii] += probeBmultcoeff*probeamp*CPAstep*CPAstep*cos(M_PI*(x*CPAx+y*CPAy)+probeBphaseoffset);
                }
        }

    // DM gain error
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            dmactgain = 1.0 + DMgainErr*(1.0-2.0*ran1());
            data.image[IDwfA].array.F[jj*size+ii]*= dmactgain;
            data.image[IDwfB].array.F[jj*size+ii]*= dmactgain;
        }


    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            data.image[IDwfA].array.F[jj*size+ii] *= data.image[IDpupa].array.F[jj*size+ii];
            data.image[IDwfB].array.F[jj*size+ii] *= data.image[IDpupa].array.F[jj*size+ii];
        }
    save_fl_fits("wfA", "!wfA.fits");
    save_fl_fits("wfB", "!wfB.fits");
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = (1.0*ii-0.5*size)/puprad;
            y = (1.0*jj-0.5*size)/puprad;
            data.image[IDpupa].array.F[jj*size+ii] *= exp(-8.0*(x*x+y*y));
        }
    save_fl_fits("pupa", "!pupa.fits");


    for(pr=0; pr<NBprobesG; pr++)
    {
        if(pr==0)
            {
                coeffA = 0.0;
                coeffB = 0.0;
            }
        else
            {
                coeffA = cos(2.0*M_PI/(NBprobesG-1)*(pr-1)); //nprobe_re[pr];
                coeffB = sin(2.0*M_PI/(NBprobesG-1)*(pr-1)); //nprobe_im[pr];
            }

        for(ii=0; ii<size*size; ii++)
            data.image[IDwf].array.F[ii] = data.image[IDwf0].array.F[ii] + coeffA*data.image[IDwfA].array.F[ii] + coeffB*data.image[IDwfB].array.F[ii];

        sprintf(fname, "!DMprobe%02ld.fits", pr);
        save_fl_fits("wf", fname);

        mk_complex_from_amph("pupa", "wf", "wfc", 0);
        permut("wfc");
        do2dfft("wfc", "imc");
        permut("imc");
        mk_amph_from_complex("imc", "ima", "imp", 0);
        delete_image_ID("imc");
        delete_image_ID("imp");
        IDa = image_ID("ima");
        for(ii=0; ii<size*size; ii++)
        {
            data.image[IDpsfC].array.F[pr*size*size+ii] = data.image[IDa].array.F[ii]*data.image[IDa].array.F[ii];
            tot1 += data.image[IDpsfC].array.F[pr*size*size+ii];
        }
        delete_image_ID("ima");
    }

    // ADD COMPANIONS
    printf("Adding companions\n");
    fflush(stdout);
    IDtmp = create_3Dimage_ID("tmp3dim", size, size, NBprobesG);
    for(ii=0; ii<size*size*NBprobesG; ii++)
        data.image[IDtmp].array.F[ii] = data.image[IDpsfC].array.F[ii];

    for(pr=0; pr<NBprobesG; pr++)
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                // companion #1
                ii1 = ii + (long) sepx;
                jj1 = jj + (long) sepy;
                if((ii1>0)&&(ii1<size)&&(jj1>0)&&(jj1<size))
                    data.image[IDpsfC].array.F[pr*size*size+jj1*size+ii1] += contrast * data.image[IDtmp].array.F[pr*size*size+jj*size+ii];
            
                // companion #2, 3.333x fainter, -10pix in x, +10 pix in y
                ii1 = ii + (long) sepx - 10;
                jj1 = jj + (long) sepy + 10;
                if((ii1>0)&&(ii1<size)&&(jj1>0)&&(jj1<size))
                    data.image[IDpsfC].array.F[pr*size*size+jj1*size+ii1] += 0.3*contrast * data.image[IDtmp].array.F[pr*size*size+jj*size+ii];
                
                // companion #3, 10x fainter, +10pix in x, +10 pix in y
                ii1 = ii + (long) sepx + 10;
                jj1 = jj + (long) sepy + 10;
                if((ii1>0)&&(ii1<size)&&(jj1>0)&&(jj1<size))
                    data.image[IDpsfC].array.F[pr*size*size+jj1*size+ii1] += 0.1*contrast * data.image[IDtmp].array.F[pr*size*size+jj*size+ii];

                // companion #4, 33.333x fainter, +10pix in x, -10 pix in y
                ii1 = ii + (long) sepx + 10;
                jj1 = jj + (long) sepy - 10;
                if((ii1>0)&&(ii1<size)&&(jj1>0)&&(jj1<size))
                    data.image[IDpsfC].array.F[pr*size*size+jj1*size+ii1] += 0.03*contrast * data.image[IDtmp].array.F[pr*size*size+jj*size+ii];

               // companion #5, 100x fainter, -10pix in x, -10 pix in y
                ii1 = ii + (long) sepx - 10;
                jj1 = jj + (long) sepy - 10;
                if((ii1>0)&&(ii1<size)&&(jj1>0)&&(jj1<size))
                    data.image[IDpsfC].array.F[pr*size*size+jj1*size+ii1] += 0.01*contrast * data.image[IDtmp].array.F[pr*size*size+jj*size+ii];

            }
    }
    delete_image_ID("tmp3dim");


    for(ii=0; ii<size*size*NBprobesG; ii++)
        data.image[IDpsfC].array.F[ii] *= totFlux/tot1;
    peak = 0.0;
    for(ii=0; ii<size*size; ii++)
        if(data.image[IDpsfC].array.F[ii]>peak)
            peak = data.image[IDpsfC].array.F[ii];


    pixscaleld = size/(2.0*puprad);
    xmin = 0.5*size + pixscaleld*CPAxmin;
    xmax = 0.5*size + pixscaleld*CPAxmax;
    ymin = 0.5*size + pixscaleld*CPAymin;
    ymax = 0.5*size + pixscaleld*CPAymax;
    xsize = xmax-xmin;
    ysize = ymax-ymin;


    // CREATE PROBE AMPLITUDE IMAGE IN FOCAL PLANE





    printf("Cropping and adding photon noise\n");
    fflush(stdout);

    ID = create_3Dimage_ID("psfCcrop", xsize, ysize, NBprobesG);
    tot1 = 0.0;
    for(pr=0; pr<NBprobesG; pr++)
        for(ii1=0; ii1<xsize; ii1++)
            for(jj1=0; jj1<ysize; jj1++)
            {
                ii = ii1+xmin;
                jj = jj1+ymin;
                data.image[ID].array.F[pr*xsize*ysize+jj1*xsize+ii1] = data.image[IDpsfC].array.F[pr*size*size+jj*size+ii];

                tot1 += data.image[ID].array.F[pr*xsize*ysize+jj1*xsize+ii1];
            }

    // CREATE PROBE AMPLITUDE IMAGE IN FOCAL PLANE
    ID1 = create_2Dimage_ID("psfprobeampC", xsize, ysize);

    
    for(ii1=0; ii1<xsize; ii1++)
        for(jj1=0; jj1<ysize; jj1++)
        {
            tot1 = 0.0;
            for(pr=CENTERprobe; pr<NBprobesG; pr++)
                tot1 += data.image[ID].array.F[pr*xsize*ysize+jj1*xsize+ii1]/peak;
            tot1 /= (NBprobesG-CENTERprobe);
            data.image[ID1].array.F[jj1*xsize+ii1] = tot1 - data.image[ID].array.F[jj1*xsize+ii1]/peak;
        }
    save_fl_fits("psfprobeampC", "!psfprobeampC.fits");

    // noise image
    IDnoise = create_3Dimage_ID("psfCcropnCn", xsize, ysize, NBprobesG);

    put_poisson_noise("psfCcrop", "psfCcropn");
    // add readout noise
    ID = image_ID("psfCcropn");
    for(ii1=0; ii1<xsize; ii1++)
        for(jj1=0; jj1<ysize; jj1++)
            data.image[ID].array.F[jj1*xsize+ii1] += RON * gauss();
            
    
    save_fl_fits("psfCcrop", "!psfCcrop.fits");
    save_fl_fits("psfCcropn", "!psfCcropn.fits");

    ID = image_ID("psfCcropn");
    ID1 =  create_3Dimage_ID("psfCcropnC", xsize, ysize, NBprobesG);
    for(pr=0; pr<NBprobesG; pr++)
        for(ii1=0; ii1<xsize; ii1++)
            for(jj1=0; jj1<ysize; jj1++)
            {
                data.image[ID1].array.F[pr*xsize*ysize+jj1*xsize+ii1] = data.image[ID].array.F[pr*xsize*ysize+jj1*xsize+ii1]/peak;
                
                val = data.image[ID].array.F[pr*xsize*ysize+jj1*xsize+ii1];
                if(val<1.0)
                    val = 1.0;
                
                data.image[IDnoise].array.F[pr*xsize*ysize+jj1*xsize+ii1] = sqrt( val + RON*RON )/peak ;  // assuming photon noise + readout noise
                if(data.image[IDnoise].array.F[pr*xsize*ysize+jj1*xsize+ii1] < CnoiseFloor)
                    data.image[IDnoise].array.F[pr*xsize*ysize+jj1*xsize+ii1]  = CnoiseFloor;
            }
    save_fl_fits("psfCcropnC", "!psfCcropnC.fits");
    printf("Saving psfCcropnCn\n");
    save_fl_fits("psfCcropnCn", "!psfCcropnCnoise.fits");
   
    for(pr=0; pr<NBprobesG; pr++)
    {
        sprintf(imname, "psfC_%03ld", pr);
        ID = create_2Dimage_ID(imname, size, size);
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                data.image[ID].array.F[jj*size+ii] = data.image[IDpsfC].array.F[pr*size*size+jj*size+ii]/peak;
            }

        sprintf(fname, "!psfC_%03ld.fits", pr);
        save_fl_fits(imname, fname);
    }

    printf("PSFs creation is complete\n");
    fflush(stdout);

    return(ID1);
}






// modegeom:
// 
// last digit
// 0  : horizontal
// 1  : vertical 
// 2  : quadrants #1
// 3  : quadrants #2
//
//
//
//
int AOsystSim_FPWFS_mkprobes(char *IDprobeA_name, char *IDprobeB_name, long dmxsize, long dmysize, double CPAmax, double CPArmin, double CPArmax, double RMSampl, long modegeom)
{
    long IDdmA, IDdmB;
    double CPAstep = 0.1;
    double x, y, r;
    double CPAx, CPAy, CPAr;
    double CPAxmin, CPAxmax;
    double CPAymin, CPAymax;
    long dmsize;
    long ii, jj;
    double rms;
    long ID;
    long imsize;
    double pha;
    
    IDdmA = create_2Dimage_ID(IDprobeA_name, dmxsize, dmysize);
    IDdmB = create_2Dimage_ID(IDprobeB_name, dmxsize, dmysize);
    
    switch (modegeom) {
    case 0 :    // mode 1: horizontal
        CPAxmin = 0.5;
        CPAxmax = CPAmax;
        CPAymin = -CPAmax;
        CPAymax = CPAmax;
        break;
    case 1 : // mode 2: vertical
        CPAxmin = -CPAmax;
        CPAxmax = CPAmax;
        CPAymin = 0.5;
        CPAymax = CPAmax;
        break;
    case 2 : // quadrants #1
        CPAxmin = 0.5;
        CPAxmax = CPAmax;
        CPAymin = 0.5;
        CPAymax = CPAmax;
        break;
    case 3 : // quadrants #1
        CPAxmin = 0.5;
        CPAxmax = CPAmax;
        CPAymin = -0.5;
        CPAymax = -CPAmax;
        break;
    default :
        printf("ERROR: mode not supported\n");
        exit(0);
        break;
    }

    dmsize = dmxsize;
    if(dmysize>dmxsize)
        dmsize = dmysize;

    for(CPAx=CPAxmin; CPAx<CPAmax; CPAx += CPAstep)
        for(CPAy=CPAymin;CPAy<CPAymax; CPAy += CPAstep)
            {
                CPAr = sqrt(CPAx*CPAx+CPAy*CPAy);
                if((CPAr>CPArmin)&&(CPAr<CPArmax))
                    {
                        pha = -1.6*CPAx;
                        for(ii=0;ii<dmxsize;ii++)
                            for(jj=0;jj<dmysize;jj++)
                                {
                                    x = 2.0*(1.0*ii-0.5*dmxsize)/dmsize;
                                    y = 2.0*(1.0*jj-0.5*dmysize)/dmsize;
                                    data.image[IDdmA].array.F[jj*dmxsize+ii] += cos(M_PI*(x*CPAx+y*CPAy)+pha);
                                    data.image[IDdmB].array.F[jj*dmxsize+ii] += cos(M_PI*(x*CPAx+y*CPAy)+M_PI/2+pha);
                                }
                    }
            }

    rms = 0.0;
    for(ii=0;ii<dmxsize*dmysize;ii++)
        rms += data.image[IDdmA].array.F[ii]*data.image[IDdmA].array.F[ii];
    rms = sqrt(rms/(dmxsize*dmysize));
    
    for(ii=0;ii<dmxsize*dmysize;ii++)
        {
            data.image[IDdmA].array.F[ii] *= RMSampl/rms;
            data.image[IDdmB].array.F[ii] *= RMSampl/rms;
        }
    
    
    // TEST
    imsize = 5*dmsize;
    ID = make_disk("pupa", imsize, imsize, 0.5*imsize, 0.5*imsize, 0.45*dmsize);
    for(ii=0;ii<imsize;ii++)
        for(jj=0;jj<imsize;jj++)
        {
            x = (1.0*ii-0.5*imsize)/(0.45*dmsize);
            y = (1.0*jj-0.5*imsize)/(0.45*dmsize);
            r = sqrt(x*x+y*y);
            data.image[ID].array.F[jj*imsize+ii] *= 1.0; //exp(-r*r*4.0);
            if(r<0.3)
                data.image[ID].array.F[jj*imsize+ii] = 0.0;
        }
    ID = create_2Dimage_ID("pupp", imsize, imsize);
    for(ii=0;ii<dmxsize;ii++)
        for(jj=0;jj<dmysize;jj++)
            data.image[ID].array.F[(jj+(imsize-dmysize)/2)*imsize+(ii+(imsize-dmxsize)/2)] = data.image[IDdmA].array.F[jj*dmxsize+ii];
    
    mk_complex_from_amph("pupa", "pupp", "pupc", 0);
    permut("pupc");
    do2dfft("pupc","focc");
    permut("focc");
    mk_amph_from_complex("focc","foca","focp", 0);
    save_fits("pupa", "!test_pupa.fits");
    save_fits("pupp", "!test_pupp_A.fits");
    save_fits("foca", "!test_foca_A.fits");
    delete_image_ID("pupc");
    delete_image_ID("focc");
    delete_image_ID("foca");
    delete_image_ID("focp");
 
    
    
   ID = image_ID("pupp");
    for(ii=0;ii<dmxsize;ii++)
        for(jj=0;jj<dmysize;jj++)
            data.image[ID].array.F[(jj+(imsize-dmysize)/2)*imsize+(ii+(imsize-dmxsize)/2)] = data.image[IDdmB].array.F[jj*dmxsize+ii];
       mk_complex_from_amph("pupa", "pupp", "pupc", 0);
    permut("pupc");
    do2dfft("pupc","focc");
    permut("focc");
    mk_amph_from_complex("focc","foca","focp", 0);
    save_fits("pupp", "!test_pupp_B.fits");
    save_fits("foca", "!test_foca_B.fits");
    delete_image_ID("pupc");
    delete_image_ID("focc");
    delete_image_ID("foca");
    delete_image_ID("focp");
 

  ID = image_ID("pupp");
    for(ii=0;ii<dmxsize;ii++)
        for(jj=0;jj<dmysize;jj++)
            data.image[ID].array.F[(jj+(imsize-dmysize)/2)*imsize+(ii+(imsize-dmxsize)/2)] = -data.image[IDdmA].array.F[jj*dmxsize+ii];
       mk_complex_from_amph("pupa", "pupp", "pupc", 0);
    permut("pupc");
    do2dfft("pupc","focc");
    permut("focc");
    mk_amph_from_complex("focc","foca","focp", 0);
    save_fits("pupp", "!test_pupp_mA.fits");
    save_fits("foca", "!test_foca_mA.fits");
    delete_image_ID("pupc");
    delete_image_ID("focc");
    delete_image_ID("foca");
    delete_image_ID("focp");
 
 
  ID = image_ID("pupp");
    for(ii=0;ii<dmxsize;ii++)
        for(jj=0;jj<dmysize;jj++)
            data.image[ID].array.F[(jj+(imsize-dmysize)/2)*imsize+(ii+(imsize-dmxsize)/2)] = -data.image[IDdmB].array.F[jj*dmxsize+ii];
       mk_complex_from_amph("pupa", "pupp", "pupc", 0);
    permut("pupc");
    do2dfft("pupc","focc");
    permut("focc");
    mk_amph_from_complex("focc","foca","focp", 0);
    save_fits("pupp", "!test_pupp_mB.fits");
    save_fits("foca", "!test_foca_mB.fits");
    delete_image_ID("pupc");
    delete_image_ID("focc");
    delete_image_ID("foca");
    delete_image_ID("focp");

   
   
 
    ID = image_ID("pupp");
    for(ii=0;ii<dmxsize;ii++)
        for(jj=0;jj<dmysize;jj++)
            data.image[ID].array.F[(jj+(imsize-dmysize)/2)*imsize+(ii+(imsize-dmxsize)/2)] = 0.0;
    mk_complex_from_amph("pupa", "pupp", "pupc", 0);
    permut("pupc");
    do2dfft("pupc","focc");
    permut("focc");
    mk_amph_from_complex("focc","foca","focp", 0);
    save_fits("pupp", "!test_pupp_00.fits");
    save_fits("foca", "!test_foca_00.fits");
    delete_image_ID("pupc");
    delete_image_ID("focc");
    delete_image_ID("foca");
    delete_image_ID("focp");

   
   
    
    return(0);
}


//
// explore optimal theoretical sensitivity for focal plane WFS
//
// optmode:
// 1: 3 free parameters: ptre, ptim, Iflux
// 2: 6 free parameters: ptre, ptim, Iflux, are, aim, e
//
int AOsystSim_FPWFS_sensitivityAnalysis(int mapmode, int mode, int optmode, int NBprobes)
{
    // conventions
    // CA : complex amplitude
    //
    // CA of point to be measured is uniformly distributed in unit circle

    char imname[200];
    double ptre, ptim; // real/imaginary components of point to be probed

    double probe_noise_prop;
    double ProbeNoise;

    double probe_mflux[100]; // measured flux (ideal, no noise)
    double probe_mflux_dre[100]; // derivative against ptre
    double probe_mflux_dim[100]; // derivative against ptim
    double probe_mflux_dI[100]; // derivative against Iflux

    double probeamp;
    int pr;
    double Iflux = 0.0; // incoherent flux
    double re, im;
    double eps = 1.0e-8;

    double re_re, re_im, im_re, im_im;
    long k, ks;
    double tmp1;
    int execmode;

    // solver
    double ptre_test, ptim_test, Iflux_test, are_test, aim_test, e_test;
    double ptre_best, ptim_best, Iflux_best, are_best, aim_best, e_best;
    double value, bvalue, value0;
    double x, y;

    double optampl;
    long ksmax = 100000000;

    double mapampl = 2.0;
    long mapsize = 41; // map size (linear)
    long mapxsize, mapysize;
    long kx, ky, kz;

    long IDmap;

    // output solution
    long IDmap_ptre;
    long IDmap_ptim;
    long IDmap_ptre_in;
    long IDmap_ptim_in;
    long IDmap_Iflux;
    long IDmap_are;
    long IDmap_aim;
    long IDmap_e;

    long mapz;
    long mapzsize;


    long st;
    long optvar;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *xvect;
    gsl_multimin_function feval_func;
    int status;
    long iter;
    long iterMax = 100000;
    double optsize;
    double bestvalue;
    double tmp;

    double FLUXph = 1.0e4; // total number of photon per cycle

    int imSimulMode = 0; // 1 if real image simulated
    long IDpsfC;

    long avecnt;
    double ave, rms;
    long ii, jj;
    long kxtest = 21;
    long kytest = 11;

    long IDprobampC = -1;

    long vID;
    double probeampl;
    double WFerr;
    double contrast;
    double totFlux;

    int initmap = 0;

    long IDmap_Iflux_ave, IDmap_Iflux_rms;
    int NewWF = 0;
    long IDpsfCnoise;

    double RON = 1.0; // detector readout noise (phe-)
    long IDmap_Iflux_rmsn;
    long IDmap_Iflux_rmsn1;
    long IDmap_CA_rms, IDmap_CA_rmsn;
    double dx, dy, val;
    double CnoiseFloor;
    

    NBprobesG = NBprobes;



    if(mapmode==2) // complex amplitude field
        imSimulMode = 1;

    mapzsize = 10;
    if((vID=variable_ID("mapzsize"))!=-1)
    {
        mapzsize = (long) (0.1+data.variable[vID].value.f);
        printf("mapzsize = %ld\n", mapzsize);
    }

    if((vID=variable_ID("NewWF"))!=-1)
    {
        NewWF = (int) (0.1+data.variable[vID].value.f);
        printf("NewWF = %d\n", NewWF);
    }

    ProbeNoise = 0.2;
    if((vID=variable_ID("ProbeNoise"))!=-1)
    {
        ProbeNoise = data.variable[vID].value.f;
        printf("ProbeNoise = %g\n", ProbeNoise);
    }

    RON = 1.0;
    if((vID=variable_ID("RON"))!=-1)
    {
        RON = data.variable[vID].value.f;
        printf("RON = %f\n", RON);
    }

    CnoiseFloor = 1.0e-5;
    if((vID=variable_ID("CnoiseFloor"))!=-1)
    {
        CnoiseFloor = data.variable[vID].value.f;
        printf("CnoiseFloor = %g\n", CnoiseFloor);
    }

    CENTERprobe = 1;
    if((vID=variable_ID("CENTERprobe"))!=-1)
    {
        CENTERprobe = (int) (0.1+data.variable[vID].value.f);
        printf("CENTERprobe = %d\n", CENTERprobe);
    }




    if(imSimulMode == 1)
    {
        probeampl = 0.0001;
        if((vID=variable_ID("probeampl"))!=-1)
        {
            probeampl = data.variable[vID].value.f;
            printf("probeamp = %f rad\n", probeampl);
        }

        contrast = 5.0e-8;
        if((vID=variable_ID("contrast"))!=-1)
        {
            contrast = data.variable[vID].value.f;
            printf("contrast = %e\n", contrast);
        }

        WFerr = 0.000;
        if((vID=variable_ID("WFerr"))!=-1)
        {
            WFerr = data.variable[vID].value.f;
            printf("WFerr = %f rad\n", WFerr);
        }

        mapmode = 2; // use existing file
    }

        
    totFlux= 1e12;
        if((vID=variable_ID("totFlux"))!=-1)
        {
            totFlux = data.variable[vID].value.f;
            FLUXph = data.variable[vID].value.f;
            printf("totFlux = %f ph\n", totFlux);
        }



    if(mapmode==2)
        CENTERprobe = 1;


    probeamp = 1.0;
    if(CENTERprobe==1) // include center probe
    {
        probe_re[pr] = 0.0;
        probe_im[pr] = 0.0;
        for(pr=1; pr<NBprobes; pr++)
        {
            probe_re[pr] = probeamp*cos(2.0*M_PI/(NBprobes-1)*(pr-1));
            probe_im[pr] = probeamp*sin(2.0*M_PI/(NBprobes-1)*(pr-1));
        }
    }
    else
    {
        for(pr=0; pr<NBprobes; pr++)
        {
            probe_re[pr] = probeamp*cos(2.0*M_PI/(NBprobes)*pr);
            probe_im[pr] = probeamp*sin(2.0*M_PI/(NBprobes)*pr);
        }
    }
    
    
    execmode = 1;





    if(mapmode==0)
    {
        mapsize = 1;
    }

    mapxsize = mapsize;
    mapysize = mapsize;







    switch (mode) {
    case 1 :    // mode 1: perfectly calibrated system -> no probe noise
        probe_noise_prop = 0.0;
        for(pr=0; pr<NBprobes; pr++)
        {
            nprobe_re[pr] = probe_re[pr];
            nprobe_im[pr] = probe_im[pr];
        }
        break;
    case 2 : // mode 2: uncorrelated probes noise
        probe_noise_prop = ProbeNoise;
        for(pr=0; pr<NBprobes; pr++)
        {
            nprobe_re[pr] = probe_re[pr] + probe_noise_prop*probeamp*(2.0*ran1()-1.0);
            nprobe_im[pr] = probe_im[pr] + probe_noise_prop*probeamp*(2.0*ran1()-1.0);
        }
        break;
    case 3 : // noise on probe axes
        probe_noise_prop = ProbeNoise;
        re_re = 1.0 + (2.0*ran1()-1.0)*probe_noise_prop;
        re_im = 0.0 + (2.0*ran1()-1.0)*probe_noise_prop;
        im_re = 0.0 + (2.0*ran1()-1.0)*probe_noise_prop;
        im_im = 1.0 + (2.0*ran1()-1.0)*probe_noise_prop;
        for(pr=0; pr<NBprobes; pr++)
        {
            nprobe_re[pr] = probe_re[pr]*re_re + probe_im[pr]*im_re;
            nprobe_im[pr] = probe_re[pr]*re_im + probe_im[pr]*im_im;
        }
        break;
    default :
        printf("ERROR: mode not supported\n");
        execmode = 0;
        break;
    }


    printf("noise mode  = %d (%f)\n", mode, probe_noise_prop);
    printf("mapmode     = %d\n", mapmode);
    printf("mapsize     = %ld   %ld   %ld\n", mapxsize, mapysize, mapzsize);
    printf("execmode    = %d\n", execmode);
    printf("mapzsize    = %ld\n", mapzsize);
    printf("imSimulMode = %d\n", imSimulMode);


    if(IDprobampC==-1)
        IDprobampC = create_2Dimage_ID("psfprobeamp", mapxsize, mapysize);
    for(ii=0;ii<mapxsize*mapysize;ii++)
        data.image[IDprobampC].array.F[ii] = 1.0;

    printf("\n\n");



    for(mapz=0; mapz<mapzsize; mapz++)
    {

        if(imSimulMode == 1)
        {
            IDpsfC = AOsystSim_FPWFS_imsimul(probeampl, 30.0, 30.0, contrast, WFerr, totFlux, probe_noise_prop, RON, CnoiseFloor); // computes data cube
            IDpsfCnoise = image_ID("psfCcropnCn");
            IDprobampC = image_ID("psfprobeampC");
            save_fl_fits("psfC", "!psfC.fits");

            ave = 0.0;
            avecnt = 0;

            for(pr=1; pr<NBprobes; pr++)
            {
                for(ii=0; ii<data.image[IDpsfC].md[0].size[0]; ii++)
                    for(jj=0; jj<data.image[IDpsfC].md[0].size[1]; jj++)
                    {
                        ave += data.image[IDpsfC].array.F[pr*data.image[IDpsfC].md[0].size[0]*data.image[IDpsfC].md[0].size[1]+jj*data.image[IDpsfC].md[0].size[0]+ii];
                        avecnt++;
                    }
            }
            ave /= avecnt;
            printf("ave = %g\n", ave);

            mapmode = 2; // use existing file
        }


        if(mapmode == 2)
        {
            mapxsize = data.image[IDpsfC].md[0].size[0];
            mapysize = data.image[IDpsfC].md[0].size[1];
        }

        if(initmap==0)
        {
            if(mapmode>0)
            {
                IDmap = create_3Dimage_ID("WFSerrmap", mapxsize, mapysize, mapzsize);
                IDmap_ptre = create_3Dimage_ID("WFSsol_ptre", mapxsize, mapysize, mapzsize);
                IDmap_ptim = create_3Dimage_ID("WFSsol_ptim", mapxsize, mapysize, mapzsize);
                IDmap_ptre_in = create_3Dimage_ID("WFSsol_ptre_in", mapxsize, mapysize, mapzsize);
                IDmap_ptim_in = create_3Dimage_ID("WFSsol_ptim_in", mapxsize, mapysize, mapzsize);
                IDmap_Iflux = create_3Dimage_ID("WFSsol_Iflux", mapxsize, mapysize, mapzsize);
                IDmap_are = create_3Dimage_ID("WFSsol_are", mapxsize, mapysize, mapzsize);
                IDmap_aim = create_3Dimage_ID("WFSsol_aim", mapxsize, mapysize, mapzsize);
                IDmap_e = create_3Dimage_ID("WFSsol_e", mapxsize, mapysize, mapzsize);
            }
            else
                IDmap = -1;
            initmap = 1;
        }





        if(mapmode==0)
        {
            printf("Preparing optimization ... \n");
            fflush(stdout);
        }
        switch (optmode) {
        case 0 :
            NBoptVar = 3;
            break;
        case 1 :
            NBoptVar = 4;
            break;
        case 2 :
            NBoptVar = 6;
            break;
        }

        xvect = gsl_vector_alloc (NBoptVar);
        ss = gsl_vector_alloc (NBoptVar);
        feval_func.n = NBoptVar;  /* number of function components */
        feval_func.f = &f_eval;
        feval_func.params = (void *) NULL;
        s = gsl_multimin_fminimizer_alloc (T, NBoptVar);


        if(execmode>0)
        {
            if((mapmode==0)&&((kx==kxtest)&&(ky==kytest)))
                for(pr=0; pr<NBprobes; pr++)
                    printf("PROBE %2d : %10lf %10lf\n", pr, (double) probe_re[pr], (double) probe_im[pr]);


            for(kx=0; kx<mapxsize; kx++)
                for(ky=0; ky<mapysize; ky++)
                {
                    printf("\r slice %ld    pixel %5ld %5ld       ", mapz, kx, ky);
                    fflush(stdout);
                    if(mapmode==0)
                    {
                        ptre = 2.0*ran1()-1.0;
                        ptim = 2.0*ran1()-1.0;
                    }
                    else
                    {
                        ptre = mapampl*(2.0*kx/(mapxsize-1)-1.0);
                        ptim = mapampl*(2.0*ky/(mapysize-1)-1.0);
                    }
                    if(mapmode==2)
                    {
                        ptre = 0.0;
                        ptim = 0.0;
                    }

                    if((mapmode==0)||((kx==kxtest)&&(ky==kytest))||((kx==kxtest+1)&&(ky==kytest)))
                        printf("\n\n ptre, ptim = %g, %g\n", ptre, ptim);
                    // compute flux
                    
                    
                  
                    for(pr=0; pr<NBprobes; pr++)
                    {
                        // ideal measurements (no noise)
                        re = probe_re[pr] - ptre;
                        im = probe_im[pr] - ptim;
                        probe_mflux[pr] = Iflux + Cflux*(re*re+im*im);


                        // noisy measurements
                        re = nprobe_re[pr] - ptre;
                        im = nprobe_im[pr] - ptim;
                        val = fast_poisson(((Iflux + Cflux*(re*re+im*im))*FLUXph)/NBprobes) + RON*gauss();
                        probe_nmflux[pr] = ( val ) / (FLUXph/NBprobes);
                        if(val<1.0)
                            val = 1.0;
                        probe_nmnoise[pr] = sqrt(val);
                        if((mapmode==0)&&(((kx==kxtest)&&(ky==kytest))||((kx==kxtest+1)&&(ky==kytest))))
                        {
                            printf("M0    NO NOISE:  probe %3d  normalized flux = %8.5lf  (%g ph)\n", pr, (double) (probe_mflux[pr]/Cflux), (double) probe_mflux[pr]);
                            printf("M0  WITH NOISE:             normalized flux = %8.5lf  (%g ph)\n",     (double) (probe_nmflux[pr]/Cflux), (double) probe_nmflux[pr]);
                            fflush(stdout);
                        }


                        if(imSimulMode == 1)
                        {
                            probe_nmflux[pr] = data.image[IDpsfC].array.F[pr*mapxsize*mapysize+ky*mapxsize+kx] / data.image[IDprobampC].array.F[ky*mapxsize+kx]; // unit : normalized contrast
                            probe_nmnoise[pr] = data.image[IDpsfCnoise].array.F[pr*mapxsize*mapysize+ky*mapxsize+kx] / data.image[IDprobampC].array.F[ky*mapxsize+kx];
                            if(((kx==kxtest)&&(ky==kytest))||((kx==kxtest+1)&&(ky==kytest)))
                                printf("PROBE %d   -> %g\n", pr, (double) probe_nmflux[pr]);
                        }
                        
                    }

                    if((mapmode==0)&&((kx==kxtest)&&(ky==kytest))||((kx==kxtest+1)&&(ky==kytest)))
                        fflush(stdout);

                    // SOLVER


                    for(optvar=0; optvar<NBoptVar; optvar++)
                        gsl_vector_set (xvect, optvar, 0.0);

                    gsl_vector_set (xvect, 0, ptre);
                    gsl_vector_set (xvect, 1, ptim);

                    if(NBoptVar>3)
                        gsl_vector_set (xvect, 3, 1.0);


                    if(NBoptVar>4)
                    {
                        gsl_vector_set (xvect, 4, 0.0);
                        gsl_vector_set (xvect, 5, 1.0);
                    }

                    /* Set initial step sizes  */
                    gsl_vector_set_all (ss, 0.1);

                    /* Initialize method and iterate */
                    iter = 0;
                    bestvalue = 100000.0;


                    gsl_multimin_fminimizer_set (s, &feval_func, xvect, ss);


                    do
                    {
                        iter++;
                        status = gsl_multimin_fminimizer_iterate(s);
                        if (status)
                            break;

                        optsize = gsl_multimin_fminimizer_size (s);
                        // printf("Iteration %ld   %g\n", iter, optsize);
                        status = gsl_multimin_test_size (optsize, 1.0/pow(10.0, 10.0-4.0*iter/iterMax)); //1.0e-10);


                        //printf ("............[%05ld] ->  %e   (%e)\n", iter, s->fval, bestvalue);
                        if (status == GSL_SUCCESS)
                        {
                            if(((kx==kxtest)&&(ky==kytest))||((kx==kxtest+1)&&(ky==kytest)))
                            {
                                printf ("  ->  %e [%5ld]  (%e)", s->fval, iter, bestvalue);
                                printf("\n");
                            }
                            ptre_best = gsl_vector_get(s->x, 0);
                            ptim_best = gsl_vector_get(s->x, 1);
                            Iflux_best = gsl_vector_get(s->x, 2);
                            if(NBoptVar>3)
                                are_best = gsl_vector_get(s->x, 3);
                            if(NBoptVar>4)
                            {
                                aim_best = gsl_vector_get(s->x, 4);
                                e_best = gsl_vector_get(s->x, 5);

                                if(e_best>1.0)
                                {
                                    e_best = 1.0/e_best;
                                    tmp = are_best;
                                    are_best = aim_best;
                                    aim_best = -tmp;
                                }
                            }



                            if(s->fval < bestvalue)
                                bestvalue = s->fval;
                        }
                    }
                    while (status == GSL_CONTINUE && iter < iterMax );
                    if(iter>iterMax-1) {
                        printf("Max number of iteration reached (optsize = %g)\n", optsize);
                    }

                    ptre_best = gsl_vector_get(s->x, 0);
                    ptim_best = gsl_vector_get(s->x, 1);
                    Iflux_best = gsl_vector_get(s->x, 2);
                    if(NBoptVar>3)
                        are_best = gsl_vector_get(s->x, 3);
                    if(NBoptVar>4)
                    {
                        aim_best = gsl_vector_get(s->x, 4);
                        e_best = gsl_vector_get(s->x, 5);
                        if(e_best>1.0)
                        {
                            e_best = 1.0/e_best;
                            tmp = are_best;
                            are_best = aim_best;
                            aim_best = -tmp;
                        }
                    }

                    if((mapmode==0)||((kx==kxtest)&&(ky==kytest))||((kx==kxtest+1)&&(ky==kytest)))
                    {
                        printf("\n\n");
                        printf("[%3ld %3ld] OPTIMAL SOLUTION  [ %6ld / %6ld  %g ]: \n", kx, ky, iter, iterMax, optsize);
                        printf("       ptre = %.18f\n", ptre_best);
                        printf("       ptim = %.18f\n", ptim_best);
                        printf("      Iflux = %.18f\n", Iflux_best);
                        if(NBoptVar>3)
                            printf("        are = %.18f\n", are_best);
                        if(NBoptVar>4)
                        {
                            printf("        aim = %.18f\n", aim_best);
                            printf("          e = %.18f\n", e_best);
                        }
                        printf("OPT VALUES :  %g  %g  %g\n", tmpvalue1, tmpvalue2, bestvalue);
                        printf("\n\n");
                    }


                    if(mapmode>0)
                    {
                        
                        data.image[IDmap].array.F[mapz*mapxsize*mapysize+ky*mapxsize+kx] = Iflux_best*sqrt(data.image[IDprobampC].array.F[ky*mapxsize+kx]);
                        data.image[IDmap_ptre].array.F[mapz*mapysize*mapxsize+ky*mapxsize+kx] = ptre_best*sqrt(data.image[IDprobampC].array.F[ky*mapxsize+kx]);
                        data.image[IDmap_ptim].array.F[mapz*mapysize*mapxsize+ky*mapxsize+kx] = ptim_best*sqrt(data.image[IDprobampC].array.F[ky*mapxsize+kx]);
                        data.image[IDmap_ptre_in].array.F[mapz*mapysize*mapxsize+ky*mapxsize+kx] = ptre*sqrt(data.image[IDprobampC].array.F[ky*mapxsize+kx]);
                        data.image[IDmap_ptim_in].array.F[mapz*mapysize*mapxsize+ky*mapxsize+kx] = ptim*sqrt(data.image[IDprobampC].array.F[ky*mapxsize+kx]);
                        data.image[IDmap_Iflux].array.F[mapz*mapysize*mapxsize+ky*mapxsize+kx] = Iflux_best*data.image[IDprobampC].array.F[ky*mapxsize+kx];                                          
                        if(NBoptVar>3)
                            data.image[IDmap_are].array.F[mapz*mapysize*mapxsize+ky*mapxsize+kx] = are_best*sqrt(data.image[IDprobampC].array.F[ky*mapxsize+kx]);
                        if(NBoptVar>4)
                        {
                            data.image[IDmap_aim].array.F[mapz*mapysize*mapxsize+ky*mapxsize+kx] = aim_best*sqrt(data.image[IDprobampC].array.F[ky*mapxsize+kx]);
                            data.image[IDmap_e].array.F[mapz*mapysize*mapxsize+ky*mapxsize+kx] = e_best;
                        }
                    }
                }
        }
        if(mapmode>0)
        {
            save_fits("WFSerrmap", "!WFSerrmap.fits");
            save_fits("WFSsol_ptre", "!WFSsol_ptre.fits");
            save_fits("WFSsol_ptim", "!WFSsol_ptim.fits");
            save_fits("WFSsol_ptre_in", "!WFSsol_ptre_in.fits");
            save_fits("WFSsol_ptim_in", "!WFSsol_ptim_in.fits");
            save_fits("WFSsol_Iflux", "!WFSsol_Iflux.fits");
            if(NBoptVar>3)
                save_fits("WFSsol_are", "!WFSsol_are.fits");
            if(NBoptVar>4)
            {
                save_fits("WFSsol_aim", "!WFSsol_aim.fits");
                save_fits("WFSsol_e", "!WFSsol_e.fits");
            }
        }


        delete_image_ID("pupa");

        if(NewWF==1)
            delete_image_ID("wf0");

        delete_image_ID("wfA");
        delete_image_ID("wfB");
        delete_image_ID("wf");
        delete_image_ID("psfC");
        delete_image_ID("wfc");
        delete_image_ID("psfCcrop");
        delete_image_ID("psfprobeampC");
        delete_image_ID("psfCcropn");
        delete_image_ID("psfCcropnC");
        for(pr=0; pr<NBprobes; pr++)
        {
            sprintf(imname, "psfC_%03d", pr);
            delete_image_ID(imname);
        }

        list_image_ID();



        if((mapz>2)&&(mapmode>0))
        {
            IDmap_Iflux_ave = image_ID("WFSsol_Iflux_ave");
            if(IDmap_Iflux_ave==-1)
                IDmap_Iflux_ave = create_2Dimage_ID("WFSsol_Iflux_ave", mapxsize, mapysize);

            IDmap_Iflux_rms = image_ID("WFSsol_Iflux_rms");
            if(IDmap_Iflux_rms==-1)
                IDmap_Iflux_rms = create_2Dimage_ID("WFSsol_Iflux_rms", mapxsize, mapysize);

            


            for(kx=0; kx<mapxsize; kx++)
                for(ky=0; ky<mapysize; ky++)
                {
                    ave = 0.0;
                    rms = 0.0;
                    for(kz=0; kz<mapz; kz++)
                    {
                        ave += data.image[IDmap_Iflux].array.F[kz*mapysize*mapxsize+ky*mapxsize+kx];
                        rms += data.image[IDmap_Iflux].array.F[kz*mapysize*mapxsize+ky*mapxsize+kx]*data.image[IDmap_Iflux].array.F[kz*mapysize*mapxsize+ky*mapxsize+kx];
                    }
                    ave /= mapz;
                    rms /= mapz;
                    rms -= ave*ave;
                    rms = sqrt(rms);
                    data.image[IDmap_Iflux_ave].array.F[ky*mapxsize+kx] = ave;
                    data.image[IDmap_Iflux_rms].array.F[ky*mapxsize+kx] = rms;
                }
            save_fits("WFSsol_Iflux_ave", "!WFSsol_Iflux_ave.fits");
            save_fits("WFSsol_Iflux_rms", "!WFSsol_Iflux_rms.fits");
        
        

        
        
            if(imSimulMode == 0) 
            {                
                // INCOHERENT COMPONENT ERROR

                IDmap_Iflux_rmsn = image_ID("WFSsol_Iflux_rmsn");
                if(IDmap_Iflux_rmsn==-1)
                    IDmap_Iflux_rmsn = create_2Dimage_ID("WFSsol_Iflux_rmsn", mapxsize, mapysize);
                
                IDmap_Iflux_rmsn1 = image_ID("WFSsol_Iflux_rmsn1");
                if(IDmap_Iflux_rmsn1==-1)
                    IDmap_Iflux_rmsn1 = create_2Dimage_ID("WFSsol_Iflux_rmsn1", mapxsize, mapysize);
                
                
                for(kx=0; kx<mapxsize; kx++)
                    for(ky=0; ky<mapysize; ky++)
                    {
                        ptre = mapampl*(2.0*kx/mapxsize-1.0);
                        ptim = mapampl*(2.0*ky/mapysize-1.0);
                        
                        data.image[IDmap_Iflux_rmsn].array.F[ky*mapxsize+kx] = data.image[IDmap_Iflux_rms].array.F[ky*mapxsize+kx] / sqrt(1.0/FLUXph);
                        data.image[IDmap_Iflux_rmsn1].array.F[ky*mapxsize+kx] = data.image[IDmap_Iflux_rms].array.F[ky*mapxsize+kx] / sqrt(1.0/FLUXph) / sqrt(ptre*ptre+ptim*ptim);
                    }
                save_fits("WFSsol_Iflux_rmsn", "!WFSsol_Iflux_rmsn.fits");
                save_fits("WFSsol_Iflux_rmsn1", "!WFSsol_Iflux_rmsn1.fits");
            
              
              
                // COHERENT COMPONENT ERROR
                
                IDmap_CA_rms = image_ID("WFSsol_CA_rms");
                if(IDmap_CA_rms==-1)
                    IDmap_CA_rms = create_2Dimage_ID("WFSsol_CA_rms", mapxsize, mapysize);

                for(kx=0; kx<mapxsize; kx++)
                for(ky=0; ky<mapysize; ky++)
                {
                    ave = 0.0;
                    rms = 0.0;
                    for(kz=0; kz<mapz; kz++)
                    {
                        dx = data.image[IDmap_ptre].array.F[kz*mapysize*mapxsize+ky*mapxsize+kx]-data.image[IDmap_ptre_in].array.F[kz*mapysize*mapxsize+ky*mapxsize+kx];
                        dy = data.image[IDmap_ptim].array.F[kz*mapysize*mapxsize+ky*mapxsize+kx]-data.image[IDmap_ptim_in].array.F[kz*mapysize*mapxsize+ky*mapxsize+kx];
                        val = dx*dx+dy*dy;
                        ave += val;
                        rms += val*val;
                    }
                    ave /= mapz;
                    rms /= mapz;
                    rms -= ave*ave;
                    rms = sqrt(rms);
                    data.image[IDmap_CA_rms].array.F[ky*mapxsize+kx] = sqrt(ave);
                }
                save_fits("WFSsol_CA_rms", "!WFSsol_CA_rms.fits");
        
                IDmap_CA_rmsn = image_ID("WFSsol_CA_rmsn");
                if(IDmap_CA_rmsn==-1)
                    IDmap_CA_rmsn = create_2Dimage_ID("WFSsol_CA_rmsn", mapxsize, mapysize);
                
                for(kx=0; kx<mapxsize; kx++)
                    for(ky=0; ky<mapysize; ky++)
                        data.image[IDmap_CA_rmsn].array.F[ky*mapxsize+kx] = data.image[IDmap_CA_rms].array.F[ky*mapxsize+kx] * sqrt(FLUXph);
                
                save_fits("WFSsol_CA_rmsn", "!WFSsol_CA_rmsn.fits");                
            }
        }

    }

    return(0);
}
















