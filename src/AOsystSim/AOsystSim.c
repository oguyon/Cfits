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

#include "OptSystProp/OptSystProp.h"
#include "AOsystSim/AOsystSim.h"

extern DATA data;

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
  
  AOsystSim_run();

  return 0;
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
	for(k=0;k<wfcbuff_size;k++)
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
				for(k1=0;k1<dmmoveNBpt;k1++)
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

	IDdmctrl = image_ID(IDdmctrl_name);

	IDdmifc = image_ID(IDdmifc_name);
	dmsizex = data.image[IDdmifc].md[0].size[0];
	dmsizey = data.image[IDdmifc].md[0].size[1];
	DMnbact = data.image[IDdmifc].md[0].size[2];

	IDdm = image_ID(IDdm_name);
	if(IDdm==-1)
	{
		IDdm = create_2Dimage_ID(IDdm_name, dmsizex, dmsizey);
	}
	for(ii=0;ii<dmsizex*dmsizey;ii++)
		data.image[IDdm].array.F[ii] = 0.0;
	
	for(dmact=0;dmact<DMnbact;dmact++)
	{
		for(ii=0;ii<dmsizex*dmsizey;ii++)
			data.image[IDdm].array.F[ii] += data.image[IDdmctrl].array.F[dmact]*data.image[IDdmifc].array.F[dmact*dmsizex*dmsizey+jj*dmsizex+ii];
	}
	
	return (0);
}







/** \brief simplified AO system simulator (closed loop part)
 *
 * creates a DM map(s) and a WF error input
 * When either DM map or WF error input changes, compute intensity outputs (images)
 *
 */

int AOsystSim_run()
{
    long arraysize = 512;
    long ii, jj;
    double puprad, dmrad;
    double dmifscale = 20.0; // scale magnification between full DM map and DM influence function
    long *imsize;
    long DMsize = 50;
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

    puprad = 0.12*arraysize;
    dmrad = 0.13*arraysize;



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
    printf("actuator pix size = %f pix\n", dmifscale*(2.0*dmrad/DMsize));
    for(ii=0; ii<arraysize; ii++)
        for(jj=0; jj<arraysize; jj++)
        {
            x = (1.0*ii-0.5*arraysize); // [pix]
            y = (1.0*jj-0.5*arraysize); // [pix]
            x /= dmifscale*(2.0*dmrad/DMsize);
            y /= dmifscale*(2.0*dmrad/DMsize);
            if((fabs(x)<0.5)&&(fabs(y)<0.5))
                data.image[IDif].array.F[jj*arraysize+ii] = 1.0;
            else
                data.image[IDif].array.F[jj*arraysize+ii] = 0.0;
        }
    // convolve dmif
    sig = 0.5*dmifscale*(2.0*dmrad/DMsize);
    gauss_filter("dmif0", "dmif", sig, (long) (2.0*sig));
    delete_image_ID("dmif0");
    IDif = image_ID("dmif");

    save_fits("dmif", "!dmif.fits");

    IDifc = create_image_ID("dmifc", 3, imsize, FLOAT, 0, 0);
    printf("\n");
    for(mx=0; mx<DMsize; mx++)
        for(my=0; my<DMsize; my++)
        {
            printf("\r actuator %2ld %2ld    ", mx, my);
            fflush(stdout);
            // actuator center coordinates
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
    //	save_fits("dmifc","!dmifc.fits");
    printf("\n");

    // INITIALIZE DM CONTROL ARRAY
    dmsizearray = (long*) malloc(sizeof(long)*2);
    dmsizearray[0] = DMsize;
    dmsizearray[1] = DMsize;
    IDdmctrl = create_image_ID("aosimdmctrl", 2, dmsizearray, FLOAT, 1, 0);
    COREMOD_MEMORY_image_set_createsem("aosimdmctrl");

    for(k=0; k<DMsize*DMsize; k++)
        data.image[IDdmctrl].array.F[k] = ran1()*2.0e-7;


    // INITIALIZE TURBULENCE SCREEN
    imsize = (long*) malloc(sizeof(long)*2);
    imsize[0] = arraysize;
    imsize[1] = arraysize;
    IDturb = create_image_ID("WFturb", 2, imsize, FLOAT, 1, 0);
    free(imsize);
    COREMOD_MEMORY_image_set_createsem("WFturb");

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

    // INPUT PUPIL
    IDpupm = make_disk("pupmask", arraysize, arraysize, 0.5*arraysize, 0.5*arraysize, puprad);
    elem = 0;
    optsystsim[0].elemtype[elem] = 1; // pupil mask
    optsystsim[0].elemarrayindex[elem] = IDpupm;
    optsystsim[0].elemZpos[elem] = 0.0;
    elem++;



    // Turbulence screen
    optsystsim[0].elemtype[elem] = 3; // reflective surface
    optsystsim[0].elemarrayindex[elem] = 0; // index
    optsystsim[0].ASPHSURFMarray[0].surfID = IDturb;
    optsystsim[0].elemZpos[elem] = 0.0;
    elem++;

    // DM 0
    optsystsim[0].elemtype[elem] = 3; // reflective surface
    optsystsim[0].elemarrayindex[elem] = 1; // index
    optsystsim[0].ASPHSURFMarray[1].surfID = IDdm0shape;
    optsystsim[0].elemZpos[elem] = 0.0;
    elem++;



    // FOCAL PLANE MASK
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

    optsystsim[0].NBelem = elem;

    optsystsim[0].SAVE = 1;

    // propagate
    OptSystProp_run(optsystsim, 0, 0, optsystsim[0].NBelem, "./conf/");
    ID = image_ID("psfi0");
    imsize = (long*) malloc(sizeof(long)*2);
    imsize[0] = data.image[ID].md[0].size[0];
    imsize[1] = data.image[ID].md[0].size[1];
    imsize[2] = data.image[ID].md[0].size[2];
    IDout = create_image_ID("aosimpsfout", 3, imsize, FLOAT, 1, 0);
    free(imsize);

    COREMOD_MEMORY_image_set_createsem("aosimpsfout");
    data.image[IDout].md[0].write = 1;
    memcpy(data.image[IDout].array.F, data.image[ID].array.F, sizeof(FLOAT)*data.image[ID].md[0].size[0]*data.image[ID].md[0].size[1]*data.image[ID].md[0].size[2]);
    data.image[IDout].md[0].cnt0++;
    data.image[IDout].md[0].write = 0;
    COREMOD_MEMORY_image_set_sempost("aosimpsfout");

    IDarray = (long*) malloc(sizeof(long)*2);
    IDarray[0] = image_ID("WFturb");
    IDarray[1] = image_ID("aosimdmctrl");



    while(1)
    {
        COREMOD_MEMORY_image_set_semwait_OR(IDarray, 2);
        AOsystSim_DMshape("aosimdmctrl", "dmifc", "dmdisp");
        OptSystProp_run(optsystsim, 0, 0, optsystsim[0].NBelem, "./conf/");
        ID = image_ID("psfi0");
        data.image[IDout].md[0].write = 1;
        memcpy(data.image[IDout].array.F, data.image[ID].array.F, sizeof(FLOAT)*data.image[ID].md[0].size[0]*data.image[ID].md[0].size[1]*data.image[ID].md[0].size[2]);
        data.image[IDout].md[0].cnt0++;
        data.image[IDout].md[0].write = 0;
        COREMOD_MEMORY_image_set_sempost("aosimpsfout");
    }


    return(0);
}





/*
     long dmxsize = 50;
    long dmysize = 50;
    long *dmsize;
    long twait = 1; // us
    long long cnt = 0;
    long long cnt0;
    float lambda = 0.7; // um


    long pupsize = 256;
    double puprad = 22.0;
    long IDpupa, IDpupp, IDpuppIR;
    long IDpyrpha, IDpyramp;
    double x, y, r;
    long ii, jj, ii1, jj1;
    double pcoeff = 0.75;
    long pupsize2;
    long ID, IDa, IDp;
    double lenssize = 1200.0;

    long wfssize = 120;
    long *wfssize_array;
    long ID_wfs, ID_wfe;
    long offset, offset1;
    long IDflat, IDdmdisp;

    long IDatmpha;
    long iioffset, jjoffset;
    long IDpuppc, IDpuppcIR;
    long IDpsf, IDpsfIR;
    long *pupsize_array;

    float alphaCorr = 0.8;
    long IDdmt, IDdmtg;

    long moffset;
    long msize;
    long usleeptime = 1000; // delay at each loop

    // Pyramid modulation 
    long modpt;
    long NBmodpt = 64;
    double modr = 0.1;
    double modPA, modx, mody;
    long IDwfscumul;


    long IDpsfnm, IDpsfnmIR;





    pupsize2 = pupsize*pupsize;


    printf("Running fake AO system simulation\n");
    



    dmsize = (long*) malloc(sizeof(long)*2);
    wfssize_array = (long*) malloc(sizeof(long)*2);
    pupsize_array = (long*) malloc(sizeof(long)*2);


    // INITIALIZE

    // create wavefront error input
    dmsize[0] = dmxsize;
    dmsize[1] = dmysize;
    ID_wfe = create_2Dimage_ID("pyrwfeinput", dmxsize, dmysize); // input WF error, to DM sampling

    // create_image_ID("dm_wfe_sim", 2, dmsize, FLOAT, 1, 0);

    // create PSF images
    pupsize_array[0] = pupsize;
    pupsize_array[1] = pupsize;
    IDpsf = create_image_ID("aosimpsf", 2, pupsize_array, FLOAT, 1, 0);
    IDpsfIR = create_image_ID("aosimpsfIR", 2, pupsize_array, FLOAT, 1, 0);

    IDpsfnm = create_image_ID("aosimpsfnm", 2, pupsize_array, FLOAT, 1, 0);
    IDpsfnmIR = create_image_ID("aosimpsfIRnm", 2, pupsize_array, FLOAT, 1, 0);



	
    IDatmpha = read_sharedmem_image("outfiltwf"); // [um]
	if(IDatmpha == -1)
		{
			printf("Running without atmospheric turbulence\n");			
		}



    IDdmdisp = read_sharedmem_image("dmdisp");
    if(IDdmdisp==-1)
		{
			printf("Please create DM channels\n");
			exit(0);
		}

    // create DM
    dmsize[0] = dmxsize;
    dmsize[1] = dmysize;


    // create WFS image
    wfssize_array[0] = wfssize;
    wfssize_array[1] = wfssize;
    wfssize_array[2] = 3; // number of slices in rolling buffer
    ID_wfs = create_image_ID("wfs_sim", 3, wfssize_array, FLOAT, 1, 0);
    data.image[ID_wfs].md[0].cnt1 = 0; // next slice to write

    IDpuppc = create_image_ID("puppcrop", 2, dmsize, FLOAT, 1, 0);



    IDpyrpha = create_2Dimage_ID("pyrpha0", pupsize, pupsize);
    IDpyramp = create_2Dimage_ID("pyramp", pupsize, pupsize);
    for(ii=0; ii<pupsize; ii++)
        for(jj=0; jj<pupsize; jj++)
        {
            x = 1.0*(ii-pupsize/2);
            y = 1.0*(jj-pupsize/2);

            data.image[IDpyrpha].array.F[jj*pupsize+ii] = pcoeff*(fabs(x)+fabs(y));
            if((fabs(x)>lenssize)||(fabs(y)>lenssize))
                data.image[IDpyramp].array.F[jj*pupsize+ii] = 0.0;
            else
                data.image[IDpyramp].array.F[jj*pupsize+ii] = 1.0;
        }
    gauss_filter("pyrpha0", "pyrpha", 15.0, 20);
    IDpyrpha = image_ID("pyrpha");
    save_fits("pyrpha","!pyrpha.fits");

    offset1 = (pupsize-dmxsize)/2;
    printf("OFFSET = %ld\n", offset);
    cnt = -1;
    printf("\n");


    if(ID=image_ID("TpupMask")==-1)
    {
        IDpupa = make_disk("pupa", pupsize, pupsize, 0.5*pupsize, 0.5*pupsize, puprad);
        for(ii=0; ii<pupsize; ii++)
            for(jj=0; jj<pupsize; jj++)
            {
                x = (1.0*ii-0.5*pupsize)/puprad;
                y = (1.0*jj-0.5*pupsize)/puprad;
                r = sqrt(x*x+y*y);
                if(r<0.3)
                    data.image[IDpupa].array.F[jj*pupsize+ii] = 0.0;
                if((fabs(0.5*x+0.5*y)<0.05)||(fabs(0.5*x-0.5*y)<0.05))
                    data.image[IDpupa].array.F[jj*pupsize+ii] = 0.0;
            }
    }
    else
    {
        msize = data.image[ID].md[0].size[0];
        IDpupa = make_disk("pupa", pupsize, pupsize, 0.5*pupsize, 0.5*pupsize, puprad);
        offset = (pupsize-msize)/2;
        for(ii=0; ii<msize; ii++)
            for(jj=0; jj<msize; jj++)
                data.image[IDpupa].array.F[(jj+offset)*pupsize+(ii+offset)] = data.image[ID].array.F[jj*msize+ii];
    }

    //  save_fits("pupa", "!Tpupa.fits");

    IDwfscumul = create_2Dimage_ID("wfscumul", pupsize, pupsize);

    while(1)
    {
        usleep(usleeptime);

        for(ii=0; ii<pupsize2; ii++)
            data.image[IDwfscumul].array.F[ii] = 0.0;

        for(modpt=0; modpt<NBmodpt; modpt++)
        {
            // IMPORT WF ERROR
            iioffset = 6;
            jjoffset = 6;
            for(ii=0; ii<dmxsize; ii++)
                for(jj=0; jj<dmysize; jj++)
                {
                    data.image[ID_wfe].array.F[jj*dmxsize+ii] = 0.0;
                    ii1 = ii+iioffset;
                    jj1 = jj+jjoffset;
                    data.image[ID_wfe].array.F[jj*dmxsize+ii] = data.image[IDatmpha].array.F[jj1*data.image[IDatmpha].md[0].size[0]+ii1]; // um

                }

            data.image[ID_wfe].md[0].cnt0++;
            data.image[ID_wfe].md[0].write = 0;

            IDpupp = create_2Dimage_ID("pupp", pupsize, pupsize);

            // make WFS PSF (without the modulation) 
            data.image[IDpuppc].md[0].write==1;
            for(ii=0; ii<dmxsize; ii++)
                for(jj=0; jj<dmysize; jj++)
                    data.image[IDpupp].array.F[(jj+offset1)*pupsize+(ii+offset1)] = data.image[ID_wfe].array.F[jj*dmxsize+ii]/lambda*M_PI*2.0;

            mk_complex_from_amph("pupa","pupp","pupc");
            permut("pupc");
            do2dfft("pupc","focc");
            delete_image_ID("pupc");
            permut("focc");
            mk_amph_from_complex("focc", "foca", "focp");

            delete_image_ID("focc");
            delete_image_ID("focp");
            IDa = image_ID("foca");


            data.image[IDpsfnm].md[0].write = 1;
            for(ii=0; ii<pupsize2; ii++)
                data.image[IDpsfnm].array.F[ii] = data.image[IDa].array.F[ii]*data.image[IDa].array.F[ii];
            data.image[IDpsfnm].md[0].cnt0++;
            data.image[IDpsfnm].md[0].write = 0;
            delete_image_ID("foca");






            IDpuppIR = create_2Dimage_ID("puppIR", pupsize, pupsize);



            data.image[IDpuppc].md[0].write==1;

            for(ii=0; ii<dmxsize; ii++)
                for(jj=0; jj<dmysize; jj++)
                    data.image[IDpuppc].array.F[jj*dmxsize+ii] = data.image[ID_wfe].array.F[jj*dmxsize+ii] - 2.0*data.image[IDdmdisp].array.F[jj*dmxsize+ii] + modr*cos(2.0*M_PI*modpt/NBmodpt)*(1.0*ii-0.5*dmxsize) + modr*sin(2.0*M_PI*modpt/NBmodpt)*(1.0*jj-0.5*dmysize);   // [um]

            for(ii=0; ii<dmxsize; ii++)
                for(jj=0; jj<dmysize; jj++)
                    data.image[IDpupp].array.F[(jj+offset1)*pupsize+(ii+offset1)] = data.image[IDpuppc].array.F[jj*dmxsize+ii]/lambda*M_PI*2.0;  // [rad]

            data.image[IDpuppc].md[0].cnt0++;
            data.image[IDpuppc].md[0].write = 0;

            for(ii=0; ii<dmxsize; ii++)
                for(jj=0; jj<dmysize; jj++)
                    data.image[IDpuppIR].array.F[(jj+offset1)*pupsize+(ii+offset1)] = (lambda/1.65)*data.image[IDpuppc].array.F[jj*dmxsize+ii]/lambda*M_PI*2.0;




            mk_complex_from_amph("pupa","pupp","pupc");
            permut("pupc");
            do2dfft("pupc","focc");
            delete_image_ID("pupc");
            permut("focc");
            mk_amph_from_complex("focc", "foca", "focp");

            delete_image_ID("focc");
            IDp = image_ID("focp");
            IDa = image_ID("foca");



            data.image[IDpsf].md[0].write = 1;
            for(ii=0; ii<pupsize2; ii++)
            {
                data.image[IDpsf].array.F[ii] = data.image[IDa].array.F[ii]*data.image[IDa].array.F[ii];
                data.image[IDa].array.F[ii] *= data.image[IDpyramp].array.F[ii];
                data.image[IDp].array.F[ii] += data.image[IDpyrpha].array.F[ii];
            }
            data.image[IDpsf].md[0].cnt0++;
            data.image[IDpsf].md[0].write = 0;


            mk_complex_from_amph("pupa","puppIR","pupc");
            permut("pupc");
            do2dfft("pupc","focc");
            delete_image_ID("pupc");
            permut("focc");
            mk_amph_from_complex("focc", "focaIR", "focpIR");

            delete_image_ID("focc");
            delete_image_ID("focpIR");
            IDa = image_ID("focaIR");


            data.image[IDpsfIR].md[0].write = 1;
            for(ii=0; ii<pupsize2; ii++)
                data.image[IDpsfIR].array.F[ii] = data.image[IDa].array.F[ii]*data.image[IDa].array.F[ii];
            data.image[IDpsfIR].md[0].cnt0++;
            data.image[IDpsfIR].md[0].write = 0;



            mk_complex_from_amph("foca","focp","focc1");
            delete_image_ID("foca");
            delete_image_ID("focp");
            permut("focc1");
            do2dfft("focc1","pupc1");
            delete_image_ID("focc1");
            permut("pupc1");
            mk_amph_from_complex("pupc1", "pupa1", "pupp1");
            delete_image_ID("pupc1");
            delete_image_ID("pupp1");

            ID = image_ID("pupa1");
            offset = (pupsize-wfssize)/2;

            for(ii=0; ii<pupsize2; ii++)
                data.image[IDwfscumul].array.F[ii] += data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
            delete_image_ID("pupa1");
        }


        data.image[ID_wfs].md[0].write = 1;
        moffset = data.image[ID_wfs].md[0].cnt1*wfssize*wfssize;
        for(ii1=0; ii1<wfssize; ii1++)
            for(jj1=0; jj1<wfssize; jj1++)
                data.image[ID_wfs].array.F[moffset+jj1*wfssize+ii1] = data.image[IDwfscumul].array.F[(jj1+offset)*pupsize+ii1+offset]; //*data.image[ID].array.F[(jj1+offset)*pupsize+ii1+offset];
        data.image[ID_wfs].md[0].cnt1++;
        if(data.image[ID_wfs].md[0].cnt1==data.image[ID_wfs].md[0].size[2])
            data.image[ID_wfs].md[0].cnt1 = 0;
        data.image[ID_wfs].md[0].cnt0++;
        data.image[ID_wfs].md[0].write = 0;


        printf("\r%10ld       ", data.image[ID_wfs].md[0].cnt0);
        fflush(stdout);
        //  cnt = cnt0;
        //	}
    }

    delete_image_ID("wfscumul");
    free(dmsize);
    free(wfssize_array);
    free(pupsize_array);


    return 0;
}
*/
