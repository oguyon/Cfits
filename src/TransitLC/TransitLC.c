#include <fitsio.h> 
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include "CLIcore.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"

#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "statistic/statistic.h"
#include "TransitLC/TransitLC.h"


extern DATA data;

#define SBUFFERSIZE 2000




// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//



int TransitLC_run_cli()
{

  if(CLI_checkarg(1, 2)==0)
    {
      TransitLC_run(data.cmdargtoken[1].val.numf);
      return 0;
    }
  else
    return 1;
}



int init_TransitLC()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "exoplanet transit light curve analysis");
  data.NBmodule++;
  
  strcpy(data.cmd[data.NBcmd].key,"tlcrun");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = TransitLC_run_cli;
  strcpy(data.cmd[data.NBcmd].info,"Run transit light curve analysis");
  strcpy(data.cmd[data.NBcmd].syntax,"<index>");
  strcpy(data.cmd[data.NBcmd].example,"tlcrun 4");
  strcpy(data.cmd[data.NBcmd].Ccall,"long TransitLC_run(long index)");
  data.NBcmd++;
 
 // add atexit functions here

  return 0;

}


//
// make light curve
//
// input IDplanet stores planets data
// size0 : NB of planet
// line0 : period (sec)
// line1 : start-transit time (sec) - this is the first transit
// line2 : transit duration (sec)
// line3 : transit depth 
//
// light curve is stored as image
// line 1 : time
// line 2 : flux
// line 3 : error
//
long TransitLC_mkLC(char *IDpl_name, char *IDout_name)
{
	long IDout;
	long NBpt;
	long IDpl;
	
	double dt = 10.0; // interval
	double tstart = 0.0;
	double tend = 3600.0*24*50.0; 
	double lctime;
	
	long NBpl;
	long pl;
	double P, t0, w, a;
	double pha;
	double tmpf;
	long ii;
	
	double noise = 1.0e-4;
	
	
	
	NBpt = (tend-tstart)/dt;
	
	printf("Initiating ligt curve  (%ld samples) ... ", NBpt);
	fflush(stdout);
	IDout = create_2Dimage_ID(IDout_name, NBpt, 3);
	for(ii=0;ii<NBpt;ii++)
	{
		lctime = tstart + dt*ii;
		data.image[IDout].array.F[ii] = lctime;
		data.image[IDout].array.F[NBpt*1 + ii] = 1.0 + noise*gauss(); // flux
		data.image[IDout].array.F[NBpt*2 + ii] = 0.1; // measurement error
	}
	printf("\n");
	fflush(stdout);
	
	
	
	
	IDpl = image_ID(IDpl_name);
	printf("Creating planetary transits (%ld planets) ...", data.image[IDpl].md[0].size[1]);
	NBpl = data.image[IDpl].md[0].size[0];
	for(pl=0;pl<NBpl;pl++)
		{
			P = data.image[IDpl].array.F[NBpl*0+pl];
			t0 = data.image[IDpl].array.F[NBpl*1+pl];
			w = data.image[IDpl].array.F[NBpl*2+pl];
			a = data.image[IDpl].array.F[NBpl*3+pl];
			
			for(ii=0;ii<NBpt;ii++)
				{
					lctime = data.image[IDout].array.F[ii];
					pha = modf((lctime-t0)/P+1.0, &tmpf);
					if(pha<w/P)
						data.image[IDout].array.F[NBpt*1 + ii] -= a;
				}
		}
	printf("\n");
	fflush(stdout);
	
	return IDout;
}


//
// output has regular sampling
//

long TransitLC_EdgeDetect(char *ID_name, char *IDout_name, double edgedt, double dt)
{
	long ID, IDout;
	long NBpt, NBpt1;
	long ii, ii1, iistart;
	double lctime, lctime1;
	
	double fluxN, fluxP;
	double fluxNcnt, fluxPcnt;
	double alpha;
	
	double tstart, tend;
	
	
	ID = image_ID(ID_name);
	NBpt = data.image[ID].md[0].size[0];
	tstart = data.image[ID].array.F[0];
	tend = data.image[ID].array.F[NBpt-1];
	
	
	NBpt1 = (tend-tstart)/dt;
	
	IDout = create_2Dimage_ID(IDout_name, NBpt1, 4);
	// line 0 : time
	// line 1 : flux
	// line 2 : convolved flux (edge detect)
	// line 3 : error
	
	iistart = 0;
	printf("\n");
	for(ii1=0; ii1<NBpt1; ii1++)
	{
		printf("\r %6ld / %6ld    ", ii1, NBpt1);
		fflush(stdout);
		
		lctime1 = tstart + dt*ii1;
		data.image[IDout].array.F[ii1] = lctime1;
		
		while(data.image[ID].array.F[iistart] < lctime1-edgedt)
			iistart++;
		
		fluxN = 0.0;
		fluxP = 0.0;
		fluxNcnt = 0.0;
		fluxPcnt = 0.0;
		ii = iistart;
		while((data.image[ID].array.F[ii] < (lctime1+edgedt)) && (ii<NBpt))
			{
				lctime = data.image[ID].array.F[ii];
				alpha = (lctime-lctime1)/edgedt;
				if(alpha<0.0)
				{
					alpha += 1.0;
					if(alpha<0.0)
						alpha = 0.0;
					fluxP += alpha*data.image[ID].array.F[NBpt*1+ii];
					fluxPcnt += alpha;
				}
				else
				{
					alpha = 1.0-alpha;
					if(alpha<0.0)
						alpha = 0.0;
					fluxN += alpha*data.image[ID].array.F[NBpt*1+ii];
					fluxNcnt += alpha;
				}				
				
				ii++;
			}
		fluxP /= fluxPcnt;
		fluxN /= fluxNcnt;
				
		data.image[IDout].array.F[ii1] = lctime1;
		data.image[IDout].array.F[NBpt1+ii1] = (fluxP+fluxN)/2.0;			
		data.image[IDout].array.F[NBpt1*2+ii1] = fluxN-fluxP;
		data.image[IDout].array.F[NBpt1*3+ii1] = 1.0/sqrt(fluxPcnt+fluxNcnt);				
	}
	
	return(IDout);
}



//
// scan for transit events
// 
long TransitLC_scanTE(char *IDin_name, char *IDout_name, double Pmin, double Pmax, double Pstep, double t0step, double wmax, double wstep)
{
	long IDin, IDout;
	long t0size, Psize, wsize;
	double t0, P, w;
	long ii, jj, kk;
	double dt; // time interval in input
	long NBpt;
	double val;
	long valcnt;
	double t1;
	long ii1;
	
	
	// output 3D image size
	t0size = (long) (Pmax/t0step);
	Psize = (long) ((Pmax-Pmin)/Pstep);
	wsize = (long) (wmax/wstep);
	
	IDin = image_ID(IDin_name);
	dt = (data.image[IDin].array.F[11] - data.image[IDin].array.F[1])/10.0;
	NBpt = data.image[IDin].md[0].size[0];
	
	printf("dt = %f sec,  NBpt = %ld\n", dt, NBpt);
	
	
	IDout = create_3Dimage_ID(IDout_name, Psize, t0size, wsize);
	printf("3D scan:   %ld x %ld x %ld\n", Psize, t0size, wsize);
	printf("\n");
		
	
	 for(ii=0;ii<Psize;ii++) // scan P [sec]
	 {
		 P = Pmin + 1.0*ii*Pstep;
		 
		 printf("\rP scan  %6ld / %6ld    %f sec   (%f - %f)     ", ii, Psize, P, Pmin, Pmax);
		fflush(stdout);
		 
		 for(jj=0; jj<t0size; jj++) // scan t0 [sec]
			{
				t0 = 1.0*jj*t0step;
				if(t0<P)
					{
						for(kk=0;kk<wsize;kk++) // scan w [sec]
							{
								w = 1.0*kk*wstep;
								
								val = 0.0;
								valcnt = 0;
								
								t1 = t0;
								ii1 = (long) (t1/dt);
								
								while(ii1<NBpt)									
									{
										val -= data.image[IDin].array.F[NBpt*2+ii1];
										valcnt++;
										t1 += P;
										ii1 = (long) (t1/dt);
									}
								t1 = t0+w;
								ii1 = (long) (t1/dt);
								while(ii1<NBpt)									
									{
										val += data.image[IDin].array.F[NBpt*2+ii1];
										valcnt++;
										t1 += P;
										ii1 = (long) (t1/dt);
									}								
							
								
								
								data.image[IDout].array.F[kk*t0size*Psize+jj*Psize+ii] = val/valcnt;
							}
					}
			}
	}
	
	printf("\n\n");
	
	return(IDout);
}





long TransitLC_run(long index)
{
	long IDpl;
	long pl;
	long NBpl = 1;
	
	// light curve
	long IDlc;
	long lcNBpt;
	
	long IDout1;
	long NBpt, ii, jj, kk;
	double Pmin, Pmax, Pstep;
	double t0step;
	double wmax, wstep;
	long IDout2;
	
	long iimax, jjmax, kkmax;
	double val, valmax, Pval, t0val, wval;
	long xsize, ysize, zsize;
	
	
	FILE *fpout;
	FILE *fp;
	
	
	double pha;
	double tmpf;
	
	
	
	IDpl = create_2Dimage_ID("lcplanets", NBpl, 4);
	
	for(pl = 0; pl < NBpl; pl++)
		{
			data.image[IDpl].array.F[0*NBpl + pl] = 3600.0*24.0*1.0; // period
			data.image[IDpl].array.F[1*NBpl + pl] = 3600.0*12.0; // t0
			data.image[IDpl].array.F[2*NBpl + pl] = 3600.0*1.0; // duration
			data.image[IDpl].array.F[3*NBpl + pl] = 1.0e-5; // depth			
		}
	
	IDlc = image_ID("tlc00");
	if(IDlc == -1)
	{
		IDlc = TransitLC_mkLC("lcplanets", "tlc00");
		save_fits("tlc00", "!tlc00.fits");
	}
	
	fp = fopen("transitcurve.txt", "w");
	lcNBpt = data.image[IDlc].md[0].size[0];
	for(ii=0;ii<lcNBpt;ii++)
		fprintf(fp, "%f %f\n", data.image[IDlc].array.F[ii], data.image[IDlc].array.F[lcNBpt*1+ii]);
	fclose(fp);
	
	
	IDout1 = image_ID("tlc00o");
	if(IDout1==-1)
	{
		printf("Edge detection ...");
		fflush(stdout);
		IDout1 = TransitLC_EdgeDetect("tlc00", "tlc00o", 3000.0, 100.0);
		save_fits("tlc00o", "!tlc00o.fits");
		printf("\n");
		fflush(stdout);
		list_image_ID();
	}
	
	NBpt = data.image[IDout1].md[0].size[0];
	fpout = fopen("outlc.dat", "w");
	for(ii=0;ii<NBpt;ii++)
		fprintf(fpout, "%ld  %18.16f %18.16f %18.16f %18.16f\n", ii, data.image[IDout1].array.F[ii], data.image[IDout1].array.F[NBpt+ii], data.image[IDout1].array.F[2*NBpt+ii], data.image[IDout1].array.F[3*NBpt+ii]);
	fclose(fpout);
	
		
	Pmin = 3600.0*23.5;
	Pmax = 3600.0*24.5;
	Pstep = 2.0;
	
	t0step = 20.0;
	wmax = 3600.0*1.2;
	wstep = 200.0;
	
	IDout2 = TransitLC_scanTE("tlc00o", "tlc3Dscan", Pmin, Pmax, Pstep, t0step, wmax, wstep);
	save_fits("tlc3Dscan", "!tlc3Dscan.fits");
	
	xsize = data.image[IDout2].md[0].size[0];
	ysize = data.image[IDout2].md[0].size[1];
	zsize = data.image[IDout2].md[0].size[2];
	valmax = 0.0;
	for(ii=0;ii<xsize;ii++)
		for(jj=0;jj<ysize;jj++)
			for(kk=1;kk<zsize;kk++)
				{
					val = data.image[IDout2].array.F[kk*xsize*ysize+jj*xsize+ii];
					if(val>valmax)
						{
							valmax = val;
							iimax = ii;
							jjmax = jj;
							kkmax = kk;
						}
				}
				
	Pval = Pmin + Pstep*iimax;
	t0val = t0step*jjmax;
	wval = wstep*kkmax;
	
	
	printf("BEST SOLUTION   %g : %f (%ld/%ld)   %f (%ld/%ld)     %f (%ld/%ld)\n", valmax, Pval, iimax, xsize, t0val, jjmax, ysize, wval, kkmax, zsize);
	
	list_image_ID();
	
	//
	// co-add individual transits
	// put raw data as well as binned in 1/10th transit width
	//
	printf("IDlc = %ld   lcNBpt = %ld\n", IDlc, lcNBpt);
	fflush(stdout);
	fp = fopen("tpts.dat", "w");
	for(ii=0;ii<lcNBpt;ii++)
		{
			pha = (data.image[IDlc].array.F[ii]-t0val)/Pval;
			pha = modf(pha, &tmpf);
			while(pha<0.0)
				pha += 1.0;
			
			
			if((pha<2.0*wval/Pval)||(pha>1.0-wval/Pval))
			{
				if(pha>0.5)
					pha -= 1.0;
				pha -= 0.5*wval/Pval;
				fprintf(fp, "%10.8f  %20.18f  %20.18f\n", pha/(wval/Pval), data.image[IDlc].array.F[ii], data.image[IDlc].array.F[lcNBpt*1+ii]);
			}
				
		}
	fclose(fp);
	
	
	

	
	return(IDlc);
}


