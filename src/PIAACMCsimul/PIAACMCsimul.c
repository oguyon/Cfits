#include <fitsio.h> 
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include "Cfits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate/WFpropagate.h"

#include "PIAACMCsimul/PIAACMCsimul.h"

extern DATA data;

#define SBUFFERSIZE 2000



long PIAACMCSIMUL_ARRAYSIZE = 4096; 
double PIAACMCSIMUL_BEAMRADPIX = 400.0;

long PIAACMCSIMUL_NBLAMBDA = 5;
double PIAACMCSIMUL_LAMBDASTART = 0.4e-6;
double PIAACMCSIMUL_LAMBDAEND = 0.9e-6;
double *lambda_array;


//
// OPTICAL CONFIGURATION
// all Zpos relative to input beam
// 
double beamrad = 0.01; // [m]
double DM1_Zpos = -0.2;
double DM2_Zpos = 0.2;
double PIAAM1_Zpos = 0.5; // [m]
double PIAAM2_Zpos = 1.5; // [m]

double PUPIL_SCALE = 0.0;



// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//

int PIAACMCsimul_run_cli()
{

  //  if(CLI_checkarg(1, 4)+CLI_checkarg(2, 3)+CLI_checkarg(3, 1)+CLI_checkarg(4, 1)+CLI_checkarg(5, 1)==0)
  //  {

  PIAACMCsimul_run();
  //data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numf, data.cmdargtoken[5].val.numf);
  return 0;

  //  }
  //else
  // return 1;
}



int init_PIAACMCsimul()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "PIAACMC system simulation");
  data.NBmodule++;
  
  strcpy(data.cmd[data.NBcmd].key,"piaacmcsimrun");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = PIAACMCsimul_run_cli;
  strcpy(data.cmd[data.NBcmd].info,"Simulate PIAACMC");
  strcpy(data.cmd[data.NBcmd].syntax,"no argument");
  strcpy(data.cmd[data.NBcmd].example,"piaacmcsimrun");
  strcpy(data.cmd[data.NBcmd].Ccall,"int PIAACMCsimul_run()");
  data.NBcmd++;
 



  PIAACMCsimul_init();
    
  // add atexit functions here
  atexit(PIAACMCsimul_free);

  return 0;

}




void PIAACMCsimul_init( void )
{
  long k;

  PUPIL_SCALE = beamrad/PIAACMCSIMUL_BEAMRADPIX; // [m/pix]

  lambda_array = (double*) malloc(sizeof(double)*PIAACMCSIMUL_NBLAMBDA);
  
  for(k=0;k<PIAACMCSIMUL_NBLAMBDA;k++)
    lambda_array[k] = PIAACMCSIMUL_LAMBDASTART + (0.5+k)*(PIAACMCSIMUL_LAMBDAEND-PIAACMCSIMUL_LAMBDASTART)/PIAACMCSIMUL_NBLAMBDA;
}


void PIAACMCsimul_free( void )
{
  free(lambda_array);
}






int PIAACMCsimu_propoagateCube(char *IDin_amp_name, char *IDin_pha_name, char *IDout_amp_name, char *IDout_pha_name, double zprop)
{
  long k;
  long ii;
  long size = PIAACMCSIMUL_ARRAYSIZE;
  long size2;
  long IDin_amp, IDin_pha;
  long IDc_in, IDc_out;
  long IDout_amp, IDout_pha;
  
  double amp, pha, re, im;

  size2 = size*size;

  IDin_amp = image_ID(IDin_amp_name);
  IDin_pha = image_ID(IDin_pha_name);
  IDc_in = create_2DCimage_ID("tmppropCin", size, size);

  IDout_amp = create_3Dimage_ID(IDout_amp_name, size, size, PIAACMCSIMUL_NBLAMBDA);
  IDout_pha = create_3Dimage_ID(IDout_pha_name, size, size, PIAACMCSIMUL_NBLAMBDA);

  for(k=0; k<PIAACMCSIMUL_NBLAMBDA; k++)
    {
      printf("k = %ld / %ld\n", k, PIAACMCSIMUL_NBLAMBDA);
      for(ii=0;ii<size2;ii++)
	{
	  amp = data.image[IDin_amp].array.F[k*size2+ii];
	  pha = data.image[IDin_pha].array.F[k*size2+ii];
	  data.image[IDc_in].array.CF[ii].re = amp*cos(pha);
	  data.image[IDc_in].array.CF[ii].im = amp*sin(pha);
	}
      Fresnel_propagate_wavefront("tmppropCin", "tmppropCout", PUPIL_SCALE, zprop, lambda_array[k]);

      IDc_out = image_ID("tmppropCout");
      for(ii=0;ii<size2;ii++)
	{
	  re = data.image[IDc_out].array.CF[ii].re;
	  im = data.image[IDc_out].array.CF[ii].im;
	  amp = sqrt(re*re+im*im);
	  pha = atan2(im,re);
	  data.image[IDout_amp].array.F[k*size2+ii] = amp;
	  data.image[IDout_pha].array.F[k*size2+ii] = pha;
	}
      delete_image_ID("tmppropCout");
    }

  return 0;
}





int PIAACMCsimul_mkPIAAMshapes_from_RadSag(char *fname, double radius_edge, char *ID_PIAAM1_name, char *ID_PIAAM2_name)
{
  FILE *fp;
  long size = PIAACMCSIMUL_ARRAYSIZE;
  long ii, jj;
  long ID_PIAAM1, ID_PIAAM2;
 
  long nbpt = 20000;
  long k;

  double x, y, r, r1;

  double *r0array;
  double *z0array;
  double *r1array;
  double *z1array;

  double alpha;
  double r00, r01;
  double val;


  r0array = (double*) malloc(sizeof(double)*nbpt);
  z0array = (double*) malloc(sizeof(double)*nbpt);
  r1array = (double*) malloc(sizeof(double)*nbpt);
  z1array = (double*) malloc(sizeof(double)*nbpt);

  fp = fopen(fname, "r");
  for(k=0;k<nbpt;k++)
    fscanf(fp,"%lf %lf %lf %lf\n", &r0array[k], &z0array[k], &r1array[k], &z1array[k]);
  fclose(fp);

  //  for(k=0;k<nbpt;k++)
  //  printf("%ld %.8lf %.8lf %.8lf %.8lf\n", k, r0array[k], z0array[k], r1array[k], z1array[k]);
  

   for(k=0;k<nbpt;k++)
     z1array[k] -= 1.0;
      
  

  ID_PIAAM1 = create_2Dimage_ID(ID_PIAAM1_name, size, size);
  ID_PIAAM2 = create_2Dimage_ID(ID_PIAAM2_name, size, size);
  
  printf("\n\n");
  for(ii=0;ii<size;ii++)
    {
      printf("\r %ld / %ld     ", ii, size);
      fflush(stdout);
      for(jj=0;jj<size;jj++)
	{
	  x = (1.0*ii-0.5*size)/PIAACMCSIMUL_BEAMRADPIX;
	  y = (1.0*jj-0.5*size)/PIAACMCSIMUL_BEAMRADPIX;
	  r = sqrt(x*x+y*y);
	  r1 = r * radius_edge;
	  if(r<2.0)
	    {
	      k = 1;
	      while((r0array[k]<r1)&&(k<nbpt-2))
		k++;
	      r00 = r0array[k-1];
	      r01 = r0array[k];
	      alpha = (r1-r00)/(r01-r00);
	      if(alpha>1.0)
		alpha = 1.0;
	      val = (1.0-alpha)*z0array[k-1] + alpha*z0array[k];
	      data.image[ID_PIAAM1].array.F[jj*size+ii] = val;
	      
	      k = 1;
	      while((r1array[k]<r1)&&(k<nbpt-2))
		k++;
	      r00 = r1array[k-1];
	      r01 = r1array[k];
	      alpha = (r1-r00)/(r01-r00);
	      if(alpha>1.0)
		alpha = 1.0;
	      val = (1.0-alpha)*z1array[k-1] + alpha*z1array[k];
	      data.image[ID_PIAAM2].array.F[jj*size+ii] = val;	
	    }
	}
    }
  printf("\n\n");

  free(r0array);
  free(z0array);
  free(r1array);
  free(z1array);

  return 0;
}





//
// Detailed simulation of PIAACMC
// 
// INPUT
// 
//
// 
int PIAACMCsimul_run()
{
  long IDx, IDy, IDr, IDPA;
  double x, y;
  long IDa, IDp;
  long size = PIAACMCSIMUL_ARRAYSIZE;
  long nblambda = PIAACMCSIMUL_NBLAMBDA;
  long size2;
  long ii, jj, k;
  long IDpiaa1z, IDpiaa2z;


  printf("PIAACMC system simulation\n");

  size2 = size*size;

  // x, y, r and PA coordinates in beam (for convenience)
  IDx = create_2Dimage_ID("xcoord", size, size);
  IDy = create_2Dimage_ID("ycoord", size, size);
  IDr = create_2Dimage_ID("rcoord", size, size);
  IDPA = create_2Dimage_ID("PAcoord", size, size);
  for(ii=0; ii<size; ii++)
    for(jj=0; jj<size; jj++)
      {
	x = (1.0*ii-0.5*size)/PIAACMCSIMUL_BEAMRADPIX;
	y = (1.0*jj-0.5*size)/PIAACMCSIMUL_BEAMRADPIX;
	data.image[IDx].array.F[jj*size+ii] = x;
	data.image[IDy].array.F[jj*size+ii] = y;
	data.image[IDr].array.F[jj*size+ii] = sqrt(x*x+y*y);
	data.image[IDPA].array.F[jj*size+ii] = atan2(y,x);	
      }
  save_fl_fits("xcoord", "!xcoord.fits");
  save_fl_fits("ycoord", "!ycoord.fits");
  save_fl_fits("rcoord", "!rcoord.fits");
  save_fl_fits("PAcoord", "!PAcoord.fits");


  // PIAA shapes
  if(file_exists("piaa1z.fits")==1)
    load_fits("piaa1z.fits", "piaa1z");
  if(file_exists("piaa2z.fits")==1)
    load_fits("piaa2z.fits", "piaa2z");
  IDpiaa1z = image_ID("piaa1z");
  IDpiaa2z = image_ID("piaa2z");
  
  if((IDpiaa1z==-1)||(IDpiaa2z==-1))
    {
      printf("CREATING PIAA SAGS\n");
      if(IDpiaa1z!=-1)
	delete_image_ID("piaa1z");
      if(IDpiaa2z!=-1)
	delete_image_ID("piaa2z");
      PIAACMCsimul_mkPIAAMshapes_from_RadSag("PIAA_Mshapes.txt", beamrad, "piaa1z", "piaa2z");
      save_fl_fits("piaa1z", "!piaa1z.fits");
      save_fl_fits("piaa2z", "!piaa2z.fits");      
      
      IDpiaa1z = image_ID("piaa1z");
      IDpiaa2z = image_ID("piaa2z");
    }
  



  // input pupil
  if(file_exists("pupa0.fits")==1)
    load_fits("pupa0.fits", "pupa0");
  if(file_exists("pupp0.fits")==1)
    load_fits("pupp0.fits", "pupp0");

  IDa = image_ID("pupa0");
  IDp = image_ID("pupp0");
  if((IDa==-1)||(IDp==-1))
    {
      printf("CREATING INPUT PUPIL\n");
      if(IDa!=-1)
	delete_image_ID("pupa0");
      if(IDp!=-1)
	delete_image_ID("pupp0");
      IDa = create_3Dimage_ID("pupa0", size, size, nblambda);
      IDp = create_3Dimage_ID("pupp0", size, size, nblambda);
      
      for(k=0;k<nblambda;k++)
	for(ii=0;ii<size2;ii++)
	  {
	    data.image[IDp].array.F[k*size2+ii] = 0.0;
	    if((data.image[IDr].array.F[ii]>0.05)&&(data.image[IDr].array.F[ii]<1.0))
	      data.image[IDa].array.F[k*size2+ii] = 1.0;
	    else
	      data.image[IDa].array.F[k*size2+ii] = 0.0;
	  }
      save_fl_fits("pupa0", "!pupa0.fits");
      save_fl_fits("pupp0", "!pupp0.fits");
    }

  
  PIAACMCsimu_propoagateCube("pupa0","pupp0", "pupa1", "pupp1", PIAAM1_Zpos); 
  delete_image_ID("pupa0");
  delete_image_ID("pupp0");


  save_fl_fits("pupa1", "!pupa1.fits");
  save_fl_fits("pupp1", "!pupp1.fits");
   

  // ADD PIAA SHAPES
  printf("Adding PIAA shapes ... ");
  fflush(stdout);
  IDp = image_ID("pupp1");
  for(k=0; k<PIAACMCSIMUL_NBLAMBDA; k++)
    {
      for(ii=0;ii<size2;ii++)
	data.image[IDp].array.F[k*size2+ii] -= 2.0 * 2.0*M_PI*data.image[IDpiaa1z].array.F[ii]/lambda_array[k];
    }
  printf("done\n");
  fflush(stdout);

  save_fl_fits("pupp1", "!pupp1b.fits");


  PIAACMCsimu_propoagateCube("pupa1", "pupp1", "pupa2", "pupp2", PIAAM2_Zpos-PIAAM1_Zpos); 
  delete_image_ID("pupa1");
  delete_image_ID("pupp1");
  save_fl_fits("pupa2", "!pupa2.fits");
  save_fl_fits("pupp2", "!pupp2.fits");

  IDp = image_ID("pupp2");
  for(k=0; k<PIAACMCSIMUL_NBLAMBDA; k++)
    {
      for(ii=0;ii<size2;ii++)
	data.image[IDp].array.F[k*size2+ii] += 2.0 * 2.0*M_PI*data.image[IDpiaa2z].array.F[ii]/lambda_array[k];
    }

  save_fl_fits("pupp2", "!pupp2b.fits");



  list_image_ID();  

  return(0);
}

