#include <fitsio.h> 
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "info/info.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate/WFpropagate.h"
#include "statistic/statistic.h"
#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "PIAACMCsimul/PIAACMCsimul.h"
#include "OptSystProp/OptSystProp.h"



extern DATA data;

#define SBUFFERSIZE 2000

OPTSYST *optsyst;
int optsystinit = 0;
long IDx, IDy, IDr, IDPA;


double LAMBDASTART = 0.7e-6;
double LAMBDAEND = 0.9e-6;

MIRRORPIAACMCDESIGN *piaacmc;
OPTSYST *optsyst;


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
 
   
  // add atexit functions here
  atexit(PIAACMCsimul_free);

  return 0;

}


void PIAACMCsimul_free( void )
{
  if(optsystinit ==1)
    {
      free(optsyst);
    }
}


























// ************************************************************************************
//
//           FOCAL PLANE MASK
//
// ************************************************************************************




long PIAACMCsimul_mkFPM_zonemap(char *IDname)
{
  long NBzones;
  long ID;
  double x, y, r;
  long ii, jj;
  long zi;
  long *sizearray;
  
  sizearray = (long*) malloc(sizeof(long)*2);
  sizearray[0] = piaacmc[0].fpmarraysize;
  sizearray[1] = piaacmc[0].fpmarraysize;
  ID = create_image_ID(IDname, 2, sizearray, LONG, 0, 0);
  free(sizearray);

  //  ID = create_2Dimage_ID(IDname, piaacmc[0].fpmarraysize, piaacmc[0].fpmarraysize);
  for(ii=0;ii<piaacmc[0].fpmarraysize;ii++)
    for(jj=0;jj<piaacmc[0].fpmarraysize;jj++)
      {
	x = (2.0*ii-1.0*piaacmc[0].fpmarraysize)/piaacmc[0].fpmarraysize;
	y = (2.0*jj-1.0*piaacmc[0].fpmarraysize)/piaacmc[0].fpmarraysize;
	r = sqrt(x*x+y*y); 
	zi = (long) ceil((1.0-r)*piaacmc[0].NBrings);
	if(zi<0.1)
	  zi = 0;
	if(zi>piaacmc[0].NBrings)
	  zi = piaacmc[0].NBrings;
	data.image[ID].array.F[jj*piaacmc[0].fpmarraysize+ii] = zi;
      }
  piaacmc[0].focmNBzone = piaacmc[0].NBrings;

  return ID;
}





// makes 1-fpm CA
// if mode = -1, make whole 1-fpm
// if mode = zone, make only 1 zone 
long PIAACMCsimul_mkFocalPlaneMask(char *IDzonemap_name, char *ID_name, int mode)
{
  long ID;
  long IDz;
  long size;
  long nblambda;
  long k;
  long ii, jj;
  double x, y, r; // in meter
  long ii1, jj1;
  double fpscale; // [m/pix]
  long zi;
  double t, re, im, amp, pha;
  long size2;

  size = optsyst[0].size;
  size2 = size*size;
  nblambda = optsyst[0].nblambda;

  IDz = image_ID(IDzonemap_name);
  ID = create_3DCimage_ID(ID_name, size, size, nblambda);
 
  // CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/CORONAGRAPHS_ARRAYSIZE

  for(k=0;k<nblambda;k++)
    {
      fpscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size/piaacmc[0].fpzfactor*optsyst[0].lambdaarray[k]*piaacmc[0].Fratio;
      printf("LAMBDA = %.5g m    SCALE = %.5g m/pix \n", optsyst[0].lambdaarray[k], fpscale);
      

      for(ii=0;ii<size;ii++)
	for(jj=0;jj<size;jj++)
	  {
	    x = (1.0*ii-size/2)*fpscale; // [m]
	    y = (1.0*jj-size/2)*fpscale; // [m]	    

	    ii1 = (long) ( (0.5 + 0.5*x/piaacmc[0].fpmRad)*piaacmc[0].fpmarraysize);
	    jj1 = (long) ( (0.5 + 0.5*y/piaacmc[0].fpmRad)*piaacmc[0].fpmarraysize);
	    if((ii1>-1)&&(ii1<piaacmc[0].fpmarraysize)&&(jj1>-1)&&(jj1<piaacmc[0].fpmarraysize))
	      {
		zi = data.image[IDz].array.F[jj1*piaacmc[0].fpmarraysize+ii1];
		t = data.image[piaacmc[0].zonezID].array.F[zi-1];
	      }
	    else
	      {
		zi = 0;
		t = 0.0;
	      }


	    if(mode == -1)
	      {
		if(zi>0.1)
		  {
		    amp = 1.0;
		    pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[k]);
		    re = amp*cos(pha);
		    im = amp*sin(pha);
		    data.image[ID].array.CF[k*size2+jj*size+ii].re = 1.0-re;
		    data.image[ID].array.CF[k*size2+jj*size+ii].im = -im;
		  }
		else
		  {
		    data.image[ID].array.CF[k*size2+jj*size+ii].re = 0.0;
		    data.image[ID].array.CF[k*size2+jj*size+ii].im = 0.0;
		  }
	      }
	    else if (mode == zi)
	      {
		amp = 1.0;
		pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[k]);
		re = amp*cos(pha);
		im = amp*sin(pha);
		data.image[ID].array.CF[k*size2+jj*size+ii].re = 1.0-re;
		data.image[ID].array.CF[k*size2+jj*size+ii].im = -im;
	      }

	  }
    }


  return(ID);
}















//
// initializes the optsyst structure to simulate reflective PIAACMC system
//
void PIAACMCsimul_init( MIRRORPIAACMCDESIGN *design, long index, double TTxld, double TTyld )
{
  long k;
  long size;
  double x, y, PA;
  long ii, jj;
  long nblambda;
  long size2;
  double beamradpix;
  long kx, ky, kxy;
  long IDpiaaz0, IDpiaaz1;
  long surf;
  long IDa;
  char fname_pupa0[200];
  long ID;
  long elem;



  optsyst = (OPTSYST*) malloc(sizeof(OPTSYST));
  optsyst[0].nblambda = 3;
  nblambda = optsyst[0].nblambda;
  for(k=0;k<optsyst[0].nblambda;k++)
    optsyst[0].lambdaarray[k] = LAMBDASTART + (0.5+k)*(LAMBDAEND-LAMBDASTART)/optsyst[0].nblambda;


  optsyst[0].beamrad = design[index].beamrad; // 8mm
  optsyst[0].size = design[index].size;
  size = optsyst[0].size;
  size2 = size*size;
  optsyst[0].pixscale = design[index].pixscale;
  optsyst[0].DFTgridpad = 2; // 0 for full DFT sampling, >0 for faster execution

  beamradpix = optsyst[0].beamrad/optsyst[0].pixscale;
  printf("BEAM RADIUS = %f pix\n", beamradpix);



  // define optical elements and locations
  optsyst[0].NB_DM = 0;
  optsyst[0].NB_asphsurfm = 2;
  optsyst[0].NB_asphsurfr = 0;
  
  optsyst[0].NBelem = 9; 


  elem = 0;
  // ------------------- elem 0: input pupil -----------------------
  optsyst[0].elemtype[elem] = 1; // pupil mask
    // input pupil
  sprintf(fname_pupa0, "pupa0_%ld.fits", size);
  if(file_exists(fname_pupa0)==1)
    load_fits(fname_pupa0, "pupa0");
  IDa = image_ID("pupa0");
  if(IDa==-1)
    {
      printf("CREATING INPUT PUPIL\n");
      if(IDa!=-1)
	delete_image_ID("pupa0");
      IDa = create_3Dimage_ID("pupa0", size, size, nblambda);
            
      for(k=0;k<nblambda;k++)
	for(ii=0;ii<size2;ii++)
	  {
	    if((data.image[IDr].array.F[ii]>0.3)&&(data.image[IDr].array.F[ii]<1.0))
	      data.image[IDa].array.F[k*size2+ii] = 1.0;
	    else
	      data.image[IDa].array.F[k*size2+ii] = 0.0;
	  }     
      sprintf(fname_pupa0, "!pupa0_%ld.fits", size);
      save_fl_fits("pupa0", fname_pupa0);      
    }
  optsyst[0].elemarrayindex[elem] = IDa;  
  optsyst[0].elemZpos[elem] = 0.0;



  elem++;





  // pointing (simulated as mirror)
  ID = create_2Dimage_ID("TTm", size, size);

  for(ii=0;ii<size;ii++)
    for(jj=0;jj<size;jj++)
      {
	x = (1.0*ii-0.5*size)/beamradpix;
	y = (1.0*jj-0.5*size)/beamradpix;
	data.image[ID].array.F[jj*size+ii] = 0.25*(TTxld*x+TTyld*y)*(LAMBDAEND+LAMBDASTART)*0.5;
      }
  save_fits("TTm","!TTm.fits");
  optsyst[0].elemtype[elem] = 3; // reflective mirror   
  optsyst[0].elemarrayindex[elem] = 0; // index
  optsyst[0].ASPHSURFMarray[0].surfID = ID; 
  optsyst[0].elemZpos[elem] = 0.0;
  elem++;

 
  IDpiaaz0 = image_ID("piaa0z");
  IDpiaaz1 = image_ID("piaa1z");


  // ------------------- elem 1: reflective PIAA M0  -----------------------  
  optsyst[0].elemtype[elem] = 3; // reflective PIAA M0   
  optsyst[0].elemarrayindex[elem] = 1; // index
  optsyst[0].ASPHSURFMarray[1].surfID = IDpiaaz0; 
  optsyst[0].elemZpos[elem] = 0.0;
  elem++;

  // ------------------- elem 2: reflective PIAA M1  -----------------------  
  optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
  optsyst[0].elemarrayindex[elem] = 2;
  optsyst[0].ASPHSURFMarray[2].surfID = IDpiaaz1; 
  optsyst[0].elemZpos[elem] = piaacmc[index].piaasep;
  elem++;


  // ------------------- elem 3: opaque mask at reflective PIAA M1  -----------------------  
  optsyst[0].elemtype[elem] = 1; // opaque mask
  ID = make_disk("piaam1mask", size, size, 0.5*size, 0.5*size, piaacmc[0].r1lim*beamradpix);
  optsyst[0].elemarrayindex[elem] = ID;
  optsyst[0].elemZpos[elem] = piaacmc[index].piaasep;
  save_fits("piaam1mask", "!piaam1mask.fits");
  elem++;
  

  // --------------------  elem 4: focal plane mask ------------------------
  optsyst[0].elemtype[elem] = 5; // focal plane mask 
  optsyst[0].elemarrayindex[elem] = 0;
  optsyst[0]. FOCMASKarray[0].fpmID = PIAACMCsimul_mkFocalPlaneMask("fpmzmap", "piaacmcfpm", -1); // this is 1-fpm
  if(1)// testing
    {
      mk_amph_from_complex("piaacmcfpm", "fpma", "fpmp");
      save_fits("fpma", "!fpma.fits");
      save_fits("fpmp", "!fpmp.fits");
      delete_image_ID("fpma");
      delete_image_ID("fpmp");
      //     exit(0);
    }
  optsyst[0]. FOCMASKarray[0].zfactor = piaacmc[0].fpzfactor;
  optsyst[0].elemZpos[elem] = piaacmc[index].piaasep; // plane from which FT is done
  elem++;

  
  // --------------------  elem 5: Lyot mask ------------------------
  optsyst[0].elemtype[elem] = 1; // Lyot mask 
  printf("CREATING LYOT MASK\n");
  IDa = create_3Dimage_ID("lmask", size, size, nblambda);
  
  for(k=0;k<nblambda;k++)
    for(ii=0;ii<size2;ii++)
      {
	if((data.image[IDr].array.F[ii]>0.21)&&(data.image[IDr].array.F[ii]<0.99))
	  data.image[IDa].array.F[k*size2+ii] = 1.0;
	else
	  data.image[IDa].array.F[k*size2+ii] = 0.0;
	  }     
  sprintf(fname_pupa0, "!LyotMask_%ld.fits", size);
  save_fl_fits("lmask", fname_pupa0);      

  optsyst[0].elemarrayindex[elem] = IDa;  
  optsyst[0].elemZpos[elem] = 0.0;
  elem++;


  // --------------------  elem 6: inv PIAA1 ------------------------
  optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
  optsyst[0].elemarrayindex[elem] = 2;
  optsyst[0].ASPHSURFMarray[2].surfID = IDpiaaz1; 
  optsyst[0].elemZpos[elem] = 0.0;
  elem++;

  // --------------------  elem 7: inv PIAA0 ------------------------
  optsyst[0].elemtype[elem] = 3; // reflective PIAA M0
  optsyst[0].elemarrayindex[elem] = 1;
  optsyst[0].ASPHSURFMarray[1].surfID = IDpiaaz0; 
  optsyst[0].elemZpos[elem] = piaacmc[index].piaasep;
  elem++;

  optsystinit = 1;
}











//
// RADIAL PIAACMC SYSTEM DESIGN (geometrical optics)
//



//
// load and fit radial apodization profile 
// modal basis is mk(r) : cos(r*k*M_PI/1.3) 
//
int PIAACMCsimul_load2DRadialApodization(char *IDapo_name, float beamradpix, float centralObs, char *IDapofit_name)
{
  long NBpts;
  long IDm;
  long sizem;
  long kmax = 10;
  long ID, IDmask, IDin;
  long ii, jj;
  long offset;
  long sizein;
  float eps = 1.0e-4;

  sizem = (long) (beamradpix*2);

  // CREATE MODES IF THEY DO NOT EXIST
  if((IDm=image_ID("APOmodesCos"))==-1)
    {
      IDm = linopt_imtools_makeCosRadModes("APOmodesCos", sizem, kmax, ApoFitCosFact*beamradpix, 1.0);  
      save_fits("APOmodesCos", "!APOmodesCos.fits");
    }
  
  // CREATE MASK AND CROP INPUT
  IDmask = create_2Dimage_ID("fitmaskapo", sizem, sizem);
  
  IDin = image_ID(IDapo_name);
  sizein = data.image[IDin].md[0].size[0];
  ID = create_2Dimage_ID("_apoincrop", sizem, sizem);
  offset = (sizein-sizem)/2;
  for(ii=0;ii<sizem;ii++)
    for(jj=0;jj<sizem;jj++)
      {
	data.image[ID].array.F[jj*sizem+ii] = data.image[IDin].array.F[(jj+offset)*sizein+(ii+offset)];
	if((data.image[ID].array.F[jj*sizem+ii]>eps)&&(ii%1==0)&&(jj%1==0))
	  data.image[IDmask].array.F[jj*sizem+ii] = 1.0;
      }


  save_fits("_apoincrop", "!_apoincrop.fits");
  save_fits("fitmaskapo", "!fitmaskapo.fits");
  linopt_imtools_image_fitModes("_apoincrop", "APOmodesCos", "fitmaskapo", 1.0e-8, IDapofit_name, 0);


  if(1) // test fit quality
    {
      linopt_imtools_image_construct("APOmodesCos", IDapofit_name, "testapofitsol");
      save_fits("testapofitsol", "!testapofitsol.fits");
      arith_image_sub("_apoincrop", "testapofitsol", "apofitres");
      arith_image_mult("apofitres", "fitmaskapo", "apofitresm");
      save_fits("apofitres", "!apofitres.fits");
      save_fits("apofitresm", "!apofitresm.fits");
      // linopt_imtools_image_fitModes("apofitres", "APOmodesCos", "fitmaskapo", 1.0e-5, "test2c", 0);
      //list_image_ID();
      info_image_stats("apofitresm", "");
    }

  delete_image_ID("_apoincrop");
  delete_image_ID("fitmaskapo");


  return 0;
}







//
// this function only works for circular PIAA
// uses radial PIAACMC design to initialize PIAA optics shapes and focal plane mask
//
int PIAACMCsimul_init_geomPIAA_rad(char *IDapofit_name)
{
  long i, ii, k;
  double *pup0;
  double *pup1;
  double *flux0cumul;
  double *flux1cumul;

  long IDcoeff;
  long nbcoeff;
  double r;
  FILE *fp;
  double total;

  // to convert r ro r1 (assymptotic outer radius on pup1)
  double coeffa = 3.0; // convergence rate from linear to assymptotic value
  double coeffa1 = 0.5; // convergence radius limit (added to 1.0)

  double r0, r1;
  double FLUX0_in = 0.0; // inside central obstruction
  double FLUX0_out = 0.0; // outside beam edge
  double FLUX1_in = 0.0; // inside central obstruction
  double FLUX1_out = 0.0; // outside beam edge
  double normcoeff;

  // inner profile adjustment
  double a, b, verr, bstep, x, eps1;
  double t0, t0cnt, value;
  int dir, odir;
  long NBistep;
  long iioffset;

  double fluxdens, F0, F1, F2;
  double dr0, dr1, ndr0, ndr1;
  double *piaar0;
  double *piaar1;

  long NBpoints;
  double *piaar00;
  double *piaar11;
  double *piaar01;
  double *piaar10;
  long cnt;
  double tmp;
  double epsilon = 0.000000000001;

  double *piaaM0z;
  double *piaaM1z;
  double r0c, r1c, dx, dy, dist, y3, r0n, slope, dz;


  pup0 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
  pup1 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
  flux0cumul = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
  flux1cumul = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);


  // STEP 1: CREATE AMPLITUDE AND CUMULATIVE INTENSITY PROFILES


  // CREATE OUTPUT AMPLITUDE APODIZATION PROFILE AND ITS CUMUL

  IDcoeff = image_ID(IDapofit_name);
  nbcoeff = data.image[IDcoeff].md[0].size[0];
  printf("%ld coefficients\n", nbcoeff);

  total = 0.0;
  for(ii=0;ii<piaacmc[0].NBradpts;ii++)
    {
      pup1[ii] = 0.0;
      r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
      if(r<1.0)
	r1 = r;
      else
	r1 = 1.0 + (r-1) / pow((1.0 + pow(1.0/coeffa1 * (r-1),coeffa)), 1.0/coeffa);

      for(k=0;k<nbcoeff;k++)
	pup1[ii] += data.image[IDcoeff].array.F[k]*cos(r1*k*M_PI/ApoFitCosFact);      
      if(r<piaacmc[0].centObs1)
	FLUX1_in += pup1[ii]*pup1[ii]*r;
      if(r>1.0)
	FLUX1_out += pup1[ii]*pup1[ii]*r;
      total += pup1[ii]*pup1[ii]*r;
      flux1cumul[ii] = total;
    }


  normcoeff = 1.0/(total-FLUX1_in-FLUX1_out);

  FLUX1_in *= normcoeff;
  FLUX1_out *= normcoeff;
  for(ii=0;ii<piaacmc[0].NBradpts;ii++)
    {
      r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
      flux1cumul[ii] *= normcoeff;
    }


  printf("outer fluxes 1: %lf %lf\n", FLUX1_in, FLUX1_out);
  

 

 // CREATE FLUX0 
  
  total = 0.0;
  for(ii=0;ii<piaacmc[0].NBradpts;ii++)
    {
      r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
      pup0[ii] = 1.0;
  
      if(r<piaacmc[0].centObs0)
	FLUX0_in += pup0[ii]*pup0[ii]*r;
      if(r>1.0)
	FLUX0_out += pup0[ii]*pup0[ii]*r;
      
      total += pup0[ii]*pup0[ii]*r;
      flux0cumul[ii] = total;      
    }
  normcoeff = 1.0/(total-FLUX0_in-FLUX0_out);

  FLUX0_in *= normcoeff;
  FLUX0_out *= normcoeff;

  printf("outer fluxes 0: %lf (%lf)    %lf\n", FLUX0_in, FLUX1_in, FLUX0_out);


  //
  // Compute inner pseudo profile 
  //
  b = 0.5;
  bstep = 0.1;
  verr = 1.0;
  NBistep = piaacmc[0].centObs0*piaacmc[0].NBradpts/piaacmc[0].r0lim;
  //  innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);

  while(fabs(verr)>1.0e-9)
    {
      t0 = 0.0;
      t0cnt = 0.0;
      for(ii=0;ii<NBistep;ii++)
	{
	  x = 1.0*ii/NBistep;
	  r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
	  a = 0.5;
	  eps1 = 1e-8;
	  if(x<eps1)
	    x = eps1;
	  if(x>1.0-eps1)
	    x = 1.0-eps1;
	  pup0[ii] =  b + (1.0-b)*(0.5+atan(-a/b/(x*x) + a/pow(x-1.0,2))/M_PI);
	  t0 += r*pup0[ii]* pup0[ii];
	  flux0cumul[ii] = t0;
	}
      
      verr = t0*normcoeff - FLUX1_in;
      
      odir = dir;
      if(verr>0.0) // too much light
	{
	  b /= (1.0+bstep);
	  dir = -1.0;
	}
      else
	{
	  b *= (1.0+bstep);
	  dir = 1.0;
	}
      if(odir*dir<0.0)
	bstep *= 0.1;
      printf(".");
      fflush(stdout);
    }
  printf("\n");
  printf("TOTAL = %f -> %g (%g %g)\n", b, t0*normcoeff, bstep, verr);


  // outer region
  b = 0.5;
  bstep = 0.1;
  verr = 1.0;
  NBistep = piaacmc[0].NBradpts*(piaacmc[0].r0lim-1.0)/piaacmc[0].r0lim;
  //  innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);
  iioffset = (long) (1.0*piaacmc[0].NBradpts/piaacmc[0].r0lim);
  NBistep = piaacmc[0].NBradpts-iioffset;
  while(fabs(verr)>1.0e-9)
    {
      t0 = 0.0;
      t0cnt = 0.0;
      for(ii=0;ii<NBistep;ii++)
	{
	  x = 1.0-1.0*ii/NBistep;
	  r = 1.0+1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
	  a = 0.5;
	  eps1 = 1e-8;
	  if(x<eps1)
	    x = eps1;
	  if(x>1.0-eps1)
	    x = 1.0-eps1;
	  pup0[ii+iioffset] =  b + (1.0-b)*(0.5+atan(-a/b/(x*x) + a/pow(x-1.0,2))/M_PI);
	  t0 += r*pup0[ii+iioffset]* pup0[ii+iioffset];
	  flux0cumul[ii+iioffset] = t0;
	}
      
      verr = t0*normcoeff - FLUX1_out;
      
      odir = dir;
      if(verr>0.0) // too much light
	{
	  b /= (1.0+bstep);
	  dir = -1.0;
	}
      else
	{
	  b *= (1.0+bstep);
	  dir = 1.0;
	}
      if(odir*dir<0.0)
	bstep *= 0.1;
      printf(".");
      fflush(stdout);
    }
  printf("\n");
  printf("TOTAL = %f -> %g (%g %g)\n", b, t0*normcoeff, bstep, verr);



  total = 0.0;
  FLUX0_in = 0.0;
  FLUX0_out = 0.0;
  for(ii=0;ii<piaacmc[0].NBradpts;ii++)
    {
      r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
      if(r<piaacmc[0].centObs0)
	FLUX0_in += pup0[ii]*pup0[ii]*r;
      if(r>1.0)
	FLUX0_out += pup0[ii]*pup0[ii]*r;
      
      total += pup0[ii]*pup0[ii]*r;
      flux0cumul[ii] = total;      
      flux0cumul[ii] *= normcoeff;;
    }
  FLUX0_in *= normcoeff;
  FLUX0_out *= normcoeff;
  
  printf("outer fluxes 0: %lf (%lf)    %lf\n", FLUX0_in, FLUX1_in, FLUX0_out);


  


  fp = fopen("pup01.prof", "w");
  for(ii=0;ii<piaacmc[0].NBradpts;ii++)
    {
      r0 = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
      r1 = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
      fprintf(fp, "%f %f %g %g %g %g\n", r0, r1, pup0[ii], pup1[ii], flux0cumul[ii], flux1cumul[ii]);
    }
  fclose(fp);

  

   
  // STEP 2: COMPUTE r0 - r1 CORRESPONDANCE

  piaar00 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r0 as a function of r0 index
  piaar11 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r1 as a function of r1 index
  piaar10 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r1 as a function of r0 index
  piaar01 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r0 as a function of r1 index

  /* computing r0 and r1 */
  /* r0 and r1 are dimensionless */

  /* first, r0 is evenly distributed on the first optic */
  for(i=0;i<piaacmc[0].NBradpts;i++)
    {
      piaar00[i] = piaacmc[0].r0lim*i/piaacmc[0].NBradpts;
      piaar11[i] = piaacmc[0].r1lim*i/piaacmc[0].NBradpts;      
    }

  i=0;
  ii=0;  
  cnt = 0;
  piaar00[0] = 0.0;
  piaar10[0] = 0.0;
  fp = fopen("test0.txt", "w");
  for(i=1;i<piaacmc[0].NBradpts;i++)
    {
      F0 = flux0cumul[i];
      while((flux1cumul[ii]<flux0cumul[i])&&(ii<piaacmc[0].NBradpts))
	ii++;
      F1 = flux1cumul[ii-1];
      F2 = flux1cumul[ii];
      
      /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
      if(fabs(F2-F1)>0.0000000001)
	fluxdens = (F2-F1)/(2.0*ii-1.0);
      else
	fluxdens = 0.0000000001;
      x = sqrt((F0-F1)/fluxdens+(1.0*ii*ii-2.0*ii+1.0)) + 1.0 - 1.0*ii;

      piaar10[i] = piaacmc[0].r1lim*(1.0*ii-1.0+x)/piaacmc[0].NBradpts;
      fprintf(fp, "%lf %lf %lf\n", piaar00[i], piaar10[i], F0);
    }
  fclose(fp);
  


  i=0;
  ii=0;  
  cnt = 0;
  piaar01[0] = 0.0;
  piaar11[0] = 0.0;
  fp = fopen("test1.txt", "w");
  for(i=1;i<piaacmc[0].NBradpts;i++)
    {
      F0 = flux1cumul[i];
      while((flux0cumul[ii]<flux1cumul[i])&&(ii<piaacmc[0].NBradpts))
	ii++;
      F1 = flux0cumul[ii-1];
      F2 = flux0cumul[ii];
      
      /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
      if(fabs(F2-F1)>0.0000000001)
	fluxdens = (F2-F1)/(2.0*ii-1.0);
      else
	fluxdens = 0.0000000001;
      x = sqrt((F0-F1)/fluxdens+(1.0*ii*ii-2.0*ii+1.0)) + 1.0 - 1.0*ii;

      piaar01[i] = piaacmc[0].r0lim*(1.0*ii-1.0+x)/piaacmc[0].NBradpts;
      fprintf(fp, "%lf %lf %lf\n", piaar11[i], piaar01[i], F0);
    }
  fclose(fp);
  




  printf("======== Compute PIAA mirror shapes ============\n");
  piaaM0z = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
  piaaM1z = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
  
  piaaM0z[0] = 0.0;
  piaaM1z[0] = piaacmc[0].piaasep;


  for(i=0;i<piaacmc[0].NBradpts-1;i++)
    {
      r0c = piaar00[i];
      r1c = piaar10[i];
      dx = (r0c-r1c)*piaacmc[0].beamrad;
      dz = piaaM1z[i]-piaaM0z[i];
      dist = sqrt(dx*dx+dz*dz);
      y3 = dist - dz;
      if(fabs(dx)>0.000000001)
	slope = y3/dx;
      else
	slope = 0.0;
      r0n = piaacmc[0].r0lim*(i+1)/piaacmc[0].NBradpts;
      piaaM0z[i+1] = piaaM0z[i] + slope*(r0n-r0c)*piaacmc[0].beamrad;
      
      if(fabs(dx)>0.000000001)
	slope = y3/dx;
      else
	slope = 0.0;      
      piaaM1z[i+1] = piaaM1z[i] + slope*(piaar10[i+1]-r1c)*piaacmc[0].beamrad;
    }      

  fp = fopen("PIAA_Mshapes.txt", "w");
  for(ii=0;ii<piaacmc[0].NBradpts;ii++)
    fprintf(fp, "%18.16f %18.16f %18.16f %18.16f\n", piaar00[ii]*piaacmc[0].beamrad, piaaM0z[ii], piaar10[ii]*piaacmc[0].beamrad, piaaM1z[ii]);
  fclose(fp);



  
  free(flux0cumul);
  free(flux1cumul);
  free(pup0);
  free(pup1);

  free(piaar0);
  free(piaar1);
  free(piaar00);
  free(piaar10);
  free(piaar01);
  free(piaar11);
  

  return(0);
}


















//
// make PIAA shapes from radial sag profile
//
int PIAACMCsimul_mkPIAAMshapes_from_RadSag(char *fname, char *ID_PIAAM0_name, char *ID_PIAAM1_name)
{
  FILE *fp;
  long size;
  long ii, jj;
  long ID_PIAAM0, ID_PIAAM1;
 
   long k;

  double x, y, r, r1;

  double *r0array;
  double *z0array;
  double *r1array;
  double *z1array;

  double alpha;
  double r00, r01;
  double val;

  double beamradpix;
  int ret;


  size = piaacmc[0].size;
  beamradpix = piaacmc[0].beamrad/piaacmc[0].pixscale;
  printf("SIZE = %ld, beamrad = %f pix, sep = %f m\n", size, beamradpix, piaacmc[0].piaasep);
  fflush(stdout);
  
 

  r0array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
  z0array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
  r1array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
  z1array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);

  fp = fopen(fname, "r");
  for(k=0;k<piaacmc[0].NBradpts;k++)
    ret = fscanf(fp,"%lf %lf %lf %lf\n", &r0array[k], &z0array[k], &r1array[k], &z1array[k]);
  fclose(fp);

  //  for(k=0;k<nbpt;k++)
  //  printf("%ld %.8lf %.8lf %.8lf %.8lf\n", k, r0array[k], z0array[k], r1array[k], z1array[k]);
   

  for(k=0;k<piaacmc[0].NBradpts;k++)
    z1array[k] -= piaacmc[0].piaasep;
  



  ID_PIAAM0 = create_2Dimage_ID(ID_PIAAM0_name, size, size);
  ID_PIAAM1 = create_2Dimage_ID(ID_PIAAM1_name, size, size);
  
  printf("\n\n");
  for(ii=0;ii<size;ii++)
    {
      printf("\r %ld / %ld     ", ii, size);
      fflush(stdout);
      for(jj=0;jj<size;jj++)
	{
	  x = (1.0*ii-0.5*size)/beamradpix;
	  y = (1.0*jj-0.5*size)/beamradpix;
	  r = sqrt(x*x+y*y)*piaacmc[0].beamrad;

	  if(r<piaacmc[0].r0lim*piaacmc[0].beamrad)
	    {
	      k = 1;
	      while((r0array[k]<r)&&(k<piaacmc[0].NBradpts-2))
		k++;
	      r00 = r0array[k-1];
	      r01 = r0array[k];
	      alpha = (r-r00)/(r01-r00);
	      if(alpha>1.0)
		alpha = 1.0;
	      val = (1.0-alpha)*z0array[k-1] + alpha*z0array[k];
	      data.image[ID_PIAAM0].array.F[jj*size+ii] = val;
	    }
	  else
	    data.image[ID_PIAAM0].array.F[jj*size+ii] = 0.0;

	  if(r<piaacmc[0].r1lim*piaacmc[0].beamrad)
	    {
	      k = 1;
	      while((r1array[k]<r)&&(k<piaacmc[0].NBradpts-2))
		k++;
	      r00 = r1array[k-1];
	      r01 = r1array[k];
	      alpha = (r-r00)/(r01-r00);
	      if(alpha>1.0)
		alpha = 1.0;
	      val = (1.0-alpha)*z1array[k-1] + alpha*z1array[k];
	      data.image[ID_PIAAM1].array.F[jj*size+ii] = -val;//-piaacmc[0].piaasep);	
	    }
	  else
	    data.image[ID_PIAAM1].array.F[jj*size+ii] = 0.0;
	}
    }
  printf("\n\n");

  free(r0array);
  free(z0array);
  free(r1array);
  free(z1array);

  return 0;
}



int PIAACMCsimul_unwrapPhase()
{
  long IDa, IDp;
  long k;
  long ii, jj;
  long ID, IDm, IDd;
  long offset;
  double tmpf;
  char imname[200];
  char imnamem[200];
  char fname[200];
  char fnamem[200];
  long IDdm;
  long Cmsize, Fmsize;
  long size;
  
  size = piaacmc[0].size;
  Cmsize = piaacmc[0].Cmsize;
  Fmsize = piaacmc[0].Fmsize;

  // unwrap plane #2
  IDa = image_ID("WFamp_003"); 
  IDp = image_ID("WFpha_003"); 
  for(k=0; k<optsyst[0].nblambda; k++)
    {
      ID = create_2Dimage_ID("_tmppha", Cmsize, Cmsize);
      IDm = create_2Dimage_ID("_maskfit", Cmsize, Cmsize);
      offset = (size-Cmsize)/2;
      for(ii=0;ii<Cmsize;ii++)
	for(jj=0;jj<Cmsize;jj++)
	  {
	    data.image[ID].array.F[jj*Cmsize+ii] = data.image[IDp].array.F[k*size*size+(jj+offset)*size+(ii+offset)];
	    data.image[IDm].array.F[jj*Cmsize+ii] = data.image[IDa].array.F[k*size*size+(jj+offset)*size+(ii+offset)];
	  }      
      save_fits("_tmppha","!_tmppha.fits");
      save_fits("_maskfit","!_maskfit.fits");
      
      sprintf(imname, "_tmpphadxdy%02ld",k);
      sprintf(fname, "!_tmpphadxdy%02ld.fits",k);
 
      sprintf(imnamem, "_tmpmaskdxdy%02ld",k);
      sprintf(fnamem, "!_tmpmaskdxdy%02ld.fits",k);
      
      IDd = create_2Dimage_ID(imname, 2*Cmsize, Cmsize);
      IDdm = create_2Dimage_ID(imnamem, 2*Cmsize, Cmsize);
      for(ii=1;ii<Cmsize;ii++)
	for(jj=1;jj<Cmsize;jj++)
	  {
	    data.image[IDd].array.F[jj*Cmsize*2+ii] = data.image[ID].array.F[jj*Cmsize+ii] - data.image[ID].array.F[jj*Cmsize+ii-1];       
	    data.image[IDd].array.F[jj*Cmsize*2+(ii+Cmsize)] = data.image[ID].array.F[jj*Cmsize+ii] - data.image[ID].array.F[(jj-1)*Cmsize+ii];
	    data.image[IDdm].array.F[jj*Cmsize*2+ii] = sqrt(data.image[IDm].array.F[jj*Cmsize+ii] * data.image[IDm].array.F[jj*Cmsize+ii-1]);       
	    data.image[IDdm].array.F[jj*Cmsize*2+(ii+Cmsize)] = sqrt(data.image[IDm].array.F[jj*Cmsize+ii] * data.image[IDm].array.F[(jj-1)*Cmsize+ii]);
	  }
      for(ii=0;ii<Cmsize*Cmsize*2;ii++)
	data.image[IDd].array.F[ii] = M_PI*modf (data.image[IDd].array.F[ii]/M_PI, &tmpf);
     
     printf("Saving %s -> %s\n", imname, fname);
     save_fits(imname, fname);
     printf("Saving %s -> %s\n", imnamem, fnamem);
     save_fits(imnamem, fnamem);
     
     linopt_imtools_image_fitModes(imname, "Cmodesdxdy", imnamem, 0.01, "wfp002Cmodescoeff", 0);
     save_fits("wfp002Cmodescoeff", "!wfp002Cmodescoeff.fits");
     linopt_imtools_image_construct("Cmodes", "wfp002Cmodescoeff", "wfp002r");
     save_fits("wfp002r", "!wfp002r.fits");
     linopt_imtools_image_construct("Cmodesdxdy", "wfp002Cmodescoeff", "wfp002rdxdy");
     save_fits("wfp002rdxdy", "!wfp002rdxdy.fits");
   }

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
  double x, y;
  long IDa, IDp;
  long size;
  long nblambda;
  long size2;
  long ii, jj, k;
  long IDpiaa1z, IDpiaa2z;
  long elem;
  long kl;

  char fname_piaa1z[200];
  char fname_piaa2z[200];
  char fname_pupa0[200];
  char fname_pupp0[200];
  char fname[200];

  long ID;
  long index;

  double proplim = 1.0e-4;
  double total;

  long NBpiaacmcdesign = 1;

  float beamradpix;
  long ID0, ID1;
  long size0, size1;
  long Cmsize, Fmsize;


  // how to measure quality
  float focscale; // l/D per pix
  float scoringIWA = 2.0;
  float scoringOWA = 5.0;
  float scoringIWAx = 1.5;
  long IDsm;
  float r;

  // What to do
  int PIAACMCSIM_CREATESMODES = 0;




  printf("PIAACMC system simulation\n");

  size = 1024;


  // allocate design struct array
  piaacmc = (MIRRORPIAACMCDESIGN*) malloc(sizeof(MIRRORPIAACMCDESIGN)*NBpiaacmcdesign);

  piaacmc[0].beamrad = 0.01; // 20mm diam
  piaacmc[0].size = size;
  piaacmc[0].pixscale = 0.00005;
  piaacmc[0].piaasep = 0.5; // [m]
  piaacmc[0].fpzfactor = 4.0;
  piaacmc[0].Fratio = 30.0;

  piaacmc[0].centObs0 = 0.3; // input central obstruction
  piaacmc[0].centObs1 = 0.2; // output central obstruction
  piaacmc[0].NBradpts = 50000;
  piaacmc[0].r0lim = 1.1425; // outer radius after extrapolation, piaa mirror 0
  piaacmc[0].r1lim = 2.0; // outer radius after extrapolation, piaa mirror 1
 

  piaacmc[0].fpmRad = 50.0e-6; // outer radius [m]
  piaacmc[0].NBrings = 20; // number of rings
  piaacmc[0].fpmarraysize = 2048;
  piaacmc[0].fpmmaterial = 4; // PMMA

  focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;





  // create modes for aspheric optical surfaces description
  beamradpix = piaacmc[0].beamrad/piaacmc[0].pixscale;
 
  // x, y, r and PA coordinates in beam (for convenience & speed)
  IDx = create_2Dimage_ID("xcoord", size, size);
  IDy = create_2Dimage_ID("ycoord", size, size);
  IDr = create_2Dimage_ID("rcoord", size, size);
  IDPA = create_2Dimage_ID("PAcoord", size, size);
  for(ii=0; ii<size; ii++)
    for(jj=0; jj<size; jj++)
      {
	x = (1.0*ii-0.5*size)/beamradpix;
	y = (1.0*jj-0.5*size)/beamradpix;
	data.image[IDx].array.F[jj*size+ii] = x;
	data.image[IDy].array.F[jj*size+ii] = y;
	data.image[IDr].array.F[jj*size+ii] = sqrt(x*x+y*y);
	data.image[IDPA].array.F[jj*size+ii] = atan2(y,x);	
      }
  save_fl_fits("xcoord", "!xcoord.fits");
  save_fl_fits("ycoord", "!ycoord.fits");
  save_fl_fits("rcoord", "!rcoord.fits");
  save_fl_fits("PAcoord", "!PAcoord.fits");






  // ==================== CREATE MODES USED TO FIT AND DESCRIBE PIAA SHAPES ===============
 

  if(PIAACMCSIM_CREATESMODES == 1)
    {
      Cmsize = (long) (beamradpix*4);
      piaacmc[0].Cmsize = Cmsize;
      linopt_imtools_makeCosRadModes("Cmodes", Cmsize, 20, beamradpix, 2.0);
      piaacmc[0].CmodesID = image_ID("Cmodes");
      save_fits("Cmodes", "!Cmodes.fits");
      piaacmc[0].NBCmodes = data.image[ piaacmc[0].CmodesID].md[0].size[2];
      
      Fmsize = (long) (beamradpix*4);
      piaacmc[0].Fmsize = Fmsize;
      linopt_imtools_makeCPAmodes("Fmodes",  Fmsize, 5.0, 0.8, beamradpix, 2.0, 1);
      piaacmc[0].FmodesID = image_ID("Fmodes");
      save_fits("Fmodes", "!Fmodes.fits");
      piaacmc[0].NBFmodes = data.image[piaacmc[0].FmodesID].md[0].size[2];
      
      
      // mode derivatives
      ID = create_3Dimage_ID("Cmodesdxdy", 2*Cmsize, Cmsize, piaacmc[0].NBCmodes);
      for(k=0;k<piaacmc[0].NBCmodes;k++)
	for(ii=1;ii<Cmsize;ii++)
	  for(jj=1;jj<Cmsize;jj++)
	    {
	      data.image[ID].array.F[k*Cmsize*Cmsize*2+jj*Cmsize*2+ii] = data.image[piaacmc[0].CmodesID].array.F[k*Cmsize*Cmsize+jj*Cmsize+ii] - data.image[piaacmc[0].CmodesID].array.F[k*Cmsize*Cmsize+jj*Cmsize+ii-1];
	      data.image[ID].array.F[k*Cmsize*Cmsize*2+jj*Cmsize*2+(ii+Cmsize)] = data.image[piaacmc[0].CmodesID].array.F[k*Cmsize*Cmsize+jj*Cmsize+ii] - data.image[piaacmc[0].CmodesID].array.F[k*Cmsize*Cmsize+(jj-1)*Cmsize+ii];
	    }
      save_fits("Cmodesdxdy", "!Cmodesdxdy.fits");
      
      ID = create_3Dimage_ID("Fmodesdxdy", 2*Fmsize, Fmsize, piaacmc[0].NBFmodes);
      for(k=0;k<piaacmc[0].NBFmodes;k++)
	for(ii=1;ii<Fmsize;ii++)
	  for(jj=1;jj<Fmsize;jj++)
	    {
	      data.image[ID].array.F[k*Fmsize*Fmsize*2+jj*Fmsize*2+ii] = data.image[piaacmc[0].FmodesID].array.F[k*Fmsize*Fmsize+jj*Fmsize+ii] - data.image[piaacmc[0].FmodesID].array.F[k*Fmsize*Fmsize+jj*Fmsize+ii-1];
	      data.image[ID].array.F[k*Fmsize*Fmsize*2+jj*Fmsize*2+(ii+Fmsize)] = data.image[piaacmc[0].FmodesID].array.F[k*Fmsize*Fmsize+jj*Fmsize+ii] - data.image[piaacmc[0].FmodesID].array.F[k*Fmsize*Fmsize+(jj-1)*Fmsize+ii];
	    }
      save_fits("Fmodesdxdy", "!Fmodesdxdy.fits");
    }
  else
    {
      piaacmc[0].CmodesID = load_fits("Cmodes.fits", "Cmodes");
      piaacmc[0].NBCmodes = data.image[ piaacmc[0].CmodesID].md[0].size[2];
      piaacmc[0].Cmsize = data.image[piaacmc[0].CmodesID].md[0].size[0];

      piaacmc[0].FmodesID = load_fits("Fmodes.fits", "Fmodes");
      piaacmc[0].NBFmodes = data.image[ piaacmc[0].FmodesID].md[0].size[2];
      piaacmc[0].Fmsize = data.image[piaacmc[0].FmodesID].md[0].size[0];      
    }





  // =================== IMPORT / CREATE PIAA SHAPES =====================


  // CREATE FROM SCRATCH 
  if (0)
    {
      piaacmc[0].piaa0CmodesID = create_2Dimage_ID("piaa0Cmodescoeff",piaacmc[0].NBCmodes, 1);
      piaacmc[0].piaa0FmodesID = create_2Dimage_ID("piaa0Fmodescoeff",piaacmc[0].NBFmodes, 1);
      piaacmc[0].piaa1CmodesID = create_2Dimage_ID("piaa1Cmodescoeff",piaacmc[0].NBCmodes, 1);
      piaacmc[0].piaa1FmodesID = create_2Dimage_ID("piaa1Fmodescoeff",piaacmc[0].NBFmodes, 1);
    }
  else if (0) // read PIAAMshape from radial sag file
    {
      if(0) // initialize PIAACMC configuration with radial system
	{
	  load_fits("Apo2Drad.fits", "apo2D");
	  PIAACMCsimul_load2DRadialApodization("apo2D", 200.0, 0.3, "outApofit");
	  
	  PIAACMCsimul_init_geomPIAA_rad("outApofit");
	  
	  PIAACMCsimul_mkPIAAMshapes_from_RadSag("PIAA_Mshapes.txt", "piaam0z", "piaam1z");
	  save_fits("piaam0z", "!piaam0z.fits");
	  save_fits("piaam1z", "!piaam1z.fits");
	}
      else
	{
	  load_fits("piaam0z.fits", "piaam0z");
	  load_fits("piaam1z.fits", "piaam1z");
	}

      
 
      // ========== FIT SHAPES ========
      
      // crop piaam0z and piaam1z to Cmodes size
      ID0 = image_ID("Cmodes");
      size0 = data.image[ID0].md[0].size[0];
      ID1 = create_2Dimage_ID("piaa0zcrop", size0, size0);
      ID = image_ID("piaam0z");
      for(ii=0;ii<size0;ii++)
	for(jj=0;jj<size0;jj++)
	  data.image[ID1].array.F[jj*size0+ii] = data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)];
      ID1 = create_2Dimage_ID("piaa1zcrop", size0, size0);
      ID = image_ID("piaam1z");
      for(ii=0;ii<size0;ii++)
	for(jj=0;jj<size0;jj++)
	  data.image[ID1].array.F[jj*size0+ii] = data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)];

      make_disk("maskd", size0, size0, 0.5*size0, 0.5*size0, beamradpix);      
      make_2Dgridpix("gridpix", size0, size0, 1, 1, 0, 0);
      arith_image_mult("maskd", "gridpix", "maskfit");
      save_fits("maskfit", "!maskfit.fits");

      printf("--------- FITTING COSINE MODES ---------\n");
      fflush(stdout);

      linopt_imtools_image_fitModes("piaa0zcrop", "Cmodes", "maskfit", 0.01, "piaa0Cmodescoeff", 0);
      save_fits("piaa0Cmodescoeff", "!piaa0Cmodescoeff.fits");
      linopt_imtools_image_fitModes("piaa1zcrop", "Cmodes", "maskfit", 0.01, "piaa1Cmodescoeff", 0);
      save_fits("piaa1Cmodescoeff", "!piaa1Cmodescoeff.fits");

      linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");
      save_fits("piaa0Cz", "!piaa0Cz.fits");
      linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");
      save_fits("piaa1Cz", "!piaa1Cz.fits");

      ID0 = image_ID("piaa0Cz");
      size0 = data.image[ID0].md[0].size[0];
      ID1 = image_ID("piaam0z");
      ID = create_2Dimage_ID("piaa0Cres", size0, size0);
      for(ii=0;ii<size0;ii++)
	for(jj=0;jj<size0;jj++)
	  data.image[ID].array.F[jj*size0+ii] = data.image[ID1].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)]-data.image[ID0].array.F[jj*size0+ii];
      save_fits("piaa0Cres", "!piaa0Cres.fits");

      ID0 = image_ID("piaa1Cz");
      size0 = data.image[ID0].md[0].size[0];
      ID1 = image_ID("piaam1z");
      ID = create_2Dimage_ID("piaa1Cres", size0, size0);
      for(ii=0;ii<size0;ii++)
	for(jj=0;jj<size0;jj++)
	  data.image[ID].array.F[jj*size0+ii] = data.image[ID1].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)]-data.image[ID0].array.F[jj*size0+ii];
      save_fits("piaa1Cres", "!piaa1Cres.fits");
       
      printf("--------- FITTING FOURIER MODES ---------\n");
      fflush(stdout);

      linopt_imtools_image_fitModes("piaa0Cres", "Fmodes", "maskfit", 0.01, "piaa0Fmodescoeff", 0);
      save_fits("piaa0Fmodescoeff", "!piaa0Fmodescoeff.fits");
      linopt_imtools_image_fitModes("piaa1Cres", "Fmodes", "maskfit", 0.01, "piaa1Fmodescoeff", 0);
      save_fits("piaa1Fmodescoeff", "!piaa1Fmodescoeff.fits");
      
      linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
      save_fits("piaa0Fz", "!piaa0Fz.fits");
      arith_image_sub("piaa0Cres", "piaa0Fz", "piaa0CFres");
      save_fits("piaa0CFres", "!piaa0CFres.fits");
      delete_image_ID("piaa0zcrop");

      linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
      save_fits("piaa1Fz", "!piaa1Fz.fits");
      arith_image_sub("piaa1Cres", "piaa1Fz", "piaa1CFres");
      save_fits("piaa1CFres", "!piaa1CFres.fits");
      delete_image_ID("piaa1zcrop");

      delete_image_ID("maskfit");

      piaacmc[0].piaa0CmodesID = image_ID("piaa0Cmodescoeff");
      piaacmc[0].piaa0FmodesID = image_ID("piaa0Fmodescoeff");
      piaacmc[0].piaa1CmodesID = image_ID("piaa1Cmodescoeff");
      piaacmc[0].piaa1FmodesID = image_ID("piaa1Fmodescoeff");


    }
  else if (1) // read from file
    {
      piaacmc[0].piaa0CmodesID = load_fits("piaa0Cmodescoeff.fits", "piaa0Cmodescoeff");
      piaacmc[0].piaa0FmodesID = load_fits("piaa0Fmodescoeff.fits", "piaa0Fmodescoeff");
      piaacmc[0].piaa1CmodesID = load_fits("piaa1Cmodescoeff.fits", "piaa1Cmodescoeff");
      piaacmc[0].piaa1FmodesID = load_fits("piaa1Fmodescoeff.fits", "piaa1Fmodescoeff");

    }


 



  // ============ construct PIAA shapes from fitting coefficients ==================
  
  // assemble piaa0z and piaa1z images
  ID0 = linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");
  ID1 = linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
  ID = create_2Dimage_ID("piaa0z", size, size);
  size0 = data.image[ID0].md[0].size[0];
  size1 = data.image[ID1].md[0].size[0];
  for(ii=0;ii<size0;ii++)
    for(jj=0;jj<size0;jj++)
      data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)] += data.image[ID0].array.F[jj*size0+ii];
  for(ii=0;ii<size1;ii++)
    for(jj=0;jj<size1;jj++)
      data.image[ID].array.F[(jj+(size-size1)/2)*size+(ii+(size-size1)/2)] += data.image[ID1].array.F[jj*size1+ii];
  save_fits("piaa0Cz", "!piaa0Cz.fits");
  save_fits("piaa0Fz", "!piaa0Fz.fits");
  save_fits("piaa0z", "!piaa0z.fits");
  
  ID0 = linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");
  ID1 = linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
  ID = create_2Dimage_ID("piaa1z", size, size);
  size0 = data.image[ID0].md[0].size[0];
  size1 = data.image[ID1].md[0].size[0];
  for(ii=0;ii<size0;ii++)
    for(jj=0;jj<size0;jj++)
      data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)] += data.image[ID0].array.F[jj*size0+ii];
  for(ii=0;ii<size1;ii++)
    for(jj=0;jj<size1;jj++)
      data.image[ID].array.F[(jj+(size-size1)/2)*size+(ii+(size-size1)/2)] += data.image[ID1].array.F[jj*size1+ii];
  save_fits("piaa1Cz", "!piaa1Cz.fits");
  save_fits("piaa1Fz", "!piaa1Fz.fits");
  save_fits("piaa1z", "!piaa1z.fits");



  // ============ MAKE FOCAL PLANE MASK ===============

  PIAACMCsimul_mkFPM_zonemap("fpmzmap");
  save_fits("fpmzmap", "!fpmzmap.fits");
  piaacmc[0].zonezID = create_2Dimage_ID("fpmzt", piaacmc[0].focmNBzone, 1);
  for(ii=0;ii<piaacmc[0].focmNBzone;ii++)
    {
      if(ii>10)
	data.image[piaacmc[0].zonezID].array.F[ii] = 1.0e-6;
    }

  // ========== initializes optical system to piaacmc design ===========
  PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);
  
  // ============ perform propagations ================
  OptSystProp_run(optsyst, 0, 0, optsyst[0].NBelem);
  
  
  // PIAACMCsimul_unwrapPhase();
  printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);




  // CREATE SCORING MASK
  printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);
  IDsm = create_3Dimage_ID("scoringmask", size, size, optsyst[0].nblambda);
  for(k=0;k<optsyst[0].nblambda;k++)
    {
      for(ii=0;ii<size;ii++)
	for(jj=0;jj<size;jj++)
	  {
	    x = (1.0*ii-0.5*size)*focscale;
	    y = (1.0*jj-0.5*size)*focscale;
	    r = sqrt(x*x+y*y);
	    if((r>scoringIWA)&&(r<scoringOWA)&&(x>scoringIWAx))
	      data.image[IDsm].array.F[k*size2+jj*size+ii] = 1.0;
	  }
    }
  save_fits("scoringmask","!scoringmask.fits");
  
  

  list_image_ID();
  free(piaacmc);



  return(0);
}

