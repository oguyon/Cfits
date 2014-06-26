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
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "info/info.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate/WFpropagate.h"
#include "statistic/statistic.h"
#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "image_filter/image_filter.h"

#include "PIAACMCsimul/PIAACMCsimul.h"
#include "OptSystProp/OptSystProp.h"


# ifdef HAVE_LIBGOMP
#include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
#endif


/**
 * @file PIAACMCsimul.c
 * @author Olivier Guyon
 */




extern DATA data;

#define SBUFFERSIZE 2000

///  Current configuration directory
char piaacmcconfdir[200];

OPTSYST *optsyst;
int optsystinit = 0;
long IDx, IDy, IDr, IDPA;


double LAMBDASTART = 0.5e-6;
double LAMBDAEND = 0.6e-6;
#define NBLAMBDA 5

MIRRORPIAACMCDESIGN *piaacmc;
OPTSYST *optsyst;


int FORCE_CREATE_Cmodes = 0;
int CREATE_Cmodes = 0;
int FORCE_CREATE_Fmodes = 0;
int CREATE_Fmodes = 0;

int FORCE_CREATE_fpmzmap = 0;
int CREATE_fpmzmap = 0;
int FORCE_CREATE_fpmzt = 0;
int CREATE_fpmzt = 0;

int FORCE_CREATE_fpmza = 0;
int CREATE_fpmza;

int FORCE_MAKE_PIAA0shape = 0;
int MAKE_PIAA0shape = 0;
int FORCE_MAKE_PIAA1shape = 0;
int MAKE_PIAA1shape = 0;

int focmMode = -1; // if != -1, compute only impulse response to corresponding zone



// declared here for speed

double evalval;
long evali; 
long evalk, evalki, evalki1, evalmz, evalii, evalii1, evalii2, evalkv;
double evalcosp, evalsinp, evalre, evalim, evalre1, evalim1, evalpha;
double evalv1;



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
  
  if(CLI_checkarg(1,2)+CLI_checkarg(2,2)==0)
    {
      PIAACMCsimul_run(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl);
      return 0;
    }
  else
    return 1;
}


/**
 * Initializes module
 */
int init_PIAACMCsimul()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "PIAACMC system simulation");
  data.NBmodule++;
  


  strcpy(data.cmd[data.NBcmd].key,"piaacmcsimrun");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = PIAACMCsimul_run_cli;
  strcpy(data.cmd[data.NBcmd].info,"Simulate PIAACMC");
  strcpy(data.cmd[data.NBcmd].syntax,"<configuration index [long]>");
  strcpy(data.cmd[data.NBcmd].example,"piaacmcsimrun");
  strcpy(data.cmd[data.NBcmd].Ccall,"int PIAACMCsimul_run(int confindex)");
  data.NBcmd++;
 
   
  // add atexit functions here
  atexit(PIAACMCsimul_free);

  return 0;

}

/**
 * Frees memory for module
 */
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

/**
 * @param[out]  IDname  Name of output image
 */

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
  ID = create_image_ID(IDname, 2, sizearray, USHORT, 0, 0);
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
	data.image[ID].array.U[jj*piaacmc[0].fpmarraysize+ii] = (unsigned short int) zi;
      }
  piaacmc[0].focmNBzone = piaacmc[0].NBrings;
  
  return ID;
}




//
// makes 1-fpm CA
// if mode = -1, make whole 1-fpm
// if mode = zone, make only 1 zone with CA = (1.0, 0.0) 
//
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
  int zi;
  double t, a, amp;
  long size2;
  double pha, re, im;
   
  size = optsyst[0].size;
  size2 = size*size;
  nblambda = optsyst[0].nblambda;
 

  IDz = image_ID(IDzonemap_name);
  ID = create_3DCimage_ID(ID_name, size, size, nblambda);
 
  // CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/CORONAGRAPHS_ARRAYSIZE

  /*  printf("%ld %ld\n", piaacmc[0].zoneaID, piaacmc[0].zonezID);
  for(k=0;k<data.image[piaacmc[0].zonezID].md[0].size[0];k++)
    {
      t = data.image[piaacmc[0].zonezID].array.D[k]; // thickness
      a = data.image[piaacmc[0].zoneaID].array.D[k]; // amplitude transmission
      amp = a;
      pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[0]);
      
      re = amp*cos(pha);
      im = amp*sin(pha);

      printf("ZONE %2ld  %12g %12g   %12g %12g    %12g %12g\n", k, a, t, amp, pha, 1.0-re, -im);
    }
  */
   


  for(k=0;k<nblambda;k++)
    {
      fpscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size/piaacmc[0].fpzfactor*optsyst[0].lambdaarray[k]*piaacmc[0].Fratio;
      printf("LAMBDA = %10.5g m    SCALE = %10.5g m/pix   size=%4ld  rad=%g\n", optsyst[0].lambdaarray[k], fpscale, size, piaacmc[0].fpmRad);
      

      for(ii=0;ii<size;ii++)
	for(jj=0;jj<size;jj++)
	  {
	    x = (1.0*ii-size/2)*fpscale; // [m]
	    y = (1.0*jj-size/2)*fpscale; // [m]	    

	    ii1 = (long) ( (0.5 + 0.5*x/piaacmc[0].fpmRad)*piaacmc[0].fpmarraysize);
	    jj1 = (long) ( (0.5 + 0.5*y/piaacmc[0].fpmRad)*piaacmc[0].fpmarraysize);
	    if((ii1>-1)&&(ii1<piaacmc[0].fpmarraysize)&&(jj1>-1)&&(jj1<piaacmc[0].fpmarraysize))
	      {
		zi = data.image[IDz].array.U[jj1*piaacmc[0].fpmarraysize+ii1];
		t = data.image[piaacmc[0].zonezID].array.D[zi-1]; // thickness
		a = data.image[piaacmc[0].zoneaID].array.D[zi-1]; // amplitude transmission
	      }
	    else
	      {
		zi = 0;
		t = 0.0;
		a = 1.0;
	      }

	    if(mode == -1)
	      {
		if(zi>0.1)
		  {
		    amp = a;
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
	    else // impusle response from single zone
	      {
		if(mode == zi)  
		  {
		    amp = 1.0;
		    pha = 0.0; //OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[k]);
		    re = amp*cos(pha);
		    im = amp*sin(pha);
		    data.image[ID].array.CF[k*size2+jj*size+ii].re = 1.0;
		    data.image[ID].array.CF[k*size2+jj*size+ii].im = 0.0;
		  }
		else
		  {
		    data.image[ID].array.CF[k*size2+jj*size+ii].re = 0.0;
		    data.image[ID].array.CF[k*size2+jj*size+ii].im = 0.0;
		  }
	      }
	    
	  }
    }

  return(ID);
}
















/// initializes the optsyst structure to simulate reflective PIAACMC system

void PIAACMCsimul_init( MIRRORPIAACMCDESIGN *design, long index, double TTxld, double TTyld )
{
  FILE *fp;
  long k, i;
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
  char fname[200];
  long IDv;

  optsyst[0].nblambda = design[index].nblambda;
  nblambda = optsyst[0].nblambda;
  
  fp = fopen("lambdalist.txt", "w");
  for(k=0;k<optsyst[0].nblambda;k++)
    {
      optsyst[0].lambdaarray[k] = LAMBDASTART + (0.5+k)*(LAMBDAEND-LAMBDASTART)/optsyst[0].nblambda;
      fprintf(fp, "%02ld %20g\n", k, optsyst[0].lambdaarray[k]);
    }
  fclose(fp);

  optsyst[0].beamrad = design[index].beamrad; // 8mm
  optsyst[0].size = design[index].size;
  size = optsyst[0].size;
  size2 = size*size;
  optsyst[0].pixscale = design[index].pixscale;
  optsyst[0].DFTgridpad = 0; // 0 for full DFT sampling, >0 for faster execution

  beamradpix = optsyst[0].beamrad/optsyst[0].pixscale;
  printf("BEAM RADIUS = %f pix\n", beamradpix);

  if((IDv=variable_ID("PIAACMC_dftgrid"))!=-1)
    optsyst[0].DFTgridpad = (long) (data.variable[IDv].value+0.001);


  // define optical elements and locations
  optsyst[0].NB_DM = 0;
  optsyst[0].NB_asphsurfm = 2;
  optsyst[0].NB_asphsurfr = 0;
  
  optsyst[0].NBelem = 100; // to be updated later 


  fp = fopen("conjugations.txt","w");

  elem = 0;
  // ------------------- elem 0: input pupil -----------------------
 
  optsyst[0].elemtype[elem] = 1; // pupil mask
    // input pupil
  sprintf(fname_pupa0, "%s/pupa0_%ld.fits", piaacmcconfdir, size);
  if(file_exists(fname_pupa0)==1)
    load_fits(fname_pupa0, "pupa0");
  IDa = image_ID("pupa0");
  if(IDa==-1)
    {
      printf("CREATING INPUT PUPIL\n");
      if(IDa!=-1)
	delete_image_ID("pupa0");
      IDa = create_3Dimage_ID("pupa0", size, size, nblambda);
      
      ID = image_ID("telpup");
      if(ID==-1)
	if(file_exists("telpup.fits")==1)
	  ID = load_fits("telpup.fits", "telpup");
      

      if(ID==-1)
	{
	  for(k=0;k<nblambda;k++)
	    for(ii=0;ii<size2;ii++)
	      {
		if((data.image[IDr].array.F[ii]>piaacmc[0].centObs0)&&(data.image[IDr].array.F[ii]<1.0))
		  data.image[IDa].array.F[k*size2+ii] = 1.0;
		else
		  data.image[IDa].array.F[k*size2+ii] = 0.0;
	      }     
	}
      else
	 for(k=0;k<nblambda;k++)
	    for(ii=0;ii<size2;ii++)
	      {
		if(data.image[ID].array.F[ii]>0.5)
		  data.image[IDa].array.F[k*size2+ii] = 1.0;
		else
		  data.image[IDa].array.F[k*size2+ii] = 0.0;
	      }     

      sprintf(fname_pupa0, "!%s/pupa0_%ld.fits", piaacmcconfdir, size);
      save_fl_fits("pupa0", fname_pupa0);      
    }
  optsyst[0].elemarrayindex[elem] = IDa;  
  optsyst[0].elemZpos[elem] = 0.0;
  fprintf(fp,"%02ld  %f    Input pupil\n", elem, optsyst[0].elemZpos[elem]);
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
  fprintf(fp,"%02ld  %f    Fold mirror used to induce pointing offsets\n", elem, optsyst[0].elemZpos[elem]);
  elem++;

 
  IDpiaaz0 = image_ID("piaa0z");
  IDpiaaz1 = image_ID("piaa1z");


  // ------------------- elem 2: reflective PIAA M0  -----------------------  
  optsyst[0].elemtype[elem] = 3; // reflective PIAA M0   
  optsyst[0].elemarrayindex[elem] = 1; // index
  optsyst[0].ASPHSURFMarray[1].surfID = IDpiaaz0; 
  optsyst[0].elemZpos[elem] = 2.0;
  fprintf(fp,"%02ld  %f    PIAAM0\n", elem, optsyst[0].elemZpos[elem]);
  elem++;

  // ------------------- elem 3: reflective PIAA M1  -----------------------  
  optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
  optsyst[0].elemarrayindex[elem] = 2;
  optsyst[0].ASPHSURFMarray[2].surfID = IDpiaaz1; 
  optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1]+design[index].piaasep;
  fprintf(fp,"%02ld  %f    PIAAM1\n", elem, optsyst[0].elemZpos[elem]);
  elem++;


  // ------------------- elem 4 opaque mask at reflective PIAA M1  -----------------------  
  optsyst[0].elemtype[elem] = 1; // opaque mask
  ID = make_disk("piaam1mask", size, size, 0.5*size, 0.5*size, design[index].r1lim*beamradpix);
  optsyst[0].elemarrayindex[elem] = ID;
  optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1];
  save_fits("piaam1mask", "!piaam1mask.fits");
  fprintf(fp,"%02ld  %f    PIAAM1 edge opaque mask\n", elem, optsyst[0].elemZpos[elem]);
  elem++;
  

  // --------------------  elem 5: focal plane mask ------------------------
  optsyst[0].elemtype[elem] = 5; // focal plane mask 
  optsyst[0].elemarrayindex[elem] = 0;
  
  optsyst[0].FOCMASKarray[0].fpmID = PIAACMCsimul_mkFocalPlaneMask("fpmzmap", "piaacmcfpm", focmMode); // if -1, this is 1-fpm; otherwise, this is impulse response from single zone
  //  save_fits("fpmza", "!TESTfpmza.fits");
  //  save_fits("fpmzt", "!TESTfpmzt.fits");
  // save_fits("fpmzmap", "!TESTfpmzmap.fits");

  if(0)// testing
    {
      sprintf(fname, "!focma_%d.fits", focmMode);
      mk_amph_from_complex("piaacmcfpm", "fpma", "fpmp");
      save_fits("fpma", fname);
      save_fits("fpmp", "!fpmp.fits");
      delete_image_ID("fpma");
      delete_image_ID("fpmp");
    }
  


  optsyst[0]. FOCMASKarray[0].zfactor = design[index].fpzfactor;
  optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1]; // plane from which FT is done
  fprintf(fp,"%02ld  %f    post-focal plane mask pupil\n", elem, optsyst[0].elemZpos[elem]);
  elem++;

  
  // --------------------  elem 6, 7: Lyot masks  ------------------------
  for(i=0;i<design[index].NBLyotStop;i++)
    {
      optsyst[0].elemtype[elem] = 1; // Lyot mask 
      optsyst[0].elemarrayindex[elem] = design[index].IDLyotStop[i];  
      printf("elem %ld  Lyot mask %ld : %ld\n", elem, i, design[index].IDLyotStop[i]);
      optsyst[0].elemZpos[elem] =  design[index].LyotStop_zpos[i];
      fprintf(fp,"%02ld  %f  Lyot Stop %ld\n", elem, optsyst[0].elemZpos[elem], i);
      elem++;
    }



  // --------------------  elem 8: inv PIAA1 ------------------------
  optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
  optsyst[0].elemarrayindex[elem] = 2;
  optsyst[0].ASPHSURFMarray[2].surfID = IDpiaaz1; 
  optsyst[0].elemZpos[elem] = 0.0;
  fprintf(fp,"%02ld  %f    invPIAA1\n", elem, optsyst[0].elemZpos[elem]);
  elem++;

  // --------------------  elem 9: inv PIAA0 ------------------------
  optsyst[0].elemtype[elem] = 3; // reflective PIAA M0
  optsyst[0].elemarrayindex[elem] = 1;
  optsyst[0].ASPHSURFMarray[1].surfID = IDpiaaz0; 
  optsyst[0].elemZpos[elem] = design[index].piaasep;
  fprintf(fp,"%02ld  %f    invPIAA0\n", elem, optsyst[0].elemZpos[elem]);
  elem++;

  // --------------------  elem 9: back end mask  ------------------------

  optsyst[0].elemtype[elem] = 1; // Lyot mask 
  ID = make_disk("outmask", size, size, 0.5*size, 0.5*size, 0.92*design[index].beamrad/design[index].pixscale);
  optsyst[0].elemarrayindex[elem] = ID;  
  optsyst[0].elemZpos[elem] =  design[index].piaasep;
  fprintf(fp,"%02ld  %f   back end mask\n", elem, optsyst[0].elemZpos[elem]);
  elem++;
  
  fclose(fp);

  optsyst[0].NBelem = elem;

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


  if(0) // test fit quality
    {
      linopt_imtools_image_construct("APOmodesCos", IDapofit_name, "testapofitsol");
      save_fits("testapofitsol", "!testapofitsol.fits");
      arith_image_sub("_apoincrop", "testapofitsol", "apofitres");
      arith_image_mult("apofitres", "fitmaskapo", "apofitresm");
      save_fits("apofitres", "!apofitres.fits");
      save_fits("apofitresm", "!apofitresm.fits");
      // linopt_imtools_image_fitModes("apofitres", "APOmodesCos", "fitmaskapo", 1.0e-5, "test2c", 0);
      info_image_stats("apofitresm", "");
    }

  delete_image_ID("_apoincrop");
  delete_image_ID("fitmaskapo");


  return 0;
}







/**
 * computes radial PIAA optics sag
 *
 * this function only works for circular PIAA
 * uses radial PIAACMC design to initialize PIAA optics shapes and focal plane mask
 */

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

  char fname[200];

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
  //  fp = fopen("test1.txt", "w");
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
      //  fprintf(fp, "%lf %lf %lf\n", piaar11[i], piaar01[i], F0);
    }
  //  fclose(fp);
  




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

  sprintf(fname, "%s/PIAA_Mshapes.txt", piaacmcconfdir);
  fp = fopen(fname, "w");
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

# ifdef HAVE_LIBGOMP
#pragma omp parallel default(shared) private(ii, jj, x, y, r, k, r00, r01, alpha, val)
  {
# endif


# ifdef HAVE_LIBGOMP      
#pragma omp for
# endif
  for(ii=0;ii<size;ii++)
    {
      //      printf("\r %ld / %ld     ", ii, size);
      //fflush(stdout);


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
# ifdef HAVE_LIBGOMP
  }
# endif



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
// 
//
// 


// transmits between rin and rout
long PIAAsimul_mkSimpleLyotStop(char *ID_name, float rin, float rout)
{
  long size;
  long size2;
  long ii, k;
  long ID;
  float r;


  size = piaacmc[0].size;
  size2 = size*size;

  ID = create_3Dimage_ID(ID_name, size, size, piaacmc[0].nblambda);      
  for(k=0;k<piaacmc[0].nblambda;k++)
    for(ii=0;ii<size2;ii++)
      {
	if((data.image[IDr].array.F[ii]<rout)&&(data.image[IDr].array.F[ii]>rin))
	  data.image[ID].array.F[k*size2+ii] = 1.0;
	else
	  data.image[ID].array.F[k*size2+ii] = 0.0;
      }   

  return(ID);
}



/**
 * @brief Creates/initializes piaacmcconf structure and directory
 * 
 * @param[in] piaacmctype  Type of system
 * @param[in] fpmradld     Focal plane mask nominal radius   
 * @param[in] centobs0     Input central obstruction
 * @param[in] centobs1     Output central obstruction
 * @param[in] load         if 1, attempt to load configuration from file
 *
 * piaacmctype:
 * - 0: if configuration does not exist, create Monochromatic idealized PIAACMC, otherwise, read configuration
 * - 1: change array size to variable PIAACMC_size 
 *
 */

int PIAAsimul_initpiaacmcconf(long piaacmctype, double fpmradld, double centobs0, double centobs1, int load)
{
  float beamradpix;
  long NBpiaacmcdesign = 1;
  long ii, jj, k, i;
  double x, y;
  long size, size0;
  long Cmsize;
  long Fmsize;
  long ID, ID0, ID1;
  long size2;
  int r;
  char command[500];
  long IDv1, IDv2;
  char fname[200];
  char name[200];

  double pha0, t0, t;
  int loaded;

  int saveconf = 0; // if 1, save conf at end of this function
  long IDv;


  if(piaacmc == NULL)
    {
      piaacmc = (MIRRORPIAACMCDESIGN*) malloc(sizeof(MIRRORPIAACMCDESIGN)*NBpiaacmcdesign);

      
      // Default Values for PIAACMC (will adopt them unless configuration file exists)
      
      piaacmc[0].nblambda = 8;
      
      //piaacmc[0].nblambda = NBLAMBDA;
      
      
      // high resolution
      //piaacmc[0].size = 4096;
      //piaacmc[0].pixscale = 0.000055;
       
      // mid resolution
      piaacmc[0].size = 2048;
      piaacmc[0].pixscale = 0.000055;	  
 
      // low resolution
      //piaacmc[0].size = 1024;
      //piaacmc[0].pixscale = 0.00011;	  
       

      // very low resolution
      //    piaacmc[0].size = 512;
      // piaacmc[0].pixscale = 0.00022;	  
 

      piaacmc[0].beamrad = 0.022; // 44 mm diam
      piaacmc[0].piaasep = 2.302606; // [m]
      piaacmc[0].fpzfactor = 16.0;
      piaacmc[0].Fratio = 80.0;
      
      piaacmc[0].centObs0 = centobs0; // input central obstruction
      piaacmc[0].centObs1 = centobs1; // output central obstruction
      piaacmc[0].NBradpts = 50000;
      piaacmc[0].r0lim = 1.1425; // outer radius after extrapolation, piaa mirror 0
      piaacmc[0].r1lim = 1.5; // outer radius after extrapolation, piaa mirror 1
      
      piaacmc[0].NBLyotStop = 2;
      for(i=0;i<10;i++)
	{
	  piaacmc[0].LyotStop_zpos[i] = 0.0;
	  piaacmc[0].IDLyotStop[i] = -1;
	}
      
      piaacmc[0].fpmaskradld = fpmradld; // to compute prolate spheroidal function
      piaacmc[0].fpmarraysize = 2048;
      
      piaacmc[0].fpmRad = 100.0e-6; // focal plane radius [m]
      piaacmc[0].NBrings = 4; // number of rings in focal plane mask
      piaacmc[0].fpmmaterial = 4;  // PMMA
      piaacmc[0].fpmaskamptransm = 1.0;
      
      
      if(piaacmctype==0) // idealized focal plane mask
	{
	  piaacmc[0].NBrings = 1;
	  piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio;      
	}
      
      
      piaacmc[0].CmodesID = -1; // Cosine radial mode
      piaacmc[0].FmodesID = -1; // Fourier 2D modes
      piaacmc[0].piaa0CmodesID = -1;
      piaacmc[0].piaa0FmodesID = -1;
      piaacmc[0].piaa1CmodesID = -1;
      piaacmc[0].piaa1FmodesID = -1;
      piaacmc[0].zonezID = -1;  // focm zone material thickness, double precision image
      piaacmc[0].zoneaID = -1;  // focm zone amplitude transmission, double precision image
    }


  if(load==1)
    {
      printf("Loading PIAACMC configuration\n");
      fflush(stdout);
      sprintf(command, "mkdir -p %s", piaacmcconfdir);
      r = system(command);
      loaded = PIAAsimul_loadpiaacmcconf(piaacmcconfdir);
      if(loaded==0)
	{
	  printf("Saving default configuration\n");
	  fflush(stdout);
	  saveconf = 1;
	}
    }

  if((IDv=variable_ID("PIAACMC_nblambda"))!=-1)
    piaacmc[0].nblambda = data.variable[IDv].value;

  if((IDv=variable_ID("PIAACMC_size"))!=-1)
    piaacmc[0].size = (long) (data.variable[IDv].value+0.01);
  
  if((IDv=variable_ID("PIAACMC_pixscale"))!=-1)
    piaacmc[0].pixscale = data.variable[IDv].value;


  // create modes for aspheric optical surfaces description
  beamradpix = piaacmc[0].beamrad/piaacmc[0].pixscale;
  size = piaacmc[0].size;
  


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



  // ==================== CREATE MODES USED TO FIT AND DESCRIBE PIAA SHAPES ===============
  CREATE_Cmodes = 0;
  sprintf(fname, "%s/Cmodes.fits", piaacmcconfdir);
  if(FORCE_CREATE_Cmodes==0)
    {
      piaacmc[0].CmodesID = image_ID("Cmodes");
      if(piaacmc[0].CmodesID==-1)
	piaacmc[0].CmodesID = load_fits(fname, "Cmodes");
      if(piaacmc[0].CmodesID==-1)
	CREATE_Cmodes = 1;
    }
  else
    CREATE_Cmodes = 1;
  if(CREATE_Cmodes == 1)
    {
      if(piaacmc[0].CmodesID!=-1)
	delete_image_ID("Cmodes");      
      Cmsize = (long) (beamradpix*4);
      piaacmc[0].Cmsize = Cmsize;
      linopt_imtools_makeCosRadModes("Cmodes", Cmsize, 40, ApoFitCosFact*beamradpix, 2.0);
      piaacmc[0].CmodesID = image_ID("Cmodes");
      save_fits("Cmodes", fname);
    }
  piaacmc[0].NBCmodes = data.image[piaacmc[0].CmodesID].md[0].size[2];
  piaacmc[0].Cmsize = data.image[piaacmc[0].CmodesID].md[0].size[0];

  CREATE_Fmodes = 0;
  sprintf(fname, "%s/Fmodes.fits", piaacmcconfdir);
  if(FORCE_CREATE_Fmodes == 0)
    {
      piaacmc[0].FmodesID = image_ID("Fmodes");
      if(piaacmc[0].FmodesID==-1)
	piaacmc[0].FmodesID = load_fits(fname, "Fmodes");
      if(piaacmc[0].FmodesID==-1)
	CREATE_Fmodes = 1;
    }
  else
    CREATE_Fmodes = 1;
  if(CREATE_Fmodes == 1)
    {
      Fmsize = (long) (beamradpix*4);
      piaacmc[0].Fmsize = Fmsize;
      linopt_imtools_makeCPAmodes("Fmodes",  Fmsize, 5.0, 0.8, beamradpix, 2.0, 1);
      piaacmc[0].FmodesID = image_ID("Fmodes");
      save_fits("Fmodes", fname);
    }
  piaacmc[0].NBFmodes = data.image[piaacmc[0].FmodesID].md[0].size[2];      
  piaacmc[0].Fmsize = data.image[piaacmc[0].FmodesID].md[0].size[0];      
  


 
  // =================== IMPORT / CREATE PIAA SHAPES =====================

 
  piaacmc[0].piaa0CmodesID = image_ID("piaa0Cmodescoeff");
  piaacmc[0].piaa0FmodesID = image_ID("piaa0Fmodescoeff");
  piaacmc[0].piaa1CmodesID = image_ID("piaa1Cmodescoeff");
  piaacmc[0].piaa1FmodesID = image_ID("piaa1Fmodescoeff");
 
  if((piaacmc[0].piaa0CmodesID==-1)||( piaacmc[0].piaa0FmodesID==-1)||(piaacmc[0].piaa1CmodesID==-1)||( piaacmc[0].piaa1FmodesID==-1))
    {
      sprintf(fname, "%s/apo2Drad.fits", piaacmcconfdir);
      if(load_fits(fname, "apo2Drad")==-1)
	{
	  printf("Creating 2D apodization for idealized circular monochromatic PIAACMC\n");
	  fflush(stdout);
	  
	  // first iteration: half size image, 2x zoom
	  IDv1 = create_variable_ID("DFTZFACTOR", 2);
	  IDv2 = create_variable_ID("PNBITER", 15);	  
	  coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix*0.5, piaacmc[0].centObs1, "apotmp1", size/2);
	  
	  // expand solution to full size
	  basic_resizeim("apotmp1", "apostart", size, size);
	  delete_image_ID("apotmp1");

	  // full size, 4x zoom
	  IDv1 = create_variable_ID("DFTZFACTOR", 4);
	  IDv2 = create_variable_ID("PNBITER", 5);	  
	  coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size);

	  // full size, 8x zoom
	  chname_image_ID("apo", "apostart");
	  IDv1 = create_variable_ID("DFTZFACTOR", 8);
	  IDv2 = create_variable_ID("PNBITER", 5);	  
	  coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size);
	  
	  // full size, 20x zoom
	  chname_image_ID("apo", "apostart");
	  IDv1 = create_variable_ID("DFTZFACTOR", 16);
	  IDv2 = create_variable_ID("PNBITER", 10);	  
	  coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size);
	  

	  chname_image_ID("apo", "apo2Drad");
	  save_fits("apo2Drad", fname);	  	  

	  if((piaacmctype==0)&&(loaded==0)) // idealized focal plane mask
	    {
	      piaacmc[0].fpmaskamptransm =  -data.variable[variable_ID("APLCmaskCtransm")].value;
	      printf("FOCAL PLANE MASK TRANSM = %f\n", piaacmc[0].fpmaskamptransm);
	      printf("Saving default configuration\n");
	      fflush(stdout);
	      saveconf = 1;
	    }
	}

      // load apodization profile and fit it a series of cosines
      PIAACMCsimul_load2DRadialApodization("apo2Drad", beamradpix, 0.3, "outApofit");                  

      // compute radial PIAA sag 
      PIAACMCsimul_init_geomPIAA_rad("outApofit");
      
      // make 2D sag maps
      sprintf(fname, "%s/PIAA_Mshapes.txt", piaacmcconfdir);
      PIAACMCsimul_mkPIAAMshapes_from_RadSag(fname, "piaam0z", "piaam1z");

      sprintf(fname, "!%s/piaam0z.fits", piaacmcconfdir);
      save_fits("piaam0z", fname);

      sprintf(fname, "!%s/piaam1z.fits", piaacmcconfdir);
      save_fits("piaam1z", fname);



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

      linopt_imtools_image_fitModes("piaa0zcrop", "Cmodes", "maskfit", 1.0e-6, "piaa0Cmodescoeff", 0);
      save_fits("piaa0Cmodescoeff", "!piaa0Cmodescoeff.fits");
      linopt_imtools_image_fitModes("piaa1zcrop", "Cmodes", "maskfit", 1.0e-6, "piaa1Cmodescoeff", 0);
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
      //      save_fits("piaa0Fmodescoeff", "!piaa0Fmodescoeff.fits");
      linopt_imtools_image_fitModes("piaa1Cres", "Fmodes", "maskfit", 0.01, "piaa1Fmodescoeff", 0);
      // save_fits("piaa1Fmodescoeff", "!piaa1Fmodescoeff.fits");
      
      //linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
      //   save_fits("piaa0Fz", "!piaa0Fz.fits");
      //arith_image_sub("piaa0Cres", "piaa0Fz", "piaa0CFres");
      //save_fits("piaa0CFres", "!piaa0CFres.fits");
      delete_image_ID("piaa0zcrop");

      //linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
      //save_fits("piaa1Fz", "!piaa1Fz.fits");
      //arith_image_sub("piaa1Cres", "piaa1Fz", "piaa1CFres");
      //save_fits("piaa1CFres", "!piaa1CFres.fits");
      delete_image_ID("piaa1zcrop");

      delete_image_ID("maskfit");



      piaacmc[0].piaa0CmodesID = image_ID("piaa0Cmodescoeff");
      piaacmc[0].piaa0FmodesID = image_ID("piaa0Fmodescoeff");
      piaacmc[0].piaa1CmodesID = image_ID("piaa1Cmodescoeff");
      piaacmc[0].piaa1FmodesID = image_ID("piaa1Fmodescoeff");


      sprintf(fname, "!%s/piaa0Cmodes.fits", piaacmcconfdir);
      save_fits(data.image[piaacmc[0].piaa0CmodesID].md[0].name, fname);

      sprintf(fname, "!%s/piaa0Fmodes.fits", piaacmcconfdir);
      save_fits(data.image[piaacmc[0].piaa0FmodesID].md[0].name, fname);

      sprintf(fname, "!%s/piaa1Cmodes.fits", piaacmcconfdir);
      save_fits(data.image[piaacmc[0].piaa1CmodesID].md[0].name, fname);

      sprintf(fname, "!%s/piaa1Fmodes.fits", piaacmcconfdir);
      save_fits(data.image[piaacmc[0].piaa1FmodesID].md[0].name, fname);
    }




  // ============ MAKE FOCAL PLANE MASK ===============

  CREATE_fpmzmap = 0;
  if(FORCE_CREATE_fpmzmap == 0)
    { 
      if(image_ID("fpmzmap")==-1)
	CREATE_fpmzmap = 1;
    }
  else
    CREATE_fpmzmap = 1;
  if(CREATE_fpmzmap == 1)
    {
      if(image_ID("fpmzmap")!=-1)
	delete_image_ID("fpmzmap");
      PIAACMCsimul_mkFPM_zonemap("fpmzmap");
      sprintf(fname, "!%s/fpmzmap.fits", piaacmcconfdir);
      save_fits("fpmzmap", fname);
    }

  CREATE_fpmzt = 0;
  if(FORCE_CREATE_fpmzt == 0)
    {
      piaacmc[0].zonezID = image_ID("fpmzt");
      if(piaacmc[0].zonezID == -1)
	{
	  sprintf(fname, "%s/fpm_zonez.fits", piaacmcconfdir);
	  piaacmc[0].zonezID = load_fits(fname, "fpmzt");	  
	  if(piaacmc[0].zonezID == -1)
	    CREATE_fpmzt = 1;
	}
    }
  else
    CREATE_fpmzt = 1;
  if(CREATE_fpmzt == 1)
    {
      piaacmc[0].zonezID = image_ID("fpmzt");
      if(piaacmc[0].zonezID!=-1)
	delete_image_ID("fpmzt");

      piaacmc[0].zonezID = create_2Dimagedouble_ID("fpmzt", piaacmc[0].focmNBzone, 1);
      t = 1.0e-9;
      if((piaacmctype==0)&&(loaded==0)) // idealized focal plane mask
	{
	  t0 = 1.0e-8;
	  pha0 = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t0, 0.5*(LAMBDASTART+LAMBDAEND));
	  t = (M_PI/pha0)*t0;
	  printf("t = %g m (%lf %g)\n", t, pha0, t0);
	}
      for(ii=0;ii<piaacmc[0].focmNBzone;ii++)
	data.image[piaacmc[0].zonezID].array.D[ii] = t;
      sprintf(fname, "!%s/fpm_zonez.fits", piaacmcconfdir);
      save_fits("fpmzt", fname);
    }



  if(FORCE_CREATE_fpmza == 0)
    {
      piaacmc[0].zoneaID = image_ID("fpmza");
      if(piaacmc[0].zoneaID == -1)
	{
	  sprintf(fname, "%s/fpm_zonea.fits", piaacmcconfdir);
	  load_fits(fname, "fpmza");
	  
	  if(piaacmc[0].zoneaID == -1)
	    CREATE_fpmza = 1;
	}
    }
  else
    CREATE_fpmza = 1;
    
  if(CREATE_fpmza == 1)
    {
      if(piaacmc[0].zoneaID != -1)
	delete_image_ID("fpmza");
      piaacmc[0].zoneaID = create_2Dimagedouble_ID("fpmza", piaacmc[0].focmNBzone, 1);
      for(ii=0;ii<piaacmc[0].focmNBzone;ii++)
	data.image[piaacmc[0].zoneaID].array.D[ii] = piaacmc[0].fpmaskamptransm;          
      sprintf(fname, "!%s/fpm_zonea.fits", piaacmcconfdir);
      save_fits("fpmza", fname);
    }
    


  // ============= MAKE LYOT STOPS =======================
  printf("LOADING/CREATING LYOT MASK\n");
  size2 = size*size;

  for(i=0;i<piaacmc[0].NBLyotStop;i++)
    {      
      sprintf(fname, "%s/LyotStop%ld.fits", piaacmcconfdir, i);
      sprintf(name, "lyotstop%ld", i);
     
      piaacmc[0].IDLyotStop[i] = image_ID(name);	      
      if(piaacmc[0].IDLyotStop[i]==-1)
	{
	  sprintf(fname, "!%s/LyotStop%ld.fits", piaacmcconfdir, i);
	  switch (i) {
	  case 0 :
	    piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, -0.01, 0.98);
	    break;
	  case 1 :
	    piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, piaacmc[0].centObs1+0.02, 1.2);
	    break;
	  default :
	    piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, piaacmc[0].centObs1+0.02, 0.98);
	    break;
	  }
	  save_fl_fits(name, fname);      
	}
    }


  if(saveconf==1)
    PIAAsimul_savepiaacmcconf(piaacmcconfdir);

  return(0);
}





int PIAACMCsimul_makePIAAshapes()
{
  long ID, ID0, ID1;
  long size, size0, size1;
  long ii, jj;

  size = piaacmc[0].size;

  // ============ construct PIAA shapes from fitting coefficients ==================

  MAKE_PIAA0shape = 0;
  if(FORCE_MAKE_PIAA0shape == 0)
    {
      ID = image_ID("piaa0z");
      if(ID==-1)
	MAKE_PIAA0shape = 1;
    }
  else
    MAKE_PIAA0shape = 1;

  if(MAKE_PIAA0shape == 1)
    {
      // assemble piaa0z and piaa1z images
      ID0 = linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");
      ID1 = linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
      ID = image_ID("piaa0z");
      if(ID==-1)
	ID = create_2Dimage_ID("piaa0z", size, size);
      size0 = data.image[ID0].md[0].size[0];
      size1 = data.image[ID1].md[0].size[0];
      for(ii=0;ii<size*size;ii++)
	data.image[ID].array.F[ii] = 0.0;
      for(ii=0;ii<size0;ii++)
	for(jj=0;jj<size0;jj++)
	  data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)] += data.image[ID0].array.F[jj*size0+ii];
      for(ii=0;ii<size1;ii++)
	for(jj=0;jj<size1;jj++)
	  data.image[ID].array.F[(jj+(size-size1)/2)*size+(ii+(size-size1)/2)] += data.image[ID1].array.F[jj*size1+ii];
      save_fits("piaa0Cz", "!piaa0Cz.fits");
      save_fits("piaa0Fz", "!piaa0Fz.fits");
      save_fits("piaa0z", "!piaa0z.fits");
      delete_image_ID("piaa0Cz");
      delete_image_ID("piaa0Fz");
    }


  MAKE_PIAA1shape = 0;
  if(FORCE_MAKE_PIAA1shape == 0)
    {
      ID = image_ID("piaa1z");
      if(ID==-1)
	MAKE_PIAA1shape = 1;
    }
  else
    MAKE_PIAA1shape = 1;

  if(MAKE_PIAA1shape == 1)
    {
      ID0 = linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");
      ID1 = linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
      ID = image_ID("piaa1z");
      if(ID==-1)
	ID = create_2Dimage_ID("piaa1z", size, size);
      for(ii=0;ii<size*size;ii++)
	data.image[ID].array.F[ii] = 0.0;
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
      delete_image_ID("piaa1Cz");
      delete_image_ID("piaa1Fz");
    }

  return 0;
}




/// returns average contrast in evaluation zone

double PIAACMCsimul_computePSF(float xld, float yld, long startelem, long endelem)
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

  long size0, size1;
  long Cmsize, Fmsize;


  // how to measure quality
  float focscale; // l/D per pix
  float scoringIWA = 1.5;
  float scoringOWA = 8.0;
  float scoringIWAx = 0.5;
  long IDsm;
  float r;

  double value;
  double avContrast;
  double peakcontrast;
  double tmpv;

  printf("PIAACMC system simulation\n");

  size = piaacmc[0].size;
  size2 = size*size;


  focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;


  // ========== initializes optical system to piaacmc design ===========
  PIAACMCsimul_init(piaacmc, 0, xld, yld);
  PIAACMCsimul_makePIAAshapes();



  // ============ perform propagations ================
  OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem);
  
  
  printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);
  



  // CREATE SCORING MASK IF IT DOES NOT EXIST
  if((IDsm=image_ID("scoringmask"))==-1)
    {  
      printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);
      IDsm = create_2Dimage_ID("scoringmask", size, size);
      for(ii=0;ii<size;ii++)
	for(jj=0;jj<size;jj++)
	  {
	    x = (1.0*ii-0.5*size)*focscale;
	    y = (1.0*jj-0.5*size)*focscale;
	    r = sqrt(x*x+y*y);
	    if((r>scoringIWA)&&(r<scoringOWA)&&(x>scoringIWAx))
	      data.image[IDsm].array.F[jj*size+ii] = 1.0;
	    if((x>scoringIWA)&&(fabs(y)<scoringIWA*0.5)&&(r<50.0))
	      data.image[IDsm].array.F[jj*size+ii] = 1.0;
	  }
      save_fits("scoringmask","!scoringmask.fits");
   
      linopt_imtools_mask_to_pixtable("scoringmask", "pixindex", "pixmult");
    }

  linopt_imtools_Image_to_vec("psfc", "pixindex", "pixmult", "imvect");
  //save_fits("imvect", "!imvect.fits");
   
 
  value = 0.0;
  peakcontrast = 0.0;
  ID = image_ID("imvect");
  for(ii=0;ii<data.image[ID].md[0].nelement;ii++)
    {
      tmpv = data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
      value += tmpv;
      if(tmpv>peakcontrast)
	peakcontrast = tmpv;
    }

  for(elem=0;elem<optsyst[0].NBelem;elem++)
    printf("    FLUX %3ld   %12.4lf %8.6lf\n", elem, optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0]);
  value = value/size/size/optsyst[0].flux[0];

  avContrast = value/(arith_image_total("scoringmask")*focscale*focscale);
  printf("Peak constrast (rough estimate)= %g\n", peakcontrast/size/size/optsyst[0].flux[0]/focscale/focscale*3.0);
  printf("Total light in scoring field = %g  -> Average contrast = %g\n", value, value/(arith_image_total("scoringmask")*focscale*focscale));

  return(avContrast);
}




int PIAAsimul_savepiaacmcconf(char *dname)
{
  char command[200];
  int r;
  FILE *fp;
  char fname[200];
  long i;


  sprintf(command, "mkdir -p %s", dname);
  r = system(command);
  
  sprintf(fname,"%s/piaacmcparams.conf", dname);
  fp = fopen(fname, "w");

  
  fprintf(fp, "%10.6f   beamrad\n", piaacmc[0].beamrad);
  fprintf(fp, "%10ld    size\n", piaacmc[0].size);
  fprintf(fp, "%10.6g   pixscale\n", piaacmc[0].pixscale);
  fprintf(fp, "%10.6f   piaasep\n", piaacmc[0].piaasep);

  fprintf(fp, "%10.6f   centObs0\n", piaacmc[0].centObs0);
  fprintf(fp, "%10.6f   centObs1\n", piaacmc[0].centObs1);
  fprintf(fp, "%10ld   NBradpts\n", piaacmc[0].NBradpts);
  fprintf(fp, "%10.6f   r0lim\n", piaacmc[0].r0lim);
  fprintf(fp, "%10.6f   r1lim\n", piaacmc[0].r1lim);

  
  

  fprintf(fp, "%10ld    NBLyotStop\n", piaacmc[0].NBLyotStop);

  printf("%10ld    NBLyotStop\n", piaacmc[0].NBLyotStop);
  fflush(stdout);


  for(i=0;i<10;i++)
    {
      if(i<piaacmc[0].NBLyotStop)
	{
	  sprintf(fname, "!%s/LyotStop%ld.fits", dname, i);
	  if(piaacmc[0].IDLyotStop[i]!=-1)
	    save_fits(data.image[piaacmc[0].IDLyotStop[i]].md[0].name, fname);
	  fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
	}
      else
	fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
    }

  sprintf(fname, "!%s/piaa0Cmodes.fits", dname);
  if(piaacmc[0].piaa0CmodesID!=-1)
    save_fits(data.image[piaacmc[0].piaa0CmodesID].md[0].name, fname);

  sprintf(fname, "!%s/piaa0Fmodes.fits", dname);
  if(piaacmc[0].piaa0FmodesID!=-1)
    save_fits(data.image[piaacmc[0].piaa0FmodesID].md[0].name, fname);

  sprintf(fname, "!%s/piaa1Cmodes.fits", dname);
  if(piaacmc[0].piaa1CmodesID!=-1)
    save_fits(data.image[piaacmc[0].piaa1CmodesID].md[0].name, fname);

  sprintf(fname, "!%s/piaa1Fmodes.fits", dname);
  if(piaacmc[0].piaa1FmodesID!=-1)
    save_fits(data.image[piaacmc[0].piaa1FmodesID].md[0].name, fname);


  fprintf(fp, "%10.6f    fpmaskradld\n", piaacmc[0].fpmaskradld);
  fprintf(fp, "%10ld    focmNBzone\n", piaacmc[0].focmNBzone);
  fprintf(fp, "%10.6f   Fratio\n", piaacmc[0].Fratio);

  sprintf(fname, "!%s/fpm_zonez.fits", dname);
  if(piaacmc[0].zonezID!=-1)
    save_fits(data.image[piaacmc[0].zonezID].md[0].name, fname);

  fprintf(fp, "%10.6f    fpmaskamptransm\n", piaacmc[0].fpmaskamptransm);

  sprintf(fname, "!%s/fpm_zonea.fits", dname);
  if(piaacmc[0].zoneaID!=-1)
    save_fits(data.image[piaacmc[0].zoneaID].md[0].name, fname);



  fprintf(fp, "%10.6f   fpzfactor\n", piaacmc[0].fpzfactor);
  fprintf(fp, "%10.6g   fpmRad\n", piaacmc[0].fpmRad);
  fprintf(fp, "%10ld    NBrings\n", piaacmc[0].NBrings);
  fprintf(fp, "%10ld    fpmarraysize \n", piaacmc[0].fpmarraysize);
  fprintf(fp, "%10d    fpmmaterial\n", piaacmc[0].fpmmaterial);


  fclose(fp);
  


  
  return(0);
}




int PIAAsimul_loadpiaacmcconf(char *dname)
{
  char command[200];
  int r;
  FILE *fp;
  char fname[200];
  char imname[200];
  long i;
  
  int tmpi;
  long tmpl;
  float tmpf;
  double tmplf;

  

  sprintf(fname,"%s/piaacmcparams.conf", dname);

  fp = fopen(fname, "r");
  if(fp==NULL)
    {
      printf("Configuration file \"%s\" does not exist (yet), using previously set configuration\n", fname);
      fflush(stdout);
      r = 0;
    }
  else
    {
      r = fscanf(fp, "%f   beamrad\n", &tmpf);
      piaacmc[0].beamrad = tmpf;
      
      r = fscanf(fp, "%ld    size\n", &tmpl);
      piaacmc[0].size = tmpl;
      
      r = fscanf(fp, "%g   pixscale\n", &tmpf);
      piaacmc[0].pixscale = tmpf;
      
      r = fscanf(fp, "%f   piaasep\n", &tmpf);
      piaacmc[0].piaasep = tmpf;
      
      r = fscanf(fp, "%lf   centObs0\n", &tmplf);
      piaacmc[0].centObs0 = tmplf;
      
      r = fscanf(fp, "%lf   centObs1\n", &tmplf);
      piaacmc[0].centObs1 = tmplf;
      
      r = fscanf(fp, "%ld   NBradpts\n", &tmpl);
      piaacmc[0].NBradpts = tmpl;
      
      r = fscanf(fp, "%lf   r0lim\n", &tmplf);
      piaacmc[0].r0lim = tmplf;
      
      r = fscanf(fp, "%lf   r1lim\n", &tmplf);
      piaacmc[0].r1lim = tmplf;
      
      
      
      
      r = fscanf(fp, "%10ld    NBLyotStop\n", &piaacmc[0].NBLyotStop);
      for(i=0;i<10;i++)
	{
	  if(i<piaacmc[0].NBLyotStop)
	    {
	      sprintf(fname, "%s/LyotStop%ld.fits", dname, i);
	      sprintf(imname, "lyotstop%ld", i);
	      piaacmc[0].IDLyotStop[i] = load_fits(fname, imname);
	      r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
	      piaacmc[0].LyotStop_zpos[i] = tmplf;
	    }
	  else
	    {
	      r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
	      piaacmc[0].LyotStop_zpos[i] = tmplf;
	      printf("LYOT STOP %ld POS : %lf\n", i, tmplf);
	    }
	}
      
      
      sprintf(fname, "%s/piaa0Cmodes.fits", dname);
      piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff");
      
      sprintf(fname, "%s/piaa0Fmodes.fits", dname);
      piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff");
      
      sprintf(fname, "%s/piaa1Cmodes.fits", dname);
      piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff");
      
      sprintf(fname, "%s/piaa1Fmodes.fits", dname);
      piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff");
      



      r = fscanf(fp, "%f    fpmaskradld\n", &tmpf);
      piaacmc[0].fpmaskradld = tmpf;

      r = fscanf(fp, "%ld    focmNBzone\n",  &tmpl);
      piaacmc[0].focmNBzone = tmpl;
      
      r = fscanf(fp, "%f   Fratio\n",      &tmpf);
      piaacmc[0].Fratio = tmpf;
      
      sprintf(fname, "%s/fpm_zonez.fits", dname);
      delete_image_ID("fpmzt");
      piaacmc[0].zonezID = load_fits(fname, "fpmzt");

      r = fscanf(fp, "%f   fpmaskamptransm\n",    &tmpf);
      piaacmc[0].fpmaskamptransm = tmpf;

      sprintf(fname, "%s/fpm_zonea.fits", dname);
      delete_image_ID("fpmza");
      piaacmc[0].zoneaID = load_fits(fname, "fpmza");
      
      r = fscanf(fp, "%f   fpzfactor\n",   &tmpf);
      piaacmc[0].fpzfactor = tmpf;
      
      r = fscanf(fp, "%f   fpmRad\n",      &tmpf);
      piaacmc[0].fpmRad = tmpf;
      
      r = fscanf(fp, "%ld    NBrings\n",     &tmpl);
      piaacmc[0].NBrings = tmpl;
      
      r = fscanf(fp, "%ld    fpmarraysize \n", &tmpl);
      piaacmc[0].fpmarraysize = tmpl;
      
      r = fscanf(fp, "%d    fpmmaterial\n",  &tmpi);
      piaacmc[0].fpmmaterial = tmpi;
      
      r = 1;

      fclose(fp);
    }

  return(r);
}


/// Make Lyot stop geometry
/// param[in] IDincoh_name   Incoherent Lyot pupil intensity response to off-axis sources
/// parampin] IDmc_name      Intensity Lyot pupil image for on-axis source
//

long PIAACMCsimul_mkLyotMask(char *IDincoh_name, char *IDmc_name, char *IDzone_name, double throughput, char *IDout_name)
{
  long ID;
  long IDmc, IDincoh, IDzone;
  double val, val1;
  double rsl;
  long iter, NBiter;
  long ii;
  long xsize, ysize;
  long IDout;
  float sigma = 4.0;
  int filter_size = 10;

  NBiter = 10;

  
  printf("IDincoh_name : %s   %ld\n", IDincoh_name, image_ID(IDincoh_name));
  printf("IDmc_name    : %s   %ld\n", IDmc_name, image_ID(IDmc_name));
  printf("IDzone_name  : %s   %ld\n", IDzone_name, image_ID(IDzone_name));

  IDincoh = gauss_filter(IDincoh_name, "incohg", sigma, filter_size);
  IDmc = gauss_filter(IDmc_name, "mcg", sigma, filter_size);
  
  printf("STEP 0\n");
  fflush(stdout);


  //  IDincoh = image_ID(IDincoh_name);
  // IDmc = image_ID(IDmc_name);
  IDzone = image_ID(IDzone_name);
  xsize = data.image[IDmc].md[0].size[0];
  ysize = data.image[IDmc].md[0].size[1];

  IDout = create_2Dimage_ID(IDout_name, xsize, ysize);

  printf("STEP 1\n");
  fflush(stdout);

  // normalize both images to 1.0
  val = 0.0;
  for(ii=0;ii<xsize*ysize;ii++)
    val += data.image[IDmc].array.F[ii];
  for(ii=0;ii<xsize*ysize;ii++)
    data.image[IDmc].array.F[ii] /= val;

  val = 0.0;
  for(ii=0;ii<xsize*ysize;ii++)
    val += data.image[IDincoh].array.F[ii];
  for(ii=0;ii<xsize*ysize;ii++)
    data.image[IDincoh].array.F[ii] /= val;

  printf("STEP 1\n");
  fflush(stdout);


  rsl = 1.0;
  for(iter=0;iter<NBiter;iter++)
    {
      val = 0.0;
      val1 = 0.0;
      
      for(ii=0;ii<xsize*ysize;ii++)
	{
	  if((data.image[IDzone].array.F[ii]>-1)&&(data.image[IDincoh].array.F[ii]/data.image[IDmc].array.F[ii]>rsl))
	    {
	      val += data.image[IDincoh].array.F[ii];
	      val1 += data.image[IDmc].array.F[ii];
	      data.image[IDout].array.F[ii] = 1.0;
	    }
	  else
	    data.image[IDout].array.F[ii] = 0.0;
	}
      printf("rsl = %f  ->  %f %f\n", rsl, val, val1);
      if(val>throughput) // too much light came through
	rsl *= 1.1;
      else
	rsl *= 0.9;
    }

  delete_image_ID("incohg");
  delete_image_ID("mcg");

  return(IDout);
}



//
/// Lyot stops positions from zmin to zmax relative to current, working back (light goes from 0 to zmax)
/// @param[in] FluxTOT  total flux in current plane
/// @param[in] FluxLim   max flux allowed from star
// 
double PIAACMCsimul_optimizeLyotStop(char *IDamp_name, char *IDpha_name, char *IDincoh_name, float zmin, float zmax, double throughput, long NBz, long NBmasks)
{
  // initial guess places Lyot stops regularly from zmin to zmax
  // light propagates from zmin to zmax
  // we start with a single mask in zmax, and work back
  //

  double ratio;

  long ID, IDa, IDp;
  long nblambda; // number of wavelengths, read from input cube
  float *zarray;
  long l;
  double zprop;
  
  char nameamp[200];
  char namepha[200];
  char nameint[200];
  char fname[200];
  long xsize, ysize;
  long ii, jj, k, m;

  float *rinarray;
  float *routarray;
  float dr = 0.02;
  double tot;


  double *totarray;
  double *tot2array;
  long IDzone;
  double x, y, r;

  double zbest, valbest, val;
  long lbest;

  long IDincoh, IDint, IDmc;
  float rsl;
  long iter;
  long NBiter = 100;
  double val1;
  long IDm;
  char name[200];
  
  FILE *fp;
  double alpha = 1.05; // norm alpha used to identify best plane

  zarray = (float*) malloc(sizeof(float)*NBz);

  rinarray = (float*) malloc(sizeof(float)*NBmasks);
  routarray = (float*) malloc(sizeof(float)*NBmasks);

  totarray = (double*) malloc(sizeof(double)*NBmasks*NBz);
  tot2array = (double*) malloc(sizeof(double)*NBmasks*NBz);


  routarray[0] = 1.0;
  rinarray[0] = 1.0 - 1.0/NBmasks;
  for(m=1;m<NBmasks;m++)
    {
      routarray[m] = rinarray[m-1];
      rinarray[m] = routarray[m] - 1.0/NBmasks;
    }
  rinarray[NBmasks-1] = 0.0;

  for(m=0;m<NBmasks;m++)
    printf("annulus %ld : %f - %f\n", m, routarray[m], rinarray[m]);


  IDa = image_ID(IDamp_name);
  IDp = image_ID(IDpha_name);
  IDincoh = image_ID(IDincoh_name);

  xsize = data.image[IDa].md[0].size[0];
  ysize = data.image[IDa].md[0].size[1];

  if(data.image[IDa].md[0].naxis==3)
    nblambda = data.image[IDa].md[0].size[2];
  else
    nblambda = 1;

  IDzone = create_2Dimage_ID("LMzonemap", xsize, ysize);
  for(ii=0;ii<xsize;ii++)
    for(jj=0;jj<ysize;jj++)
      {
	data.image[IDzone].array.F[jj*xsize+ii] = -2;
	x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
	y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
	r = sqrt(x*x+y*y);
	for(m=0;m<NBmasks;m++)
	  if((r>rinarray[m])&&(r<routarray[m]))
	    data.image[IDzone].array.F[jj*xsize+ii] = m;
      }
  save_fits("LMzonemap", "!LMzonemap.fits");
  // initialize zarray
  for(l=0;l<NBz;l++)
    zarray[l] = zmin + (zmax-zmin)*l/(NBz-1);
  
  save_fits(nameint, fname);
  

  ID = create_3Dimage_ID("LMintC", xsize, ysize, NBz);
  for(l=0;l<NBz;l++)
    {
      sprintf(nameamp, "LMPamp%02ld", l);
      sprintf(namepha, "LMPpha%02ld", l);
      zprop = zarray[l];
      OptSystProp_propagateCube(optsyst, 0, IDamp_name, IDpha_name, nameamp, namepha, zprop);

      // collapse broadband intensities
      //  sprintf(nameint, "LMPint%02ld", l);
      IDa = image_ID(nameamp);
      //ID = create_2Dimage_ID(nameint, xsize, ysize);
 
      
      for(m=0;m<NBmasks;m++)
	{
	  tot2array[l*NBmasks+m] = 0.0;
	  totarray[l*NBmasks+m] = 0.0;
	}

      for(ii=0;ii<xsize*ysize;ii++)	
	{
	  m = (long) (data.image[IDzone].array.F[ii]+0.1);
	  for(k=0;k<nblambda;k++)
	    {
	      data.image[ID].array.F[l*xsize*ysize+ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii];      
	      if((m>-1)&&(m<NBmasks))
		{
		  totarray[l*NBmasks+m] += data.image[ID].array.F[l*xsize*ysize+ii];
		  tot2array[l*NBmasks+m] += pow(data.image[ID].array.F[l*xsize*ysize+ii], alpha);
		}
	    }
	}
      //sprintf(fname, "!LMPint%02ld.fits", l);
      save_fits(nameint, fname);

      delete_image_ID(nameamp);
      delete_image_ID(namepha);
    }

  save_fits("LMintC", "!LMintC.fits");
  IDmc = create_2Dimage_ID("Lcomb", xsize, ysize);

  fp = fopen("LyotMasks_zpos.txt", "w");
  IDint = image_ID("LMintC");
  for(m=0;m<NBmasks;m++)
    {
      valbest = 0.0;
      lbest = 0;
      zbest = 0.0;
      for(l=0;l<NBz;l++)
	{
	  val =  tot2array[l*NBmasks+m]/pow(totarray[l*NBmasks+m], alpha);
	  printf("MASK %ld   z= %f  ->  %g   ( %g %g) \n", m, zarray[l], val, tot2array[l*NBmasks+m], totarray[l*NBmasks+m]);
	  if(val>valbest)
	    {
	      valbest = val;
	      zbest = zarray[l];
	      lbest = l;
	    }
	}
      printf(" ==========  MASK %ld   BEST CONJUGATION : %ld %f (%g)\n", m, lbest, zbest, valbest);
      piaacmc[0].LyotStop_zpos[m] = zbest; // relative to starting plane
      fprintf(fp, "%02ld %f\n", lbest, zbest);
      //     sprintf(nameint, "LMPint%02ld", lbest);
      //      IDint = image_ID(nameint);

      for(ii=0;ii<xsize*ysize;ii++)
	if(m==data.image[IDzone].array.F[ii])
	  data.image[IDmc].array.F[ii] = data.image[IDint].array.F[lbest*xsize*ysize+ii];
    }
  fclose(fp);

  save_fits("Lcomb", "!Lcomb.fits");
  
  ID = PIAACMCsimul_mkLyotMask(IDincoh_name, "Lcomb", "LMzonemap", throughput, "LMask");
  delete_image_ID("Lcomb");

  for(m=0;m<NBmasks;m++)
    {
      sprintf(name, "optLM%02ld", m);
      IDm = create_2Dimage_ID(name, xsize, ysize);
      for(ii=0;ii<xsize;ii++)
	for(jj=0;jj<ysize;jj++)
	  {
	    x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
	    y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
	    r = sqrt(x*x+y*y);
	    
	    if((r>rinarray[m]-dr)&&(r<routarray[m]+dr))
	      data.image[IDm].array.F[jj*xsize+ii] = data.image[ID].array.F[jj*xsize+ii];
	    else
	      data.image[IDm].array.F[jj*xsize+ii] = 1.0;
	  }
      if(m==0)
	for(ii=0;ii<xsize;ii++)
	  for(jj=0;jj<ysize;jj++)
	    {
	      x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
	      y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
	      r = sqrt(x*x+y*y);
	      if(r>1.0)
		data.image[IDm].array.F[jj*xsize+ii] = 0.0;
	    }

      sprintf(fname, "!optLM%02ld.fits", m);
      save_fits(name, fname);
    }

 
  free(totarray);
  free(tot2array);
  free(rinarray);
  free(routarray);
  free(zarray);
  
  delete_image_ID("LMzonemap");

  return(ratio);
}


///
/// solves for focal plane mask solution using pre-computed zone responses
///
/// @param[in] fpmresp_array   Mask zones responses, float array
/// @param[in] zonez_array     zone thicknesses, double array
/// @param[in] dphadz_array    for each lambda, pha = thickness x dphadt_array[lambdaindex]      
/// @param[in] outtmp_array    output temp array
///
/// written to be fast, no checking of array sizes
/// all arrays pre-allocated outside this function
///
double PIAACMCsimul_achromFPMsol_eval(float *fpmresp_array, double *zonez_array, double *dphadz_array, float *outtmp_array, long vsize, long nbz, long nbl)
{
   
  // axis 0: eval pts (ii)   size = data.image[IDfpmresp].md[0].size[0] -> vsize
  // axis 1: zones (mz)      size = data.image[piaacmc[0].zonezID].md[0].size[0]+1 = nbz+1
  // axis 3: lambda (k)      size = piaacmc[0].nblambda -> nbl
  //
  // indexing :  k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*vsize + mz*vsize + ii

	
  for(evalk=0;evalk<nbl;evalk++)
    {
      evalki = evalk*(nbz+1)*vsize;
      
      // outer zone
      for(evalii=0;evalii<vsize;evalii++)
	outtmp_array[evalk*vsize+evalii] = fpmresp_array[evalk*(nbz+1)*vsize+evalii];
      
      // inner zones
      for(evalmz=0;evalmz<nbz;evalmz++) // zone
	{
	  // phase offset for this zone
	  evalpha = zonez_array[evalmz]*dphadz_array[evalk];
	  evalcosp = cos(evalpha);
	  evalsinp = sin(evalpha);
	  evalki1 = evalki + (evalmz+1)*vsize;
	  evalkv = evalk*vsize;
	  for(evalii=0;evalii<vsize/2;evalii++)
	    {
	      evalii1 = 2*evalii;
	      evalii2 = 2*evalii+1;
	      evalre = fpmresp_array[evalki1 + evalii1];
	      evalim = fpmresp_array[evalki1 + evalii2];
	      evalre1 = evalre*evalcosp - evalim*evalsinp;
	      evalim1 = evalre*evalsinp + evalim*evalcosp;
	      outtmp_array[evalkv + evalii1] += evalre1;
	      outtmp_array[evalkv + evalii2] += evalim1;
	    }
	}
    }
  evalval = 0.0;
  for(evalii=0;evalii<vsize*(nbz+1)*nbl;evalii++)
    {
      evalv1 = outtmp_array[evalii];
      evalval += evalv1*evalv1;
    }
  evalval /= vsize*(nbz+1)*nbl;
  
  return evalval;
}

/// derivative against parameter n
//double PIAACMCsimul_achromFPMsol_eval_d(float *fpmresp_array, double *zonez_array, double *dphadz_array, float *outtmp_array, long vsize, long nbz, long nbl, long n)





/**
 *
 * @brief Main simulation routine
 *
 * @param[in] confindex  PIAACMC configuration index pointing to the input/output directory number
 * @param[in] mode       Type of operation to be performed
 * 
 * mode values:
 * - 0 : compute on-axis propagation for specified configuration. If configuration does not exist, create idealized monochromatic PIAACMC (intended to create a new index) and compute on-axis propagation
 * - 1 : optimize Lyot stop(s) locations
 * - 2 : optimize focal plane mask transmission for idealized monochromatic PIAACMC
 * - 3 : TO BE WRITTEN
 * - 4 : linear optimization around current design, free parameters = PIAA optics shapes (runs for infinite number of iterations, track progress by looking at val.opt file)
 */

int PIAACMCsimul_run(long confindex, long mode)
{
  long NBparam;
  FILE *fp;
  char command[500];

  int paramtype[1000]; // FLOAT or DOUBLE
  double *paramval[1000]; // array of pointers, double
  float *paramvalf[1000]; // array of pointers, float
  double paramrefval[1000]; 

  double paramdelta[1000]; 
  double parammaxstep[1000]; // maximum single iteration step
  double parammin[1000]; // minimum value
  double parammax[1000]; // maximum value

  double paramdeltaval[1000];

  double valref, valbest, val0;
  double parambest[1000]; // for scanning
  double paramref[1000];
  long i, ii, jj;
  long IDv, ID, IDref, IDa;
  long IDmodes;
  long xsize, ysize;
  long k;

  long iter;
  long NBiter = 1000;

  long IDfpmresp, IDref1;
  double t, a, dpha, amp;
  int zi;

  char fname[200];

  long IDm, ID1D, ID1Dref;
  long size1Dvec;

  
  // OPTIMIZATION PARAMETERS
  int REGPIAASHAPES = 0;
  float piaa0C_regcoeff = 0.0e-7; // regularization coeff
  float piaa1C_regcoeff = 0.0e-7; // regularization coeff
  
  float piaa0C_regcoeff_alpha = 1.0; // regularization coeff power
  float piaa1C_regcoeff_alpha = 1.0; // regularization coeff power
  int r;

  double val, v1, pha, cosp, sinp, re, im, re1, im1;


  double fpmradld = 0.9;
  double centobs0 = 0.3;
  double centobs1 = 0.2;

  double range, stepsize;
  int loopOK;
  long NBls, ls, ls1;
  double lstransm;
  long mz;

  int LINOPT = 0; // 1 if perform linear optimization
  long vsize;
  long ii1, ii2, ki, kv, ki1;

  double scangain;
  double scanstepgain = 0.01;
  int linscanOK;
  double valold;
  double bestgain;
  double linoptgainarray[100];
  double linoptvalarray[100];
  int linoptlimflagarray[100];
  long NBlinoptgain;
  long kmax;


  float *fpmresp_array;
  double *zonez_array;
  double *dphadz_array;
  float *outtmp_array;
  double valtest;

  piaacmc = NULL;

  if(optsyst==NULL)
    optsyst = (OPTSYST*) malloc(sizeof(OPTSYST));


  sprintf(piaacmcconfdir, "piaacmcconf%03ld", confindex);

  
  switch (mode) {

  case 0 :  // Run existing config for on-axis point source. If new, create centrally obscured idealized PIAACMC
    if((IDv=variable_ID("PIAACMC_centobs0"))!=-1)
      centobs0 = data.variable[IDv].value;
    if((IDv=variable_ID("PIAACMC_centobs1"))!=-1)
      centobs1 = data.variable[IDv].value;
    if((IDv=variable_ID("PIAACMC_fpmradld"))!=-1)
      fpmradld = data.variable[IDv].value;
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    PIAACMCsimul_makePIAAshapes();      
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 
    valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
    printf("valref = %g\n", valref);
    break;
      

  case 1 : // optimize Lyot stop positions
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    PIAACMCsimul_makePIAAshapes();      
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 

    // initialization
    if((IDv=variable_ID("PIAACMC_lsoptrange"))!=-1)
      range = data.variable[IDv].value;
    else
      range = 2.0;
    stepsize = range/3.0;
    for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
      paramref[ls] = piaacmc[0].LyotStop_zpos[ls];
    NBiter = 4;
    
    

    fp = fopen("result.log", "w");
    fclose(fp);


    if(0) // systematic 3D grid search, slow
      {
	for(iter=0;iter<NBiter;iter++)
	  {
	    for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
	      {
		piaacmc[0].LyotStop_zpos[ls] = paramref[ls]-range;
		parambest[ls] = piaacmc[0].LyotStop_zpos[ls];
	      }
	    loopOK = 1;	
	    valbest = 1.0;
	    
	    while(loopOK==1)
	      {
		val = PIAACMCsimul_computePSF(0.0, 0.0, 6, optsyst[0].NBelem);
		
		if(val<valbest)
		  {
		    for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
		      parambest[ls] = piaacmc[0].LyotStop_zpos[ls];
		    valbest = val;
		  }
		
		fp = fopen("result.log", "a");
		for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
		  fprintf(fp," %lf", piaacmc[0].LyotStop_zpos[ls]);
		fprintf(fp, " %g\n", val);
		fclose(fp);
		
		ls = 0;
		piaacmc[0].LyotStop_zpos[ls] += stepsize;
		while((piaacmc[0].LyotStop_zpos[ls]>paramref[ls]+range+0.001*stepsize)&&(ls<piaacmc[0].NBLyotStop))
		  {
		    piaacmc[0].LyotStop_zpos[ls] = paramref[ls]-range;
		    ls++;
		    if(ls<piaacmc[0].NBLyotStop)
		      piaacmc[0].LyotStop_zpos[ls] += stepsize;
		  }
		if(ls==piaacmc[0].NBLyotStop)
		  loopOK = 0;
	      }
	    
	    printf("BEST SOLUTION :  ");
	    for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
	      {
		paramref[ls] = parambest[ls];
		printf(" %lf", parambest[ls]);
	      }
	    printf(" %g\n", valbest);
	    
	    
	    fp = fopen("result.log", "a");
	    fprintf(fp, "\n");
	    fclose(fp);
	    
	    range *= 0.3;
	    stepsize = range/3.0;
	  }
      }
    else // 1 axis at a time, much faster
      {
	stepsize = range/5.0;
	for(iter=0;iter<NBiter;iter++)
	  {	    
	    for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
	      {
		piaacmc[0].LyotStop_zpos[ls] = paramref[ls]-range;
		parambest[ls] = piaacmc[0].LyotStop_zpos[ls];
		
		loopOK = 1;	
		valbest = 1.0;
		
		while(piaacmc[0].LyotStop_zpos[ls]<paramref[ls]+range)
		  {
		    val = PIAACMCsimul_computePSF(0.0, 0.0, 6, optsyst[0].NBelem);
		    
		    if(val<valbest)
		      {
			parambest[ls] = piaacmc[0].LyotStop_zpos[ls];
			valbest = val;
		      }
		
		    fp = fopen("result.log", "a");
		    for(ls1=0;ls1<piaacmc[0].NBLyotStop;ls1++)
		      fprintf(fp," %lf", piaacmc[0].LyotStop_zpos[ls1]);
		    fprintf(fp, " %g\n", val);
		    fclose(fp);
		    
		    piaacmc[0].LyotStop_zpos[ls] += stepsize;
		  }
		printf("BEST SOLUTION :  ");
		paramref[ls] = parambest[ls];
		piaacmc[0].LyotStop_zpos[ls] = paramref[ls];
		printf(" %lf", parambest[ls]);
		printf(" %g\n", valbest);
	      }
	 
	    fp = fopen("result.log", "a");
	    fprintf(fp, "\n");
	    fclose(fp);
	    
	    range *= 0.3;
	    stepsize = range/3.0;	    
	  }
      }






    
    for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
      piaacmc[0].LyotStop_zpos[ls] = parambest[ls];    
    PIAAsimul_savepiaacmcconf(piaacmcconfdir);
    break;



  case 2 : // optimize focal plane mask transmission for monochromatic idealized PIAACMC
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    PIAACMCsimul_makePIAAshapes();      
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 

    // initialization
    range = 0.1;
    stepsize = range/3.0;
    paramref[0] = piaacmc[0].fpmaskamptransm;
    NBiter = 4;
    
    fp = fopen("result.log", "w");
    fclose(fp);

    for(iter=0;iter<NBiter;iter++)
      {
	piaacmc[0].fpmaskamptransm = paramref[0]-range;
	parambest[0] = piaacmc[0].fpmaskamptransm;

	loopOK = 1;	
	valbest = 1.0;
	
	while(loopOK==1)
	  {
	    FORCE_CREATE_fpmza = 1;
	    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
 	    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
	    
	    if(val<valbest)
	      {
		parambest[0] = piaacmc[0].fpmaskamptransm;
		valbest = val;
	      }
	    
	    fp = fopen("result.log", "a");
	    fprintf(fp," %lf", piaacmc[0].fpmaskamptransm);
	    fprintf(fp, " %g  %ld %g %g\n", val, iter, range, stepsize);
	    fclose(fp);
	    
	    ls = 0;
	    piaacmc[0].fpmaskamptransm += stepsize;
	    if(piaacmc[0].fpmaskamptransm>paramref[0]+range+0.001*stepsize)
	      loopOK = 0;
	  }
	
	printf("BEST SOLUTION :  ");

	paramref[0] = parambest[0];
	printf(" %lf", parambest[0]);
	
	printf(" %g\n", valbest);

	
	fp = fopen("result.log", "a");
	fprintf(fp, "\n");
	fclose(fp);

	range *= 0.3;
	stepsize = range/3.0;
      }
    
    piaacmc[0].fpmaskamptransm = parambest[0];    
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
    PIAAsimul_savepiaacmcconf(piaacmcconfdir);
    FORCE_CREATE_fpmza = 0;
    break;


  case 3 : // calibrate
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    PIAACMCsimul_makePIAAshapes();      
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 

    paramref[0] = piaacmc[0].fpmaskamptransm;

    piaacmc[0].fpmaskamptransm = -1.0;
    FORCE_CREATE_fpmza = 1;
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
    
    // restore original configuration
    piaacmc[0].fpmaskamptransm = paramref[0];    
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
    PIAAsimul_savepiaacmcconf(piaacmcconfdir);
    FORCE_CREATE_fpmza = 0;

    break;

  case 4 : // optimize PIAA optics shapes
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    LINOPT = 1; // perform linear optimization
    if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
      NBiter = (long) data.variable[IDv].value+0.01;
    else
      NBiter = 1000; 

    kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
    if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
      kmax = (long) data.variable[IDv].value+0.01;
 
    if(kmax>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
      kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];

    

    NBparam = 0;
    for(k=0;k<kmax;k++)
      {
	paramtype[NBparam] = FLOAT;
	paramvalf[NBparam] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
	paramdelta[NBparam] = 1.0e-9;
	parammaxstep[NBparam] = 1.0e-8;
	parammin[NBparam] = -1.0e-7;
	parammax[NBparam] = 1.0e-7;
	NBparam++;
      }
    
    for(k=0;k<kmax;k++)
      {
	paramtype[NBparam] = FLOAT;
	paramvalf[NBparam] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
	paramdelta[NBparam] = 1.0e-9;
	parammaxstep[NBparam] = 1.0e-8;
	parammin[NBparam] = -1.0e-7;
	parammax[NBparam] = 1.0e-7;
	NBparam++;
      }
    FORCE_MAKE_PIAA0shape = 1;
    FORCE_MAKE_PIAA1shape = 1;
    break;


    
  case 5 : // optimize Lyot stops shapes and positions
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    PIAACMCsimul_makePIAAshapes();      
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 
    
    NBls = 3;
    if((IDv=variable_ID("PIAACMC_nblstop"))!=-1)
      NBls = (long) data.variable[IDv].value+0.01;

    lstransm = 0.85;
    if((IDv=variable_ID("PIAACMC_lstransm"))!=-1)
      lstransm = (long) data.variable[IDv].value;
  

    PIAACMCsimul_computePSF(3.0, 0.0, 0, optsyst[0].NBelem);
    IDa = image_ID("WFamp_005");
    xsize = data.image[IDa].md[0].size[0];
    ysize = data.image[IDa].md[0].size[1];
    
    ID = create_2Dimage_ID("OAincoh", xsize, ysize);
    for(ii=0;ii<xsize*ysize;ii++)
      for(k=0;k<optsyst[0].nblambda;k++)
	data.image[ID].array.F[ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii]/4;
    
    PIAACMCsimul_computePSF(-3.0, 0.0, 0, optsyst[0].NBelem);
    IDa = image_ID("WFamp_005");
    for(ii=0;ii<xsize*ysize;ii++)
      for(k=0;k<optsyst[0].nblambda;k++)
	data.image[ID].array.F[ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii]/4;
    
    PIAACMCsimul_computePSF(0.0, 3.0, 0, optsyst[0].NBelem);
    IDa = image_ID("WFamp_005");
    for(ii=0;ii<xsize*ysize;ii++)
      for(k=0;k<optsyst[0].nblambda;k++)
	data.image[ID].array.F[ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii]/4;
    
    PIAACMCsimul_computePSF(0.0, -3.0, 0, optsyst[0].NBelem);
    IDa = image_ID("WFamp_005");
    for(ii=0;ii<xsize*ysize;ii++)
      for(k=0;k<optsyst[0].nblambda;k++)
	data.image[ID].array.F[ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii]/4;
    
    save_fits("OAincoh", "!OAincoh.fits");
    
    PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
    PIAACMCsimul_optimizeLyotStop("WFamp_005", "WFpha_005", "OAincoh", -20.0, 5.0, lstransm, 200, NBls);
    delete_image_ID("OAincoh");

    piaacmc[0].NBLyotStop = NBls;

    for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
      piaacmc[0].LyotStop_zpos[ls] += optsyst[0].elemZpos[5];

    PIAAsimul_savepiaacmcconf(piaacmcconfdir);
    for(ls=0;ls<piaacmc[0].NBLyotStop;ls++)
      {	
	sprintf(command, "cp optLM%02ld.fits ./%s/LyotStop%ld.fits", ls, piaacmcconfdir, ls);
	r = system(command);
      }
    break;


  case 6: // test off-axis performance
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    PIAACMCsimul_makePIAAshapes();      
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 
    PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
    break;
 

  case 10 : // setup polychromatic mask optimization
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 
    piaacmc[0].fpmaskamptransm = 1.0;
    piaacmc[0].NBrings = 16;
    piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio * 2.0; // 2 l/D radius at central lambda      
    FORCE_CREATE_fpmzmap = 1;
    FORCE_CREATE_fpmza = 1;
    FORCE_CREATE_fpmzt = 1;
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
    PIAAsimul_savepiaacmcconf(piaacmcconfdir);
    break;

  case 11 : // Compute polychromatic response to zones, store result in FPMresp
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    PIAACMCsimul_makePIAAshapes();      
        
    focmMode = data.image[piaacmc[0].zonezID].md[0].size[0]+10;
    optsyst[0].FOCMASKarray[0].mode = 1; // 1-fpm
    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
    printf("val = %g\n", val);
    ID = image_ID("imvect");


    // axis 0: eval pts (ii) - size = data.image[ID].md[0].size[0]
    // axis 1: zones (mz) - size = data.image[piaacmc[0].zonezID].md[0].size[0]+1
    // axis 3: lambda (k) - size = piaacmc[0].nblambda
    //
    // indexing :  k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + mz*data.image[ID].md[0].size[0] + ii
    IDfpmresp = image_ID("FPMresp");
    if(IDfpmresp==-1)
      IDfpmresp = create_3Dimage_ID("FPMresp", data.image[ID].md[0].size[0], data.image[piaacmc[0].zonezID].md[0].size[0]+1, piaacmc[0].nblambda);
    

    //    for(k=0;k<piaacmc[0].nblambda;k++)
    //    for(k=0;k<1;k++)
    //{


    for(k=0;k<piaacmc[0].nblambda;k++)
      for(ii=0;ii<data.image[ID].md[0].size[0];ii++)
	data.image[IDfpmresp].array.F[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + ii] = data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];	      
    save_fits("FPMresp","!FPMresp.fits");
    
    for(mz=1;mz<data.image[piaacmc[0].zonezID].md[0].size[0]+1;mz++)
      {	
	focmMode = mz;
	optsyst[0].FOCMASKarray[0].mode = 0; 
	val = PIAACMCsimul_computePSF(0.0, 0.0, 4, optsyst[0].NBelem);
	
	 for(k=0;k<piaacmc[0].nblambda;k++)
	   for(ii=0;ii<data.image[ID].md[0].size[0];ii++)
	     {
	       data.image[IDfpmresp].array.F[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + mz*data.image[ID].md[0].size[0]+ii] = data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];	
	       data.image[IDfpmresp].array.F[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + ii] -= data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];
	     }
	 save_fits("FPMresp","!FPMresp.fits");
	 sprintf(fname, "!psfi_z%02ld.fits", mz);
	 save_fits("psfi", fname);
	}
    
	 /*	 sprintf(fname, "!piaacmcfpm_z%02ld.fits", mz);
	 printf("SAVING NOW --------------- %ld ---------------->  %s\n", image_ID("piaacmcfpm"), fname);
	 mk_amph_from_complex("piaacmcfpm", "fpma", "fpmp");
	 save_fits("fpma", fname);
	 delete_image_ID("fpma");
	 delete_image_ID("fpmp");	
	 */
     
  
    focmMode = -1;
    break;


  case 12 : // search for best mask solution using FPMresp, random search 
    PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
    PIAACMCsimul_makePIAAshapes();      
    optsyst[0].FOCMASKarray[0].mode = 1; // 1-fpm
    for(k=0;k<data.image[piaacmc[0].zonezID].md[0].size[0];k++)
      data.image[piaacmc[0].zonezID].array.D[k] = 0.005e-6*k;

    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
    save_fits("imvect", "!imvect.fits");
    // imvect is 2D image,  2xnb eval pts   x   nb lambda

    mk_amph_from_complex("piaacmcfpm", "fpma", "fpmp");
    save_fits("fpma", "!piaacmcfpma.fits");
    save_fits("fpmp", "!piaacmcfpmp.fits");
    delete_image_ID("fpma");
    delete_image_ID("fpmp");	
    
    // axis 0: eval pts (ii) - size = data.image[IDfpmresp].md[0].size[0] = vsize
    // axis 1: zones (mz) - size = data.image[piaacmc[0].zonezID].md[0].size[0]+1
    // axis 3: lambda (k) - size = piaacmc[0].nblambda
    //
    // indexing :  k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*vsize + mz*vsize + ii

    ID = image_ID("imvect");
    
    val = 0.0;
    for(ii=0;ii<data.image[ID].md[0].nelement;ii++)
      {
	v1 = data.image[ID].array.F[ii];
	val += v1*v1;
      }
    val /= data.image[ID].md[0].nelement;
    val0 = val;
    valbest = val0;
    printf("val0 = %lf\n", val0);
    
    IDfpmresp = load_fits("FPMresp.fits", "FPMresp");
    vsize = data.image[IDfpmresp].md[0].size[0]; // number of eval pts x2 
    ID = create_2Dimage_ID("imvect1", vsize, piaacmc[0].nblambda);
    // measured speed: 
    //
    // [nblambda=5]
    // 4.35 kHz on single thread (without omp)
    // 9.09 kHz with omp, 8 threads
    // -> better to launch multiple instances 
    //
    // [nblambda=8]
    // 2.78 kHz  
    // 13.89 kHz with omp, 8 threads
    // -> x5 speedup

    // allocate arrays for fast routine
    
   
    fpmresp_array = data.image[IDfpmresp].array.F;
    zonez_array = data.image[piaacmc[0].zonezID].array.D;
    dphadz_array = (double*) malloc(sizeof(double)*piaacmc[0].nblambda);
    for(k=0;k<piaacmc[0].nblambda;k++)
      dphadz_array[k] = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, 1.0, optsyst[0].lambdaarray[k]);
    outtmp_array = (float*) malloc(sizeof(float)*vsize*piaacmc[0].nblambda);
	
    for(i=0;i<10000000;i++)
      {
	for(k=0;k<data.image[piaacmc[0].zonezID].md[0].size[0];k++)
	  data.image[piaacmc[0].zonezID].array.D[k] = 1.0e-6*(1.0-2.0*ran1());
	
	// test fast routine
	valtest = PIAACMCsimul_achromFPMsol_eval(fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);


	for(k=0;k<piaacmc[0].nblambda;k++)
	  {
	    ki = k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*vsize;

	    // outer zone
	    for(ii=0;ii<vsize;ii++)
	      data.image[ID].array.F[k*vsize+ii] = data.image[IDfpmresp].array.F[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*vsize+ii];
	    
	    // inner zones
	    for(mz=0;mz<data.image[piaacmc[0].zonezID].md[0].size[0];mz++) // zone
	      {
		// phase offset for this zone
		t = data.image[piaacmc[0].zonezID].array.D[mz];
		pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[k]);
		cosp = cos(pha);
		sinp = sin(pha);
		ki1 = ki + (mz+1)*vsize;
		kv = k*vsize;
		for(ii=0;ii<vsize/2;ii++)
		  {
		    ii1 = 2*ii;
		    ii2 = 2*ii+1;
		    re = data.image[IDfpmresp].array.F[ki1 + ii1];
		    im = data.image[IDfpmresp].array.F[ki1 + ii2];
		    re1 = re*cosp - im*sinp;
		    im1 = re*sinp + im*cosp;
		    data.image[ID].array.F[kv + ii1] += re1;
		    data.image[ID].array.F[kv + ii2] += im1;
		  }
	      }
	  }
	val = 0.0;
	for(ii=0;ii<data.image[ID].md[0].nelement;ii++)
	  {
	    v1 = data.image[ID].array.F[ii];
	    val += v1*v1;
	  }
	val /= data.image[ID].md[0].nelement;
	
	if(val<valbest)
	  {
	    printf("%10ld  best value = %20lf  (%20lf)\n", i, val, val0);
	    valbest = val;
	  }
	//	else
	//printf("%20lf %20lf\n", val, val0);
      }
    save_fits("imvect1", "!imvect1.fits");
    break;


 

  default :
    printERROR(__FILE__,__func__,__LINE__, "mode not recognized");
    break;
  }
  




  if(LINOPT == 1) // linear optimization
    { 

      // Compute Reference 
      PIAACMCsimul_makePIAAshapes();      
      optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 
      valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
      
      printf("Reference = %g\n", valref);
      chname_image_ID("imvect", "vecDHref");
      ID = image_ID("vecDHref");
      xsize = data.image[ID].md[0].size[0];
      ysize = data.image[ID].md[0].size[1];
      save_fits("vecDHref", "!vecDHref.fits");
      size1Dvec = data.image[ID].md[0].nelement;
      if(REGPIAASHAPES==1)
	{
	  size1Dvec += data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
	  size1Dvec += data.image[piaacmc[0].piaa1CmodesID].md[0].size[0];
	}
      
      // re-package vector into 1D array and add regularization terms
      IDm = create_2Dimage_ID("DHmask", size1Dvec, 1);
      ID1Dref = create_2Dimage_ID("vecDHref1D", size1Dvec, 1);
      
      ID = image_ID("vecDHref");
      for(ii=0;ii<data.image[ID].md[0].nelement; ii++)
	{
	  data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii]; 
	  data.image[IDm].array.F[ii] = 1.0;
	}
      if(REGPIAASHAPES == 1)
	{
	  ID = piaacmc[0].piaa0CmodesID;
	  for(jj=0;jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];jj++)
	    {
	      data.image[ID1Dref].array.F[ii] = piaa0C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa0C_regcoeff_alpha);
	      data.image[IDm].array.F[ii] = 1.0;
	      ii++;
	    }
	  
	  ID = piaacmc[0].piaa1CmodesID;
	  for(jj=0;jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0];jj++)
	    {
	      data.image[ID1Dref].array.F[ii] = piaa1C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa1C_regcoeff_alpha);
	      data.image[IDm].array.F[ii] = 1.0;
	      ii++;
	    }
	}
      delete_image_ID("vecDHref");
      
      fp = fopen("val.opt", "w");
      fclose(fp);
      
      
      //
      // LINEAR OPTIMIZATION AROUND CURRENT POINT
      //
      
      for(iter=0;iter<NBiter;iter++)
	{
	  IDmodes = create_3Dimage_ID("DHmodes", size1Dvec, 1, NBparam);
	  
	  for(i=0;i<NBparam;i++)
	    {
	      optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 
	      if(paramtype[i]==FLOAT)
		{
		  *(paramvalf[i]) += (float) paramdelta[i];
		  val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
		}
	      else
		{
		  *(paramval[i]) += paramdelta[i];
		  val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
		}
	      
	      sprintf(fname,"!imvect_%02ld.fits", i);
	      save_fits("imvect", fname);
	      ID = image_ID("imvect");
	      
	      
	      
	      // re-package vector into 1D array and add regularization terms
	      ID1D = create_2Dimage_ID("imvect1D", size1Dvec, 1);
		    
	      for(ii=0;ii<data.image[ID].md[0].nelement; ii++)
		data.image[ID1D].array.F[ii] = data.image[ID].array.F[ii]; 
	      
	      if(REGPIAASHAPES==1)
		{		    
		  ID = piaacmc[0].piaa0CmodesID;
		  for(jj=0;jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];jj++)
		    {
		      data.image[ID1D].array.F[ii] = piaa0C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa0C_regcoeff_alpha);
		      ii++;
		    }
		  
		  ID = piaacmc[0].piaa1CmodesID;
		  for(jj=0;jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0];jj++)
		    {
		      data.image[ID1D].array.F[ii] = piaa1C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa1C_regcoeff_alpha);
		      ii++;
		    }
		}
	      delete_image_ID("imvect");
	      
	      
		    
	      if(paramtype[i]==FLOAT)
		*(paramvalf[i]) -= (float) paramdelta[i];
	      else
		*(paramval[i]) -= paramdelta[i];
	      
	      
	      
	      for(ii=0;ii<data.image[ID1D].md[0].nelement;ii++)
		data.image[IDmodes].array.F[i*data.image[ID1D].md[0].nelement+ii] = (data.image[ID1D].array.F[ii] - data.image[ID1Dref].array.F[ii]);
	      
	      
	      printf("%3ld %g %g\n", i, val, valref);
	      
	      ID = create_2Dimage_ID("DHmodes2D", size1Dvec, NBparam);
	      for(ii=0;ii<data.image[IDmodes].md[0].nelement;ii++)
		data.image[ID].array.F[ii] = data.image[IDmodes].array.F[ii];
	      save_fits("DHmodes2D", "!DHmodes.fits");	      
	      delete_image_ID("DHmodes2D");
	    }
	    
	  linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 1.0e-5, "optcoeff", 0);
	  
	  
	  // delete_image_ID("vecDHref1D");
	  
	  ID = image_ID("optcoeff");

	  // do linear scan 
	  linscanOK = 1;
	  scangain = 0.0; //scanstepgain;
	  val = 100000000000.0; // big number
	  bestgain = 0.0;
	  k = 0;
	  while(linscanOK==1)
	    {
	      linoptlimflagarray[k] = 0;
	      for(i=0;i<NBparam;i++)
		{
		  paramdeltaval[i] = -scangain*data.image[ID].array.F[i]*paramdelta[i];
		  if(paramdeltaval[i]<-parammaxstep[i])
		    {
		      paramdeltaval[i] = -parammaxstep[i];
		      linoptlimflagarray[k] = 1;
		    }
		  if(paramdeltaval[i]>parammaxstep[i])
		    {
		      paramdeltaval[i] = parammaxstep[i];
		      linoptlimflagarray[k] = 1;
		    }		  
		  if(paramtype[i]==FLOAT)
		    *(paramvalf[i]) += (float) paramdeltaval[i];
		  else
		    *(paramval[i]) += paramdeltaval[i];
		}
	      valold = val;

	      val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);

	      for(i=0;i<NBparam;i++)
		{
		  if(paramtype[i]==FLOAT)
		    *(paramvalf[i]) -= (float) paramdeltaval[i];
		  else
		    *(paramval[i]) -= paramdeltaval[i];
		}
	      
	      linoptgainarray[k] = scangain;
	      linoptvalarray[k] = val;
	      k++;

	      
	      if(val<valold)
		{
		  linscanOK = 1; 
		  bestgain = scangain;
		  scangain += scanstepgain;	      
		}
	      else
		linscanOK = 0;
	      
	      if(k>90)
		linscanOK = 0;
	      scangain += scanstepgain;	      
	    }
	  NBlinoptgain = k;
	  for(k=0;k<NBlinoptgain;k++)
	    printf("%2ld %10f %20g %d\n", k, linoptgainarray[k], linoptvalarray[k], linoptlimflagarray[k]);


	  for(i=0;i<NBparam;i++)
	    {
	      paramdeltaval[i] = -bestgain*data.image[ID].array.F[i]*paramdelta[i];
	      if(paramdeltaval[i]<-parammaxstep[i])
		paramdeltaval[i] = -parammaxstep[i];
	      if(paramdeltaval[i]>parammaxstep[i])
		paramdeltaval[i] = parammaxstep[i];
	      
	      if(paramtype[i]==FLOAT)
		*(paramvalf[i]) += (float) paramdeltaval[i];
	      else
		*(paramval[i]) += paramdeltaval[i];	       	      
	    }
	  valold = val;
	  
	  val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
	  printf("gain: %lf -> val = %20g\n", bestgain, val);
	  
	  

	  r = system("cp psfi.fits psfi_ref.fits");
	  printf(" %g -> %g\n", valref, val);
	  
	  delete_image_ID("DHmodes");
	    
	  
	  ID1Dref = image_ID("vecDHref1D"); //create_2Dimage_ID("vecDHref1D", size1Dvec, 1);      
	  ID = image_ID("imvect");
	  for(ii=0;ii<data.image[ID].md[0].nelement; ii++)
	    data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii]; 
	  
	  
	  if(REGPIAASHAPES==1)
	    {
	      ID = piaacmc[0].piaa0CmodesID;
	      for(jj=0;jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];jj++)
		{
		  data.image[ID1Dref].array.F[ii] = piaa0C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa0C_regcoeff_alpha);
		  ii++;
		}
	      
	      ID = piaacmc[0].piaa1CmodesID;
	      for(jj=0;jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0];jj++)
		{
		  data.image[ID1Dref].array.F[ii] = piaa1C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa1C_regcoeff_alpha);
		  ii++;
		}
	    }
	  delete_image_ID("imvect");
	  
	  printf("Writing results\n");
	  fflush(stdout);
	  
	  ID1Dref = image_ID("vecDHref1D");
	  
	  ID = image_ID("optcoeff");
	  fp = fopen("param.opt", "w");
	  for(i=0;i<NBparam;i++)
	    {
	      if(paramtype[i]==FLOAT)
		fprintf(fp, "%5ld %20g %20g %20g %20g\n", i, *(paramvalf[i]), data.image[ID].array.F[i], data.image[ID].array.F[i]*paramdelta[i], paramdeltaval[i]);
	      else
		fprintf(fp, "%5ld %20g %20g %20g %20g\n", i, *(paramval[i]), data.image[ID].array.F[i], data.image[ID].array.F[i]*paramdelta[i], paramdeltaval[i]);
	      
	    }      
	  fclose(fp);
	  delete_image_ID("optcoeff");
	  
	  fp = fopen("val.opt", "a");
	  fprintf(fp, "%5ld %20g %20g\n", iter, val, valref);
	  fclose(fp);
	  
	  
	  PIAAsimul_savepiaacmcconf("piaacmclinopt"); // staging area
	  sprintf(command, "rsync -au --progress ./piaacmclinopt/* ./%s/", piaacmcconfdir);
	  r = system(command); 	  
	}   
    }
  




  //  PIAAsimul_savepiaacmcconf("piaacmc0");
  //  PIAAsimul_loadpiaacmcconf("piaacmc0");
  // PIAAsimul_savepiaacmcconf("piaacmc1");
  //exit(0);


  
  if(0) // Lyot mask #0 position
    {
      paramtype[NBparam] = DOUBLE;
      paramval[NBparam] = &piaacmc[0].LyotStop_zpos[0];
      paramdelta[NBparam] = 0.05;
      parammaxstep[NBparam] = 0.05;
      parammin[NBparam] = 0.0;
      parammax[NBparam] = 2.5;
      NBparam++;
    }

  if(0) // Lyot mask #1 position
    {
      paramtype[NBparam] = DOUBLE;
      paramval[NBparam] = &piaacmc[0].LyotStop_zpos[1];
      paramdelta[NBparam] = 0.05;
      parammaxstep[NBparam] = 0.05;
      parammin[NBparam] = -0.5;
      parammax[NBparam] = 0.5;
      NBparam++;
    }

  if(0) // Focal plane mask radius
    {
      paramtype[NBparam] = DOUBLE;
      paramval[NBparam] = &piaacmc[0].fpmRad;
      paramdelta[NBparam] = 1.0e-6;
      parammaxstep[NBparam] = 5.0e-6;
      parammin[NBparam] = 1.0e-6;
      parammax[NBparam] = 1.0e-4;
      NBparam++;
    }

  if(0) // Focal plane material thickness
    {
      for(k=0;k<data.image[piaacmc[0].zonezID].md[0].size[0];k++)
	{	  
	  paramtype[NBparam] = DOUBLE;
	  paramval[NBparam] = &data.image[piaacmc[0].zonezID].array.D[k];
	  paramdelta[NBparam] = 1.0e-9;
	  parammaxstep[NBparam] = 1.0e-7;
	  parammin[NBparam] = -1.0e-5;
	  parammax[NBparam] = 1.0e-5;
	  NBparam++;
	}
    }

  if(0) // Focal plane material transmission
    {
      for(k=0;k<data.image[piaacmc[0].zoneaID].md[0].size[0];k++)
	{	
	  paramtype[NBparam] = DOUBLE;
  	  paramval[NBparam] = &data.image[piaacmc[0].zoneaID].array.D[k];
	  paramdelta[NBparam] = 1.0e-4;
	  parammaxstep[NBparam] = 5.0e-2;
	  parammin[NBparam] = 1.0e-5;
	  parammax[NBparam] = 0.99;
	  NBparam++;
	}
    }

  
 

  // FOCAL PLANE MASK THICKNESS OPTIMIZATION
  if(0)
    {
      if(1)
	{
	  IDfpmresp = create_3Dimage_ID("FPMresp", xsize, ysize, data.image[piaacmc[0].zonezID].md[0].size[0]+1);
	  
	  focmMode = data.image[piaacmc[0].zonezID].md[0].size[0]+10;
	  optsyst[0].FOCMASKarray[0].mode = 1; // 1-fpm
	  val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
	  printf("val = %g\n", val);
	  ID = image_ID("imvect");
	  for(ii=0;ii<data.image[ID].md[0].nelement;ii++)
	    data.image[IDfpmresp].array.F[ii] = data.image[ID].array.F[ii];	      
	  
	  for(k=1;k<data.image[piaacmc[0].zonezID].md[0].size[0]+1;k++)
	    {
	      focmMode = k;
	      optsyst[0].FOCMASKarray[0].mode = 0; 
	      val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
	      ID = image_ID("imvect");
	      for(ii=0;ii<data.image[ID].md[0].nelement;ii++)
		{
		  data.image[IDfpmresp].array.F[k*data.image[ID].md[0].nelement+ii] = data.image[ID].array.F[ii];	
		  data.image[IDfpmresp].array.F[ii] -= data.image[ID].array.F[ii];
		}
	      save_fits("FPMresp","!FPMresp.fits");
	    }
	  focmMode = -1;
	}
      else
	IDfpmresp = load_fits("FPMresp.fits", "FPMresp");
      

      if(1)
	{
	  IDref1 = create_2Dimage_ID("vecDHref1", xsize, ysize);
	  for(ii=0;ii<xsize*ysize;ii++)
	    data.image[IDref1].array.F[ii] = data.image[IDfpmresp].array.F[ii];
	  
	  
	  for(zi=1;zi<data.image[piaacmc[0].zonezID].md[0].size[0]+1;zi++)
	    {	  
	      t = data.image[piaacmc[0].zonezID].array.D[zi-1];	  
	      a = data.image[piaacmc[0].zoneaID].array.D[zi-1];
	      
	      for(k=0;k<optsyst[0].nblambda;k++)
		{
		  dpha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[k]);	  
		  printf("ZONE %3d    t = %10.5g   lambda = %10.5g   dpha = %10.5g\n", zi, t, optsyst[0].lambdaarray[k], dpha);
		  fflush(stdout);
		  for(ii=0;ii<xsize;ii+=2)
		    {
		      re = data.image[IDfpmresp].array.F[zi*xsize*ysize+k*xsize+ii];
		      im = data.image[IDfpmresp].array.F[zi*xsize*ysize+k*xsize+ii+1];
		      amp = sqrt(re*re+im*im);
		      pha = atan2(im,re);
		      
		      amp *= a;
		      pha += dpha;
		      
		      data.image[IDref1].array.F[k*xsize+ii] += amp*cos(pha);
		      data.image[IDref1].array.F[k*xsize+ii+1] += amp*sin(pha);	  
		    }
		}
	    }
	  save_fits("vecDHref1", "!vecDHref1.fits");
	}
    }



  //
  // RANDOM SCAN AROUND CURRENT POINT
  //
  if(0)
    {
     optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm 
     if(paramtype[i]==FLOAT)
       paramrefval[i] = *(paramvalf[i]);
     else
       paramrefval[i] = *(paramval[i]);

     val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
     fp = fopen("param_rand.opt", "w");
     fprintf(fp, "%15g ", val);
     for(i=0;i<NBparam;i++)	   
       {
	 if(paramtype[i]==FLOAT)
	   fprintf(fp, " %8f", *(paramvalf[i]));
	 else
	   fprintf(fp, " %8lf", *(paramval[i]));
       }
     fprintf(fp,"\n");
     fclose(fp);

     while(1)
       {
	 for(i=0;i<NBparam;i++)
	   {
	     if(paramtype[i]==FLOAT)	   
	       *(paramvalf[i]) = (float) (parammin[i]+ran1()*(parammax[i]-parammin[i]));
	     else
	       *(paramval[i]) = (double) (parammin[i]+ran1()*(parammax[i]-parammin[i]));	 
	   }
	 val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
	 fp = fopen("param_rand.opt", "a");
	 fprintf(fp, "%15g ", val);
	 for(i=0;i<NBparam;i++)	   
	   {
	     if(paramtype[i]==FLOAT)
	       fprintf(fp, " %8f", *(paramvalf[i]));
	     else
	       fprintf(fp, " %8lf", *(paramval[i]));
	   }
	 fprintf(fp,"\n");
	 fclose(fp);
	}
    }
  



  free(piaacmc);

  return 0;
}



