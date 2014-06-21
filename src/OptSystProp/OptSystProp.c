#include <fitsio.h> 
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate/WFpropagate.h"
#include "statistic/statistic.h"
#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"

#include "OptSystProp/OptSystProp.h"

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



int init_OptSystProp()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "Optical propagation through system");
  data.NBmodule++;
   
   
  // add atexit functions here

  return 0;

}










int OptSystProp_propagateCube(OPTSYST *optsyst, long index, char *IDin_amp_name, char *IDin_pha_name, char *IDout_amp_name, char *IDout_pha_name, double zprop)
{
  int kl;
  long ii;
  long size;
  long size2;
  long IDin_amp, IDin_pha;
  long IDc_in, IDc_out;
  long IDout_amp, IDout_pha;

  double amp, pha, re, im;



  printf("propagating by %lf m\n", zprop);

  IDin_amp = image_ID(IDin_amp_name);
  IDin_pha = image_ID(IDin_pha_name);
  size = data.image[IDin_amp].md[0].size[0];
  size2 = size*size;
  IDc_in = create_2DCimage_ID("tmppropCin", size, size);

  IDout_amp = create_3Dimage_ID(IDout_amp_name, size, size, optsyst[index].nblambda);
  IDout_pha = create_3Dimage_ID(IDout_pha_name, size, size, optsyst[index].nblambda);

  for(kl=0; kl<optsyst[0].nblambda; kl++)
    {
      printf("kl = %d / %d  %f\n", kl, optsyst[index].nblambda, optsyst[index].lambdaarray[kl]);
      for(ii=0;ii<size2;ii++)
	{
	  amp = data.image[IDin_amp].array.F[kl*size2+ii];
	  pha = data.image[IDin_pha].array.F[kl*size2+ii];
	  data.image[IDc_in].array.CF[ii].re = amp*cos(pha);
	  data.image[IDc_in].array.CF[ii].im = amp*sin(pha);
	}
      Fresnel_propagate_wavefront("tmppropCin", "tmppropCout", optsyst[index].pixscale, zprop, optsyst[index].lambdaarray[kl]);

      IDc_out = image_ID("tmppropCout");
      for(ii=0;ii<size2;ii++)
	{
	  re = data.image[IDc_out].array.CF[ii].re;
	  im = data.image[IDc_out].array.CF[ii].im;
	  amp = sqrt(re*re+im*im);
	  pha = atan2(im,re);
	  data.image[IDout_amp].array.F[kl*size2+ii] = amp;
	  data.image[IDout_pha].array.F[kl*size2+ii] = pha;
	}
      delete_image_ID("tmppropCout");
    }

  return 0;
}




int OptSystProp_run(OPTSYST *optsyst, long index, long elemstart, long elemend)
{
  long IDx, IDy, IDr, IDPA;
  double x, y;
  long IDa, IDp;
  long size;
  long nblambda;
  long size2;
  long ii, jj, k;
  long elem;
  long kl;

  char fname_pupa0[200];
  char fname_pupp0[200];
  char fname[200];

  long ID;
  double proplim = 1.0e-4;
  double total;


  float beamradpix;
  long ID0, ID1, ID2;
  long size0, size1;
  long i, j;

  char imnameamp_in[200];
  char imnamepha_in[200];
  char imnameamp_out[200];
  char imnamepha_out[200];
  long emax;
  double propdist;

  long gsize, offset, ii1, jj1;
  float re, im;
  long IDre, IDim, IDre1, IDim1;
  float *convkern;
  float val, tot;
  float u, t;

  long elemstart1 = 0;
  int elemOK;

  size = optsyst[0].size;
  size2 = size*size;
  nblambda = optsyst[0].nblambda;


  
  // create base complex amplitude
  IDa = image_ID("WFamp");
  if(IDa==-1)
    IDa = create_3Dimage_ID("WFamp", size, size, nblambda);

  IDp = image_ID("WFpha");
  if(IDp==-1)
    IDp = create_3Dimage_ID("WFpha", size, size, nblambda);

  for(ii=0;ii<size2;ii++)
    for(kl=0;kl<nblambda;kl++)
      data.image[IDa].array.F[size2*kl+ii] = 1.0;
  


  elemstart1 = 0;
  elemOK = 1;
  while(elemOK==1)
    {
      if(elemstart1==0)
	{
	  sprintf(imnameamp_in, "WFamp");
	  sprintf(imnamepha_in, "WFpha");
	}
      else
	{
	  sprintf(imnameamp_in, "WFamp_%03ld", elemstart1-1);
	  sprintf(imnamepha_in, "WFpha_%03ld", elemstart1-1);
	}      
      if(((ID1=image_ID(imnameamp_in))!=-1)&&((ID2=image_ID(imnamepha_in))!=-1)&&(elemstart1<elemstart+1))
	{
	  elemstart1++;
	  elemOK = 1;
	}
      else
	elemOK = 0;

      printf("%ld/%ld %d    %s %ld   %s %ld\n", elemstart1, elemstart, elemOK, imnameamp_in, ID1, imnamepha_in, ID2);
    }
  elemstart1--;

  printf("STARTING AT ELEMENT %ld\n", elemstart1);

  emax = elemend;
  if(emax>optsyst[0].NBelem)
    emax = optsyst[0].NBelem;

  for(elem=elemstart1; elem<emax; elem++)
    {
      if(elem==0)
	{
	  sprintf(imnameamp_in, "WFamp");
	  sprintf(imnamepha_in, "WFpha");
	  sprintf(imnameamp_out, "WFamp_000");
	  sprintf(imnamepha_out, "WFpha_000");
	  propdist = optsyst[0].elemZpos[0];
	}
      else
	{
	  sprintf(imnameamp_in, "WFamp_%03ld", elem-1);
	  sprintf(imnamepha_in, "WFpha_%03ld", elem-1);
	  sprintf(imnameamp_out, "WFamp_%03ld", elem);
	  sprintf(imnamepha_out, "WFpha_%03ld", elem);  
	  propdist = optsyst[0].elemZpos[elem]-optsyst[0].elemZpos[elem-1];
	}
      
      
      if(image_ID(imnameamp_out)!=-1)
	delete_image_ID(imnameamp_out);
      
      if(image_ID(imnamepha_out)!=-1)
	delete_image_ID(imnamepha_out);


      if( fabs(propdist)>proplim )
	{
	  printf("Propagating to element %ld  (%lf m)\n", elem,  propdist);
	  OptSystProp_propagateCube(optsyst, 0, imnameamp_in, imnamepha_in, imnameamp_out, imnamepha_out, propdist); 
	}
      else
	{
	  copy_image_ID(imnameamp_in, imnameamp_out);
	  copy_image_ID(imnamepha_in, imnamepha_out);	  
	}
      IDa = image_ID(imnameamp_out);
      IDp = image_ID(imnamepha_out);

      printf("Applying element %ld\n", elem);
      fflush(stdout);

      if(optsyst[0].elemtype[elem]==1)   // OPAQUE MASK
	{
	  printf("============= Opaque mask =================\n");
	  fflush(stdout);
	  ID = optsyst[0].elemarrayindex[elem];

	  if(ID == -1)
	    {
	      printf("ERROR: ID = -1, missing mask image\n");
	      exit(0);
	    }

	  if(data.image[ID].md[0].size[2] != nblambda)
	    {
	      for(ii=0;ii<size2;ii++)
		for(kl=0;kl<nblambda;kl++)
		  data.image[IDa].array.F[size2*kl+ii] *= data.image[ID].array.F[ii];
	    }
	  else
	    {
	      for(ii=0;ii<size2*nblambda;ii++)
		data.image[IDa].array.F[ii] *= data.image[ID].array.F[ii];	      
	    }
	}

      if(optsyst[0].elemtype[elem]==3)  // MIRROR
	{
	  printf("============= Mirror =======================\n");
	  fflush(stdout);
	  ID = optsyst[0].ASPHSURFMarray[optsyst[0].elemarrayindex[elem]].surfID;
	  printf("%ld mirror ID = %ld\n", optsyst[0].elemarrayindex[elem], ID);
	  fflush(stdout);
	  for(ii=0;ii<size2;ii++)
	    for(kl=0;kl<nblambda;kl++)
	      data.image[IDp].array.F[size2*kl+ii] -= 4.0*M_PI*data.image[ID].array.F[ii]/optsyst[0].lambdaarray[kl];
	}


      if(optsyst[0].elemtype[elem]==5)  // FOCAL PLANE MASK - MASK INPUT IS 1-MASK FOR EFFICIENT DFT
	{
	  printf("============= Focal Plane Mask ==============\n");
	  fflush(stdout);
	  // uses 1-fpm 
	  
	  // test
	  //	  save_fits(imnameamp_out, "!TESTamp.fits");
	  //exit(0);

	  ID = mk_complex_from_amph(imnameamp_out, imnamepha_out, "_WFctmp");
	  delete_image_ID(imnameamp_out);
	  delete_image_ID(imnamepha_out);

	  
	  if(optsyst[0].DFTgridpad>0)
	    {
	      // RESAMPLE ON A SPARSE GRID TO SPEED UP DFT
	      IDre = create_3Dimage_ID("dftgridre", size, size, nblambda);
	      IDim = create_3Dimage_ID("dftgridim", size, size, nblambda);
	      gsize = 2*optsyst[0].DFTgridpad+1; // grid size, odd number - this is the space between pixels
	      offset = optsyst[0].DFTgridpad; // offset from box edge to active pixel
	      
	      ID = image_ID("_WFctmp");
	      for(kl=0;kl<nblambda;kl++)
		for(ii=0;ii<size;ii++)
		  for(jj=0;jj<size;jj++)
		    {
		      re = data.image[ID].array.CF[size2*kl+jj*size+ii].re;
		      im = data.image[ID].array.CF[size2*kl+jj*size+ii].im;
		      ii1 = offset + ((long) (ii/gsize))*gsize; 
		      jj1 = offset + ((long) (jj/gsize))*gsize; 
		      if((ii1<size)&&(jj1<size))
			{
			  data.image[IDre].array.F[size2*kl+jj1*size+ii1] += re;
			  data.image[IDim].array.F[size2*kl+jj1*size+ii1] += im;
			}
		    }
	      // save_fits("dftgridre", "!dftgridre.fits");
	      // save_fits("dftgridim", "!dftgridim.fits");
	      mk_complex_from_reim("dftgridre", "dftgridim", "_WFctmpc");
	      delete_image_ID("dftgridre");
	      delete_image_ID("dftgridim");
	      
	      
	      i = optsyst[0].elemarrayindex[elem];
	      ID = optsyst[0].FOCMASKarray[i].fpmID;
	      printf("focm : %s\n", data.image[ID].md[0].name);
	      fflush(stdout);
	      fft_DFTinsertFPM("_WFctmpc", data.image[ID].md[0].name, optsyst[0].FOCMASKarray[i].zfactor, "_WFcout");
	      delete_image_ID("_WFctmpc");

	      //mk_reim_from_complex("_WFcout", "_twfre", "_twfim");
	      //save_fits("_twfre", "!_twfre.fits");
	      //save_fits("_twfim", "!_twfim.fits");
	      //
	      // INTERPOLATE SPARSE RESULT ON CONTINUOUS GRID
	      //
	      convkern = (float*) malloc(sizeof(float)*(2*gsize+1)*(2*gsize+1));
	      tot = 0.0;
	      for(i=0;i<2*gsize+1;i++)
		for(j=0;j<2*gsize+1;j++)
		  {
		    u = fabs(1.0*(i-gsize)/gsize);
		    t = fabs(1.0*(j-gsize)/gsize);
		    val = (1.0-u)*(1.0-t);
		    convkern[j*(2*gsize+1)+i] = val;
		    //printf("   %d %d %f\n", i, j, val);
		    tot += val;
		  }
	      for(i=0;i<(2*gsize+1)*(2*gsize+1);i++)
		convkern[i] *= gsize*gsize/tot;   
	      
	      
	      ID = image_ID("_WFcout");
	      IDre1 = create_3Dimage_ID("dftgridre1", size, size, nblambda);
	      IDim1 = create_3Dimage_ID("dftgridim1", size, size, nblambda);
	      for(kl=0;kl<nblambda;kl++)
		for(ii1=offset+gsize;ii1<size-gsize;ii1+=gsize)
		  for(jj1=offset+gsize;jj1<size-gsize;jj1+=gsize)
		    {
		      re = data.image[ID].array.CF[size2*kl+jj1*size+ii1].re;
		      im = data.image[ID].array.CF[size2*kl+jj1*size+ii1].im;
		      
		      for(i=0;i<2*gsize+1;i++)
			for(j=0;j<2*gsize+1;j++)
			  {
			    ii = ii1+(i-gsize);
			    jj = jj1+(j-gsize);
			    data.image[IDre1].array.F[size2*kl+jj*size+ii] += re*convkern[j*(2*gsize+1)+i];
			    data.image[IDim1].array.F[size2*kl+jj*size+ii] += im*convkern[j*(2*gsize+1)+i];
			  }
		    }
	      //save_fits("dftgridre1", "!dftgridre1.fits");
	      //save_fits("dftgridim1", "!dftgridim1.fits");	  	  
	      free(convkern);
	      delete_image_ID("_WFcout");
	      mk_complex_from_reim("dftgridre1", "dftgridim1", "_WFcout");
	      delete_image_ID("dftgridre1");
	      delete_image_ID("dftgridim1");
	    }
	  else
	    {
	      i = optsyst[0].elemarrayindex[elem];
	      ID = optsyst[0].FOCMASKarray[i].fpmID;
	      printf("focm : %s\n", data.image[ID].md[0].name);
	      fflush(stdout);
	      fft_DFTinsertFPM("_WFctmp", data.image[ID].md[0].name, optsyst[0].FOCMASKarray[i].zfactor, "_WFcout");
	    }

	  i = optsyst[0].elemarrayindex[elem];

	  if(optsyst[0].FOCMASKarray[i].mode == 1)
	    {
	      arith_image_sub_inplace("_WFctmp", "_WFcout");
	      mk_amph_from_complex("_WFctmp", imnameamp_out, imnamepha_out);	  
	    }
	  else
	    mk_amph_from_complex("_WFcout", imnameamp_out, imnamepha_out);	  


	  delete_image_ID("_WFctmp");
	  delete_image_ID("_WFcout");	  
	  //  delete_image_ID("dftgrid");
	}

      IDa = image_ID(imnameamp_out);
      optsyst[0].flux[elem] = 0.0;
      for(kl=0;kl<nblambda;kl++)
	for(ii=0;ii<size2;ii++)
	  optsyst[0].flux[elem] += data.image[IDa].array.F[kl*size2+ii]*data.image[IDa].array.F[kl*size2+ii];

      printf("Element %ld   Flux = %lf\n", elem, optsyst[0].flux[elem]/nblambda);
      

      printf("Saving intermediate plane [%ld] ... ", elem);
      fflush(stdout);

      sprintf(fname, "!WFamp_%03ld.fits", elem);
      save_fits(imnameamp_out, fname);
      sprintf(fname, "!WFpha_%03ld.fits", elem);
      save_fits(imnamepha_out, fname);
 
      printf("done\n");
      fflush(stdout);
    }
  
  if(elem==optsyst[0].NBelem) // Compute final focal plane image
    {
      printf("COMPUTING FINAL IMAGE AS FFT OF %ld\n", elem-1);
      mk_complex_from_amph(imnameamp_out, imnamepha_out, "_WFctmp");
      permut("_WFctmp");
      do2dfft("_WFctmp","psfc");
      delete_image_ID("_WFctmp");
      permut("psfc");
      mk_amph_from_complex("psfc", "psfa", "psfp");
      save_fits("psfa","!psfa.fits");
      save_fits("psfp","!psfp.fits");
      
      ID = image_ID("psfa");
      for(ii=0;ii<size2*nblambda;ii++)
	data.image[ID].array.F[ii] /= sqrt(size2*optsyst[0].flux[0]/nblambda);
      
      arith_image_mult("psfa", "psfa", "psfi");
      total = arith_image_total("psfi")/nblambda;
      printf("TOTAL = %lf\n", total);
      
      save_fits("psfi", "!psfi.fits");
    }

 
  return(0);
}

