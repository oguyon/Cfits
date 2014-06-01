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

#include "AOsystSim/AOsystSim.h"

extern DATA data;


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


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




int AOsystSim_run()
{
  long dmxsize = 50;
  long dmysize = 50;
  long *dmsize;
  long twait = 1; // us
  long long cnt = 0;
  long long cnt0;
  long ID_dm;
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
  double lenssize = 120.0;

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
 
  long usleeptime = 100; // delay at each loop

  pupsize2 = pupsize*pupsize;


  printf("Running fake AO system simulation\n");

  dmsize = (long*) malloc(sizeof(long)*2);
  wfssize_array = (long*) malloc(sizeof(long)*2);
  pupsize_array = (long*) malloc(sizeof(long)*2);
   
  
  // INITIALIZE
  
  // create wavefront error input
  dmsize[0] = dmxsize;
  dmsize[1] = dmysize;
  ID_wfe = read_sharedmem_image("dmdisp1"); // turbulence channel



  IDflat = read_sharedmem_image("dmdisp0"); // flat
  data.image[IDflat].md[0].write = 1;
  for(ii=0;ii<dmxsize;ii++)
    for(jj=0;jj<dmysize;jj++)
      data.image[IDflat].array.F[jj*dmxsize+ii] = 5.5;
  data.image[IDflat].md[0].cnt0++; 
  data.image[IDflat].md[0].write = 0;


  // create_image_ID("dm_wfe_sim", 2, dmsize, FLOAT, 1, 0);  
 
  // create PSF images
  pupsize_array[0] = pupsize;
  pupsize_array[1] = pupsize;
  IDpsf = create_image_ID("aosimpsf", 2, pupsize_array, FLOAT, 1, 0);  
  IDpsfIR = create_image_ID("aosimpsfIR", 2, pupsize_array, FLOAT, 1, 0);  



  IDatmpha = read_sharedmem_image("atmwfspha");
  

  

  IDdmdisp = read_sharedmem_image("dmdisp");

  // create DM
  dmsize[0] = dmxsize;
  dmsize[1] = dmysize;
  ID_dm = read_sharedmem_image("dmdisp2"); // correction
  //create_image_ID("dm_sim", 2, dmsize, FLOAT, 1, 0);


  // create WFS image
  wfssize_array[0] = wfssize;
  wfssize_array[1] = wfssize;
  wfssize_array[2] = 3; // number of slices in rolling buffer
  ID_wfs = create_image_ID("wfs_sim", 3, wfssize_array, FLOAT, 1, 0);
  data.image[ID_wfs].md[0].cnt1 = 0; // next slice to write

  IDpuppc = create_image_ID("puppcrop", 2, dmsize, FLOAT, 1, 0);


  
  IDpyrpha = create_2Dimage_ID("pyrpha0", pupsize, pupsize);
  IDpyramp = create_2Dimage_ID("pyramp", pupsize, pupsize);
  for(ii=0;ii<pupsize;ii++)
    for(jj=0;jj<pupsize;jj++)
      {
	x = 1.0*(ii-pupsize/2);
	y = 1.0*(jj-pupsize/2);

	data.image[IDpyrpha].array.F[jj*pupsize+ii] = pcoeff*(fabs(x)+fabs(y));
	if((fabs(x)>lenssize)||(fabs(y)>lenssize))
	  data.image[IDpyramp].array.F[jj*pupsize+ii] = 0.0;
	else
	  data.image[IDpyramp].array.F[jj*pupsize+ii] = 1.0;
      }
  gauss_filter("pyrpha0", "pyrpha", 5.0, 8);
  IDpyrpha = image_ID("pyrpha");
  save_fits("pyrpha","!pyrpha.fits");

  offset1 = (pupsize-dmxsize)/2;
  printf("OFFSET = %ld\n", offset);
  cnt = -1;
  printf("\n");
  while(1)
    {
      usleep(usleeptime);
      
      // IMPORT WF ERROR
      iioffset = 6;
      jjoffset = 6;
      IDdmt = create_2Dimage_ID("dmdispt", dmxsize, dmysize);
      for(ii=0;ii<dmxsize;ii++)
	for(jj=0;jj<dmysize;jj++)
	  {
	    data.image[IDdmt].array.F[jj*dmxsize+ii] = 0.0;
	    ii1 = ii+iioffset;
	    jj1 = jj+jjoffset;
	    data.image[IDdmt].array.F[jj*dmxsize+ii] = data.image[IDatmpha].array.F[jj1*data.image[IDatmpha].md[0].size[0]+ii1]/2.0/M_PI*lambda/2.0; // um
	    
	  }
      gauss_filter("dmdispt", "dmdisptg", 4.0, 5);
      IDdmtg = image_ID("dmdisptg");
      data.image[ID_wfe].md[0].write = 1;
      for(ii=0;ii<dmxsize;ii++)
	for(jj=0;jj<dmysize;jj++)
	  data.image[ID_wfe].array.F[jj*dmxsize+ii] = data.image[IDdmt].array.F[jj*dmxsize+ii] - alphaCorr*data.image[IDdmtg].array.F[jj*dmxsize+ii];
	        
      data.image[ID_wfe].md[0].cnt0++;
      data.image[ID_wfe].md[0].write = 0;
      delete_image_ID("dmdispt");
      delete_image_ID("dmdisptg");



      // cnt0 = data.image[ID_dm].md[0].cnt0;
      // if(cnt0!=cnt)
      //{
	  //	  printf("COMPUTING WFS IMAGE\n");
	  //fflush(stdout);
	  
	  IDpupa = make_disk("pupa", pupsize, pupsize, 0.5*pupsize, 0.5*pupsize, puprad);
	  for(ii=0;ii<pupsize;ii++)
	    for(jj=0;jj<pupsize;jj++)
	      {
		x = (1.0*ii-0.5*pupsize)/puprad;
		y = (1.0*jj-0.5*pupsize)/puprad;
		r = sqrt(x*x+y*y);
		if(r<0.3)
		  data.image[IDpupa].array.F[jj*pupsize+ii] = 0.0;
		if((fabs(0.5*x+0.5*y)<0.05)||(fabs(0.5*x-0.5*y)<0.05))
		  data.image[IDpupa].array.F[jj*pupsize+ii] = 0.0;
	      }
	  save_fits("pupa", "!pupa.fits");
	  
	  IDpupp = create_2Dimage_ID("pupp", pupsize, pupsize);
	  IDpuppIR = create_2Dimage_ID("puppIR", pupsize, pupsize);

	  
	  while((data.image[IDdmdisp].md[0].write==1)||(data.image[IDpuppc].md[0].write==1))
	    usleep(100);

	  data.image[IDpuppc].md[0].write==1;

	  for(ii=0;ii<dmxsize*dmysize;ii++)
	    data.image[IDpuppc].array.F[ii] = data.image[IDdmdisp].array.F[ii];
	    
	  for(ii=0;ii<dmxsize;ii++)
	    for(jj=0;jj<dmysize;jj++)
	      data.image[IDpupp].array.F[(jj+offset1)*pupsize+(ii+offset1)] = data.image[IDpuppc].array.F[jj*dmxsize+ii]/lambda*2.0*M_PI*2.0;
	  
	  data.image[IDpuppc].md[0].cnt0++;
	  data.image[IDpuppc].md[0].write = 0;
	  
	  for(ii=0;ii<dmxsize;ii++)
	    for(jj=0;jj<dmysize;jj++)
	      data.image[IDpuppIR].array.F[(jj+offset1)*pupsize+(ii+offset1)] = (0.65/1.65)*data.image[IDpuppc].array.F[jj*dmxsize+ii]/lambda*2.0*M_PI*2.0;
	  



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
	  for(ii=0;ii<pupsize2;ii++)
	    {
	      data.image[IDpsf].array.F[ii] = data.image[IDa].array.F[ii]*data.image[IDa].array.F[ii];
	      data.image[IDa].array.F[ii] *= data.image[IDpyramp].array.F[ii];
	      data.image[IDp].array.F[ii] += data.image[IDpyrpha].array.F[ii];
	    }
	  data.image[IDpsf].md[0].cnt0++;
	  data.image[IDpsf].md[0].write = 0;


	  mk_complex_from_amph("pupa","puppIR","pupc");
	  delete_image_ID("pupa");
	  permut("pupc");
	  do2dfft("pupc","focc");
	  delete_image_ID("pupc");
	  permut("focc");
	  mk_amph_from_complex("focc", "focaIR", "focpIR");
	  
	  delete_image_ID("focc");
	  delete_image_ID("focpIR");
	  IDa = image_ID("focaIR");
	  

	  data.image[IDpsfIR].md[0].write = 1;
	  for(ii=0;ii<pupsize2;ii++)
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
	  

	  data.image[ID_wfs].md[0].write = 1;
	  moffset = data.image[ID_wfs].md[0].cnt1*wfssize*wfssize;
	  for(ii1=0;ii1<wfssize;ii1++)
	    for(jj1=0;jj1<wfssize;jj1++)	     
	      data.image[ID_wfs].array.F[moffset+jj1*wfssize+ii1] = data.image[ID].array.F[(jj1+offset)*pupsize+ii1+offset]*data.image[ID].array.F[(jj1+offset)*pupsize+ii1+offset];
	  data.image[ID_wfs].md[0].cnt1++;
	  if(data.image[ID_wfs].md[0].cnt1==data.image[ID_wfs].md[0].size[2])
	    data.image[ID_wfs].md[0].cnt1 = 0;
	  data.image[ID_wfs].md[0].cnt0++;
	  data.image[ID_wfs].md[0].write = 0;





	  delete_image_ID("pupa1");



	  printf("\r%10ld       ", data.image[ID_wfs].md[0].cnt0);
	  fflush(stdout);
	  //  cnt = cnt0;
	  //	}
    }
      
  free(dmsize);
  free(wfssize_array);
  free(pupsize_array);


  return 0;
}
