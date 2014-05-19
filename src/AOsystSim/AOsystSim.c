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
 
  long pupsize = 256;
  double puprad = 22.0;
  long IDpupa, IDpupp;
  long IDpyrpha, IDpyramp;
  double x, y;
  long ii, jj, ii1, jj1;
  double pcoeff = 0.75;
  long pupsize2;
  long ID, IDa, IDp;
  double lenssize = 70.0;

  long wfssize = 120;
  long *wfssize_array;
  long ID_wfs, ID_wfe;
  long offset, offset1;

  pupsize2 = pupsize*pupsize;


  printf("Running fake AO system simulation\n");

  dmsize = (long*) malloc(sizeof(long)*2);
  wfssize_array = (long*) malloc(sizeof(long)*2);
 
  
  // INITIALIZE
  
  // create wavefront error input
  dmsize[0] = dmxsize;
  dmsize[1] = dmysize;
  ID_wfe = create_image_ID("dm_wfe_sim", 2, dmsize, FLOAT, 1, 0);



  // create DM
  dmsize[0] = dmxsize;
  dmsize[1] = dmysize;
  ID_dm = create_image_ID("dm_sim", 2, dmsize, FLOAT, 1, 0);


  // create WFS image
  wfssize_array[0] = wfssize;
  wfssize_array[1] = wfssize;
  ID_wfs = create_image_ID("wfs_sim", 2, wfssize_array, FLOAT, 1, 0);


  
  IDpyrpha = create_2Dimage_ID("pyrpha", pupsize, pupsize);
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


  offset1 = (pupsize-dmxsize)/2;
  printf("OFFSET = %ld\n", offset);
  cnt = -1;
  printf("\n");
  while(1)
    {
      usleep(10000); // 100 Hz 
      // cnt0 = data.image[ID_dm].md[0].cnt0;
      // if(cnt0!=cnt)
      //{
	  //	  printf("COMPUTING WFS IMAGE\n");
	  //fflush(stdout);
	  
	  IDpupa = make_disk("pupa", pupsize, pupsize, 0.5*pupsize, 0.5*pupsize, puprad);
	  IDpupp = create_2Dimage_ID("pupp", pupsize, pupsize);

	  for(ii=0;ii<dmxsize;ii++)
	    for(jj=0;jj<dmysize;jj++)
	      data.image[IDpupp].array.F[(jj+offset1)*pupsize+(ii+offset1)] = data.image[ID_dm].array.F[jj*dmxsize+ii]+data.image[ID_wfe].array.F[jj*dmxsize+ii];


	  //	  save_fl_fits("pupp","!pupp.fits");

	  mk_complex_from_amph("pupa","pupp","pupc");
	  delete_image_ID("pupa");
	  delete_image_ID("pupp");
	  permut("pupc");
	  do2dfft("pupc","focc");
	  delete_image_ID("pupc");
	  permut("focc");
	  mk_amph_from_complex("focc", "foca", "focp");
	  delete_image_ID("focc");
	  IDp = image_ID("focp");
	  IDa = image_ID("foca");
	  
	  
	  for(ii=0;ii<pupsize2;ii++)
	    {
	      data.image[IDa].array.F[ii] *= data.image[IDpyramp].array.F[ii];
	      data.image[IDp].array.F[ii] += data.image[IDpyrpha].array.F[ii];
	    }
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
	  for(ii1=0;ii1<wfssize;ii1++)
	    for(jj1=0;jj1<wfssize;jj1++)
	      data.image[ID_wfs].array.F[jj1*wfssize+ii1] = data.image[ID].array.F[(jj1+offset)*pupsize+ii1+offset]*data.image[ID].array.F[(jj1+offset)*pupsize+ii1+offset];
	  data.image[ID_wfs].md[0].write = 0;
	  data.image[ID_wfs].md[0].cnt0++;
	  delete_image_ID("pupa1");


	  printf("\r%10ld       ", data.image[ID_wfs].md[0].cnt0);
	  fflush(stdout);
	  //  cnt = cnt0;
	  //	}
    }
      
  free(dmsize);
  free(wfssize_array);

  return 0;
}
