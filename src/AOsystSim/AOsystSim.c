#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "Cfits.h"
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
  long ID_dm;
  double a;
 
  long pupsize = 256;
  double puprad = 22.0;
  long IDpupa, IDpupp;
  long IDpyrpha;
  double x, y;
  long ii, jj, ii1, jj1;
  double pcoeff = 1.0;
  long pupsize2;
  long ID;
  double lenssize = 30.0;

  long wfssize = 120;
  long *wfssize_array;
  long ID_wfs;
  long offset;

  pupsize2 = pupsize*pupsize;


  printf("Running fake AO system simulation\n");

  dmsize = (long*) malloc(sizeof(long)*2);
  wfssize_array = (long*) malloc(sizeof(long)*2);
 
  
  // INITIALIZE
  
  // create DM
  dmsize[0] = dmxsize;
  dmsize[1] = dmysize;
  ID_dm = create_image_ID("dm_sim", 2, dmsize, FLOAT, 1, 0);


  // create WFS image
  wfssize_array[0] = wfssize;
  wfssize_array[1] = wfssize;
  ID_wfs = create_image_ID("wfs_sim", 2, wfssize_array, FLOAT, 1, 0);



  IDpyrpha = create_2Dimage_ID("pyrpha", pupsize, pupsize);
  for(ii=0;ii<pupsize;ii++)
    for(jj=0;jj<pupsize;jj++)
      {
	x = 1.0*(ii-pupsize/2);
	y = 1.0*(jj-pupsize/2);
	data.image[IDpyrpha].array.F[jj*pupsize+ii] = pcoeff*(fabs(x)+fabs(y));
	if(fabs(x)>lenssize)
	  data.image[IDpyrpha].array.F[jj*pupsize+ii] = 0.0;
	if(fabs(y)>lenssize)
	  data.image[IDpyrpha].array.F[jj*pupsize+ii] = 0.0;
      }

  IDpupa = make_disk("pupa", pupsize, pupsize, 0.5*pupsize, 0.5*pupsize, puprad);
  IDpupp = create_2Dimage_ID("pupp", pupsize, pupsize);
  mk_complex_from_amph("pupa","pupp","pupc");
  permut("pupc");
  do2dfft("pupc","focc");
  permut("focc");
  mk_amph_from_complex("focc", "foca", "focp");
  ID = image_ID("focp");
  for(ii=0;ii<pupsize2;ii++)
    data.image[ID].array.F[ii] += data.image[IDpyrpha].array.F[ii];
  mk_complex_from_amph("foca","focp","focc1");
  permut("focc1");
  do2dfft("focc1","pupc1");
  permut("pupc1");
  mk_amph_from_complex("pupc1", "pupa1", "pupp1");

  ID = image_ID("pupa1");
  offset = (pupsize-wfssize)/2;
  for(ii1=0;ii1<wfssize;ii1++)
    for(jj1=0;jj1<wfssize;jj1++)
      {
	data.image[ID_wfs].array.F[jj1*wfssize+ii1] = data.image[ID].array.F[(jj+offset)*pupsize+ii+offset]*data.image[ID].array.F[(jj+offset)*pupsize+ii+offset];
      }

  free(dmsize);
  free(wfssize);

  return 0;
}
