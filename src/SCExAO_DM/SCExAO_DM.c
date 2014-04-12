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
#include "SCExAO_DM/SCExAO_DM.h"

extern DATA data;


#define DMSTROKE100 0.7 // um displacement for 100V
#define MAXVOLT 150
#define NBact 2500

// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


int SCExAO_DM_CombineChannels_cli()
{
  
  SCExAO_DM_CombineChannels();

  return 0;
}




int init_SCExAO_DM()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "SCExAO DM operation");
  data.NBmodule++;


  strcpy(data.cmd[data.NBcmd].key,"scexaoDMcomb");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = SCExAO_DM_CombineChannels_cli;
  strcpy(data.cmd[data.NBcmd].info,"combine channels");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"scexaoDMcomb");
  strcpy(data.cmd[data.NBcmd].Ccall,"int SCExAO_DM_CombineChannels()");
  data.NBcmd++;


 // add atexit functions here

  return 0;
}





int SCExAO_DM_disp2V(long IDdisp, long IDvolt)
{
  long ii;
  float volt;

  data.image[IDvolt].md[0].write = 1;
  for(ii=0;ii<NBact;ii++)
    {
      volt = 100.0*sqrt(data.image[IDdisp].array.F[ii]/DMSTROKE100);
      if(volt>MAXVOLT)
	volt = MAXVOLT;
      data.image[IDvolt].array.U[ii] = (unsigned short int) (volt/300.0*65536.0);
    }

  data.image[IDvolt].md[0].write = 0;
  data.image[IDvolt].md[0].cnt0++;

  

  return 0;
}



int SCExAO_DM_CombineChannels()
{
  long naxis = 2;
  long xsize = 50;
  long ysize = 50;
  long *size;
  long ch;
  char name[200];
  long *IDch;
  long NBch = 8;
  long cnt = 0;
  long long cntsumold;
  long long cntsum;

  long IDdisp;
  long IDvolt;

  size = (long*) malloc(sizeof(long)*naxis);
  IDch = (long*) malloc(sizeof(long)*NBch);
  size[0] = xsize;
  size[1] = ysize;

  printf("Initialize channels\n");  

  for(ch=0;ch<NBch;ch++)
    {
      sprintf(name, "dmdisp%ld", ch);
      printf("Channel %ld \n", ch);
      IDch[ch] = create_image_ID(name, naxis, size, FLOAT, 1, 10);
    }

  IDdisp = create_image_ID("dmdisp", naxis, size, FLOAT, 1, 10);
  IDvolt = create_image_ID("dmvolt", naxis, size, USHORT, 1, 10);

  cntsumold = 0;
  
  while(1)
    {
      usleep(1);
      cntsum = 0;
      for(ch=0;ch<NBch;ch++)
	cntsum += data.image[IDch[ch]].md[0].cnt0;
      
      if(cntsum != cntsumold)
	{
	  printf("NEW DM SHAPE %ld   %ld %ld\n", cnt, (long) cntsum, (long) cntsumold);
	  fflush(stdout);
	  cnt++;

	  copy_image_ID("dmdisp0","dmdisptmp");
	  for(ch=1;ch<NBch;ch++)
	    {
	      sprintf(name, "dmdisp%ld", ch);
	      arith_image_add_inplace("dmdisptmp",name);
	    }
	  copy_image_ID("dmdisptmp","dmdisp");	 
	  SCExAO_DM_disp2V(IDdisp, IDvolt);
	  cntsumold = cntsum;
	}
    }
    
  free(size);
  free(IDch);

  return 0;
}
