#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <err.h>
#include <fcntl.h>
#include <sched.h>
#include <ncurses.h>

#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "image_filter/image_filter.h"


#include "img_reduce/img_reduce.h"

# ifdef _OPENMP
# include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
# endif


/** image analysis/reduction routines for astronomy
 * 
 * 
 */



extern DATA data;


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//



int IMG_REDUCUE_cubesimplestat_cli()
{
    if(CLI_checkarg(1,4)==0)
        IMG_REDUCUE_cubesimplestat(data.cmdargtoken[1].val.string);
    else
        return 1;

    return(0);
}




int init_img_reduce()
{

  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "image analysis for astronomy: basic routines");
  data.NBmodule++;


  strcpy(data.cmd[data.NBcmd].key,"cubesimplestat");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = IMG_REDUCUE_cubesimplestat_cli;
  strcpy(data.cmd[data.NBcmd].info,"simple data cube stats");
  strcpy(data.cmd[data.NBcmd].syntax,"<image>");
  strcpy(data.cmd[data.NBcmd].example,"cubesimplestat");
  strcpy(data.cmd[data.NBcmd].Ccall,"long IMG_REDUCUE_cubesimplestat(char *IDin_name)");
  data.NBcmd++;

  // add atexit functions here
  
  
  return 0;
}




/** compute ave, RMS
 * 
 */


long IMG_REDUCUE_cubesimplestat(char *IDin_name)
{
    long IDin;
    long xsize, ysize, zsize;
    long xysize;
    long ii, jj, kk;
    long offset;
    double tmpf;
    
    long IDave, IDrms;

    IDin = image_ID(IDin_name);

    xsize = data.image[IDin].md[0].size[0];
    ysize = data.image[IDin].md[0].size[1];
    ysize = data.image[IDin].md[0].size[2];

    xysize = xsize*ysize;


    IDave = create_2Dimage_ID("c_ave", xsize, ysize);
    IDrms = create_2Dimage_ID("c_ave", xsize, ysize);

    for(kk=jj; kk<zsize; kk++)
    {
        offset = kk*xysize;
        for(ii=0; ii<xysize; ii++)
        {
            tmpf = data.image[IDin].array.F[offset+ii];
            data.image[IDave].array.F[ii] += tmpf;
            data.image[IDrms].array.F[ii] += tmpf*tmpf;
        }
    }

    return(IDin);
}





/** this is the main routine to pre-process a cube stream of images 
 * 
 * 
 * Optional inputs:
 * 	calib_darkim  (can be single frame or cube)
 *  calib_badpix
 *  calib_flat 
 * 
 * 
 */

int IMG_REDUCE_cubeprocess(char *IDin_name)
{
	long IDin;
	long xsize, ysize, zsize;
	long xysize;
	long ii, jj;
	
	long IDflat;
	long IDbadpix;
	long IDdark;
	
	
	IDin = image_ID(IDin_name);
	
	xsize = data.image[IDin].md[0].size[0];
	ysize = data.image[IDin].md[0].size[1];
	ysize = data.image[IDin].md[0].size[2];
	
	
	
	
	return(0);
}
