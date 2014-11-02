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

#include <ncurses.h>

#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

 
# ifdef _OPENMP
# include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
# endif

#define STYPESIZE 10
#define SBUFFERSIZE 1000


// MEMORY MONITOR 
FILE *listim_scr_fpo;
FILE *listim_scr_fpi;
SCREEN *listim_scr; // for memory monitoring
int MEM_MONITOR = 0; // 1 if memory monitor is on
int listim_scr_wrow;
int listim_scr_wcol;


extern DATA data;


char errmsg[SBUFFERSIZE];



/** data logging of shared memory image stream
 * 
 */

struct savethreadmsg {
    char iname[100];
    char fname[200];
	int partial; // 1 if partial cube
    long cubesize; // if partial cube, this is the size of the cube
};
long tret; // thread return value






// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//

int create_image_cli()
{
  long *imsize;
  long naxis = 0;
  long i;

  
  if(CLI_checkarg(1,3)+CLI_checkarg(2,2)==0)
    {
      naxis = 0;
      imsize = (long*) malloc(sizeof(long)*5);
      i = 2;
      while(data.cmdargtoken[i].type==2)
	{
	  imsize[naxis] = data.cmdargtoken[i].val.numl;
	  naxis++;
	  i++;
	}
      switch(data.precision){
      case 0:
	create_image_ID(data.cmdargtoken[1].val.string, naxis, imsize, FLOAT, data.SHARED_DFT, data.NBKEWORD_DFT);
	break;
      case 1:
	create_image_ID(data.cmdargtoken[1].val.string, naxis, imsize, DOUBLE, data.SHARED_DFT, data.NBKEWORD_DFT);
	break;
      }
      free(imsize);
    }
  else
    return 1;
}



int image_write_keyword_L_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,2)+CLI_checkarg(4,3)==0)
    image_write_keyword_L(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.string);
  else
    return 1;
}



int image_list_keywords_cli()
{
  if(CLI_checkarg(1,4)==0)
    image_list_keywords(data.cmdargtoken[1].val.string);
  else
    return 1;
}



int read_sharedmem_image_cli()
{
  if(CLI_checkarg(1,3)==0)
    read_sharedmem_image(data.cmdargtoken[1].val.string);
  else
    return 1;
}


int create_image_shared_cli() // default precision
{
  long *imsize;
  long naxis = 0;
  long i;

  
  if(CLI_checkarg(1,3)+CLI_checkarg(2,2)==0)
    {
      naxis = 0;
      imsize = (long*) malloc(sizeof(long)*5);
      i = 2;
      while(data.cmdargtoken[i].type==2)
	{
	  imsize[naxis] = data.cmdargtoken[i].val.numl;
	  naxis++;
	  i++;
	}
      switch(data.precision){
      case 0:
	create_image_ID(data.cmdargtoken[1].val.string, naxis, imsize, FLOAT, 1, data.NBKEWORD_DFT);
	break;
      case 1:
	create_image_ID(data.cmdargtoken[1].val.string, naxis, imsize, DOUBLE, 1, data.NBKEWORD_DFT);
	break;
      }
      free(imsize);
    }
  else
    return 1;
}


int create_ushort_image_shared_cli() // default precision
{
  long *imsize;
  long naxis = 0;
  long i;

  
  if(CLI_checkarg(1,3)+CLI_checkarg(2,2)==0)
    {
      naxis = 0;
      imsize = (long*) malloc(sizeof(long)*5);
      i = 2;
      while(data.cmdargtoken[i].type==2)
	{
	  imsize[naxis] = data.cmdargtoken[i].val.numl;
	  naxis++;
	  i++;
	}
      create_image_ID(data.cmdargtoken[1].val.string, naxis, imsize, USHORT, 1, data.NBKEWORD_DFT);
	
      free(imsize);
    }
  else
    return 1;
}



int create_2Dimage_float()
{
  long *imsize;

  // CHECK ARGS
//  printf("CREATING IMAGE\n");
  imsize = (long*) malloc(sizeof(long)*2);

  imsize[0] = data.cmdargtoken[2].val.numl;
  imsize[1] = data.cmdargtoken[3].val.numl;

  create_image_ID(data.cmdargtoken[1].val.string, 2, imsize, FLOAT, data.SHARED_DFT, data.NBKEWORD_DFT);

  free(imsize);

  return 0;
}

int create_3Dimage_float()
{
  long *imsize;

  // CHECK ARGS
//  printf("CREATING 3D IMAGE\n");
  imsize = (long*) malloc(sizeof(long)*3);

  imsize[0] = data.cmdargtoken[2].val.numl;
  imsize[1] = data.cmdargtoken[3].val.numl;
  imsize[2] = data.cmdargtoken[4].val.numl;

  create_image_ID(data.cmdargtoken[1].val.string, 3, imsize, FLOAT, data.SHARED_DFT, data.NBKEWORD_DFT);

  free(imsize);

  return 0;
}


int memory_monitor_cli()
{
  memory_monitor(data.cmdargtoken[1].val.string);
  return 0;
}

int list_variable_ID_file_cli()
{
 if(CLI_checkarg(1,3)==0)
    list_variable_ID_file(data.cmdargtoken[1].val.string);
  else
    return 1;
}


int delete_image_ID_cli()
{
  long i = 1;
  printf("%ld : %d\n", i, data.cmdargtoken[i].type);
  while(data.cmdargtoken[i].type != 0)
    {
      if(data.cmdargtoken[i].type == 4)
	delete_image_ID(data.cmdargtoken[i].val.string);
      else
	printf("Image %s does not exist\n", data.cmdargtoken[i].val.string);
      i++;
    }

  return 0;
}

int copy_image_ID_cli()
{
  if(data.cmdargtoken[1].type != 4)
    {
      printf("Image %s does not exist\n", data.cmdargtoken[1].val.string);
      return -1;
    }
  
  copy_image_ID(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);

  return 0;
}

int chname_image_ID_cli()
{
  if(data.cmdargtoken[1].type != 4)
    {
      printf("Image %s does not exist\n", data.cmdargtoken[1].val.string);
      return -1;
    }
  
  chname_image_ID(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);

  return 0;
}


int mk_complex_from_reim_cli()
{
  if(data.cmdargtoken[1].type != 4)
    {
      printf("Image %s does not exist\n", data.cmdargtoken[1].val.string);
      return -1;
    }
  if(data.cmdargtoken[2].type != 4)
    {
      printf("Image %s does not exist\n", data.cmdargtoken[2].val.string);
      return -1;
    }

  mk_complex_from_reim(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string);

  return 0;
}


int mk_complex_from_amph_cli()
{
  if(data.cmdargtoken[1].type != 4)
    {
      printf("Image %s does not exist\n", data.cmdargtoken[1].val.string);
      return -1;
    }
  if(data.cmdargtoken[2].type != 4)
    {
      printf("Image %s does not exist\n", data.cmdargtoken[2].val.string);
      return -1;
    }

  mk_complex_from_amph(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string);

  return 0;
}


int mk_reim_from_complex_cli()
{
  if(data.cmdargtoken[1].type != 4)
    {
      printf("Image %s does not exist\n", data.cmdargtoken[1].val.string);
      return -1;
    }

  mk_reim_from_complex(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string);

  return 0;
}

int mk_amph_from_complex_cli()
{
  if(data.cmdargtoken[1].type != 4)
    {
      printf("Image %s does not exist\n", data.cmdargtoken[1].val.string);
      return -1;
    }

  mk_amph_from_complex(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string);

  return 0;
}



//long COREMOD_MEMORY_cp2shm(char *IDname, char *IDshmname);
int COREMOD_MEMORY_cp2shm_cli()
{
	if(CLI_checkarg(1,4)+CLI_checkarg(2,3)==0)
    {
		COREMOD_MEMORY_cp2shm(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);
      return 0;
    }
  else
    return 1; 
}



int COREMOD_MEMORY_sharedMem_2Dim_log_cli()
{
	if(CLI_checkarg(1,3)+CLI_checkarg(2,2)+CLI_checkarg(3,3)==0)
    {
		COREMOD_MEMORY_sharedMem_2Dim_log(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.string);
      return 0;
    }
  else
    return 1; 
}









int init_COREMOD_memory()
{
  strcpy(data.module[data.NBmodule].name,__FILE__);
  strcpy(data.module[data.NBmodule].info,"memory management for images and variables");
  data.NBmodule++;

  strcpy(data.cmd[data.NBcmd].key,"listim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = list_image_ID;
  strcpy(data.cmd[data.NBcmd].info,"list images in memory");
  strcpy(data.cmd[data.NBcmd].syntax,"no argument");
  strcpy(data.cmd[data.NBcmd].example,"listim");
  strcpy(data.cmd[data.NBcmd].Ccall,"int list_image_ID()");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"mmon");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = memory_monitor_cli;
  strcpy(data.cmd[data.NBcmd].info,"Monitor memory content");
  strcpy(data.cmd[data.NBcmd].syntax,"terminal tty name");
  strcpy(data.cmd[data.NBcmd].example,"mmon /dev/pts/4");
  strcpy(data.cmd[data.NBcmd].Ccall,"int memory_monitor(char *ttyname)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"listvar");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = list_variable_ID;
  strcpy(data.cmd[data.NBcmd].info,"list variables in memory");
  strcpy(data.cmd[data.NBcmd].syntax,"no argument");
  strcpy(data.cmd[data.NBcmd].example,"listvar");
  strcpy(data.cmd[data.NBcmd].Ccall,"int list_variable_ID()");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"listvarf");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = list_variable_ID_file_cli;
  strcpy(data.cmd[data.NBcmd].info,"list variables in memory, write to file");
  strcpy(data.cmd[data.NBcmd].syntax,"<file name>");
  strcpy(data.cmd[data.NBcmd].example,"listvarf var.txt");
  strcpy(data.cmd[data.NBcmd].Ccall,"int list_variable_ID_file()");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"creaim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = create_image_cli;
  strcpy(data.cmd[data.NBcmd].info,"create image, default precision");
  strcpy(data.cmd[data.NBcmd].syntax,"<name> <xsize> <ysize> <opt: zsize>");
  strcpy(data.cmd[data.NBcmd].example,"creaim imname 512 512");
  strcpy(data.cmd[data.NBcmd].Ccall,"long create_image_ID(char *name, long naxis, long *size, int atype, 0, 10)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"imwritekwL");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = image_write_keyword_L_cli;
  strcpy(data.cmd[data.NBcmd].info,"write long type keyword");
  strcpy(data.cmd[data.NBcmd].syntax,"<imname> <kname> <value [long]> <comment>");
  strcpy(data.cmd[data.NBcmd].example,"imwritekwL im1 kw2 34 my_keyword_comment");
  strcpy(data.cmd[data.NBcmd].Ccall,"long image_write_keyword_L(char *IDname, char *kname, long value, char *comment)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"imlistkw");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = image_list_keywords_cli;
  strcpy(data.cmd[data.NBcmd].info,"list image keywords");
  strcpy(data.cmd[data.NBcmd].syntax,"<imname>");
  strcpy(data.cmd[data.NBcmd].example,"imlistkw im1");
  strcpy(data.cmd[data.NBcmd].Ccall,"long image_list_keywords(char *IDname)");
  data.NBcmd++;
 
   strcpy(data.cmd[data.NBcmd].key,"readshmim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = read_sharedmem_image_cli;
  strcpy(data.cmd[data.NBcmd].info,"read shared memory image");
  strcpy(data.cmd[data.NBcmd].syntax,"<name>");
  strcpy(data.cmd[data.NBcmd].example,"readshmim im1");
  strcpy(data.cmd[data.NBcmd].Ccall,"read_sharedmem_image(char *name)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"creaimshm");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = create_image_shared_cli;
  strcpy(data.cmd[data.NBcmd].info,"create image in shared mem, default precision");
  strcpy(data.cmd[data.NBcmd].syntax,"<name> <xsize> <ysize> <opt: zsize>");
  strcpy(data.cmd[data.NBcmd].example,"creaimshm imname 512 512");
  strcpy(data.cmd[data.NBcmd].Ccall,"long create_image_ID(char *name, long naxis, long *size, int atype, 0, 10)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"creaushortimshm");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = create_ushort_image_shared_cli;
  strcpy(data.cmd[data.NBcmd].info,"create unsigned short image in shared mem");
  strcpy(data.cmd[data.NBcmd].syntax,"<name> <xsize> <ysize> <opt: zsize>");
  strcpy(data.cmd[data.NBcmd].example,"creaushortimshm imname 512 512");
  strcpy(data.cmd[data.NBcmd].Ccall,"long create_image_ID(char *name, long naxis, long *size, USHORT, 0, 10)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"crea3dim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = create_3Dimage_float;
  strcpy(data.cmd[data.NBcmd].info,"creates 3D image, single precision");
  strcpy(data.cmd[data.NBcmd].syntax,"<name> <xsize> <ysize> <zsize>");
  strcpy(data.cmd[data.NBcmd].example,"crea3dim imname 512 512 100");
  strcpy(data.cmd[data.NBcmd].Ccall,"long create_image_ID(char *name, long naxis, long *size, FLOAT, 0, 10)");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"rm");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = delete_image_ID_cli;
  strcpy(data.cmd[data.NBcmd].info,"remove image(s)");
  strcpy(data.cmd[data.NBcmd].syntax,"list of images");
  strcpy(data.cmd[data.NBcmd].example,"rm im1 im4");
  strcpy(data.cmd[data.NBcmd].Ccall,"int delete_image_ID(char* imname)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"cp");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = copy_image_ID_cli;
  strcpy(data.cmd[data.NBcmd].info,"copy image");
  strcpy(data.cmd[data.NBcmd].syntax,"source, dest");
  strcpy(data.cmd[data.NBcmd].example,"cp im1 im4");
  strcpy(data.cmd[data.NBcmd].Ccall,"long copy_image_ID(char *name, char *newname)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"mv");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = chname_image_ID_cli;
  strcpy(data.cmd[data.NBcmd].info,"change image name");
  strcpy(data.cmd[data.NBcmd].syntax,"source, dest");
  strcpy(data.cmd[data.NBcmd].example,"mv im1 im4");
  strcpy(data.cmd[data.NBcmd].Ccall,"long chname_image_ID(char *name, char *newname)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"ri2c");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = mk_complex_from_reim_cli;
  strcpy(data.cmd[data.NBcmd].info,"real, imaginary -> complex");
  strcpy(data.cmd[data.NBcmd].syntax,"real imaginary complex");
  strcpy(data.cmd[data.NBcmd].example,"ri2c imr imi imc");
  strcpy(data.cmd[data.NBcmd].Ccall,"int mk_complex_from_reim(char *re_name, char *im_name, char *out_name)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"ap2c");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = mk_complex_from_amph_cli;
  strcpy(data.cmd[data.NBcmd].info,"ampl, pha -> complex");
  strcpy(data.cmd[data.NBcmd].syntax,"ampl pha complex");
  strcpy(data.cmd[data.NBcmd].example,"ap2c ima imp imc");
  strcpy(data.cmd[data.NBcmd].Ccall,"int mk_complex_from_amph(char *re_name, char *im_name, char *out_name)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"c2ri");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = mk_reim_from_complex_cli;
  strcpy(data.cmd[data.NBcmd].info,"complex -> real, imaginary");
  strcpy(data.cmd[data.NBcmd].syntax,"complex real imaginary");
  strcpy(data.cmd[data.NBcmd].example,"c2ri imc imr imi");
  strcpy(data.cmd[data.NBcmd].Ccall,"int mk_reim_from_complex(char *re_name, char *im_name, char *out_name)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"c2ap");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = mk_amph_from_complex_cli;
  strcpy(data.cmd[data.NBcmd].info,"complex -> ampl, pha");
  strcpy(data.cmd[data.NBcmd].syntax,"complex ampl pha");
  strcpy(data.cmd[data.NBcmd].example,"c2ap imc ima imp");
  strcpy(data.cmd[data.NBcmd].Ccall,"int mk_amph_from_complex(char *re_name, char *im_name, char *out_name)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"rmall");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = clearall;
  strcpy(data.cmd[data.NBcmd].info,"remove all images");
  strcpy(data.cmd[data.NBcmd].syntax,"no argument");
  strcpy(data.cmd[data.NBcmd].example,"rmall");
  strcpy(data.cmd[data.NBcmd].Ccall,"int clearall()");
  data.NBcmd++;
 
 
  strcpy(data.cmd[data.NBcmd].key,"imcp2shm");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = COREMOD_MEMORY_cp2shm_cli;
  strcpy(data.cmd[data.NBcmd].info,"copy image ot shared memory");
  strcpy(data.cmd[data.NBcmd].syntax,"<image> <shared mem image>");
  strcpy(data.cmd[data.NBcmd].example,"imcp2shm im1 ims1");
  strcpy(data.cmd[data.NBcmd].Ccall,"long COREMOD_MEMORY_cp2shm(char *IDname, char *IDshmname)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"shmimstreamlog");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = COREMOD_MEMORY_sharedMem_2Dim_log_cli;
  strcpy(data.cmd[data.NBcmd].info,"logs shared memory stream (run in current directory)");
  strcpy(data.cmd[data.NBcmd].syntax,"<shm image> <cubesize [long]> <logdir>");
  strcpy(data.cmd[data.NBcmd].example,"shmimstreamlog wfscamim 10000 /media/data");
  strcpy(data.cmd[data.NBcmd].Ccall,"long COREMOD_MEMORY_sharedMem_2Dim_log(char *IDname, long zsize, char *logdir)");
  data.NBcmd++;
 

 
  // add atexit functions here

  return 0;
}























long compute_nb_image()
{
  long i;
  long total=0;

  for(i=0;i<data.NB_MAX_IMAGE;i++)
    {
      if(data.image[i].used==1)
	total += 1;
    }
  return(total);
}

long compute_nb_variable()
{
  long i;
  long total=0;

  for(i=0;i<data.NB_MAX_VARIABLE;i++)
    {
      if(data.variable[i].used==1)
	total += 1;
    }
  return(total);
}

long long compute_image_memory()
{
  long i;
  long long total=0;

	printf("\n");
  for(i=0;i<data.NB_MAX_IMAGE;i++)
    {
if(data.image[i].used==1)
{
		printf("image %ld  %s", i, data.image[i].md[0].name);
		fflush(stdout);
		
 		printf("used=%d ", data.image[i].used);
		fflush(stdout);

		printf("nelement=%ld ", data.image[i].md[0].nelement);
		fflush(stdout);
	
		printf("atype=%d ", data.image[i].md[0].atype);
		fflush(stdout);
		
		printf("typesize=%d ", TYPESIZE[data.image[i].md[0].atype]);
		fflush(stdout);
		
		printf("\n");
		fflush(stdout);
	}
     if(data.image[i].used==1)
		total += data.image[i].md[0].nelement*TYPESIZE[data.image[i].md[0].atype];
    }

  return(total);
}

long compute_variable_memory()
{
  long i;
  long total=0;

  for(i=0;i<data.NB_MAX_VARIABLE;i++)
    {
      total += sizeof(VARIABLE);
      if(data.variable[i].used==1)
	{
	  total += 0;
	}
    }
  return(total);
}

long image_ID(char *name) /* ID number corresponding to a name */
{
  long i,ID;
  int found;
  long tmp = 0;
  struct timespec timenow;

  i = 0;
  found = 0;
  while(found == 0)
    {
      if(data.image[i].used == 1)
	{
	  if((strncmp(name,data.image[i].md[0].name,strlen(name))==0)&&(data.image[i].md[0].name[strlen(name)]=='\0'))
	    {
	      found = 1;
	      tmp = i;
	      clock_gettime(CLOCK_REALTIME, &timenow);
	      data.image[i].md[0].last_access = 1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec;	      
	    }
	}
      i++;
      if(i == data.NB_MAX_IMAGE)
	{
	  found = 1;
	  tmp = -1;
	}}
  ID = tmp;
  
  return(tmp);    
}

long image_ID_noaccessupdate(char *name) /* ID number corresponding to a name */
{
  long i,ID;
  int found;
  long tmp = 0;

  i = 0;
  found = 0;
  while(found == 0)
    {
      if(data.image[i].used == 1)
	{
	  if((strncmp(name,data.image[i].md[0].name,strlen(name))==0)&&(data.image[i].md[0].name[strlen(name)]=='\0'))
	    {
	      found = 1;
	      tmp = i;
	    }
	}
      i++;
      if(i == data.NB_MAX_IMAGE)
	{
	  found = 1;
	  tmp = -1;
	}}
  ID = tmp;
  
  return(tmp);    
}

long variable_ID(char *name) /* ID number corresponding to a name */
{
  long i,ID;
  int found;
  long tmp = -1;

  i = 0;
  found = 0;
  while(found == 0)
    {
      if(data.variable[i].used == 1)
	{
	  if((strncmp(name,data.variable[i].name,strlen(name))==0)&&(data.variable[i].name[strlen(name)]=='\0'))
	    {
	      found = 1;
	      tmp = i;
	    }
	}
      i++;
      if(i == data.NB_MAX_VARIABLE)
	{
	  found = 1;
	  tmp = -1;
	}
    }
  ID = tmp;
  
  /*  if(tmp==-1) printf("error : no variable named \"%s\" in memory\n", name);*/
  return(tmp);    
}


long next_avail_image_ID() /* next available ID number */
{
  long i;
  long ID = -1;
  int found = 0;

# ifdef _OPENMP
#pragma omp critical
  {
#endif
  for (i=0;i<data.NB_MAX_IMAGE;i++)
    {
      if((data.image[i].used == 0)&&(found == 0))
	{
	  ID = i;
	  found = 1;
	}
    }
# ifdef _OPENMP
  }
# endif

  if(ID==-1)
    ID = data.NB_MAX_IMAGE;
  
  return(ID);
}

long next_avail_variable_ID() /* next available ID number */
{
  long i;
  long ID = -1;
  int found = 0;

  for (i=0;i<data.NB_MAX_VARIABLE;i++)
    {
      if((data.variable[i].used == 0)&&(found == 0))
	{
	  ID = i;
	  found = 1;
	}
    }
  if(ID==-1)
    {
      ID = data.NB_MAX_VARIABLE;
    }
  return(ID);
}



int delete_image_ID(char* imname) /* deletes an ID */
{
    long ID;
    char command[200];
    int r;

    ID = image_ID(imname);

    if (ID!=-1)
    {
        data.image[ID].used = 0;

        if(data.image[ID].md[0].shared == 1)
        {
            if (munmap(data.image[ID].md, data.image[ID].memsize) == -1) {
                printf("unmapping ID %ld : %p  %ld\n", ID, data.image[ID].md, data.image[ID].memsize);
                perror("Error un-mmapping the file");
            }
            close(data.image[ID].shmfd);
            data.image[ID].md = NULL;
            data.image[ID].kw = NULL;
            data.image[ID].shmfd = -1;
            data.image[ID].memsize = 0;

            sprintf(command, "rm %s/%s.im.shm", SHAREDMEMDIR, imname);
            r = system(command);
        }
        else
        {

            if(data.image[ID].md[0].atype==CHAR)
            {
                if(data.image[ID].array.C == NULL)
                {
                    printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
                    exit(0);
                }
                free(data.image[ID].array.C);
                data.image[ID].array.C = NULL;
            }
            if(data.image[ID].md[0].atype==INT)
            {
                if(data.image[ID].array.I == NULL)
                {
                    printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
                    exit(0);
                }
                free(data.image[ID].array.I);
                data.image[ID].array.I = NULL;
            }
            if(data.image[ID].md[0].atype==FLOAT)
            {
                if(data.image[ID].array.F == NULL)
                {
                    printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
                    exit(0);
                }
                free(data.image[ID].array.F);
                data.image[ID].array.F = NULL;
            }
            if(data.image[ID].md[0].atype==DOUBLE)
            {
                if(data.image[ID].array.D == NULL)
                {
                    printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
                    exit(0);
                }
                free(data.image[ID].array.D);
                data.image[ID].array.D = NULL;
            }
            if(data.image[ID].md[0].atype==COMPLEX_FLOAT)
            {
                if(data.image[ID].array.CF == NULL)
                {
                    printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
                    exit(0);
                }
                free(data.image[ID].array.CF);
                data.image[ID].array.CF = NULL;
            }
            if(data.image[ID].md[0].atype==COMPLEX_DOUBLE)
            {
                if(data.image[ID].array.CD == NULL)
                {
                    printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
                    exit(0);
                }
                free(data.image[ID].array.CD);
                data.image[ID].array.CD = NULL;
            }


            if(data.image[ID].md == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
                exit(0);
            }
            free(data.image[ID].md);
            data.image[ID].md = NULL;
        

        free(data.image[ID].kw);
        data.image[ID].kw = NULL;
		}
        /*      free(data.image[ID].size);*/
        //      data.image[ID].md[0].last_access = 0;
    }
    else
        fprintf(stderr,"%c[%d;%dm WARNING: image %s does not exist [ %s  %s  %d ] %c[%d;m\n", (char) 27, 1, 31, imname, __FILE__, __func__, __LINE__, (char) 27, 0);


    if(MEM_MONITOR == 1)
        list_image_ID_ncurses();

    return(0);
}


// delete all images with a prefix
int delete_image_ID_prefix(char *prefix)
{
    long i;

    for (i=0; i<data.NB_MAX_IMAGE; i++)
    {
        if(data.image[i].used==1)
            if((strncmp(prefix,data.image[i].md[0].name,strlen(prefix)))==0)
            {
                printf("deleting image %s\n",data.image[i].md[0].name);
                delete_image_ID(data.image[i].md[0].name);
            }
    }
    return(0);
}


int delete_variable_ID(char* varname) /* deletes a variable ID */
{
    long ID;

    ID = variable_ID(varname);
    if (ID!=-1)
    {
        data.variable[ID].used = 0;
        /*      free(data.variable[ID].name);*/
    }
    else
        fprintf(stderr,"%c[%d;%dm WARNING: variable %s does not exist [ %s  %s  %d ] %c[%d;m\n", (char) 27, 1, 31, varname, __FILE__, __func__, __LINE__, (char) 27, 0);

    return(0);
}




/* creates an image ID */
/* all images should be created by this function */
long create_image_ID(char *name, long naxis, long *size, int atype, int shared, int NBkw)
{
    long ID;
    long i,ii;
    time_t lt;
    long nelement;
    struct timespec timenow;


    size_t sharedsize = 0; // shared memory size in bytes
    int SM_fd; // shared memory file descriptor
    char SM_fname[200];
    int result;
    IMAGE_METADATA *map;
    char *mapv; // pointed cast in bytes

    int kw;
    char comment[80];
    char kname[16];

//	printf("NBkw = %ld\n", NBkw);

    ID = -1;
    if(image_ID(name)==-1)
    {
        ID = next_avail_image_ID();
        data.image[ID].used = 1;

        nelement = 1;
        for(i=0; i<naxis; i++)
            nelement*=size[i];

        // compute total size to be allocated
        if(shared==1)
        {
            sharedsize = sizeof(IMAGE_METADATA);

            if(atype==CHAR)
                sharedsize += nelement*sizeof(char);
            if(atype==INT)
                sharedsize += nelement*sizeof(int);
            if(atype==FLOAT)
                sharedsize += nelement*sizeof(float);
            if(atype==DOUBLE)
                sharedsize += nelement*sizeof(double);
            if(atype==COMPLEX_FLOAT)
                sharedsize += nelement*2*sizeof(float);
            if(atype==COMPLEX_DOUBLE)
                sharedsize += nelement*2*sizeof(double);
            if(atype==USHORT)
                sharedsize += nelement*sizeof(unsigned short int);
            if(atype==LONG)
                sharedsize += nelement*sizeof(long);


            sharedsize += NBkw*sizeof(IMAGE_KEYWORD);


            sprintf(SM_fname, "%s/%s.im.shm", SHAREDMEMDIR, name);
            SM_fd = open(SM_fname, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);

            if (SM_fd == -1) {
                perror("Error opening file for writing");
                exit(0);
            }
            data.image[ID].shmfd = SM_fd;
            data.image[ID].memsize = sharedsize;

            result = lseek(SM_fd, sharedsize-1, SEEK_SET);
            if (result == -1) {
                close(SM_fd);
                perror("Error calling lseek() to 'stretch' the file");
                exit(0);
            }

            result = write(SM_fd, "", 1);
            if (result != 1) {
                close(SM_fd);
                perror("Error writing last byte of the file");
                exit(0);
            }

            map = (IMAGE_METADATA*) mmap(0, sharedsize, PROT_READ | PROT_WRITE, MAP_SHARED, SM_fd, 0);
            if (map == MAP_FAILED) {
                close(SM_fd);
                perror("Error mmapping the file");
                exit(0);
            }

            data.image[ID].md = (IMAGE_METADATA*) map;
            data.image[ID].md[0].shared = 1;            
        }
        else
        {
            data.image[ID].md = (IMAGE_METADATA*) malloc(sizeof(IMAGE_METADATA));
            data.image[ID].md[0].shared = 0;
            data.image[ID].kw = (IMAGE_KEYWORD*) malloc(sizeof(IMAGE_KEYWORD)*NBkw);
        }

        data.image[ID].md[0].atype = atype;
        data.image[ID].md[0].naxis = naxis;
        strcpy(data.image[ID].md[0].name, name);
        for(i=0; i<naxis; i++)
            data.image[ID].md[0].size[i] = size[i];
        data.image[ID].md[0].NBkw = NBkw;




        if(atype==CHAR)
        {
            if(shared==1)
                data.image[ID].array.C = (char*) (map + sizeof(IMAGE));
            else
                data.image[ID].array.C = (char*) calloc ((size_t) nelement, sizeof(char));


            if(data.image[ID].array.C == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
                fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
                fprintf(stderr,"Image name = %s\n",name);
                fprintf(stderr,"Image size = ");
                fprintf(stderr,"%ld",size[0]);
                for(i=1; i<naxis; i++)
                    fprintf(stderr,"x%ld",size[i]);
                fprintf(stderr,"\n");
                fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(char));
                fprintf(stderr," %c[%d;m",(char) 27, 0);
                list_image_ID();
                exit(0);
            }
        }
        if(atype==INT)
        {
            if(shared==1)
                data.image[ID].array.I = (int*) (map + sizeof(IMAGE));
            else
                data.image[ID].array.I = (int*) calloc ((size_t) nelement,sizeof(int));

            if(data.image[ID].array.I == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
                fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
                fprintf(stderr,"Image name = %s\n",name);
                fprintf(stderr,"Image size = ");
                fprintf(stderr,"%ld",size[0]);
                for(i=1; i<naxis; i++)
                    fprintf(stderr,"x%ld",size[i]);
                fprintf(stderr,"\n");
                fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(int));
                fprintf(stderr," %c[%d;m",(char) 27, 0);
                list_image_ID();
                exit(0);
            }
        }
        if(atype==LONG)
        {
            if(shared==1)
                data.image[ID].array.L = (long*) (map + sizeof(IMAGE));
            else
                data.image[ID].array.L = (long*) calloc ((size_t) nelement,sizeof(long));

            if(data.image[ID].array.L == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
                fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
                fprintf(stderr,"Image name = %s\n",name);
                fprintf(stderr,"Image size = ");
                fprintf(stderr,"%ld",size[0]);
                for(i=1; i<naxis; i++)
                    fprintf(stderr,"x%ld",size[i]);
                fprintf(stderr,"\n");
                fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(long));
                fprintf(stderr," %c[%d;m",(char) 27, 0);
                list_image_ID();
                exit(0);
            }
        }
        if(atype==FLOAT)	{
            if(shared==1)
            {
                mapv = (char*) map;
                mapv += sizeof(IMAGE_METADATA);
                data.image[ID].array.F = (float*) (mapv);
                memset(data.image[ID].array.F, '\0', nelement*sizeof(float));
                mapv += sizeof(float)*nelement;
                data.image[ID].kw = (IMAGE_KEYWORD*) (mapv);            
			}
            else
                data.image[ID].array.F = (float*) calloc ((size_t) nelement, sizeof(float));

            if(data.image[ID].array.F == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
                fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
                fprintf(stderr,"Image name = %s\n",name);
                fprintf(stderr,"Image size = ");
                fprintf(stderr,"%ld",size[0]);
                for(i=1; i<naxis; i++)
                    fprintf(stderr,"x%ld",size[i]);
                fprintf(stderr,"\n");
                fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(float));
                fprintf(stderr," %c[%d;m",(char) 27, 0);
                list_image_ID();
                exit(0);
            }
        }
        if(atype==DOUBLE)
        {
            if(shared==1)
                data.image[ID].array.D = (double*) (map + sizeof(IMAGE));
            else
                data.image[ID].array.D = (double*) calloc ((size_t) nelement,sizeof(double));
            if(data.image[ID].array.D == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
                fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
                fprintf(stderr,"Image name = %s\n",name);
                fprintf(stderr,"Image size = ");
                fprintf(stderr,"%ld",size[0]);
                for(i=1; i<naxis; i++)
                    fprintf(stderr,"x%ld",size[i]);
                fprintf(stderr,"\n");
                fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(double));
                fprintf(stderr," %c[%d;m",(char) 27, 0);
                list_image_ID();
                exit(0);
            }
        }
        if(atype==COMPLEX_FLOAT)
        {
            if(shared==1)
                data.image[ID].array.CF = (complex_float*) (map + sizeof(IMAGE));
            else
                data.image[ID].array.CF = (complex_float*) calloc ((size_t) nelement,sizeof(complex_float));
            for(ii=0; ii<nelement; ii++)
            {
                data.image[ID].array.CF[ii].re = 0.0;
                data.image[ID].array.CF[ii].im = 0.0;
            }
            if(data.image[ID].array.CF == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
                fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
                fprintf(stderr,"Image name = %s\n",name);
                fprintf(stderr,"Image size = ");
                fprintf(stderr,"%ld",size[0]);
                for(i=1; i<naxis; i++)
                    fprintf(stderr,"x%ld",size[i]);
                fprintf(stderr,"\n");
                fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(float)*2);
                fprintf(stderr," %c[%d;m",(char) 27, 0);
                list_image_ID();
                exit(0);
            }
        }
        if(atype==COMPLEX_DOUBLE)
        {
            if(shared==1)
                data.image[ID].array.CD = (complex_double*) (map + sizeof(IMAGE));
            else
                data.image[ID].array.CD = (complex_double*) calloc ((size_t) nelement,sizeof(complex_double));
            if(data.image[ID].array.CD == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
                fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
                fprintf(stderr,"Image name = %s\n",name);
                fprintf(stderr,"Image size = ");
                fprintf(stderr,"%ld",size[0]);
                for(i=1; i<naxis; i++)
                    fprintf(stderr,"x%ld",size[i]);
                fprintf(stderr,"\n");
                fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(double)*2);
                fprintf(stderr," %c[%d;m",(char) 27, 0);
                list_image_ID();
                exit(0);
            }
        }
        if(atype==USHORT)
        {
            if(shared==1)
            {
                mapv = (char*) map;
                mapv += sizeof(IMAGE_METADATA);
                data.image[ID].array.U = (unsigned short*) (mapv);
                memset(data.image[ID].array.U, '\0', nelement*sizeof(unsigned short));
                mapv += sizeof(unsigned short)*nelement;
                data.image[ID].kw = (IMAGE_KEYWORD*) (mapv);
            }
            else
                data.image[ID].array.U = (unsigned short*) calloc ((size_t) nelement, sizeof(unsigned short));

            if(data.image[ID].array.U == NULL)
            {
                printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
                fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
                fprintf(stderr,"Image name = %s\n",name);
                fprintf(stderr,"Image size = ");
                fprintf(stderr,"%ld",size[0]);
                for(i=1; i<naxis; i++)
                    fprintf(stderr,"x%ld",size[i]);
                fprintf(stderr,"\n");
                fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(double)*2);
                fprintf(stderr," %c[%d;m",(char) 27, 0);
                list_image_ID();
                exit(0);
            }
        }

        clock_gettime(CLOCK_REALTIME, &timenow);
        data.image[ID].md[0].last_access = 1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec;
        data.image[ID].md[0].creation_time = data.image[ID].md[0].last_access;
		data.image[ID].md[0].write = 0;
		data.image[ID].md[0].cnt0 = 0;
		data.image[ID].md[0].cnt1 = 0;
        data.image[ID].md[0].nelement = nelement;
    }
    else
    {
        //      printf("Cannot create image : name \"%s\" already in use\n",name);
        ID = image_ID(name);

        if(data.image[ID].md[0].atype != atype)
        {
            fprintf(stderr,"%c[%d;%dm ERROR: [ %s %s %d ] %c[%d;m\n", (char) 27, 1, 31, __FILE__, __func__, __LINE__, (char) 27, 0);
            fprintf(stderr,"%c[%d;%dm Pre-existing image \"%s\" has wrong type %c[%d;m\n", (char) 27, 1, 31,name, (char) 27, 0);
            exit(0);
        }
        if(data.image[ID].md[0].naxis != naxis)
        {
            fprintf(stderr,"%c[%d;%dm ERROR: [ %s %s %d ] %c[%d;m\n", (char) 27, 1, 31, __FILE__, __func__, __LINE__, (char) 27, 0);
            fprintf(stderr,"%c[%d;%dm Pre-existing image \"%s\" has wrong naxis %c[%d;m\n", (char) 27, 1, 31,name, (char) 27, 0);
            exit(0);
        }

        for(i=0; i<naxis; i++)
            if(data.image[ID].md[0].size[i] != size[i])
            {
                fprintf(stderr,"%c[%d;%dm ERROR: [ %s %s %d ] %c[%d;m\n", (char) 27, 1, 31, __FILE__, __func__, __LINE__, (char) 27, 0);
                fprintf(stderr,"%c[%d;%dm Pre-existing image \"%s\" has wrong size %c[%d;m\n", (char) 27, 1, 31,name, (char) 27, 0);
                exit(0);
            }
    }


    // initialize keywords (test)
    for(kw=0; kw<data.image[ID].md[0].NBkw; kw++)
		data.image[ID].kw[kw].type = 'N';
/*      {
        sprintf(kname, "KEY%d", kw);
        strcpy(data.image[ID].kw[kw].name, kname);
        data.image[ID].kw[kw].type = 'D';
        data.image[ID].kw[kw].value.f.numf = 1.0*kw;
        sprintf(comment, "this is keyword %d", kw);
        strcpy(data.image[ID].kw[kw].comment, comment);
      }
    */


    if(MEM_MONITOR == 1)
        list_image_ID_ncurses();

    return(ID);
}




long image_write_keyword_L(char *IDname, char *kname, long value, char *comment)
{
	long ID;
	long kw, NBkw, kw0;	
	
	ID = image_ID(IDname);
	NBkw = data.image[ID].md[0].NBkw;
	
	kw = 0;
	while((data.image[ID].kw[kw].type!='N')&&(kw<NBkw))
		kw++;
	kw0 = kw;
	
	if(kw0==NBkw)
		{
			printf("ERROR: no available keyword entry\n");
			exit(0);
		}
	else
	{
        strcpy(data.image[ID].kw[kw].name, kname);
        data.image[ID].kw[kw].type = 'L';
        data.image[ID].kw[kw].value.numl = value;
        strcpy(data.image[ID].kw[kw].comment, comment);
	}
	
	return(kw0);
}

long image_write_keyword_D(char *IDname, char *kname, double value, char *comment)
{
	long ID;
	long kw, NBkw, kw0;	
	
	ID = image_ID(IDname);
	NBkw = data.image[ID].md[0].NBkw;
	
	kw = 0;
	while((data.image[ID].kw[kw].type!='N')&&(kw<NBkw))
		kw++;
	kw0 = kw;
	
	if(kw0==NBkw)
		{
			printf("ERROR: no available keyword entry\n");
			exit(0);
		}
	else
	{
        strcpy(data.image[ID].kw[kw].name, kname);
        data.image[ID].kw[kw].type = 'D';
        data.image[ID].kw[kw].value.numf = value;
        strcpy(data.image[ID].kw[kw].comment, comment);
	}
	
	return(kw0);
}

long image_write_keyword_S(char *IDname, char *kname, char *value, char *comment)
{
	long ID;
	long kw, NBkw, kw0;	
	
	ID = image_ID(IDname);
	NBkw = data.image[ID].md[0].NBkw;
	
	kw = 0;
	while((data.image[ID].kw[kw].type!='N')&&(kw<NBkw))
		kw++;
	kw0 = kw;
	
	if(kw0==NBkw)
		{
			printf("ERROR: no available keyword entry\n");
			exit(0);
		}
	else
	{
        strcpy(data.image[ID].kw[kw].name, kname);
        data.image[ID].kw[kw].type = 'D';
        strcpy(data.image[ID].kw[kw].value.valstr, value);
        strcpy(data.image[ID].kw[kw].comment, comment);
	}
	
	return(kw0);
}
 
 


long image_list_keywords(char *IDname)
{
	long ID;
	long kw, NBkw;	
	
	ID = image_ID(IDname);
	
	for(kw=0; kw<data.image[ID].md[0].NBkw; kw++)
		{
			if(data.image[ID].kw[kw].type=='L')
				printf("%18s  %20ld %s\n", data.image[ID].kw[kw].name, data.image[ID].kw[kw].value.numl, data.image[ID].kw[kw].comment);
			if(data.image[ID].kw[kw].type=='D')
				printf("%18s  %20lf %s\n", data.image[ID].kw[kw].name, data.image[ID].kw[kw].value.numf, data.image[ID].kw[kw].comment);
			if(data.image[ID].kw[kw].type=='S')
				printf("%18s  %20s %s\n", data.image[ID].kw[kw].name, data.image[ID].kw[kw].value.valstr, data.image[ID].kw[kw].comment);			
		//	if(data.image[ID].kw[kw].type=='N')
			//	printf("-------------\n");
		}
		
	return(ID);
}


long image_read_keyword_D(char *IDname, char *kname, double *val)
{
    long ID;
    long kw, NBkw, kw0;

    ID = image_ID(IDname);
    kw0 = -1;
    for(kw=0; kw<data.image[ID].md[0].NBkw; kw++)
    {
        if((data.image[ID].kw[kw].type=='D')&&(strncmp(kname, data.image[ID].kw[kw].name ,strlen(kname))==0))
        {
            kw0 = kw;
			*val = data.image[ID].kw[kw].value.numf;
        }
    }

    return(kw0);
}


long image_read_keyword_L(char *IDname, char *kname, long *val)
{
    long ID;
    long kw, NBkw, kw0;

    ID = image_ID(IDname);
    kw0 = -1;
    for(kw=0; kw<data.image[ID].md[0].NBkw; kw++)
    {
        if((data.image[ID].kw[kw].type=='L')&&(strncmp(kname, data.image[ID].kw[kw].name ,strlen(kname))==0))
        {
            kw0 = kw;
			*val = data.image[ID].kw[kw].value.numl;
        }
    }

    return(kw0);
}






long read_sharedmem_image(char *name)
{
  long ID = -1;
  int SM_fd;
  struct stat file_stat;
  char SM_fname[200];
  IMAGE_METADATA *map;
  char *mapv;
  int atype;
  int kw;

 
  ID = next_avail_image_ID();
  data.image[ID].used = 1;

  sprintf(SM_fname, "%s/%s.im.shm", SHAREDMEMDIR, name); 	
  printf("Importing mmap file \"%s\"\n",SM_fname);
  
  SM_fd = open(SM_fname, O_RDWR);
  if(SM_fd==-1)
    {
      ID = -1;
      printf("Cannot import file\n");
    }
  else
    {
      fstat(SM_fd, &file_stat);
      printf("File %s size: %zd\n", SM_fname, file_stat.st_size);
      
      map = (IMAGE_METADATA*) mmap(0, file_stat.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, SM_fd, 0);
      if (map == MAP_FAILED) {
	close(SM_fd);
	perror("Error mmapping the file");
	exit(0);
      }
      
      data.image[ID].memsize = file_stat.st_size;
      data.image[ID].shmfd = SM_fd;
      
      data.image[ID].md = map;
      atype = data.image[ID].md[0].atype;
      data.image[ID].md[0].shared = 1;
      
      printf("image size = %ld %ld\n", data.image[ID].md[0].size[0], data.image[ID].md[0].size[1]);
      fflush(stdout);
      
      mapv = (char*) map;
      mapv += sizeof(IMAGE_METADATA);
      
      if(atype==CHAR)
	{
	  data.image[ID].array.C = (char*) mapv;
	  mapv += sizeof(char)*data.image[ID].md[0].nelement;
	}
      if(atype==INT)
	{
	  data.image[ID].array.I = (int*) mapv;
	  mapv += sizeof(int)*data.image[ID].md[0].nelement;
	}
      if(atype==FLOAT)
	{
	  data.image[ID].array.F = (float*) mapv;
	  mapv += sizeof(float)*data.image[ID].md[0].nelement;
	}      
      if(atype==DOUBLE)
	{
	  data.image[ID].array.D = (double*) mapv;
	  mapv += sizeof(double)*data.image[ID].md[0].nelement;
	}      
      if(atype==COMPLEX_FLOAT)
	{
	  data.image[ID].array.CF = (complex_float*) mapv;
	  mapv += sizeof(complex_float)*data.image[ID].md[0].nelement;
	}      
      if(atype==COMPLEX_DOUBLE)
	{
	  data.image[ID].array.CD = (complex_double*) mapv;
	  mapv += sizeof(complex_double)*data.image[ID].md[0].nelement;
	}      
      if(atype==USHORT)
	{
	  data.image[ID].array.U = (unsigned short*) mapv;
	  mapv += sizeof(unsigned short)*data.image[ID].md[0].nelement;
	}


      data.image[ID].kw = (IMAGE_KEYWORD*) (mapv);
      
      
      for(kw=0; kw<data.image[ID].md[0].NBkw; kw++)
	{
	  if(data.image[ID].kw[kw].type == 'L')
	    printf("%d  %s %ld %s\n", kw, data.image[ID].kw[kw].name, data.image[ID].kw[kw].value.numl, data.image[ID].kw[kw].comment);
	  if(data.image[ID].kw[kw].type == 'D')
	    printf("%d  %s %lf %s\n", kw, data.image[ID].kw[kw].name, data.image[ID].kw[kw].value.numf, data.image[ID].kw[kw].comment);
	  if(data.image[ID].kw[kw].type == 'S')
	    printf("%d  %s %s %s\n", kw, data.image[ID].kw[kw].name, data.image[ID].kw[kw].value.valstr, data.image[ID].kw[kw].comment);      
	}
      
      if(MEM_MONITOR == 1)
	list_image_ID_ncurses();
    }

  return(ID);
}





long create_1Dimage_ID(char *ID_name, long xsize)
{
  long ID = -1;
  long naxis = 1;
  long naxes[1];

  naxes[0]=xsize;

  if(data.precision == 0)
    ID = create_image_ID(ID_name, naxis, naxes, 3, data.SHARED_DFT, data.NBKEWORD_DFT); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name, naxis, naxes, 4, data.SHARED_DFT, data.NBKEWORD_DFT); // double precision
 
  return(ID);
}

long create_1DCimage_ID(char *ID_name, long xsize)
{
  long ID = -1;
  long naxis=1;
  long naxes[1];

  naxes[0]=xsize;

  if(data.precision == 0)
    ID = create_image_ID(ID_name,naxis,naxes,5, data.SHARED_DFT, data.NBKEWORD_DFT); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,6, data.SHARED_DFT, data.NBKEWORD_DFT); // double precision

  return(ID);
}

long create_2Dimage_ID(char *ID_name, long xsize, long ysize)
{
  long ID = -1;
  long naxis=2;
  long naxes[2];

  naxes[0]=xsize;
  naxes[1]=ysize;
  
 // printf("Creating 2D image %s, %ld x %ld [%d %d]", ID_name, xsize, ysize, data.SHARED_DFT, data.NBKEWORD_DFT);
 // fflush(stdout);

  if(data.precision == 0)
    ID = create_image_ID(ID_name, naxis, naxes, 3, data.SHARED_DFT, data.NBKEWORD_DFT); // single precision
  else if (data.precision == 1)
    ID = create_image_ID(ID_name, naxis, naxes, 4, data.SHARED_DFT, data.NBKEWORD_DFT); // double precision
  else 
    {
      printf("Default precision (%d) invalid value: assuming single precision\n", data.precision);
      ID = create_image_ID(ID_name, naxis, naxes, 3, data.SHARED_DFT, data.NBKEWORD_DFT); // single precision
    }
    
//  printf("\n");
 // fflush(stdout);
    


  return(ID);
}

long create_2Dimagedouble_ID(char *ID_name, long xsize, long ysize)
{
  long ID = -1;
  long naxis=2;
  long naxes[2];

  naxes[0]=xsize;
  naxes[1]=ysize;

  ID = create_image_ID(ID_name,naxis,naxes,4, data.SHARED_DFT, data.NBKEWORD_DFT);

  return(ID);
}


/* 2D complex image */
long create_2DCimage_ID(char *ID_name, long xsize, long ysize)
{
  long ID = -1;
  long naxis=2;
  long naxes[2];

    naxes[0]=xsize;
  naxes[1]=ysize;

  if(data.precision == 0)
    ID = create_image_ID(ID_name,naxis,naxes, 5, data.SHARED_DFT, data.NBKEWORD_DFT); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes, 6, data.SHARED_DFT, data.NBKEWORD_DFT); // double precision

  return(ID);
}


/* 3D image, single precision */
long create_3Dimage_ID_float(char *ID_name, long xsize, long ysize, long zsize)
{
  long ID = -1;
  long naxis=3;
  long naxes[3];

  naxes[0] = xsize;
  naxes[1] = ysize;
  naxes[2] = zsize;
  
  //  printf("CREATING 3D IMAGE: %s %ld %ld %ld\n", ID_name, xsize, ysize, zsize);
  //  fflush(stdout);

  ID = create_image_ID(ID_name,naxis,naxes,3, data.SHARED_DFT, data.NBKEWORD_DFT); // single precision

  //  printf("IMAGE CREATED WITH ID = %ld\n",ID);
  //  fflush(stdout);

  return(ID);
}


/* 3D image, double precision */
long create_3Dimage_ID_double(char *ID_name, long xsize, long ysize, long zsize)
{
  long ID;
  long naxis=3;
  long naxes[3];

  naxes[0] = xsize;
  naxes[1] = ysize;
  naxes[2] = zsize;

  ID = create_image_ID(ID_name,naxis,naxes,4, data.SHARED_DFT, data.NBKEWORD_DFT); // double precision
  
  return(ID);
}



/* 3D image, default precision */
long create_3Dimage_ID(char *ID_name, long xsize, long ysize, long zsize)
{
  long ID = -1;
  long naxis=3;
  long naxes[3];

  naxes[0]=xsize;
  naxes[1]=ysize;
  naxes[2]=zsize;

  if(data.precision == 0)
    ID = create_image_ID(ID_name,naxis,naxes,3, data.SHARED_DFT, data.NBKEWORD_DFT); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,4, data.SHARED_DFT, data.NBKEWORD_DFT); // double precision

  return(ID);
}

/* 3D complex image */
long create_3DCimage_ID(char *ID_name, long xsize, long ysize, long zsize)
{
  long ID = -1;
  long naxis=3;
  long naxes[3];

  naxes[0]=xsize;
  naxes[1]=ysize;
  naxes[2]=zsize;

  if(data.precision == 0)
    ID = create_image_ID(ID_name,naxis,naxes,5, data.SHARED_DFT, data.NBKEWORD_DFT); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,6, data.SHARED_DFT, data.NBKEWORD_DFT); // double precision

  return(ID);
}


long copy_image_ID(char *name, char *newname)
{
  long ID, IDout;
  long naxis;
  long *size = NULL;
  int atype;
  long nelement;
  long i;
  
  ID = image_ID(name);
  naxis = data.image[ID].md[0].naxis;

  size = (long*) malloc(sizeof(long)*naxis);
  if(size==NULL)
    {
      printERROR(__FILE__,__func__,__LINE__,"malloc error");
      exit(0);
    }

  for(i=0;i<naxis;i++)
    size[i] = data.image[ID].md[0].size[i];
  atype  = data.image[ID].md[0].atype;

  nelement = data.image[ID].md[0].nelement;

  IDout = image_ID(newname);
  if(IDout==-1)
    {
      create_image_ID(newname,naxis,size,atype, data.SHARED_DFT, data.NBKEWORD_DFT);
      IDout = image_ID(newname);
    }
  else
    {
      // verify newname has the right size and type
      if(data.image[ID].md[0].nelement != data.image[IDout].md[0].nelement)
	{
	  fprintf(stderr,"ERROR [copy_image_ID]: images %s and %s do not have the same size\n",name,newname);
	  exit(0);
	}
      if(data.image[ID].md[0].atype!=data.image[IDout].md[0].atype)
	{
	  fprintf(stderr,"ERROR [copy_image_ID]: images %s and %s do not have the same type\n",name,newname);	
	  exit(0);
	}
    }
  data.image[IDout].md[0].write = 1;
  
  if(atype==CHAR)
    memcpy (data.image[IDout].array.C,data.image[ID].array.C, sizeof(char)*nelement);
  
  if(atype==INT)
    memcpy (data.image[IDout].array.I,data.image[ID].array.I, sizeof(int)*nelement);
  
  if(atype==FLOAT)
    memcpy (data.image[IDout].array.F,data.image[ID].array.F, sizeof(float)*nelement);
  
  if(atype==DOUBLE)
    memcpy (data.image[IDout].array.D,data.image[ID].array.D, sizeof(double)*nelement);
  
  if(atype==COMPLEX_FLOAT)
    memcpy (data.image[IDout].array.CF,data.image[ID].array.CF, sizeof(float)*2*nelement);
  
  if(atype==COMPLEX_DOUBLE)
    memcpy (data.image[IDout].array.CD,data.image[ID].array.CD, sizeof(double)*2*nelement);
  
  if(atype==USHORT)
    memcpy (data.image[IDout].array.U, data.image[ID].array.U, sizeof(double)*nelement);

  data.image[IDout].md[0].write = 0;
  data.image[IDout].md[0].cnt0++;

  free(size);

  return(IDout);
}


/* creates floating point variable */
long create_variable_ID(char *name, double value)
{
  long ID;
  long i1,i2;

  ID = -1;
  i1 = image_ID(name);
  i2 = variable_ID(name);

  if(i1!=-1)
    {
      printf("ERROR: cannot create variable \"%s\": name already used as an image\n",name);
    }
  else
    {
      if(i2!=-1)
	{
	  //	  printf("Warning : variable name \"%s\" is already in use\n",name);
	  ID = i2;
	}
      else
	ID = next_avail_variable_ID();

      data.variable[ID].used = 1;
	data.variable[ID].type = 0; /** floating point double */
      strcpy(data.variable[ID].name,name);
      data.variable[ID].value.f = value;

    }

  return(ID);
}


/* creates long variable */
long create_variable_long_ID(char *name, long value)
{
  long ID;
  long i1,i2;

  ID = -1;
  i1 = image_ID(name);
  i2 = variable_ID(name);

  if(i1!=-1)
    {
      printf("ERROR: cannot create variable \"%s\": name already used as an image\n",name);
    }
  else
    {
      if(i2!=-1)
	{
	  //	  printf("Warning : variable name \"%s\" is already in use\n",name);
	  ID = i2;
	}
      else
	ID = next_avail_variable_ID();

      data.variable[ID].used = 1;
	data.variable[ID].type = 1; /** long */
      strcpy(data.variable[ID].name,name);
      data.variable[ID].value.l = value;

    }

  return(ID);
}



/* creates long variable */
long create_variable_string_ID(char *name, char *value)
{
  long ID;
  long i1,i2;

  ID = -1;
  i1 = image_ID(name);
  i2 = variable_ID(name);

  if(i1!=-1)
    {
      printf("ERROR: cannot create variable \"%s\": name already used as an image\n",name);
    }
  else
    {
      if(i2!=-1)
	{
	  //	  printf("Warning : variable name \"%s\" is already in use\n",name);
	  ID = i2;
	}
      else
	ID = next_avail_variable_ID();

      data.variable[ID].used = 1;
	data.variable[ID].type = 2; /** string */
      strcpy(data.variable[ID].name,name);
      strcpy(data.variable[ID].value.s, value);
    }

  return(ID);
}











int init_list_image_ID_ncurses(char *termttyname)
{
  int wrow, wcol;

  listim_scr_fpi = fopen(termttyname, "r");
  listim_scr_fpo = fopen(termttyname, "w");
  listim_scr = newterm(NULL, listim_scr_fpo, listim_scr_fpi);
  
  getmaxyx(stdscr, listim_scr_wrow, listim_scr_wcol);
  start_color();
  init_pair(1, COLOR_BLACK, COLOR_WHITE); 
  init_pair(2, COLOR_BLACK, COLOR_RED);
  init_pair(3, COLOR_GREEN, COLOR_BLACK);
  init_pair(4, COLOR_RED, COLOR_BLACK);
  init_pair(5, COLOR_BLACK, COLOR_GREEN);
  init_pair(6, COLOR_CYAN, COLOR_BLACK);
  init_pair(7, COLOR_MAGENTA, COLOR_BLACK);
  init_pair(8, COLOR_BLACK, COLOR_MAGENTA);
  init_pair(9, COLOR_YELLOW, COLOR_BLACK);
  
  return 0;
}

int list_image_ID_ncurses()
{
  char str[500];
  long i, j;
  long long tmp_long;
  char type[STYPESIZE];
  int atype;
  int n;
  long long sizeb, sizeKb, sizeMb, sizeGb;

  struct timespec timenow;
  double timediff;

  clock_gettime(CLOCK_REALTIME, &timenow);

  set_term(listim_scr);

  clear();					


  sizeb = compute_image_memory();


  printw("INDEX    NAME         SIZE                    TYPE        SIZE  [percent]    LAST ACCESS\n");
  printw("\n");

  for (i=0;i<data.NB_MAX_IMAGE;i++)
    {
      if(data.image[i].used==1) 
	{
	  atype = data.image[i].md[0].atype;
	  tmp_long = ((long long) (data.image[i].md[0].nelement)) * TYPESIZE[atype];
	  
	  if(data.image[i].md[0].shared == 1)
	    printw("%4ldS", i);
	  else
	    printw("%4ld ", i);

	  if(data.image[i].md[0].shared == 1)
	    attron(A_BOLD | COLOR_PAIR(9));
	  else
	    attron(A_BOLD | COLOR_PAIR(6));
	  sprintf(str, "%10s ", data.image[i].md[0].name);
	  printw(str);

	  if(data.image[i].md[0].shared == 1)
	    attroff(A_BOLD | COLOR_PAIR(9));
	  else
	    attroff(A_BOLD | COLOR_PAIR(6));

	  sprintf(str, "[ %6ld",data.image[i].md[0].size[0]);
	  
	  for(j=1;j<data.image[i].md[0].naxis;j++)
	    {
	      sprintf(str, "%s x %6ld", str, data.image[i].md[0].size[j]);
	    }
	  sprintf(str, "%s]", str);
	  
	  printw("%-28s", str);
	  
	  attron(COLOR_PAIR(3));
	  n = 0;
	  if(atype==CHAR)
	    n = snprintf(type,STYPESIZE,"CHAR   ");	
	  if(atype==INT)
	    n = snprintf(type,STYPESIZE,"INT    ");
	  if(atype==LONG)
	    n = snprintf(type,STYPESIZE,"LONG   ");
	  if(atype==FLOAT)
	    n = snprintf(type,STYPESIZE,"FLOAT  ");
	  if(atype==DOUBLE)
	    n = snprintf(type,STYPESIZE,"DOUBLE ");
	  if(atype==COMPLEX_FLOAT)
	    n = snprintf(type,STYPESIZE,"CFLOAT ");
	  if(atype==COMPLEX_DOUBLE)
	    n = snprintf(type,STYPESIZE,"CDOUBLE");
	  printw("%7s ", type);
	  
	  attroff(COLOR_PAIR(3));
	  
	  if(n >= STYPESIZE)
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  
	  printw("%10ld Kb %6.2f   ", (long) (tmp_long/1024), (float) (100.0*tmp_long/sizeb));
	  
	  timediff = (1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec) - data.image[i].md[0].last_access;
	  
	  if(timediff<0.01)
	    {
	      attron(COLOR_PAIR(4));
	      printw("%15.9f\n", timediff);
	      attroff(COLOR_PAIR(4));
	    }
	  else
	    printw("%15.9f\n", timediff);
	}
      else
	{
	  printw("\n");
	}
    }

  sizeGb = 0;
  sizeMb = 0;
  sizeKb = 0;
  sizeb = compute_image_memory();
  
  if(sizeb>1024-1)
    {
      sizeKb = sizeb/1024;
      sizeb = sizeb-1024*sizeKb;
    }
  if(sizeKb>1024-1)
    {
      sizeMb = sizeKb/1024;
      sizeKb = sizeKb-1024*sizeMb;
    }
  if(sizeMb>1024-1)
    {
      sizeGb = sizeMb/1024;
      sizeMb = sizeMb-1024*sizeGb;
    }
  
  //attron(A_BOLD);

  sprintf(str, "%ld image(s)      ", compute_nb_image());
  if(sizeGb>0)
    sprintf(str, "%s %ld GB", str, (long) (sizeGb));

  if(sizeMb>0)
    sprintf(str, "%s %ld MB", str, (long) (sizeMb));
  
  if(sizeKb>0)
    sprintf(str, "%s %ld KB", str, (long) (sizeKb));
  
  if(sizeb>0)
    sprintf(str, "%s %ld B", str, (long) (sizeb));
  
  mvprintw(listim_scr_wrow-1, 0, "%s\n", str);
  //  attroff(A_BOLD);

  refresh();

      
  return 0;
}

void close_list_image_ID_ncurses( void )
{
  printf("Closing monitor cleanly\n");
  set_term(listim_scr);
  endwin();
  fclose(listim_scr_fpo);
  fclose(listim_scr_fpi);
  MEM_MONITOR = 0;
}


int memory_monitor(char *termttyname)
{
  if(data.Debug>0)
    printf("starting memory_monitor on \"%s\"\n", termttyname);

  MEM_MONITOR = 1;
  init_list_image_ID_ncurses(termttyname);
  list_image_ID_ncurses();
  atexit(close_list_image_ID_ncurses);

}







int list_image_ID_ofp(FILE *fo)
{
  long i,j;
  long long tmp_long;
  char type[STYPESIZE];
  int atype;
  int n;
  long long sizeb, sizeKb, sizeMb, sizeGb;
  char str[500];
  struct timespec timenow;
  double timediff;

	printf("--- test 00 ---\n");
	fflush(stdout);

	sizeb = compute_image_memory();

	printf("--- test 01 ---\n");
	fflush(stdout);

  clock_gettime(CLOCK_REALTIME, &timenow);

  fprintf(fo, "\n");
  fprintf(fo, "INDEX    NAME         SIZE                    TYPE        SIZE  [percent]    LAST ACCESS\n");
	fflush(stdout);
  fprintf(fo, "\n");

  for (i=0;i<data.NB_MAX_IMAGE;i++)
    if(data.image[i].used==1) 
      {
	  atype = data.image[i].md[0].atype;
	  tmp_long = ((long long) (data.image[i].md[0].nelement)) * TYPESIZE[atype];
	  
	  if(data.image[i].md[0].shared==1)
	    fprintf(fo, "%4ld %c[%d;%dm%10s%c[%d;m ",i, (char) 27, 1, 34, data.image[i].md[0].name, (char) 27, 0);
	  else
	    fprintf(fo, "%4ld %c[%d;%dm%10s%c[%d;m ",i, (char) 27, 1, 33, data.image[i].md[0].name, (char) 27, 0);
//fprintf(fo, "%s", str);
	  
	  sprintf(str, "[ %6ld",data.image[i].md[0].size[0]);
	  
	  for(j=1;j<data.image[i].md[0].naxis;j++)
	    {
	      sprintf(str, "%s x %6ld", str, data.image[i].md[0].size[j]);
	    }
	  sprintf(str, "%s]", str);
	  
	  fprintf(fo, "%-28s", str);
	  

	  n = 0;
	  if(atype==CHAR)
	    n = snprintf(type,STYPESIZE,"CHAR");	
	  if(atype==INT)
	    n = snprintf(type,STYPESIZE,"INT");
	  if(atype==LONG)
	    n = snprintf(type,STYPESIZE,"LONG");
	  if(atype==FLOAT)
	    n = snprintf(type,STYPESIZE,"FLOAT");
	  if(atype==DOUBLE)
	    n = snprintf(type,STYPESIZE,"DOUBLE");
	  if(atype==COMPLEX_FLOAT)
	    n = snprintf(type,STYPESIZE,"CFLOAT");
	  if(atype==COMPLEX_DOUBLE)
	    n = snprintf(type,STYPESIZE,"CDOUBLE");
	  if(atype==USHORT)
	    n = snprintf(type,STYPESIZE,"USHORT");
	  fprintf(fo, "%7s ", type);
	  
	  
	  if(n >= STYPESIZE)
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  
	  fprintf(fo, "%10ld Kb %6.2f   ", (long) (tmp_long/1024), (float) (100.0*tmp_long/sizeb));
	  
	  timediff = (1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec) - data.image[i].md[0].last_access;
	  
	  fprintf(fo, "%15.9f\n", timediff);
	  fflush(stdout);
	}
  fprintf(fo, "\n");

 
  sizeGb = 0;
  sizeMb = 0;
  sizeKb = 0;
  sizeb = compute_image_memory();
  
  if(sizeb>1024-1)
    {
      sizeKb = sizeb/1024;
      sizeb = sizeb-1024*sizeKb;
    }
  if(sizeKb>1024-1)
    {
      sizeMb = sizeKb/1024;
      sizeKb = sizeKb-1024*sizeMb;
    }
  if(sizeMb>1024-1)
    {
      sizeGb = sizeMb/1024;
      sizeMb = sizeMb-1024*sizeGb;
    }

  fprintf(fo, "%ld image(s)   ",compute_nb_image());
  if(sizeGb>0)
    fprintf(fo, " %ld Gb",(long) (sizeGb));
  if(sizeMb>0)
    fprintf(fo, " %ld Mb",(long) (sizeMb));
  if(sizeKb>0)
    fprintf(fo, " %ld Kb",(long) (sizeKb));
   if(sizeb>0)
     fprintf(fo, " %ld",(long) (sizeb));
   fprintf(fo, "\n\n");

   fflush(fo);
   
   return(0);
}

int list_image_ID_ofp_simple(FILE *fo)
{
  long i,j;
  long long tmp_long;
  int atype;

  for (i=0;i<data.NB_MAX_IMAGE;i++)
    if(data.image[i].used==1) 
      {
	  atype = data.image[i].md[0].atype;
	  tmp_long = ((long long) (data.image[i].md[0].nelement)) * TYPESIZE[atype];
	  
		fprintf(fo, "%20s %d %ld %d %4ld", data.image[i].md[0].name, atype, data.image[i].md[0].naxis, data.image[i].md[0].shared, data.image[i].md[0].size[0]);
	  
	  for(j=1;j<data.image[i].md[0].naxis;j++)
	      fprintf(fo, " %4ld", data.image[i].md[0].size[j]);
	  fprintf(fo, "\n");
	}
  fprintf(fo, "\n");
   
   return(0);
}



int list_image_ID()
{
  list_image_ID_ofp(stdout);
  return 0;
}



/* list all images in memory
   output is written in ASCII file
   only basic info is listed 
   image name
   number of axis
   size
   type
 */

int list_image_ID_file(char *fname)
{
  FILE *fp;
  long i,j;
  int atype;
  char type[STYPESIZE];
  int n;
  
  fp = fopen(fname,"w");
  if(fp == NULL)
    {
      n = snprintf(errmsg,SBUFFERSIZE,"Cannot create file %s",fname);
      if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
      printERROR(__FILE__,__func__,__LINE__,errmsg);
      exit(0);
    }

  for (i=0;i<data.NB_MAX_IMAGE;i++)
    if(data.image[i].used == 1) 
      {
	atype = data.image[i].md[0].atype;
	fprintf(fp,"%ld %s",i, data.image[i].md[0].name);
	fprintf(fp," %ld",data.image[i].md[0].naxis);
	for(j=0;j<data.image[i].md[0].naxis;j++)
	  fprintf(fp," %ld",data.image[i].md[0].size[j]);
	
	n = 0;
	if(atype==CHAR)
	  n = snprintf(type,STYPESIZE,"CHAR");
	if(atype==INT)
	  n = snprintf(type,STYPESIZE,"INT");
	if(atype==LONG)
	  n = snprintf(type,STYPESIZE,"LONG");
	if(atype==FLOAT)
	  n = snprintf(type,STYPESIZE,"FLOAT");
	if(atype==DOUBLE)
	  n = snprintf(type,STYPESIZE,"DOUBLE");
	if(atype==COMPLEX_FLOAT)
	  n = snprintf(type,STYPESIZE,"CFLOAT");
	if(atype==COMPLEX_DOUBLE)
	  n = snprintf(type,STYPESIZE,"CDOUBLE");
	
	if(n >= STYPESIZE) 
	  printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

	fprintf(fp," %s\n",type);
      }
  fclose(fp);

  return(0);
}


int list_variable_ID()
{
  long i;

  for (i=0;i<data.NB_MAX_VARIABLE;i++)
    if(data.variable[i].used == 1) 
      printf("%4ld %16s %25.18g\n",i, data.variable[i].name,data.variable[i].value.f);
  
  return(0);
}

int list_variable_ID_file(char *fname)
{
  long i;
	FILE *fp;

	fp = fopen(fname, "w");
  for (i=0;i<data.NB_MAX_VARIABLE;i++)
    if(data.variable[i].used == 1) 
      fprintf(fp, "%s=%.18g\n",data.variable[i].name, data.variable[i].value.f);

  fclose(fp);
  
  return(0);
}



long chname_image_ID(char *ID_name, char *new_name)
{
  long ID;
  
  ID=-1;
  if((image_ID(new_name)==-1)&&(variable_ID(new_name)==-1))
    {  
      ID = image_ID(ID_name);
      strcpy(data.image[ID].md[0].name, new_name);
      //      if ( Debug > 0 ) { printf("change image name %s -> %s\n",ID_name,new_name);}
    }
  else
    printf("Cannot change name %s -> %s : new name already in use\n", ID_name, new_name);

  if(MEM_MONITOR==1)
    list_image_ID_ncurses();

  return(ID);
}

int mk_complex_from_reim(char *re_name, char *im_name, char *out_name)
{
  long IDre,IDim,IDout;
  long *naxes = NULL;
  long naxis;
  long nelement;
  long ii;
  long i;
  int n;
  int atype_re, atype_im, atype_out;

  IDre = image_ID(re_name);
  IDim = image_ID(im_name);
  
  atype_re = data.image[IDre].md[0].atype;
  atype_im = data.image[IDim].md[0].atype;
  naxis = data.image[IDre].md[0].naxis;

  naxes = (long *) malloc(sizeof(long)*naxis);
  if(naxes==NULL)
    {
      printERROR(__FILE__,__func__,__LINE__,"malloc error");
      exit(0);
    }

  for(i=0;i<naxis;i++)
    naxes[i] = data.image[IDre].md[0].size[i];
  nelement = data.image[IDre].md[0].nelement;


  if((atype_re==FLOAT)&&(atype_im==FLOAT))
    {
      atype_out = COMPLEX_FLOAT;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out, data.SHARED_DFT, data.NBKEWORD_DFT);
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CF[ii].re = data.image[IDre].array.F[ii];
	  data.image[IDout].array.CF[ii].im = data.image[IDim].array.F[ii];
	}  
    }
  else if((atype_re==FLOAT)&&(atype_im==DOUBLE))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out, data.SHARED_DFT, data.NBKEWORD_DFT);
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CD[ii].re = data.image[IDre].array.F[ii];
	  data.image[IDout].array.CD[ii].im = data.image[IDim].array.D[ii];
	}  
    }
  else if((atype_re==DOUBLE)&&(atype_im==FLOAT))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out, data.SHARED_DFT, data.NBKEWORD_DFT);
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CD[ii].re = data.image[IDre].array.D[ii];
	  data.image[IDout].array.CD[ii].im = data.image[IDim].array.F[ii];
	}   
    }
  else if((atype_re==DOUBLE)&&(atype_im==DOUBLE))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out, data.SHARED_DFT, data.NBKEWORD_DFT);
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CD[ii].re = data.image[IDre].array.D[ii];
	  data.image[IDout].array.CD[ii].im = data.image[IDim].array.D[ii];
	}   
    }
  else
    {
      n = snprintf(errmsg,SBUFFERSIZE,"Wrong image type(s)\n");
      if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

      printERROR(__FILE__,__func__,__LINE__,errmsg);
      exit(0);
    }
  // Note: openMP doesn't help here
  
  free(naxes);

  return(0);
}

int mk_complex_from_amph(char *am_name, char *ph_name, char *out_name)
{
  long IDam,IDph,IDout;
  long naxes[3];
  long naxis;
  long nelement;
  long ii;
  long i;
  int atype_am, atype_ph, atype_out;
  int n;

  IDam = image_ID(am_name);
  IDph = image_ID(ph_name);
  atype_am = data.image[IDam].md[0].atype;
  atype_ph = data.image[IDph].md[0].atype;

  naxis = data.image[IDam].md[0].naxis;
  for(i=0;i<naxis;i++)
    naxes[i] = data.image[IDam].md[0].size[i];
  nelement = data.image[IDam].md[0].nelement;

  if((atype_am==FLOAT)&&(atype_ph==FLOAT))
    {
      atype_out = COMPLEX_FLOAT;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out, data.SHARED_DFT, data.NBKEWORD_DFT);
  # ifdef _OPENMP
  #pragma omp parallel if (nelement>OMP_NELEMENT_LIMIT)
  {
  #pragma omp for
  # endif
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CF[ii].re = data.image[IDam].array.F[ii]*((float) cos(data.image[IDph].array.F[ii]));
	  data.image[IDout].array.CF[ii].im = data.image[IDam].array.F[ii]*((float) sin(data.image[IDph].array.F[ii]));
	}  
  # ifdef _OPENMP
  }
  # endif
    }
  else if((atype_am==FLOAT)&&(atype_ph==DOUBLE))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out, data.SHARED_DFT, data.NBKEWORD_DFT);
  # ifdef _OPENMP
  #pragma omp parallel if (nelement>OMP_NELEMENT_LIMIT)
  {
  #pragma omp for
  # endif
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CD[ii].re = data.image[IDam].array.F[ii]*cos(data.image[IDph].array.D[ii]);
	  data.image[IDout].array.CD[ii].im = data.image[IDam].array.F[ii]*sin(data.image[IDph].array.D[ii]);
	}  
  # ifdef _OPENMP
  }
  # endif
    }
  else if((atype_am==DOUBLE)&&(atype_ph==FLOAT))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out, data.SHARED_DFT, data.NBKEWORD_DFT);
  # ifdef _OPENMP
  #pragma omp parallel if (nelement>OMP_NELEMENT_LIMIT)
  {
  #pragma omp for
  # endif
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CD[ii].re = data.image[IDam].array.D[ii]*cos(data.image[IDph].array.F[ii]);
	  data.image[IDout].array.CD[ii].im = data.image[IDam].array.D[ii]*sin(data.image[IDph].array.F[ii]);
	}   
  # ifdef _OPENMP
  }
  # endif
    }
  else if((atype_am==DOUBLE)&&(atype_ph==DOUBLE))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out, data.SHARED_DFT, data.NBKEWORD_DFT);
  # ifdef _OPENMP
  #pragma omp parallel if (nelement>OMP_NELEMENT_LIMIT)
  {
  #pragma omp for
  # endif
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CD[ii].re = data.image[IDam].array.D[ii]*cos(data.image[IDph].array.D[ii]);
	  data.image[IDout].array.CD[ii].im = data.image[IDam].array.D[ii]*sin(data.image[IDph].array.D[ii]);
	}   
  # ifdef _OPENMP
  }
  # endif
    }
  else
    {
      n = snprintf(errmsg,SBUFFERSIZE,"Wrong image type(s)\n");
      if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
      
      printERROR(__FILE__,__func__,__LINE__,errmsg);
      exit(0);
    }

  return(0);
}

int mk_reim_from_complex(char *in_name, char *re_name, char *im_name)
{
  long IDre,IDim,IDin;
  long naxes[3];
  long naxis;
  long nelement;
  long ii;
  long i;
  long atype;
  int n;

  IDin = image_ID(in_name);
  atype = data.image[IDin].md[0].atype;
  naxis = data.image[IDin].md[0].naxis;
  for(i=0;i<naxis;i++)
    naxes[i] = data.image[IDin].md[0].size[i];
  nelement = data.image[IDin].md[0].nelement;

  if(atype == COMPLEX_FLOAT) // single precision
    {
      IDre = create_image_ID(re_name,naxis,naxes,FLOAT, data.SHARED_DFT, data.NBKEWORD_DFT);
      IDim = create_image_ID(im_name,naxis,naxes,FLOAT, data.SHARED_DFT, data.NBKEWORD_DFT);
  
  # ifdef _OPENMP
  #pragma omp parallel if (nelement>OMP_NELEMENT_LIMIT)
  {
  #pragma omp for
  # endif
  for(ii=0;ii<nelement;ii++)
    {
      data.image[IDre].array.F[ii] = data.image[IDin].array.CF[ii].re;
      data.image[IDim].array.F[ii] = data.image[IDin].array.CF[ii].im;
    }
  # ifdef _OPENMP
  }
  # endif
    }
  else if(atype==COMPLEX_DOUBLE) // double precision
    {
      IDre = create_image_ID(re_name,naxis,naxes,DOUBLE, data.SHARED_DFT, data.NBKEWORD_DFT);
      IDim = create_image_ID(im_name,naxis,naxes,DOUBLE, data.SHARED_DFT, data.NBKEWORD_DFT);
  
  # ifdef _OPENMP
  #pragma omp parallel if (nelement>OMP_NELEMENT_LIMIT)
  {
  #pragma omp for
  # endif
  for(ii=0;ii<nelement;ii++)
    {
      data.image[IDre].array.D[ii] = data.image[IDin].array.CD[ii].re;
      data.image[IDim].array.D[ii] = data.image[IDin].array.CD[ii].im;
    }
  # ifdef _OPENMP
  }
  # endif
    }
   else
    {
      n = snprintf(errmsg,SBUFFERSIZE,"Wrong image type(s)\n");
      if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

      printERROR(__FILE__,__func__,__LINE__,errmsg);
      exit(0);
    }
 

  return(0);
}

int mk_amph_from_complex(char *in_name, char *am_name, char *ph_name)
{
  long IDam,IDph,IDin;
  long naxes[3];
  long naxis;
  long nelement;
  long ii;
  long i;
  float amp_f,pha_f;
  double amp_d,pha_d;
  int atype;
  int n;

  IDin = image_ID(in_name);
  atype = data.image[IDin].md[0].atype;
  naxis = data.image[IDin].md[0].naxis;

  for(i=0;i<naxis;i++)
    naxes[i] = data.image[IDin].md[0].size[i];
  nelement = data.image[IDin].md[0].nelement;

  if(atype==COMPLEX_FLOAT) // single precision
    {
      IDam = create_image_ID(am_name,naxis,naxes,FLOAT, data.SHARED_DFT, data.NBKEWORD_DFT);
      IDph = create_image_ID(ph_name,naxis,naxes,FLOAT, data.SHARED_DFT, data.NBKEWORD_DFT);
      
  # ifdef _OPENMP
  #pragma omp parallel if (nelement>OMP_NELEMENT_LIMIT) private(ii,amp_f,pha_f)
  {
  #pragma omp for
  # endif
  for(ii=0;ii<nelement;ii++)
    {
      amp_f = (float) sqrt(data.image[IDin].array.CF[ii].re*data.image[IDin].array.CF[ii].re + data.image[IDin].array.CF[ii].im*data.image[IDin].array.CF[ii].im);
      pha_f = (float) atan2(data.image[IDin].array.CF[ii].im,data.image[IDin].array.CF[ii].re);
      data.image[IDam].array.F[ii] = amp_f;
      data.image[IDph].array.F[ii] = pha_f;
    }
  # ifdef _OPENMP
  }
  # endif
    }
  else if(atype==COMPLEX_DOUBLE) // double precision
    {
      IDam = create_image_ID(am_name,naxis,naxes,DOUBLE, data.SHARED_DFT, data.NBKEWORD_DFT);
      IDph = create_image_ID(ph_name,naxis,naxes,DOUBLE, data.SHARED_DFT, data.NBKEWORD_DFT);
      
  # ifdef _OPENMP
  #pragma omp parallel if (nelement>OMP_NELEMENT_LIMIT) private(ii,amp_d,pha_d)
  {
  #pragma omp for
  # endif
  for(ii=0;ii<nelement;ii++)
    {
      amp_d = sqrt(data.image[IDin].array.CD[ii].re*data.image[IDin].array.CD[ii].re + data.image[IDin].array.CD[ii].im*data.image[IDin].array.CD[ii].im);
      pha_d = atan2(data.image[IDin].array.CD[ii].im,data.image[IDin].array.CD[ii].re);
      data.image[IDam].array.D[ii] = amp_d;
      data.image[IDph].array.D[ii] = pha_d;
    }
  # ifdef _OPENMP
  }
  # endif
    }
  else
    {
      n = snprintf(errmsg,SBUFFERSIZE,"Wrong image type(s)\n");
      if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
      printERROR(__FILE__,__func__,__LINE__,errmsg);
      exit(0);
    }

  return(0);
}

int mk_reim_from_amph(char *am_name, char *ph_name, char *re_out_name, char *im_out_name)
{
  mk_complex_from_amph(am_name,ph_name,"Ctmp");
  mk_reim_from_complex("Ctmp",re_out_name,im_out_name);
  delete_image_ID("Ctmp");

  return(0);
}

int mk_amph_from_reim(char *re_name, char *im_name, char *am_out_name, char *ph_out_name)
{
  mk_complex_from_reim(re_name,im_name,"Ctmp");
  mk_amph_from_complex("Ctmp",am_out_name, ph_out_name);
  delete_image_ID("Ctmp");

  return(0);
}

int clearall()
{
  long ID;

  for(ID=0;ID<data.NB_MAX_IMAGE;ID++)
    {
      if(data.image[ID].used==1)
	delete_image_ID(data.image[ID].md[0].name);
    }
  for(ID=0;ID<data.NB_MAX_VARIABLE;ID++)
    {
      if(data.variable[ID].used==1)
	delete_variable_ID(data.variable[ID].name);
    }

  return(0);
}

int check_2Dsize(char *ID_name, long xsize, long ysize)
{
  int value;
  long ID;

  value = 1;
  ID=image_ID(ID_name);
  if(data.image[ID].md[0].naxis!=2)
    value=0;
  if(value==1)
    {
      if(data.image[ID].md[0].size[0]!=xsize)
	value = 0;
      if(data.image[ID].md[0].size[1]!=ysize)
	value = 0;
    }
  
  return(value);
}

int check_3Dsize(char *ID_name, long xsize, long ysize, long zsize)
{
  int value;
  long ID;

  value = 1;
  ID=image_ID(ID_name);
  if(data.image[ID].md[0].naxis!=3)
    {
      /*      printf("Wrong naxis : %ld - should be 3\n",data.image[ID].md[0].naxis);*/
      value = 0;
    }
  if(value==1)
    {
      if(data.image[ID].md[0].size[0]!=xsize)
	{
	  /*	  printf("Wrong xsize : %ld - should be %ld\n",data.image[ID].md[0].size[0],xsize);*/
	  value = 0;
	}
      if(data.image[ID].md[0].size[1]!=ysize)
	{
	  /*	  printf("Wrong ysize : %ld - should be %ld\n",data.image[ID].md[0].size[1],ysize);*/
	  value = 0;
	}
      if(data.image[ID].md[0].size[2]!=zsize)
	{
	  /*	  printf("Wrong zsize : %ld - should be %ld\n",data.image[ID].md[0].size[2],zsize);*/
	  value = 0;
	}
    }
  /*  printf("CHECK = %d\n",value);*/

  return(value);
}

int rotate_cube(char *ID_name, char *ID_out_name, int orientation)
{
    /* 0 is from x axis */
    /* 1 is from y axis */
    long ID,IDout;
    long xsize,ysize,zsize;
    long xsize1,ysize1,zsize1;
    long ii,jj,kk;
    int atype;
    int n;

    ID = image_ID(ID_name);
    atype = data.image[ID].md[0].atype;

    if(data.image[ID].md[0].naxis!=3)
    {
        n = snprintf(errmsg,SBUFFERSIZE,"Wrong naxis : %ld - should be 3\n",data.image[ID].md[0].naxis);
        if(n >= SBUFFERSIZE)
            printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

        printERROR(__FILE__,__func__,__LINE__,errmsg);
        exit(0);
    }
    xsize = data.image[ID].md[0].size[0];
    ysize = data.image[ID].md[0].size[1];
    zsize = data.image[ID].md[0].size[2];

    if(atype==FLOAT) // single precision
    {
        if(orientation==0)
        {
            xsize1 = zsize;
            ysize1 = ysize;
            zsize1 = xsize;
            IDout = create_3Dimage_ID_float(ID_out_name,xsize1,ysize1,zsize1);
            for(ii=0; ii<xsize1; ii++)
                for(jj=0; jj<ysize1; jj++)
                    for(kk=0; kk<zsize1; kk++)
                        data.image[IDout].array.F[kk*ysize1*xsize1+jj*xsize1+ii] = data.image[ID].array.F[ii*xsize*ysize+jj*xsize+kk];
        }
        else
        {
            xsize1 = xsize;
            ysize1 = zsize;
            zsize1 = ysize;
            IDout = create_3Dimage_ID_float(ID_out_name,xsize1,ysize1,zsize1);
            for(ii=0; ii<xsize1; ii++)
                for(jj=0; jj<ysize1; jj++)
                    for(kk=0; kk<zsize1; kk++)
                        data.image[IDout].array.F[kk*ysize1*xsize1+jj*xsize1+ii] = data.image[ID].array.F[jj*xsize*ysize+kk*xsize+ii];
        }
    }
    else if(atype==DOUBLE)
    {
        if(orientation==0)
        {
            xsize1 = zsize;
            ysize1 = ysize;
            zsize1 = xsize;
            IDout = create_3Dimage_ID_double(ID_out_name,xsize1,ysize1,zsize1);
            for(ii=0; ii<xsize1; ii++)
                for(jj=0; jj<ysize1; jj++)
                    for(kk=0; kk<zsize1; kk++)
                        data.image[IDout].array.D[kk*ysize1*xsize1+jj*xsize1+ii] = data.image[ID].array.D[ii*xsize*ysize+jj*xsize+kk];
        }
        else
        {
            xsize1 = xsize;
            ysize1 = zsize;
            zsize1 = ysize;
            IDout = create_3Dimage_ID_double(ID_out_name,xsize1,ysize1,zsize1);
            for(ii=0; ii<xsize1; ii++)
                for(jj=0; jj<ysize1; jj++)
                    for(kk=0; kk<zsize1; kk++)
                        data.image[IDout].array.D[kk*ysize1*xsize1+jj*xsize1+ii] = data.image[ID].array.D[jj*xsize*ysize+kk*xsize+ii];
        }
    }
    else
    {
        n = snprintf(errmsg,SBUFFERSIZE,"Wrong image type(s)\n");
        if(n >= SBUFFERSIZE)
            printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

        printERROR(__FILE__,__func__,__LINE__,errmsg);
        exit(0);
    }

    return(0);
}



/// ---------------------------------------- LOGGING FUNCTIONS --------------------------------



/// creates logshimconf shared memory and loads it
LOGSHIM_CONF* COREMOD_MEMORY_logshim_create_SHMconf(char *logshimname)
{
    int SM_fd;
    size_t sharedsize = 0; // shared memory size in bytes
	char SM_fname[200];
	int result;
	LOGSHIM_CONF *map;
	
	sharedsize = sizeof(LOGSHIM_CONF);
	
    sprintf(SM_fname, "%s/%s.logshimconf.shm", SHAREDMEMDIR, logshimname);
    SM_fd = open(SM_fname, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);

    if (SM_fd == -1) {
        perror("Error opening file for writing");
        exit(0);
    }
  
    result = lseek(SM_fd, sharedsize-1, SEEK_SET);
    if (result == -1) {
        close(SM_fd);
        perror("Error calling lseek() to 'stretch' the file");
        exit(0);
    }

    result = write(SM_fd, "", 1);
    if (result != 1) {
        close(SM_fd);
        perror("Error writing last byte of the file");
        exit(0);
    }

    map = (LOGSHIM_CONF*) mmap(0, sharedsize, PROT_READ | PROT_WRITE, MAP_SHARED, SM_fd, 0);
    if (map == MAP_FAILED) {
        close(SM_fd);
        perror("Error mmapping the file");
        exit(0);
    }

	map[0].on = 0;
	map[0].cnt = 0;
	map[0].filecnt = 0;
	map[0].interval = 1;
	map[0].logexit = 0;
	strcpy(map[0].fname, SM_fname);
	
    return(map);
}








void *save_fits_function( void *ptr )
{
    long ID;
    struct savethreadmsg *tmsg = malloc(sizeof(struct savethreadmsg));
    long *imsizearray;
    long xsize, ysize;
    int atype;
    long IDc;
    long framesize; // in bytes
    char *ptr0; // source
    char *ptr1; // destination


    imsizearray = (long*) malloc(sizeof(long)*3);

    tmsg = (struct savethreadmsg*) ptr;
   // printf("THREAD : SAVING  %s -> %s \n", tmsg->iname, tmsg->fname);
//fflush(stdout);
    if(tmsg->partial==0) // full image
        save_fits(tmsg->iname, tmsg->fname);
    else
    {
  //      printf("Saving partial image (name = %s   zsize = %ld)\n", tmsg->iname, tmsg->cubesize);
	
	//	list_image_ID();
 
        ID = image_ID(tmsg->iname);
        atype = data.image[ID].md[0].atype;
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];

		//printf("step00\n");
		//fflush(stdout);
		
        imsizearray[0] = xsize;
        imsizearray[1] = ysize;
        imsizearray[2] = tmsg->cubesize;

	

        IDc = create_image_ID("tmpsavecube", 3, imsizearray, atype, 0, 1);
        
       // list_image_ID();
        
        switch ( atype ) {
        case CHAR:
            framesize = sizeof(char)*xsize*ysize;
            break;
        case INT:
            framesize = sizeof(int)*xsize*ysize;
            break;
        case FLOAT:
            framesize = sizeof(float)*xsize*ysize;
            break;
        case DOUBLE:
            framesize = sizeof(double)*xsize*ysize;
            break;
        case USHORT:
            framesize = sizeof(unsigned short)*xsize*ysize;
            break;

        default:
            printf("ERROR: WRONG DATA TYPE\n");
            exit(0);
            break;
        }

        ptr0 = (char*) data.image[ID].array.F;  // source
        ptr1 = (char*) data.image[IDc].array.F; // destination


        memcpy((void *) ptr1, (void *) ptr0, framesize*tmsg->cubesize);
		save_fits("tmpsavecube", tmsg->fname);
        delete_image_ID("tmpsavecube");
    }


//    printf(" DONE\n");
//fflush(stdout);

    ID = image_ID(tmsg->iname);
    tret = ID;
    free(imsizearray);
    pthread_exit(&tret);
}


/** copy an image to shared memory 
 * 
 * 
 */
long COREMOD_MEMORY_cp2shm(char *IDname, char *IDshmname)
{
	long ID;
	long IDshm;
	long atype;
	long naxis;
	long *sizearray;
	char *ptr1;
	char *ptr2;
	long k;
	
	
	ID = image_ID(IDname);
	naxis = data.image[ID].md[0].naxis;
	
	sizearray = (long*) malloc(sizeof(long)*naxis);
	atype = data.image[ID].md[0].atype;
	for(k=0;k<naxis;k++)
		sizearray[k] = data.image[ID].md[0].size[k];

				
	IDshm = create_image_ID(IDshmname, naxis, sizearray, atype, 1, 0);
	free(sizearray);

	//data.image[IDshm].md[0].nelement = data.image[ID].md[0].nelement;
	printf("======= %ld %ld ============\n", data.image[ID].md[0].nelement, data.image[IDshm].md[0].nelement);

	switch (atype) {
		case FLOAT :
			ptr1 = (char*) data.image[ID].array.F;
			ptr2 = (char*) data.image[IDshm].array.F;
			memcpy(ptr2, ptr1, sizeof(float)*data.image[ID].md[0].nelement);
			break;
		case DOUBLE :
			ptr1 = (char*) data.image[ID].array.D;
			ptr2 = (char*) data.image[IDshm].array.D;
			memcpy(ptr2, ptr1, sizeof(float)*data.image[ID].md[0].nelement);
			break;
		case USHORT :
			ptr1 = (char*) data.image[ID].array.U;
			ptr2 = (char*) data.image[IDshm].array.U;
			memcpy(ptr2, ptr1, sizeof(unsigned short)*data.image[ID].md[0].nelement);
			break;
		default :
			printf("data type not supported\n");
			break;
	}
	
	return(0);
}


/** logs a shared memory stream onto disk
 *
 * uses data cube buffer to store frames
 *
 */
long COREMOD_MEMORY_sharedMem_2Dim_log(char *IDname, long zsize, char *logdir)
{
    long ID;
    long xsize, ysize;
    long ii;
    long IDb, IDb0, IDb1;
    long index = 0;
    long cnt = -1;
    int buffer;
    int atype;
    long *imsizearray;
    char fname[200];
    char iname[200];
    time_t t;
    struct tm *uttime;
    struct timespec *thetime = (struct timespec *)malloc(sizeof(struct timespec));
    long kw;

    char *ptr0; // source
    char *ptr1; // destination
    long framesize; // in bytes

    FILE *fp;
    char fname_asciilog[200];

    pthread_t thread_savefits;
    int tOK = 0;
    int iret_savefits;
    //	char tmessage[500];
    struct savethreadmsg *tmsg = malloc(sizeof(struct savethreadmsg));

    long fnb = 0;
    long NBfiles = -1; // run forever
    long long cntwait;
    long waitdelayus = 50;
    long long cntwaitlim = 200000; // 10 sec
    int wOK;
    int noframe;

	int is3Dcube = 0; // this is a rolling buffer


	LOGSHIM_CONF* logshimconf;


	logshimconf = COREMOD_MEMORY_logshim_create_SHMconf(IDname);

	
	logshimconf[0].on = 1;
	logshimconf[0].cnt = 0;
	logshimconf[0].filecnt = 0;
	logshimconf[0].logexit = 0;
	logshimconf[0].interval = 1;
	
	

    imsizearray = (long*) malloc(sizeof(long)*3);
	
	

    read_sharedmem_image(IDname);
    ID = image_ID(IDname);
    atype = data.image[ID].md[0].atype;
    xsize = data.image[ID].md[0].size[0];
    ysize = data.image[ID].md[0].size[1];
	
	if(data.image[ID].md[0].naxis==3)
		is3Dcube = 1;

    /** create the 2 buffers */

    imsizearray[0] = xsize;
    imsizearray[1] = ysize;
    imsizearray[2] = zsize;

    IDb0 = create_image_ID("logbuff0", 3, imsizearray, atype, 0, 1);
    IDb1 = create_image_ID("logbuff1", 3, imsizearray, atype, 0, 1);

    IDb = IDb0;

    switch ( atype ) {
    case CHAR:
        framesize = sizeof(char)*xsize*ysize;
        break;
    case INT:
        framesize = sizeof(int)*xsize*ysize;
        break;
    case FLOAT:
        framesize = sizeof(float)*xsize*ysize;
        break;
    case DOUBLE:
        framesize = sizeof(double)*xsize*ysize;
        break;
    case USHORT:
        framesize = sizeof(unsigned short)*xsize*ysize;
        break;

    default:
        printf("ERROR: WRONG DATA TYPE\n");
        exit(0);
        break;
    }

    cnt = data.image[ID].md[0].cnt0;

    buffer = 0;
    index = 0;
    while( (logshimconf[0].filecnt != NBfiles) && (logshimconf[0].logexit==0) )
    {
        cntwait = 0;
        noframe = 0;
        wOK = 1;
        // printf("Entering wait loop   index = %ld %d\n", index, noframe);
        while((cnt==data.image[ID].md[0].cnt0)&&(wOK==1))
        {
            usleep(10);
            cntwait++;
            if(cntwait>cntwaitlim) // save current cube
            {
                strcpy(tmsg->iname, iname);
                strcpy(tmsg->fname, fname);
                tmsg->partial = 1; // partial cube
                tmsg->cubesize = index;
                wOK=0;
                if(index==0)
                    noframe = 1;
                else
                    noframe = 0;
                //		printf("Exiting wait loop - %ld %d\n", index, noframe);
            }
        }
        //		printf("Entering main loop   index = %ld %d  [%d]\n", index, noframe, wOK);


        if(index==0)
        {
            /// measure time
            t = time(NULL);
            uttime = gmtime(&t);
            clock_gettime(CLOCK_REALTIME, thetime);
            sprintf(fname,"!%s/%s_%02d:%02d:%02d.%09ld.fits", logdir, IDname, uttime->tm_hour, uttime->tm_min, uttime->tm_sec, thetime->tv_nsec);
            sprintf(fname_asciilog,"%s/%s_%02d:%02d:%02d.%09ld.txt", logdir, IDname, uttime->tm_hour, uttime->tm_min, uttime->tm_sec, thetime->tv_nsec);
        }


        if(wOK==1) // normal step: a frame has arrived
        {
            /// measure time
            t = time(NULL);
            uttime = gmtime(&t);
            clock_gettime(CLOCK_REALTIME, thetime);

            if(index==0)
                fp = fopen(fname_asciilog, "w");

            ptr0 = (char*) data.image[ID].array.F;
            if(is3Dcube==1)
				ptr0 += framesize*data.image[ID].md[0].cnt1;
            
            ptr1 = (char*) data.image[IDb].array.F;
            ptr1 += framesize*index;

            memcpy((void *) ptr1, (void *) ptr0, framesize);

            fprintf(fp, "%02d:%02d:%02d.%09ld ", uttime->tm_hour, uttime->tm_min, uttime->tm_sec, thetime->tv_nsec);

            for(kw=0; kw<data.image[ID].md[0].NBkw; kw++)
            {
                switch (data.image[ID].kw[kw].type) {
                case 'L' :
                    fprintf(fp, " %ld", data.image[ID].kw[kw].value.numl);
                    break;
                case 'D' :
                    fprintf(fp, " %f", data.image[ID].kw[kw].value.numf);
                    break;
                }
            }
            fprintf(fp, "\n");

            index++;
        }


        /// cases:
        /// index>zsize-1  buffer full
        /// wOK==0 && index>0
        if(  (index>zsize-1)  ||  ((wOK==0)&&(index>0)) )
        {
            /// save image
            sprintf(iname, "logbuff%d", buffer);
            //          printf("Saving %s -> %s\n", iname, fname);
            //         fflush(stdout);
            if(wOK==1) // image has arrived
            {
                strcpy(tmsg->iname, iname);
                strcpy(tmsg->fname, fname);
                tmsg->partial = 0; // full cube
            }

            fclose(fp);

            if(tOK == 1)
            {
                //           printf("WAITING FOR SAVE THREAD TO COMPLETE ...");
                //          fflush(stdout);
                pthread_join(thread_savefits, (void**)&thread_savefits);
                //          printf("OK\n");
                //          fflush(stdout);
            }
			
			strcpy(tmsg->iname, iname);
            iret_savefits = pthread_create( &thread_savefits, NULL, save_fits_function, tmsg);
            
            logshimconf[0].cnt ++;
            
            tOK = 1;
            if(iret_savefits)
            {
                fprintf(stderr,"Error - pthread_create() return code: %d\n", iret_savefits);
                exit(EXIT_FAILURE);
            }

            index = 0;
            buffer++;
            if(buffer==2)
                buffer = 0;
            //            printf("[%ld -> %d]", cnt, buffer);
            //           fflush(stdout);
            if(buffer==0)
                IDb = IDb0;
            else
                IDb = IDb1;

			logshimconf[0].filecnt ++;            
        }


        cnt = data.image[ID].md[0].cnt0;
    }

    free(imsizearray);

    return(0);
}









