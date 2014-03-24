#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <ncurses.h>

#include "../Cfits.h"
#include "../00CORE/00CORE.h"
#include "../COREMOD_memory/COREMOD_memory.h"
#include "../COREMOD_iofits/COREMOD_iofits.h"

 
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

/*void print_sys_mem_info()
{
   for solaris 
  printf("total computer RAM  %ld x %ld = %ld Mb\n",sysconf(_SC_PHYS_PAGES),sysconf(_SC_PAGESIZE), sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGESIZE)/1024/1024);
  printf("free RAM  %ld x %ld = %ld Mb\n",sysconf(_SC_AVPHYS_PAGES),sysconf(_SC_PAGESIZE), sysconf(_SC_AVPHYS_PAGES)*sysconf(_SC_PAGESIZE)/1024/1024);
}
*/



/*
* Forward references 
*/

int list_image_ID();
int list_variable_ID();
long create_image_ID(char *name, long naxis, long *size, int atype);
int memory_monitor(char *termttyname);
int list_image_ID_ncurses();
void close_list_image_ID_ncurses( void );
long copy_image_ID(char *name, char *newname);


// create image, default precision
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
	create_image_ID(data.cmdargtoken[1].val.string, naxis, imsize, FLOAT);
	break;
      case 1:
	create_image_ID(data.cmdargtoken[1].val.string, naxis, imsize, DOUBLE);
	break;
      }
      free(imsize);
    }
  else
    return 1;
}



int create_2Dimage_float()
{
  long *imsize;

  // CHECK ARGS
  printf("CREATING IMAGE\n");
  imsize = (long*) malloc(sizeof(long)*2);

  imsize[0] = data.cmdargtoken[2].val.numl;
  imsize[1] = data.cmdargtoken[3].val.numl;

  create_image_ID(data.cmdargtoken[1].val.string, 2, imsize, FLOAT);

  free(imsize);

  return 0;
}

int create_3Dimage_float()
{
  long *imsize;

  // CHECK ARGS
  printf("CREATING 3D IMAGE\n");
  imsize = (long*) malloc(sizeof(long)*3);

  imsize[0] = data.cmdargtoken[2].val.numl;
  imsize[1] = data.cmdargtoken[3].val.numl;
  imsize[2] = data.cmdargtoken[4].val.numl;

  create_image_ID(data.cmdargtoken[1].val.string, 3, imsize, FLOAT);

  free(imsize);

  return 0;
}


int memory_monitor_cli()
{
  memory_monitor(data.cmdargtoken[1].val.string);
  return 0;
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
 
  strcpy(data.cmd[data.NBcmd].key,"creaim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = create_image_cli;
  strcpy(data.cmd[data.NBcmd].info,"create image, default precision");
  strcpy(data.cmd[data.NBcmd].syntax,"<name> <xsize> <ysize> <opt: zsize>");
  strcpy(data.cmd[data.NBcmd].example,"creaim imname 512 512");
  strcpy(data.cmd[data.NBcmd].Ccall,"long create_image_ID(char *name, long naxis, long *size, int atype)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"crea3dim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = create_3Dimage_float;
  strcpy(data.cmd[data.NBcmd].info,"creates 3D image, single precision");
  strcpy(data.cmd[data.NBcmd].syntax,"<name> <xsize> <ysize> <zsize>");
  strcpy(data.cmd[data.NBcmd].example,"crea3dim imname 512 512 100");
  strcpy(data.cmd[data.NBcmd].Ccall,"long create_image_ID(char *name, long naxis, long *size, FLOAT)");
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
  strcpy(data.cmd[data.NBcmd].Ccall,"int mk_complex_from_reim(char *re_name, char *im_name, char *out_name)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"c2ap");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = mk_amph_from_complex_cli;
  strcpy(data.cmd[data.NBcmd].info,"complex -> ampl, pha");
  strcpy(data.cmd[data.NBcmd].syntax,"complex ampl pha");
  strcpy(data.cmd[data.NBcmd].example,"c2ap imc ima imp");
  strcpy(data.cmd[data.NBcmd].Ccall,"int mk_complex_from_amph(char *re_name, char *im_name, char *out_name)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"rmall");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = clearall;
  strcpy(data.cmd[data.NBcmd].info,"remove all images");
  strcpy(data.cmd[data.NBcmd].syntax,"no argument");
  strcpy(data.cmd[data.NBcmd].example,"rmall");
  strcpy(data.cmd[data.NBcmd].Ccall,"int clearall()");
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

  for(i=0;i<data.NB_MAX_IMAGE;i++)
    {
      if(data.image[i].used==1)
	total += data.image[i].nelement*TYPESIZE[data.image[i].atype];
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
	  if((strncmp(name,data.image[i].name,strlen(name))==0)&&(data.image[i].name[strlen(name)]=='\0'))
	    {
	      found = 1;
	      tmp = i;
	      clock_gettime(CLOCK_REALTIME, &timenow);
	      data.image[i].last_access = 1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec;	      
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
	  if((strncmp(name,data.image[i].name,strlen(name))==0)&&(data.image[i].name[strlen(name)]=='\0'))
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

  ID = image_ID(imname);
 
  if (ID!=-1)
    {
      data.image[ID].used = 0;
      
      if(data.image[ID].atype==CHAR)
	{
	  if(data.image[ID].array.C == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
	      exit(0);
	    }
	  free(data.image[ID].array.C);
	  data.image[ID].array.C = NULL;
	}	 
      if(data.image[ID].atype==INT)
	{
	  if(data.image[ID].array.I == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
	      exit(0);
	    }
	  free(data.image[ID].array.I);
	  data.image[ID].array.I = NULL;
	}
      if(data.image[ID].atype==FLOAT)
	{
	  if(data.image[ID].array.F == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
	      exit(0);
	    }
	  free(data.image[ID].array.F);
	  data.image[ID].array.F = NULL;
	}
      if(data.image[ID].atype==DOUBLE)
	{
	  if(data.image[ID].array.D == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
	      exit(0);
	    }
	  free(data.image[ID].array.D);
	  data.image[ID].array.D = NULL;
	}
      if(data.image[ID].atype==COMPLEX_FLOAT)
	{
	  if(data.image[ID].array.CF == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
	      exit(0);
	    }
	  free(data.image[ID].array.CF);
	  data.image[ID].array.CF = NULL;
	}
       if(data.image[ID].atype==COMPLEX_DOUBLE)
	{
	  if(data.image[ID].array.CD == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"data array pointer is null\n");
	      exit(0);
	    }
	  free(data.image[ID].array.CD);
	  data.image[ID].array.CD = NULL;
	}

       /*      free(data.image[ID].size);*/
       //      data.image[ID].last_access = 0;
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

  for (i=0;i<data.NB_MAX_IMAGE;i++)
    {
      
      if((data.image[i].used==1)&&((strncmp(prefix,data.image[i].name,strlen(prefix)))==0))
	{
	  printf("deleting image %s\n",data.image[i].name);
	  delete_image_ID(data.image[i].name);
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
/* all Cfits images should be created by this function */
long create_image_ID(char *name, long naxis, long *size, int atype)
{
  long ID;
  long i,ii;
  time_t lt;
  long nelement;
  struct timespec timenow;

  ID = -1;
  if(image_ID(name)==-1)
    {
      ID = next_avail_image_ID();
      data.image[ID].used = 1;
      data.image[ID].atype = atype;
      data.image[ID].naxis = naxis;
      strcpy(data.image[ID].name,name);
      

      for(i=0;i<naxis;i++)
	data.image[ID].size[i] = size[i];
	
      nelement = 1;
      for(i=0;i<naxis;i++)
	nelement*=size[i];


      if(atype==CHAR)
	{
	  data.image[ID].array.C = (char*) calloc ((size_t) nelement,sizeof(char));
	  if(data.image[ID].array.C == NULL)
	    {    
	      printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
	      fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
	      fprintf(stderr,"Image name = %s\n",name);
	      fprintf(stderr,"Image size = ");
	      fprintf(stderr,"%ld",size[0]);
	      for(i=1;i<naxis;i++)
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
	  data.image[ID].array.I = (int*) calloc ((size_t) nelement,sizeof(int));
	  if(data.image[ID].array.I == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
	      fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
	      fprintf(stderr,"Image name = %s\n",name);
	      fprintf(stderr,"Image size = ");
	      fprintf(stderr,"%ld",size[0]);
	      for(i=1;i<naxis;i++)
		fprintf(stderr,"x%ld",size[i]);
	      fprintf(stderr,"\n");
	      fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(int));
	      fprintf(stderr," %c[%d;m",(char) 27, 0);
	      list_image_ID();
	      exit(0);   
	    }
	}
      if(atype==FLOAT)
	{
	  data.image[ID].array.F = (float*) calloc ((size_t) nelement, sizeof(float));
	  if(data.image[ID].array.F == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
	      fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
	      fprintf(stderr,"Image name = %s\n",name);
	      fprintf(stderr,"Image size = ");
	      fprintf(stderr,"%ld",size[0]);
	      for(i=1;i<naxis;i++)
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
	  data.image[ID].array.D = (double*) calloc ((size_t) nelement,sizeof(double));
	  if(data.image[ID].array.D == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
	      fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
	      fprintf(stderr,"Image name = %s\n",name);
	      fprintf(stderr,"Image size = ");
	      fprintf(stderr,"%ld",size[0]);
	      for(i=1;i<naxis;i++)
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
	  
	  data.image[ID].array.CF = (complex_float*) calloc ((size_t) nelement,sizeof(complex_float));
	  for(ii=0;ii<nelement;ii++)
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
	      for(i=1;i<naxis;i++)
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
	  data.image[ID].array.CD = (complex_double*) calloc ((size_t) nelement,sizeof(complex_double));
	  if(data.image[ID].array.CD == NULL)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
	      fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
	      fprintf(stderr,"Image name = %s\n",name);
	      fprintf(stderr,"Image size = ");
	      fprintf(stderr,"%ld",size[0]);
	      for(i=1;i<naxis;i++)
		fprintf(stderr,"x%ld",size[i]);
	      fprintf(stderr,"\n");
	      fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n",nelement,1.0/1024/1024*nelement*sizeof(double)*2);
	      fprintf(stderr," %c[%d;m",(char) 27, 0);
	      list_image_ID();
	      exit(0);   
	    }
	}
    
      clock_gettime(CLOCK_REALTIME, &timenow);
      data.image[ID].last_access = 1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec;
      data.image[ID].creation_time = data.image[ID].last_access;

      data.image[ID].nelement = nelement;
    }
  else
    {
      //      printf("Cannot create image : name \"%s\" already in use\n",name);
      ID = image_ID(name);

      if(data.image[ID].atype != atype)
	{
	  fprintf(stderr,"%c[%d;%dm ERROR: [ %s %s %d ] %c[%d;m\n", (char) 27, 1, 31, __FILE__, __func__, __LINE__, (char) 27, 0);
	  fprintf(stderr,"%c[%d;%dm Pre-existing image \"%s\" has wrong type %c[%d;m\n", (char) 27, 1, 31,name, (char) 27, 0);
	  exit(0);
	}
      if(data.image[ID].naxis != naxis)
	{
	  fprintf(stderr,"%c[%d;%dm ERROR: [ %s %s %d ] %c[%d;m\n", (char) 27, 1, 31, __FILE__, __func__, __LINE__, (char) 27, 0);
	  fprintf(stderr,"%c[%d;%dm Pre-existing image \"%s\" has wrong naxis %c[%d;m\n", (char) 27, 1, 31,name, (char) 27, 0);
	  exit(0);
	}
      
      for(i=0;i<naxis;i++)
	if(data.image[ID].size[i] != size[i])
	  {
	    fprintf(stderr,"%c[%d;%dm ERROR: [ %s %s %d ] %c[%d;m\n", (char) 27, 1, 31, __FILE__, __func__, __LINE__, (char) 27, 0);
	    fprintf(stderr,"%c[%d;%dm Pre-existing image \"%s\" has wrong size %c[%d;m\n", (char) 27, 1, 31,name, (char) 27, 0);
	    exit(0);
	  }
    }

  if(MEM_MONITOR == 1)
    list_image_ID_ncurses();
  
  return(ID);
}

long create_1Dimage_ID(char *ID_name, long xsize)
{
  long ID = -1;
  long naxis = 1;
  long naxes[1];

  naxes[0]=xsize;

  if(data.precision == 0)
    ID = create_image_ID(ID_name,naxis,naxes,3); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,4); // double precision
 
  return(ID);
}

long create_1DCimage_ID(char *ID_name, long xsize)
{
  long ID = -1;
  long naxis=1;
  long naxes[1];

  naxes[0]=xsize;

  if(data.precision == 0)
    ID = create_image_ID(ID_name,naxis,naxes,5); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,6); // double precision

  return(ID);
}

long create_2Dimage_ID(char *ID_name, long xsize, long ysize)
{
  long ID = -1;
  long naxis=2;
  long naxes[2];

  naxes[0]=xsize;
  naxes[1]=ysize;

  if(data.precision == 0)
    ID = create_image_ID(ID_name,naxis,naxes,3); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,4); // double precision

  return(ID);
}

long create_2Dimagedouble_ID(char *ID_name, long xsize, long ysize)
{
  long ID = -1;
  long naxis=2;
  long naxes[2];

  naxes[0]=xsize;
  naxes[1]=ysize;

  ID = create_image_ID(ID_name,naxis,naxes,4);

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
    ID = create_image_ID(ID_name,naxis,naxes,5); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,6); // double precision

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

  ID = create_image_ID(ID_name,naxis,naxes,3); // single precision

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

  ID = create_image_ID(ID_name,naxis,naxes,4); // double precision
  
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
    ID = create_image_ID(ID_name,naxis,naxes,3); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,4); // double precision

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
    ID = create_image_ID(ID_name,naxis,naxes,5); // single precision
  if(data.precision == 1)
    ID = create_image_ID(ID_name,naxis,naxes,6); // double precision

  return(ID);
}


long copy_image_ID(char *name, char *newname)
{
  long ID,IDout;
  long naxis;
  long *size = NULL;
  int atype;
  long nelement;
  long i;
  
  ID = image_ID(name);
  naxis = data.image[ID].naxis;

  size = (long*) malloc(sizeof(long)*naxis);
  if(size==NULL)
    {
      printERROR(__FILE__,__func__,__LINE__,"malloc error");
      exit(0);
    }

  for(i=0;i<naxis;i++)
    size[i] = data.image[ID].size[i];
  atype  = data.image[ID].atype;

  nelement = data.image[ID].nelement;

  IDout = image_ID(newname);
  if(IDout==-1)
    {
      create_image_ID(newname,naxis,size,atype);
      IDout = image_ID(newname);
    }
  else
    {
      // verify newname has the right size and type
      if(data.image[ID].nelement!=data.image[IDout].nelement)
	{
	  fprintf(stderr,"ERROR [copy_image_ID]: images %s and %s do not have the same size\n",name,newname);
	  exit(0);
	}
      if(data.image[ID].atype!=data.image[IDout].atype)
	{
	  fprintf(stderr,"ERROR [copy_image_ID]: images %s and %s do not have the same type\n",name,newname);	
	  exit(0);
	}
    }

  if(atype==CHAR)
    memcpy (data.image[IDout].array.C,data.image[ID].array.C,sizeof(char)*nelement);
  
  if(atype==INT)
    memcpy (data.image[IDout].array.I,data.image[ID].array.I,sizeof(int)*nelement);
  
  if(atype==FLOAT)
    memcpy (data.image[IDout].array.F,data.image[ID].array.F,sizeof(float)*nelement);
  
  if(atype==DOUBLE)
    memcpy (data.image[IDout].array.D,data.image[ID].array.D,sizeof(double)*nelement);
  
  if(atype==COMPLEX_FLOAT)
    memcpy (data.image[IDout].array.CF,data.image[ID].array.CF,sizeof(float)*2*nelement);
  
  if(atype==COMPLEX_DOUBLE)
    memcpy (data.image[IDout].array.CD,data.image[ID].array.CD,sizeof(double)*2*nelement);
  
  free(size);

  return(IDout);
}


/* creates an ID */
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
      strcpy(data.variable[ID].name,name);
      data.variable[ID].value = value;

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
	  atype = data.image[i].atype;
	  tmp_long = ((long long) (data.image[i].nelement)) * TYPESIZE[atype];
	  
	  printw("%4ld ", i);
	  attron(A_BOLD | COLOR_PAIR(6));
	  sprintf(str, "%10s ", data.image[i].name);
	  printw(str);
	  attroff(A_BOLD | COLOR_PAIR(6));
	  
	  sprintf(str, "[ %6ld",data.image[i].size[0]);
	  
	  for(j=1;j<data.image[i].naxis;j++)
	    {
	      sprintf(str, "%s x %6ld", str, data.image[i].size[j]);
	    }
	  sprintf(str, "%s]", str);
	  
	  printw("%-28s", str);
	  
	  attron(COLOR_PAIR(3));
	  n = 0;
	  if(atype==CHAR)
	    n = snprintf(type,STYPESIZE,"CHAR");	
	  if(atype==INT)
	    n = snprintf(type,STYPESIZE,"INT");
	  if(atype==FLOAT)
	    n = snprintf(type,STYPESIZE,"FLOAT");
	  if(atype==DOUBLE)
	    n = snprintf(type,STYPESIZE,"DOUBLE");
	  if(atype==COMPLEX_FLOAT)
	    n = snprintf(type,STYPESIZE,"CFLOAT");
	  if(atype==COMPLEX_DOUBLE)
	    n = snprintf(type,STYPESIZE,"CDOUBLE");
	  printw("%7s ", type);
	  
	  attroff(COLOR_PAIR(3));
	  
	  if(n >= STYPESIZE)
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  
	  printw("%10ld Kb %6.2f   ", (long) (tmp_long/1024), (float) (100.0*tmp_long/sizeb));
	  
	  timediff = (1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec) - data.image[i].last_access;
	  
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

  sizeb = compute_image_memory();

  clock_gettime(CLOCK_REALTIME, &timenow);

  fprintf(fo, "\n");
  fprintf(fo, "INDEX    NAME         SIZE                    TYPE        SIZE  [percent]    LAST ACCESS\n");
  fprintf(fo, "\n");

  for (i=0;i<data.NB_MAX_IMAGE;i++)
    if(data.image[i].used==1) 
      {
	  atype = data.image[i].atype;
	  tmp_long = ((long long) (data.image[i].nelement)) * TYPESIZE[atype];
	  
	  fprintf(fo, "%4ld %c[%d;%dm%10s%c[%d;m ",i, (char) 27, 1, 33, data.image[i].name, (char) 27, 0);
	  //fprintf(fo, "%s", str);
	  
	  sprintf(str, "[ %6ld",data.image[i].size[0]);
	  
	  for(j=1;j<data.image[i].naxis;j++)
	    {
	      sprintf(str, "%s x %6ld", str, data.image[i].size[j]);
	    }
	  sprintf(str, "%s]", str);
	  
	  fprintf(fo, "%-28s", str);
	  

	  n = 0;
	  if(atype==CHAR)
	    n = snprintf(type,STYPESIZE,"CHAR");	
	  if(atype==INT)
	    n = snprintf(type,STYPESIZE,"INT");
	  if(atype==FLOAT)
	    n = snprintf(type,STYPESIZE,"FLOAT");
	  if(atype==DOUBLE)
	    n = snprintf(type,STYPESIZE,"DOUBLE");
	  if(atype==COMPLEX_FLOAT)
	    n = snprintf(type,STYPESIZE,"CFLOAT");
	  if(atype==COMPLEX_DOUBLE)
	    n = snprintf(type,STYPESIZE,"CDOUBLE");
	  fprintf(fo, "%7s ", type);
	  
	  
	  if(n >= STYPESIZE)
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  
	  fprintf(fo, "%10ld Kb %6.2f   ", (long) (tmp_long/1024), (float) (100.0*tmp_long/sizeb));
	  
	  timediff = (1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec) - data.image[i].last_access;
	  
	  fprintf(fo, "%15.9f\n", timediff);
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
	atype = data.image[i].atype;
	fprintf(fp,"%ld %s",i, data.image[i].name);
	fprintf(fp," %ld",data.image[i].naxis);
	for(j=0;j<data.image[i].naxis;j++)
	  fprintf(fp," %ld",data.image[i].size[j]);
	
	n = 0;
	if(atype==CHAR)
	  n = snprintf(type,STYPESIZE,"CHAR");
	if(atype==INT)
	  n = snprintf(type,STYPESIZE,"INT");
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
      printf("%4ld %10s %25.18g\n",i, data.variable[i].name,data.variable[i].value);
  
  return(0);
}

long chname_image_ID(char *ID_name, char *new_name)
{
  long ID;
  
  ID=-1;
  if((image_ID(new_name)==-1)&&(variable_ID(new_name)==-1))
    {  
      ID = image_ID(ID_name);
      strcpy(data.image[ID].name,new_name);
      //      if ( Debug > 0 ) { printf("change image name %s -> %s\n",ID_name,new_name);}
    }
  else
    printf("Cannot change name : new name already in use\n");

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
  
  atype_re = data.image[IDre].atype;
  atype_im = data.image[IDim].atype;
  naxis = data.image[IDre].naxis;

  naxes = (long *) malloc(sizeof(long)*naxis);
  if(naxes==NULL)
    {
      printERROR(__FILE__,__func__,__LINE__,"malloc error");
      exit(0);
    }

  for(i=0;i<naxis;i++)
    naxes[i] = data.image[IDre].size[i];
  nelement = data.image[IDre].nelement;


  if((atype_re==FLOAT)&&(atype_im==FLOAT))
    {
      atype_out = COMPLEX_FLOAT;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out);
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CF[ii].re = data.image[IDre].array.F[ii];
	  data.image[IDout].array.CF[ii].im = data.image[IDim].array.F[ii];
	}  
    }
  else if((atype_re==FLOAT)&&(atype_im==DOUBLE))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out);
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CD[ii].re = data.image[IDre].array.F[ii];
	  data.image[IDout].array.CD[ii].im = data.image[IDim].array.D[ii];
	}  
    }
  else if((atype_re==DOUBLE)&&(atype_im==FLOAT))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out);
      for(ii=0;ii<nelement;ii++)
	{
	  data.image[IDout].array.CD[ii].re = data.image[IDre].array.D[ii];
	  data.image[IDout].array.CD[ii].im = data.image[IDim].array.F[ii];
	}   
    }
  else if((atype_re==DOUBLE)&&(atype_im==DOUBLE))
    {
      atype_out = COMPLEX_DOUBLE;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out);
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
  atype_am = data.image[IDam].atype;
  atype_ph = data.image[IDph].atype;

  naxis = data.image[IDam].naxis;
  for(i=0;i<naxis;i++)
    naxes[i] = data.image[IDam].size[i];
  nelement = data.image[IDam].nelement;

  if((atype_am==FLOAT)&&(atype_ph==FLOAT))
    {
      atype_out = COMPLEX_FLOAT;
      IDout = create_image_ID(out_name,naxis,naxes,atype_out);
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
      IDout = create_image_ID(out_name,naxis,naxes,atype_out);
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
      IDout = create_image_ID(out_name,naxis,naxes,atype_out);
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
      IDout = create_image_ID(out_name,naxis,naxes,atype_out);
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
  atype = data.image[IDin].atype;
  naxis = data.image[IDin].naxis;
  for(i=0;i<naxis;i++)
    naxes[i] = data.image[IDin].size[i];
  nelement = data.image[IDin].nelement;

  if(atype == COMPLEX_FLOAT) // single precision
    {
      IDre = create_image_ID(re_name,naxis,naxes,FLOAT);
      IDim = create_image_ID(im_name,naxis,naxes,FLOAT);
  
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
      IDre = create_image_ID(re_name,naxis,naxes,DOUBLE);
      IDim = create_image_ID(im_name,naxis,naxes,DOUBLE);
  
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
  atype = data.image[IDin].atype;
  naxis = data.image[IDin].naxis;

  for(i=0;i<naxis;i++)
    naxes[i] = data.image[IDin].size[i];
  nelement = data.image[IDin].nelement;

  if(atype==COMPLEX_FLOAT) // single precision
    {
      IDam = create_image_ID(am_name,naxis,naxes,FLOAT);
      IDph = create_image_ID(ph_name,naxis,naxes,FLOAT);
      
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
      IDam = create_image_ID(am_name,naxis,naxes,DOUBLE);
      IDph = create_image_ID(ph_name,naxis,naxes,DOUBLE);
      
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
	delete_image_ID(data.image[ID].name);
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
  if(data.image[ID].naxis!=2)
    value=0;
  if(value==1)
    {
      if(data.image[ID].size[0]!=xsize)
	value = 0;
      if(data.image[ID].size[1]!=ysize)
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
  if(data.image[ID].naxis!=3)
    {
      /*      printf("Wrong naxis : %ld - should be 3\n",data.image[ID].naxis);*/
      value = 0;
    }
  if(value==1)
    {
      if(data.image[ID].size[0]!=xsize)
	{
	  /*	  printf("Wrong xsize : %ld - should be %ld\n",data.image[ID].size[0],xsize);*/
	  value = 0;
	}
      if(data.image[ID].size[1]!=ysize)
	{
	  /*	  printf("Wrong ysize : %ld - should be %ld\n",data.image[ID].size[1],ysize);*/
	  value = 0;
	}
      if(data.image[ID].size[2]!=zsize)
	{
	  /*	  printf("Wrong zsize : %ld - should be %ld\n",data.image[ID].size[2],zsize);*/
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
  atype = data.image[ID].atype;

  if(data.image[ID].naxis!=3)
    {
      n = snprintf(errmsg,SBUFFERSIZE,"Wrong naxis : %ld - should be 3\n",data.image[ID].naxis);
      if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

      printERROR(__FILE__,__func__,__LINE__,errmsg);
      exit(0);
    }
  xsize = data.image[ID].size[0];
  ysize = data.image[ID].size[1];
  zsize = data.image[ID].size[2];

  if(atype==FLOAT) // single precision
    {
      if(orientation==0)
	{
	  xsize1 = zsize;
	  ysize1 = ysize;
	  zsize1 = xsize;
	  IDout = create_3Dimage_ID_float(ID_out_name,xsize1,ysize1,zsize1);
	  for(ii=0;ii<xsize1;ii++)
	    for(jj=0;jj<ysize1;jj++)
	      for(kk=0;kk<zsize1;kk++)
		data.image[IDout].array.F[kk*ysize1*xsize1+jj*xsize1+ii] = data.image[ID].array.F[ii*xsize*ysize+jj*xsize+kk];
	}
      else
	{
	  xsize1 = xsize;
	  ysize1 = zsize;
	  zsize1 = ysize;
	  IDout = create_3Dimage_ID_float(ID_out_name,xsize1,ysize1,zsize1);
	  for(ii=0;ii<xsize1;ii++)
	    for(jj=0;jj<ysize1;jj++)
	      for(kk=0;kk<zsize1;kk++)
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
	  for(ii=0;ii<xsize1;ii++)
	    for(jj=0;jj<ysize1;jj++)
	      for(kk=0;kk<zsize1;kk++)
		data.image[IDout].array.D[kk*ysize1*xsize1+jj*xsize1+ii] = data.image[ID].array.D[ii*xsize*ysize+jj*xsize+kk];
	}
      else
	{
	  xsize1 = xsize;
	  ysize1 = zsize;
	  zsize1 = ysize;
	  IDout = create_3Dimage_ID_double(ID_out_name,xsize1,ysize1,zsize1);
	  for(ii=0;ii<xsize1;ii++)
	    for(jj=0;jj<ysize1;jj++)
	      for(kk=0;kk<zsize1;kk++)
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
