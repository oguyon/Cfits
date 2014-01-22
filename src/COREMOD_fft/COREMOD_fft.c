#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
#endif

#include <Cfits.h>
#include <00CORE/00CORE.h>
#include <COREMOD_fft/COREMOD_fft.h>
#include <COREMOD_memory/COREMOD_memory.h>
#include <COREMOD_iofits/COREMOD_iofits.h>
#include <COREMOD_arith/COREMOD_arith.h>
#include <COREMOD_tools/COREMOD_tools.h>


#define PI 3.14159265358979323846264338328

#define SWAP(x,y)  tmp=(x);x=(y);y=tmp;
#define CSWAP(x,y)  tmp=(x.re);x.re=(y.re);y.re=tmp;tmp=(x.im);x.im=(y.im);y.im=tmp;

#define SBUFFERSIZE 1000

#define FFTWOPTMODE FFTW_ESTIMATE
//#define FFTWOPTMODE FFTW_MEASURE
//#define FFTWOPTMODE FFTW_PATIENT
//#define FFTWOPTMODE FFTW_EXHAUSTIVE


int NB_FFTW_THREADS = 1;

extern DATA data;




int fft_setNthreads(int nt)
{
  printf("set number of thread to %d (FFTWMT)\n",nt);
# ifdef FFTWMT
  fftwf_cleanup_threads();
  fftwf_cleanup();

  //  printf("Multi-threaded fft enabled, max threads = %d\n",nt);
  fftwf_init_threads();
  fftwf_plan_with_nthreads(nt);
# endif
  

  import_wisdom();
  
  return(0);
}

int import_wisdom()
{
  FILE *fp;
  char wisdom_file[SBUFFERSIZE];  
  char warnmessg[SBUFFERSIZE];
  int n;

  # ifdef FFTWMT
  if(sizeof(PRECISION)==sizeof(float))
    n = snprintf(wisdom_file,SBUFFERSIZE,"%s/fftwf_mt_wisdom.dat", CONFIGDIR);
  else
    n = snprintf(wisdom_file,SBUFFERSIZE,"%s/fftw_mt_wisdom.dat", CONFIGDIR);
  # endif
    
  # ifndef FFTWMT
  if(sizeof(PRECISION)==sizeof(float))
    n = snprintf(wisdom_file,SBUFFERSIZE,"%s/fftwf_wisdom.dat", CONFIGDIR);
  else
    n = snprintf(wisdom_file,SBUFFERSIZE,"%s/fftw_wisdom.dat", CONFIGDIR);
  # endif
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  
  if((fp = fopen(wisdom_file,"r"))==NULL)
    {
      n = snprintf(warnmessg,SBUFFERSIZE,"No wisdom file in %s\n FFTs will not be optimized, and may run slower than if a wisdom file is used\n type \"initfft\" to create the wisdom file (this will take time)", wisdom_file);
      if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
      printWARNING(__FILE__,__func__,__LINE__,warnmessg);      
    } else  { 
    if(sizeof(PRECISION)==sizeof(float))
      {
	if (fftwf_import_wisdom_from_file(fp)==0)
	  printERROR(__FILE__,__func__,__LINE__,"Error reading wisdom");
      }
    else
      {
	if (fftw_import_wisdom_from_file(fp)==0)
	  printERROR(__FILE__,__func__,__LINE__,"Error reading wisdom");
      }
    fclose(fp); 
  }
  
  return(0);
}

int export_wisdom() 
{
  FILE *fp;
  char wisdom_file[SBUFFERSIZE];
  char errmessg[SBUFFERSIZE];
  int n;

  # ifdef FFTWMT
  if(sizeof(PRECISION)==sizeof(float))
    n = snprintf(wisdom_file,SBUFFERSIZE,"%s/fftwf_mt_wisdom.dat", CONFIGDIR);
  else
    n = snprintf(wisdom_file,SBUFFERSIZE,"%s/fftw_mt_wisdom.dat", CONFIGDIR);
  # endif
    
  # ifndef FFTWMT
  if(sizeof(PRECISION)==sizeof(float))
    n = snprintf(wisdom_file,SBUFFERSIZE,"%s/fftwf_wisdom.dat", CONFIGDIR);
  else
    n = snprintf(wisdom_file,SBUFFERSIZE,"%s/fftw_wisdom.dat", CONFIGDIR);
  # endif
  if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  if((fp = fopen(wisdom_file, "w"))==NULL)
    {
      n = snprintf(errmessg,SBUFFERSIZE,"Error creating wisdom file \"%s\"",wisdom_file);
      if(n >= SBUFFERSIZE) 
	printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
      printERROR(__FILE__,__func__,__LINE__,errmessg);
      exit(0);
    }

  if(sizeof(PRECISION)==sizeof(float))
    fftwf_export_wisdom_to_file(fp);
  else
    fftw_export_wisdom_to_file(fp);
  
  fclose(fp);

  return(0);
}


/*^-----------------------------------------------------------------------------
| 
| 
| 
|
| 
|
|
+-----------------------------------------------------------------------------*/
int init_fftw_plans(int mode)
{
  int n;
  int size;
  fftwf_complex *inf = NULL;
  fftwf_complex *outf = NULL;
  float *rinf = NULL;
  unsigned int plan_mode;


  printf("Optimization of FFTW\n");
  printf("The optimization is done for 2D complex to complex FFTs, with size equal to 2^n x 2^n\n");
  printf("You can kill the optimization anytime, and resume later where it previously stopped.\nAfter each size is optimized, the result is saved\n");
  printf("It might be a good idea to run this overnight or when your computer is not busy\n");

  fflush(stdout);

  size = 1;

  //  plan_mode = FFTWOPTMODE;
  plan_mode = FFTW_EXHAUSTIVE;
  
  for(n=0;n<14;n++)
    {
      if(mode==0)
	{
	  printf("Optimizing 2D FFTs - size = %d\n",size);
	  fflush(stdout);
	}
      rinf = (float*) fftwf_malloc(size*size*sizeof(float));
      inf = (fftwf_complex*) fftwf_malloc(size*size*sizeof(fftwf_complex));
      outf = (fftwf_complex*) fftwf_malloc(size*size*sizeof(fftwf_complex));

      fftwf_plan_dft_2d(size, size, inf, outf, FFTW_FORWARD, plan_mode);
      fftwf_plan_dft_2d(size, size, inf, outf, FFTW_BACKWARD, plan_mode);
      fftwf_plan_dft_r2c_2d(size, size, rinf, outf, plan_mode);
      
      fftwf_free(inf);
      fftwf_free(rinf);
      fftwf_free(outf);
      
      size*=2;
      if(mode==0)
	export_wisdom();
    }
  size = 1;
  for(n=0;n<15;n++)
    {
      if(mode==0)
	{
	  printf("Optimizing 1D FFTs - size = %d\n",size);
	  fflush(stdout);
	}
      rinf = (float*) fftwf_malloc(size*sizeof(float));
      inf = (fftwf_complex*) fftwf_malloc(size*sizeof(fftwf_complex));
      outf = (fftwf_complex*) fftwf_malloc(size*sizeof(fftwf_complex));

      fftwf_plan_dft_1d(size, inf, outf, FFTW_FORWARD, plan_mode);
      fftwf_plan_dft_1d(size, inf, outf, FFTW_BACKWARD, plan_mode);
      fftwf_plan_dft_r2c_1d(size, rinf, outf, plan_mode);
      
      fftwf_free(inf);
      fftwf_free(rinf);
      fftwf_free(outf);
      
      size*=2;
      if(mode==0)
	export_wisdom();
    }
  

  export_wisdom();
  
  return(0);
}

int init_fftw_plans0()
{
  init_fftw_plans(0);

  return(0);
}

int permut(char *ID_name)
{
  PRECISION tmp;
  long *naxes;   
  int ID;
  long xhalf,yhalf;
  long ii,jj,kk;
  long i;
  long naxis;
  int atype;
  int OK=0;

  ID = image_ID(ID_name);
  naxis = data.image[ID].naxis;
  naxes = (long*) malloc(naxis*sizeof(long));
  for(i=0;i<naxis;i++)
    naxes[i] = data.image[ID].size[i];
  atype = data.image[ID].atype;

  tmp=0;

  if(atype==Dtype)
    {
      if(naxis==1)
	{
	  OK=1;
	  xhalf = (long) (naxes[0]/2);
	  for (ii = 0; ii < xhalf; ii++)
	    SWAP(data.image[ID].arrayD[ii],data.image[ID].arrayD[ii+xhalf])
	}
      if(naxis==2)
	{
	  OK=1;
	  xhalf = (long) (naxes[0]/2);
	  yhalf = (long) (naxes[1]/2);
	  for (jj = 0; jj < yhalf; jj++) 
	    for (ii = 0; ii < xhalf; ii++){
	      SWAP(data.image[ID].arrayD[jj*naxes[0]+ii],data.image[ID].arrayD[(jj+yhalf)*naxes[0]+(ii+xhalf)])
		}
	  for (jj = yhalf; jj < naxes[1]; jj++) 
	    for (ii = 0; ii < xhalf; ii++){
	      SWAP(data.image[ID].arrayD[jj*naxes[0]+ii],data.image[ID].arrayD[(jj-yhalf)*naxes[0]+(ii+xhalf)])
		}
	}
      if(naxis==3)
	{
	  OK=1;
	  xhalf = (long) (naxes[0]/2);
	  yhalf = (long) (naxes[1]/2);
	  for (jj = 0; jj < yhalf; jj++) 
	    for (ii = 0; ii < xhalf; ii++)
	      {
		for(kk=0;kk<naxes[2];kk++)
		  SWAP(data.image[ID].arrayD[kk*naxes[0]*naxes[1]+jj*naxes[0]+ii],data.image[ID].arrayD[kk*naxes[0]*naxes[1]+(jj+yhalf)*naxes[0]+(ii+xhalf)])
	      }
	  for (jj = yhalf; jj < naxes[1]; jj++) 
	    for (ii = 0; ii < xhalf; ii++)
	      {
		for(kk=0;kk<naxes[2];kk++)
		  SWAP(data.image[ID].arrayD[kk*naxes[0]*naxes[1]+jj*naxes[0]+ii],data.image[ID].arrayD[kk*naxes[0]*naxes[1]+(jj-yhalf)*naxes[0]+(ii+xhalf)])
	      } 
	}
    }
  
  if(atype==CDtype)
    {
      if(naxis==1)
	{
	  OK=1;
	  xhalf = (long) (naxes[0]/2);
	  for (ii = 0; ii < xhalf; ii++)
	    CSWAP(data.image[ID].arrayCD[ii],data.image[ID].arrayCD[ii+xhalf])
	}
      if(naxis==2)
	{
	  OK=1;
	  xhalf = (long) (naxes[0]/2);
	  yhalf = (long) (naxes[1]/2);
	  for (jj = 0; jj < yhalf; jj++) 
	    for (ii = 0; ii < xhalf; ii++){
	      CSWAP(data.image[ID].arrayCD[jj*naxes[0]+ii],data.image[ID].arrayCD[(jj+yhalf)*naxes[0]+(ii+xhalf)])
		}
	  for (jj = yhalf; jj < naxes[1]; jj++) 
	    for (ii = 0; ii < xhalf; ii++){
	      CSWAP(data.image[ID].arrayCD[jj*naxes[0]+ii],data.image[ID].arrayCD[(jj-yhalf)*naxes[0]+(ii+xhalf)])
		}
	}
      if(naxis==3)
	{
	  OK=1;
	  xhalf = (long) (naxes[0]/2);
	  yhalf = (long) (naxes[1]/2);
	  for (jj = 0; jj < yhalf; jj++) 
	    for (ii = 0; ii < xhalf; ii++)
	      {
		for(kk=0;kk<naxes[2];kk++)
		  CSWAP(data.image[ID].arrayCD[kk*naxes[0]*naxes[1]+jj*naxes[0]+ii],data.image[ID].arrayCD[kk*naxes[0]*naxes[1]+(jj+yhalf)*naxes[0]+(ii+xhalf)])
	      }
	  for (jj = yhalf; jj < naxes[1]; jj++) 
	    for (ii = 0; ii < xhalf; ii++)
	      {
		for(kk=0;kk<naxes[2];kk++)
		  CSWAP(data.image[ID].arrayCD[kk*naxes[0]*naxes[1]+jj*naxes[0]+ii],data.image[ID].arrayCD[kk*naxes[0]*naxes[1]+(jj-yhalf)*naxes[0]+(ii+xhalf)])
	      } 
	}
    }

  if(OK==0)
    printf("Error : data format not supported by permut\n");

  free(naxes);
  
  return(0);
}

int array_index(long size)
{
  int i;
  
  switch (size) {
  case 1: 
    i=0;
    break;
  case 2:
    i=1;
    break;
  case 4:
    i=2;
    break;
  case 8: 
    i=3;
    break;
  case 16:
    i=4;
    break;
  case 32:
    i=5;
    break;
  case 64: 
    i=6;
    break;
  case 128:
    i=7;
    break;
  case 256:
    i=8;
    break;
  case 512: 
    i=9;
    break;
  case 1024:
    i=10;
    break;
  case 2048:
    i=11;
    break;
  case 4096: 
    i=12;
    break;
  case 8192:
    i=13;
    break;
  case 16384:
    i=14;
    break;
  default:
    i=100;
  }

  return(i);
}





/* 1d complex fft */
int do1dfft(char *in_name, char *out_name)
{
  int *naxes;
  long *naxesl;
  long naxis;
  long IDin,IDout;
  long i;
  int OK=0;
  fftwf_plan plan;

  IDin=image_ID(in_name);
  naxis=data.image[IDin].naxis;
  naxes = (int *) malloc(naxis*sizeof(int));
  naxesl = (long *) malloc(naxis*sizeof(long));
  for (i=0;i<naxis;i++)
    {
      naxesl[i]= (long) data.image[IDin].size[i];
      naxes[i]= (int) data.image[IDin].size[i];
    }

  IDout = create_image_ID(out_name,naxis,naxesl,CDtype);

  if(naxis==1)
    {
      if(array_index(naxes[0])!=100)
	{
	  OK=1;
	  plan = fftwf_plan_dft_1d(naxes[0], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, -1,FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}
      else
	{
	  OK=1;
	  plan = fftwf_plan_dft_1d(naxes[0], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, -1,FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}
    }

  if(naxis==2)
    {
      if((naxes[1]==1)&&(array_index(naxes[0])!=100))
	{
	  OK=1;
	  plan = fftwf_plan_dft_1d(naxes[0], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, -1,FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}
      else
	{
	  OK=1;
	  plan = fftwf_plan_dft_1d(naxes[0], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, -1,FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}	
    }

  if(OK==0)
    {
      printf("Error : image dimension not appropriate for FFT\n");
    }
  free(naxes);
  free(naxesl);

  return(0);
}


/* 1d complex fft */
int do1drfft(char *in_name, char *out_name)
{
  int *naxes;
  long *naxesl;
  long *naxestmp;
  long naxis;
  long IDin,IDout,IDtmp;
  long i;
  int OK=0;
  long ii,jj;
  fftwf_plan plan;

  char ffttmpname[SBUFFERSIZE];
  int n;

  IDin=image_ID(in_name);
  naxis=data.image[IDin].naxis;
  naxes = (int *) malloc(naxis*sizeof(int));
  naxesl = (long *) malloc(naxis*sizeof(long));
  naxestmp = (long *) malloc(naxis*sizeof(long));

  for (i=0;i<naxis;i++)
    {
      naxesl[i]= (long) data.image[IDin].size[i];
      naxes[i]= (int) data.image[IDin].size[i];
      naxestmp[i]=data.image[IDin].size[i];
      if(i==0)
	naxestmp[i]=data.image[IDin].size[i]/2+1;
    }


  n = snprintf(ffttmpname,SBUFFERSIZE,"_ffttmpname_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
  IDtmp = create_image_ID(ffttmpname,naxis,naxestmp,CDtype);


  IDout = create_image_ID(out_name,naxis,naxesl,CDtype);

  if(naxis==2)
    {
      if((naxes[1]==1)&&(array_index(naxes[0])!=100))
	{
	  OK=1;
	  plan = fftwf_plan_dft_r2c_1d(naxes[0], data.image[IDin].arrayD, (fftwf_complex*) data.image[IDout].arrayCD, FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}
      else
	{
	  OK=1;
	  plan = fftwf_plan_dft_r2c_1d(naxes[0], data.image[IDin].arrayD, (fftwf_complex*) data.image[IDout].arrayCD, FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}	
      
      for(ii=0;ii<naxes[0]/2+1;ii++)
	for(jj=0;jj<naxes[1];jj++)
	  {
	    data.image[IDout].arrayCD[jj*naxes[0]+ii] = data.image[IDtmp].arrayCD[jj*naxestmp[0]+ii];
	  }
      for(ii=1;ii<naxes[0]/2+1;ii++)
	{
	  jj=0;
	  data.image[IDout].arrayCD[jj*naxes[0]+(naxes[0]-ii)].re = data.image[IDtmp].arrayCD[jj*naxestmp[0]+ii].re;
	  data.image[IDout].arrayCD[jj*naxes[0]+(naxes[0]-ii)].im = -data.image[IDtmp].arrayCD[jj*naxestmp[0]+ii].im;
	  for(jj=1;jj<naxes[1];jj++)
	    {
	      data.image[IDout].arrayCD[jj*naxes[0]+(naxes[0]-ii)].re = data.image[IDtmp].arrayCD[(naxes[1]-jj)*naxestmp[0]+ii].re;
	      data.image[IDout].arrayCD[jj*naxes[0]+(naxes[0]-ii)].im = -data.image[IDtmp].arrayCD[(naxes[1]-jj)*naxestmp[0]+ii].im;
	    }
	}
    }

  if(OK==0)
    {
      printf("Error : image dimension not appropriate for FFT\n");
    }
  free(naxes);
  free(naxesl);
  free(naxestmp);
  delete_image_ID(ffttmpname);
  
  return(0);
}

/* 1d inverse complex fft */
int do1dffti(char *in_name, char *out_name)
{
  int *naxes;
  long *naxesl;
  long naxis;
  long IDin,IDout;
  long i;
  int OK=0;
  fftwf_plan plan;

  IDin=image_ID(in_name);
  naxis=data.image[IDin].naxis;
  naxes = (int *) malloc(naxis*sizeof(int));
  naxesl = (long *) malloc(naxis*sizeof(long));
  for (i=0;i<naxis;i++)
    {
      naxesl[i]= (long) data.image[IDin].size[i];
      naxes[i]= (int) data.image[IDin].size[i];
    }
  IDout = create_image_ID(out_name,naxis,naxesl,CDtype);

  if(naxis==1)
    {
      if(array_index(naxes[0])!=100)
	{
	  OK=1;
	  plan = fftwf_plan_dft_1d(naxes[0], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, 1,FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}
      else
	{
	  OK=1;
	  plan = fftwf_plan_dft_1d(naxes[0], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, 1,FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}
    }

  if(naxis==2)
    {
      if((naxes[1]==1)&&(array_index(naxes[0])!=100))
	{
	  OK=1;
	  plan = fftwf_plan_dft_1d(naxes[0], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, 1,FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}
      else
	{
	  OK=1;
	  plan = fftwf_plan_dft_1d(naxes[0], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, 1,FFTWOPTMODE);
	  fftwf_execute(plan);
	  fftwf_destroy_plan(plan);
	}	
    }

  if(OK==0)
    {
      printf("Error : image dimension not appropriate for FFT\n");
    }
  free(naxes);

  return(0);
}

/* 2d complex fft */
int do2dfft(char *in_name, char *out_name)
{
  int *naxes;
  long *naxesl;
  long naxis;
  long IDin,IDout;
  long i;
  int OK=0;
  fftwf_plan plan;
  long tmp1;
  long IDcpy;

  char ffttmpcpyname[SBUFFERSIZE];
  int n;
  long nextID;


  IDin = image_ID(in_name);
  naxis = data.image[IDin].naxis;
  naxes = (int *) malloc(naxis*sizeof(int));
  naxesl = (long *) malloc(naxis*sizeof(long));

  for (i=0;i<naxis;i++)
    {
      naxesl[i]= (long) data.image[IDin].size[i];
      naxes[i]= (int) data.image[IDin].size[i];
    }


  IDout = create_image_ID(out_name,naxis,naxesl,CDtype);
  
  // need to swap first 2 axis for fftw
  if(naxis>1)
    {
      tmp1 = naxes[0];
      naxes[0] = naxes[1];
      naxes[1] = tmp1;
    }
  if(naxis==2)
    {
      OK=1;
      plan = fftwf_plan_dft_2d(naxes[0],naxes[1], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, -1,FFTWOPTMODE);
      
      if(plan==NULL)
	{
	  //	  if ( Debug > 2)  
	  fprintf(stdout,"New FFT size [do2dfft %d x %d]: optimizing ...",naxes[1],naxes[0]);
	  fflush(stdout);	  

	  n = snprintf(ffttmpcpyname,SBUFFERSIZE,"_ffttmpcpyname_%d",(int) getpid());
	  if(n >= SBUFFERSIZE) 
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  copy_image_ID(in_name,ffttmpcpyname);

	  plan = fftwf_plan_dft_2d(naxes[0],naxes[1], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, -1,FFTWOPTMODE);
	  copy_image_ID(ffttmpcpyname,in_name);
	  delete_image_ID(ffttmpcpyname);
	  export_wisdom();
	  fprintf(stdout,"\n");
	}

      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
    }

  if(naxis==3)
    {
      OK=1;
     
      plan = fftwf_plan_many_dft(2,naxes,naxes[2],(fftwf_complex*) data.image[IDin].arrayCD,NULL,1,naxes[0]*naxes[1],(fftwf_complex*) data.image[IDout].arrayCD,NULL,1,naxes[0]*naxes[1],-1,FFTWOPTMODE);
      if(plan==NULL)
	{
	  //if ( Debug > 2) 
	  fprintf(stdout,"New FFT size [do2dfft %d x %d x %d]: optimizing ...",naxes[1],naxes[0],naxes[2]);
	  fflush(stdout);

	  n = snprintf(ffttmpcpyname,SBUFFERSIZE,"_ffttmpcpyname_%d",(int) getpid());
	  if(n >= SBUFFERSIZE) 
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  copy_image_ID(in_name,ffttmpcpyname);

	  plan = fftwf_plan_many_dft(2,naxes,naxes[2],(fftwf_complex*) data.image[IDin].arrayCD,NULL,1,naxes[0]*naxes[1],(fftwf_complex*) data.image[IDout].arrayCD,NULL,1,naxes[0]*naxes[1],-1,FFTWOPTMODE);
	  copy_image_ID(ffttmpcpyname,in_name);
	  delete_image_ID(ffttmpcpyname);
	  export_wisdom();
	  fprintf(stdout,"\n");
	}

      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
    }

  if(OK==0)
    printf("Error : image dimension not appropriate for FFT\n");
    
  free(naxes);

  return(0);
}


/* 2d inverse complex fft */
int do2dffti(char *in_name, char *out_name)
{
  long *naxesl;
  int *naxes;
  long naxis;
  long IDin,IDout;
  long i;
  int OK=0;
  fftwf_plan plan;
  long tmp1;
  
  char ffttmpcpyname[SBUFFERSIZE];
  int n;

  IDin=image_ID(in_name);
  naxis=data.image[IDin].naxis;
  naxes = (int *) malloc(naxis*sizeof(int));
  naxesl = (long *) malloc(naxis*sizeof(long));
  for (i=0;i<naxis;i++)
    {
      naxes[i]= (int) data.image[IDin].size[i];
      naxesl[i]= (long) data.image[IDin].size[i];
    }


  IDout = create_image_ID(out_name,naxis,naxesl,CDtype);

  // need to swap first 2 axis for fftw
  if(naxis>1)
    {
      tmp1 = naxes[0];
      naxes[0] = naxes[1];
      naxes[1] = tmp1;
    }


  //  IDout = create_image_ID(out_name,naxis,naxesl,CDtype);

  if(naxis==2)
    {
      OK=1;
      plan = fftwf_plan_dft_2d(naxes[0],naxes[1], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, 1,FFTWOPTMODE);
      if(plan==NULL)
	{
	  //	  if ( Debug > 2) 
	  fprintf(stdout,"New FFT size [do2dffti %d x %d]: optimizing ...",naxes[1],naxes[0]);
	  fflush(stdout);

	  n = snprintf(ffttmpcpyname,SBUFFERSIZE,"_ffttmpcpyname_%d",(int) getpid());
	  if(n >= SBUFFERSIZE) 
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  copy_image_ID(in_name,ffttmpcpyname);
	  
	  plan = fftwf_plan_dft_2d(naxes[0],naxes[1], (fftwf_complex*) data.image[IDin].arrayCD, (fftwf_complex*) data.image[IDout].arrayCD, 1,FFTWOPTMODE);
	  copy_image_ID(ffttmpcpyname,in_name);
	  delete_image_ID(ffttmpcpyname);
	  export_wisdom();
	  fprintf(stdout,"\n");	  
	}
      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
    }
  if(naxis==3)
    {
      OK=1;
      plan = fftwf_plan_many_dft(2,naxes,naxes[2],(fftwf_complex*) data.image[IDin].arrayCD,NULL,1,naxes[0]*naxes[1],(fftwf_complex*) data.image[IDout].arrayCD,NULL,1,naxes[0]*naxes[1],1,FFTWOPTMODE);
      if(plan==NULL)
	{
	  // if ( Debug > 2) 
	  fprintf(stdout,"New FFT size [do2dffti %d x %d x %d]: optimizing ...",naxes[1],naxes[0],naxes[2]);
	  fflush(stdout);

	  n = snprintf(ffttmpcpyname,SBUFFERSIZE,"_ffttmpcpyname_%d",(int) getpid());
	  if(n >= SBUFFERSIZE) 
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  copy_image_ID(in_name,ffttmpcpyname);
	    
	  plan = fftwf_plan_many_dft(2,naxes,naxes[2],(fftwf_complex*) data.image[IDin].arrayCD,NULL,1,naxes[0]*naxes[1],(fftwf_complex*) data.image[IDout].arrayCD,NULL,1,naxes[0]*naxes[1],1,FFTWOPTMODE);
	  copy_image_ID(ffttmpcpyname,in_name);
	  delete_image_ID(ffttmpcpyname);
	  export_wisdom();
	  fprintf(stdout,"\n");
	}
      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
    }

  if(OK==0)
    printf("Error (do2dffti): image dimension not appropriate for FFT: naxis = %ld\n",naxis);

  free(naxes);

  return(0);
}



  /* inv = 0 for direct fft and 1 for inverse fft */
  /* direct = focal plane -> pupil plane  equ. fft2d(..,..,..,1) */
  /* inverse = pupil plane -> focal plane equ. fft2d(..,..,..,0) */
  /* options :  -reim  takes real/imaginary input and creates real/imaginary output
                 -inv  for inverse fft (inv=1) */
int pupfft(char *ID_name_ampl, char *ID_name_pha, char *ID_name_ampl_out, char *ID_name_pha_out, char *options)
{
  int reim;
  int inv;
  
  char Ctmpname[SBUFFERSIZE];
  char C1tmpname[SBUFFERSIZE];
  int n;

  reim = 0;
  inv = 0;
  
  if(strstr(options,"-reim")!=NULL)
    {
      /*	printf("taking real / imaginary input/output\n");*/
      reim = 1;
    }
  
  if(strstr(options,"-inv")!=NULL)
    {
      /*printf("doing the inverse Fourier transform\n");*/
      inv = 1;
    }


  n = snprintf(Ctmpname,SBUFFERSIZE,"_Ctmp_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
  if (reim==0)
    {
      mk_complex_from_amph(ID_name_ampl,ID_name_pha,Ctmpname);
    }   
  else
    {
      mk_complex_from_reim(ID_name_ampl,ID_name_pha,Ctmpname);
    }
  
    permut(Ctmpname);

    n = snprintf(C1tmpname,SBUFFERSIZE,"_C1tmp_%d",(int) getpid());
    if(n >= SBUFFERSIZE) 
      printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
    if(inv==0)
      do2dfft(Ctmpname,C1tmpname); /* equ. fft2d(..,1) */
    else
      do2dffti(Ctmpname,C1tmpname); /* equ. fft2d(..,0) */

    delete_image_ID(Ctmpname);

    if (reim==0)
      {
	/* if this line is removed, the program crashes... why ??? */
	/*	list_image_ID(data); */
	mk_amph_from_complex(C1tmpname,ID_name_ampl_out,ID_name_pha_out);
      }
    else
      {
	mk_reim_from_complex(C1tmpname,ID_name_ampl_out,ID_name_pha_out);
      }
    
    delete_image_ID(C1tmpname);

    permut(ID_name_ampl_out);
    permut(ID_name_pha_out);

    return(0);
  }


/* real fft : real to complex */
int do2drfft(char *in_name, char *out_name)
{
  int *naxes;
  long *naxesl;
  long *naxestmp;

  long naxis;
  long IDin,IDout,IDtmp;
  long i;
  int OK=0;
  long idist;
  long ii,jj,kk;
  fftwf_plan plan;
  long tmp1;

  char ffttmpname[SBUFFERSIZE];
  char ffttmpcpyname[SBUFFERSIZE];
  int n;

  IDin = image_ID(in_name);
  naxis = data.image[IDin].naxis;
  naxes = (int *) malloc(naxis*sizeof(int));
  naxesl = (long *) malloc(naxis*sizeof(long));
  naxestmp = (long *) malloc(naxis*sizeof(long));

  for (i=0;i<naxis;i++)
    {
      naxes[i] = (int) data.image[IDin].size[i];
      naxesl[i] = (long) data.image[IDin].size[i];
      naxestmp[i]=data.image[IDin].size[i];
      if(i==0)
	naxestmp[i] = data.image[IDin].size[i]/2+1;
    }

  n = snprintf(ffttmpname,SBUFFERSIZE,"_ffttmp_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
  IDtmp = create_image_ID(ffttmpname,naxis,naxestmp,CDtype);

  IDout = create_image_ID(out_name,naxis,naxesl,CDtype);

  if(naxis==2)
    {
      OK=1;
      plan = fftwf_plan_dft_r2c_2d(naxes[1],naxes[0],data.image[IDin].arrayD,(fftwf_complex*) data.image[IDtmp].arrayCD,FFTWOPTMODE);
      if(plan==NULL)
	{
	  // if ( Debug > 2) 
	  fprintf(stdout,"New FFT size [do2drfft %d x %d]: optimizing ...",naxes[1],naxes[0]);
	  fflush(stdout);
	  
	  n = snprintf(ffttmpcpyname,SBUFFERSIZE,"_ffttmpcpy_%d",(int) getpid());
	  if(n >= SBUFFERSIZE) 
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  copy_image_ID(in_name,ffttmpcpyname);

	  plan = fftwf_plan_dft_r2c_2d(naxes[1],naxes[0],data.image[IDin].arrayD,(fftwf_complex*) data.image[IDtmp].arrayCD,FFTWOPTMODE);
	  copy_image_ID(ffttmpcpyname,in_name);
	  delete_image_ID(ffttmpcpyname);
	  export_wisdom();
	  fprintf(stdout,"\n");
	}
      fftwf_execute(plan);
      fftwf_destroy_plan(plan);

      for(ii=0;ii<naxes[0]/2+1;ii++)
	for(jj=0;jj<naxes[1];jj++)
	  data.image[IDout].arrayCD[jj*naxes[0]+ii] = data.image[IDtmp].arrayCD[jj*naxestmp[0]+ii];
	  
      for(ii=1;ii<naxes[0]/2+1;ii++)
	{
	  jj=0;
	  data.image[IDout].arrayCD[jj*naxes[0]+(naxes[0]-ii)].re = data.image[IDtmp].arrayCD[jj*naxestmp[0]+ii].re;
	  data.image[IDout].arrayCD[jj*naxes[0]+(naxes[0]-ii)].im = -data.image[IDtmp].arrayCD[jj*naxestmp[0]+ii].im;
	  for(jj=1;jj<naxes[1];jj++)
	    {
	      data.image[IDout].arrayCD[jj*naxes[0]+(naxes[0]-ii)].re = data.image[IDtmp].arrayCD[(naxes[1]-jj)*naxestmp[0]+ii].re;
	      data.image[IDout].arrayCD[jj*naxes[0]+(naxes[0]-ii)].im = -data.image[IDtmp].arrayCD[(naxes[1]-jj)*naxestmp[0]+ii].im;
	    }
	}
    }
  if(naxis==3)
    {
      OK=1;
      idist = naxes[0]*naxes[1];
      
      // swapping first 2 axis
      tmp1 = naxes[0];
      naxes[0] = naxes[1];
      naxes[1] = tmp1;
      
      plan = fftwf_plan_many_dft_r2c(2,naxes,naxes[2],data.image[IDin].arrayD,NULL,1,naxes[0]*naxes[1],(fftwf_complex*) data.image[IDout].arrayCD,NULL,1,naxes[0]*naxes[1],FFTWOPTMODE);
      if(plan==NULL)
	{
	  //	  if ( Debug > 2) fprintf(stdout,"New FFT size [do2drfft %d x %d x %d]: optimizing ...",naxes[1],naxes[0],naxes[2]);
	  fflush(stdout);

	  n = snprintf(ffttmpcpyname,SBUFFERSIZE,"_ffttmpcpy_%d",(int) getpid());
	  if(n >= SBUFFERSIZE) 
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  copy_image_ID(in_name,ffttmpcpyname);

	  plan = fftwf_plan_many_dft_r2c(2,naxes,naxes[2],data.image[IDin].arrayD,NULL,1,naxes[0]*naxes[1],(fftwf_complex*) data.image[IDout].arrayCD,NULL,1,naxes[0]*naxes[1],FFTWOPTMODE);
	  copy_image_ID(ffttmpcpyname,in_name);
	  delete_image_ID(ffttmpcpyname);
	  export_wisdom();
	  fprintf(stdout,"\n");
	}

      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
      
      // unswapping first 2 axis
      tmp1 = naxes[0];
      naxes[0] = naxes[1];
      naxes[1] = tmp1;
      
      for(ii=0;ii<naxes[0]/2+1;ii++)
	for(jj=0;jj<naxes[1];jj++)
	  for(kk=0;kk<naxes[2];kk++)
	    {
	      data.image[IDout].arrayCD[naxes[0]*naxes[1]*kk+jj*naxes[0]+ii]=data.image[IDtmp].arrayCD[naxestmp[0]*naxestmp[1]*kk+jj*naxestmp[0]+ii];
	      if(ii!=0)
		data.image[IDout].arrayCD[naxes[0]*naxes[1]*kk+jj*naxes[0]+(naxes[0]-ii)]=data.image[IDtmp].arrayCD[naxestmp[0]*naxestmp[1]*kk+jj*naxestmp[0]+ii];
	    }
    }
  
  if(OK==0)
    printf("Error : image dimension not appropriate for FFT\n");
    
  delete_image_ID(ffttmpname);

  free(naxestmp);
  free(naxesl);
  free(naxes);

  return(0);
}

/* inverse real fft : complex to real */
int do2drffti(char *in_name, char *out_name)
{
  int *naxes;
  long *naxesl;
  long naxis;
  long IDin,IDout;
  long i;
  int OK=0;
  long idist;
  fftwf_plan plan;
  long tmp1;

  char ffttmpcpyname[SBUFFERSIZE];
  int n;

  IDin = image_ID(in_name);
  naxis = data.image[IDin].naxis;
  naxes = (int *) malloc(naxis*sizeof(int));
  naxesl = (long *) malloc(naxis*sizeof(long));

  for (i=0;i<naxis;i++)
    {
      naxes[i]=data.image[IDin].size[i];
      naxesl[i]=data.image[IDin].size[i];
    }
  IDout = create_image_ID(out_name,naxis,naxesl,Dtype);

  if(naxis==2)
    {
      OK=1;
      plan = fftwf_plan_dft_r2c_2d(naxes[1],naxes[0],data.image[IDin].arrayD,(fftwf_complex*) data.image[IDout].arrayCD,FFTWOPTMODE);
      if(plan==NULL)
	{
	  //if ( Debug > 2) 
	  fprintf(stdout,"New FFT size [do2drffti %d x %d]: optimizing ...",naxes[0],naxes[1]);
	  fflush(stdout);

	  n = snprintf(ffttmpcpyname,SBUFFERSIZE,"_ffttmpcpy_%d",(int) getpid());
	  copy_image_ID(in_name,ffttmpcpyname);

	  plan = fftwf_plan_dft_r2c_2d(naxes[1],naxes[0],data.image[IDin].arrayD,(fftwf_complex*) data.image[IDout].arrayCD,FFTWOPTMODE);
 	  copy_image_ID(ffttmpcpyname,in_name);
	  delete_image_ID(ffttmpcpyname);
	  export_wisdom();
	  fprintf(stdout,"\n");
	}

      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
      
    }
  if(naxis==3)
    {
      OK=1;
      idist = naxes[0]*naxes[1];
      
      // swapping first 2 axis
      tmp1 = naxes[0];
      naxes[0] = naxes[1];
      naxes[1] = tmp1;
      
      plan = fftwf_plan_many_dft_r2c(2,naxes,naxes[2],data.image[IDin].arrayD,NULL,1,naxes[0]*naxes[1],(fftwf_complex*) data.image[IDout].arrayCD,NULL,1,naxes[0]*naxes[1],FFTWOPTMODE);
      if(plan==NULL)
	{
	  //if ( Debug > 2) 
	  fprintf(stdout,"New FFT size [do2drffti %d x %d x %d]: optimizing ...",naxes[1],naxes[0],naxes[2]);
	  fflush(stdout);

	  n = snprintf(ffttmpcpyname,SBUFFERSIZE,"_ffttmpcpy_%d",(int) getpid());
	  if(n >= SBUFFERSIZE) 
	    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
	  copy_image_ID(in_name,ffttmpcpyname);

	  plan = fftwf_plan_many_dft_r2c(2,naxes,naxes[2],data.image[IDin].arrayD,NULL,1,naxes[0]*naxes[1],(fftwf_complex*) data.image[IDout].arrayCD,NULL,1,naxes[0]*naxes[1],FFTWOPTMODE);
 	  copy_image_ID(ffttmpcpyname,in_name);
	  delete_image_ID(ffttmpcpyname);
	  export_wisdom();
	  fprintf(stdout,"\n");
	}

      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
    }
  
  if(OK==0)
    printf("Error : image dimension not appropriate for FFT\n");
    
  free(naxes);
  free(naxesl);

  return(0);
}


long fft_correlation(char *ID_name1, char *ID_name2, char *ID_nameout)
{
  long ID1,ID2,IDout;
  long nelement;
  
  char ft1name[SBUFFERSIZE];
  char ft2name[SBUFFERSIZE];
  char fta1name[SBUFFERSIZE];
  char fta2name[SBUFFERSIZE];
  char ftp1name[SBUFFERSIZE];
  char ftp2name[SBUFFERSIZE];
  char fta12name[SBUFFERSIZE];
  char ftp12name[SBUFFERSIZE];
  char fftname[SBUFFERSIZE];
  char fft1name[SBUFFERSIZE];
  char fft1pname[SBUFFERSIZE];
  int n;
  
  ID1 = image_ID(ID_name1);
  nelement = data.image[ID1].nelement;  


  n = snprintf(ft1name,SBUFFERSIZE,"_ft1_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
  do2drfft(ID_name1,ft1name);
  
  n = snprintf(ft2name,SBUFFERSIZE,"_ft2_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
  do2drfft(ID_name2,ft2name);

  n = snprintf(fta1name,SBUFFERSIZE,"_%s_a_%d",ft1name,(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
  
  n = snprintf(ftp1name,SBUFFERSIZE,"_%s_p_%d",ft1name,(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  n = snprintf(fta2name,SBUFFERSIZE,"_%s_a_%d",ft2name,(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  n = snprintf(ftp2name,SBUFFERSIZE,"_%s_p_%d",ft2name,(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  n = snprintf(fta12name,SBUFFERSIZE,"_%s_12a_%d",ft1name,(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  n = snprintf(ftp12name,SBUFFERSIZE,"_%s_12p_%d",ft1name,(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  mk_amph_from_complex(ft1name,fta1name,ftp1name);
  mk_amph_from_complex(ft2name,fta2name,ftp2name);
  delete_image_ID(ft1name);
  delete_image_ID(ft2name);

  arith_image_mult(fta1name, fta2name, fta12name);
  arith_image_sub(ftp1name, ftp2name, ftp12name);
  delete_image_ID(fta1name);
  delete_image_ID(fta2name);
  delete_image_ID(ftp1name);
  delete_image_ID(ftp2name);
  
  arith_image_cstmult_inplace(fta12name,1.0/sqrt(nelement)/(1.0*nelement));
  
  n = snprintf(fftname,SBUFFERSIZE,"_fft_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  mk_complex_from_amph(fta12name,ftp12name,fftname);
  delete_image_ID(fta12name);
  delete_image_ID(ftp12name);
 
  n = snprintf(fft1name,SBUFFERSIZE,"_fft1_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  do2dfft(fftname,fft1name);
  delete_image_ID(fftname);

  n = snprintf(fft1pname,SBUFFERSIZE,"_fft1p_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  mk_amph_from_complex(fft1name,ID_nameout,fft1pname);
  permut(ID_nameout);
  delete_image_ID(fft1name);
  delete_image_ID(fft1pname);


  IDout = image_ID(ID_nameout);

  return(IDout);
}


int autocorrelation(char *ID_name, char *ID_out)
{
  long ID;
  long nelement;

  char atmp1name[SBUFFERSIZE];
  char aampname[SBUFFERSIZE];
  char aphaname[SBUFFERSIZE];
  char sqaampname[SBUFFERSIZE];
  char sqaamp1name[SBUFFERSIZE];
  int n;

  ID = image_ID(ID_name);
  nelement = data.image[ID].nelement;

  n = snprintf(atmp1name,SBUFFERSIZE,"_atmp1_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");
 
  do2drfft(ID_name,atmp1name);
  n = snprintf(aampname,SBUFFERSIZE,"_aamp_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  n = snprintf(aphaname,SBUFFERSIZE,"_apha_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  mk_amph_from_complex(atmp1name,aampname,aphaname);
 
  n = snprintf(sqaampname,SBUFFERSIZE,"_sqaamp_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  arith_image_mult(aampname,aampname,sqaampname);
  delete_image_ID(aampname);
  delete_image_ID(aphaname);
  delete_image_ID(atmp1name);

  n = snprintf(sqaamp1name,SBUFFERSIZE,"_sqaamp1_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  arith_image_cstmult(sqaampname,1.0/sqrt(nelement)/(1.0*nelement),sqaamp1name);
  delete_image_ID(sqaampname);
  do2drfft(sqaamp1name,atmp1name);
  mk_reim_from_complex(atmp1name,ID_out,aphaname);
  delete_image_ID(sqaamp1name);
  delete_image_ID(atmp1name);
  delete_image_ID(aphaname);
 
  return(0);
}

int fftczoom(char *ID_name, char *ID_out, long factor)
{
  long ID,ID1;
  long naxes[2];
  long ii,jj;
  PRECISION coeff;

  char tmpzname[SBUFFERSIZE];
  char tmpz1name[SBUFFERSIZE];
  int n;

  ID=image_ID(ID_name);

  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];

  coeff = 1.0/(factor*factor*naxes[0]*naxes[1]);
  permut(ID_name);

  n = snprintf(tmpzname,SBUFFERSIZE,"_tmpz_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  do2dfft(ID_name,tmpzname);


  permut(ID_name);
  permut(tmpzname);
  ID = image_ID(tmpzname);

  n = snprintf(tmpz1name,SBUFFERSIZE,"_tmpz1_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  ID1 = create_2DCimage_ID(tmpz1name,factor*naxes[0],factor*naxes[1]);

  for(ii=0;ii<naxes[0];ii++)
     for(jj=0;jj<naxes[1];jj++)
       {
	 data.image[ID1].arrayCD[(jj+factor*naxes[1]/2-naxes[1]/2)*naxes[0]*factor+(ii+factor*naxes[0]/2-naxes[0]/2)].re = data.image[ID].arrayCD[jj*naxes[0]+ii].re*coeff;
	 data.image[ID1].arrayCD[(jj+factor*naxes[1]/2-naxes[1]/2)*naxes[0]*factor+(ii+factor*naxes[0]/2-naxes[0]/2)].im = data.image[ID].arrayCD[jj*naxes[0]+ii].im*coeff;
       }
  delete_image_ID(tmpzname);

  permut(tmpz1name);
  do2dffti(tmpz1name,ID_out);
  permut(ID_out);
  delete_image_ID(tmpz1name);

  return(0);
}

int fftzoom(char *ID_name, char *ID_out, long factor)
{
  long ID,ID1;
  long naxes[2];
  long ii,jj;
  PRECISION coeff;

  char tmpzname[SBUFFERSIZE];
  char tmpz1name[SBUFFERSIZE];
  char tmpz2name[SBUFFERSIZE];
  char tbename[SBUFFERSIZE];
  int n;

  ID = image_ID(ID_name);

  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];

  coeff = 1.0/(factor*factor*naxes[0]*naxes[1]);
  permut(ID_name);

  n = snprintf(tmpzname,SBUFFERSIZE,"_tmpz_%d",(int) getpid());
    if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  do2drfft(ID_name,tmpzname);

  permut(ID_name);
  permut(tmpzname);
  ID = image_ID(tmpzname);

  n = snprintf(tmpz1name,SBUFFERSIZE,"_tmpz1_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  ID1 = create_2DCimage_ID(tmpz1name,factor*naxes[0],factor*naxes[1]);

  for(ii=0;ii<naxes[0];ii++)
     for(jj=0;jj<naxes[1];jj++)
       {
	 data.image[ID1].arrayCD[(jj+factor*naxes[1]/2-naxes[1]/2)*naxes[0]*factor+(ii+factor*naxes[0]/2-naxes[0]/2)].re = data.image[ID].arrayCD[jj*naxes[0]+ii].re*coeff;
	 data.image[ID1].arrayCD[(jj+factor*naxes[1]/2-naxes[1]/2)*naxes[0]*factor+(ii+factor*naxes[0]/2-naxes[0]/2)].im = data.image[ID].arrayCD[jj*naxes[0]+ii].im*coeff;
       }
  delete_image_ID(tmpzname);

  permut(tmpz1name);

  n = snprintf(tmpz2name,SBUFFERSIZE,"_tmpz2_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  do2dffti(tmpz1name,tmpz2name);

  permut(tmpz2name);
  delete_image_ID(tmpz1name);

  n = snprintf(tbename,SBUFFERSIZE,"_tbe_%d",(int) getpid());
  if(n >= SBUFFERSIZE) 
    printERROR(__FILE__,__func__,__LINE__,"Attempted to write string buffer with too many characters");

  mk_reim_from_complex(tmpz2name,ID_out,tbename);

  delete_image_ID(tbename);
  delete_image_ID(tmpz2name);

  return(0);
}

int test_fftspeed(int nmax)
{
  int n;
  long size;
  int nbiter,iter;

  struct timespec tS0;
  struct timespec tS1;
  struct timespec tS2;
  double ti0,ti1,ti2;
  double dt1;
  struct timeval tv;
  int nb_threads=1;
  int nb_threads_max = 8;
  
  /*  printf("%ld ticks per second\n",CLOCKS_PER_SEC);*/
  nbiter = 10000;
  size=2;

  printf("Testing complex FFT, nxn pix\n");
  
  printf("size(pix)");
# ifdef FFTWMT
  for(nb_threads=1;nb_threads<nb_threads_max;nb_threads++)	
    printf("%13d",nb_threads);
# endif
  printf("\n");
      

  size = 2;
  for(n=0;n<nmax;n++)
    {
      printf("%9ld",size);
# ifdef FFTWMT
      for(nb_threads=1;nb_threads<nb_threads_max;nb_threads++)
	{
	  fft_setNthreads(nb_threads);
# endif

	  #if _POSIX_TIMERS > 0
	  clock_gettime(CLOCK_REALTIME, &tS0);  
          #else
	  gettimeofday(&tv, NULL);
	  tS0.tv_sec = tv.tv_sec;
	  tS0.tv_nsec = tv.tv_usec*1000;
          #endif



	  //	  clock_gettime(CLOCK_REALTIME, &tS0);
	  for(iter=0;iter<nbiter;iter++)
	    {  
	      create_2DCimage_ID("tmp",size,size);
	      do2dfft("tmp","tmpf");
	      delete_image_ID("tmp");
	      delete_image_ID("tmpf");
	    }

          #if _POSIX_TIMERS > 0
	  clock_gettime(CLOCK_REALTIME, &tS1);  
          #else
	  gettimeofday(&tv, NULL);
	  tS1.tv_sec = tv.tv_sec;
	  tS1.tv_nsec = tv.tv_usec*1000;
          #endif
	  //	  clock_gettime(CLOCK_REALTIME, &tS1);
	  
	  for(iter=0;iter<nbiter;iter++)
	    {  
	      create_2DCimage_ID("tmp",size,size);
	      delete_image_ID("tmp");
	    }

          #if _POSIX_TIMERS > 0
	  clock_gettime(CLOCK_REALTIME, &tS2);  
          #else
	  gettimeofday(&tv, NULL);
	  tS2.tv_sec = tv.tv_sec;
	  tS2.tv_nsec = tv.tv_usec*1000;
          #endif
	  //clock_gettime(CLOCK_REALTIME, &tS2);
	  
	  
	  ti0 = 1.0*tS0.tv_sec+0.000000001*tS0.tv_nsec;
	  ti1 = 1.0*tS1.tv_sec+0.000000001*tS1.tv_nsec;
	  ti2 = 1.0*tS2.tv_sec+0.000000001*tS2.tv_nsec;
	  dt1 = 1.0*(ti1-ti0)-1.0*(ti2-ti1);

	  dt1 /= nbiter;
	  
	  printf("%10.3f ms",dt1*1000.0);
	  //printf("Complex FFT %ldx%ld [%d threads] : %f ms  [%ld]\n",size,size,nb_threads,dt1*1000.0,nbiter);
	  fflush(stdout);
# ifdef FFTWMT
	}
# endif
      printf("\n");
      nbiter = 0.1/dt1;
      if(nbiter<2)
	nbiter = 2;
      size = size*2;
    }

  return(0);
}






/* ----------------- CUSTOM DFT ------------- */

// Zfactor is zoom factor
// dir = -1 for FT, 1 for inverse FT

long COREMOD_fft_DFT( char *IDin_name, char *IDinmask_name, char *IDout_name, char *IDoutmask_name, double Zfactor, int dir)
{
  long IDin;
  long IDout;
  long IDinmask;
  long IDoutmask;

  long NBptsin;
  long NBptsout;

  long xsize, ysize;
  long ii, jj, k, kout;
  double val;
  double re, im, pha;

  long *iiinarray;
  long *jjinarray;
  double *xinarray;
  double *yinarray;
  double *valinamp;
  double *valinpha;

  long *iioutarray;
  long *jjoutarray;
  double *xoutarray;
  double *youtarray;

  IDin = image_ID(IDin_name);

  IDinmask = image_ID(IDinmask_name);
  xsize = data.image[IDinmask].size[0];
  ysize = data.image[IDinmask].size[1];
  NBptsin = 0;
  for(ii=0; ii<xsize; ii++)
    for(jj=0; jj<ysize; jj++)
      {
	val = data.image[IDinmask].arrayD[jj*xsize+ii];
	if(val>0.5)
	  NBptsin ++;
      }

  printf("DFT (%f):  %ld input points -> ", Zfactor, NBptsin);

  iiinarray = (long *) malloc(sizeof(long)*NBptsin);
  jjinarray = (long *) malloc(sizeof(long)*NBptsin);
  xinarray = (double *) malloc(sizeof(double)*NBptsin);
  yinarray = (double *) malloc(sizeof(double)*NBptsin);
  valinamp = (double *) malloc(sizeof(double)*NBptsin);
  valinpha = (double *) malloc(sizeof(double)*NBptsin);
  k = 0;


  for(ii=0; ii<xsize; ii++)
    for(jj=0; jj<ysize; jj++)
      {
	val = data.image[IDinmask].arrayD[jj*xsize+ii];
	if(val>0.5)
	  {
	    iiinarray[k] = ii;
	    jjinarray[k] = jj;
	    xinarray[k] = 1.0*ii/xsize-0.5;
	    yinarray[k] = 1.0*jj/xsize-0.5;
	    re = data.image[IDin].arrayCD[jj*xsize+ii].re;
	    im = data.image[IDin].arrayCD[jj*xsize+ii].im;
	    valinamp[k] = sqrt(re*re+im*im);
	    valinpha[k] = atan2(im,re);
	    k++;
	  }
      }

  


  IDoutmask = image_ID(IDoutmask_name);
  NBptsout = 0;
  for(ii=0; ii<xsize; ii++)
    for(jj=0; jj<ysize; jj++)
      {
	val = data.image[IDoutmask].arrayD[jj*xsize+ii];
	if(val>0.5)
	  NBptsout ++;
      }
  
  printf("%ld output points\n", NBptsout);
  
  iioutarray = (long *) malloc(sizeof(long)*NBptsout);
  jjoutarray = (long *) malloc(sizeof(long)*NBptsout);
  xoutarray = (double *) malloc(sizeof(double)*NBptsout);
  youtarray = (double *) malloc(sizeof(double)*NBptsout);
  kout = 0;
  for(ii=0; ii<xsize; ii++)
    for(jj=0; jj<ysize; jj++)
      {
	val = data.image[IDoutmask].arrayD[jj*xsize+ii];
	if(val>0.5)
	  {
	    iioutarray[kout] = ii;
	    jjoutarray[kout] = jj;
	    xoutarray[kout] = (1.0/Zfactor) * (1.0*ii/xsize-0.5) * xsize;
	    youtarray[kout] = (1.0/Zfactor) * (1.0*jj/ysize-0.5) * ysize;
	    kout++;
	  }
      }

  IDout = create_2DCimage_ID(IDout_name, xsize, ysize);


  # ifdef _OPENMP
#pragma omp parallel default(shared) private(kout, k, pha, re, im)
  {
  # endif

  # ifdef _OPENMP
      #pragma omp for
      # endif
  for(kout=0; kout<NBptsout; kout++)
    {
      re = 0.0;
      im = 0.0;
      for(k=0; k<NBptsin; k++)
	{
	  pha = valinpha[k] + 2.0*dir*M_PI*(xinarray[k]*xoutarray[kout] + yinarray[k]*youtarray[kout]);
	  re += valinamp[k]*cos(pha);
	  im += valinamp[k]*sin(pha);
	}
      data.image[IDout].arrayCD[jjoutarray[kout]*xsize+iioutarray[kout]].re = re/Zfactor;
      data.image[IDout].arrayCD[jjoutarray[kout]*xsize+iioutarray[kout]].im = im/Zfactor;
    }
 # ifdef _OPENMP
  }
  # endif



  free(iiinarray);
  free(jjinarray);
  free(xinarray);
  free(yinarray);
  free(valinamp);
  free(valinpha);


  free(iioutarray);
  free(jjoutarray);
  free(xoutarray);
  free(youtarray);

  return(IDout);
}


//
// pupil convolution by complex focal plane mask of limited support
// typically used with fpmz = zoomed copy of 1-fpm
// high resolution focal plane mask using DFT
// zoom factor
//
//
//
long COREMOD_fft_DFTinsertFPM( char *pupin_name, char *fpmz_name, double zfactor, char *pupout_name)
{
  double eps = 1.0e-16;
  long ID, ID1;
  long IDpupin_mask;
  long IDfpmz;
  long IDfpmz_mask;
  long xsize, ysize;
  long IDin, IDout;
  long ii, jj;
  double re, im, rein, imin, amp, pha, ampin, phain, amp2;
  double x, y, r;
  double total = 0;
  
  int FORCE_IMZERO = 0;
  double imresidual = 0.0;
  double tx, ty, tcx, tcy;

  if(variable_ID("_FORCE_IMZERO")!=-1)
    {
      FORCE_IMZERO = 1;
      printf("---------------FORCING IMAGINARY PART TO ZERO-------------\n");
    }
  

  printf("zfactor = %f\n", zfactor);

  IDin = image_ID(pupin_name);
  xsize = data.image[ID].size[0];
  ysize = data.image[ID].size[1];

  IDpupin_mask = create_2Dimage_ID("_pupinmask", xsize, ysize);
  for(ii=0; ii<xsize*ysize; ii++)
    {
      re = data.image[IDin].arrayCD[ii].re;
      im = data.image[IDin].arrayCD[ii].im;
      amp2 = re*re+im*im;
      if(amp2>eps)
	data.image[IDpupin_mask].arrayD[ii] = 1.0;
      else
	data.image[IDpupin_mask].arrayD[ii] = 0.0;
    }

  IDfpmz = image_ID(fpmz_name);
  IDfpmz_mask = create_2Dimage_ID("_fpmzmask", xsize, ysize);
  for(ii=0; ii<xsize*ysize; ii++)
    {
      re = data.image[IDfpmz].arrayCD[ii].re;
      im = data.image[IDfpmz].arrayCD[ii].im;
      amp2 = re*re+im*im;
      if(amp2>eps)
	data.image[IDfpmz_mask].arrayD[ii] = 1.0;
      else
	data.image[IDfpmz_mask].arrayD[ii] = 0.0;
    }

  COREMOD_fft_DFT( pupin_name, "_pupinmask", "_foc0", "_fpmzmask", zfactor, -1);

  ID = image_ID("_foc0");
  total = 0.0;
  tx = 0.0;
  ty = 0.0;
  tcx = 0.0;
  tcy = 0.0;
  for(ii=0; ii<xsize; ii++)
    for(jj=0; jj<ysize; jj++)
      {
	x = 1.0*ii-0.5*xsize;
	y = 1.0*jj-0.5*ysize;
	re = data.image[IDfpmz].arrayCD[jj*xsize+ii].re;
	im = data.image[IDfpmz].arrayCD[jj*xsize+ii].im;
	amp = sqrt(re*re+im*im);
	pha = atan2(im,re);
	
	rein = data.image[ID].arrayCD[jj*xsize+ii].re;
	imin = data.image[ID].arrayCD[jj*xsize+ii].im;
	ampin = sqrt(rein*rein+imin*imin);
	phain = atan2(imin, rein);
	
	ampin *= amp;
	total += ampin*ampin;
	phain += pha;
	
	data.image[ID].arrayCD[jj*xsize+ii].re = ampin*cos(phain);
	data.image[ID].arrayCD[jj*xsize+ii].im = ampin*sin(phain);      

	tx += x*ampin*sin(phain)*ampin;
	ty += y*ampin*sin(phain)*ampin;
	tcx += x*x*ampin*ampin;
	tcy += y*y*ampin*ampin;
      }
  printf("TX TY = %.18lf %.18lf", tx/tcx, ty/tcy);
  if(FORCE_IMZERO==1) // Remove tip-tilt in focal plane mask imaginary part
    {
      tx = 0.0;
      ty = 0.0;
      for(ii=0; ii<xsize; ii++)
	for(jj=0; jj<ysize; jj++)
	  {
	    x = 1.0*ii-0.5*xsize;
	    y = 1.0*jj-0.5*ysize;
	    
	    re = data.image[ID].arrayCD[jj*xsize+ii].re;
	    im = data.image[ID].arrayCD[jj*xsize+ii].im;
	    amp = sqrt(re*re+im*im);
	    
	    data.image[ID].arrayCD[jj*xsize+ii].im -= amp*(x*tx/tcx + y*ty/tcy);
	    tx += x*data.image[ID].arrayCD[jj*xsize+ii].im*amp;
	    ty += y*data.image[ID].arrayCD[jj*xsize+ii].im*amp;
	  }
      printf("  ->   %.18lf %.18lf", tx/tcx, ty/tcy);
    
      mk_amph_from_complex("_foc0","_foc0_amp","_foc0_pha");
      save_fl_fits("_foc0_amp", "!_foc_amp.fits");
      save_fl_fits("_foc0_pha", "!_foc_pha.fits");
      delete_image_ID("_foc0_amp");
      delete_image_ID("_foc0_pha");
    }
  printf("\n");


  data.CFITSVARRAY[0] = (PRECISION) total;

  /*  if(FORCE_IMZERO==1) // Remove tip-tilt in focal plane mask imaginary part
    {
      imresidual = 0.0;
      ID = image_ID("_foc0");
      ID1 = create_2Dimage_ID("imresidual", xsize, ysize);
      for(ii=0; ii<xsize*ysize; ii++)
	{
	  data.image[ID1].arrayD[ii] = data.image[ID].arrayCD[ii].im;
	  imresidual += data.image[ID].arrayCD[ii].im*data.image[ID].arrayCD[ii].im;
	  data.image[ID].arrayCD[ii].im = 0.0;
	}
      printf("IM RESIDUAL = %lf\n", imresidual);
      save_fl_fits("imresidual", "!imresidual.fits");
      delete_image_ID("imresidual");
    }
  */

  if(1) // TEST
    {
      mk_amph_from_complex("_foc0", "tmp_foc0_a", "tmp_foc0_p");
      save_fl_fits("tmp_foc0_a", "!_DFT_foca");
      save_fl_fits("tmp_foc0_p", "!_DFT_focp");
      delete_image_ID("tmp_foc0_a");
      delete_image_ID("tmp_foc0_p");
    }


  /* for(ii=0; ii<xsize; ii++)
    for(jj=0; jj<ysize; jj++)
      {
	x = 1.0*ii-xsize/2;
	y = 1.0*jj-ysize/2;
	r = sqrt(x*x+y*y);
	if(r<150.0)
	  data.image[IDpupin_mask].arrayD[jj*xsize+ii] = 1.0;
	  }*/

  COREMOD_fft_DFT( "_foc0", "_fpmzmask", pupout_name, "_pupinmask", zfactor, 1);

  IDout = image_ID(pupout_name);
  for(ii=0; ii<xsize*ysize; ii++)
    {
      data.image[IDout].arrayCD[ii].re /= xsize*ysize;
      data.image[IDout].arrayCD[ii].im /= xsize*ysize;
    }  

  delete_image_ID("_foc0");

  delete_image_ID("_pupinmask");
  delete_image_ID("_fpmzmask");
  
  return(IDout);
}


//
// pupil convolution by real focal plane mask of limited support
// typically used with fpmz = zoomed copy of 1-fpm
// high resolution focal plane mask using DFT
// zoom factor
//
//
//
long COREMOD_fft_DFTinsertFPM_re( char *pupin_name, char *fpmz_name, double zfactor, char *pupout_name)
{
  double eps = 1.0e-10;
  long ID;
  long IDpupin_mask;
  long IDfpmz;
  long IDfpmz_mask;
  long xsize, ysize;
  long IDin, IDout;
  long ii, jj;
  double re, im, rein, imin, amp, pha, ampin, phain, amp2;
  double x, y, r;
  double total = 0;
  
  IDin = image_ID(pupin_name);
  xsize = data.image[ID].size[0];
  ysize = data.image[ID].size[1];


  printf("zfactor = %f\n", zfactor);

  IDpupin_mask = create_2Dimage_ID("_pupinmask", xsize, ysize);
  for(ii=0; ii<xsize*ysize; ii++)
    {
      re = data.image[IDin].arrayCD[ii].re;
      im = data.image[IDin].arrayCD[ii].im;
      amp2 = re*re+im*im;
      if(amp2>eps)
	data.image[IDpupin_mask].arrayD[ii] = 1.0;
      else
	data.image[IDpupin_mask].arrayD[ii] = 0.0;
    }

  IDfpmz = image_ID(fpmz_name);
  IDfpmz_mask = create_2Dimage_ID("_fpmzmask", xsize, ysize);
  for(ii=0; ii<xsize*ysize; ii++)
    {
      amp = fabs(data.image[IDfpmz].arrayD[ii]);
      if(amp>eps)
	data.image[IDfpmz_mask].arrayD[ii] = 1.0;
      else
	data.image[IDfpmz_mask].arrayD[ii] = 0.0;
    }

  COREMOD_fft_DFT( pupin_name, "_pupinmask", "_foc0", "_fpmzmask", zfactor, -1);

  ID = image_ID("_foc0");
  total = 0.0;
  for(ii=0; ii<xsize*ysize; ii++)
    {
      amp = data.image[IDfpmz].arrayD[ii];

      rein = data.image[ID].arrayCD[ii].re;
      imin = data.image[ID].arrayCD[ii].im;
      ampin = sqrt(rein*rein+imin*imin);
      phain = atan2(imin, rein);

      ampin *= amp;
      total += ampin*ampin;

      data.image[ID].arrayCD[ii].re = ampin*cos(phain);
      data.image[ID].arrayCD[ii].im = ampin*sin(phain);      
    }

  data.CFITSVARRAY[0] = (PRECISION) total;


  if(0) // TEST
    {
      mk_amph_from_complex("_foc0", "tmp_foc0_a", "tmp_foc0_p");
      save_fl_fits("tmp_foc0_a", "!_DFT_foca");
      save_fl_fits("tmp_foc0_p", "!_DFT_focp");
      delete_image_ID("tmp_foc0_a");
      delete_image_ID("tmp_foc0_p");
    }

  /* for(ii=0; ii<xsize; ii++)
    for(jj=0; jj<ysize; jj++)
      {
	x = 1.0*ii-xsize/2;
	y = 1.0*jj-ysize/2;
	r = sqrt(x*x+y*y);
	if(r<150.0)
	  data.image[IDpupin_mask].arrayD[jj*xsize+ii] = 1.0;
	  }*/

  COREMOD_fft_DFT( "_foc0", "_fpmzmask", pupout_name, "_pupinmask", zfactor, 1);

  IDout = image_ID(pupout_name);
  for(ii=0; ii<xsize*ysize; ii++)
    {
      data.image[IDout].arrayCD[ii].re /= xsize*ysize;
      data.image[IDout].arrayCD[ii].im /= xsize*ysize;
    }  

  delete_image_ID("_foc0");

  delete_image_ID("_pupinmask");
  delete_image_ID("_fpmzmask");
  
  return(IDout);
}
