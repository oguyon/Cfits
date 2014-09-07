#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sched.h>


#include <sys/types.h>
#include <sys/file.h>
#include <sys/mman.h>

#include <assert.h>

#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <device_types.h>
#include <pthread.h>




#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"

#include "cudacomp/cudacomp.h"

# ifdef _OPENMP
# include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
# endif




extern DATA data;




// data passed to each thread
typedef struct
{
  int thread_no;
  long numl0;
} THDATA;

double vala = 0.0;

void *compute_function( void *ptr );

/** Need to install process with setuid.  Then, so you aren't running privileged all the time */
uid_t euid_real;
uid_t euid_called;
uid_t suid;



int M = 2000; // number of DM modes
int N = 14400; // number of WFelements


#define TESTMODE 0


// computer memory (host)
float *cMat;
float **cMat_part;
float *wfsVec;
float **wfsVec_part;
float *dmVec;
float *dmVecTMP;
float **dmVec_part;

// GPU memory (device)
int bufferindex = 0; // index for double buffer
float **d_cMat;
float **d_wfsVec;
float **d_dmVec;

// splitting limits
long *Nsize;
long *Noffset;

int NBstreams;
cudaStream_t *stream;
cublasHandle_t *handle;
cudaError_t error;
cublasStatus_t stat;
float alpha = 1.0f;
float beta  = 0.0f;


int orientation = 0;





 
 



// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


//int CUDACOMP_init_cli()
//{
//    if(CLI_checkarg(1,2)==0)
	//CUDACOMP_init();
  //  else
    //    return 1;
//}




int init_CUDACOMP()
{
  strcpy(data.module[data.NBmodule].name,__FILE__);
  strcpy(data.module[data.NBmodule].info,"memory management for images and variables");
  data.NBmodule++;

  strcpy(data.cmd[data.NBcmd].key,"cudacompinit");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = CUDACOMP_init;
  strcpy(data.cmd[data.NBcmd].info,"init CUDA comp");
  strcpy(data.cmd[data.NBcmd].syntax,"no argument");
  strcpy(data.cmd[data.NBcmd].example,"cudacompinit");
  strcpy(data.cmd[data.NBcmd].Ccall,"int CUDACOMP_init()");
  data.NBcmd++;
 
  // add atexit functions here

  return 0;
}





int CUDACOMP_init()
{
  int device;
   int deviceCount;
  struct cudaDeviceProp deviceProp;

    cudaGetDeviceCount(&deviceCount);
    printf("%d devices found\n", deviceCount);
    printf("\n");
    for (device = 0; device < deviceCount; ++device) {
        cudaGetDeviceProperties(&deviceProp, device);
        printf("Device %d [ %20s ]  has compute capability %d.%d.\n",
               device, deviceProp.name, deviceProp.major, deviceProp.minor);
        printf("  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n", (float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
        printf("  (%2d) Multiprocessors\n", deviceProp.multiProcessorCount);
        printf("  GPU Clock rate:                                %.0f MHz (%0.2f GHz)\n", deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);
        printf("\n");
    }



    return(0);
}



