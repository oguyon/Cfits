#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef __MACH__
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct mach_timespec *t){
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    uint64_t time;
    time = mach_absolute_time();
    double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
    double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1e9);
    t->tv_sec = seconds;
    t->tv_nsec = nseconds;
    return 0;
}
#else
#include <time.h>
#endif



#include <math.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <err.h>
#include <fcntl.h>
#include <sched.h>
#include <ncurses.h>
#include <semaphore.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <pthread.h>


#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "linopt_imtools/linopt_imtools.h"
#include "AOloopControl/AOloopControl.h"
#include "image_filter/image_filter.h"
#include "info/info.h"
#include "ZernikePolyn/ZernikePolyn.h"
#include "linopt_imtools/linopt_imtools.h"
#include "image_gen/image_gen.h"
#include "statistic/statistic.h"

#include "FPAOloopControl/FPAOloopControl.h"

#ifdef HAVE_CUDA
#include "cudacomp/cudacomp.h"
#endif




# ifdef _OPENMP
# include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000 
# endif








extern DATA data;


#define FPAOconfname "/tmp/FPAOconf.shm"
FPAOLOOPCONTROL_CONF *FPAOconf; // configuration - this can be an array
long NB_FPAOloopcontrol = 1;
int FPAOconf_init = 0;

int *FPAOconf_fd; 





// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string, not existing image
// 4: existing image
// 5: string 




int FPAOloopControl_initMem_cli()
{
    if(CLI_checkarg(1,4)+CLI_checkarg(2,4)+CLI_checkarg(3,4)+CLI_checkarg(4,4)==0)
        FPAOloopControl_initMem(0, data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string, data.cmdargtoken[4].val.string, 0);
    else
        return 1;     
}




int init_FPAOloopControl()
{


    strcpy(data.cmd[data.NBcmd].key,"fpaolinitmem");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = FPAOloopControl_initMem_cli;
    strcpy(data.cmd[data.NBcmd].info,"Initialize focal plane AO memory");
    strcpy(data.cmd[data.NBcmd].syntax,"<dmRM stream> <dmC stream> <FPim stream> <FPim dark stream>");
    strcpy(data.cmd[data.NBcmd].example,"fpaolinitmem dmRM dmC fpim fpimdark");
    strcpy(data.cmd[data.NBcmd].Ccall,"long FPAOloopControl_initMem(long loop, char *IDdmRM_name, char *IDdmC_name, char *IDfpim_name, char *IDfpim_dark_name, int mode)");
    data.NBcmd++;


    strcpy(data.module[data.NBmodule].name, __FILE__);
    strcpy(data.module[data.NBmodule].info, "FP AO loop control");
    data.NBmodule++;


    return 0;
}







/*** mode = 0 or 1. if mode == 1, simply connect */


long FPAOloopControl_initMem(long loop, char *IDdmRM_name, char *IDdmC_name, char *IDfpim_name, char *IDfpim_dark_name, int mode)
{
    int SM_fd;
    struct stat file_stat;
    int create = 0;
    int result;
    long *sizearray;
    char cntname[200];
    int k;
    FILE *fp;
    // FILE *fp1; // testing
    int tmpi;
    int ret;
    char fname[200];


    SM_fd = open(FPAOconfname, O_RDWR);
    if(SM_fd==-1)
    {
        printf("Cannot import file \"%s\" -> creating file\n", FPAOconfname);
        create = 1;
    }
    else
    {
        fstat(SM_fd, &file_stat);
        printf("File %s size: %zd\n", FPAOconfname, file_stat.st_size);
        if(file_stat.st_size!=sizeof(FPAOLOOPCONTROL_CONF)*NB_FPAOloopcontrol)
        {
            printf("File \"%s\" size is wrong -> recreating file\n", FPAOconfname);
            create = 1;
            close(SM_fd);
        }
    }

    if(create==1)
    {
        SM_fd = open(FPAOconfname, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);

        if (SM_fd == -1) {
            perror("Error opening file for writing");
            exit(0);
        }

        result = lseek(SM_fd, sizeof(FPAOLOOPCONTROL_CONF)*NB_FPAOloopcontrol-1, SEEK_SET);
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
    }


	
	
    FPAOconf = (FPAOLOOPCONTROL_CONF*) mmap(0, sizeof(FPAOLOOPCONTROL_CONF)*NB_FPAOloopcontrol, PROT_READ | PROT_WRITE, MAP_SHARED, SM_fd, 0);
    if (FPAOconf == MAP_FAILED) {
        close(SM_fd);
        perror("Error mmapping the file");
        exit(0);
    }
    
   
	// DM streams
	
	FPAOconf[loop].dmRM_ID = image_ID(IDdmRM_name);
	FPAOconf[loop].dmC_ID = image_ID(IDdmC_name);
	FPAOconf[loop].dmxsize = data.image[FPAOconf[loop].dmC_ID].md[0].size[0];
	FPAOconf[loop].dmysize = data.image[FPAOconf[loop].dmC_ID].md[0].size[1];
	
	// Focal plane image stream
	
	FPAOconf[loop].FPim_ID = image_ID(IDfpim_name);
	FPAOconf[loop].FPim_dark_ID = image_ID(IDfpim_dark_name);
	FPAOconf[loop].fpimxsize = data.image[FPAOconf[loop].FPim_ID].md[0].size[0];
	FPAOconf[loop].fpimysize = data.image[FPAOconf[loop].FPim_ID].md[0].size[1];



	// Calibration
	
	FPAOconf[loop].IDrespMat_act2amp = -1;
	FPAOconf[loop].IDrespMat_act2pha = -1;

	FPAOconf_init = 1;


    return(0);
}





