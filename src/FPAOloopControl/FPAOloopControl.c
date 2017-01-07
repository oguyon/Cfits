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



long NB_FPAOloopcontrol = 1;
long FPLOOPNUMBER = 0; // current loop index
int FPAOloopcontrol_meminit = 0;
int FPAOlooploadconf_init = 0;

#define FPAOconfname "/tmp/FPAOconf.shm"
FPAOLOOPCONTROL_CONF *FPAOconf; // configuration - this can be an array










long FPaoconfID_wfsim = -1;
int FPWFSatype;
long FPaoconfID_wfsdark = -1;

long FPaoconfID_dmC = -1;
long FPaoconfID_dmRM = -1;




int FPAO_loadcreateshm_log = 0; // 1 if results should be logged in ASCII file
FILE *FPAO_loadcreateshm_fplog;








// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string, not existing image
// 4: existing image
// 5: string 





int FPAOloopControl_loadconfigure_cli()
{
  if(CLI_checkarg(1,2)==0)
    {
      FPAOloopControl_loadconfigure(data.cmdargtoken[1].val.numl, 1, 10);
      return 0;
    }
  else
    return 1;
}




int init_FPAOloopControl()
{


    strcpy(data.cmd[data.NBcmd].key,"FPaolloadconf");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = FPAOloopControl_loadconfigure_cli;
    strcpy(data.cmd[data.NBcmd].info,"load FPAO loop configuration");
    strcpy(data.cmd[data.NBcmd].syntax,"<loop #>");
    strcpy(data.cmd[data.NBcmd].example,"FPaolloadconf 1");
    strcpy(data.cmd[data.NBcmd].Ccall,"int FPAOloopControl_loadconfigure(long loopnb, 1, 10)");
    data.NBcmd++;


    strcpy(data.module[data.NBmodule].name, __FILE__);
    strcpy(data.module[data.NBmodule].info, "FP AO loop control");
    data.NBmodule++;


    return 0;
}







/*** mode = 0 or 1. if mode == 1, simply connect */


long FPAOloopControl_InitializeMemory(int mode)
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
	int loop;
	
	
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
    
   
   
   for(loop=0;loop<NB_FPAOloopcontrol; loop++)
	{
	// DM streams
	FPAOconf[loop].dmxsize = 0;
	FPAOconf[loop].dmysize = 0;
	
	// Focal plane image stream
	FPAOconf[loop].sizexWFS = 0;
	FPAOconf[loop].sizeyWFS = 0;

	}

	// Calibration
	FPAOloopcontrol_meminit = 1;


    return(0);
}






int FPAOloopControl_loadconfigure(long loop, int mode, int level)
{
   FILE *fp;
    char content[200];
    char name[200];
    char name1[200];
    char fname[200];
    long ID;
    long *sizearray;
    int vOK;
    int kw;
    long k;
    int r;
    int sizeOK;
    char command[500];
    int CreateSMim;
    long ID1tmp, ID2tmp;
    long ii;
    long kk, tmpl;
    char testdirname[200];
    
    FILE *fplog; // human-readable log of load sequence

    if((fplog=fopen("loadconf.log","w"))==NULL)
    {
        printf("ERROR: file loadconf.log missing\n");
        exit(0);
    }
    FPAO_loadcreateshm_log = 1;
    FPAO_loadcreateshm_fplog = fplog;


    if(FPAOloopcontrol_meminit==0)
        FPAOloopControl_InitializeMemory(0);



    // printf("mode = %d\n", mode); // not used yet


    // Name definitions for shared memory

    sprintf(name, "FPaol%ld_dmC", loop);
    printf("FP loop DM control file name : %s\n", name);
    strcpy(FPAOconf[loop].dmCname, name);

    sprintf(name, "FPaol%ld_dmRM", loop);
    printf("FP loop DM RM file name : %s\n", name);
    strcpy(FPAOconf[loop].dmRMname, name);

    sprintf(name, "FPaol%ld_wfsim", loop);
    printf("FP loop WFS file name: %s\n", name);
    strcpy(FPAOconf[loop].WFSname, name);



    sizearray = (long*) malloc(sizeof(long)*3);


    // READ LOOP NAME

    if((fp=fopen("./conf/conf_LOOPNAME.txt","r"))==NULL)
    {
        printf("ERROR: file ./conf/conf_LOOPNAME.txt missing\n");
        exit(0);
    }
    r = fscanf(fp, "%s", content);
    printf("loop name : %s\n", content);
    fprintf(fplog, "FPAOconf[%ld].name = %s\n", loop, FPAOconf[loop].name);
    fclose(fp);
    fflush(stdout);
    strcpy(FPAOconf[loop].name, content);

    
   
   
	
	if((fp=fopen("./conf/conf_hardwlatency.txt", "r"))==NULL)
    {
        printf("WARNING: file ./conf/conf_hardwlatency.txt missing\n");
    }
    else
    {
        r = fscanf(fp, "%f", &FPAOconf[loop].hardwlatency);
        printf("hardwlatency : %f\n", FPAOconf[loop].hardwlatency);
        fclose(fp);
        fflush(stdout);
        fprintf(fplog, "AOconf[%ld].hardwlatency = %f\n", loop, FPAOconf[loop].hardwlatency);
   }
	
	FPAOconf[loop].hardwlatency_frame = FPAOconf[loop].hardwlatency * FPAOconf[loop].loopfrequ;
	






	if((fp=fopen("./conf/conf_loopfrequ.txt","r"))==NULL)
    {
        printf("WARNING: file ./conf/conf_loopfrequ.txt missing\n");
        printf("Using default loop speed\n");
        fprintf(fplog, "WARNING: file ./conf/conf_loopfrequ.txt missing. Using default loop speed\n");
        FPAOconf[loop].loopfrequ = 2000.0;
    }
    else
    {
        r = fscanf(fp, "%s", content);
        printf("loopfrequ : %f\n", atof(content));
        fclose(fp);
        fflush(stdout);
        FPAOconf[loop].loopfrequ = atof(content);
        fprintf(fplog, "FPAOconf[%ld].loopfrequ = %f\n", loop, FPAOconf[loop].loopfrequ);
    }




    // Connect to WFS camera
    // This is where the size of the WFS is fixed
    FPaoconfID_wfsim = read_sharedmem_image(FPAOconf[loop].WFSname);
    if(FPaoconfID_wfsim == -1)
        fprintf(fplog, "ERROR : cannot read shared memory stream %s\n", FPAOconf[loop].WFSname);
    else
        fprintf(fplog, "stream %s loaded as ID = %ld\n", FPAOconf[loop].WFSname, FPaoconfID_wfsim);


    FPAOconf[loop].sizexWFS = data.image[FPaoconfID_wfsim].md[0].size[0];
    FPAOconf[loop].sizeyWFS = data.image[FPaoconfID_wfsim].md[0].size[1];
    FPAOconf[loop].sizeWFS = FPAOconf[loop].sizexWFS*FPAOconf[loop].sizeyWFS;

    fprintf(fplog, "FPAO WFS stream size = %ld x %ld\n", FPAOconf[loop].sizexWFS, FPAOconf[loop].sizeyWFS);




    // The AOloopControl_xDloadcreate_shmim functions work as follows:
    // If file already loaded, use it (we assume it's already been properly loaded)
    // If not, attempt to read it from shared memory
    // If not available in shared memory, create it in shared memory
    // if "fname" exists, attempt to load it into the shared memory image

    sprintf(name, "FPaol%ld_wfsdark", loop);
    sprintf(fname, "./conf/FPaol%ld_wfsdark.fits", loop);
    FPaoconfID_wfsdark = AOloopControl_2Dloadcreate_shmim(name, fname, FPAOconf[loop].sizexWFS, FPAOconf[loop].sizeyWFS);
    
    


    // Connect to DM
    // Here the DM size is fixed
    //


    FPaoconfID_dmC = image_ID(FPAOconf[loop].dmCname);
    if(FPaoconfID_dmC==-1)
    {
        printf("connect to %s\n", FPAOconf[loop].dmCname);
        FPaoconfID_dmC = read_sharedmem_image(FPAOconf[loop].dmCname);
        if(FPaoconfID_dmC==-1)
        {
            printf("ERROR: cannot connect to shared memory %s\n", FPAOconf[loop].dmCname);
            exit(0);
        }
    }
    FPAOconf[loop].dmxsize = data.image[FPaoconfID_dmC].md[0].size[0];
    FPAOconf[loop].dmysize = data.image[FPaoconfID_dmC].md[0].size[1];
    FPAOconf[loop].dmsize = FPAOconf[loop].dmxsize*FPAOconf[loop].dmysize;
    
    fprintf(fplog, "Connected to DM %s, size = %ld x %ld\n", FPAOconf[loop].dmCname, FPAOconf[loop].dmxsize, FPAOconf[loop].dmysize);



    FPaoconfID_dmRM = image_ID(FPAOconf[loop].dmRMname);
    if(FPaoconfID_dmRM==-1)
    {
        printf("connect to %s\n", FPAOconf[loop].dmRMname);
        FPaoconfID_dmRM = read_sharedmem_image(FPAOconf[loop].dmRMname);
        if(FPaoconfID_dmRM==-1)
        {
            printf("ERROR: cannot connect to shared memory %s\n", FPAOconf[loop].dmRMname);
            exit(0);
        }
    }
    fprintf(fplog, "stream %s loaded as ID = %ld\n", FPAOconf[loop].dmRMname, FPaoconfID_dmRM);



    list_image_ID();

	FPAOlooploadconf_init = 1;
    
    FPAO_loadcreateshm_log = 0;
    fclose(fplog);


    return(0);

}



















long FPAOloopControl_acquireRM_level1(float ampl)
{
	// pokes X and Y patterns
	// X pattens:
	// period = 1 act
	// period = 2 act
	
	
	
	
	return 0;
}


