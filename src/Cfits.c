/*^****************************************************************************
* FILE: Cfits.c : Cfits main
*
*
*
*
*
******************************************************************************/
#include <string.h>
#include <Cfits.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>

#include <fftw3.h>
#include <readline/readline.h>  
#include <readline/history.h>  


# ifdef _OPENMP
# include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
# endif



#include <gsl/gsl_rng.h> // for random numbers
#include <fitsio.h>

#include "00CORE/00CORE.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_fft/COREMOD_fft.h"


/*-----------------------------------------
*       Globals
*/
int Debug   = 0;
int Journal = 0;
int Verbose = 0;
int LIST_IMFILES = 0;
DATA data;
PIAAconf piaaconf;
ZERNIKE Zernike;
int NORMAL_EXIT = 0;
PRECISION CFITSVARRAY[SZ_CFITSVARRAY];
long CFITSVARRAY_LONG[SZ_CFITSVARRAY];
int ECHO;

/*-----------------------------------------
*       Forward References
*/
int user_function();
void fnExit1 (void);
void cfits_init();

int re_alloc();
int command_line( int argc, char **argv);



/*^-----------------------------------------------------------------------------
| 
| 
| 
|
|
+-----------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  // char cmd1[1000];
  // char cmd[1000];
  //  int executed;
  long i;
  int quiet=0;
  long tmplong;
  //  FILE *fp;
  const gsl_rng_type * rndgenType;
  PRECISION v1;
  char prompt[200];
  char promptname[100];
  char *line;
  int terminate = 0;

  TYPESIZE[0] = 0;
  TYPESIZE[1] = sizeof(char);
  TYPESIZE[2] = sizeof(int);
  TYPESIZE[3] = sizeof(float);
  TYPESIZE[4] = sizeof(double);
  TYPESIZE[5] = 2*sizeof(float);
  TYPESIZE[6] = 2*sizeof(double);
  
  atexit(fnExit1);

  data.overwrite = 0;
  data.precision = 0; // float is default precision

  // Get command-line options
  command_line( argc, argv );
  
  // 
  if( Verbose ) {
    fprintf(stdout, "%s: compiled %s %s\n",__FILE__,__DATE__,__TIME__);
    // printf("TIME = %ld\n",time(NULL));
  }

  CfitsPID = getpid();

  sprintf(promptname, "%s", argv[0]);
  sprintf(prompt,"%c[%d;%dm%s >%c[%dm ",0x1B, 1, 36, promptname, 0x1B, 0);

  printf("type \"help\" for instructions\n");
    
# ifdef _OPENMP
  printf("Running with openMP, max threads = %d  (defined by environment variable OMP_NUM_THREADS)\n",omp_get_max_threads());
# endif
  
# ifdef FFTWMT
  printf("Multi-threaded fft enabled, max threads = %d\n",omp_get_max_threads());
  fftwf_init_threads();
  fftwf_plan_with_nthreads(omp_get_max_threads());
# endif


  //    sprintf(CfitsDocDir,"%s",DOCDIR);
  //   sprintf(CfitsSrcDir,"%s",SOURCEDIR);
  //  sprintf(CfitsBuildFile,"%s",__FILE__);
  //  sprintf(CfitsBuildDate,"%s",__DATE__);
  //  sprintf(CfitsBuildTime,"%s",__TIME__);

    // Initialize random-number generator
    //
    //rndgenType = gsl_rng_ranlxs2; // best algorithm but slow
    //rndgenType = gsl_rng_ranlxs0; // not quite as good, slower
    rndgenType = gsl_rng_rand; // not as good but ~10x faster fast
    data.rndgen = gsl_rng_alloc (rndgenType);
    gsl_rng_set (data.rndgen,time(NULL));

    // warm up
    for(i=0;i<10;i++)
      v1 = gsl_rng_uniform (data.rndgen);
    
    // 
    Zernike.init = 0;
    Zernike.ZERMAX = 5000;
  
    // FFTW init
       // load fftw wisdom
    import_wisdom();

    if(sizeof(PRECISION)==sizeof(float))
      fftwf_set_timelimit(1000.0);
    else
      fftw_set_timelimit(1000.0);
    

    /*--------------------------------------------------          
    |  Check command-line arguements  
    +-------------------------------------------------*/
   

    for(i=1;i<argc;i++)
        if(strcmp(argv[i],"-quiet")==0)
            quiet=1;

    ECHO = 0;
    for(i=1;i<argc;i++)
    if(strcmp(argv[i],"-echo")==0)
        ECHO = 1;


    /* Initialize data control block
    */
    cfits_init();

    tmplong = data.NB_MAX_VARIABLE;
    C_ERRNO = 0; // initialize C error variable to 0 (no error)
    terminate = 0;
    while(terminate==0) {

      if( fopen( "STOPCFITS","r" ) != NULL ) {
	fprintf(stdout, "STOPCFITS FILE FOUND. Exiting...\n");
	exit(1);
        }
      
      if(LIST_IMFILES==1) {
	list_image_ID_file("imlist.txt");
      }

      /* Keep the number of image addresses available   
       *  NB_IMAGES_BUFFER above the number of used images 
       *
       *  Keep the number of variables addresses available 
       *  NB_VARIABLES_BUFFER above the number of used variables 
       */
      if( re_alloc() != 0 )
	{
	  fprintf(stderr,"%c[%d;%dm ERROR [ FILE: %s   FUNCTION: %s   LINE: %d ]  %c[%d;m\n", (char) 27, 1, 31, __FILE__, __func__, __LINE__, (char) 27, 0);
	  fprintf(stderr,"%c[%d;%dm Cfits memory re-allocation failed  %c[%d;m\n", (char) 27, 1, 31, (char) 27, 0);
	  exit(0);
	}
      
      compute_image_memory(data);
      compute_nb_image(data);
      
      // -------------------------------------------------------------
      //                 get user input
      // -------------------------------------------------------------
      line = readline (prompt);

          if(!strncmp(line,"help",4))
	{
	  printf("\n");
	  printf("----------- COMMANDS ---------\n");
	  printf("\n");
	  printf("\n");
	  printf("   exit                       : exit\n");
	  printf("\n");	      
	}
      else if ((!strcmp(line, "quit"))||(!strcmp(line, "exit")))
	{
	  terminate = 1;
	  printf("Closing PID %ld (prompt process)\n", (long) getpid());
	}
      add_history(line);		  
    }
   
    
    // Free 
    free(data.image);
    free(data.variable);
    gsl_rng_free (data.rndgen);

    # ifdef FFTWMT
    fftwf_cleanup_threads();
    # endif

    # ifndef FFTWMT
    fftwf_cleanup();
    # endif


    NORMAL_EXIT = 1;
    
    if(LIST_IMFILES==1) {
      if(system("rm imlist.txt")==-1)
	{
	  printERROR(__FILE__,__func__,__LINE__,"system() error");
	  exit(0);
	}
    }
    

    return(0);
}


/*^-----------------------------------------------------------------------------
| void cfits_init : Initialization the "data" structure  
| 
| 
|
|
+-----------------------------------------------------------------------------*/
void cfits_init()
{

  long tmplong;
  int i;
  

  
  /* initialization of the data structure 
   */
  data.quiet           = 1;
  data.NB_MAX_IMAGE    = 100;
  data.NB_MAX_VARIABLE = 100;
  data.INVRANDMAX      = 1.0/RAND_MAX;
  
  // Allocate data.image
  data.image           = (IMAGE *) calloc(data.NB_MAX_IMAGE, sizeof(IMAGE));
  if(data.image==NULL)  {
    printERROR(__FILE__,__func__,__LINE__,"Allocation of data.image has failed - exiting program");
    exit(0);
  }
  
  
  // Allocate data.variable
  data.variable = (VARIABLE *) calloc(data.NB_MAX_VARIABLE, sizeof(VARIABLE));
  if(data.variable==NULL)  {
    printERROR(__FILE__,__func__,__LINE__,"Allocation of data.variable has failed - exiting program");       
    exit(0);
  }
  
  data.image[0].used   = 0;
  tmplong              = data.NB_MAX_VARIABLE;
  data.NB_MAX_VARIABLE = data.NB_MAX_VARIABLE + NB_VARIABLES_BUFFER_REALLOC ;
  
  // 
  data.variable = (VARIABLE *) realloc(data.variable, data.NB_MAX_VARIABLE*sizeof(VARIABLE));
  for(i=tmplong;i<data.NB_MAX_VARIABLE;i++)
    data.variable[i].used = 0;
  
  
  
  tmplong = data.NB_MAX_VARIABLE;
  if (data.variable == NULL)   {
    printERROR(__FILE__,__func__,__LINE__,"Reallocation of data.variable has failed - exiting program");
    exit(0);
  }
  
 
    
  create_variable_ID("_PI",3.14159265358979323846264338328);
  create_variable_ID("_e",exp(1));
  create_variable_ID("_gamma",0.5772156649);
  create_variable_ID("_c",299792458.0);
  create_variable_ID("_h",6.626075540e-34);
  create_variable_ID("_k",1.38065812e-23);
  create_variable_ID("_pc",3.0856776e16);
  create_variable_ID("_ly",9.460730472e15);
  create_variable_ID("_AU",1.4959787066e11);

  srand(time(NULL));




}

/*^-----------------------------------------------------------------------------
| 
| 
| 
|
|
+-----------------------------------------------------------------------------*/
int user_function()
{
    printf("-");
    fflush(stdout);
    printf("-");
    fflush(stdout);
  
    return(0);
}

/*^-----------------------------------------------------------------------------
| 
| 
| 
|
|
+-----------------------------------------------------------------------------*/
void fnExit1 (void)
{
  //  
}



/*^-----------------------------------------------------------------------------
| 
|  re_alloc    : keep the number of images addresses available   
| 		 NB_IMAGES_BUFFER above the number of used images 
|
|                keep the number of variables addresses available 
|                NB_VARIABLES_BUFFER above the number of used variables 
|
| NOTE:  this should probably be renamed and put in the module/memory/memory.c 
|
+-----------------------------------------------------------------------------*/
int re_alloc()
{
  int i;
  long  tmplong;


  /* keeps the number of images addresses available   
   *  NB_IMAGES_BUFFER above the number of used images 
   */
  if((compute_nb_image(data)+NB_IMAGES_BUFFER)>data.NB_MAX_IMAGE)
    {
      tmplong = data.NB_MAX_IMAGE;
      data.NB_MAX_IMAGE = data.NB_MAX_IMAGE + NB_IMAGES_BUFFER_REALLOC;
      data.image = (IMAGE *) realloc(data.image, data.NB_MAX_IMAGE*sizeof(IMAGE));
      if(data.image==NULL)   {
	printERROR(__FILE__,__func__,__LINE__,"Reallocation of data.image has failed - exiting program");
	return -1;      //  exit(0);
      }
      
      for(i=tmplong;i<data.NB_MAX_IMAGE;i++)   {
	data.image[i].used = 0;
      }
    }
  
  /* keeps the number of variables addresses available 
   *  NB_VARIABLES_BUFFER above the number of used variables 
   */
  if((compute_nb_variable(data)+NB_VARIABLES_BUFFER)>data.NB_MAX_VARIABLE)
    {
      data.NB_MAX_VARIABLE = data.NB_MAX_VARIABLE + NB_VARIABLES_BUFFER_REALLOC;
      data.variable = (VARIABLE *) realloc(data.variable, data.NB_MAX_VARIABLE*sizeof(VARIABLE));
      if (data.variable==NULL)
	{ 
	  printERROR(__FILE__,__func__,__LINE__,"Reallocation of data.variable has failed - exiting program");
	  return -1;   // exit(0);
	}
    }
  
  return 0;
}















/*^-----------------------------------------------------------------------------
| static PF 
| command_line  : parse unix command line options. 
|
|   int argc    :
|   char **argv :  
|
|   TO DO : allow option values. eg: debug=3 
+-----------------------------------------------------------------------------*/
int command_line( int argc, char **argv)
{
  FILE *fp;
  int rc;
  int i;
  char startup_info[1024];
  float f1;
  struct tm *ptr;
  time_t tm;
  char command[1000];
  
    NORMAL_EXIT=1;

    for (i=1; i<argc; i++) {



        //  set debug

        if (!strcasecmp(argv[i],"-debug") ) {
            Debug = 2;
            strcat( startup_info, "        DEBUG=2 Option: '-debug'\n"); 
        }
    }


    for (i=1; i<argc; i++) {

      if ( Verbose )
	fprintf(stdout, "argv[%d]:[%s]\n", i, argv[i]);

        /* We handled the -debug op above. 
        *  but this prevents logic from complaining if we see it
        *  it again.
        */
        if (!strcasecmp(argv[i],"-debug") ); 
        else if ( !strcasecmp(argv[i],"-i") ) {
	  fprintf(stdout, "\n");
	  fprintf(stdout, "--------------- GENERAL ----------------------\n");
	  fprintf(stdout, "%s VERSION   %s\n",  argv[0], PACKAGE_VERSION    ); 
	  fprintf(stdout, "%s BUILT   %s %s\n", __FILE__,__DATE__,__TIME__);
	  fprintf(stdout, "\n");
	  fprintf(stdout, "--------------- SETTINGS ---------------------\n");
	  if(data.precision==0)
	    fprintf(stdout, "Default precision upon startup : float\n");
	  if(data.precision==1)
	    fprintf(stdout, "Default precision upon startup : double\n");
	  fprintf(stdout, "\n");
	  fprintf(stdout, "--------------- LIBRARIES --------------------\n");
	  fprintf(stdout, "READLINE : version %x\n",RL_READLINE_VERSION);
# ifdef _OPENMP
	  fprintf(stdout, "OPENMP   : Compiled by an OpenMP-compliant implementation.\n");
# endif
	  fprintf(stdout, "CFITSIO  : version %f\n",fits_get_version(&f1));
	  fprintf(stdout, "\n");
    	  fprintf(stdout, "--------------- DIRECTORIES ------------------\n");

	  //	  fprintf(stdout, "Config directory:        %s\n", CONFIGDIR);
	  // fprintf(stdout, "Object directory:        %s\n", OBJDIR);
	  //fprintf(stdout, "Source directory:        %s\n", SOURCEDIR);
	  //fprintf(stdout, "Documentation directory: %s\n", DOCDIR);
	  fprintf(stdout, "\n");
	  exit(0);
        }
        else if ( !strcasecmp(argv[i],"-verbose") ) {
	    Verbose = 1; 
        }
	else if ( !strcasecmp(argv[i],"-o") ) {
	  data.overwrite = 1;
	  printf("OVERWRITE = 1: USE WITH CAUTION - WILL OVERWRITE EXISTING FITS FILES\n"); 
	}
	else if ( !strcasecmp(argv[i],"-j") ) {
	  Journal = 1;
	  tm = time(NULL);
	  ptr = localtime(&tm);
	  fp = fopen("cfits_cmdlog.txt","a");
	  fprintf(fp,"# NEW SESSION ------------ %s\n",asctime(ptr));
	  fclose(fp);
	}
        else if ( !strcasecmp(argv[i],"-quiet") ) {
            data.quiet           = 1;
  	    Debug   = 0;
	    Verbose = 0; 
        }
        else if ( !strcasecmp(argv[i],"-h") ) {

	  //	  sprintf(command,"more %s/help.txt",DOCDIR);
	  if(system(command)==-1)
	    {
	      printERROR(__FILE__,__func__,__LINE__,"system() error");
	      exit(0);
	    }
	  exit(0);
        }	
	else if  ( !strcasecmp(argv[i],"-l") ) {
	  printf("List of images on file imlist.txt\n");
	  if(system("rm imlist.txt")==-1) // just in case last Cfits run crashed and left a imlist.txt file
	    {
	      printERROR(__FILE__,__func__,__LINE__,"system() error");
	      exit(0);
	    }
	  LIST_IMFILES = 1;
	}
	else {
          fprintf(stderr, "%s : Unrecognized option: %s\n", argv[0], argv[i]);
          exit(0);
        }
    }

    NORMAL_EXIT=0;
    return 0;

}


/*^-----------------------------------------------------------------------------
| 
| 
| 
|
|
+-----------------------------------------------------------------------------*/
