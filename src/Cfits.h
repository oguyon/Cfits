#ifndef _CFITS_H
#define _CFITS_H


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>

#include <time.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>	// for random numbers


#define PI 3.14159265358979323846264338328

// Size of array CFITVARRAY
#define SZ_CFITSVARRAY 1000

// important directories and info
pid_t CfitsPID;
char CfitsDocDir[200];		// location of Cfits documentation
char CfitsSrcDir[200];		// location of Cfits source
char CfitsBuildFile[200];	// file name for Cfits source
char CfitsBuildDate[200];
char CfitsBuildTime[200];

int C_ERRNO;			// C errno (from errno.h)

/* #define DEBUG */
#define CFITSEXIT  printf("Cfits abnormally terminated, File \"%s\", line %d\n",__FILE__,__LINE__);exit(0)

#ifdef DEBUG
#define nmalloc(f,type,n) f = (type*) malloc(sizeof(type)*n);if(f==NULL){printf("ERROR: pointer \"" #f "\" allocation failed\n");exit(0);}else{printf("\nMALLOC: \""#f "\" allocated\n");}
#define nfree(f) free(f);printf("\nMALLOC: \""#f"\" freed\n");
#else
#define nmalloc(f,type,n) f = (type*) malloc(sizeof(type)*n);if(f==NULL){printf("ERROR: pointer \"" #f "\" allocation failed\n");exit(0);}
#define nfree(f) free(f);
#endif

#define TEST_ALLOC(f) if(f==NULL){printf("ERROR: pointer \"" #f "\" allocation failed\n");exit(0);}


#define SHAREDMEMDIR "/tmp"



#define NB_ARG_MAX                 20





/*^-----------------------------------------------------------------------------
| commands available through the CLI
+-----------------------------------------------------------------------------*/



typedef struct {
  char key[100];        // command keyword                 
  char module[50];      // module name
  int (* fp) ();         // command function pointer        
  char info   [1000];     // short description/help         
  char syntax [1000];   // command syntax
  char example[1000];  // command example
  char Ccall[1000];
} CMD;

typedef struct {
  char name[50];   // module name
  char info[1000]; // short description
} MODULE;




/* ---------------------------------------------------------- */
/*                                                            */
/*                                                            */
/*       COMMAND LINE ARGs / TOKENS                           */
/*                                                            */
/*                                                            */
/* ---------------------------------------------------------- */


// The command line is parsed and 

// cmdargtoken type
// 0 : unsolved
// 1 : floating point (double precision)
// 2 : long 
// 3 : string 
// 4 : existing image
// 5 : command
typedef struct
{
  int type; 
  union
  {
    double numf;
    long numl;
    char string[200];
  } val;
} CMDARGTOKEN;



int CLI_checkarg(int argnum, int argtype);







typedef struct
{
  float re;
  float im;
} complex_float;

typedef struct
{
  double re;
  double im;
} complex_double;



#define CHAR 1
#define INT 2
#define FLOAT 3
#define DOUBLE 4
#define COMPLEX_FLOAT 5
#define COMPLEX_DOUBLE 6
#define USHORT 7

int TYPESIZE[8];


//typedef float PRECISION;
//typedef complex_float CPRECISION;
//#define arrayD array.F
//#define arrayCD array.CF
#define Dtype 3
#define CDtype 5



typedef struct
{
  char name[16];
  char type; // N: unused, L: long, D: double, S: 16-char string 
  union {
    long numl;
    double numf;
    char valstr[16];
  } value;
  char comment[80];
} IMAGE_KEYWORD;



typedef struct
{
  char name[80];               // image name

  long naxis;                   // number of axis
  long size[3];                 // image size 
  long nelement;		// number of elements in image
  int atype;			// data type code   

  double creation_time;	        // creation time (since program start)
  double last_access;		// last time the image was accessed  (since program start)
  struct timespec wtime;

  int shared; // 1 if in shared memory

  int write;                 // 1 if image is being written  
  int status;
  long cnt0;                 // counter (if image is updated)
  long cnt1;
  
  long NBkw; // number of keywords

} IMAGE_METADATA;



typedef struct			/* structure used to store data arrays */
{
  int used;
  int shmfd; // if shared memory, file descriptor
  size_t memsize; // total size in memory if shared 

  IMAGE_METADATA *md;

  union
  {
    char *C;
    int *I;
    float *F;
    double *D;
    complex_float *CF;
    complex_double *CD;
    unsigned short *U;
  } array;                      // pointer to data array
 
  IMAGE_KEYWORD *kw;

} IMAGE;






typedef struct
{
  int used;
  char name[20];
  double value;
} VARIABLE;

typedef struct
{
  int used;
  char name[20];
  long value;
} VARIABLELONG;







// THIS IS WHERE EVERYTHING THAT NEEDS TO BE WIDELY ACCESSIBLE GETS STORED
typedef struct
{
  int Debug;
  int quiet;
  int overwrite;		// automatically overwrite FITS files
  double INVRANDMAX;
  gsl_rng *rndgen;		// random number generator  
  int precision;		// default precision: 0 for float, 1 for double

  // Command Line Interface (CLI) INPUT
  long NBcmd;
  long NB_MAX_COMMAND;
  CMD *cmd;
  int parseerror; // 1 if error, 0 otherwise
  long cmdNBarg;  // number of arguments in last command line
  CMDARGTOKEN cmdargtoken[NB_ARG_MAX];
  long cmdindex; // when command is found in command line, holds index of command
  long calctmp_imindex; // used to create temporary images
  int CMDexecuted; // 0 if command has not been executed, 1 otherwise
  long NBmodule;
  long NB_MAX_MODULE;
  MODULE *module;


  // shared memory default
  int SHARED_DFT;

  // Number of keyword per iamge default
  int NBKEWORD_DFT;

  // images, variables
  long NB_MAX_IMAGE;
  IMAGE *image;

  long NB_MAX_VARIABLE;
  VARIABLE *variable;

  float FLOATARRAY[1000];	// array to store temporary variables
  double DOUBLEARRAY[1000];    // for convenience
} DATA;




#define MAX_NB_FRAMES 500
#define MAX_NB_FRAMENAME_CHAR 500
#define MAX_NB_EXCLUSIONS 40






//extern int ECHO;
//extern PRECISION CFITSVARRAY[SZ_CFITSVARRAY];
//extern long CFITSVARRAY_LONG[SZ_CFITSVARRAY];

#endif
