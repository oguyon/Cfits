#ifndef _CFITS_H
#define _CFITS_H


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>

#include <time.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h> // for random numbers


#define PI 3.14159265358979323846264338328

// Size of array CFITVARRAY
#define SZ_CFITSVARRAY 1000

// important directories and info
pid_t CfitsPID;
char CfitsDocDir[200]; // location of Cfits documentation
char CfitsSrcDir[200]; // location of Cfits source
char CfitsBuildFile[200]; // file name for Cfits source
char CfitsBuildDate[200];
char CfitsBuildTime[200];

int C_ERRNO; // C errno (from errno.h)

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

/* maximum number of characters in the array name */
#define max_nb_charact 450

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

int TYPESIZE[7];


typedef float PRECISION;
typedef complex_float CPRECISION;
#define arrayD array.F
#define arrayCD array.CF
#define Dtype 3
#define CDtype 5


/*
typedef double PRECISION;
typedef complex_double CPRECISION;
#define arrayD array.D
#define arrayCD array.CD
#define Dtype 4
#define CDtype 6
*/


typedef struct   /* structure used to store data arrays */
{
  int used; /* 0 if unused, 1 if used */
  char name[100];
  long naxis;
  long size[3];
  long nelement; // number of elements in image
  int atype; // data type code
  union {
    char *C;
    int *I;
    float *F;
    double *D;
    complex_float *CF;
    complex_double *CD;
  } array;
  time_t creation_time; /* creation time */
  time_t last_access;  /* last time the image was accessed */
} IMAGE;


typedef struct
{
  int used;
  char name[20];
  double value;
} VARIABLE;

typedef struct
{
  int quiet;
  int overwrite; // automatically overwrite FITS files
  long NB_MAX_IMAGE;
  int precision; // default precision: 0 for float, 1 for double
  IMAGE *image;
  long NB_MAX_VARIABLE;
  VARIABLE *variable;
  PRECISION INVRANDMAX;
  gsl_rng *rndgen; // random number generator  
  PRECISION CFITSVARRAY[1000]; // array to store temporary variables
} DATA;

#define MAX_NB_FRAMES 500
#define MAX_NB_FRAMENAME_CHAR 500
#define MAX_NB_EXCLUSIONS 40

typedef struct /* structure to store Zernike coefficients */
{
  int init;
  long ZERMAX;
  long *Zer_n;
  long *Zer_m;
  double *R_array;
} ZERNIKE;

typedef struct /* structure used to store the log file in img_reduce */
{
  int OBS_YEAR;
  int OBS_MONTH;
  int OBS_DAY;
  PRECISION SITE_LONG; /* in degrees */
  PRECISION SITE_LAT; /* in degrees */
  PRECISION target_RA; /* in deg */
  PRECISION target_DEC; /* in deg */
  PRECISION detector_saturation; /* in ADU */
  PRECISION pixel_scale; /* in arcsec per pixel */
  PRECISION detector_readout_noise; /* in e- */
  PRECISION detector_gain; /* in e-/ADU */
  PRECISION magnitude_zeropt; /* magn of 1 ADU/s */
  char FLAT_file_name[MAX_NB_FRAMENAME_CHAR];
  char BAD_PIX_file_name[MAX_NB_FRAMENAME_CHAR];
  int SKY;
  char SKY_file_name[MAX_NB_FRAMENAME_CHAR];
  int TIMELOG;
  char TIMELOG_file_name[MAX_NB_FRAMENAME_CHAR];
  int nb_exclusion_files;
  char EXCLUSION_file_name[MAX_NB_EXCLUSIONS][MAX_NB_FRAMENAME_CHAR];
  int nbFRAMES;
  char FRAME_file_name[MAX_NB_FRAMENAME_CHAR][MAX_NB_FRAMES];
  PRECISION FRAME_etime[MAX_NB_FRAMES];
  PRECISION FRAME_UT[MAX_NB_FRAMES];
  PRECISION FRAME_FWHM[MAX_NB_FRAMES];
  PRECISION FRAME_cx[MAX_NB_FRAMES];
  PRECISION FRAME_cy[MAX_NB_FRAMES];
  PRECISION FRAME_orientation[MAX_NB_FRAMES];
  int FRAME_field_rotator[MAX_NB_FRAMES];
  int FRAME_isobject[MAX_NB_FRAMES];
  int FRAME_ispsf[MAX_NB_FRAMES];
  int FRAME_issky[MAX_NB_FRAMES];
  int FRAME_exclusion[MAX_NB_FRAMES];
  int FRAME_exclusion_centering_mode[MAX_NB_FRAMES]; /* 0 for absolute, 1 for relative to PSF center */
  PRECISION FRAME_exclusion_x[MAX_NB_FRAMES];
  PRECISION FRAME_exclusion_y[MAX_NB_FRAMES];
  PRECISION FRAME_exclusion_r[MAX_NB_FRAMES];
  int FRAME_setzeroh[MAX_NB_FRAMES];
  int FRAME_setzerov[MAX_NB_FRAMES];
  int FRAME_hvexclusion[MAX_NB_FRAMES];

} LOGFILE;

typedef struct
{
  int direction; /* 0 for apodization, 1 for de-apodization */
  long PUPFILESIZE; /* number of points in the entrance pupil radial profile file */
  long NPUPFILESIZE; /* number of points in the PIAAsimul internal pupil radial profile */
  double *pup_amp_profile; /* internal pupil radial profile */
  double epsilon;  /* below this value, the exit pupil is considered "empty" */

  long size;
  PRECISION lambda;
  double scale;

  double beam_radius;
  double extrapolate_dist;
  double sphere_radius;
  double parabola_distance;
  double focal;
  int FOCAL;
  double screendist;
  double offaxisdist;

  double M0_pos_err_z; /* z position error of M0 in m */
  double M0_pos_err_x; /* x position error of M0 in m */
  double M0_pos_err_y; /* y position error of M0 in m */
  double M0_pos_err_xtilt; /* xtilt position error of M0 in rad */
  double M0_pos_err_ytilt; /* ytilt position error of M0 in rad */
  double M0_pos_PA; /* PA of the mirror in rad */
  double M1_pos_err_z; /* z position error of M1 in m */
  double M1_pos_err_x; /* x position error of M1 in m */
  double M1_pos_err_y; /* x position error of M1 in m */
  double M1_pos_err_xtilt; /* xtilt position error of M1 in rad */
  double M1_pos_err_ytilt; /* ytilt position error of M1 in rad */
  double M1_pos_PA; /* PA of the mirror in rad */

  double *r0;
  double *r1;
  double *r1fr0; /* r1 as a function of r0 */
  double *r0fr1; /* r0 as a function of r1 */
  long NBr1fr0_pts; /* number of points for the 2 above functions */
  long NBpoints;
  double *M0_z; /* in m */
  double *M1_z; /* in m */
  long NBradsamppts; /* regular sampling */
  double *M0r_z; /*in m, regular sampling */
  double *M1r_z; /*in m, regular sampling */

  long M0radfit_order;
  long M1radfit_order;
  double *M0radfit_comp;
  double *M1radfit_comp;

  PRECISION Lyot_pixel_scale; /* m per pixel */
  PRECISION Lyot_beam_radius; /* in m */
  long Lyot_pupil_size; /* in pixels (image size) */
  PRECISION Lyot_fmask_radius; /* in pixels */

  PRECISION oneDpupcoeff;
  PRECISION oneDpsfcoeff;
  PRECISION oneDPSFSIZE;

  int amp_noremap; /* no remapping for amplitude */

  long ZERMAX; /* max number of Zernikes in the off-axis residual fit */
  int M0_ADD_ZERN; /* M0 Zernikes included */
  int M1_ADD_ZERN; /* M1 Zernikes included */
  PRECISION *M0_ZER; /* M0 Zernikes added */
  PRECISION *M1_ZER; /* M1 Zernikes added */

  int M0_ZERSHAPE; /* 1 if M0 is entirely defined by Zernikes */
  int M1_ZERSHAPE; /* 1 if M1 is entirely defined by Zernikes */
  PRECISION *M0ZER;
  PRECISION *M1ZER; 
  
  int M0_ADD_SURF; /* M0 additionnal surface included - image name M0adds */
  int M1_ADD_SURF; /* M1 additionnal surface included - image name M1adds */
  long M0_ADD_SURF_ID;
  long M1_ADD_SURF_ID;
  long M_ADD_SURF_SIZE; /* size of above 2 images */
  PRECISION M_ADD_SURF_RADIUS; /* radius of pupil in the above 2 images */

  /* raytracing memory structures */
  long NB1ray;
  double *i_in; /* intensity of the light ray */
  double *x_in;
  double *y_in;
  double *x_out1;
  double *y_out1;
  double *x_out2;
  double *y_out2;
  double *OPD_out;
  int *ray_OK;

  int alignerror;

  int init_cosapo1D; /* 1 if intialized */
  int cosapo1D; /* 1 if a 1D cosine apodization is multiplied to the default apodization */
  double cosapo1D_period; /* period as a fraction of pupil size */
  double cosapo1D_coeff;
  double *x1fx0;
  double *x0fx1;
  double *xfactor;
  long apo1D_NBpts;

} PIAAconf;

extern int ECHO;
extern PRECISION CFITSVARRAY[SZ_CFITSVARRAY];
extern long CFITSVARRAY_LONG[SZ_CFITSVARRAY];

#endif
