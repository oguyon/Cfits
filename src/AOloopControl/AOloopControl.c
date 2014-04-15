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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>




#include "Cfits.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "AOloopControl/AOloopControl.h"


#include "ZernikePolyn/ZernikePolyn.h"




extern DATA data;


#define NB_AOloopcontrol 10 // max number of loops
long LOOPNUMBER = 0; // current loop index
int AOloopcontrol_meminit = 0;
int AOlooploadconf_init = 0;

#define AOconfname "/tmp/AOconf.shm"
AOLOOPCONTROL_CONF *AOconf; // configuration - this can be an array
int *AOconf_loaded = 0;
int *AOconf_fd; 


int AOloopControl_camimage_extract2D_sharedmem_loop(char *in_name, char *out_name, long size_x, long size_y, long xstart, long ystart);
int Measure_Resp_Matrix(long loop, long NbAve, float amp, long nbloop);
int compute_ControlMatrix(long loop, long NB_MODE_REMOVED, char *ID_Rmatrix_name, char *ID_Cmatrix_name, char *ID_VTmatrix_name);
int Average_cam_frames(long loop, long NbAve);
int AOloopControl_loadconfigure(long loopnb, char *fname);
int AOloopControl_run();

int AOloopControl_setLoopNumber(long loop);
int AOloopControl_loopkill();
int AOloopControl_loopon();
int AOloopControl_loopoff();
int AOloopControl_setgain(float gain);
int AOloopControl_setmaxlimit(float maxlimit);

// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//



int AOloopControl_camimage_extract2D_sharedmem_loop_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,2)+CLI_checkarg(4,2)+CLI_checkarg(5,2)+CLI_checkarg(6,2)==0)
    {
      AOloopControl_camimage_extract2D_sharedmem_loop(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string , data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl, data.cmdargtoken[5].val.numl, data.cmdargtoken[6].val.numl);
      return 0;
    }
  else
    return 1;
}

int AOloopControl_loadconfigure_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,3)==0)
    {
      AOloopControl_loadconfigure(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.string);
      return 0;
    }
  else
    return 1;
}

int AOloopControl_setLoopNumber_cli()
{
  if(CLI_checkarg(1,2)==0)
    {
      AOloopControl_setLoopNumber(data.cmdargtoken[1].val.numl);
      return 0;
    }
  else
    return 1;
}


int AOloopControl_setgain_cli()
{
  if(CLI_checkarg(1,1)==0)
    {
      AOloopControl_setgain(data.cmdargtoken[1].val.numf);
      return 0;
    }
  else
    return 1;
}

int AOloopControl_setmaxlimit_cli()
{
  if(CLI_checkarg(1,1)==0)
    {
      AOloopControl_setmaxlimit(data.cmdargtoken[1].val.numf);
      return 0;
    }
  else
    return 1;
}

int AOloopControl_Measure_Resp_Matrix_cli()
{
  if(CLI_checkarg(1,2)+CLI_checkarg(2,1)+CLI_checkarg(3,2)==0)
    {
      Measure_Resp_Matrix(LOOPNUMBER, data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numl);
      return 0;
    }
  else
    return 1;
}




int init_AOloopControl()
{
  LOOPNUMBER = 0;
  
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "AO loop control");
  data.NBmodule++;

  strcpy(data.cmd[data.NBcmd].key,"cropshim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_camimage_extract2D_sharedmem_loop_cli;
  strcpy(data.cmd[data.NBcmd].info,"crop shared mem image");
  strcpy(data.cmd[data.NBcmd].syntax,"<input image> <output image>");
  strcpy(data.cmd[data.NBcmd].example,"cropshim imin imout");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_camimage_extract2D_sharedmem_loop(char *in_name, char *out_name, long size_x, long size_y, long xstart, long ystart)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolloadconf");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loadconfigure_cli;
  strcpy(data.cmd[data.NBcmd].info,"load AO loop configuration from file");
  strcpy(data.cmd[data.NBcmd].syntax,"<loop #> <conf file>");
  strcpy(data.cmd[data.NBcmd].example,"AOlooploadconf 1 aoloop1.conf");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loadconfigure(long loopnb, char *fname)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolacqresp");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_Measure_Resp_Matrix_cli;
  strcpy(data.cmd[data.NBcmd].info,"acquire AO response matrix and WFS reference");
  strcpy(data.cmd[data.NBcmd].syntax,"<ave#> <ampl> <nbloop>");
  strcpy(data.cmd[data.NBcmd].example,"aolacqresp 50 0.1 5");
  strcpy(data.cmd[data.NBcmd].Ccall,"int Measure_Resp_Matrix(long loop, long NbAve, float amp, long nbloop)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolrun");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_run;
  strcpy(data.cmd[data.NBcmd].info,"run AO loop");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"AOlooprun");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_run()");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolkill");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopkill;
  strcpy(data.cmd[data.NBcmd].info,"kill AO loop");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aolkill");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setLoopNumber()");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"aolnb");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setLoopNumber_cli;
  strcpy(data.cmd[data.NBcmd].info,"set AO loop #");
  strcpy(data.cmd[data.NBcmd].syntax,"<loop nb>");
  strcpy(data.cmd[data.NBcmd].example,"AOloopnb 0");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setLoopNumber(long loop)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"aolon");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopon;
  strcpy(data.cmd[data.NBcmd].info,"turn loop on");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aolon");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loopon()");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aoloff");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_loopoff;
  strcpy(data.cmd[data.NBcmd].info,"turn loop off");
  strcpy(data.cmd[data.NBcmd].syntax,"no arg");
  strcpy(data.cmd[data.NBcmd].example,"aoloff");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_loopoff()");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aolsetgain");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setgain_cli;
  strcpy(data.cmd[data.NBcmd].info,"set gain");
  strcpy(data.cmd[data.NBcmd].syntax,"<gain value>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetgain");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setgain(float gain)");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"aolsetmaxlim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AOloopControl_setmaxlimit_cli;
  strcpy(data.cmd[data.NBcmd].info,"set max limit for AO mode correction");
  strcpy(data.cmd[data.NBcmd].syntax,"<limit value>");
  strcpy(data.cmd[data.NBcmd].example,"aolsetmaxlim");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_setmaxlimit(float maxlimit)");
  data.NBcmd++;
  


  // add atexit functions here
  // atexit((void*) SCEXAO_DM_unloadconf);
  
  return 0;
}





//
// every time im_name changes (counter increments), crop it to out_name in shared memory
//
int AOloopControl_camimage_extract2D_sharedmem_loop(char *in_name, char *out_name, long size_x, long size_y, long xstart, long ystart)
{
  long iiin,jjin, iiout, jjout;
  long IDin, IDout;
  int atype;
  long *sizeout;
  long long cnt0;
  
  sizeout = (long*) malloc(sizeof(long)*2);
  sizeout[0] = size_x;
  sizeout[1] = size_y;

  IDin = image_ID(in_name);
  atype = data.image[IDin].md[0].atype;

  // Create shared memory output image 
  IDout = create_image_ID(out_name, 2, sizeout, atype, 1, 0);
  
  cnt0 = -1;

  switch (atype) {
  case USHORT :
    while(1)
      {
	usleep(100);
	if(data.image[IDin].md[0].cnt0!=cnt0)
	  {
	    cnt0 = data.image[IDin].md[0].cnt0;
	    for(iiout=0; iiout<size_x; iiout++)
	      for(jjout=0; jjout<size_y; jjout++)
		{
		  iiin = xstart + iiout;
		  jjin = ystart + jjout;
		  data.image[IDout].array.U[jjout*size_x+iiout] = data.image[IDin].array.U[jjin*data.image[IDin].md[0].size[0]+iiin];
		}
	    data.image[IDout].md[0].cnt0 = cnt0;
	  }
      }
    break;
  case FLOAT :
    while(1)
      {
	usleep(100);
	if(data.image[IDin].md[0].cnt0!=cnt0)
	  {
	    cnt0 = data.image[IDin].md[0].cnt0;
	    for(iiout=0; iiout<size_x; iiout++)
	      for(jjout=0; jjout<size_y; jjout++)
		{
		  iiin = xstart + iiout;
		  jjin = ystart + jjout;
		  data.image[IDout].array.F[jjout*size_x+iiout] = data.image[IDin].array.F[jjin*data.image[IDin].md[0].size[0]+iiin];
		}
	    data.image[IDout].md[0].cnt0 = cnt0;
	  }
      }
    break;
  default :
     printf("ERROR: DATA TYPE NOT SUPPORTED\n");
      exit(0);
    break;
  }
  free(sizeout);

  return(0);
}




//
// Computes control matrix
// Conventions:
//   m: number of actuators (= NB_MODES)
//   n: number of sensors  (= # of pixels)
int compute_ControlMatrix(long loop, long NB_MODE_REMOVED, char *ID_Rmatrix_name, char *ID_Cmatrix_name, char *ID_VTmatrix_name) /* works even for m != n */
{
  FILE *fp;
  long ii1, jj1, k, ii;
  gsl_matrix *matrix_D; /* this is the response matrix */
  gsl_matrix *matrix_Ds; /* this is the pseudo inverse of D */ 
  gsl_matrix *matrix_Dtra;
  gsl_matrix *matrix_DtraD;
  gsl_matrix *matrix_DtraDinv;
  gsl_matrix *matrix_DtraD_evec;
  gsl_matrix *matrix1;
  gsl_matrix *matrix2;
  gsl_vector *matrix_DtraD_eval;
  gsl_eigen_symmv_workspace *w;
 
  gsl_matrix *matrix_save;

  long m;
  long n;
  long ID_Rmatrix, ID_Cmatrix, ID_VTmatrix;
  long *arraysizetmp;



  arraysizetmp = (long*) malloc(sizeof(long)*3);


  ID_Rmatrix = image_ID(ID_Rmatrix_name);

  n = data.image[ID_Rmatrix].md[0].size[0]*data.image[ID_Rmatrix].md[0].size[1]; //AOconf[loop].NBDMmodes;
  m = data.image[ID_Rmatrix].md[0].size[2]; //AOconf[loop].sizeWFS;

  /* in this procedure, m=number of actuators/modes, n=number of WFS elements */
  //  long m = smao[0].NBmode;
  // long n = smao[0].NBwfselem;

  printf("m = %ld actuators (modes), n = %ld sensors\n", m, n);
  fflush(stdout);

  matrix_DtraD_eval = gsl_vector_alloc (m); 
  matrix_D = gsl_matrix_alloc (n,m);
  matrix_Ds = gsl_matrix_alloc (m,n);
  matrix_Dtra = gsl_matrix_alloc (m,n);
  matrix_DtraD = gsl_matrix_alloc (m,m); 
  matrix_DtraDinv = gsl_matrix_alloc (m,m); 
  matrix_DtraD_evec = gsl_matrix_alloc (m,m);
  


  /* write matrix_D */
  for(k=0;k<m;k++)
    for(ii=0;ii<n;ii++)
      gsl_matrix_set (matrix_D, ii, k, data.image[ID_Rmatrix].array.F[k*n+ii]);

  /* compute DtraD */
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, matrix_D, matrix_D, 0.0, matrix_DtraD);


  /* compute the inverse of DtraD */

  /* first, compute the eigenvalues and eigenvectors */
  w =   gsl_eigen_symmv_alloc (m);
  matrix_save = gsl_matrix_alloc (m,m);
  gsl_matrix_memcpy(matrix_save, matrix_DtraD);
  gsl_eigen_symmv (matrix_save, matrix_DtraD_eval, matrix_DtraD_evec, w);
  gsl_matrix_free(matrix_save);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort (matrix_DtraD_eval, matrix_DtraD_evec, GSL_EIGEN_SORT_ABS_DESC);

  printf("Eigenvalues\n");
  fflush(stdout);

  // Write eigenvalues
  if((fp=fopen("eigenv.dat","w"))==NULL)
    {
      printf("ERROR: cannot create file \"eigenv.dat\"\n");
      exit(0);
    }
  for(k=0; k<m; k++)
    fprintf(fp,"%ld %g\n", k, gsl_vector_get(matrix_DtraD_eval,k));
  fclose(fp);
  
  for(k=0; k<m; k++)
    printf("Mode %ld eigenvalue = %g\n", k, gsl_vector_get(matrix_DtraD_eval,k));
  
  
  // Write rotation matrix to go from DM modes to eigenmodes
  arraysizetmp[0] = m;
  arraysizetmp[1] = m;
  ID_VTmatrix = create_image_ID(ID_VTmatrix_name, 2, arraysizetmp, FLOAT, 0, 0);
  for(ii=0;ii<m;ii++) // modes
    for(k=0;k<m;k++) // modes
      data.image[ID_VTmatrix].array.F[k*m+ii] = (float) gsl_matrix_get( matrix_DtraD_evec, k, ii);
  
  
  /* second, build the "inverse" of the diagonal matrix of eigenvalues (matrix1) */
  matrix1 = gsl_matrix_alloc (m, m);
  for(ii1=0; ii1<m; ii1++)
    for(jj1=0; jj1<m; jj1++)
      {
	if(ii1==jj1)
	  {
	    if((m-ii1-1)<NB_MODE_REMOVED)
	      gsl_matrix_set(matrix1, ii1, jj1, 0.0);
	    else
	      gsl_matrix_set(matrix1, ii1, jj1, 1.0/gsl_vector_get(matrix_DtraD_eval,ii1));
	  }
	else
	  gsl_matrix_set(matrix1, ii1, jj1, 0.0);
      }
 
  printf("Compute inverse\n");
  fflush(stdout);

  /* third, compute the "inverse" of DtraD */
  matrix2 = gsl_matrix_alloc (m, m);
  //  printf("step 0\n");
  // fflush(stdout);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, matrix_DtraD_evec, matrix1, 0.0, matrix2);
  // printf("step 1\n");
  //fflush(stdout);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, matrix2, matrix_DtraD_evec, 0.0, matrix_DtraDinv);
  //printf("step 2\n");
  //fflush(stdout);
  gsl_matrix_free(matrix1);
  gsl_matrix_free(matrix2);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, matrix_DtraDinv, matrix_D, 0.0, matrix_Ds);

  printf("Write result\n");
  fflush(stdout);

  arraysizetmp[0] = AOconf[loop].sizexWFS;
  arraysizetmp[1] = AOconf[loop].sizexWFS;
  arraysizetmp[2] = m;
  ID_Cmatrix = create_image_ID(ID_Cmatrix_name, 3, arraysizetmp, FLOAT, 0, 0);


  /* write result */
  for(ii=0;ii<n;ii++) // sensors
    for(k=0;k<m;k++) // actuator modes
      data.image[ID_Cmatrix].array.F[k*n+ii] = (float) gsl_matrix_get(matrix_Ds, k, ii);
      

  gsl_vector_free(matrix_DtraD_eval);
  gsl_matrix_free(matrix_D);
  gsl_matrix_free(matrix_Ds);
  gsl_matrix_free(matrix_Dtra);
  gsl_matrix_free(matrix_DtraD);
  gsl_matrix_free(matrix_DtraDinv);
  gsl_matrix_free(matrix_DtraD_evec);

  free(arraysizetmp);

  return(ID_Cmatrix);
}










int AOloopControl_InitializeMemory()
{
  int SM_fd;
  struct stat file_stat;
  int create = 0;
  int result;
  long loop;

  SM_fd = open(AOconfname, O_RDWR);
  if(SM_fd==-1)
    {
      printf("Cannot import file -> creating file\n");
      create = 1;
    }
  else
    {
      fstat(SM_fd, &file_stat);
      printf("File %s size: %zd\n", AOconfname, file_stat.st_size);
      if(file_stat.st_size!=sizeof(AOLOOPCONTROL_CONF)*NB_AOloopcontrol)
	{
	  printf("File size wrong -> recreating file\n");
	  create = 1;
	  close(SM_fd);
	}
    }
  
  if(create==1)
    {
      SM_fd = open(AOconfname, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
	 
      if (SM_fd == -1) {
	perror("Error opening file for writing");
	exit(0);
      }

      result = lseek(SM_fd, sizeof(AOLOOPCONTROL_CONF)*NB_AOloopcontrol-1, SEEK_SET);
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

  AOconf = (AOLOOPCONTROL_CONF*) mmap(0, sizeof(AOLOOPCONTROL_CONF)*NB_AOloopcontrol, PROT_READ | PROT_WRITE, MAP_SHARED, SM_fd, 0);
  if (AOconf == MAP_FAILED) {
    close(SM_fd);
    perror("Error mmapping the file");
    exit(0);
  }

  if(create==1)
    {
      for(loop=0; loop<NB_AOloopcontrol; loop++)
	{
	  AOconf[loop].init = 0;
	  AOconf[loop].on = 0;
	  AOconf[loop].cnt = 0;	  
	  AOconf[loop].maxlimit = 0.1;
	  AOconf[loop].gain = 0.0;
	}
    }
  else
    {
      for(loop=0; loop<NB_AOloopcontrol; loop++)
	if(AOconf[loop].init == 1)
	  {
	    
	    printf("----- Loop %ld   (%s) ----------\n", loop, AOconf[loop].name);
	    printf("  WFS:  %s  [%ld]  %ld x %ld\n", AOconf[loop].WFSname, AOconf[loop].ID_WFS, AOconf[loop].sizexWFS, AOconf[loop].sizeyWFS);
	    printf("   DM:  %s  [%ld]  %ld x %ld\n", AOconf[loop].DMname, AOconf[loop].ID_DM, AOconf[loop].sizexDM, AOconf[loop].sizeyDM);
	  }
    }
  
  AOloopcontrol_meminit = 1;


  return 0;
}






int Average_cam_frames(long loop, long NbAve)
{
  long imcnt;
  long ii;
  double total;
  char name[200];
  int atype;

  atype = data.image[AOconf[loop].ID_WFS].md[0].atype;


  if(NbAve>1)
    for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
      data.image[AOconf[loop].ID_WFS1].array.F[ii] = 0.0;



  switch (atype) {
  case FLOAT :
    imcnt = 0;
    while(imcnt<NbAve)
      {
	usleep(100);	  
	if(AOconf[loop].WFScnt!=data.image[AOconf[loop].ID_WFS].md[0].cnt0)
	  {
	    AOconf[loop].WFScnt = data.image[AOconf[loop].ID_WFS].md[0].cnt0;
	    if(NbAve>1)
	      {
		for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
		  data.image[AOconf[loop].ID_WFS1].array.F[ii] += data.image[AOconf[loop].ID_WFS].array.F[ii];
	      }
	    else
	      {
		for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
		  data.image[AOconf[loop].ID_WFS1].array.F[ii] = data.image[AOconf[loop].ID_WFS].array.F[ii];
	      }
	    imcnt++;
	  }      
      }
    break;
  case USHORT :
    imcnt = 0;
    while(imcnt<NbAve)
      {
	usleep(100);	  
	if(AOconf[loop].WFScnt!=data.image[AOconf[loop].ID_WFS].md[0].cnt0)
	  {
	    AOconf[loop].WFScnt = data.image[AOconf[loop].ID_WFS].md[0].cnt0;
	    if(NbAve>1)
	      {
		for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
		  data.image[AOconf[loop].ID_WFS1].array.F[ii] += data.image[AOconf[loop].ID_WFS].array.U[ii];
	      }
	    else
	      {
		for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
		  data.image[AOconf[loop].ID_WFS1].array.F[ii] = data.image[AOconf[loop].ID_WFS].array.U[ii];
	      }
	    imcnt++;
	  }
      }         
    break;
  default :
    printf("ERROR: DATA TYPE NOT SUPPORTED\n");
    exit(0);
    break;
  }
  

  if(NbAve>1)
    for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
      data.image[AOconf[loop].ID_WFS1].array.F[ii] /= NbAve;

  // SUBTRACT DARK

  // Normalize  
  sprintf(name, "imWFS1_%ld", loop);
  total = arith_image_total(name);
  
  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
    data.image[AOconf[loop].ID_WFS1].array.F[ii] /= total;

  total = arith_image_total(name);
  
  return(0);
}






long AOloopControl_MakeDMModes(long loop, long NBmodes, char *IDname)
{
  long ID;
  long IDtmp;
  long ii, jj;
  double x, y;
  long size = AOconf[loop].sizexDM;
  long m;
  float rpix;
  long size2;

  size2 = size*size;
  rpix = 0.5*size;

  ID = create_3Dimage_ID_float(IDname, size, size, NBmodes);
  
  for(m=0;m<NBmodes;m++)
    {
      IDtmp = mk_zer("zertmp", size, m+1, rpix);
      for(ii=0;ii<size2;ii++)
	data.image[ID].array.F[size2*m+ii] = data.image[IDtmp].array.F[ii];
      delete_image_ID("zertmp");
    }

  return(ID);
}



//
// read loop configuration file
//
int AOloopControl_loadconfigure(long loop, char *config_fname)
{
  FILE *fp;
  char content[200];
  char name[200];
  char fname[200];
  long ID;
  long *sizearray;
  int OK;

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  sizearray = (long*) malloc(sizeof(long)*3);


  if(read_config_parameter(config_fname, "loopname", content)==0)
    exit(0);
  printf("loop name : %s\n", content);
  strcpy(AOconf[loop].name, content);

  // Connect to WFS camera
  if(read_config_parameter(config_fname, "WFSname", content)==0)
    exit(0);
  printf("WFS file name : %s\n", content);
  strcpy(AOconf[loop].WFSname, content);
  AOconf[loop].ID_WFS = read_sharedmem_image(content);
  AOconf[loop].sizexWFS = data.image[AOconf[loop].ID_WFS].md[0].size[0];
  AOconf[loop].sizeyWFS = data.image[AOconf[loop].ID_WFS].md[0].size[1];
  AOconf[loop].sizeWFS = AOconf[loop].sizexWFS*AOconf[loop].sizeyWFS;
  


  // Create WFS1 image memory (averaged, dark-subtracted)
  sprintf(name, "imWFS1_%ld", loop);
  sizearray[0] =  AOconf[loop].sizexWFS;
  sizearray[1] =  AOconf[loop].sizeyWFS;
  AOconf[loop].ID_WFS1 = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);

  // Create WFS1 image memory (averaged, dark-subtracted)
  sprintf(name, "imWFS2_%ld", loop);
  sizearray[0] =  AOconf[loop].sizexWFS;
  sizearray[1] =  AOconf[loop].sizeyWFS;
  AOconf[loop].ID_WFS2 = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);



  // Create WFSref image memory (averaged, dark-subtracted)
  sprintf(name, "refWFS_%ld", loop);
  sizearray[0] =  AOconf[loop].sizexWFS;
  sizearray[1] =  AOconf[loop].sizeyWFS;
  AOconf[loop].ID_refWFS = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);

  if(read_config_parameter(config_fname, "WFSrefim", content)==0)
    exit(0);
  printf("WFS ref file name : %s\n", content);
  OK = 0;
  AOconf[loop].init_refWFS = 0;
  if(is_fits_file(content)==1)
    if((ID=load_fits(content, "tmprefwfs"))!=-1)
      {
	printf("Verifying file\n");
	OK = 1;
	if(data.image[ID].md[0].naxis!=2)
	  {
	    printf("refWFS has wrong dimension\n");
	    OK = 0;
	  }
	if(data.image[ID].md[0].atype!=FLOAT)
	  {
	    printf("refWFS has wrong type\n");
	    OK = 0;
	  }
	if(OK==1)
	  {
	    if(data.image[ID].md[0].size[0]!=AOconf[loop].sizexWFS)
	      {
		printf("refWFS has wrong x size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexWFS);
		OK = 0;
	    }
	    if(data.image[ID].md[0].size[1]!=AOconf[loop].sizeyWFS)
	      {
		printf("refWFS has wrong y size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexWFS);
		OK = 0;
	      }
	  }
	if(OK==1)
	  {
	    sprintf(name, "refWFS_%ld", loop);
	    sizearray[0] =  AOconf[loop].sizexWFS;
	    sizearray[1] =  AOconf[loop].sizeyWFS;
	    AOconf[loop].ID_refWFS = create_image_ID(name, 2, sizearray, FLOAT, 1, 0);
	    list_image_ID();
	    copy_image_ID("tmprefwfs", name);
	    AOconf[loop].init_refWFS = 1;
	  }
	delete_image_ID("tmprefwfs");
      }  



  // Connect to DM 
  if(read_config_parameter(config_fname, "DMname", content)==0)
    exit(0);
  printf("DM file name : %s\n", content);
  strcpy(AOconf[loop].DMname, content);
  AOconf[loop].ID_DM = read_sharedmem_image(content);
  AOconf[loop].sizexDM = data.image[AOconf[loop].ID_DM].md[0].size[0];
  AOconf[loop].sizeyDM = data.image[AOconf[loop].ID_DM].md[0].size[1];
  AOconf[loop].sizeDM = AOconf[loop].sizexDM*AOconf[loop].sizeyDM;

 

  // Load DM modes (will exit if not successful)
  sprintf(fname,"DMmodes_%ld.fits", loop);
  sprintf(name,"DMmodes_%ld", loop);
  if(read_config_parameter(config_fname, "DMmodes", content)==0)
    exit(0);
  printf("DM modes file name : %s\n", content);

  OK = 0;
  if(is_fits_file(content)==1)
    if((ID=load_fits(content, "tmpdmmodes"))!=-1)
      {
	printf("Verifying file\n");
	OK = 1;
	if(data.image[ID].md[0].naxis != 3)
	  {
	    printf("DM modes has wrong dimension\n");
	    OK = 0;
	  }
	if(data.image[ID].md[0].atype != FLOAT)
	  {
	    printf("DM modes has wrong type\n");
	    OK = 0;
	  }
	if(OK==1)
	  {
	    if(data.image[ID].md[0].size[0]!=AOconf[loop].sizexDM)
	      {
		printf("DM modes has wrong x size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexDM);
		OK = 0;
	    }
	    if(data.image[ID].md[0].size[1]!=AOconf[loop].sizeyDM)
	      {
		printf("DM modes has wrong y size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexDM);
		OK = 0;
	      }
	  }
	if(OK==1)
	  {
	    AOconf[loop].NBDMmodes = data.image[ID].md[0].size[2];
	    sprintf(name, "DMmodes_%ld", loop);
	    sizearray[0] =  AOconf[loop].sizexDM;
	    sizearray[1] =  AOconf[loop].sizeyDM;
	    sizearray[2] =  AOconf[loop].NBDMmodes;
	    AOconf[loop].ID_DMmodes = create_image_ID(name, 3, sizearray, FLOAT, 1, 0);
	    list_image_ID();
	    copy_image_ID("tmpdmmodes", name);
	  }
	delete_image_ID("tmpdmmodes");
      }
  if(OK == 0)
    {
      printf("\n");
      printf("========== ERROR: NEED DM MODES TO START AO LOOP ===========\n");	
      printf("\n");
      exit(0);
    }
	


  // AOconf[loop].ID_DMmodes = AOloopControl_MakeDMModes(loop, 5, name);

  printf("%ld modes\n", AOconf[loop].NBDMmodes);








  // Create modal command vector memory
  sprintf(name, "DMmode_cmd_%ld", loop);
  sizearray[0] =  AOconf[loop].NBDMmodes;
  AOconf[loop].ID_cmd_modes = create_image_ID(name, 1, sizearray, FLOAT, 1, 0);
  sprintf(name, "DMmode_cmd1_%ld", loop);
  sizearray[0] =  AOconf[loop].NBDMmodes;
  AOconf[loop].ID_cmd1_modes = create_image_ID(name, 1, sizearray, FLOAT, 1, 0);







  // Create Response Matrix shared memory
  sprintf(name, "RespM_%ld", loop);
  sizearray[0] =  AOconf[loop].sizexWFS;
  sizearray[1] =  AOconf[loop].sizeyWFS;
  sizearray[2] =  AOconf[loop].NBDMmodes;
  AOconf[loop].ID_respM = create_image_ID(name, 3, sizearray, FLOAT, 1, 0);
  

  // Load response matrix if possible
  AOconf[loop].init_RM = 0;
  if(read_config_parameter(config_fname, "RespMatrix", content)==0)
    exit(0);
  printf("Response Matrix file name : %s\n", content);
  if((ID=load_fits(content, "tmprespM"))!=-1)
    {
      // check size is OK
      OK = 1;
      if(data.image[ID].md[0].naxis!=3)
	{
	  printf("Response matrix has wrong dimension\n");
	  OK = 0;
	}
      if(data.image[ID].md[0].atype!=FLOAT)
	{
	  printf("Response matrix has wrong type\n");
	  OK = 0;
	}
     if(OK==1)
	{
	  if(data.image[ID].md[0].size[0]!=AOconf[loop].sizexWFS)
	    {
	      printf("Response matrix has wrong x size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexWFS);
	      OK = 0;
	    }
	  if(data.image[ID].md[0].size[1]!=AOconf[loop].sizeyWFS)
	    {
	      printf("Response matrix has wrong y size\n");
	      OK = 0;
	    }
	  if(data.image[ID].md[0].size[2]!=AOconf[loop].NBDMmodes)
	    {
	      printf("Response matrix has wrong z size\n");
	      OK = 0;
	    }	  
	}
      if(OK==1)
	{
	  AOconf[loop].init_RM = 1;
	  sprintf(name, "RespM_%ld", loop);
	  copy_image_ID("tmprespM", name);
	}
      delete_image_ID("tmprespM");
    }

  


  // Create Control Matrix shared memory
  sprintf(name, "ContrM_%ld", loop);
  sizearray[0] =  AOconf[loop].sizexWFS;
  sizearray[1] =  AOconf[loop].sizeyWFS;
  sizearray[2] =  AOconf[loop].NBDMmodes;
  AOconf[loop].ID_contrM = create_image_ID(name, 3, sizearray, FLOAT, 1, 0);
  

  // Load response matrix if possible
  AOconf[loop].init_CM = 0;
    if(read_config_parameter(config_fname, "ContrMatrix", content)==0)
    exit(0);
  printf("Control Matrix file name : %s\n", content);
  if((ID=load_fits(content, "tmpcontrM"))!=-1)
    {
      // check size is OK
      OK = 1;
      if(data.image[ID].md[0].naxis!=3)
	{
	  printf("Control matrix has wrong dimension\n");
	  OK = 0;
	}
      if(data.image[ID].md[0].atype!=FLOAT)
	{
	  printf("Control matrix has wrong type\n");
	  OK = 0;
	}
     if(OK==1)
	{
	  if(data.image[ID].md[0].size[0]!=AOconf[loop].sizexWFS)
	    {
	      printf("Control matrix has wrong x size : is %ld, should be %ld\n", data.image[ID].md[0].size[0], AOconf[loop].sizexWFS);
	      OK = 0;
	    }
	  if(data.image[ID].md[0].size[1]!=AOconf[loop].sizeyWFS)
	    {
	      printf("Control matrix has wrong y size\n");
	      OK = 0;
	    }
	  if(data.image[ID].md[0].size[2]!=AOconf[loop].NBDMmodes)
	    {
	      printf("Control matrix has wrong z size\n");
	      OK = 0;
	    }	  
	}
      if(OK==1)
	{
	  AOconf[loop].init_CM = 1;
	  sprintf(name, "ContrM_%ld", loop);
	  copy_image_ID("tmpcontrM", name);
	}
      delete_image_ID("tmpcontrM");
    }
 
  



  list_image_ID();



  
  free(sizearray);

  AOconf[loop].init = 1;


  printf("   init_WFSref    %d\n", AOconf[loop].init_refWFS);
  printf("   init_RM        %d\n", AOconf[loop].init_RM);
  printf("   init_CM        %d\n", AOconf[loop].init_CM);

  AOconf[loop].init = 1;

  return(0);
}




int set_DM_modes(long loop)
{
  long k;
  long i, j;

  data.image[AOconf[loop].ID_DM].md[0].write = 1;

  for(j=0;j<AOconf[loop].sizeDM;j++)
    data.image[AOconf[loop].ID_DM].array.F[j] = 0.0;

  for(k=0; k < AOconf[loop].NBDMmodes; k++)
    {
      for(i=0;i<AOconf[loop].sizeDM;i++)
	data.image[AOconf[loop].ID_DM].array.F[i] += data.image[AOconf[loop].ID_cmd_modes].array.F[k] * data.image[AOconf[loop].ID_DMmodes].array.F[k*AOconf[loop].sizeDM+i];
    }

  data.image[AOconf[loop].ID_DM].md[0].cnt0++;
  data.image[AOconf[loop].ID_DM].md[0].write = 0;


  AOconf[loop].DMupdatecnt ++;

  return(0);
}




// measures response matrix AND reference 

int Measure_Resp_Matrix(long loop, long NbAve, float amp, long nbloop)
{
  long NBloops;
  long kloop;
  long delayus = 100; // delay in us
  long ii, i, imax;
  int Verbose = 1;
  long k1, k, k2;
  char fname[200];
  char name0[200];
  char name[200];

  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();
  
  printf("SETTING UP...\n");
  sprintf(fname, "AOloop%ld.conf", LOOPNUMBER);
  AOloopControl_loadconfigure(LOOPNUMBER, fname);
  
  NBloops = nbloop;

  // initialize RM to zero
  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
    for(k=0;k<AOconf[loop].NBDMmodes;k++)
      data.image[AOconf[loop].ID_respM].array.F[k*AOconf[loop].sizeWFS+ii] = 0.0;
  

  // initialize reference to zero
  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
    data.image[AOconf[loop].ID_refWFS].array.F[ii] = 0.0;


  printf("\n");
  printf("Testing (in measure_resp_matrix function) :,  NBloops = %ld, NBmode = %ld\n",  NBloops, AOconf[loop].NBDMmodes);



  arith_image_zero(AOconf[loop].WFSname);
  for(k2 = 0; k2 < AOconf[loop].NBDMmodes; k2++)
    data.image[AOconf[loop].ID_cmd_modes].array.F[k2] = 0.0;
  
  for (kloop = 0; kloop < NBloops; kloop++) // number of loops
    {
      if(Verbose)
	{
	  printf("\n Loop %ld / %ld (%f)\n", kloop, NBloops, amp);
	  fflush(stdout);
	}

      for(k1 = 0; k1 < AOconf[loop].NBDMmodes; k1++)
	{	 
	  for(k2 = 0; k2 < AOconf[loop].NBDMmodes; k2++)
	    data.image[AOconf[loop].ID_cmd_modes].array.F[k2] = 0.0;
	 
	  // positive 
	  data.image[AOconf[loop].ID_cmd_modes].array.F[k1] = amp;
	  set_DM_modes(loop);

	  usleep(delayus);
	  Average_cam_frames(loop, NbAve);
	 
	  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	    {
	      data.image[AOconf[loop].ID_respM].array.F[k1*AOconf[loop].sizeWFS+ii] += data.image[AOconf[loop].ID_WFS1].array.F[ii];
	      data.image[AOconf[loop].ID_refWFS].array.F[ii] += data.image[AOconf[loop].ID_WFS1].array.F[ii];
	    }

	  // negative
	  data.image[AOconf[loop].ID_cmd_modes].array.F[k1] = 0.0-amp;
	  set_DM_modes(loop);
	  usleep(delayus);
	  Average_cam_frames(loop, NbAve);
	 
	  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
	    {
	      data.image[AOconf[loop].ID_respM].array.F[k1*AOconf[loop].sizeWFS+ii] -= data.image[AOconf[loop].ID_WFS1].array.F[ii];
	      data.image[AOconf[loop].ID_refWFS].array.F[ii] += data.image[AOconf[loop].ID_WFS1].array.F[ii];
	    }
	  
	}
    }

  set_DM_modes(loop);

  printf("Acquisition done, compiling results...");
  fflush(stdout);

  arith_image_zero(AOconf[loop].WFSname);
  for(k2 = 0; k2 < AOconf[loop].NBDMmodes; k2++)
    data.image[AOconf[loop].ID_cmd_modes].array.F[k2] = 0.0;


  
  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
    for(k1=0;k1<AOconf[loop].NBDMmodes;k1++)
      data.image[AOconf[loop].ID_respM].array.F[k1*AOconf[loop].sizeWFS+ii] /= (NBloops*2.0*amp);

  for(ii=0;ii<AOconf[loop].sizeWFS;ii++)
    data.image[AOconf[loop].ID_refWFS].array.F[ii] /= (NBloops*2.0*AOconf[loop].NBDMmodes);

  printf("\n");
  fflush(stdout);

  printf("Computing Control matrix(es) ... ");
  fflush(stdout);


  imax = 5;
  if(imax<AOconf[loop].NBDMmodes-1)
    imax = AOconf[loop].NBDMmodes-1;

  for(i=0;i<imax;i++)
    {
      printf("[%ld] ", i);
      fflush(stdout);
      sprintf(name0, "RespM_%ld", loop);
      sprintf(name, "cmat%ld", i);
      compute_ControlMatrix(LOOPNUMBER, i, name0, name, "evecM");
    }

  printf("Done\n");

  return(0);
}




int ControlMatrixMultiply( float *cm_array, float *imarray, long m, long n, float *outvect)
{
  long i;

  cblas_sgemv (CblasRowMajor, CblasNoTrans, m, n, 1.0, cm_array, n, imarray, 1, 0.0, outvect, 1);


  return(0);
}




int AOcompute(long loop)
{
  float total = 0.0;
  long k, k1, k2;
  long ii;
  
  long m, n;
  

  // get dark-subtracted image
  Average_cam_frames(loop, 1);
  
  AOconf[loop].status = 4;  // 4: REMOVING REF
  

  for(ii=0; ii<AOconf[loop].sizeWFS; ii++)
    data.image[AOconf[loop].ID_WFS2].array.F[ii] = data.image[AOconf[loop].ID_WFS1].array.F[ii] - data.image[AOconf[loop].ID_refWFS].array.F[ii];

  AOconf[loop].status = 5; // MULTIPLYING BY CONTROL MATRIX -> MODE VALUES


  ControlMatrixMultiply( data.image[AOconf[loop].ID_contrM].array.F, data.image[AOconf[loop].ID_WFS2].array.F, AOconf[loop].NBDMmodes, AOconf[loop].sizeWFS, data.image[AOconf[loop].ID_cmd1_modes].array.F);
  data.image[AOconf[loop].ID_WFS2].md[0].cnt0 ++;

 
  if(1)
    {
      printf("GSL MATRIX MULT  :  ");
      for(k=0; k<AOconf[loop].NBDMmodes; k++)
	printf(" %6.3f ",data.image[AOconf[loop].ID_cmd1_modes].array.F[k]);
      printf("\n");       

       
      printf("Conventional mat mult %d %d\n", M, N);
      for(m=0; m<AOconf[loop].NBDMmodes; m++)
	{
	  data.image[AOconf[loop].ID_cmd1_modes].array.F[m] = 0.0;
	  for(n=0; n<AOconf[loop].sizeWFS; n++)
	    data.image[AOconf[loop].ID_cmd1_modes].array.F[m] += cMat[n*AOconf[loop].NBDMmodes+m]*wfsVec[n];
	}	  

      printf("CONV MATRIX MULT :  ");
      for(k=0; k<AOconf[loop].NBDMmodes; k++)
	printf(" %6.3f ",data.image[AOconf[loop].ID_cmd1_modes].array.F[k]);
      printf("\n");          
    }

  AOconf[loop].RMSmodes = 0;
  for(k=0; k<AOconf[loop].NBDMmodes; k++)
    AOconf[loop].RMSmodes += data.image[AOconf[loop].ID_cmd1_modes].array.F[k]*data.image[AOconf[loop].ID_cmd1_modes].array.F[k];
    

  AOconf[loop].status = 6; //  MULTIPLYING BY GAINS
   

  for(k=0; k<AOconf[loop].NBDMmodes; k++)
    {
      data.image[AOconf[loop].ID_cmd_modes].array.F[k] -= AOconf[loop].gain * data.image[AOconf[loop].ID_cmd1_modes].array.F[k];
      
      if(data.image[AOconf[loop].ID_cmd_modes].array.F[k] < -AOconf[loop].maxlimit)
	data.image[AOconf[loop].ID_cmd_modes].array.F[k] = -AOconf[loop].maxlimit;
      
      if(data.image[AOconf[loop].ID_cmd_modes].array.F[k] > AOconf[loop].maxlimit)
	data.image[AOconf[loop].ID_cmd_modes].array.F[k] = AOconf[loop].maxlimit;
    }	


  /*
    for(k=0; k<smao[0].NBmode; k++)
      {
	smao[0].mon_mode_residual[k] += smao[0].cmmode_array1[k];
	smao[0].mon_mode_residual_cnt[k]++;
	
	smao[0].mon_mode_residualRMS[k] += smao[0].cmmode_array1[k]*smao[0].cmmode_array1[k];
	smao[0].mon_mode_residualRMS_cnt[k]++;
	
	smao[0].mon_mode_applied[k] += smao[0].command_modes[k];
	smao[0].mon_mode_applied_cnt[k]++;
      }
  */

    AOconf[loop].status = 7; // fill in covariance matrix
    /*    if(smao[0].covM_on == 1)
      {
	for(k1=0; k1<smao[0].NBmode; k1++)
	  for(k2=0; k2<k1+1; k2++)
	    smao[0].covM[k2*smao[0].NBmode+k1] += smao[0].cmmode_array1[k1] * smao[0].cmmode_array1[k2];
	smao[0].covMcnt++;
	}*/

  return(0);
}


int AOloopControl_run()
{
  char fname[200];
  long loop;
  int OK;

  loop = LOOPNUMBER;
  
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  printf("SETTING UP...\n");
  sprintf(fname, "AOloop%ld.conf", LOOPNUMBER);
  AOloopControl_loadconfigure(LOOPNUMBER, fname);
 

  OK = 1;
  if(AOconf[loop].init_refWFS==0)
    {
      printf("ERROR: CANNOT RUN LOOP WITHOUT WFS REFERENCE\n");
      OK = 0;
    }
  if(AOconf[loop].init_CM==0)
    {
      printf("ERROR: CANNOT RUN LOOP WITHOUT CONTROL MATRIX\n");
      OK = 0;
    }
  
  if(OK==1)
    {
      AOconf[loop].kill = 0;
      AOconf[loop].on = 0;
      printf("\n");
      while( AOconf[loop].kill == 0)
	{
	  printf(" WAITING                    \r");
	  fflush(stdout);
	  usleep(100);
	  
	  while(AOconf[loop].on == 1)
	    {
	      printf("LOOP IS RUNNING  %llu  %g      Gain = %f \r", AOconf[loop].cnt, AOconf[loop].RMSmodes, AOconf[loop].gain);
	      fflush(stdout);
	      usleep(10000);
	      AOcompute(loop);
	      set_DM_modes(loop);
	      AOconf[loop].cnt++;
	    }
	}
    }

  /*  
      printf("Acquiring image\n");
  Average_cam_frames(LOOPNUMBER, 10);
  save_fl_fits("imWFS1_0","!imave.fits");
  */
  
  /*  printf("Acquiring response matrix\n");
  Measure_Resp_Matrix(LOOPNUMBER, 10, 1.0, 3);
  
  save_fl_fits("RespM_0", "!respm0.fits");
  save_fl_fits("refWFS_0", "!refwfs0.fits");

  printf("Computing control matrix\n");
  compute_ControlMatrix(LOOPNUMBER, 1, "RespM_0", "ContrM_0", "evecM");
  save_fl_fits("ContrM_0", "!ContrM_0.fits");
  save_fl_fits("evecM", "!evecM.fits");
  */

  

  list_image_ID();

  return(0);
}





int AOloopControl_setLoopNumber(long loop)
{
  printf("LOOPNUMBER = %ld\n", loop);
  LOOPNUMBER = loop;

  return 0;
}

int AOloopControl_loopkill()
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].kill = 1;

  return 0;
}

int AOloopControl_loopon()
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].on = 1;

  return 0;
}

int AOloopControl_loopoff()
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].on = 0;

  return 0;
}

int AOloopControl_setgain(float gain)
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].gain = gain;

  return 0;
}

int AOloopControl_setmaxlimit(float maxlimit)
{
  if(AOloopcontrol_meminit==0)
    AOloopControl_InitializeMemory();

  AOconf[LOOPNUMBER].maxlimit = maxlimit;

  return 0;
}
