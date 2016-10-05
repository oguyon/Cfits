#include <fitsio.h> 
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <semaphore.h>
#include <sched.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#ifdef __MACH__
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct timespec *t){
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    uint64_t time;
    time = mach_absolute_time();
    double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
    double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1e9);
    t->tv_sec = seconds;
    t->tv_nsec = nseconds;
    return 0;, optres[ii]
}
#else
#include <time.h>
#endif


#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "statistic/statistic.h"
#include "info/info.h"
#include "linopt_imtools/linopt_imtools.h"

#include "linARfilterPred/linARfilterPred.h"

extern DATA data;




// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string (not image)
// 4: existing image
// 5: string


int LINARFILTERPRED_SelectBlock_cli()
{

	if(CLI_checkarg(1,4)+CLI_checkarg(2,4)+CLI_checkarg(3,2)+CLI_checkarg(4,3)==0)
		{
			LINARFILTERPRED_SelectBlock(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.string);
		}
	else
		return(1);
}


int LINARFILTERPRED_Build_LinPredictor_cli()
{
	if(CLI_checkarg(1,4)+CLI_checkarg(2,2)+CLI_checkarg(3,1)+CLI_checkarg(4,1)+CLI_checkarg(5,1)+CLI_checkarg(6,3)==0)
		LINARFILTERPRED_Build_LinPredictor(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numf, data.cmdargtoken[5].val.numf, data.cmdargtoken[6].val.string, 1);
	else
       return 1;

  return(0);
}

int LINARFILTERPRED_Apply_LinPredictor_cli()
{
	if(CLI_checkarg(1,4)+CLI_checkarg(2,4)+CLI_checkarg(3,1)+CLI_checkarg(4,3)==0)
		LINARFILTERPRED_Apply_LinPredictor(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.string);
	else
		return 1;

  return(0);
}

int LINARFILTERPRED_ScanGain_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,1)+CLI_checkarg(3,1)==0)
    LINARFILTERPRED_ScanGain(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numf);
	else
       return 1;

  return(0);
}







int init_linARfilterPred()
{
    strcpy(data.module[data.NBmodule].name, __FILE__);
    strcpy(data.module[data.NBmodule].info, "linear auto-regressive predictive filters");
    data.NBmodule++;


    strcpy(data.cmd[data.NBcmd].key,"mselblock");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = LINARFILTERPRED_SelectBlock_cli;
    strcpy(data.cmd[data.NBcmd].info,"select modes belonging to a block");
    strcpy(data.cmd[data.NBcmd].syntax,"<input mode values> <block map> <selected block> <output>");
    strcpy(data.cmd[data.NBcmd].example,"mselblock modevals blockmap 23 blk23modevals");
    strcpy(data.cmd[data.NBcmd].Ccall,"long LINARFILTERPRED_SelectBlock(char *IDin_name, char *IDblknb_name, long blkNB, char *IDout_name)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"mkARpfilt");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = LINARFILTERPRED_Build_LinPredictor_cli;
    strcpy(data.cmd[data.NBcmd].info,"Make linear auto-regressive filter");
    strcpy(data.cmd[data.NBcmd].syntax,"<input data> <PForder> <PFlag> <SVDeps> <regularization param> <output filters>");
    strcpy(data.cmd[data.NBcmd].example,"mkARpfilt indata 5 2.4 0.0001 0.0 outPF");
    strcpy(data.cmd[data.NBcmd].Ccall,"int LINARFILTERPRED_Build_LinPredictor(char *IDin_name, long PForder, float PFlag, double SVDeps, double RegLambda, char *IDoutPF, int outMode)");
    data.NBcmd++;


    strcpy(data.cmd[data.NBcmd].key,"applyARpfilt");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = LINARFILTERPRED_Apply_LinPredictor_cli;
    strcpy(data.cmd[data.NBcmd].info,"Apply linear auto-regressive filter");
    strcpy(data.cmd[data.NBcmd].syntax,"<input data> <predictor> <PFlag> <prediction>");
    strcpy(data.cmd[data.NBcmd].example,"applyARpfilt indata Pfilt 2.4 outPF");
    strcpy(data.cmd[data.NBcmd].Ccall,"long LINARFILTERPRED_Apply_LinPredictor(char *IDfilt_name, char *IDin_name, float PFlag, char *IDout_name)");
    data.NBcmd++;

	
	strcpy(data.cmd[data.NBcmd].key,"mscangain");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = LINARFILTERPRED_ScanGain_cli;
    strcpy(data.cmd[data.NBcmd].info,"scan gain");
    strcpy(data.cmd[data.NBcmd].syntax,"<input mode values> <multiplicative factor (leak)> <latency [frame]>");
    strcpy(data.cmd[data.NBcmd].example,"mscangain olwfsmeas 0.98 2.65");
    strcpy(data.cmd[data.NBcmd].Ccall,"LINARFILTERPRED_ScanGain(char* IDin_name, float multfact, float framelag)");
    data.NBcmd++;



    // add atexit functions here

    return 0;
}













// select block on first dimension 
long LINARFILTERPRED_SelectBlock(char *IDin_name, char *IDblknb_name, long blkNB, char *IDout_name)
{
	long IDin, IDblknb;
	long naxis, axis;
	long m;
	long NBmodes1;
	long *sizearray;
	long xsize, ysize, zsize;
	long cnt;
	long ii, jj, kk;
	long IDout;
	char imname[200];
	long mmax;
	
	printf("Selecting block %ld ...\n", blkNB);
	fflush(stdout);
	
	IDin = image_ID(IDin_name);
	IDblknb = image_ID(IDblknb_name);
	naxis = data.image[IDin].md[0].naxis;
	
	if(data.image[IDin].md[0].size[0] != data.image[IDblknb].md[0].size[0])
		{			
			printf("WARNING: block index file and telemetry have different sizes\n");
			fflush(stdout);
			mmax = data.image[IDin].md[0].size[0];
			if(data.image[IDblknb].md[0].size[0]<mmax)
				mmax = data.image[IDblknb].md[0].size[0];
		}
	
	
	NBmodes1 = 0;
	for(m=0;m<mmax;m++)
		{
			if(data.image[IDblknb].array.U[m] == blkNB)
				NBmodes1++;
		}
	
	sizearray = (long*) malloc(sizeof(long)*naxis);
	
	for(axis=0;axis<naxis;axis++)
		sizearray[axis] = data.image[IDin].md[0].size[axis];
	sizearray[0] = NBmodes1;
	
	IDout = create_image_ID(IDout_name, naxis, sizearray, FLOAT, 0, 0);
	
	xsize = data.image[IDin].md[0].size[0];
	if(naxis>1)
		ysize = data.image[IDin].md[0].size[1];
	else
		ysize = 1;
	if(naxis>2)
		zsize = data.image[IDin].md[0].size[2];
	else
		zsize = 1;
	
	
	
	cnt = 0;

	for(jj=0;jj<ysize;jj++)
		for(kk=0;kk<zsize;kk++)
			for(ii=0;ii<mmax;ii++)
				if(data.image[IDblknb].array.U[ii] == blkNB)
						{
							//printf("%ld / %ld   cnt = %8ld / %ld\n", ii, xsize, cnt, NBmodes1*ysize*zsize);
							//fflush(stdout);
							data.image[IDout].array.F[cnt] = data.image[IDin].array.F[kk*xsize*ysize + jj*ysize + ii];
							cnt++;
						}		
	
	free(sizearray);
	
		
	return(IDout);
}






//
// IDin_name is a 2D or 3D image
//
// optional: inmask selects input pixels to be used
//           outmask selects output pixel(s) to be used
//
// Note: if atmospheric wavefronts, data should be piston-free 
//
// outMode
//	0: do not write individual filters
//	1: write individual filters
// (note: output filter cube always written)
//
long LINARFILTERPRED_Build_LinPredictor(char *IDin_name, long PForder, float PFlag, double SVDeps, double RegLambda, char *IDoutPF_name, int outMode)
{
    long IDin, IDmatA, IDout, IDinmask;
	long nbspl; // number of samples
	long NBpix;
	long NBmvec, NBmvec1;
	long mvecsize;
	long xsize, ysize;
	long ii, jj;
	long *pixarray_x;
	long *pixarray_y;
	long *pixarray_xy;
	double *ave_inarray;
	int REG = 0;  // 1 if regularization
	long m, m1, pix, k0, dt;
	int Save = 1;
	long xysize;
	long IDmatC;
	int use_magma = 1; // use MAGMA library if available
	int magmacomp = 0;

	long IDfiltC;
	float *valfarray;
	float alpha;
	long PFpix;
	char filtname[200];
	char filtfname[200];
	long ID_Pfilt;
	float val;
	long ind1;
	long IDoutmask;
	int ret;
	long IDoutPF;
	
	
	int DC_MODE = 0; // 1 if average value of each mode is removed
	
	
	// =========== SELECT INPUT VALUES =======================
	
	IDin = image_ID(IDin_name);

	switch (data.image[IDin].md[0].naxis) {
		
		case 2 :
		nbspl = data.image[IDin].md[0].size[1];
		xsize = data.image[IDin].md[0].size[0];
		ysize = 1;
		break;
		
		case 3 :
		nbspl = data.image[IDin].md[0].size[2];
		xsize = data.image[IDin].md[0].size[0];
		ysize = data.image[IDin].md[0].size[1];
		break;
		
		default :
		printf("Invalid image size\n");
		break;
	}
	xysize = xsize*ysize;
	
	pixarray_x = (long*) malloc(sizeof(long)*xsize*ysize);
	pixarray_y = (long*) malloc(sizeof(long)*xsize*ysize);
	pixarray_xy = (long*) malloc(sizeof(long)*xsize*ysize);
	ave_inarray = (double*) malloc(sizeof(double)*xsize*ysize);
	
	IDinmask = image_ID("inmask");
	if(IDinmask==-1)
		{
			NBpix = 0; //xsize*ysize;
			
			for(ii=0;ii<xsize;ii++)
				for(jj=0;jj<ysize;jj++)
					{
					pixarray_x[NBpix] = ii;
					pixarray_y[NBpix] = jj;
					pixarray_xy[NBpix] = jj*xsize+ii;
					NBpix++;
				}
		}
	else
	{
	NBpix = 0;
	for(ii=0;ii<xsize;ii++)
		for(jj=0;jj<ysize;jj++)
			if(data.image[IDinmask].array.F[jj*xsize+ii] > 0.5)
				{
					pixarray_x[NBpix] = ii;
					pixarray_y[NBpix] = jj;
					pixarray_xy[NBpix] = jj*xsize+ii;
					NBpix++;
				}
	}	
	printf("NBpix = %ld\n", NBpix);
	
	
	// ===================== BUILD DATA MATRIX ============================
	// build data matrix
	NBmvec = nbspl - PForder - (int) (PFlag) - 1;
	mvecsize = NBpix * PForder; // size of each sample vector for AR filter, excluding regularization
	
	if(REG==0) // no regularization
		{
			printf("NBmvec   = %ld  -> %ld \n", NBmvec, NBmvec);
			NBmvec1 = NBmvec;
			IDmatA = create_2Dimage_ID("PFmatD", NBmvec, mvecsize);
		}
    else // with regularization
		{
			printf("NBmvec   = %ld  -> %ld \n", NBmvec, NBmvec + mvecsize);
			NBmvec1 = NBmvec + mvecsize;
			IDmatA = create_2Dimage_ID("PFmatD", NBmvec + mvecsize, mvecsize);
		}
		
	IDmatA = image_ID("PFmatD");
 	
 	
 	// each column (ii = cst) is a measurement
    // m index is measurement
    // dt*NBpix+pix index is pixel
   
    printf("mvecsize = %ld  (%ld x %ld)\n", mvecsize, PForder, NBpix);
	printf("NBpix = %ld\n", NBpix);
	printf("NBmvec1 = %ld\n", NBmvec1);
	printf("PForder = %ld\n", PForder);

	printf("xysize = %ld\n", xysize);
	printf("IDin = %ld\n\n", IDin);
	list_image_ID();
	
	
	if(DC_MODE == 1) // remove average
	{
		for(pix=0; pix<NBpix; pix++)
			{
				ave_inarray[pix] = 0.0;
				for(m=0; m<nbspl; m++)
					ave_inarray[pix] += data.image[IDin].array.F[m*xysize+pixarray_xy[pix]];
				ave_inarray[pix] /= nbspl;
			}
	}
	else
	{
		for(pix=0; pix<NBpix; pix++)
			ave_inarray[pix] = 0.0;
	}


	
	for(m=0; m<NBmvec1; m++)
	{
		k0 = m + PForder-1; // dt=0 index
		for(pix=0; pix<NBpix; pix++)
			for(dt=0; dt<PForder; dt++)		
				data.image[IDmatA].array.F[(NBpix*dt+pix)*NBmvec1 + m] = data.image[IDin].array.F[(k0-dt)*xysize + pixarray_xy[pix]] - ave_inarray[pix];
	}
	free(ave_inarray);
	


	if(REG==1)
		{
			for(m=0; m<mvecsize; m++)
				{
					m1 = NBmvec + m;
					data.image[IDmatA].array.F[(m)*NBmvec1+(NBmvec+m)] = RegLambda;
				}
		}

	
	if(Save == 1)
		save_fits("PFmatD", "!PFmatD.fits");
	list_image_ID();

	
	
	// ===================== COMPUTE RECONSTRUCTION MATRIX ============================
	printf("Compute reconstruction matrix\n");
	fflush(stdout);


	
#ifdef HAVE_MAGMA
		CUDACOMP_magma_compute_SVDpseudoInverse("PFmatD", "PFmatC", SVDeps, 100000, "PF_VTmat");
	#else
		linopt_compute_SVDpseudoInverse("PFmatD", "PFmatC", SVDeps, 100000, "PF_VTmat");
	#endif 
	 
	if(Save==1)
        save_fits("PFmatC", "!PFmatC.fits");
    IDmatC = image_ID("PFmatC");

	
	
	// ===================== COMPUTE FILTERS ============================
	printf("Compute filters\n");
	fflush(stdout);
	
	
	
	
	ret = system("mkdir -p pixfilters");
	
	// FILTER MATRIX
	// axis 0 [ii] : input mode
	// axis 1 [jj] : reconstructed mode
	// axis 2 [kk] : time step
	IDoutPF = create_3Dimage_ID(IDoutPF_name, NBpix, NBpix, PForder);	
	
	IDoutmask = image_ID("outmask");
	
	valfarray = (float*) malloc(sizeof(float)*NBmvec);
	alpha = PFlag - ((long) PFlag);
	for(PFpix=0; PFpix<NBpix; PFpix++) // PFpix is the pixel for which the filter is created (axis 1 in cube, jj)
	{
		// INDIVIDUAL FILTERS
		sprintf(filtname, "PFfilt_%06ld_%03ld_%03ld", pixarray_xy[PFpix], pixarray_x[PFpix], pixarray_y[PFpix]);			
		sprintf(filtfname, "!./pixfilters/PFfilt_%06ld_%03ld_%03ld.fits", pixarray_xy[PFpix], pixarray_x[PFpix], pixarray_y[PFpix]);	
		ID_Pfilt = create_3Dimage_ID(filtname, xsize, ysize, PForder);
		
		
		// fill in valfarray
		
		for(m=0; m<NBmvec; m++)
			{
				k0 = m + PForder -1;
				k0 += (long) PFlag;
				
				valfarray[m] = (1.0-alpha)*data.image[IDin].array.F[(k0)*xysize + pixarray_xy[PFpix]] + alpha*data.image[IDin].array.F[(k0+1)*xysize + pixarray_xy[PFpix]];
			}
		
		
		for(pix=0; pix<NBpix; pix++)
			{
				for(dt=0; dt<PForder; dt++)		
					{
						val = 0.0;
						ind1 = (NBpix*dt+pix)*NBmvec1;
						for(m=0; m<NBmvec; m++)
							val += data.image[IDmatC].array.F[ind1+m] * valfarray[m];

						data.image[ID_Pfilt].array.F[xysize*dt + pixarray_xy[pix]] =  val;
						data.image[IDoutPF].array.F[dt*NBpix*NBpix  + PFpix*NBpix + pix] = val;
					}
			}
		save_fits(filtname, filtfname);	
	}
	
	free(valfarray);
 	free(pixarray_x);
 	free(pixarray_y);
 	free(pixarray_xy);
	
    return(IDoutPF);
}


//
// 
// out : prediction
//
// ADDITIONAL OUTPUTS:
// outf : time-shifted measurement
//


long LINARFILTERPRED_Apply_LinPredictor(char *IDfilt_name, char *IDin_name, float PFlag, char *IDout_name)
{
	long IDout;
	long IDin;
	long IDfilt;
	long xsize, ysize, xysize;
	long nbspl;
	long PForder;
	long step;
	long kk;
	float alpha;
	long PFlagl;
	long ii, iip;
	float valp, valf;
	
	long IDoutf;
	
		
	
	IDin = image_ID(IDin_name);
	IDfilt = image_ID(IDfilt_name);
	
	switch (data.image[IDin].md[0].naxis) {
		
		case 2 :
		nbspl = data.image[IDin].md[0].size[1];
		xsize = data.image[IDin].md[0].size[0];
		ysize = 1;
		IDout = create_2Dimage_ID(IDout_name, xsize, nbspl);
		IDoutf = create_2Dimage_ID("outf", xsize, nbspl);
		break;
		
		case 3 :
		nbspl = data.image[IDin].md[0].size[2];
		xsize = data.image[IDin].md[0].size[0];
		ysize = data.image[IDin].md[0].size[1];
		IDout = create_3Dimage_ID(IDout_name, xsize, ysize, nbspl);
		IDoutf = create_3Dimage_ID("outf", xsize, ysize, nbspl);
		break;
		
		default :
		printf("Invalid image size\n");
		break;
	}
	xysize = xsize*ysize;
	
	PForder = data.image[IDfilt].md[0].size[2];
	
	if((data.image[IDfilt].md[0].size[0]!=xysize)||(data.image[IDfilt].md[0].size[1]!=xysize))
		{
			printf("ERROR: filter \"%s\" size is incorrect\n", IDfilt_name);
			exit(0);
		}
	
	alpha = PFlag - ((long) PFlag);
	PFlagl = (long) PFlag;
	
	for(kk=PForder;kk<nbspl;kk++)
	{
		for(iip=0;iip<xysize;iip++) // predicted variable
			{
				valp = 0.0; // prediction 
				for(step=0;step<PForder;step++)
					{
						for(ii=0;ii<xsize*ysize;ii++)
							valp += data.image[IDfilt].array.F[xysize*xysize*step+iip*xysize+ii]*data.image[IDin].array.F[(kk-step)*xysize + ii];
					}
				data.image[IDout].array.F[kk*xysize+iip] = valp;
			
			
				valf = 0.0;
				if(kk+PFlag+1<nbspl)
					valf = (1.0-alpha) * data.image[IDin].array.F[(kk+PFlagl)*xysize+iip] + alpha * data.image[IDin].array.F[(kk+PFlagl+1)*xysize+iip];
				data.image[IDoutf].array.F[kk*xysize+iip] = valf;
			}
	}
	
	
	return(IDout);
}



//
// IDin_name is a 2 or 3D image, open-loop disturbance
// last axis is time (step)
// this optimization asssumes no correlation in noise
//
float LINARFILTERPRED_ScanGain(char* IDin_name, float multfact, float framelag)
{
	float gain;
	float gainmax = 1.1;
	float residual;
	float optgainblock;
	float residualblock;
	float residualblock0;
	float gainstep = 0.01;
	long IDin;
	
	long nbstep;
	long step, step0, step1;
	
	long framelag0;
	long framelag1;
	float alpha;
	
	float *actval_array; // actuator value
	float actval;
	
	long nbvar;
	long axis, naxis;
	
	double *errval;
	double errvaltot;
	long cnt;
	
	FILE *fp;
	char fname[200];
	float mval;
	long ii;
	float tmpv;
	
	int TEST = 0;
	float TESTperiod = 20.0;
	
	// results
	float *optgain;
	float *optres;
	float *res0;
	int optinit = 0;
	
	
	if(framelag<1.00000001)
		{
			printf("ERROR: framelag should be be > 1\n");
			exit(0);
		}
	
	IDin = image_ID(IDin_name);
	naxis = data.image[IDin].md[0].naxis;

	nbvar = 1;
	for(axis=0;axis<naxis-1;axis++)
		nbvar *= data.image[IDin].md[0].size[axis];
	errval = (double*) malloc(sizeof(double)*nbvar);
	
	nbstep = data.image[IDin].md[0].size[naxis-1];

	framelag0 = (long) framelag;
	framelag1 = framelag0+1;
	alpha = framelag-framelag0;

	printf("alpha = %f    nbvar = %ld\n", alpha, nbvar);
	
	list_image_ID();
	if(TEST==1)
		{
			for(ii=0;ii<nbvar;ii++)
			for(step=0;step<nbstep;step++)
				data.image[IDin].array.F[step*nbvar+ii] = 1.0*sin(2.0*M_PI*step/TESTperiod);
		}
	
	
	actval_array = (float*) malloc(sizeof(float)*nbstep);
	
	
	optgain = (float*) malloc(sizeof(float)*nbvar);
	optres = (float*) malloc(sizeof(float)*nbvar);
	res0 = (float*) malloc(sizeof(float)*nbvar);
	
	sprintf(fname, "gainscan.txt");
	
	gain = 0.2;
	ii = 0;
	fp = fopen(fname, "w");
	residualblock = 1.0e20;
	optgainblock = 0.0;
	for(gain=0;gain<gainmax;gain+=gainstep)
		{
			fprintf(fp, "%5.3f", gain);
		
			errvaltot = 0.0;
			for(ii=0;ii<nbvar;ii++)
			{
				errval[ii] = 0.0;
				cnt = 0.0;
				for(step=0;step<framelag1+2;step++)
					actval_array[step] = 0.0;
				for(step=framelag1; step<nbstep; step++)
					{
						step0 = step - framelag0;
						step1 = step - framelag1;

						actval = (1.0-alpha)*actval_array[step0] + alpha*actval_array[step1];
						mval = ((1.0-alpha)*data.image[IDin].array.F[step0*nbvar+ii] + alpha*data.image[IDin].array.F[step1*nbvar+ii]) - actval;
						actval_array[step] = multfact*(actval_array[step-1] + gain * mval);						
						tmpv = data.image[IDin].array.F[step*nbvar+ii] - actval_array[step];
						errval[ii] += tmpv*tmpv;
						cnt++;		
					}
				errval[ii] = sqrt(errval[ii]/cnt);
				fprintf(fp, " %10f", errval[ii]);
				errvaltot += errval[ii]*errval[ii];
			
				if(optinit==0)
				{
					optgain[ii] = gain;
					optres[ii] = errval[ii];
					res0[ii] = errval[ii];
				}
				else
				{
					if(errval[ii]<optres[ii])
						{
							optres[ii] = errval[ii];
							optgain[ii] = gain;
						}
				}
			}
			
			if(optinit==0)
				residualblock0 = errvaltot;
			
			optinit = 1;
			fprintf(fp, "%10f\n", errvaltot);	
			
			if(errvaltot < residualblock)
			{
				residualblock = errvaltot;
				optgainblock = gain;
			}
			
		}
	fclose(fp);
	
	free(actval_array);
	free(errval);
	
	for(ii=0;ii<nbvar;ii++)
		printf("MODE %4ld    optimal gain = %5.2f     residual = %.6f -> %.6f \n", ii, optgain[ii], res0[ii], optres[ii]);
	
	printf("\noptimal block gain = %f     residual = %.6f -> %.6f\n\n", optgainblock, sqrt(residualblock0), sqrt(residualblock));
	printf("RMS per mode = %f -> %f\n", sqrt(residualblock0/nbvar), sqrt(residualblock/nbvar));
	
	free(optgain);
	free(optres);
	free(res0);
	
	return(optgainblock);
}
















