#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



#include <sys/stat.h>

#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_tools/COREMOD_tools.h"

#include "info/info.h"
#include "fft/fft.h"
#include "WFpropagate/WFpropagate.h"
#include "ZernikePolyn/ZernikePolyn.h"
#include "image_gen/image_gen.h"
#include "coronagraphs/coronagraphs.h"
#include "statistic/statistic.h"

#include "coronagraphs/coronagraphs.h"


#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>



#define SWAP(x,y)  tmp=(x);x=(y);y=tmp;
#define PI 3.14159265358979323846264338328



extern DATA data;


#define CORONAGRAPHSDATALOCAL "/CfitsDev/data/coronagraphs"

static int useDFT = 1;
static double DFTZFACTOR = 8.0; // zoom factor for DFTs of focal plane masks



#define CORONAGRAPHS_PSCALE 0.006 //0.05 /* pupil scale [m per pixel] */



#define CORONAGRAPHS_TDIAM 2.4 /* telescope diameter, in meter */
#define CORONAGRAPHS_LAMBDA 0.0000005 /* m */


double CORONAGRAPHS_PIXSCALE = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/CORONAGRAPHS_ARRAYSIZE; // l/d per pixel
double COROTMP1;
double MASKSIZELD;

/* this module can simulate various coronagraphs 

1 : AIC
2 : phase mask
3 : 4 quadrants
4 : PIAA/CPA hybrid
5 : PIAAC/CPA hybrid (fpmask = 5.5)
6 : 8 order BL, circular
7 : 8 order BL, linear
8 : AIC/PIAAC combination
9 : CPA
10 : APLC
11 : APLC, 2 step
12 : APLC, 3 step
13 : APLC, 4 step
14 : APLC, 5 step
15 : APLC, 0 step (= CPA)
16 : ODC
19 : STRIPC
20 : SUMXY
21 : PIAAC/CPA hybrid (fpmask = 5.0)

phase / amplitude aberrations, if any, are in "corpha" and "coramp" images
*/

static int FPMASKSIZE_ERROR = 0; /* 1 if focal plane mask error is included (for chromaticity for example) */
static double FPMASK_FACTOR = 1.0;  // if focal plane mask error included, multiplicative coefficient to the focal plane mask size - currently only for MULTI APLC
static double FPMASK_FACTOR1 = 1.0;
static double FPMASK_FACTOR2 = 1.1;
static double FPM_TRANSM1 = 0.0;
static double FPM_TRANSM2 = 1.0;

/* arrays and numbers that need to be in memory to use various coronagraphs */
/* 8 order BL Lyot */
#define BL8MASK_NBSTEP 100000 /* number of points over which the BL mask is computed */
#define BL8MASK_STEP 0.01  /* lambda/d per step */
static double BL8MASK[BL8MASK_NBSTEP];
static double BL8MASK_m = 1.0;
static double BL8MASK_l = 3.0;
static double BL8MASK_eps = 0.6;
static int BL8MODE = 1; /* 0 = circular, 1 = linear */

/* 4 order BL Lyot */
static double BL4MASK_eps = 0.21;

/* Optical Differenciation Coronagraph */
static double ODC_GAUSS = 10.0; /* in lambda/d */
static double ODCMASK_eps = 0.85;

/* SHEAR4 */
static double SHEAR4_SHEAR = 0.1;

/* R&R Phase mask */
static double RRPM_RADIUS = 0.53773; /* phase mask radius in lambda/d */
static double RRPM_P2 = -0.3207; /* term of order 2 in polynomial pupil apodization */
static double RRPM_P3 = 0.00952; /* term of order 3 in polynomial pupil apodization */
static double RRPM_P4 = 0.000658; /* term of order 4 in polynomial pupil apodization */
static double RRPM_P5 = 0.001127; /* term of order 5 in polynomial pupil apodization */
static int RRPM_PIAA = 1;

/* CPA */
#define CPAAPO_NBPOINTS 1930
static double CPAAPO[CPAAPO_NBPOINTS];
#define CPAAPO_FNAME "pup_10.411256_2.0"
#define CPAFPMASKRAD 0.0 /*4.4*/ /* in l/d */
#define CPAPPMASKRAD 1.0 /* in pupil unit */

/* PIAA and PIAAC (both hybrids) */
static double PIAACENTOBS = 0.0;
static long PIAAAPO_NBPOINTS = 2042;
static double *PIAAAPO;
static double *PIAA_HYBRID_CPAAPO;
static double PIAAAPO2D[CORONAGRAPHS_ARRAYSIZE*CORONAGRAPHS_ARRAYSIZE];
static char PIAAAPO_FNAME[400]="pup_10.411256_2.0";
static char PIAAAPODIZE_2DAPOFNAME[400]="APLCapo_4.200.1024.ref.gz";

/*#define PIAAAPO_FNAME "Apod04.dat.2"*/
static int PIAAFPMASK = 0;
static int PIAALOWFS = 1;
static double PIAAFPMASKRAD = 5.5; /* in l/d */ /* 4.25 is the critical size to maintain 1e-10 contrast at all separations  default=4.6*/
#define PIAAPPMASKRAD 1.0 /* in pupil unit */
#define PIAACPPMASKRAD1 1.0 /* in pupil unit */
#define PIAA_HYBRID_CST 0.00000001 /* relative to peak */
static double *piaaconfpup_amp_profile;
static double *piaaconfr0;
static double *piaaconfr1;
#define piaaconfNPUPFILESIZE 200000
static long piaaconfNBpoints;
#define piaaconfNBr1fr0_pts 20000
static double *piaaconfr0fr1;
static double *piaaconfr1fr0;

// PIAA optics shapes
static double *piaaconfM0;
static double *piaaconfM1;

static int piaaconfdirection = 0;
static double PIAAFLUXFACTOR;
static double PIAAOVERSIZE = 1.025;
static int initPIAA=0;
static int AUTOPIAACMASK = 1;
static double APLC_CentOBS0 = 0.0;
static double APLC_CentOBS1 = 0.0;

static double PIAAextfactor0 = 1.1;
static double PIAAextfactor1 = 1.1;

/* AIC / PIAA hybrid */
#define PIAAAIC_FIELDLIMIT 15.0 /* in l/d */

/* APLC */
static int APLC_PIAA = 0; /* 1 if the apodization is performed by PIAA */
static long NB_APLC_STEP = 0; /* number of APLC steps */
static double APLC_FPMASKsize = 4.2; /* APLC focal plane mask radius in l/d */
static int APLC_FLIP = 0; // output AIC in APLC PIAA
static int APLC_PMASK = 0; // 1 if focal plane mask in APLC has phase
static double FPMASK_transm_error = 0.0;
static double FPMASK_size_error = 0.0;
/* STRIPC */
static double STRIPCOFFSET = 0.4;

/* OVC */
static long OVC_CHARGE = 2;


// RAW apodization
static long aporawN;
static double *aporaw_r; // radius
static double *aporaw_v; // value
static int fitapoINIT = 0;
static long fitapoN = 5; // number of terms of form a exp(b x^c)
static double *fitapo_a;
static double *fitapo_b;
static double *fitapo_c;
static double *fitapo_c1;


static long LOOPCNT = 0;
static double optval0;



static long IDprol_init;
static long IDprol_ffrac;
static long IDprol_transm;
static long IDprol_peak;
static long IDprol_fitapo_a, IDprol_fitapo_b, IDprol_fitapo_c;
static long IDprol_fitfit;

static double fitapo_minc = 0.1;
static double APLCapo_CO_START = 0.0;
static double APLCapo_CO_END = 0.5;
static double APLCapo_CO_STEP = 0.001;
static double APLCapo_FPMRAD_START = 0.8;
static double APLCapo_FPMRAD_END = 5.0;
static double APLCapo_FPMRAD_STEP = 0.001;






// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//

int_fast8_t coronagraph_make_2Dprolate_cli()
{
  if(CLI_checkarg(1,1)+CLI_checkarg(2,1)+CLI_checkarg(3,1)+CLI_checkarg(4,3)+CLI_checkarg(5,2)==0)
    {
      coronagraph_make_2Dprolate(data.cmdargtoken[1].val.numf, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.string, data.cmdargtoken[5].val.numl, "NULLim");
      return 0;
    }
  else
    return 1;
}

int_fast8_t coronagraph_make_2Dprolateld_cli()
{
  if(CLI_checkarg(1,1)+CLI_checkarg(2,1)+CLI_checkarg(3,1)+CLI_checkarg(4,3)+CLI_checkarg(5,2)==0)
    {
      coronagraph_make_2Dprolateld(data.cmdargtoken[1].val.numf, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.string, data.cmdargtoken[5].val.numl, "NULL");
      return 0;
    }
  else
    return 1;
}

int_fast8_t coronagraph_update_2Dprolate_cli()
{
  if(CLI_checkarg(1,1)+CLI_checkarg(2,1)+CLI_checkarg(3,1)+CLI_checkarg(4,1)==0)
    {
      coronagraph_update_2Dprolate(data.cmdargtoken[1].val.numf, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numf);
      return 0;
    }
  else
    return 1;
}


int_fast8_t  coronagraph_simulPSF_cli()
{
  if(CLI_checkarg(1,1)+CLI_checkarg(2,1)+CLI_checkarg(3,3)+CLI_checkarg(4,2)==0)
    {
      coronagraph_simulPSF(data.cmdargtoken[1].val.numf, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.string, data.cmdargtoken[4].val.numl, "");
      return 0;
    }
  else
    return 1;
}


int_fast8_t CORONAGRAPHS_scanPIAACMC_centObs_perf_cli()
{
  if(CLI_checkarg(1,1)==0)
    {
      CORONAGRAPHS_scanPIAACMC_centObs_perf(data.cmdargtoken[1].val.numf);
      return 0;
    }
  else
    return 1;
}




int_fast8_t init_coronagraphs()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "coronagraph routines");
  data.NBmodule++;

  strcpy(data.cmd[data.NBcmd].key,"cormk2Dprolate");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = coronagraph_make_2Dprolate_cli;
  strcpy(data.cmd[data.NBcmd].info,"make 2D prolate");
  strcpy(data.cmd[data.NBcmd].syntax,"<focal plane mask radius> <central obstruction> <image name>");
  strcpy(data.cmd[data.NBcmd].example,"cormk2Dprolate 1.20 0.1 prolim");
  strcpy(data.cmd[data.NBcmd].Ccall,"coronagraph_make_2Dprolate(double masksize, double centralObs, const char *outname)"); 
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"cormk2Dprolateld");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = coronagraph_make_2Dprolateld_cli;
  strcpy(data.cmd[data.NBcmd].info,"make 2D prolate, l/D unit mask size");
  strcpy(data.cmd[data.NBcmd].syntax,"<focal plane mask radius [l/D]> <central obstruction> <image name> <array size>");
  strcpy(data.cmd[data.NBcmd].example,"cormk2Dprolateld 1.20 0.1 prolim 2048");
  strcpy(data.cmd[data.NBcmd].Ccall,"coronagraph_make_2Dprolateld(double masksizeld, double centralObs, const char *outname, long size)"); 
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"corup2Dprolate");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = coronagraph_update_2Dprolate_cli;
  strcpy(data.cmd[data.NBcmd].info,"update 2D prolate");
  strcpy(data.cmd[data.NBcmd].syntax,"<focal plane mask radius [l/D]> <central obstruction> <DFT zoom factor>");
  strcpy(data.cmd[data.NBcmd].example,"cormk2Dprolateld 1.20 0.1 4.0");
  strcpy(data.cmd[data.NBcmd].Ccall,"coronagraph_update_2Dprolate(double masksizeld, double centralObs, double zfactor)"); 
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"corprolcomp");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = coronagraph_APLCapo_compile;
  strcpy(data.cmd[data.NBcmd].info,"compile prolate functions data");
  strcpy(data.cmd[data.NBcmd].syntax,"no argument");
  strcpy(data.cmd[data.NBcmd].example,"corprolcomp");
  strcpy(data.cmd[data.NBcmd].Ccall,"int coronagraph_APLCapo_compile()"); 
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"corsimpsf");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = coronagraph_simulPSF_cli;
  strcpy(data.cmd[data.NBcmd].info,"create coronagraph PSF");
  strcpy(data.cmd[data.NBcmd].syntax,"<x [l/D]> <y [l/D]> <output image name> <coronagraph index>");
  strcpy(data.cmd[data.NBcmd].example,"corsimpsf 1.0 0.1 50");
  strcpy(data.cmd[data.NBcmd].Ccall,"coronagraph_simulPSF(double xld, double yld, const char *psfname, long coronagraph_type, const char *options)"); 
  data.NBcmd++;
  
  
  strcpy(data.cmd[data.NBcmd].key,"coroscanpiaacmcperf");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = CORONAGRAPHS_scanPIAACMC_centObs_perf_cli;
  strcpy(data.cmd[data.NBcmd].info,"scan PIAACMC perf");
  strcpy(data.cmd[data.NBcmd].syntax,"<input central obstruction>");
  strcpy(data.cmd[data.NBcmd].example,"coroscanpiaacmcperf 0.2");
  strcpy(data.cmd[data.NBcmd].Ccall,"int CORONAGRAPHS_scanPIAACMC_centObs_perf( double obs0input )"); 
  data.NBcmd++;
  

 // add atexit functions here

  return 0;

}






//
// if pupmask does not exist, do not use it
//
double coronagraph_make_2Dprolate(double fpmradpix, double beamradpix, double centralObs, const char *outname, long size, const char *pupmask_name)
{
    FILE *fp;
    long size2;
    long IDfpmask;
    long IDpupa0,IDpupp0,IDpupa0m,IDpupp0m;
    long ii,jj,ii1;
    double peak;
    long IDprol,ID;
    long double total,total0,total2,v1;
    long iter;
    long NBiter = 40;
    double transm;
    double x,y,r2;
    long IDr,IDi,IDrcp,IDicp; //,ID1,ID2;
    long IDprolr,IDproli,IDprolp;
    char fname[200];
	long IDpupmask;

    int CentralObstructionFlag = 0;
    long IDpupa0co;

    size2 = size*size;

    if(centralObs > 0.001)
        CentralObstructionFlag = 1;


    ID = variable_ID("PNBITER");
    if(ID!=-1)
    {
        NBiter = (long) (1.0*data.variable[ID].value.f+0.01);
    }

    v1 = 0.0;



    MASKSIZELD = (2.0*beamradpix/size)*fpmradpix;
    printf("MASK RADIUS = %f pix = %f l/D\n", fpmradpix, MASKSIZELD);
    printf("PUPIL RADIUS = %f\n", beamradpix);
    printf("SIZE = %ld\n", size);

    IDproli = create_2Dimage_ID("proli", size, size);
    IDprolr = create_2Dimage_ID("prolr", size, size);
    IDprolp = create_2Dimage_ID("prolp", size, size);


    IDfpmask = make_subpixdisk("FPmask", size, size, size/2, size/2, fpmradpix);
    for(ii=0; ii<size2; ii++)
    {
        data.image[IDfpmask].array.F[ii] = 1.0*data.image[IDfpmask].array.F[ii];
        //      if(data.image[IDfpmask].array.F[ii]<0.01)
        //	data.image[IDfpmask].array.F[ii] = 0.0;
    }

    //  save_fl_fits("FPmask","!FPmask.tmp.fits");
    // exit(0);

    IDpupp0 = create_2Dimage_ID("pupp0", size, size);

    IDpupa0 = image_ID("pupa0");
    if(IDpupa0!=-1)
    {
        if(data.image[IDpupa0].md[0].size[0]!=size)
        {
            printf("ERROR: pupa0 should be %ld x %ld\n", size, size);
            exit(0);
        }
    }
    
    if(IDpupa0==-1)
    {
        IDpupa0 = make_disk("pupa0", size, size, size/2, size/2, beamradpix);
        if(CentralObstructionFlag == 1)
        {
            IDpupa0co = make_disk("pupa0co", size, size, size/2, size/2, beamradpix*centralObs);
            for(ii=0; ii<size*size; ii++)
                data.image[IDpupa0].array.F[ii] -= data.image[IDpupa0co].array.F[ii];
            delete_image_ID("pupa0co");
        }
        
        IDpupmask = image_ID(pupmask_name);
        if(IDpupmask != -1)
			{
				for(ii=0; ii<size*size; ii++)
					data.image[IDpupa0].array.F[ii] *= data.image[IDpupmask].array.F[ii];
			}
        save_fits("pupa0", "!test_pupa0.fits");
    }
    else
		save_fits("pupa0", "!test__pupa0.fits");


    total0 = 0.0;
    for(ii=0; ii<size2; ii++)
    {
        v1 = data.image[IDpupa0].array.F[ii];
        total0 += v1*v1;
    }


    //
    // INITIALIZE PROLATE
    //
    // IDprol : prolate
    //

    // total0 = arith_image_total("pupa0");
    IDpupa0m = create_2Dimage_ID("pupa0m",size,size);
    IDpupp0m = create_2Dimage_ID("pupp0m",size,size);

    IDprol = create_2Dimage_ID(outname,size,size);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            if(data.image[IDpupa0].array.F[jj*size+ii] > 0.01)
            {
                x = 1.0*ii-size/2;
                y = 1.0*jj-size/2;
                r2 = (x*x+y*y)/beamradpix/beamradpix;
                data.image[IDprolr].array.F[jj*size+ii] = exp(-r2*0.73/5.0*fpmradpix);
                data.image[IDproli].array.F[jj*size+ii] = 0.0;
                data.image[IDprol].array.F[jj*size+ii] = data.image[IDprolr].array.F[jj*size+ii];
            }
            else
            {
                data.image[IDprolr].array.F[jj*size+ii] = 0.0;
                data.image[IDproli].array.F[jj*size+ii] = 0.0;
                data.image[IDprol].array.F[jj*size+ii] = 0.0;
            }
        }


    ID = image_ID("apostart");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDprolr].array.F[ii] = data.image[ID].array.F[ii];
            data.image[IDprol].array.F[ii] = data.image[ID].array.F[ii];
        }
    }
    delete_image_ID("apostart");



    //
    // START LOOP
    //

    for(iter=0; iter<NBiter; iter++)
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDpupa0m].array.F[ii] =  data.image[IDpupa0].array.F[ii]*data.image[IDprol].array.F[ii];
            data.image[IDpupp0m].array.F[ii] =  data.image[IDpupp0].array.F[ii]+0.0*data.image[IDprolp].array.F[ii];
        }

        mk_complex_from_amph("pupa0m", "pupp0", "pc1", 0);
        permut("pc1");
        do2dfft("pc1","fc1");
        permut("fc1");
        mk_amph_from_complex("fc1", "fa1", "fp1", 0);

        execute_arith("fa1m=fa1*FPmask");

        total = 0.0;
        for(ii=0; ii<size2; ii++)
            total += data.image[IDprol].array.F[ii]*data.image[IDprol].array.F[ii];
        transm = total/total0;

        ID = image_ID("fa1m");
        total = 0.0;
        for(ii=0; ii<size2; ii++)
            total += data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
        total /= 1.0*size2*total0;

        ID = image_ID("fa1");
        total2 = 0.0;
        for(ii=0; ii<size2; ii++)
            total2 += data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
        total2 /= 1.0*size2*total0;

        printf("%ld/%ld   FPmask throughput = %g    prolate throughput = %.18g\n",iter,NBiter,(double) (total/total2),(double) transm);

        ID = image_ID("fp1");
        //      for(ii=0;ii<size2;ii++)
        //	data.image[ID].array.F[ii] = 0.0;

        mk_complex_from_amph("fa1m", "fp1", "fc2", 0);
        permut("fc2");
        do2dfft("fc2","pc2");
        permut("pc2");



        mk_reim_from_complex("pc2", "pr2", "pi2", 0);

        IDr = image_ID("pr2");
        copy_image_ID("pr2", "pr2cp", 0);
        IDrcp = image_ID("pr2cp");
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
                data.image[IDr].array.F[jj*size+ii] = 0.0;

        for(ii=1; ii<size; ii++)
            for(jj=1; jj<size; jj++)
                data.image[IDr].array.F[jj*size+ii] = data.image[IDrcp].array.F[(size-jj)*size+(size-ii)]/size2;
        peak = data.image[IDr].array.F[(size/2)*size+size/2];
        delete_image_ID("pr2cp");

        IDi = image_ID("pi2");
        copy_image_ID("pi2", "pi2cp", 0);
        IDicp = image_ID("pi2cp");
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
                data.image[IDi].array.F[jj*size+ii] = 0.0;

        for(ii=1; ii<size; ii++)
            for(jj=1; jj<size; jj++)
                data.image[IDi].array.F[jj*size+ii] = data.image[IDicp].array.F[(size-jj)*size+(size-ii)]/size2;

        printf("PEAK =  %f ", peak);
        peak = 0.0;
        for(ii=0; ii<size2; ii++)
            if(data.image[IDpupa0].array.F[ii] > 0.01)
                if(fabs(data.image[IDr].array.F[ii])>fabs(peak))
                    peak = data.image[IDr].array.F[ii];
        printf(" %f\n", peak);



        // write info file
        // col 1 is mask size in SYSTEM l/D
        // col 2 is fraction of light going through FPmask
        // col 3 is prolate throughput if not done by PIAA
        // col 4 is peak
        // ideal complex transmission for FPM is  MASKCAMULT = -(1.0-peak)/peak;
        printf("------- MASK SIZE = %f l/D ------\n", MASKSIZELD);
        sprintf(fname, "%s/APLCapo.%.3f.%.3f.info", data.SAVEDIR, MASKSIZELD, centralObs);
        fp = fopen(fname, "w");
        fprintf(fp,"%f %.18f %.18f %.18f %.18f\n", MASKSIZELD, (double) (total/total2), (double) transm, peak, centralObs);
        fclose(fp);
        printf("IDEAL COMPLEX TRANSMISSION = %g\n", -(1.0-peak)/peak);

        //      save_fl_fits("pr2","!pr2.tmp.fits");
        // save_fl_fits("pi2","!pi2.tmp.fits");

        //      printf("peak = %f\n",peak);


        for(ii=0; ii<size2; ii++)
            if(data.image[IDpupa0].array.F[ii] > 0.01)
            {
                data.image[IDprolr].array.F[ii] = 1.0*data.image[IDr].array.F[ii];
                data.image[IDproli].array.F[ii] = 1.0*data.image[IDi].array.F[ii];
            }

        //      for(ii=0;ii<size2;ii++)
        //	data.image[IDprol].array.F[ii] = data.image[IDprolr].array.F[ii];
        for(ii=1; ii<size-1; ii++)
            for(jj=1; jj<size-1; jj++)
            {
                ii1 = jj*size+ii;
                data.image[IDprol].array.F[jj*size+ii] = sqrt(data.image[IDprolr].array.F[ii1]*data.image[IDprolr].array.F[ii1]+data.image[IDproli].array.F[ii1]*data.image[IDproli].array.F[ii1]);
                data.image[IDprolp].array.F[jj*size+ii] = -atan2(data.image[IDproli].array.F[ii1],data.image[IDprolr].array.F[ii1]);
            }

        peak = 0.0;
        for(ii=0; ii<size2; ii++)
            if(data.image[IDpupa0].array.F[ii] > 0.01)
                if(fabs(data.image[IDprol].array.F[ii])>fabs(peak))
                    peak = data.image[IDprol].array.F[ii];

        //      printf("peak = %f\n",peak);
        for(ii=0; ii<size2; ii++)
        {
            if(data.image[IDpupa0].array.F[ii] > 0.001)
                data.image[IDprol].array.F[ii] /= peak;
            else
                data.image[IDprol].array.F[ii] = 0.0;
        }


        if (0) // TEST
        {
            mk_amph_from_complex("pc1", "pc1_a", "pc1_p", 0);
            save_fl_fits("pc1_a", "!pc1_a.1.fits");
            save_fl_fits("pc1_p", "!pc1_p.1.fits");
            save_fl_fits("pr2", "!pr2.1.fits");
            save_fl_fits("pi2", "!pi2.1.fits");
            delete_image_ID("pc1_a");
            delete_image_ID("pc1_p");
            save_fl_fits(outname, "!proltmp.1.fits");
        }


        delete_image_ID("pc1");
        delete_image_ID("fc1");
        delete_image_ID("fa1");
        delete_image_ID("fp1");
        delete_image_ID("fa1m");
        delete_image_ID("fc2");
        delete_image_ID("pc2");
        delete_image_ID("pr2");
        delete_image_ID("pi2");
    }

    for(ii=0; ii<size2; ii++)
        data.image[IDfpmask].array.F[ii] = 1.0-data.image[IDfpmask].array.F[ii]*(1.0+(1.0-peak)/peak);


    //  save_fl_fits("FPmask","!FPmask.fits");

    delete_image_ID("FPmask");
    delete_image_ID("pupp0");
    delete_image_ID("pupa0");
    delete_image_ID("pupa0m");

    delete_image_ID("proli");
    delete_image_ID("prolr");
    delete_image_ID("prolp");
    delete_image_ID("pupp0m");


    return(transm);
}



//
// High definition prolate computation
// Uses DFT
//
double coronagraph_make_2Dprolate_DFT(double fpmradpix, double beamradpix, double centralObs, const char *outname, long size, const char *pupmask_name)
{
    FILE *fp;
    long size2;
    long IDfpmask, IDfpmaskc, IDfpmaskz;
    long IDpupa0,IDpupp0,IDpupa0m,IDpupp0m;
    long ii,jj,ii1;
    double peak;
    long IDprol,ID;
    long double total,total0,total2,v1;
    long iter;
    long NBiter = 10;
    double transm;
    double x,y,r2;
    long IDr,IDi,IDrcp,IDicp; //,ID1,ID2;
    long IDprolr,IDproli,IDprolp;
    char fname[200];
    char vname[200];
	long IDpupmask;
	
    int CentralObstructionFlag = 0;
    long IDpupa0co;

    double re, im;

    // focal plane mask shape perturbations
    long k;
    long fpmshape_n = 20;
    double *fpmshape_ra;
    double *fpmshape_ka;
    double *fpmshape_pa;

    long IDv;


    size2 = size * size;

    fpmshape_ra = (double*) malloc(sizeof(double)*fpmshape_n);
    fpmshape_ka = (double*) malloc(sizeof(double)*fpmshape_n);
    fpmshape_pa = (double*) malloc(sizeof(double)*fpmshape_n);


    if((ID=variable_ID("DFTZFACTOR"))!=-1)
        DFTZFACTOR = data.variable[ID].value.f;

    if(centralObs > 0.001)
        CentralObstructionFlag = 1;

    // SET NUMBER OF ITERATION TO VALUE OTHER THAN DEFAULT IF VARIABLE PNBITER EXISTS
    ID = variable_ID("PNBITER");
    if(ID!=-1)
        NBiter = (long) (1.0*data.variable[ID].value.f+0.01);


    v1 = 0.0;
    MASKSIZELD = (2.0*beamradpix/size) * fpmradpix;
    printf("MASK SIZE = %f pixel = %f l/D\n", fpmradpix, MASKSIZELD);
    printf("PUPIL RADIUS = %f\n", beamradpix);
    printf("SIZE = %ld\n", size);

    IDproli = create_2Dimage_ID("proli", size, size);
    IDprolr = create_2Dimage_ID("prolr", size, size);
    IDprolp = create_2Dimage_ID("prolp", size, size);


    // CREATE FOCAL PLANE MASK

    fpmshape_ra[0] = 0.0;
    fpmshape_ka[0] = 3.0;
    fpmshape_pa[0] = M_PI/2.0;

    for(k=0; k<fpmshape_n; k++)
    {
        fpmshape_ra[k] = 0.0;
        fpmshape_ka[k] = 0.0;
        fpmshape_pa[k] = 0.0;

        sprintf(vname, "FPMSHAPE_%ld_r", k);
        if((ID=variable_ID(vname))!=-1)
            fpmshape_ra[k] = data.variable[ID].value.f;


        sprintf(vname, "FPMSHAPE_%ld_k", k);
        if((ID=variable_ID(vname))!=-1)
            fpmshape_ka[k] = data.variable[ID].value.f;

        sprintf(vname, "FPMSHAPE_%ld_p", k);
        if((ID=variable_ID(vname))!=-1)
            fpmshape_pa[k] = data.variable[ID].value.f;

        printf("FPM SHAPE TERM %02ld : %f %f %f\n", k,  fpmshape_ra[k],  fpmshape_ka[k],  fpmshape_pa[k]);
    }

    IDfpmask = make_subpixdisk("FPmask", size, size, size/2, size/2, fpmradpix);
    IDfpmaskz = make_subpixdisk_perturb("FPmaskz", size, size, size/2, size/2, fpmradpix*DFTZFACTOR, fpmshape_n, fpmshape_ra, fpmshape_ka, fpmshape_pa);
    IDfpmaskc = create_2DCimage_ID("_fpmz", size, size);
    for(ii=0; ii<size2; ii++)
    {
        data.image[IDfpmask].array.F[ii] = 1.0*data.image[IDfpmaskz].array.F[ii];
        data.image[IDfpmaskc].array.CF[ii].re = data.image[IDfpmaskz].array.F[ii];
        data.image[IDfpmaskc].array.CF[ii].im = 0.0;
        //      if(data.image[IDfpmask].array.F[ii]<0.01)
        //	data.image[IDfpmask].array.F[ii] = 0.0;
    }

    sprintf(fname, "!%s/FPmask.tmp.fits", data.SAVEDIR);
    save_fl_fits("FPmask",fname);

    free(fpmshape_ra);
    free(fpmshape_ka);
    free(fpmshape_pa);

    //  exit(0);

    if(0)
    {
        mk_amph_from_complex("_fpmz", "_fpmza", "_fpmzp", 0);
        save_fl_fits("_fpmza", "!_fpmza.fits");
        save_fl_fits("_fpmzp", "!_fpmzp.fits");
        delete_image_ID("_fpmza");
        delete_image_ID("_fpmzp");
    }
    // exit(0);

    IDpupp0 = create_2Dimage_ID("pupp0", size, size);

    IDpupa0 = image_ID("pupa0");
    if(IDpupa0!=-1)
    {
        if(data.image[IDpupa0].md[0].size[0]!=size)
        {
            printf("ERROR: pupa0 should be %ld x %ld\n",size,size);
            exit(0);
        }
    }
    if(IDpupa0==-1)
    {
        //      printf("Creating pupil, radius = %f\n", beamradpix);
        IDpupa0 = make_disk("pupa0", size, size, size/2, size/2, beamradpix);
        
        if(CentralObstructionFlag == 1)
        {
            //  printf("Central obstruction = %f\n", beamradpix*centralObs);
            IDpupa0co = make_disk("pupa0co", size, size, size/2, size/2, beamradpix*centralObs);
            for(ii=0; ii<size*size; ii++)
                data.image[IDpupa0].array.F[ii] -= data.image[IDpupa0co].array.F[ii];
            delete_image_ID("pupa0co");
        }
        
     //   save_fits("pupa0", "!test_pupa00.fits"); //TEST
        IDpupmask = image_ID(pupmask_name);
        if(IDpupmask != -1)
			{
				for(ii=0; ii<size*size; ii++)
					data.image[IDpupa0].array.F[ii] *= data.image[IDpupmask].array.F[ii];
			}
      //  save_fits("pupa0", "!test_pupa0.fits"); //TEST
    //    save_fits(pupmask_name, "!test_pupamask.fits");//TEST
	//printf("Test wait\n");
	//list_image_ID();
	//	sleep(100.0); // TEST
    }
	//else
	//	 save_fits("pupa0", "!test__pupa0.fits");//TEST
		
    //  save_fl_fits("pupa0", "!pupa00.fits");

    total0 = 0.0;
    for(ii=0; ii<size2; ii++)
    {
        v1 = data.image[IDpupa0].array.F[ii];
        total0 += v1*v1;
    }


    //
    // INITIALIZE PROLATE
    //
    // IDprol : prolate
    //

    // total0 = arith_image_total("pupa0");
    IDpupa0m = create_2Dimage_ID("pupa0m",size,size);
    IDpupp0m = create_2Dimage_ID("pupp0m",size,size);

    IDprol = create_2Dimage_ID(outname, size, size);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            if(data.image[IDpupa0].array.F[jj*size+ii] > 0.0001)
            {
                x = 1.0*ii-size/2;
                y = 1.0*jj-size/2;
                r2 = (x*x+y*y)/beamradpix/beamradpix;
                data.image[IDprolr].array.F[jj*size+ii] = exp(-r2*0.73/5.0*fpmradpix);
                data.image[IDproli].array.F[jj*size+ii] = 0.0;
                data.image[IDprol].array.F[jj*size+ii] = data.image[IDprolr].array.F[jj*size+ii];
            }
            else
            {
                data.image[IDprolr].array.F[jj*size+ii] = 0.0;
                data.image[IDproli].array.F[jj*size+ii] = 0.0;
                data.image[IDprol].array.F[jj*size+ii] = 0.0;
            }
        }


    ID = image_ID("apostart");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDprolr].array.F[ii] = data.image[ID].array.F[ii];
            data.image[IDprol].array.F[ii] = data.image[ID].array.F[ii];
        }
    }
    delete_image_ID("apostart");



    //
    // START LOOP
    //

    for(iter=0; iter<NBiter; iter++)
    {
        printf("[%ld/%ld - %ld] ", iter, NBiter, size2);
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDpupa0m].array.F[ii] =  data.image[IDpupa0].array.F[ii]*data.image[IDprol].array.F[ii];
            data.image[IDpupp0m].array.F[ii] =  data.image[IDpupp0].array.F[ii]+0.0*data.image[IDprolp].array.F[ii];
        }

        mk_complex_from_amph("pupa0m", "pupp0", "pc1", 0);
	//	save_fits("pupa0m", "!test_pupa0m.fits");//TEST
	//	sleep(100.0);//TEST

        //
        // insert focal plane mask
        // PUPIL -> PUPIL
        // pc1 -> pc3
        //
        // save_fl_fits("pupa0", "!pupa0.fits");
        //save_fl_fits("pupa0m", "!pupa0m.fits");
        //save_fl_fits("pupp0", "!pupp0.fits");
        //printf("DFTZFACTOR = %f\n", DFTZFACTOR);

        fft_DFTinsertFPM( "pc1", "_fpmz", DFTZFACTOR, "pc3");


        //      exit(0);


        // MEASURE OVERALL PROLATE TRANSMISSION -> transm
        total = 0.0;
        for(ii=0; ii<size2; ii++)
            total += data.image[IDprol].array.F[ii]*data.image[IDprol].array.F[ii];
        transm = total/total0;

        // MEASURE TOTAL FRACTION OF LIGHT WITHIN MASK -> total
        total = data.DOUBLEARRAY[0];
        total /= 1.0*size2*total0;


        // TOTAL LIGHT IN IMAGE -> total2
        ID = image_ID("pc1");
        total2 = 0.0;
        for(ii=0; ii<size2; ii++)
        {
            data.image[ID].array.CF[ii].re *= size;
            data.image[ID].array.CF[ii].im *= size;
            re = data.image[ID].array.CF[ii].re;
            im = data.image[ID].array.CF[ii].im;
            total2 += re*re + im*im;
        }
        total2 /= 1.0*size2*total0;



        printf("%ld/%ld   FPmask throughput = %g    prolate throughput = %.18g\n", iter, NBiter, (double) (total/total2), (double) transm);

        mk_reim_from_complex("pc3", "pr2", "pi2", 0);



        IDr = image_ID("pr2");
        IDi = image_ID("pi2");





        peak = 0.0;
        for(ii=0; ii<size2; ii++)
            if(data.image[IDpupa0].array.F[ii] > 0.0001)
                if(fabs(data.image[IDr].array.F[ii])>fabs(peak))
                    peak = data.image[IDr].array.F[ii];
        printf("%f\n", peak);



        // write info file
        // col 1 is mask size in SYSTEM l/D
        // col 2 is fraction of light going through FPmask
        // col 3 is prolate throughput if not done by PIAA
        // col 4 is peak
        // ideal complex transmission for FPM is  MASKCAMULT = -(1.0-peak)/peak;
        printf("------- MASK SIZE = %f l/D ------\n", MASKSIZELD);
        sprintf(fname, "%s/APLCapo.%.3f.%.3f.info", data.SAVEDIR, MASKSIZELD, centralObs);
        fp = fopen(fname, "w");
        fprintf(fp,"%f %.18f %.18f %.18f\n", MASKSIZELD, (double) (total/total2), (double) transm, peak);
        fclose(fp);
        printf("IDEAL COMPLEX TRANSMISSION = %g\n",-(1.0-peak)/peak);

        IDv = create_variable_ID("APLCmaskCtransm", -(1.0-peak)/peak);

        if(0) // TEST
        {
            save_fl_fits("pr2","!_prol_pr2.fits");
            save_fl_fits("pi2","!_prol_pi2.fits");
        }

        for(ii=0; ii<size2; ii++)
            if(data.image[IDpupa0].array.F[ii] > 0.0001)
            {
                data.image[IDprolr].array.F[ii] = 1.0*data.image[IDr].array.F[ii];
                data.image[IDproli].array.F[ii] = 1.0*data.image[IDi].array.F[ii];
            }

        for(ii=0; ii<size2; ii++)
            data.image[IDprol].array.F[ii] = sqrt(data.image[IDprolr].array.F[ii]*data.image[IDprolr].array.F[ii] + data.image[IDproli].array.F[ii]*data.image[IDproli].array.F[ii]);


        peak = 0.0;
        for(ii=0; ii<size2; ii++)
            if(data.image[IDpupa0].array.F[ii] > 0.0001)
                if(fabs(data.image[IDprol].array.F[ii])>fabs(peak))
                    peak = data.image[IDprol].array.F[ii];


        for(ii=0; ii<size2; ii++)
        {
            if(data.image[IDpupa0].array.F[ii] > 0.001)
                data.image[IDprol].array.F[ii] /= peak;
            else
                data.image[IDprol].array.F[ii] = 0.0;
        }


        if (0) // TEST
        {
            mk_amph_from_complex("pc1", "pc1_a", "pc1_p", 0);
            save_fl_fits("pc1_a", "!pc1_a.fits");
            save_fl_fits("pc1_p", "!pc1_p.fits");
            mk_amph_from_complex("pc3", "pc3_a", "pc3_p", 0);
            save_fl_fits("pc3_a", "!pc3_a.fits");
            save_fl_fits("pc3_p", "!pc3_p.fits");
            save_fl_fits("pr2", "!pr2.fits");
            save_fl_fits("pi2", "!pi2.fits");
            delete_image_ID("pc1_a");
            delete_image_ID("pc1_p");
            delete_image_ID("pc3_a");
            delete_image_ID("pc3_p");
            save_fl_fits(outname, "!proltmp.fits");
            //	  exit(0);
        }



        delete_image_ID("pc1");
        delete_image_ID("pc3");
        delete_image_ID("pr2");
        delete_image_ID("pi2");

    }

    for(ii=0; ii<size2; ii++)
        data.image[IDfpmask].array.F[ii] = 1.0-data.image[IDfpmask].array.F[ii]*(1.0+(1.0-peak)/peak);

    for(ii=0; ii<size2; ii++)
        data.image[IDfpmaskz].array.F[ii] = 1.0-data.image[IDfpmaskz].array.F[ii]*(1.0+(1.0-peak)/peak);

    //  save_fl_fits("FPmaskz","!FPmask.fits");
    delete_image_ID("FPmaskz");

    delete_image_ID("FPmask");
    delete_image_ID("_fpmz");
    delete_image_ID("pupp0");
    delete_image_ID("pupa0");
    delete_image_ID("pupa0m");
    delete_image_ID("proli");
    delete_image_ID("prolr");
    delete_image_ID("prolp");
    delete_image_ID("pupp0m");

    return(transm);
}








int coronagraph_make_2Dprolateld(double masksizeld, double beamradpix, double centralObs, const char *outname, long size, const char *pupmask_name)
{
    double fpmradpix;

  //  printf("Making 2D prolate (DFT=%d)...  \n", useDFT);
  //  sleep(5.0);


    fpmradpix = masksizeld/(2.0*beamradpix/size);

    if(useDFT==0)
        coronagraph_make_2Dprolate(fpmradpix, beamradpix, centralObs, outname,  size, pupmask_name);
    else
        coronagraph_make_2Dprolate_DFT(fpmradpix, beamradpix, centralObs, outname, size, pupmask_name);

    return (0);
}



// works for non circular pupil shape
//
// input is pupaCS
int coronagraph_make_2Dprolate_CS(double masksize, double beamradpix, const char *outname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDout;

    copy_image_ID("pupaCS", "pupa0", 0);
    coronagraph_make_2Dprolate_DFT(masksize, beamradpix, 0.0, "outapo", size, "NULLim");
    IDout = image_ID("outapo");
    list_image_ID();

    // extend outapo into radial profile and 2D image
    // this is done by finding best fit ??
    // remap entrance pupil according to apodization profile


    return(0);
}



int coronagraph_update_2Dprolate(double masksizeld, double beamradpix, double centralObs, double zfactor)
{
    /* masksizeld is in lambda/d */
    FILE *fp;
    char fname[500];
    char fnameprev[500];
    char fnameinfo[500];
    char fnameinfoprev[500];
    char fnameprof[500];
    char fname1[500];
    char fname1prev[500];
    char fnamenew[500];

    char command[500];
    long size = CORONAGRAPHS_ARRAYSIZE;
    double FFcoeff = 0.0;
    long iter;
    long ID, ID0, ID1, ID2;
    long ii;
    long iter1;
    long NBiter1 = 1;
    double transm;
    double tmp1,err;
    int OK1;

    int CentralObstructionFlag = 0;
    long IDpupa0co;

    int result = 0;

    float masksizestep = 0.01;



    DFTZFACTOR = zfactor;
    if((ID=variable_ID("DFTZFACTOR"))!=-1)
        DFTZFACTOR = data.variable[ID].value.f;

    if(centralObs > 0.001)
    {
        CentralObstructionFlag = 1;
        printf("Central obstruction = %f\n", centralObs);
    }


    MASKSIZELD = masksizeld;
    OK1 = 0;
    /*  while(1)
      {
        if(OK1==1)
    {
      MASKSIZELD = 0.0;
      while(MASKSIZELD<0.05)
        MASKSIZELD = 0.8+0.1*((long) (23.0*ran1()));
        }*/

    for(MASKSIZELD=0.8; MASKSIZELD<5.01; MASKSIZELD+=masksizestep)
    {
        OK1=1;
        printf("--------- MASK SIZE = %f --------\n", MASKSIZELD);

        if(CentralObstructionFlag == 0)
        {
            sprintf(fname,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
            sprintf(fnameinfo,"%s/APLCapo/APLCapo_%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
        }
        else
        {
            sprintf(fname,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
            sprintf(fnameinfo,"%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
        }


        if(CentralObstructionFlag == 0)
        {
            sprintf(fnameprev,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, size);
            sprintf(fnameinfoprev,"%s/APLCapo/APLCapo_%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, size);
        }
        else
        {
            sprintf(fnameprev,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, centralObs, size);
            sprintf(fnameinfoprev,"%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, centralObs, size);
        }

        sprintf(fnamenew, "%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.new1", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);

        //     if((file_exists(fname)==0)||(file_exists(fnameinfo)==0)) // don't compute if file is already here
        if(file_exists(fnamenew)==0)
        {
            printf("Mask size = %f l/D\n", MASKSIZELD);
            sprintf(fnamenew, "%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.new1", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
            if(   ((file_exists(fnameprev)==1)&&(file_exists(fnameinfoprev)==1))  &&  ((file_exists(fname)==0)||(file_exists(fnameinfo)==0))  )
            {
                printf("COPY FROM PREVIOUS MASK SIZE\n");
                if(CentralObstructionFlag == 0)
                {
                    sprintf(fname,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
                    sprintf(fname1,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
                    sprintf(fnameinfo,"%s/APLCapo/APLCapo_%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
                    sprintf(fnameprev,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, size);
                    sprintf(fname1prev,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, size);
                    sprintf(fnameinfoprev,"%s/APLCapo/APLCapo_%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, size);
                }
                else
                {
                    sprintf(fname,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
                    sprintf(fname1,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
                    sprintf(fnameinfo,"%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
                    sprintf(fnameprev,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, centralObs, size);
                    sprintf(fname1prev,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, centralObs, size);
                    sprintf(fnameinfoprev,"%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD-masksizestep, centralObs, size);
                }


                sprintf(command, "cp %s %s", fnameprev, fname);
                printf("executing: %s\n", command);
                result = system(command);
                sprintf(command, "cp %s %s", fname1prev, fname1);
                printf("executing: %s\n", command);
                result = system(command);
                sprintf(command, "cp %s %s", fnameinfoprev, fnameinfo);
                printf("executing: %s\n", command);
                result = system(command);
            }

            ID0 = create_2Dimage_ID("apo", size, size);
            ID1 = create_2Dimage_ID("apo1", size, size);
            ID2 = create_2Dimage_ID("apo2", size, size);
            for(ii=0; ii<size*size; ii++)
                data.image[ID2].array.F[ii] = 0.0;

            fp = fopen("prol_conv.log","w");
            fclose(fp);
            iter1 = 0;
            while(iter1<NBiter1)
            {
                iter1++;
                iter = 0;
                while((iter!=2)&&(iter1<NBiter1+1))
                {
                    printf(" ---- %f lambda/d  [%ld / %ld] ----\n", MASKSIZELD, iter1, NBiter1);
                    if(CentralObstructionFlag == 0)
                    {
                        sprintf(fname,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
                        sprintf(fname1,"!%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
                        sprintf(fnameinfo,"%s/APLCapo/APLCapo_%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
                    }
                    else
                    {
                        sprintf(fname,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
                        sprintf(fname1,"!%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
                        sprintf(fnameinfo,"%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
                    }

                    load_fits(fname, "apostart", 1);
                    // save_fl_fits("apostart", "!apostart.fits");
                    //	  exit(0);
                    printf(".... STARTING OPTIMIZATION [%s] (%f %f %f)... \n", fname, MASKSIZELD, CORONAGRAPHS_PIXSCALE, centralObs);

                    transm = coronagraph_make_2Dprolate_DFT(MASKSIZELD/CORONAGRAPHS_PIXSCALE, beamradpix, centralObs, "out", size, "NULLim");

                    save_fl_fits("out",fname1);
                    ID = image_ID("out");
                    err = 0.0;
                    for(ii=0; ii<size*size; ii++)
                    {
                        data.image[ID1].array.F[ii] = data.image[ID2].array.F[ii];
                        data.image[ID2].array.F[ii] = data.image[ID].array.F[ii];
                        tmp1 = data.image[ID1].array.F[ii]-data.image[ID2].array.F[ii];
                        err += tmp1*tmp1;
                    }
                    if(iter1==0)
                        err = 0.0;
                    if(CentralObstructionFlag == 0)
                        sprintf(fnameprof,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.prof", CORONAGRAPHSDATALOCAL, MASKSIZELD, size);
                    else
                        sprintf(fnameprof,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.prof", CORONAGRAPHSDATALOCAL, MASKSIZELD, centralObs, size);
                    profile("out","out.prof", size/2, size/2, 1.0, size/2);
                    sprintf(command,"mv out.prof %s", fnameprof);
                    result = system(command);
                    sprintf(command,"cp APLCapo.%.3f.%.3f.info %s", MASKSIZELD, centralObs, fnameinfo);
                    printf("COMMAND = \"%s\"\n",command);
                    result = system(command);


                    delete_image_ID("apostart");

                    if(iter==1)
                    {
                        fp = fopen("prol_conv.log","a");
                        fprintf(fp,"%ld %.18f %.18g\n", iter1, transm, err);
                        fclose(fp);
                    }
                    iter ++;
                    iter1 ++;
                }

                for(ii=0; ii<size*size; ii++)
                    data.image[ID0].array.F[ii] = data.image[ID2].array.F[ii] + FFcoeff*(data.image[ID2].array.F[ii]-data.image[ID1].array.F[ii]);
                save_fl_fits("apo",fname1);
            }




            delete_image_ID("apo");
            delete_image_ID("apo1");
            delete_image_ID("apo2");

            sprintf(command, "touch %s", fnamenew);
            result = system(command);
        }
    }
    return(0);
}




double coronagraph_apofit_eval_r(double r)
{
  double value = 0.0;
  long n;
  double rp;
  double eps = 1.0e-8;
  double eps1 = 1.0e-4;

  double a = 8.0; // for extrapolation beyond r=1



  if(r>1.0+eps)
    {
      rp = r-1.0;
      rp = rp/(1.0 + a*rp*exp(-1.0/(a*rp)));
      rp += 1.0;
    }
  else
    rp = r;

  if(rp<eps1)
    rp = eps1;
  
  value = 0.0;

  for(n=0; n<fitapoN; n++)   
    value += fitapo_a[n] * exp(fitapo_b[n] * pow(rp, fitapo_c1[n]));
      
  if(isfinite(value)==0)
    {
      printf("ERROR : infinite value in coronagraph_apofit_eval_r\n");

      value = 0.0;
      for(n=0; n<fitapoN; n++)
	{
	  printf("        :::: %ld  %f %f %f :::::\n", n, fitapo_a[n], fitapo_b[n], fitapo_c[n]);	  
	}

    

      for(n=0; n<fitapoN; n++)   
	{
	  value += fitapo_a[n] * exp(fitapo_b[n] * pow(rp, fitapo_c[n])); 	  
	  printf("---- %ld %f   %f %f %f -> %f %f -> %f\n", n, rp, fitapo_a[n], fitapo_b[n], fitapo_c[n],  pow(rp, fitapo_c[n]),  exp(fitapo_b[n] * pow(rp, fitapo_c[n])), value);
	  fflush(stdout);
	}
    }
  

  return(value);
}



double  coronagraph_apofit_eval()
{
    double value;
    long i;
    double r, v;
    long n;
    double tmpv;


    value = 0.0;


    for(n=0; n<fitapoN; n++)
    {
        // penalize negative and small values of fitao_c
        value += 100.0*pow( 0.5*(  atan(-(fitapo_c[n]-0.5)*10.0/(M_PI/2))/(M_PI/2)+1.0  ), 8.0);
        fitapo_c1[n] = fitapo_c[n];
        if(fitapo_c1[n]<fitapo_minc)
            fitapo_c1[n] = fitapo_minc;

        // penalize large values of fitao_a and fitapo_b
        value += exp(-1.0/(fitapo_a[n]*fitapo_a[n]/200.0/200.0));
        value += exp(-1.0/(fitapo_b[n]*fitapo_b[n]/200.0/200.0));
    }




    // don't fit the first and last point as they may be corrupted
    for(i=1; i<aporawN-1; i++)
    {
        r = aporaw_r[i];
        v = coronagraph_apofit_eval_r(r);
        tmpv = aporaw_v[i] - v;
        value += tmpv*tmpv;
    }
    value = sqrt(value/aporawN);





    if(isfinite(value)==0)
    {
        printf("ERROR: INFINITE VALUE in coronagraph_apofit_eval ----\n");
        fflush(stdout);

        for(n=0; n<fitapoN; n++)
            printf("   :::: %ld  %f %f %f :::::\n", n, fitapo_a[n], fitapo_b[n], fitapo_c[n]);

        printf("aporawN = %ld\n", aporawN);
        value = 0.0;
        for(i=0; i<aporawN; i++)
        {
            r = aporaw_r[i];
            v = coronagraph_apofit_eval_r(r);
            tmpv = aporaw_v[i] - v;
            value += tmpv*tmpv;
            printf("%ld %f %f %f %f\n", i, r, v, tmpv, value);
            fflush(stdout);
            // if(isfinite(value)==0)
            //  exit(0);
        }
        value = sqrt(value/aporawN);
        exit(0);
    }

    return(value);
}





double coronagraph_apofit_f_evalmask (const gsl_vector *v, void *params)
{
  double *p = (double *)params;
  double value; 
  long zindex;
  long n;

  for(n=0; n<fitapoN; n++)
    {
      fitapo_a[n] = gsl_vector_get(v, 3*n);
      fitapo_b[n] = gsl_vector_get(v, 3*n+1);
      fitapo_c[n] = gsl_vector_get(v, 3*n+2);      
    }

  value = coronagraph_apofit_eval();
  
  if(LOOPCNT==0)
    {
      optval0 = value;
    }
  LOOPCNT++;
  
  //printf("[[%e]]", value);
  // fflush(stdout);

  return (value);

}



double coronagraph_apofit(const char *fnameout)
{
    FILE *fp;
    int ret;
    char fname[500];
    double value;
    long i;
    double r, v;
    long n;
    double tmpv;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    int status;
    gsl_multimin_function fpmeval_func;
    long NBoptVar;
    long iterMax;
    long ITERMAX = 200;
    long iter1;
    double eps1;
    double bestvalue;
    FILE *fpbest;
    double KEEPlimit, KEEPlimit1;
    int KEEP;
    long KEEPcnt;
    int OK = 0;
    long iter;
    double size;


    double aRANGE = 10.0;
    double bRANGE = 5.0;
    double cRANGE = 5.0;
    double cMin = 0.1;

    double *fitapo_a_best;
    double *fitapo_b_best;
    double *fitapo_c_best;

    double v1, v2, v3;

    printf("Fitting apodization profile with analytical function\n");
    fflush(stdout);

    value = coronagraph_apofit_eval();

    printf("Preparing optimization ... \n");
    fflush(stdout);
    NBoptVar = 3*fitapoN;
    x = gsl_vector_alloc (NBoptVar);
    ss = gsl_vector_alloc (NBoptVar);
    fpmeval_func.n = NBoptVar;  /* number of function components */
    fpmeval_func.f = &coronagraph_apofit_f_evalmask;
    fpmeval_func.params = (void *) NULL;
    s = gsl_multimin_fminimizer_alloc (T, NBoptVar);

    fitapo_a_best = (double*) malloc(sizeof(double)*fitapoN);
    fitapo_b_best = (double*) malloc(sizeof(double)*fitapoN);
    fitapo_c_best = (double*) malloc(sizeof(double)*fitapoN);



    // read previous best solution (if it exits)
    //  sprintf(fname, "apofitsol_params.txt");
    printf("Reading previous best solution [%s] ... ", fnameout);
    fflush(stdout);
    fpbest = fopen(fnameout, "r");
    OK = 0; // use previous solution ?
    if(fpbest!=NULL)
    {
        OK = 1;
        ret = fscanf(fpbest, "%lf", &bestvalue);
        printf("PREVIOUS BEST SOLUTION = %g\n", bestvalue);
        fflush(stdout);

        for(n=0; n<fitapoN; n++)
        {
            ret = fscanf(fpbest," %lf %lf %lf", &fitapo_a_best[n], &fitapo_b_best[n], &fitapo_c_best[n]);
            fitapo_a[n] = fitapo_a_best[n];
            fitapo_b[n] = fitapo_b_best[n];
            fitapo_c[n] = fitapo_c_best[n];
            if(fitapo_c[n]<0.0)
                OK = 0;
        }
        fclose(fpbest);
    }

    if(OK==0)
    {
        printf(" Using std starting point\n");
        fflush(stdout);
        bestvalue = 100.0;
        fitapo_a_best[0] = 1.0;
        fitapo_b_best[0] = -1.0;
        fitapo_c_best[0] = 1.0;
        for(n=1; n<fitapoN; n++)
        {
            fitapo_a_best[n] = 0.0;
            fitapo_b_best[n] = -1.0;
            fitapo_c_best[n] = 1.0;
        }
    }




    iterMax = 50000;
    printf("----------- STARTING DOWNHILL SIMPLEX ------------- [%ld]\n",ITERMAX);
    OK = 0; // did it improve ?

    KEEPlimit = bestvalue;
    KEEPlimit1 = 3.0*bestvalue;
    KEEP = 0;
    KEEPcnt = 0;
    eps1 = 0.99999999;
    for(iter1=0; iter1<ITERMAX; iter1++) // DOWNHILL SIMPLEX METHOD
    {
        printf("%05ld ", iter1);

        printf("[KL %e %e]  ", KEEPlimit, KEEPlimit1);
        fflush(stdout);

        if((iter1==0)||(OK==1)) // set starting point = best point if 1st iteration or loop is still making progress
        {
            KEEPcnt = 0;
            printf("BEST  ");
            fflush(stdout);
            for(n=0; n<fitapoN; n++)
            {
                gsl_vector_set (x, 3*n, fitapo_a_best[n]);
                gsl_vector_set (x, 3*n+1, fitapo_b_best[n]);
                gsl_vector_set (x, 3*n+2, fitapo_c_best[n]);
            }

            for(n=0; n<fitapoN; n++)
            {
                fitapo_a[n] = gsl_vector_get(x, 3*n);
                fitapo_b[n] = gsl_vector_get(x, 3*n+1);
                fitapo_c[n] = gsl_vector_get(x, 3*n+2);
            }
        }
        else
        {
            if(KEEP==0)  // random starting point
            {
                KEEPcnt = 0;
                printf("RAND  ");
                fflush(stdout);
                for(n=0; n<fitapoN; n++)
                {
                    v1 = (1.0-2.0*ran1())*aRANGE;
                    v2 = (1.0-2.0*ran1())*bRANGE;
                    v3 = 0.5+(2.0*ran1())*cRANGE;
                    gsl_vector_set (x, 3*n, v1);
                    gsl_vector_set (x, 3*n+1, v2);
                    gsl_vector_set (x, 3*n+2, v3);
                    printf("%lf %lf %lf\n", v1, v2, v3);
                }
                for(n=0; n<fitapoN; n++)
                {
                    fitapo_a[n] = gsl_vector_get(x, 3*n);
                    fitapo_b[n] = gsl_vector_get(x, 3*n+1);
                    fitapo_c[n] = gsl_vector_get(x, 3*n+2);
                }
            }
            else
            {
                printf("KEEP  ");
                fflush(stdout);
                KEEPcnt++;
                KEEPlimit = (0.94*KEEPlimit + 0.05*bestvalue)*0.7 + 0.3*s->fval;
                for(n=0; n<fitapoN; n++)
                {
                    gsl_vector_set (x, 3*n, fitapo_a[n]);
                    gsl_vector_set (x, 3*n+1, fitapo_b[n]);
                    gsl_vector_set (x, 3*n+2, fitapo_c[n]);
                }
                for(n=0; n<fitapoN; n++)
                {
                    fitapo_a[n] = gsl_vector_get(x, 3*n);
                    fitapo_b[n] = gsl_vector_get(x, 3*n+1);
                    fitapo_c[n] = gsl_vector_get(x, 3*n+2);
                }
            }
        }

        /* Set initial step sizes to 1e-8 */
        // printf("ss");
        //fflush(stdout);
        gsl_vector_set_all (ss, 1.0e-9);

        /* Initialize method and iterate */
        iter = 0;


        //      for(n=0; n<fitapoN; n++)
        //	printf(":::: %ld  %f %f %f :::::\n", n, fitapo_a[n], fitapo_b[n], fitapo_c[n]);


        printf("I");
        fflush(stdout);
        gsl_multimin_fminimizer_set (s, &fpmeval_func, x, ss);
        printf("I");
        fflush(stdout);

        OK = 0;
        LOOPCNT = 0;
        do
        {
            iter++;

            //  printf("S");
            //fflush(stdout);
            status = gsl_multimin_fminimizer_iterate(s);
            if (status)
                break;
            //printf("*");
            //fflush(stdout);


            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 1e-11);


            //	  printf ("............[%05ld] %e  %e  (%e) ", iter, s->fval, size, bestvalue);
            //fflush(stdout);


            if ((status == GSL_SUCCESS)||(s->fval < bestvalue * eps1))
            {
                if(status == GSL_SUCCESS)
                {
                    printf (" %e ->  %e [%5ld]  (%e)", optval0, s->fval, iter, bestvalue);
                    if(OK==1)
                        printf("  *");
                    printf("\n");
                }
                if(KEEP==0)
                {
                    if(s->fval < KEEPlimit1)
                        KEEPlimit1 = 0.8*KEEPlimit1 + 0.2*s->fval;
                    else
                        KEEPlimit1 = 0.98*KEEPlimit1 + 0.02*s->fval;
                }

                if(s->fval < bestvalue * eps1) //if(s->f<bestvalue)
                {
                    OK = 1;

                    for(n=0; n<fitapoN; n++)
                    {
                        fitapo_a_best[n] = fitapo_a[n];
                        fitapo_b_best[n] = fitapo_b[n];
                        fitapo_c_best[n] = fitapo_c[n];
                    }

                    //		  sprintf(fname,"apofitsol_params.txt");
                    fpbest = fopen(fnameout, "w");
                    bestvalue = s->fval;
                    fprintf(fpbest,"%.20g", bestvalue);

                    for(n=0; n<fitapoN; n++)
                        fprintf(fpbest," %15.10f %15.10f %15.10f", fitapo_a_best[n], fitapo_b_best[n], fitapo_c_best[n]);
                    fprintf(fpbest, "\n");
                    n = 0;
                    fprintf(fpbest, "f(x) = %15.10f * exp(%15.10f * x**%15.10f)", fitapo_a_best[n], fitapo_b_best[n], fitapo_c_best[n]);
                    for(n=1; n<fitapoN; n++)
                        fprintf(fpbest,"+ %15.10f * exp(%15.10f * x**%15.10f)", fitapo_a_best[n], fitapo_b_best[n], fitapo_c_best[n]);
                    fprintf(fpbest, "\n");
                    fclose(fpbest);
                }
                if(KEEP==0)
                {
                    if(s->fval<KEEPlimit1)
                        KEEP = 1;
                    else
                        KEEP = 0;
                    KEEPlimit = KEEPlimit1;
                }
                else
                {
                    if(s->fval<KEEPlimit)
                        KEEP = 1;
                    else
                        KEEP = 0;
                    if(KEEPcnt > 50)
                        KEEP = 0;
                }
            }
            //	  printf("\n");
            //fflush(stdout);
        }
        while (status == GSL_CONTINUE && iter < iterMax );
        if(iter>iterMax-1) {
            printf("Max number of iteration reached... starting over\n");
            KEEP = 0;
        }
    }
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);


    for(n=0; n<fitapoN; n++)
    {
        fitapo_a[n] = fitapo_a_best[n];
        fitapo_b[n] = fitapo_b_best[n];
        fitapo_c[n] = fitapo_c_best[n];
    }

    free(fitapo_a_best);
    free(fitapo_b_best);
    free(fitapo_c_best);



    fp = fopen("apofitsol.txt", "w");
    for(i=0; i<aporawN; i++)
    {
        r = aporaw_r[i];
        v = 0.0;
        for(n=0; n<fitapoN; n++)
        {
            v += fitapo_a[n] * exp(fitapo_b[n] * pow(r,fitapo_c[n]));
        }
        fprintf(fp, "%11.8f %11.8f %11.8f\n", r, aporaw_v[i], v);
    }
    fclose(fp);



    return(bestvalue);
}




//
// compile prolate apodization parameters in a single set of FITS files
// each FITS file has 2 dimensions: x = central obstruction, y = focal plane mask radius
//
//
int_fast8_t coronagraph_APLCapo_compile()
{
    FILE *fp;
    int ret;

    struct stat buffer;


    double FITVLIMIT = 2.0e-5;

    double co = 0.0; // central obstruction
    double fpmrad = 0.8; // focal plane mask radius
    long ii, jj;
    long iisize, jjsize;

    char fnameprolinfo[500];
    char fnameout[500];
    double tmpf0, tmpf1, tmpf2, tmpf3;

    int fitOK;
    char fnameapofit[500];
    char fnameapoprof[500];
    long n;
    double fitv;
    FILE *fpfit;
    char command[500];
    int result = 0;
    long cnt, i;
    double tmpf;
    double eps0 = 1.0e-8;



    // CREATE IMAGES
    iisize = (long) ((APLCapo_CO_END-APLCapo_CO_START)/APLCapo_CO_STEP);
    jjsize = (long) ((APLCapo_FPMRAD_END-APLCapo_FPMRAD_START)/APLCapo_FPMRAD_STEP);

    IDprol_init = create_2Dimage_ID("pinit", iisize, jjsize);
    IDprol_ffrac = create_2Dimage_ID("pffrac", iisize, jjsize);
    IDprol_transm = create_2Dimage_ID("ptransm", iisize, jjsize);
    IDprol_peak = create_2Dimage_ID("ppeak", iisize, jjsize);

    IDprol_fitapo_a = create_3Dimage_ID("fitapoa", iisize, jjsize, fitapoN);
    IDprol_fitapo_b = create_3Dimage_ID("fitapob", iisize, jjsize, fitapoN);
    IDprol_fitapo_c = create_3Dimage_ID("fitapoc", iisize, jjsize, fitapoN);
    IDprol_fitfit = create_2Dimage_ID("fitfit", iisize, jjsize);

    //  for(n=0; n<fitapoN; n++)
    // value += fitapo_a[n] * exp(fitapo_b[n] * pow(rp,fitapo_c[n]));


    printf("\n\n");

    for(ii=0; ii<iisize; ii++)
    {
        for(jj=0; jj<jjsize; jj++)
            if(ran1()<1.01)
            {

                co = APLCapo_CO_START + APLCapo_CO_STEP*ii;
                fpmrad = APLCapo_FPMRAD_START + APLCapo_FPMRAD_STEP*jj;

                // .info file
                //
                // col 1 is mask size in SYSTEM l/D
                // col 2 is fraction of light going through FPmask
                // col 3 is prolate throughput if not done by PIAA
                // col 4 is peak
                // ideal complex transmission for FPM is  MASKCAMULT = -(1.0-peak)/peak;


                // load info file for apodization
                if(co < 0.0001)
                {
                    sprintf(fnameprolinfo, "%s/APLCapo/APLCapo_%5.3f.4096.ref.info", CORONAGRAPHSDATALOCAL, fpmrad);
                    sprintf(fnameapoprof, "%s/APLCapo/raw/APLCapo_%.3f.4096.ref.prof", CORONAGRAPHSDATALOCAL, fpmrad);
                }
                else
                {
                    sprintf(fnameprolinfo, "%s/APLCapo/APLCapo_%5.3f.%5.3f.4096.ref.info", CORONAGRAPHSDATALOCAL, fpmrad, co);
                    sprintf(fnameapoprof, "%s/APLCapo/raw/APLCapo_%.3f.%.3f.4096.ref.prof", CORONAGRAPHSDATALOCAL, fpmrad, co);
                }


                if((stat(fnameprolinfo, &buffer)==0)&&(stat(fnameapoprof, &buffer)==0)) // if file exists
                {
                    printf("\r");
                    printf("%5ld/%5ld   %5ld/%5ld  %f %f  ", ii, iisize, jj, jjsize, co, fpmrad);


                    fp = fopen(fnameprolinfo, "r");
                    ret = fscanf(fp, "%lf %lf %lf %lf\n", &tmpf0, &tmpf1, &tmpf2, &tmpf3);
                    fclose(fp);
                    //	    fpmasktransm =  -(1.0-tmpf3)/tmpf3;
                    //data.image[IDfpmasktransm].array.F[fpmaskrad_index*obs1arraysize + obs1_index] = fpmasktransm;

                    data.image[IDprol_init].array.F[jj*iisize+ii] = 1.0;
                    data.image[IDprol_ffrac].array.F[jj*iisize+ii] = tmpf1;
                    data.image[IDprol_transm].array.F[jj*iisize+ii] = tmpf2;
                    data.image[IDprol_peak].array.F[jj*iisize+ii] = tmpf3;



                    sprintf(fnameapofit, "%s/APLCapo/APLCapo_%.3f.%.3f.4096.ref.prof.fitparams", CORONAGRAPHSDATALOCAL, fpmrad, co);


                    fitapo_a = (double*) malloc(sizeof(double)*fitapoN);
                    fitapo_b = (double*) malloc(sizeof(double)*fitapoN);
                    fitapo_c = (double*) malloc(sizeof(double)*fitapoN);
                    fitapo_c1 = (double*) malloc(sizeof(double)*fitapoN);


                    fitapo_a[0] = 1.0;
                    fitapo_b[0] = -1.0;
                    fitapo_c[0] = 2.0;
                    for(n=0; n<fitapoN; n++)
                    {
                        fitapo_a[n] = 0.0;
                        fitapo_b[n] = -1.0;
                        fitapo_c[n] = 2.0+n;
                    }

                    fitOK = 0;
                    if(stat(fnameapofit, &buffer)==0) // if file exists
                    {
                        fp = fopen(fnameapofit, "r");
                        if(fp==NULL)
                        {
                            printf("ERROR: file \"%s\" does not exist\n", fnameapofit);
                            exit(0);
                        }

                        ret = fscanf(fp, "%lf", &fitv);
                        printf("FIT SOLUTION RESIDUAL = %g ", fitv);
                        fflush(stdout);
                        for(n=0; n<fitapoN; n++)
                            ret = fscanf(fp," %lf %lf %lf", &fitapo_a[n], &fitapo_b[n], &fitapo_c[n]);
                        fclose(fp);
                        fitOK = 1;
                        if(fitv>FITVLIMIT)
                            fitOK = 0;
                    }
                    else
                    {
                        printf("                          ");
                        fflush(stdout);
                    }

                    if(fitOK==0) // if file does not exist, or fit is not sufficiently good
                    {
                        if(co>0.0001)
                            sprintf(fnameapoprof, "%s/APLCapo/raw/APLCapo_%.3f.%.3f.4096.ref.prof", CORONAGRAPHSDATALOCAL, fpmrad, co);
                        else
                            sprintf(fnameapoprof, "%s/APLCapo/raw/APLCapo_%.3f.4096.ref.prof", CORONAGRAPHSDATALOCAL, fpmrad);


                        sprintf(command,"awk '{print $1,$2}' %s > apotmp.prof", fnameapoprof);
                        result = system(command);

                        if((fp=fopen("apotmp.prof","r"))==NULL)
                        {
                            printf("ERROR : cannot open file \"apotmp.prof\"\n");
                            exit(0);
                        }
                        else
                        {
                            printf("FILE \"apotmp.prof\" opened \n");
                            fflush(stdout);
                        }

                        // FIT APODIZATION WITH ANALYTICAL FUNCTION
                        aporaw_r = (double*) malloc(sizeof(double)*PIAAAPO_NBPOINTS);
                        aporaw_v = (double*) malloc(sizeof(double)*PIAAAPO_NBPOINTS);

                        cnt = 0;
                        for(i=0; i<PIAAAPO_NBPOINTS; i++)
                        {
                            result = fscanf(fp,"%lf %lf\n", &tmpf0, &tmpf);
                            if(tmpf>eps0)
                            {
                                aporaw_r[cnt] = (double) (1.0*tmpf0/200.0);
                                aporaw_v[cnt] = (double) tmpf;
                                cnt++;
                            }
                        }
                        fclose(fp);
                        aporawN = cnt;


                        fitv = coronagraph_apofit(fnameapofit);
                        printf("fit value = %g\n", fitv);

                        fpfit = fopen(fnameapofit, "r");
                        ret = fscanf(fpfit, "%lf", &fitv);
                        printf("FIT SOLUTION RESIDUAL = %g\n", fitv);
                        fflush(stdout);
                        for(n=0; n<fitapoN; n++)
                            ret = fscanf(fpfit," %lf %lf %lf", &fitapo_a[n], &fitapo_b[n], &fitapo_c[n]);
                        fclose(fpfit);


                        free(aporaw_r);
                        free(aporaw_v);
                    }

                    quick_sort3_double(fitapo_a, fitapo_b, fitapo_c, fitapoN);

                    for(n=0; n<fitapoN; n++)
                    {
                        data.image[IDprol_fitapo_a].array.F[n*iisize*jjsize+jj*iisize+ii] = fitapo_a[n];
                        data.image[IDprol_fitapo_b].array.F[n*iisize*jjsize+jj*iisize+ii] = fitapo_b[n];
                        data.image[IDprol_fitapo_c].array.F[n*iisize*jjsize+jj*iisize+ii] = fitapo_c[n];
                    }
                    data.image[IDprol_fitfit].array.F[jj*iisize+ii] = fitv;

                    free(fitapo_a);
                    free(fitapo_b);
                    free(fitapo_c);
                    free(fitapo_c1);
                }
                else
                    data.image[IDprol_init].array.F[jj*iisize+ii] = 0.0;


            }
    }

    sprintf(fnameout, "!%s/APLCapo_init.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("pinit", fnameout);

    sprintf(fnameout, "!%s/APLCapo_ffrac.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("pffrac", fnameout);

    sprintf(fnameout, "!%s/APLCapo_transm.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("ptransm", fnameout);

    sprintf(fnameout, "!%s/APLCapo_peak.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("ppeak", fnameout);


    sprintf(fnameout, "!%s/APLCapo_fitapo_a.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("fitapoa", fnameout);
    sprintf(fnameout, "!%s/APLCapo_fitapo_b.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("fitapob", fnameout);
    sprintf(fnameout, "!%s/APLCapo_fitapo_c.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("fitapoc", fnameout);
    sprintf(fnameout, "!%s/APLCapo_fitfit.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("fitfit", fnameout);




    delete_image_ID("pffrac");
    delete_image_ID("ptransm");
    delete_image_ID("ppeak");
    delete_image_ID("fitfit");

    delete_image_ID("fitapoa");
    delete_image_ID("fitapob");
    delete_image_ID("fitapoc");
    printf("\n");

    return(0);
}









//
// this function only works for circular symmetry PIAA
// it uses an analytical fit to apodized pupil function
//
//
int coronagraph_init_PIAA()
{
    FILE *fp;
    FILE *fpfit;
    int ret;
    int result = 0;
    long size = CORONAGRAPHS_ARRAYSIZE;
    long i,ii,jj;
    double *r0, *r1;
    long NBpoints0 = piaaconfNPUPFILESIZE;
    double *pup;
    double *pupsum,*pup0sum;
    // long tmpl;
    long NBpoints;
    double tmpf,tmpf0;
    double F0,F1,F2,fluxdens,x,total, total1, total2, tmp;
    long NB_profile_pts;
    double dr0,dr1,ndr0;
    double psr0,psr1;
    long *cntarray;
    double rad,v1,v2;
    long factor = 3;
    long ii0, ii1,jj1,ri1;
    double y,r,ri,alpha,value;
    char fname[200];
    long icentobs;
    double epsilon = 0.000000000001;
    char command[400];

    int file_write = 1; // 1 if writing output files
    FILE *fpout;
    FILE *fp1;

    double eps0 = 1.0e-8;
    double lim0, lim1; // inner and outer edges in output profile, estimated from profile itself
    double FLUXIN = 0.0001; // fraction of total flux to insert in central obstruction

    long i0, i1;
    double d2;

    long n;
    long cnt;
    double fitv;

    double a, b, eps1;
    long IDt, IDtt;

    double total_ir = 0.0;
    double total_or = 0.0;
    double normfactor = 1.0;

    long NBistep;
    double t0, t0cnt;
    double dir, odir;
    double verr;
    double bstep;
    double *innerprof_cumul;

    double dx, dz, dist, slope, y3, r0c, r1c, r0n, r1n;
    double piaaconfdist;
    double piaaconfbeamradius;


    long tmpstep = 1;



    printf("PIAA init    %ld points\n", PIAAAPO_NBPOINTS);
    printf("piaaconfNPUPFILESIZE = %ld\n", (long) piaaconfNPUPFILESIZE);
    printf("NBpoints0 = %ld\n", NBpoints0);
    printf("APLC_CentOBS0 = %f\n", APLC_CentOBS0);
    printf("APLC_CentOBS1 = %f\n", APLC_CentOBS1);
    fflush(stdout);

    PIAAextfactor1 = 2.0; // shape extrapolation oversize factor at output





    //
    // STEP 1 : GET APODIZATION RADIAL PROFILE
    //
    //
    //
    //
    //
    //
    //

    // allocate memory for profile
    PIAAAPO = (double*) malloc(sizeof(double)*PIAAAPO_NBPOINTS);
    PIAA_HYBRID_CPAAPO = (double*) malloc(sizeof(double)*PIAAAPO_NBPOINTS);






    //
    // preprocess apo profile file
    //

    /*
    sprintf(command,"awk '{print $1,$2}' %s/APLCapo/raw/%s > apotmp.prof", CORONAGRAPHSDATALOCAL, PIAAAPO_FNAME);
    result = system(command);

    if((fp=fopen("apotmp.prof","r"))==NULL)
      {
        printf("ERROR : cannot open file apotmp.prof\n");
        exit(0);
      }
    else
      {
        printf("FILE apotmp.prof opened \n");
        fflush(stdout);
      }
    */



    // FIT APODIZATION WITH ANALYTICAL FUNCTION
    /*
    aporaw_r = (double*) malloc(sizeof(double)*PIAAAPO_NBPOINTS);
    aporaw_v = (double*) malloc(sizeof(double)*PIAAAPO_NBPOINTS);

    cnt = 0;
    for(i=0;i<PIAAAPO_NBPOINTS;i++)
      {
        result = fscanf(fp,"%lf %lf\n", &tmpf0, &tmpf);
        if(tmpf>eps0)
    {
      aporaw_r[cnt] = (double) (1.0*tmpf0/200.0);
      aporaw_v[cnt] = (double) tmpf;
      cnt++;
    }
      }
    fclose(fp);
    aporawN = cnt;

    fitapo_a = (double*) malloc(sizeof(double)*fitapoN);
    fitapo_b = (double*) malloc(sizeof(double)*fitapoN);
    fitapo_c = (double*) malloc(sizeof(double)*fitapoN);
    fitapo_c1 = (double*) malloc(sizeof(double)*fitapoN);


    fitapo_a[0] = 1.0;
    fitapo_b[0] = -1.0;
    fitapo_c[0] = 2.0;
    fitapo_c1[0] = 2.0;
    for(n=0; n<fitapoN; n++)
      {
        fitapo_a[n] = 0.0;
        fitapo_b[n] = -1.0;
        fitapo_c[n] = 2.0+n;
        fitapo_c1[n] = 2.0+n;
      }




    sprintf(fname, "%s/APLCapo/%s.fitparams", CORONAGRAPHSDATALOCAL, PIAAAPO_FNAME);
    printf("Reading previous best fit solution [%s] ... ", fname);
    fflush(stdout);
    fpfit = fopen(fname, "r");
    if(fpfit==NULL)
      {
        fitv = coronagraph_apofit(fname);
        printf("fit value = %g\n", fitv);
        fpfit = fopen(fname, "r");
      }

    ret = fscanf(fpfit, "%lf", &fitv);
    printf("FIT SOLUTION RESIDUAL = %g\n", fitv);
    fflush(stdout);
    for(n=0; n<fitapoN; n++)
      {
        ret = fscanf(fpfit," %lf %lf %lf", &fitapo_a[n], &fitapo_b[n], &fitapo_c[n]);
        fitapo_c1[n] = fitapo_c[n];
        if(fitapo_c1[n]<fitapo_minc)
    fitapo_c1[n] = fitapo_minc;

      }
    fclose(fpfit);

    n = 0;
    printf("f(x) = %15.10f * exp(%15.10f * x**%15.10f)", fitapo_a[n], fitapo_b[n], fitapo_c1[n]);
    for(n=1; n<fitapoN; n++)
      printf("+ %15.10f * exp(%15.10f * x**%15.10f)", fitapo_a[n], fitapo_b[n], fitapo_c1[n]);
    printf("\n");



    free(aporaw_r);
    free(aporaw_v);

    */


    /*
    lim0 = 0.0;
    lim1 = -1.0;
    if((fp=fopen("apotmp.prof","r"))==NULL)
      {
        printf("ERROR : cannot open file apotmp.prof\n");
        exit(0);
      }
    else
      {
        printf("FILE apotmp.prof opened \n");
        fflush(stdout);
      }
    for(i=0;i<PIAAAPO_NBPOINTS;i++)
      {
        result = fscanf(fp,"%lf %lf\n", &tmpf0, &tmpf);
        PIAAAPO[i] = tmpf;
        if(tmpf>eps0)
    {
      lim1 = tmpf0;
      i1 = i;
    }
        if((tmpf<eps0)&&(lim1<0.0))
    {
      lim0 = tmpf0;
      i0 = i;
    }
          }
    fclose(fp);
    printf("PROFILE INNER EDGE = %f\n", lim0);
    printf("PROFILE OUTER EDGE = %f\n", lim1);
    printf("CENTRAL OBS = %f\n", lim0/lim1);

    */






    // HYBRID: SEPARATE APODIZATION INTO HYBRID + PIAA
    /* tmp = PIAA_HYBRID_CST*PIAAAPO[0];
    for(i=0;i<PIAAAPO_NBPOINTS;i++)
      {
        PIAA_HYBRID_CPAAPO[i] = PIAAAPO[i]/(tmp+PIAAAPO[i]);
        PIAAAPO[i] += tmp;
      }
    */

    // CENTRAL OBSTRUCTION
    //
    // input central obstruction = APLC_CentOBS0
    // output central obs = lim0/lim1
    //
    // filling in central obstruction:
    //
    // sb0 = FLUXIN/APLC_CentOBS0^2
    // sb1 = FLUXIN/(lim0/lim1)^2
    //

    // Flux inside central obstruction for output apodization function
    //






    /*  PIAAFLUXFACTOR = total1/total;*/
    printf("factor = %ld\n", factor);

    total = 0.0;
    total1 = 0.0;
    total_or = 0.0; // outer region
    total_ir = 0.0; // inner region
    IDt = create_2Dimage_ID("apo2Dim", size, size);
    IDtt = create_2Dimage_ID("apo2Dimt", size, size); // truncated


    for(ii=0; ii<factor*size; ii+=tmpstep)
    {
        ii1 = (long) ((1.0*ii-0.5)/factor+0.5);

        //      printf("%ld/%ld %ld %ld\n", ii, size, factor, ii1);
        //fflush(stdout);

        for(jj=0; jj<factor*size; jj+=tmpstep)
        {
            jj1 = (long) ((1.0*jj-0.5)/factor+0.5);
            x = 1.0/factor*(ii-factor*size/2);
            y = 1.0/factor*(jj-factor*size/2);
            r = sqrt(x*x+y*y); /* in pixel */
            r /= CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
            //	 printf("%ld %ld   %ld %ld    %f\n", ii, jj, ii1, jj1, r);
            // fflush(stdout);
            value = 0.0;
            if(r<PIAAextfactor1)
            {
                value = coronagraph_apofit_eval_r(r);
                if(r>1.0) // OUTER RING
                    total_or += value*value;
                else
                {
                    if(r<APLC_CentOBS1) // INSIDE CENTRAL OBSTRUCTION
                        total_ir += value*value;
                    else // WITHIN PUPIL
                    {
                        total += value*value;
                        total1 += 1.0;
                    }
                }
            }
            if((ii1>0)&&(jj1>0)&&(ii1<size)&&(jj1<size))
            {
                PIAAAPO2D[jj1*size+ii1] += value/factor/factor;
                data.image[IDt].array.F[jj1*size+ii1] =  PIAAAPO2D[jj1*size+ii1];

                if((r>APLC_CentOBS1)&&(r<1.0))
                    data.image[IDtt].array.F[jj1*size+ii1] =  PIAAAPO2D[jj1*size+ii1];
                else
                    data.image[IDtt].array.F[jj1*size+ii1] = 0.0;
            }
        }
    }
    normfactor = sqrt(total1/total);
    printf("NORMALIZATION FACTOR = %f\n", normfactor);
    fflush(stdout);

    for(ii=0; ii<size*size; ii++)
        data.image[IDt].array.F[ii] *= normfactor;


    printf("%g %g %g\n", total_ir/total, total/total, total_or/total);

    total_ir /= total;
    printf("Light within central obstruction (0->%f) = %g\n", APLC_CentOBS1, total_ir);
    total_or /= total;
    printf("Light in outer ring (1->%f) = %g\n", PIAAextfactor1, total_or);

    save_fl_fits("apo2Dim", "!apo2Dim.fits");
    save_fl_fits("apo2Dimt", "!apo2Dimt.fits");
    FLUXIN = total_ir;

    //  system("rm apo2Dim.prof");
    //profile("apo2Dim", "apo2Dim.prof", 2048, 2048, 1.0, 500);


    printf("GOAL = %g\n", total_ir/(APLC_CentOBS0*APLC_CentOBS0/(1.0-APLC_CentOBS0*APLC_CentOBS0)));



    //  exit(0);


    //
    // STEP 2: EXTEND PROFILE INWARD AND OUTWARD TO ENSURE OPTICS SHAPE SMOOTHNESS
    //


    // Compute inner pseudo profile
    b = 0.5;
    bstep = 0.1;
    verr = 1.0;
    NBistep = 1000000;
    innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);
	dir = 1.0; // initial direction
    while(fabs(verr)>1.0e-9)
    {
        t0 = 0.0;
        t0cnt = 0.0;
        for(i=0; i<NBistep; i++)
        {
            x = 1.0*i/NBistep;
            a = 0.5;
            eps1 = 1e-8;
            if(x<eps1)
                x = eps1;
            if(x>1.0-eps1)
                x = 1.0-eps1;
            value =  b + (1.0-b)*(0.5+atan(-a/b/(x*x) + a/pow(x-1.0,2))/M_PI);
            t0 += x*value*value;
            innerprof_cumul[i] = t0;
            t0cnt += x;
        }

        verr = t0/t0cnt - total_ir/(APLC_CentOBS0*APLC_CentOBS0/(1.0-APLC_CentOBS0*APLC_CentOBS0));

        odir = dir;
        if(verr>0.0) // too much light
        {
            b /= (1.0+bstep);
            dir = -1.0;
        }
        else
        {
            b *= (1.0+bstep);
            dir = 1.0;
        }
        if(odir*dir<0.0)
            bstep *= 0.1;
        printf(".");
        fflush(stdout);
    }
    printf("\n");
    printf("TOTAL = %f -> %g (%g %g)\n", b, t0/t0cnt, bstep, verr);

    for(i=0; i<NBistep; i++)
        innerprof_cumul[i] *= FLUXIN/innerprof_cumul[NBistep-1];
    // printf("LAST PT : %g\n", innerprof_cumul[NBistep-1]);



    // compute oversizing of INPUT beam
    PIAAextfactor0 = sqrt((total_ir + 1.0 + total_or) / (1.0 + total_ir));
    printf("PIAAextfactor0 = %f  \n", PIAAextfactor0 );
    piaaconfpup_amp_profile = (double*) malloc(sizeof(double)*piaaconfNPUPFILESIZE);

    jj = 0;
    for(ii=0; ii<piaaconfNPUPFILESIZE; ii++)
    {
        r = PIAAextfactor1*ii/piaaconfNPUPFILESIZE;
        piaaconfpup_amp_profile[ii] = normfactor*coronagraph_apofit_eval_r(r);
    }



    if (piaaconfr0 != NULL)
        free(piaaconfr0);
    if (piaaconfr1 != NULL)
        free(piaaconfr1);

    r0 = (double*) malloc(sizeof(double)*NBpoints0);
    r1 = (double*) malloc(sizeof(double)*NBpoints0);
    pup = (double*) malloc(sizeof(double)*piaaconfNPUPFILESIZE);
    pupsum = (double*) malloc(sizeof(double)*piaaconfNPUPFILESIZE);
    pup0sum = (double*) malloc(sizeof(double)*NBpoints0);

    /* computing r0 and r1 */
    /* r0 and r1 are dimensionless */

    /* first, r0 is evenly distributed on the first optic */
    for(i=0; i<NBpoints0; i++)
        r0[i] = PIAAextfactor0*i/(NBpoints0-1);

    if(APLC_CentOBS0>0.01) // Input central obstruction
        PIAACENTOBS = APLC_CentOBS0;

    //  icentobs = 1.0*PIAACENTOBS*NBpoints0;



    for(i=0; i<NBpoints0; i++)
    {
        if(r0[i]<APLC_CentOBS0)
        {
            x = r0[i]/APLC_CentOBS0;
            i1 = (long) (x*NBistep);
            if(i1>NBistep-1)
                i1 = NBistep-1;
            pup0sum[i] = innerprof_cumul[i1];
        }
        else
            pup0sum[i] = FLUXIN + (r0[i]*r0[i]-PIAACENTOBS*PIAACENTOBS)/(1.0-PIAACENTOBS*PIAACENTOBS);
        //      else
        //pup0sum[i] = FLUXIN*(r0[i]*r0[i]/(PIAACENTOBS*PIAACENTOBS));
    }
    for(i=0; i<NBpoints0; i++)
        pup0sum[i] /= 1.0+total_ir+total_or;


    /* computing r1 from r0 */
    NB_profile_pts = piaaconfNPUPFILESIZE;
    total=0.0;
    for(ii=0; ii<NB_profile_pts; ii++)
    {
        total += piaaconfpup_amp_profile[ii]*piaaconfpup_amp_profile[ii]*ii;
        pupsum[ii] = total;
    }


    i=0;
    ii=0;
    cnt = 0;
    r0[0] = 0.0;
    r1[0] = 0.0;
    for(i=1; i<NBpoints0; i++)
    {
        F0 = pup0sum[i];
        while(((pupsum[ii]/pupsum[(long) (NB_profile_pts)-1])<pup0sum[i])&&(ii<piaaconfNPUPFILESIZE))
            ii++;
        F1 = pupsum[ii-1]/pupsum[(long) (NB_profile_pts)-1];
        F2 = pupsum[ii]/pupsum[(long) (NB_profile_pts)-1];

        /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
        if(fabs(F2-F1)>0.0000000001)
            fluxdens = (F2-F1)/(2.0*ii-1.0);
        else
            fluxdens = 0.0000000001;
        x = sqrt((F0-F1)/fluxdens+(1.0*ii*ii-2.0*ii+1.0)) + 1.0 - 1.0*ii;

        r1[i] = PIAAextfactor1*(1.0*ii-1.0+x)/(NB_profile_pts-1);
        //      printf("%f %ld/%ld %ld %.18f %.18f %.18f %.18f %.18f %.18f\n", r0[i] , i, icentobs, ii, F0, F1, F2, x, fluxdens, r1[i]);
    }


    free(pup0sum);

    printf("Computing number of points [%ld]...", NBpoints0);
    fflush(stdout);
    i = 0;
    ii = 0;
    tmp = 0.0;
    //  fp1 = fopen("test_prof.txt", "w");
    while(i<NBpoints0-1)
    {
        dr0 = r0[i+1]-r0[i];
        dr1 = r1[i+1]-r1[i];
        if(dr1<epsilon)
        {
            dr1 = epsilon;
            /*	  printf("*");*/
        }

        ndr0 = dr0*dr0/dr1;
        if(dr0<ndr0)
            ndr0 = dr0;

        tmp += ndr0;
        ii++;
        if(tmp>r0[i])
            i++;
        // 1: i = index i, linear with r0
        // 2: ii =
        // 3: dr0 = derivative of r0 against i (should be constant)
        // 4: dr1 = derivative of r1 against i
        // 5: ndr0
        // 6: r0[i] = radius in input beam
        // 7: r1[i] = radius in output beam
        // 8: r0[i+1]
        // 9: r1[i+1]
        //  fprintf(fp1, "%ld %ld %g %g %g %g %g %g %g\n", i, ii, dr0, dr1, ndr0, r0[i], r1[i], r0[i+1], r1[i+1]);
    }
    //  fclose(fp);



    NBpoints = ii;
    printf("\n%ld Points\n", NBpoints);
    fflush(stdout);



    piaaconfr0 = (double*) malloc(sizeof(double)*NBpoints);
    if(piaaconfr0==NULL)
    {
        printf("COULD NOT ALLOCATE MEMORY FOR piaaconfr0\n");
        exit(0);
    }

    piaaconfr1 = (double*) malloc(sizeof(double)*NBpoints);
    if(piaaconfr0==NULL)
    {
        printf("COULD NOT ALLOCATE MEMORY FOR piaaconfr1\n");
        exit(0);
    }

    printf("Filling up nr0 ...");
    fflush(stdout);
    /* fill up nr0 */
    i = 0;
    ii = 0;
    piaaconfr0[ii] = 0.0;


    cnt = 0;
    while((i<NBpoints0-1)&&(ii<NBpoints-1))
    {
        dr0 = r0[i+1] - r0[i];
        dr1 = r1[i+1] - r1[i];

        ndr0 = dr0*dr0/dr1;
        /* ndr0 is big at the center of the beam, small at the edge */
        /* at the center of the beam, ndr0>dr0, so we want to use dr0 sampling */
        /* at the edge, ndr0<dr0, we want to use ndr0 */
        /* when we do that, we are always properly sampled both in r0 space and r1 space */
        if(dr0<ndr0) /* equivalent to (dr0>dr1) test */
            ndr0 = dr0;
        if(ndr0<0.0)
            ndr0 = 0.0;

        piaaconfr0[ii+1] = piaaconfr0[ii] + ndr0;
        ii++;
        if(piaaconfr0[ii]>r0[i])
            i++;
    }
    free(r0);
    free(r1);




    printf("----- %ld    %ld %ld ----\n", ii, NBpoints, NBpoints0);

    pup0sum = (double*) malloc(sizeof(double)*NBpoints);

    for(i=0; i<NBpoints; i++)
    {
        if(piaaconfr0[i]<APLC_CentOBS0)
        {
            x = piaaconfr0[i]/APLC_CentOBS0;
            i1 = (long) (x*NBistep);
            if(i1>NBistep-1)
                i1 = NBistep-1;
            pup0sum[i] = innerprof_cumul[i1];
        }
        else
            pup0sum[i] = FLUXIN + (piaaconfr0[i]*piaaconfr0[i]-PIAACENTOBS*PIAACENTOBS)/(1.0-PIAACENTOBS*PIAACENTOBS);
        //      else
        //pup0sum[i] = FLUXIN*(piaaconfr0[i]*piaaconfr0[i]/(PIAACENTOBS*PIAACENTOBS));
    }
    for(i=0; i<NBpoints; i++)
        pup0sum[i] /= 1.0+total_ir+total_or;

    free(innerprof_cumul);


    /* compute nr1 */
    total=0.0;
    for(ii=0; ii<NB_profile_pts; ii++)
    {
        total += 2.0*PI*piaaconfpup_amp_profile[ii]*piaaconfpup_amp_profile[ii]*ii; /* new */
        pupsum[ii]=total;
    }
    i=0;
    ii=0;

    piaaconfr1[0] = 0.0;
    for(i=1; i<NBpoints; i++)
    {
        F0 = pup0sum[i]; // fraction of total light inside input radius in input pupil
        while(((pupsum[ii]/pupsum[NB_profile_pts-1])<pup0sum[i])&&(ii<NB_profile_pts))
            ii++;
        if(ii!=0)
            F1 = pupsum[ii-1]/pupsum[NB_profile_pts-1];
        else
            F1 = 0.0;
        F2 = pupsum[ii]/pupsum[NB_profile_pts-1];
        fluxdens  = (F2-F1)/(2.0*PI*ii);
        /* compute pseudo radii psr0 and psr1 */
        psr0 = sqrt(1.0*ii*ii-1.0*ii);
        psr1 = sqrt(1.0*ii*ii+1.0*ii);
        piaaconfr1[i] = PIAAextfactor1*sqrt((F0-F1)/fluxdens/PI+psr0*psr0)/(NB_profile_pts-1);
    }




    piaaconfNBpoints = NBpoints;

    free(pup0sum);
    free(pupsum);
    free(pup);






    printf("compute r1 as a function of r0\n");
    fflush(stdout);

    /* compute r1 as a function of r0 */
    piaaconfr1fr0 = (double*) malloc(sizeof(double)*piaaconfNBr1fr0_pts);
    cntarray = (long*) malloc(sizeof(long)*piaaconfNBr1fr0_pts);

    for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
    {
        piaaconfr1fr0[ii] = 0.0;
        cntarray[ii] = 0;
    }

    for(ii=0; ii<NBpoints; ii++)
    {
        jj = (long) ((piaaconfr0[ii]/PIAAextfactor0)*piaaconfNBr1fr0_pts+0.5);
        if((jj>-1)&&(jj<piaaconfNBr1fr0_pts))
        {
            piaaconfr1fr0[jj] += piaaconfr1[ii];
            cntarray[jj]++;
        }
    }

    for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
        piaaconfr1fr0[ii] /= cntarray[ii];


    for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
        if(cntarray[ii]==0)
        {
            ii0 = ii;
            ii1 = ii;
            while(cntarray[ii0]==0)
                ii0--;
            while(cntarray[ii1]==0)
                ii1++;
            x = (1.0*ii-ii0)/(ii1-ii0);
            piaaconfr1fr0[ii] = (1.0-x)*piaaconfr1fr0[ii0] + x*piaaconfr1fr0[ii1];
        }


    if(file_write == 1)
    {
        fpout = fopen("PIAA_r0r1.txt","w");
        for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
            fprintf(fpout, "%18.16f %18.16f\n", PIAAextfactor0*ii/piaaconfNBr1fr0_pts, piaaconfr1fr0[ii]);
        fclose(fpout);
    }

    free(cntarray);


    printf("Compute PIAA mirror shapes\n");
    piaaconfbeamradius = 0.01;
    piaaconfdist = 1.0;
    piaaconfM0 = (double*) malloc(sizeof(double)*piaaconfNBr1fr0_pts);
    piaaconfM1 = (double*) malloc(sizeof(double)*piaaconfNBr1fr0_pts);

    piaaconfM0[0] = 0.0;
    piaaconfM1[0] = piaaconfdist;


    for(i=0; i<piaaconfNBr1fr0_pts-1; i++)
    {
        r0c = PIAAextfactor0*i/piaaconfNBr1fr0_pts;
        r1c = piaaconfr1fr0[i];
        dx = (r0c-r1c)*piaaconfbeamradius;
        dz = piaaconfM1[i]-piaaconfM0[i];
        dist = sqrt(dx*dx+dz*dz);
        y3 = dist - dz;
        if(fabs(dx)>0.000000001)
            slope = y3/dx;
        else
            slope = 0.0;
        r0n = PIAAextfactor0*(i+1)/piaaconfNBr1fr0_pts;
        piaaconfM0[i+1] = piaaconfM0[i] + slope*(r0n-r0c)*piaaconfbeamradius;

        if(fabs(dx)>0.000000001)
            slope = y3/dx;
        else
            slope = 0.0;
        piaaconfM1[i+1] = piaaconfM1[i] + slope*(piaaconfr1fr0[i+1]-r1c)*piaaconfbeamradius;
    }

    fp = fopen("PIAA_Mshapes.txt", "w");
    for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
        fprintf(fp, "%18.16f %18.16f %18.16f %18.16f\n", PIAAextfactor0*ii/piaaconfNBr1fr0_pts*piaaconfbeamradius, piaaconfM0[ii], piaaconfr1fr0[ii]*piaaconfbeamradius, piaaconfM1[ii]);
    fclose(fp);



    printf("compute r0 as a function of r1\n");
    fflush(stdout);
    /* compute r0 as a function of r1 */
    piaaconfr0fr1 = (double*) malloc(sizeof(double)*piaaconfNBr1fr0_pts);
    cntarray = (long*) malloc(sizeof(long)*piaaconfNBr1fr0_pts);

    for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
    {
        piaaconfr0fr1[ii] = 0.0;
        cntarray[ii] = 0;
    }

    for(ii=0; ii<NBpoints; ii++)
    {
        jj = (long) ((piaaconfr1[ii]/PIAAextfactor1)*piaaconfNBr1fr0_pts+0.5);
        if((jj>-1)&&(jj<piaaconfNBr1fr0_pts))
        {
            piaaconfr0fr1[jj] += piaaconfr0[ii];
            cntarray[jj]++;
        }
    }


    for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
        piaaconfr0fr1[ii] /= cntarray[ii];

    for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
        if(cntarray[ii]==0)
        {
            ii0 = ii;
            ii1 = ii;
            while(cntarray[ii0]==0)
                ii0--;
            while(cntarray[ii1]==0)
                ii1++;
            x = (1.0*ii-ii0)/(ii1-ii0);
            piaaconfr0fr1[ii] = (1.0-x)*piaaconfr0fr1[ii0] + x*piaaconfr0fr1[ii1];
        }

    if(file_write == 1)
    {
        fpout = fopen("PIAA_r1r0.txt", "w");
        for(ii=0; ii<piaaconfNBr1fr0_pts; ii++)
            fprintf(fpout, "%18.16f %18.16f\n", PIAAextfactor1*ii/piaaconfNBr1fr0_pts, piaaconfr0fr1[ii]);
        fclose(fpout);
    }

    free(cntarray);

    printf("init PIAA is done\n");
    fflush(stdout);

    return(NBpoints);
}


int coronagraphs_free_PIAA()
{

    free(fitapo_a);
    free(fitapo_b);
    free(fitapo_c);
    free(fitapo_c1);

    return(0);
}



int coronagraph_telescope_pupil_Subaru_inside1(double x, double y)
{
    int inside;
    double r,PA;
    double centobs = 0.285;
    double SPIDER_offset = -0.11;
    double SPIDER_width = 0.0275; /* fraction of pupil diameter */
    long nbarms = 4;
    double x1[4];
    double y1[4];
    double x2[4];
    double y2[4];
    double eps=1e-20;
    long arm;
    double tmpx,tmpy,x12,x13,y12,y13,asq,bsq,absq,dist;

    inside = 1;
    r = sqrt(x*x+y*y);
    if(r>1.0)
        inside = 0;
    if(r<centobs)
        inside = 0;

    PA=0.775;
    x1[0] = -sin(PA)*SPIDER_offset;
    x2[0] = cos(PA);
    y1[0] = cos(PA)*SPIDER_offset;
    y2[0] = sin(PA);

    x1[1] = x1[0];
    x2[1] = x2[0];
    y1[1] = -y1[0];
    y2[1] = -y2[0];

    x1[2] = -x1[0];
    x2[2] = -x2[0];
    y1[2] = y1[0];
    y2[2] = y2[0];

    x1[3] = -x1[0];
    x2[3] = -x2[0];
    y1[3] = -y1[0];
    y2[3] = -y2[0];

    PA = 0.0;
    for(arm=0; arm<nbarms; arm++)
    {
        tmpx = x1[arm]*cos(PA)+y1[arm]*sin(PA);
        tmpy = -x1[arm]*sin(PA)+y1[arm]*cos(PA);
        x1[arm] = tmpx;
        y1[arm] = tmpy;
        tmpx = x2[arm]*cos(PA)+y2[arm]*sin(PA);
        tmpy = -x2[arm]*sin(PA)+y2[arm]*cos(PA);
        x2[arm] = tmpx;
        y2[arm] = tmpy;
    }

    for(arm=0; arm<nbarms; arm++)
    {
        asq = (x1[arm]-x2[arm])*(x1[arm]-x2[arm])+(y1[arm]-y2[arm])*(y1[arm]-y2[arm]);
        x12 = x2[arm]-x1[arm];
        y12 = y2[arm]-y1[arm];


        bsq = (x1[arm]-x)*(x1[arm]-x)+(y1[arm]-y)*(y1[arm]-y);
        x13 = x-x1[arm];
        y13 = y-y1[arm];
        absq = (x12*x13+y12*y13)*(x12*x13+y12*y13);
        if(bsq-absq/asq<eps)
            dist = 0.0;
        else
            dist = sqrt(bsq-absq/asq);

        if(x12*x13+y12*y13>0)
        {
            if(dist<SPIDER_width)
                inside = 0;
        }
    }

    return(inside);
}


int coronagraphs_make_SUBARU_pupil()
{
    long ii,jj,i;
    // double epsilon = 0.0000000001;;
    long Ssize = 4096;
    double rad = 201.0;
    long ID;
    double x,y;
    //  double PA,r0,r1;
    //  long size = 1024;

    PIAACENTOBS = 0.278;
    // APLC_FPMASKsize = 1.6;
    //  PIAAAPO_NBPOINTS = 80;
    //  sprintf(PIAAAPO_FNAME,"APLCapo_%.3f.%ld.ref.prof",APLC_FPMASKsize,size);

    //  coronagraph_init_PIAA();
    ID = create_2Dimage_ID("subpup",Ssize,Ssize);
    for(ii=0; ii<Ssize; ii++)
        for(jj=0; jj<Ssize; jj++)
        {
            x = 1.0*ii-Ssize/2;
            y = 1.0*jj-Ssize/2;
            x /= rad;
            y /= rad;
            //PA = atan2(y,x);
            //r1 = sqrt(x*x+y*y);

            /*	if(r1<1.0)
              {
                if(r1<epsilon)
                  {
            	i = 1;
            	r0 = 1.0*i/piaaconfNBr1fr0_pts;
                  }
                else
                  {
            	i = (long) (r1*piaaconfNBr1fr0_pts+0.5);
                  }

                if(i>piaaconfNBr1fr0_pts-2)
                  i = piaaconfNBr1fr0_pts-2;
                r0 = piaaconfr0fr1[i];
              }
            else
              r0 = r1;
            */
            /*	r0 = r1;*/
            //x = r0*cos(PA);
            //y = r0*sin(PA);
            data.image[ID].array.F[jj*Ssize+ii] = coronagraph_telescope_pupil_Subaru_inside1(x,y);
        }

    return(0);
}




// compute 1st order PIAA remapping perturbation as a function of Zernike values
// shape is applied on PIAA optic #2
// ratio: distance between optics in pupil radius unit
// z displacement is in pupil radius unit
int coronagraph_PIAAperturbation(double *zarray, long *zindex, long NBzern, double ratio)
{
    double trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    long IDzmap; // M1 sag map
    long IDdx; // x displacement
    long IDdy; // y displacement
    long IDi; // intensity multiplication
    long ii, jj;
    double eps, epst;

    double x0, y0, r0, PA0;

    double p0v0, p0v1, p0v2;
    double p1v0, p1v1, p1v2;
    double p2v0, p2v1, p2v2;

    double p0x0, p0y0, p0r0, p0PA0;
    double p0x1, p0y1, p0r1, p0PA1;
    double p0x2, p0y2, p0r2, p0PA2;
    double p0dx, p0dy;
    double p0x0p, p0y0p;

    double p1x0, p1y0, p1r0, p1PA0;
    double p1x1, p1y1, p1r1, p1PA1;
    double p1x2, p1y2, p1r2, p1PA2;
    double p1dx, p1dy;
    double p1x0p, p1y0p;

    double p2x0, p2y0, p2r0, p2PA0;
    double p2x1, p2y1, p2r1, p2PA1;
    double p2x2, p2y2, p2r2, p2PA2;
    double p2dx, p2dy;
    double p2x0p, p2y0p;


    double Tpre, Tpost;


    long k;
    long size = CORONAGRAPHS_ARRAYSIZE;
    double dx, dy, dx1, dx2, dy1, dy2, ddx, ddy;

    double x1, y1, x2, y2;


    eps = 1.0e-4; // step to compute derivative [pupil radius]
    epst = 1.0e-3; // triangle radius

    printf("RADIUS = %f\n", trad_pix);
    IDzmap = create_2Dimage_ID("PIAApert_z", size, size);
    IDdx = create_2Dimage_ID("PIAApert_dx", size, size);
    IDdy = create_2Dimage_ID("PIAApert_dy", size, size);
    IDi = create_2Dimage_ID("PIAApert_i", size, size);

    zernike_init();

    //   2
    // 3 0 1
    //   4

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {

            x0 = 1.0*(ii-size/2)/trad_pix;
            y0 = 1.0*(jj-size/2)/trad_pix;
            r0 = sqrt(x0*x0 + y0*y0);
            PA0 = atan2(y0,x0);


            // defines 3 points and measures where they land

            if(r0< 1.1)
            {
                // point #0

                p0x0 = x0 + epst*cos(0.0);
                p0y0 = y0 + epst*sin(0.0);
                p0r0 = sqrt(p0x0*p0x0 + p0y0*p0y0);
                p0PA0 = atan2(p0y0, p0x0);

                p0x1 = p0x0 + eps;
                p0y1 = p0y0;
                p0r1 = sqrt(p0x1*p0x1 + p0y1*p0y1);
                p0PA1 = atan2(p0y1, p0x1);

                p0x2 = p0x0;
                p0y2 = p0y0 + eps;
                p0r2 = sqrt(p0x2*p0x2 + p0y2*p0y2);
                p0PA2 = atan2(p0y2, p0x2);

                p0v0 = 0.0;
                p0v1 = 0.0;
                p0v2 = 0.0;
                for(k=0; k<NBzern; k++)
                {
                    p0v0 += zarray[k]*Zernike_value(zindex[k], p0r0, p0PA0);
                    p0v1 += zarray[k]*Zernike_value(zindex[k], p0r1, p0PA1);
                    p0v2 += zarray[k]*Zernike_value(zindex[k], p0r2, p0PA2);
                }
                p0dx = -2.0*ratio * (p0v1-p0v0)/(eps); // second term is unitless, result is in pupil radius
                p0dy = -2.0*ratio * (p0v2-p0v0)/(eps);




                // point #1

                p1x0 = x0 + epst*cos(2.0*M_PI/3.0);
                p1y0 = y0 + epst*sin(2.0*M_PI/3.0);
                p1r0 = sqrt(p1x0*p1x0 + p1y0*p1y0);
                p1PA0 = atan2(p1y0, p1x0);

                p1x1 = p1x0 + eps;
                p1y1 = p1y0;
                p1r1 = sqrt(p1x1*p1x1 + p1y1*p1y1);
                p1PA1 = atan2(p1y1, p1x1);

                p1x2 = p1x0;
                p1y2 = p1y0 + eps;
                p1r2 = sqrt(p1x2*p1x2 + p1y2*p1y2);
                p1PA2 = atan2(p1y2, p1x2);

                p1v0 = 0.0;
                p1v1 = 0.0;
                p1v2 = 0.0;

                for(k=0; k<NBzern; k++)
                {
                    p1v0 += zarray[k]*Zernike_value(zindex[k], p1r0, p1PA0);
                    p1v1 += zarray[k]*Zernike_value(zindex[k], p1r1, p1PA1);
                    p1v2 += zarray[k]*Zernike_value(zindex[k], p1r2, p1PA2);
                }
                p1dx = -2.0*ratio * (p1v1-p1v0)/(eps); // second term is unitless, result is in pupil radius
                p1dy = -2.0*ratio * (p1v2-p1v0)/(eps);



                // point #2

                p2x0 = x0 + epst*cos(4.0*M_PI/3.0);
                p2y0 = y0 + epst*sin(4.0*M_PI/3.0);
                p2r0 = sqrt(p2x0*p2x0 + p2y0*p2y0);
                p2PA0 = atan2(p2y0, p2x0);

                p2x1 = p2x0 + eps;
                p2y1 = p2y0;
                p2r1 = sqrt(p2x1*p2x1 + p2y1*p2y1);
                p2PA1 = atan2(p2y1, p2x1);

                p2x2 = p2x0;
                p2y2 = p2y0 + eps;
                p2r2 = sqrt(p2x2*p2x2 + p2y2*p2y2);
                p2PA2 = atan2(p2y2, p2x2);

                p2v0 = 0.0;
                p2v1 = 0.0;
                p2v2 = 0.0;

                for(k=0; k<NBzern; k++)
                {
                    p2v0 += zarray[k]*Zernike_value(zindex[k], p2r0, p2PA0);
                    p2v1 += zarray[k]*Zernike_value(zindex[k], p2r1, p2PA1);
                    p2v2 += zarray[k]*Zernike_value(zindex[k], p2r2, p2PA2);
                }
                p2dx = -2.0*ratio * (p2v1-p2v0)/(eps); // second term is unitless, result is in pupil radius
                p2dy = -2.0*ratio * (p2v2-p2v0)/(eps);



                data.image[IDzmap].array.F[jj*size+ii] = (p0v0+p1v0+p2v0)/3.0;


                data.image[IDdx].array.F[jj*size+ii] = (p0dx+p1dx+p2dx)/3.0;
                data.image[IDdy].array.F[jj*size+ii] = (p0dy+p1dy+p2dy)/3.0;


                p0x0p = p0x0 + p0dx;
                p0y0p = p0y0 + p0dy;

                p1x0p = p1x0 + p1dx;
                p1y0p = p1y0 + p1dy;

                p2x0p = p2x0 + p2dx;
                p2y0p = p2y0 + p2dy;

                Tpre = 0.5*fabs((p0x0-p2x0)*(p1y0-p0y0)-(p0x0-p1x0)*(p2y0-p0y0));
                Tpost = 0.5*fabs((p0x0p-p2x0p)*(p1y0p-p0y0p)-(p0x0p-p1x0p)*(p2y0p-p0y0p));

                x1 = p1x0p - p0x0p;
                y1 = p1y0p - p0y0p;
                x2 = p2x0p - p0x0p;
                y2 = p2y0p - p0y0p;

                data.image[IDi].array.F[jj*size+ii] = Tpre/Tpost;
            }
        }
    save_fl_fits("PIAApert_z", "!PIAApert_z.fits");
    save_fl_fits("PIAApert_dx", "!PIAApert_dx.fits");
    save_fl_fits("PIAApert_dy", "!PIAApert_dy.fits");
    save_fl_fits("PIAApert_i", "!PIAApert_i.fits");

    return(0);
}




int coronagraphs_PIAA_apodize_beam(const char *ampl1, const char *opd1, const char *ampl2, const char *opd2)
{
    /* apodize a beam with PIAA */
    double PA;
    double r1,r2;
    long ii,jj,ii2,jj2, iin, jjn;
    long i;
    double x,y;
    long IDa1,IDo1,IDa2,IDo2;
    double factor,r;
    long *cntarray;
    double total,totalcnt;
    double epsilon = 0.0000000001;
    long size = CORONAGRAPHS_ARRAYSIZE;
    double trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    double t,u;
    double v00, v01, v10, v11;
    double eps;
    double coeff;

    int USE_2DPROL = 1;  /* if set to 1, the 2D apodized pupil will be used to make the apodization */
    int ID_2DPROL;
    char fname_2DPROL[400];
    double totalProl = 0.0;
    long totalProlcnt = 0;

    long size2;


    //  long IDr;
    //  FILE *fptest;

    size2 = size*size;
    //  save_fl_fits(ampl1, "!test_ampl1.fits");
    //exit(0);

    if(USE_2DPROL == 1)
    {
        //      sprintf(fname_2DPROL, "%s/APLCapo/raw/%s", CORONAGRAPHSDATALOCAL, PIAAAPODIZE_2DAPOFNAME);
        sprintf(fname_2DPROL, "apo2Dimt.fits"); // generated by init_PIAA

        ID_2DPROL = image_ID("prol2d");
        if(ID_2DPROL == -1)
            ID_2DPROL = load_fits(fname_2DPROL, "prol2d", 1);
        //     save_fl_fits("prol2d","!test_prol2d.fits");
        // exit(0);


        for(ii=0; ii<size2; ii++)
        {
            if( data.image[ID_2DPROL].array.F[ii]>0.00001)
            {
                totalProl += data.image[ID_2DPROL].array.F[ii]*data.image[ID_2DPROL].array.F[ii];
                totalProlcnt ++;
            }
        }

        for(ii=0; ii<size2; ii++)
            data.image[ID_2DPROL].array.F[ii] /= sqrt(totalProl/totalProlcnt);
    }


    IDa1 = image_ID(ampl1);
    IDo1 = image_ID(opd1);

    IDa2 = create_2Dimage_ID(ampl2,size,size);
    IDo2 = create_2Dimage_ID(opd2,size,size);

    cntarray = (long*) malloc(sizeof(long)*size*size);

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
            cntarray[jj*size+ii] = 0;

    total = 0.0;
    totalcnt = 0.0;

    // fptest = fopen("testfile.txt","w");

    //  IDr = create_2Dimage_ID("rtest", size, size);

    for(ii2=0; ii2<size; ii2++)
        for(jj2=0; jj2<size; jj2++)
        {
            x = 1.0*(ii2-size/2)/trad_pix;
            y = 1.0*(jj2-size/2)/trad_pix;
            r2 = sqrt(x*x+y*y);
            PA=atan2(y,x);
            r = r2;

            if(r2<1.0)
            {
                if(piaaconfdirection==0)
                {
                    if(r2<epsilon)
                    {
                        i = 1;
                        r2 = 1.0*i/piaaconfNBr1fr0_pts*PIAAextfactor1;
                    }
                    else
                    {
                        /* compute r1 */
                        i = (long) (r2*piaaconfNBr1fr0_pts/PIAAextfactor1+0.5);
                    }

                    if(i>piaaconfNBr1fr0_pts-2)
                        i = piaaconfNBr1fr0_pts-2;

                    r1 = piaaconfr0fr1[i];
                    if(i<piaaconfNBr1fr0_pts-2)
                        factor = piaaconfr0fr1[i+2]-piaaconfr0fr1[i];
                    else
                        factor = piaaconfr0fr1[i]-piaaconfr0fr1[i-2];
                }
                else
                {
                    if(r2<epsilon)
                    {
                        i = 1;
                        r2 = 1.0*i/piaaconfNBr1fr0_pts*PIAAextfactor0;
                    }
                    else
                    {
                        /* compute r1 */
                        i = (long) (r2*piaaconfNBr1fr0_pts/PIAAextfactor0+0.5);
                    }

                    if(i>piaaconfNBr1fr0_pts-2)
                        i = piaaconfNBr1fr0_pts-2;

                    r1 = piaaconfr1fr0[i];
                    //	fprintf(fptest, "%ld %f\n", i, r1);
                    if(i<piaaconfNBr1fr0_pts-2)
                        factor = piaaconfr1fr0[i+2]-piaaconfr1fr0[i];
                    else
                        factor = piaaconfr1fr0[i]-piaaconfr1fr0[i-2];
                }

                factor *= r1/r2;

                //	    data.image[IDr].array.F[jj2*size+ii2] = r1;

                ii = (long) (0.5*size + r1*trad_pix*cos(PA));
                t = 0.5*size + r1*trad_pix*cos(PA) - ii;
                jj = (long) (0.5*size + r1*trad_pix*sin(PA));
                u = 0.5*size + r1*trad_pix*sin(PA) - jj;

                iin = (long) (0.5*size + r1*trad_pix*cos(PA) + 0.5);
                jjn = (long) (0.5*size + r1*trad_pix*sin(PA) + 0.5);

                //	    printf("PIAAAPO2D = %f\n", PIAAAPO2D[size/2*size+size/2] );

                if((ii>-1)&&(jj>-1)&&(ii<size-1)&&(jj<size-1))
                {
                    /*	data.image[IDa2].array.F[jj2*size+ii2] = data.image[IDa1].array.F[jj*size+ii]*factor;*/
                    i = (long) (r*PIAAAPO_NBPOINTS);
                    if(i<PIAAAPO_NBPOINTS)
                    {
                        if(piaaconfdirection==0)
                        {
                            if(USE_2DPROL == 1)
                                data.image[IDa2].array.F[jj2*size+ii2] = data.image[ID_2DPROL].array.F[jj2*size+ii2] * data.image[IDa1].array.F[jjn*size+iin];
                            else
                                data.image[IDa2].array.F[jj2*size+ii2] = data.image[IDa1].array.F[jj*size+ii]*PIAAAPO2D[jj2*size+ii2]*PIAA_HYBRID_CPAAPO[i];

                        }
                        else
                            data.image[IDa2].array.F[jj2*size+ii2] = data.image[IDa1].array.F[jj*size+ii]/(PIAAAPO2D[jj*size+ii]+epsilon);
                        if(r>PIAACPPMASKRAD1)
                            data.image[IDa2].array.F[jj2*size+ii2] = 0.0;

                        if(1) // bilinear interpolation
                        {
                            data.image[IDo2].array.F[jj2*size+ii2] = (1.0-t)*(1.0-u)*data.image[IDo1].array.F[jj*size+ii]+t*(1.0-u)*data.image[IDo1].array.F[jj*size+ii+1]+t*u*data.image[IDo1].array.F[(jj+1)*size+ii+1]+(1.0-t)*u*data.image[IDo1].array.F[(jj+1)*size+ii];
                        }
                        else
                        {
                            v00 = data.image[IDo1].array.F[jj*size+ii];
                            v01 = data.image[IDo1].array.F[jj*size+ii+1];
                            v10 = data.image[IDo1].array.F[(jj+1)*size+ii];
                            v11 = data.image[IDo1].array.F[(jj+1)*size+ii+1];
                            data.image[IDo2].array.F[jj2*size+ii2] = 0.0;
                            coeff = 0.0;
                            eps = 1.0e-5;
                            if(v00>eps)
                            {
                                data.image[IDo2].array.F[jj2*size+ii2] += (1.0-t)*(1.0-u)*v00;
                                coeff += (1.0-t)*(1.0-u);
                            }
                            if(v01>eps)
                            {
                                data.image[IDo2].array.F[jj2*size+ii2] += t*(1.0-u)*v01;
                                coeff += t*(1.0-u);
                            }
                            if(v10>eps)
                            {
                                data.image[IDo2].array.F[jj2*size+ii2] += (1.0-t)*u*v10;
                                coeff += (1.0-t)*u;
                            }
                            if(v11>eps)
                            {
                                data.image[IDo2].array.F[jj2*size+ii2] += t*u*v11;
                                coeff += t*u;
                            }
                            if(coeff>eps)
                                data.image[IDo2].array.F[jj2*size+ii2] /= coeff;
                        }
                        total += data.image[IDo2].array.F[jj2*size+ii2];
                        totalcnt += 1.0;
                    }
                }
                else
                    data.image[IDo2].array.F[jj2*size+ii2] = 0.0;
            }
        }
    //  fclose(fptest);
    total /= totalcnt;

    /*  for(ii2=0;ii2<size;ii2++)
      for(jj2=0;jj2<size;jj2++)
        {
    x = 1.0*(ii2-size/2)/trad_pix;
    y = 1.0*(jj2-size/2)/trad_pix;
    r2 = sqrt(x*x+y*y);
    PA=atan2(y,x);

    if(r2<1.0)
      data.image[IDo2].array.F[jj2*size+ii2] -= total;
      }*/

    free(cntarray);

    // TEST

    /*  save_fl_fits(ampl1, "!test_ampl1.fits");
    save_fl_fits(ampl2, "!test_ampl2.fits");
    save_fl_fits(opd2, "!test_opd2.fits");
    save_fl_fits("rtest", "!test_r.fits");

    delete_image_ID("rtest");
    exit(0);
    */

    return(0);
}




int coronagraph_init_CPA()
{
    FILE *fp;
    int result = 0;
    long i;
    long size = CORONAGRAPHS_ARRAYSIZE;
    double x,y,r;
    long tmpl;
    double tmpf;
    long ID;
    long factor = 9;
    long ii,jj,ii1,jj1;
    /*  double trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;*/
    double ri,value,alpha;
    long ri1;
    char fname[200];
    char fname1[200];


    sprintf(fname,"%s/%s", CORONAGRAPHSDATALOCAL,CPAAPO_FNAME);
    if((fp=fopen(fname,"r"))==NULL)
    {
        printf("ERROR : cannot open file \"%s\"\n",fname);
        exit(0);
    }
    for(i=0; i<CPAAPO_NBPOINTS; i++)
    {
        result = fscanf(fp,"%ld %lf\n", &tmpl, &tmpf);
        CPAAPO[i] = tmpf;
    }

    fclose(fp);

    for(i=1; i<CPAAPO_NBPOINTS; i++)
        CPAAPO[i] /=  CPAAPO[0];
    CPAAPO[0] = 1.0;

    ID = create_2Dimage_ID("cpapupref",size,size);
    for(ii=0; ii<factor*size; ii++)
        for(jj=0; jj<factor*size; jj++)
        {
            ii1 = (long) ((1.0*ii-0.5)/factor+0.5);
            jj1 = (long) ((1.0*jj-0.5)/factor+0.5);
            x = 1.0/factor*(ii-factor*size/2);
            y = 1.0/factor*(jj-factor*size/2);
            r = sqrt(x*x+y*y); /* in pixel */
            r /= CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
            ri = r*CPAAPO_NBPOINTS;
            ri1 = (long) (ri);
            alpha = ri-ri1;
            if(ri<CPAAPO_NBPOINTS-1)
                value = (1.0-alpha)*CPAAPO[ri1]+alpha*CPAAPO[ri1+1];
            else
                value = CPAAPO[CPAAPO_NBPOINTS-1];
            if((ii1>0)&&(jj1>0)&&(ii1<size)&&(jj1<size))
                data.image[ID].array.F[jj1*size+ii1] += value/factor/factor;
        }
    sprintf(fname1,"!%s/cpapup.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    save_fl_fits("cpapupref",fname1);

    return(0);
}



int coronagraph_init_ODC()
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    int factor = 4;
    long sizeODC1 = (size/2)*factor;
    double *ODC1m;
    double x;
    long ii,ii1,jj;
    double max=0.0;
    long ID;
    char fname1[200];
    long index;
    double tmp;

    ODC1m = (double*) malloc(sizeof(double)*sizeODC1);

    for(ii=0; ii<sizeODC1; ii++)
    {
        x = CORONAGRAPHS_PIXSCALE*ii;
        ODC1m[ii] = x*exp(-x*x/ODC_GAUSS/ODC_GAUSS);
        if(ODC1m[ii]>max)
        {
            max = ODC1m[ii];
            printf("x = %f\n",x);
        }
    }

    for(ii=0; ii<sizeODC1; ii++)
        ODC1m[ii] /= max;


    ID = create_2Dimage_ID("odcmref",size,size);
    for(ii=0; ii<factor*size; ii++)
    {
        ii1 = (long) ((1.0*ii-0.5)/factor+0.5);
        x = 1.0/factor*(ii-factor*size/2);
        if(ii>factor*size/2)
        {
            index = ii-factor*size/2;
            tmp = ODC1m[index];
        }
        else
        {
            index = factor*size/2-ii;
            tmp = -ODC1m[index];
        }

        if((ii1>0)&&(ii1<size))
            for(jj=0; jj<size; jj++)
                data.image[ID].array.F[jj*size+ii1] += tmp/factor;
    }

    free(ODC1m);

    sprintf(fname1,"!%s/odcm.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    save_fl_fits("odcmref",fname1);

    return(0);
}



int coronagraph_init_BL8()
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    double x,y,r;
    long i;
    long ii,jj;
    long ii1,jj1;
    double tmp1,tmp2,tmp3;
    double x2,x3;
    double max;
    long ID;
    long factor = 9;
    char fname1[200];

    tmp1 = (BL8MASK_l-BL8MASK_m)/BL8MASK_l;
    max = 0.0;
    BL8MASK[0] = 0.0;
    for(i=1; i<BL8MASK_NBSTEP; i++)
    {
        x = BL8MASK_STEP*i;
        x2 = PI*x*BL8MASK_eps/BL8MASK_l;
        tmp2 = pow(sin(x2)/x2,BL8MASK_l);
        x3 = PI*x*BL8MASK_eps/BL8MASK_m;
        tmp3 = BL8MASK_m/BL8MASK_l*pow(sin(x3)/x3,BL8MASK_m);
        BL8MASK[i] = tmp1-tmp2+tmp3;
        if(BL8MASK[i]>max)
            max = BL8MASK[i];
    }
    for(i=0; i<BL8MASK_NBSTEP; i++)
        BL8MASK[i] /= max;


    ID = create_2Dimage_ID("bl8mref",size,size);
    for(ii=0; ii<factor*size; ii++)
        for(jj=0; jj<factor*size; jj++)
        {
            ii1 = (long) ((1.0*ii-0.5)/factor+0.5);
            jj1 = (long) ((1.0*jj-0.5)/factor+0.5);
            x = 1.0/factor*(ii-factor*size/2);
            y = 1.0/factor*(jj-factor*size/2);
            r = sqrt(x*x+y*y);
            r *= CORONAGRAPHS_PIXSCALE;
            i = (long) (r/BL8MASK_STEP);
            if((ii1>0)&&(jj1>0)&&(ii1<size)&&(jj1<size))
            {
                if(i>BL8MASK_NBSTEP)
                    data.image[ID].array.F[jj1*size+ii1] += 1.0/factor/factor;
                else
                    data.image[ID].array.F[jj1*size+ii1] += BL8MASK[i]/factor/factor;
            }
        }
    sprintf(fname1,"!%s/bl8m.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    save_fl_fits("bl8mref",fname1);

    ID = create_2Dimage_ID("bl8mrefl",size,size);
    for(ii=0; ii<factor*size; ii++)
        for(jj=0; jj<factor*size; jj++)
        {
            ii1 = (long) ((1.0*ii-0.5)/factor+0.5);
            jj1 = (long) ((1.0*jj-0.5)/factor+0.5);
            x = 1.0/factor*(ii-factor*size/2);
            y = 1.0/factor*(jj-factor*size/2);
            r = sqrt(x*x);
            r *= CORONAGRAPHS_PIXSCALE;
            i = (long) (r/BL8MASK_STEP);
            if((ii1>0)&&(jj1>0)&&(ii1<size)&&(jj1<size))
            {
                if(i>BL8MASK_NBSTEP)
                    data.image[ID].array.F[jj1*size+ii1] += 1.0/factor/factor;
                else
                    data.image[ID].array.F[jj1*size+ii1] += BL8MASK[i]/factor/factor;
            }
        }
    sprintf(fname1,"!%s/bl8ml.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    save_fl_fits("bl8mrefl",fname1);

    return(0);
}


int coronagraph_init_BL4()
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    double x,y,r,rx;
    long ii,jj;
    long ii1,jj1;
    double tmp1,tmp2;
    long ID;
    long factor = 9;
    char fname1[200];


    ID = create_2Dimage_ID("bl4mref",size,size);
    for(ii=0; ii<factor*size; ii++)
        for(jj=0; jj<factor*size; jj++)
        {
            ii1 = (long) ((1.0*ii-0.5)/factor+0.5);
            jj1 = (long) ((1.0*jj-0.5)/factor+0.5);
            x = 1.0/factor*(ii-factor*size/2);
            y = 1.0/factor*(jj-factor*size/2);
            rx = sqrt(x*x);
            r = sqrt(x*x+y*y);
            rx *= CORONAGRAPHS_PIXSCALE;
            r *= CORONAGRAPHS_PIXSCALE;
            tmp1 = sin(rx*PI*BL4MASK_eps/2);
            tmp2 = exp(-(r/50.0)*(r/50.0)/2);
            if((ii1>0)&&(jj1>0)&&(ii1<size)&&(jj1<size))
                data.image[ID].array.F[jj1*size+ii1] += tmp1*tmp1*tmp2/factor/factor;
        }
    sprintf(fname1,"!%s/bl4m.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    save_fl_fits("bl4mref",fname1);

    return(0);
}


int coronagraph_init_RRPM()
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long size2 = size*size;
    double x,y,r;
    long ii,jj;
    long ID;
    double trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    char fname1[200];

    ID = make_subpixdisk("rrpmreffm",size,size,size/2,size/2,RRPM_RADIUS/CORONAGRAPHS_PIXSCALE);
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] = 1.0-2.0*data.image[ID].array.F[ii];
    sprintf(fname1,"!%s/rrpmref_fm.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    save_fl_fits("rrpmreffm",fname1);

    ID = create_2Dimage_ID("rrpmrefpm",size,size);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = 1.0*(ii-size/2);
            y = 1.0*(jj-size/2);
            r = sqrt(x*x+y*y)/trad_pix;
            data.image[ID].array.F[jj*size+ii] = 1.0+RRPM_P2*r*r+RRPM_P3*pow(r,3.0)+RRPM_P4*pow(r,4.0)+RRPM_P5*pow(r,5.0);
            if(data.image[ID].array.F[jj*size+ii]>1.0)
                data.image[ID].array.F[jj*size+ii] = 1.0;
        }
    sprintf(fname1,"!%s/rrpmref_pm.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    save_fl_fits("rrpmrefpm",fname1);

    return(0);
}


int coronagraph_init_OVC(long charge)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    double x,y,theta,re,im;
    long ii,jj;
    long i,j;
    long cnt;
    long ssamp=3;
    long ID1,ID2;
    char fname1[200];
    double eps=1e-8;
    char fnamea[200];
    char fnamep[200];

    printf("Initializing OVC charge %ld...\n",charge);
    fflush(stdout);

    sprintf(fnamea,"ovcreffma%ld",charge);
    sprintf(fnamep,"ovcreffmp%ld",charge);
    ID1 = create_2Dimage_ID(fnamep,size,size);
    ID2 = create_2Dimage_ID(fnamea,size,size);


    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            re = 0.0;
            im = 0.0;
            cnt = 0;
            for(i=-ssamp; i<ssamp+1; i++)
                for(j=-ssamp; j<ssamp+1; j++)
                {
                    x = 1.0*(ii-size/2)+0.5*i/(0.5+ssamp);
                    y = 1.0*(jj-size/2)+0.5*j/(0.5+ssamp);
                    theta = 1.0*charge*atan2(y,x);
                    if((x*x+y*y)>eps)
                    {
                        re += cos(theta);
                        im += sin(theta);
                    }
                    cnt ++;
                }
            re /= cnt;
            im /= cnt;
            data.image[ID2].array.F[jj*size+ii] = sqrt(re*re+im*im);
            data.image[ID1].array.F[jj*size+ii] = atan2(im,re);
        }

    sprintf(fname1,"!%s/ovcref_fmp%ld.ref.%ld",CORONAGRAPHSDATALOCAL,charge,size);
    save_fl_fits(fnamep,fname1);
    sprintf(fname1,"!%s/ovcref_fma%ld.ref.%ld",CORONAGRAPHSDATALOCAL,charge,size);
    save_fl_fits(fnamea,fname1);

    printf("done\n");
    fflush(stdout);

    return(0);
}



int coronagraph_simul_SHEAR4(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1,IDc;
    long ID;
    double trad_pix;
    long ii,jj;
    double a1,a2,a3,a4,p1,p2,p3,p4,r1,r2,r3,r4,i1,i2,i3,i4;
    long index1x,index1y,index1;
    long index2x,index2y,index2;
    long index3x,index3y,index3;
    long index4x,index4y,index4;
    double total;
    long size2=size*size;
    double epsilon = 1e-10;
    long iistart,iiend;
    long shearpix;

    /*  printf("Pixel scale = %f l/d per pixel\n",CORONAGRAPHS_PIXSCALE);*/

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    shearpix = (long) (SHEAR4_SHEAR*trad_pix);

    IDa1 = make_disk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    IDc = create_2DCimage_ID("pc2",size,size);

    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }


    iistart = (long) (0.5*size-trad_pix-shearpix);
    iiend = (long) (0.5*size+trad_pix+shearpix);
    if(iistart<1)
        iistart = 1;
    if(iiend>size-1)
        iiend = size-1;


    for(jj=iistart; jj<iiend; jj++)
        for(ii=iistart; ii<iiend; ii++)
        {
            index1x = ii-shearpix;
            index1y = jj+shearpix;
            index1 = index1y*size+index1x;

            index2x = ii+shearpix;
            index2y = jj+shearpix;
            index2 = index2y*size+index2x;

            index3x = ii-shearpix;
            index3y = jj-shearpix;
            index3 = index3y*size+index3x;

            index4x = ii+shearpix;
            index4y = jj-shearpix;
            index4 = index4y*size+index4x;

            a1 = data.image[IDa1].array.F[index1];
            p1 = data.image[IDp1].array.F[index1];
            r1 = a1*cos(p1);
            i1 = a1*sin(p1);

            a2 = data.image[IDa1].array.F[index2];
            p2 = data.image[IDp1].array.F[index2];
            r2 = a2*cos(p2);
            i2 = a2*sin(p2);

            a3 = data.image[IDa1].array.F[index3];
            p3 = data.image[IDp1].array.F[index3];
            r3 = a3*cos(p3);
            i3 = a3*sin(p3);

            a4 = data.image[IDa1].array.F[index4];
            p4 = data.image[IDp1].array.F[index4];
            r4 = a4*cos(p4);
            i4 = a4*sin(p4);

            if(a1*a2*a3*a4>epsilon)
            {
                data.image[IDc].array.CF[jj*size+ii].re = 0.25*(r1-r2-r3+r4);
                data.image[IDc].array.CF[jj*size+ii].im = 0.25*(i1-i2-i3+i4);
            }
        }

    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    mk_amph_from_complex("fc2", "fa2", "fp2", 0);
    arith_image_mult("fa2","fa2",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fa2");
    delete_image_ID("fp2");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2;

    return(0);
}


int coronagraph_simul_DICC(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1,IDc;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    long ID_DICC1,ID_DICCX,ID_DICCY,ID_DICCX2,ID_DICCY2,ID_DICCXY;
    long index;
    double x,y;
    double v1r,vxr,vyr,vx2r,vy2r,vxyr;
    double v1i,vxi,vyi,vx2i,vy2i,vxyi;
    double total1,totalx,totaly,totalx2,totaly2,totalxy;
    double amp,pha,re,im;
    long iter;
    long NBiter = 2;
    double totalx2s,totaly2s,totalx2y2;
    long IDaref;

    /*  printf("Pixel scale = %f l/d per pixel\n",CORONAGRAPHS_PIXSCALE);*/

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;

    IDaref = make_disk("paref",size,size,0.5*size,0.5*size,trad_pix);
    IDa1 = make_disk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    IDc = create_2DCimage_ID("pc2",size,size);

    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    ID_DICC1 = image_ID("dicc1");
    if(ID_DICC1==-1)
    {
        printf("Initializing modes...");
        fflush(stdout);
        ID_DICC1 = create_2Dimage_ID("dicc1",size,size);
        ID_DICCX = create_2Dimage_ID("diccx",size,size);
        ID_DICCY = create_2Dimage_ID("diccy",size,size);
        ID_DICCX2 = create_2Dimage_ID("diccx2",size,size);
        ID_DICCY2 = create_2Dimage_ID("diccy2",size,size);
        ID_DICCXY = create_2Dimage_ID("diccxy",size,size);
        total1 = 0.0;
        totalx = 0.0;
        totaly = 0.0;
        totalx2 = 0.0;
        totaly2 = 0.0;
        totalxy = 0.0;
        totalx2s = 0.0;
        totaly2s = 0.0;
        totalx2y2 = 0.0;

        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                if(data.image[IDaref].array.F[jj*size+ii] > 0.5)
                {
                    x = 1.0*(ii-size/2);
                    y = 1.0*(jj-size/2);
                    totalx2s += x*x;
                    totaly2s += y*y;
                }
            }
        totalx2s /= total;
        totaly2s /= total;

        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                if(data.image[IDaref].array.F[jj*size+ii] > 0.5)
                {
                    x = 1.0*(ii-size/2);
                    y = 1.0*(jj-size/2);
                    data.image[ID_DICC1].array.F[jj*size+ii] = 1.0;
                    total1 += 1.0;
                    data.image[ID_DICCX].array.F[jj*size+ii] = x;
                    totalx += x*x;
                    data.image[ID_DICCY].array.F[jj*size+ii] = y;
                    totaly += y*y;
                    data.image[ID_DICCX2].array.F[jj*size+ii] = x*x-totalx2s;
                    totalx2 += data.image[ID_DICCX2].array.F[jj*size+ii]*data.image[ID_DICCX2].array.F[jj*size+ii];
                    data.image[ID_DICCY2].array.F[jj*size+ii] = y*y-totaly2s;
                    totaly2 += data.image[ID_DICCY2].array.F[jj*size+ii]*data.image[ID_DICCY2].array.F[jj*size+ii];
                    data.image[ID_DICCXY].array.F[jj*size+ii] = x*y;
                    totalxy += x*y*x*y;
                }
            }
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                if(data.image[IDaref].array.F[jj*size+ii] > 0.5)
                {
                    data.image[ID_DICC1].array.F[jj*size+ii] /= sqrt(total1);
                    data.image[ID_DICCX].array.F[jj*size+ii] /= sqrt(totalx);
                    data.image[ID_DICCY].array.F[jj*size+ii] /= sqrt(totaly);
                    data.image[ID_DICCX2].array.F[jj*size+ii] /= sqrt(totalx2);
                    data.image[ID_DICCY2].array.F[jj*size+ii] /= sqrt(totaly2);
                    data.image[ID_DICCXY].array.F[jj*size+ii] /= sqrt(totalxy);
                }
            }

        totalx2y2 = 0.0;
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                if(data.image[IDaref].array.F[jj*size+ii] > 0.5)
                {
                    totalx2y2 += data.image[ID_DICCX2].array.F[jj*size+ii]*data.image[ID_DICCY2].array.F[jj*size+ii];
                }
            }
        totaly2 = 0.0;
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                if(data.image[IDaref].array.F[jj*size+ii] > 0.5)
                {
                    data.image[ID_DICCY2].array.F[jj*size+ii] -= totalx2y2*data.image[ID_DICCX2].array.F[jj*size+ii];
                    totaly2 += data.image[ID_DICCY2].array.F[jj*size+ii]*data.image[ID_DICCY2].array.F[jj*size+ii];
                }
            }

        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                if(data.image[IDaref].array.F[jj*size+ii] > 0.5)
                {
                    data.image[ID_DICCY2].array.F[jj*size+ii] /= sqrt(totaly2);
                }
            }


        /*  save_fl_fits("dicc1","!dicc1");
        save_fl_fits("diccx","!diccx");
        save_fl_fits("diccy","!diccy");
        save_fl_fits("diccx2","!diccx2");
        save_fl_fits("diccy2","!diccy2");
        save_fl_fits("diccxy","!diccxy");*/
        printf(" done\n");
        fflush(stdout);
    }
    ID_DICCX = image_ID("diccx");
    ID_DICCY = image_ID("diccy");
    ID_DICCX2 = image_ID("diccx2");
    ID_DICCY2 = image_ID("diccy2");
    ID_DICCXY = image_ID("diccxy");

    for(iter=0; iter<NBiter; iter++)
    {
        v1r = 0.0;
        vxr = 0.0;
        vyr = 0.0;
        vx2r = 0.0;
        vy2r = 0.0;
        vxyr = 0.0;
        v1i = 0.0;
        vxi = 0.0;
        vyi = 0.0;
        vx2i = 0.0;
        vy2i = 0.0;
        vxyi = 0.0;

        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                index = jj*size+ii;
                if(data.image[IDaref].array.F[index] > 0.5)
                {
                    amp = data.image[IDa1].array.F[index];
                    pha = data.image[IDp1].array.F[index];
                    re = amp*cos(pha);
                    im = amp*sin(pha);
                    v1r += re*data.image[ID_DICC1].array.F[index];
                    v1i += im*data.image[ID_DICC1].array.F[index];
                    vxr += re*data.image[ID_DICCX].array.F[index];
                    vxi += im*data.image[ID_DICCX].array.F[index];
                    vyr += re*data.image[ID_DICCY].array.F[index];
                    vyi += im*data.image[ID_DICCY].array.F[index];
                    vx2r += re*data.image[ID_DICCX2].array.F[index];
                    vx2i += im*data.image[ID_DICCX2].array.F[index];
                    vy2r += re*data.image[ID_DICCY2].array.F[index];
                    vy2i += im*data.image[ID_DICCY2].array.F[index];
                    vxyr += re*data.image[ID_DICCXY].array.F[index];
                    vxyi += im*data.image[ID_DICCXY].array.F[index];
                }
            }

        /*      printf("modes : %g %g %g %g %g %g\n",sqrt(v1r*v1r+v1i*v1i)/sqrt(total),sqrt(vxr*vxr+vxi*vxi)/sqrt(total),sqrt(vyr*vyr+vyi*vyi)/sqrt(total),sqrt(vx2r*vx2r+vx2i*vx2i)/sqrt(total),sqrt(vy2r*vy2r+vy2i*vy2i)/sqrt(total),sqrt(vxyr*vxyr+vxyi*vxyi)/sqrt(total));*/

        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                index = jj*size+ii;
                if(data.image[IDaref].array.F[index] > 0.5)
                {
                    amp = data.image[IDa1].array.F[index];
                    pha = data.image[IDp1].array.F[index];
                    re = amp*cos(pha);
                    im = amp*sin(pha);
                    re -= v1r*data.image[ID_DICC1].array.F[index];
                    im -= v1i*data.image[ID_DICC1].array.F[index];
                    re -= vxr*data.image[ID_DICCX].array.F[index];
                    im -= vxi*data.image[ID_DICCX].array.F[index];
                    re -= vyr*data.image[ID_DICCY].array.F[index];
                    im -= vyi*data.image[ID_DICCY].array.F[index];
                    re -= vx2r*data.image[ID_DICCX2].array.F[index];
                    im -= vx2i*data.image[ID_DICCX2].array.F[index];
                    re -= vy2r*data.image[ID_DICCY2].array.F[index];
                    im -= vy2i*data.image[ID_DICCY2].array.F[index];
                    re -= vxyr*data.image[ID_DICCXY].array.F[index];
                    im -= vxyi*data.image[ID_DICCXY].array.F[index];
                    data.image[IDa1].array.F[index] = sqrt(re*re+im*im);
                    data.image[IDp1].array.F[index] = atan2(im,re);

                }
            }
    }
    delete_image_ID("paref");
    for(index=0; index<size2; index++)
    {
        amp = data.image[IDa1].array.F[index];
        pha = data.image[IDp1].array.F[index];
        re = amp*cos(pha);
        im = amp*sin(pha);
        data.image[IDc].array.CF[index].re = re;
        data.image[IDc].array.CF[index].im = im;
    }

    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    mk_amph_from_complex("fc2", "fa2", "fp2", 0);
    arith_image_mult("fa2","fa2",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fa2");
    delete_image_ID("fp2");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2;


    return(0);
}


int coronagraph_simul_AIC(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1,IDc;
    long ID;
    double trad_pix;
    long ii,jj;
    double a11,a12,p11,p12,r11,r12,i11,i12;
    double total;
    long size2=size*size;
    long index1,index2,index1a,index2a;
    double epsilon = 1e-10;
    long iistart,iiend;

    /*  printf("Pixel scale = %f l/d per pixel\n",CORONAGRAPHS_PIXSCALE);*/

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;

    IDa1 = make_disk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    IDc = create_2DCimage_ID("pc2",size,size);

    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    iistart = (long) (0.5*size-trad_pix-2.0);
    iiend = (long) (0.5*size+trad_pix+2.0);
    if(iistart<1)
        iistart = 1;
    if(iiend>size-1)
        iiend = size-1;

    for(jj=iistart; jj<iiend; jj++)
    {
        index1a = jj*size;
        index2a = (size-jj)*size+size;
        for(ii=iistart; ii<iiend; ii++)
        {
            index1 = index1a+ii;
            index2 = index2a-ii;
            a11 = data.image[IDa1].array.F[index1];
            if(a11>epsilon)
            {
                p11 = data.image[IDp1].array.F[index1];
                a12 = data.image[IDa1].array.F[index2];
                p12 = data.image[IDp1].array.F[index2];
                r11 = a11*cos(p11);
                r12 = a12*cos(p12);
                i11 = a11*sin(p11);
                i12 = a12*sin(p12);
                data.image[IDc].array.CF[index1].re = 0.5*(r11-r12);
                data.image[IDc].array.CF[index1].im = 0.5*(i11-i12);
            }
        }
    }
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    mk_amph_from_complex("fc2", "fa2", "fp2", 0);
    arith_image_mult("fa2","fa2",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fa2");
    delete_image_ID("fp2");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2;

    return(0);
}




int coronagraph_simul_4QPM(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    long IDrefr,IDrefi;
    double x,y;
    int init4q = 0;
    char fname1[200];

    x = xld;
    y = yld;

    sprintf(fname1,"%s/4qre.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    IDrefr = image_ID("ref4qre");
    if(IDrefr==-1)
    {
        if(file_exists(fname1)==1)
            IDrefr = load_fits(fname1,"ref4qre", 1);
        else
            init4q = 1;
    }
    sprintf(fname1,"%s/4qim.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    IDrefi = image_ID("ref4qim");
    if(IDrefi==-1)
    {
        if(file_exists(fname1)==1)
            IDrefi = load_fits(fname1,"ref4qim", 1);
        else
            init4q = 1;
    }

    if(init4q==1)
    {
        printf("initialize 4 quadrant coronagraph\n");
        x = 0.0;
        y = 0.0;
    }

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*x/trad_pix,PI*y/trad_pix);

    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    mk_complex_from_amph("pa1", "pp1", "pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    ID = image_ID("fc2");
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            if((ii>size/2)&&(jj>size/2))
            {
                data.image[ID].array.CF[jj*size+ii].re *= -1.0;
                data.image[ID].array.CF[jj*size+ii].im *= -1.0;
            }
            if((ii<size/2)&&(jj<size/2))
            {
                data.image[ID].array.CF[jj*size+ii].re *= -1.0;
                data.image[ID].array.CF[jj*size+ii].im *= -1.0;
            }
        }
    for(ii=0; ii<size; ii++)
    {
        data.image[ID].array.CF[size/2*size+ii].re = 0.0;
        data.image[ID].array.CF[size/2*size+ii].im = 0.0;
        data.image[ID].array.CF[ii*size+size/2].re = 0.0;
        data.image[ID].array.CF[ii*size+size/2].im = 0.0;
    }
    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        if(data.image[IDa1].array.F[ii]<0.99999)
        {
            data.image[ID].array.CF[ii].re = 0.0;
            data.image[ID].array.CF[ii].im = 0.0;
        }

    if(init4q==1)
    {
        IDrefr = create_2Dimage_ID("ref4qre",size,size);
        IDrefi = create_2Dimage_ID("ref4qim",size,size);
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDrefr].array.F[ii] = data.image[ID].array.CF[ii].re;
            data.image[IDrefi].array.F[ii] = data.image[ID].array.CF[ii].im;
        }
        sprintf(fname1,"!%s/4qre.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        save_fl_fits("ref4qre",fname1);
        sprintf(fname1,"!%s/4qim.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        save_fl_fits("ref4qim",fname1);
        printf("4QPM initialized, please rerun program.\n");
        exit(0);
    }
    else
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[ID].array.CF[ii].re -= data.image[IDrefr].array.F[ii];
            data.image[ID].array.CF[ii].im -= data.image[IDrefi].array.F[ii];
        }
    }
    permut("pc3");
    do2dfft("pc3","fc4");
    permut("fc4");
    mk_amph_from_complex("fc4", "fa4", "fp4", 0);

    arith_image_mult("fa4","fa4",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    delete_image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    return(0);
}

int coronagraph_simul_ODC(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    long IDrefm;
    double r,tmp1,tmp2;
    char refname[200];
    char refname1[200];
    double x,y;
    int initodc=0;
    char fname[200];
    char fname1[200];
    char fnamer[200];
    char fnamei[200];
    char fname1r[200];
    char fname1i[200];
    long IDrefr,IDrefi;

    x = xld;
    y = yld;

    sprintf(fname1,"%s/odcre.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    sprintf(fname,"refocdre");
    IDrefr = image_ID(fname);

    if(IDrefr==-1)
    {
        if(file_exists(fname1)==1)
            IDrefr = load_fits(fname1,fname, 1);
        else
            initodc = 1;
    }

    sprintf(fname1,"%s/odcim.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    sprintf(fname,"refodcim");
    IDrefi = image_ID(fname);

    if(IDrefi==-1)
    {
        if(file_exists(fname1)==1)
            IDrefi = load_fits(fname1,fname, 1);
        else
            initodc = 1;
    }

    if(initodc==1)
    {
        printf("initialize ODC coronagraph\n");
        x = 0.0;
        y = 0.0;
    }


    sprintf(refname,"%s/odcm.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    sprintf(refname1,"odcmref");
    IDrefm = image_ID(refname1);

    if(IDrefm==-1)
    {
        if(file_exists(refname)==1)
            IDrefm = load_fits(refname,refname1, 1);
        else
            coronagraph_init_ODC();
    }
    IDrefm = image_ID(refname1);

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    mk_complex_from_amph("pa1", "pp1", "pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    ID = image_ID("fc2");
    for(ii=0; ii<size2; ii++)
    {
        data.image[ID].array.CF[ii].re *= data.image[IDrefm].array.F[ii];
        data.image[ID].array.CF[ii].im *= data.image[IDrefm].array.F[ii];
    }

    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
            if(r>ODCMASK_eps)
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }


    if(initodc==1)
    {
        sprintf(fnamer,"refodcre");
        sprintf(fnamei,"refodcim");
        sprintf(fname1r,"!%s/odcre.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        sprintf(fname1i,"!%s/odcim.ref.%ld", CORONAGRAPHSDATALOCAL, size);


        IDrefr = create_2Dimage_ID(fnamer,size,size);
        IDrefi = create_2Dimage_ID(fnamei,size,size);
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDrefr].array.F[ii] = data.image[ID].array.CF[ii].re;
            data.image[IDrefi].array.F[ii] = data.image[ID].array.CF[ii].im;
        }
        save_fl_fits(fnamer,fname1r);
        save_fl_fits(fnamei,fname1i);
        printf("ODC initialized, please rerun program.\n");
        exit(0);
    }
    else
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[ID].array.CF[ii].re -= data.image[IDrefr].array.F[ii];
            data.image[ID].array.CF[ii].im -= data.image[IDrefi].array.F[ii];
        }
    }


    /*  mk_amph_from_complex("pc3", "fa", "fp", 0);
    save_fl_fits("fa","!fa");
    exit(0);*/


    permut("pc3");
    do2dfft("pc3", "fc4");
    permut("fc4");
    mk_amph_from_complex("fc4", "fa4", "fp4", 0);

    arith_image_mult("fa4", "fa4", psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    delete_image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    return(0);
}


int coronagraph_simul_BL8(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    long IDrefm;
    double r,tmp1,tmp2;
    double iioffset;
    double r1,r2,tmp1b;
    char refname[200];
    char refname1[200];
    double x,y;
    int initbl8=0;
    char fname[200];
    char fname1[200];
    char fnamer[200];
    char fnamei[200];
    char fname1r[200];
    char fname1i[200];
    long IDrefr,IDrefi;

    x = xld;
    y = yld;

    if(BL8MODE==0)
    {
        sprintf(fname1,"%s/bl8re.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        sprintf(fname,"refbl8re");
    }
    else
    {
        sprintf(fname1,"%s/bl8lre.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        sprintf(fname,"refbl8lre");
    }
    IDrefr = image_ID(fname);
    if(IDrefr==-1)
    {
        if(file_exists(fname1)==1)
            IDrefr = load_fits(fname1,fname, 1);
        else
            initbl8 = 1;
    }

    if(BL8MODE==0)
    {
        sprintf(fname1,"%s/bl8im.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        sprintf(fname,"refbl8im");
    }
    else
    {
        sprintf(fname1,"%s/bl8lim.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        sprintf(fname,"refbl8lim");
    }
    IDrefi = image_ID(fname);
    if(IDrefi==-1)
    {
        if(file_exists(fname1)==1)
            IDrefi = load_fits(fname1,fname, 1);
        else
            initbl8 = 1;
    }


    if(initbl8==1)
    {
        printf("initialize BL8 coronagraph\n");
        x = 0.0;
        y = 0.0;
    }

    if(BL8MODE==0)
    {
        sprintf(refname,"%s/bl8m.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        sprintf(refname1,"bl8mref");
    }
    else
    {
        sprintf(refname,"%s/bl8ml.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        sprintf(refname1,"bl8mrefl");
    }

    IDrefm = image_ID(refname1);

    if(IDrefm==-1)
    {
        if(file_exists(refname)==1)
            IDrefm = load_fits(refname,refname1, 1);
        else
            coronagraph_init_BL8();
    }
    IDrefm = image_ID(refname1);

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    mk_complex_from_amph("pa1", "pp1", "pc2", 0);
    permut("pc2");
    do2dfft("pc2", "fc2");
    permut("fc2");
    ID = image_ID("fc2");
    for(ii=0; ii<size2; ii++)
    {
        data.image[ID].array.CF[ii].re *= data.image[IDrefm].array.F[ii];
        data.image[ID].array.CF[ii].im *= data.image[IDrefm].array.F[ii];
    }

    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");
    iioffset = (long) (BL8MASK_eps*trad_pix);
    if(BL8MODE==0)
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                tmp1 = 1.0*ii-size/2;
                tmp2 = 1.0*jj-size/2;
                r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
                if(r>1.0-BL8MASK_eps)
                {
                    data.image[ID].array.CF[jj*size+ii].re = 0.0;
                    data.image[ID].array.CF[jj*size+ii].im = 0.0;
                }
            }
    }
    else
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                tmp1 = (1.0*ii-size/2)-iioffset;
                tmp1b = (1.0*ii-size/2)+iioffset;
                tmp2 = 1.0*jj-size/2;
                r1 = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
                r2 = sqrt(tmp1b*tmp1b+tmp2*tmp2)/trad_pix;
                if((r1>1.0)||(r2>1.0))
                {
                    data.image[ID].array.CF[jj*size+ii].re = 0.0;
                    data.image[ID].array.CF[jj*size+ii].im = 0.0;
                }
            }
    }

    if(initbl8==1)
    {
        if(BL8MODE==0)
        {
            sprintf(fnamer,"refbl8re");
            sprintf(fnamei,"refbl8im");
            sprintf(fname1r,"!%s/bl8re.ref.%ld", CORONAGRAPHSDATALOCAL, size);
            sprintf(fname1i,"!%s/bl8im.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        }
        else
        {
            sprintf(fnamer,"refbl8lre");
            sprintf(fnamei,"refbl8lim");
            sprintf(fname1r,"!%s/bl8lre.ref.%ld", CORONAGRAPHSDATALOCAL, size);
            sprintf(fname1i,"!%s/bl8lim.ref.%ld", CORONAGRAPHSDATALOCAL, size);
        }

        IDrefr = create_2Dimage_ID(fnamer,size,size);
        IDrefi = create_2Dimage_ID(fnamei,size,size);
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDrefr].array.F[ii] = data.image[ID].array.CF[ii].re;
            data.image[IDrefi].array.F[ii] = data.image[ID].array.CF[ii].im;
        }
        save_fl_fits(fnamer,fname1r);
        save_fl_fits(fnamei,fname1i);
        printf("BL8 initialized, please rerun program.\n");
        exit(0);
    }
    else
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[ID].array.CF[ii].re -= data.image[IDrefr].array.F[ii];
            data.image[ID].array.CF[ii].im -= data.image[IDrefi].array.F[ii];
        }
    }


    /*  mk_amph_from_complex("pc3", "fa", "fp", 0);
    save_fl_fits("fa", "!fa");
    exit(0);*/


    permut("pc3");
    do2dfft("pc3","fc4");
    permut("fc4");
    mk_amph_from_complex("fc4", "fa4", "fp4", 0);

    arith_image_mult("fa4", "fa4", psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    delete_image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    return(0);
}


int coronagraph_simul_BL4(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    long IDrefm;
    double tmp1,tmp2;
    double iioffset;
    double r1,r2,tmp1b;
    char refname[200];
    char refname1[200];
    double x,y;
    int initbl4=0;
    char fname[200];
    char fname1[200];
    char fnamer[200];
    char fnamei[200];
    char fname1r[200];
    char fname1i[200];
    long IDrefr,IDrefi;

    x = xld;
    y = yld;


    sprintf(fname1,"%s/bl4re.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    sprintf(fname,"refbl4re");

    IDrefr = image_ID(fname);
    if(IDrefr==-1)
    {
        if(file_exists(fname1)==1)
            IDrefr = load_fits(fname1,fname, 1);
        else
            initbl4 = 1;
    }

    sprintf(fname1,"%s/bl4im.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    sprintf(fname,"refbl4im");

    IDrefi = image_ID(fname);
    if(IDrefi==-1)
    {
        if(file_exists(fname1)==1)
            IDrefi = load_fits(fname1,fname, 1);
        else
            initbl4 = 1;
    }

    if(initbl4==1)
    {
        printf("initialize BL4 coronagraph\n");
        x = 0.0;
        y = 0.0;
    }

    sprintf(refname,"%s/bl4m.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    sprintf(refname1,"bl4mref");

    IDrefm = image_ID(refname1);

    if(IDrefm==-1)
    {
        if(file_exists(refname)==1)
            IDrefm = load_fits(refname,refname1, 1);
        else
            coronagraph_init_BL4();
    }
    IDrefm = image_ID(refname1);

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    mk_complex_from_amph("pa1", "pp1", "pc2", 0);
    permut("pc2");
    do2dfft("pc2", "fc2");
    permut("fc2");
    ID = image_ID("fc2");
    for(ii=0; ii<size2; ii++)
    {
        data.image[ID].array.CF[ii].re *= data.image[IDrefm].array.F[ii];
        data.image[ID].array.CF[ii].im *= data.image[IDrefm].array.F[ii];
    }

    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");
    iioffset = (long) (BL4MASK_eps*trad_pix+0.2);
    iioffset += 1;

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = (1.0*ii-size/2)-iioffset;
            tmp1b = (1.0*ii-size/2)+iioffset;
            tmp2 = 1.0*jj-size/2;
            r1 = (sqrt(tmp1*tmp1+tmp2*tmp2)+1.0)/trad_pix;
            r2 = (sqrt(tmp1b*tmp1b+tmp2*tmp2)+1.0)/trad_pix;
            if((r1>1.0)||(r2>1.0))
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }


    if(initbl4==1)
    {
        sprintf(fnamer,"refbl4re");
        sprintf(fnamei,"refbl4im");
        sprintf(fname1r,"!%s/bl4re.ref.%ld",CORONAGRAPHSDATALOCAL,size);
        sprintf(fname1i,"!%s/bl4im.ref.%ld",CORONAGRAPHSDATALOCAL,size);

        IDrefr = create_2Dimage_ID(fnamer,size,size);
        IDrefi = create_2Dimage_ID(fnamei,size,size);
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDrefr].array.F[ii] = data.image[ID].array.CF[ii].re;
            data.image[IDrefi].array.F[ii] = data.image[ID].array.CF[ii].im;
        }
        save_fl_fits(fnamer,fname1r);
        save_fl_fits(fnamei,fname1i);
        printf("BL4 initialized, please rerun program.\n");
        exit(0);
    }
    else
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[ID].array.CF[ii].re -= data.image[IDrefr].array.F[ii];
            data.image[ID].array.CF[ii].im -= data.image[IDrefi].array.F[ii];
        }
    }


    /*  mk_amph_from_complex("pc3", "fa", "fp", 0);
    save_fl_fits("fa", "!fa");
    exit(0);*/


    permut("pc3");
    do2dfft("pc3","fc4");
    permut("fc4");
    mk_amph_from_complex("fc4", "fa4", "fp4", 0);

    arith_image_mult("fa4", "fa4", psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    delete_image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    return(0);
}


int coronagraph_simul_RRPM(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1,IDa2,IDp2;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    long IDref_fm,IDref_pm;
    double r,tmp1,tmp2;
    char refname[200];
    char refname1[200];
    long IDrefr,IDrefi;
    double x,y;
    int initrrpm = 0;
    char fname1[200];
    double v1,v2,factor;

    sprintf(PIAAAPO_FNAME,"rrpmref_pm.ref.%ld.prof",size);
    sprintf(PIAAAPODIZE_2DAPOFNAME,"rrpmref_pm.ref.%ld",size);

    sprintf(refname,"%s/rrpmref_fm.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    sprintf(refname1,"rrpmreffm");
    IDref_fm = image_ID(refname1);
    if(IDref_fm==-1)
    {
        if(file_exists(refname)==1)
            IDref_fm = load_fits(refname,refname1, 1);
        else
            coronagraph_init_RRPM();
    }
    IDref_fm = image_ID(refname1);

    sprintf(refname,"%s/rrpmref_pm.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    sprintf(refname1,"rrpmrefpm");
    IDref_pm = image_ID(refname1);
    if(IDref_pm==-1)
    {
        if(file_exists(refname)==1)
            IDref_pm = load_fits(refname,refname1, 1);
        else
            coronagraph_init_RRPM();
    }
    IDref_pm = image_ID(refname1);


    x = xld;
    y = yld;

    sprintf(fname1,"%s/rrpmre.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    IDrefr = image_ID("refrrpmre");
    if(IDrefr==-1)
    {
        if(file_exists(fname1)==1)
            IDrefr = load_fits(fname1,"refrrpmre", 1);
        else
            initrrpm = 1;
    }
    sprintf(fname1,"%s/rrpmim.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    IDrefi = image_ID("refrrpmim");
    if(IDrefi==-1)
    {
        if(file_exists(fname1)==1)
            IDrefi = load_fits(fname1,"refrrpmim", 1);
        else
            initrrpm = 1;
    }

    if(initrrpm==1)
    {
        printf("initialize RRPM coronagraph\n");
        x = 0.0;
        y = 0.0;
    }

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*x/trad_pix,PI*y/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }


    IDa2 = create_2Dimage_ID("pa2",size,size);
    for(ii=0; ii<size2; ii++)
        data.image[IDa2].array.F[ii] = data.image[IDa1].array.F[ii] * data.image[IDref_pm].array.F[ii];

    if(RRPM_PIAA==0)
    {
        IDa2 = create_2Dimage_ID("pa2",size,size);
        for(ii=0; ii<size2; ii++)
            data.image[IDa2].array.F[ii] = data.image[IDa1].array.F[ii] * data.image[IDref_pm].array.F[ii];
        factor = 1.0;
    }
    else
    {
        PIAAAPO_NBPOINTS = 80;
        if(initPIAA==0)
        {
            coronagraph_init_PIAA();
            initPIAA = 1;
        }
        coronagraphs_PIAA_apodize_beam("pa1","pp1","pa2","pp2");
        IDa2 = image_ID("pa2");
        IDp2 = image_ID("pp2");
        //      save_fl_fits("pa2","!pa2");
        //      save_fl_fits("pp2","!pp2");
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] = data.image[IDp2].array.F[ii];
        tmp1 = 0.0;
        tmp2 = 0.0;
        for(ii=0; ii<size2; ii++)
        {
            v1 = data.image[IDa1].array.F[ii];
            v2 = data.image[IDref_pm].array.F[ii];
            tmp1 += v1*v1;
            if(v1>0.5)
                tmp2 += v2*v2;
        }
        factor = tmp1/tmp2;
        // printf("factor = %f\n",factor);

        for(ii=0; ii<size2; ii++)
            data.image[IDa2].array.F[ii] = data.image[IDa1].array.F[ii] * data.image[IDref_pm].array.F[ii];
        delete_image_ID("pp2");
    }

    mk_complex_from_amph("pa2", "pp1", "pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    ID = image_ID("fc2");
    for(ii=0; ii<size2; ii++)
    {
        data.image[ID].array.CF[ii].re *= data.image[IDref_fm].array.F[ii];
        data.image[ID].array.CF[ii].im *= data.image[IDref_fm].array.F[ii];
    }

    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
            if(data.image[IDa1].array.F[jj*size+ii]<0.999)
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }

    if(initrrpm==1)
    {
        IDrefr = create_2Dimage_ID("refrrpmre",size,size);
        IDrefi = create_2Dimage_ID("refrrpmim",size,size);
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDrefr].array.F[ii] = data.image[ID].array.CF[ii].re;
            data.image[IDrefi].array.F[ii] = data.image[ID].array.CF[ii].im;
        }
        sprintf(fname1,"!%s/rrpmre.ref.%ld",CORONAGRAPHSDATALOCAL,size);
        save_fl_fits("refrrpmre",fname1);
        sprintf(fname1,"!%s/rrpmim.ref.%ld",CORONAGRAPHSDATALOCAL,size);
        save_fl_fits("refrrpmim",fname1);
        printf("RRPM initialized, please rerun program.\n");
        exit(0);
    }
    else
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[ID].array.CF[ii].re -= data.image[IDrefr].array.F[ii];
            data.image[ID].array.CF[ii].im -= data.image[IDrefi].array.F[ii];
        }
    }

    /*  mk_amph_from_complex("pc3", "fa", "fp", 0);
    save_fl_fits("fa", "!fa");
    exit(0);*/


    permut("pc3");
    do2dfft("pc3","fc4");
    permut("fc4");
    mk_amph_from_complex("fc4", "fa4", "fp4", 0);

    arith_image_mult("fa4","fa4",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pa2");
    delete_image_ID("pp1");
    delete_image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2/factor;

    return(0);
}

/* optical vortex coronagraph */
int coronagraph_simul_OVC(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    long IDref_fmp,IDref_fma;
    double r,tmp1,tmp2;
    char refname[200];
    char refname1[200];
    long IDrefr,IDrefi;
    double x,y;
    int initovc = 0;
    char fname1[200];
    double re,im,amp,pha;

    sprintf(refname,"%s/ovcref_fmp%ld.ref.%ld",CORONAGRAPHSDATALOCAL,OVC_CHARGE,size);
    sprintf(refname1,"ovcreffmp");
    IDref_fmp = image_ID(refname1);
    if(IDref_fmp==-1)
    {
        if(file_exists(refname)==1)
            IDref_fmp = load_fits(refname,refname1, 1);
        else
            coronagraph_init_OVC(OVC_CHARGE);
    }
    IDref_fmp = image_ID(refname1);
    //  save_fl_fits("ovcreffmp","!ovcreffmp");

    sprintf(refname,"%s/ovcref_fma%ld.ref.%ld",CORONAGRAPHSDATALOCAL,OVC_CHARGE,size);
    sprintf(refname1,"ovcreffma");
    IDref_fma = image_ID(refname1);
    if(IDref_fma==-1)
    {
        if(file_exists(refname)==1)
            IDref_fma = load_fits(refname,refname1, 1);
        else
            coronagraph_init_OVC(OVC_CHARGE);
    }
    IDref_fma = image_ID(refname1);


    x = xld;
    y = yld;

    sprintf(fname1,"%s/ovcre%ld.ref.%ld",CORONAGRAPHSDATALOCAL,OVC_CHARGE,size);
    IDrefr = image_ID("refovcre");
    if(IDrefr==-1)
    {
        if(file_exists(fname1)==1)
            IDrefr = load_fits(fname1,"refovcre", 1);
        else
            initovc = 1;
    }
    sprintf(fname1,"%s/ovcim%ld.ref.%ld",CORONAGRAPHSDATALOCAL,OVC_CHARGE,size);
    IDrefi = image_ID("refovcim");
    if(IDrefi==-1)
    {
        if(file_exists(fname1)==1)
            IDrefi = load_fits(fname1,"refovcim", 1);
        else
            initovc = 1;
    }

    if(initovc==1)
    {
        printf("initialize OVC coronagraph charge %ld\n",OVC_CHARGE);
        x = 0.0;
        y = 0.0;
    }

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    /*  for(ii=0;ii<size;ii++)
      for(jj=0;jj<size;jj++)
        {
    x1 = 1.0*ii-size/2;
    y1 = 1.0*jj-size/2;
    if(sqrt(x1*x1+y1*y1)<trad_pix*0.2)
        data.image[IDa1].array.F[jj*size+ii] = 0.0;
        }
    */

    total = arith_image_total("pa1");
    ID = image_ID("coramp");
    if(ID!=-1)
    {
        printf("---------- USING AMPLITUDE FILE ----------\n");
        for(ii=0; ii<size2; ii++)
            data.image[IDa1].array.F[ii] = data.image[ID].array.F[ii];
    }


    IDp1 = make_slopexy("pp1",size,size,PI*x/trad_pix,PI*y/trad_pix);
    ID = image_ID("corpha");
    if(ID!=-1)
    {
        printf("---------- USING PHASE FILE ----------\n");
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    mk_complex_from_amph("pa1","pp1","pc2", 0);
    if(0==0)
    {
        save_fl_fits("pa1","!pa1");
        save_fl_fits("pp1","!pp1");
    }


    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    ID = image_ID("fc2");


    if(0==0)
    {
        mk_amph_from_complex("fc2", "fa2", "fp2", 0);
        save_fl_fits("fa2","!fa2");
        save_fl_fits("fp2","!fp2");
    }

    for(ii=0; ii<size2; ii++)
    {
        re = data.image[ID].array.CF[ii].re;
        im = data.image[ID].array.CF[ii].im;
        amp = sqrt(re*re+im*im);
        pha = atan2(im,re);
        amp *= data.image[IDref_fma].array.F[ii];
        pha += data.image[IDref_fmp].array.F[ii];
        data.image[ID].array.CF[ii].re = amp*cos(pha);
        data.image[ID].array.CF[ii].im = amp*sin(pha);
    }

    /*  mk_amph_from_complex("fc2","fa2","fp2");
    save_fl_fits("fa2","!fa2");
    save_fl_fits("fp2","!fp2");
    delete_image_ID("fa2");
    delete_image_ID("fp2");*/


    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");

    if(0==0)
    {
        mk_amph_from_complex("pc3", "pc3a", "pc3p", 0);
        save_fl_fits("pc3a", "!pc3a");
        save_fl_fits("pc3p", "!pc3p");
        delete_image_ID("pc3a");
        delete_image_ID("pc3p");
    }

    /*  for(ii=0;ii<size;ii++)
      for(jj=0;jj<size;jj++)
        {
    tmp1 = 1.0*ii-size/2;
    tmp2 = 1.0*jj-size/2;
    r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
    if(data.image[IDa1].array.F[jj*size+ii]<0.999)
      {
        data.image[ID].array.CF[jj*size+ii].re = 0.0;
        data.image[ID].array.CF[jj*size+ii].im = 0.0;
      }
        }
    */



    if(initovc==1)
    {
        IDrefr = create_2Dimage_ID("refovcre",size,size);
        IDrefi = create_2Dimage_ID("refovcim",size,size);
        for(ii=0; ii<size2; ii++)
        {
            data.image[IDrefr].array.F[ii] = data.image[ID].array.CF[ii].re;
            data.image[IDrefi].array.F[ii] = data.image[ID].array.CF[ii].im;
        }
        sprintf(fname1,"!%s/ovcre%ld.ref.%ld",CORONAGRAPHSDATALOCAL,OVC_CHARGE,size);
        save_fl_fits("refovcre",fname1);
        sprintf(fname1,"!%s/ovcim%ld.ref.%ld",CORONAGRAPHSDATALOCAL,OVC_CHARGE,size);
        save_fl_fits("refovcim",fname1);
        printf("OVC charge %ld initialized, please rerun program.\n",OVC_CHARGE);
        exit(0);
    }
    else
    {
        for(ii=0; ii<size2; ii++)
        {
            data.image[ID].array.CF[ii].re -= data.image[IDrefr].array.F[ii];
            data.image[ID].array.CF[ii].im -= data.image[IDrefi].array.F[ii];
        }
    }

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
            //	if(data.image[IDa1].array.F[jj*size+ii]<0.99999)
            if(r>0.96)
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }




    /*  mk_amph_from_complex("pc3", "fa", "fp", 0);
    save_fl_fits("fa", "!fa");
    exit(0);*/


    permut("pc3");
    do2dfft("pc3", "fc4");
    permut("fc4");
    mk_amph_from_complex("fc4", "fa4", "fp4", 0);

    arith_image_mult("fa4", "fa4", psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    delete_image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    return(0);
}


int coronagraph_simul_CPA(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1,IDa2;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    long IDref;
    double r,tmp1,tmp2;
    char refname[200];
    char refname1[200];


    sprintf(refname,"%s/cpapup.ref.%ld",CORONAGRAPHSDATALOCAL,size);
    sprintf(refname1,"cpapupref");
    IDref = image_ID(refname1);
    if(IDref==-1)
    {
        if(file_exists(refname)==1)
            IDref = load_fits(refname,refname1, 1);
        else
            coronagraph_init_CPA();
    }
    IDref = image_ID(refname1);

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }


    IDa2 = create_2Dimage_ID("pa2",size,size);
    for(ii=0; ii<size2; ii++)
        data.image[IDa2].array.F[ii] = data.image[IDa1].array.F[ii] * data.image[IDref].array.F[ii];

    /*save_fl_fits("pa2","!pa2");*/
    mk_complex_from_amph("pa2", "pp1", "pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");

    ID = image_ID("fc2");
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)*CORONAGRAPHS_PIXSCALE;
            if(r<CPAFPMASKRAD)
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }


    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
            if(data.image[IDa1].array.F[jj*size+ii]>CPAPPMASKRAD)
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }


    /*  mk_amph_from_complex("pc3","fa","fp", 0);
    save_fl_fits("fa","!fa");
    exit(0);*/


    permut("pc3");
    do2dfft("pc3","fc4");
    permut("fc4");
    mk_amph_from_complex("fc4", "fa4", "fp4", 0);

    arith_image_mult("fa4", "fa4", psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pa2");
    delete_image_ID("pp1");
    delete_image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    return(0);
}


int coronagraph_simul_PPA(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1,IDa2;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    double x;
    double a = 2.0;
    double epsilon = 0.001;
    double re,im;
    long IDrefa,IDrefp;

    double *pha_array;
    long *index_array;
    long NBsample = 10000;
    long IDcnt;
    double phax,phay,pha;
    long kx,ky,k;
    long IDrefr,IDrefi;
    long index;

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    IDrefa = image_ID("pparefa");
    IDrefp = image_ID("pparefp");

    if((IDrefa==-1)||(IDrefp==-1))
    {
        pha_array = (double*) malloc(NBsample*sizeof(double));
        index_array = (long*) malloc(NBsample*sizeof(long));
        for(k=0; k<NBsample; k++)
        {
            x = 2.0*k/NBsample-1.0;
            pha_array[k] = a*log(((1.0+epsilon) + x)/((1.0+epsilon) - x));
            index_array[k] = (long) (x*(trad_pix-1.0)/sqrt(2.0)+size/2);
        }


        IDrefr = create_2Dimage_ID("pparefr",size,size);
        IDrefi = create_2Dimage_ID("pparefi",size,size);
        IDrefa = create_2Dimage_ID("pparefa",size,size);
        IDrefp = create_2Dimage_ID("pparefp",size,size);
        IDcnt = create_2Dimage_ID("cntref",size,size);

        for(kx=0; kx<NBsample; kx++)
        {
            phax = pha_array[kx];
            ii = index_array[kx];
            for(ky=0; ky<NBsample; ky++)
            {
                phay = pha_array[ky];
                jj = index_array[ky];
                pha = phax+phay;
                index = jj*size+ii;
                data.image[IDrefr].array.F[index] += cos(pha);
                data.image[IDrefi].array.F[index] += sin(pha);
                data.image[IDcnt].array.F[index] += 1.0;
            }
        }

        for(ii=0; ii<size2; ii++)
            if(data.image[IDcnt].array.F[ii]>0.5)
            {
                re = data.image[IDrefr].array.F[ii]/data.image[IDcnt].array.F[ii];
                im = data.image[IDrefi].array.F[ii]/data.image[IDcnt].array.F[ii];
                data.image[IDrefa].array.F[ii] = sqrt(re*re+im*im);
                data.image[IDrefp].array.F[ii] = atan2(im,re);
            }
        delete_image_ID("pparefr");
        delete_image_ID("pparefi");
        delete_image_ID("cntref");
        free(pha_array);
        free(index_array);
    }



    IDa2 = create_2Dimage_ID("pa2",size,size);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            data.image[IDp1].array.F[jj*size+ii] += data.image[IDrefp].array.F[jj*size+ii];
            data.image[IDa2].array.F[jj*size+ii] = data.image[IDa1].array.F[jj*size+ii]*data.image[IDrefa].array.F[jj*size+ii];
        }


    /*save_fl_fits("pa2","!pa2");
      save_fl_fits("pp1","!pp1");*/
    mk_complex_from_amph("pa2","pp1","pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");


    mk_amph_from_complex("fc2","fa2","fp2", 0);

    arith_image_mult("fa2","fa2",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fa2");
    delete_image_ID("fp2");
    delete_image_ID("pa1");
    delete_image_ID("pa2");
    delete_image_ID("pp1");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2;

    return(0);
}


int coronagraph_simul_NOCORO(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID;
    double trad_pix;
    long ii;
    double total;
    long size2=size*size;

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }


    /*save_fl_fits("pa2","!pa2");
      save_fl_fits("pp1","!pp1");*/
    mk_complex_from_amph("pa1","pp1","pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");


    mk_amph_from_complex("fc2","fa2","fp2", 0);

    arith_image_mult("fa2","fa2",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fa2");
    delete_image_ID("fp2");
    delete_image_ID("pa1");
    delete_image_ID("pa2");
    delete_image_ID("pp1");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2;

    return(0);
}


int coronagraph_simul_PIAA(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID,IDm,IDp;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    double r,tmp1,tmp2;
    double re,im,amp,pha,phap,pham;


    if(initPIAA==0)
    {
        coronagraph_init_PIAA();
        initPIAA = 1;
    }


    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix*PIAAOVERSIZE);
    total = arith_image_total("pa1")/PIAAOVERSIZE/PIAAOVERSIZE;
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }


    coronagraphs_PIAA_apodize_beam("pa1","pp1","pa2","pp2");

    mk_complex_from_amph("pa2","pp2","pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");

    ID = image_ID("fc2");


    if(PIAAFPMASK==1)
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                tmp1 = 1.0*ii-size/2;
                tmp2 = 1.0*jj-size/2;
                r = sqrt(tmp1*tmp1+tmp2*tmp2)*CORONAGRAPHS_PIXSCALE;
                if(r<PIAAFPMASKRAD)
                {
                    data.image[ID].array.CF[jj*size+ii].re = 0.0;
                    data.image[ID].array.CF[jj*size+ii].im = 0.0;
                }
            }
    }
    if(PIAALOWFS==1)
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                tmp1 = 1.0*ii-size/2;
                tmp2 = 1.0*jj-size/2;
                r = sqrt(tmp1*tmp1+tmp2*tmp2)*CORONAGRAPHS_PIXSCALE;
                if((r>PIAAFPMASKRAD*1.2)||(r<1.2*PIAAFPMASKRAD/4))
                {
                    data.image[ID].array.CF[jj*size+ii].re = 0.0;
                    data.image[ID].array.CF[jj*size+ii].im = 0.0;
                }
            }
    }


    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
            if(r>PIAAPPMASKRAD)
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }


    if(PIAALOWFS==1)
    {
        copy_image_ID("pc3", "pc3p", 0);
        copy_image_ID("pc3", "pc3m", 0);
        IDp = image_ID("pc3p");
        IDm = image_ID("pc3m");
        ID = image_ID("pc3");
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                tmp1 = 1.0*ii-size/2;
                tmp2 = 1.0*jj-size/2;
                r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
                if(r<PIAAPPMASKRAD)
                {
                    re = data.image[ID].array.CF[jj*size+ii].re;
                    im = data.image[ID].array.CF[jj*size+ii].im;
                    amp = sqrt(re*re+im*im);
                    pha = atan2(im,re);
                    phap = pha+r*r*6.0;
                    pham = pha-r*r*6.0;
                    data.image[IDp].array.CF[jj*size+ii].re = amp*cos(phap);
                    data.image[IDp].array.CF[jj*size+ii].im = amp*sin(phap);
                    data.image[IDm].array.CF[jj*size+ii].re = amp*cos(pham);
                    data.image[IDm].array.CF[jj*size+ii].im = amp*sin(pham);
                }
            }

        permut("pc3p");
        do2dfft("pc3p","fc4p");
        permut("fc4p");
        mk_amph_from_complex("fc4p","fa4p","fp4p", 0);
        arith_image_mult("fa4p","fa4p","psfWFSp");
        ID = image_ID("psfWFSp");
        for(ii=0; ii<size2; ii++)
            data.image[ID].array.F[ii] /= total*size2*size2*size2;
        delete_image_ID("pc3p");
        delete_image_ID("fc4p");
        delete_image_ID("fa4p");
        delete_image_ID("fp4p");

        permut("pc3m");
        do2dfft("pc3m","fc4m");
        permut("fc4m");
        mk_amph_from_complex("fc4m","fa4m","fp4m", 0);
        arith_image_mult("fa4m","fa4m","psfWFSm");
        ID = image_ID("psfWFSm");
        for(ii=0; ii<size2; ii++)
            data.image[ID].array.F[ii] /= total*size2*size2*size2;
        delete_image_ID("pc3m");
        delete_image_ID("fc4m");
        delete_image_ID("fa4m");
        delete_image_ID("fp4m");
    }

    /* mk_amph_from_complex("pc3","fa","fp", 0);
    save_fl_fits("fa","!fa");
    exit(0);*/


    permut("pc3");
    do2dfft("pc3","fc4");
    permut("fc4");
    mk_amph_from_complex("fc4","fa4","fp4", 0);


    arith_image_mult("fa4","fa4",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pa2");
    delete_image_ID("pp1");
    delete_image_ID("pp2");
    delete_image_ID("pc3");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    return(0);
}



// This routine works for PSFs within 20 l/D from the center
int coronagraph_simul_PIAAC(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    double r,tmp1,tmp2,tmp3;
    int CHEAT = 0; // if separation > CHEATLIM, fix up output beam phase sampling problem
    // if CHEAT=1, final complex amplitude from what hits the focal plane mask is subtracted to the original complex amplitude
    //  double CHEATLIM = 6.5;
    long ID1,ID2,ID3,ID4,IDb,IDc;
    double re1,im1,re2,im2,re,im;
    long index1,index2;

    CHEAT = 0;


    if(initPIAA==0)
    {
        coronagraph_init_PIAA();
        initPIAA = 1;
    }


    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix*PIAAOVERSIZE);
    total = arith_image_total("pa1")/PIAAOVERSIZE/PIAAOVERSIZE;

    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    /*  make_subpixdisk("pa12",size,size,0.55*size,0.52*size,0.2*trad_pix*PIAAOVERSIZE);
    execute_arith("pp1=pp1+pa12");
    IDp1 = image_ID("pp1");*/
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }


    piaaconfdirection = 0;
    coronagraphs_PIAA_apodize_beam("pa1","pp1","pa2","pp2");

    //  save_fl_fits("pp2","!pp2");
    //  exit(0);

    mk_complex_from_amph("pa2","pp2","pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    ID = image_ID("fc2");


    /*  if(sqrt(xld*xld+yld*yld)>CHEATLIM)
      CHEAT = 1;
    else
    CHEAT = 0;*/

    //  CHEAT = 0;

    if(CHEAT == 1)
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                tmp1 = 1.0*ii-size/2;
                tmp2 = 1.0*jj-size/2;
                r = sqrt(tmp1*tmp1+tmp2*tmp2)*CORONAGRAPHS_PIXSCALE;
                if(r>PIAAFPMASKRAD)
                {
                    data.image[ID].array.CF[jj*size+ii].re = 0.0;
                    data.image[ID].array.CF[jj*size+ii].im = 0.0;
                }
            }
    }
    else
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                tmp1 = 1.0*ii-size/2;
                tmp2 = 1.0*jj-size/2;
                r = sqrt(tmp1*tmp1+tmp2*tmp2)*CORONAGRAPHS_PIXSCALE;
                if(r<PIAAFPMASKRAD)
                {
                    data.image[ID].array.CF[jj*size+ii].re = 0.0;
                    data.image[ID].array.CF[jj*size+ii].im = 0.0;
                }
            }
    }

    /*  mk_amph_from_complex("fc2","tfa2","tfp2", 0);
    save_fl_fits("tfa2","!tfa2");
    delete_image_ID("tfa2");
    delete_image_ID("tfp2");
    */

    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");
    mk_amph_from_complex("pc3","pa3t","pp3t", 0);
    //  save_fl_fits("pa3t","!pa3t");
    delete_image_ID("pa3t");
    delete_image_ID("pp3t");

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
            if(r>PIAAPPMASKRAD)
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }

    mk_amph_from_complex("pc3","pa3","pp3", 0);
    /*save_fl_fits("pa3","!pa3");*/
    copy_image_ID("pp3", "pp3b", 0);
    IDb = image_ID("pp3b");
    for(ii=0; ii<size*size; ii++)
    {
        tmp1 = data.image[IDb].array.F[ii];
        tmp1 += 2.0*PI/3.0;
        if(tmp1>PI)
            tmp1 -= 2.0*PI;
        data.image[IDb].array.F[ii] = tmp1;
    }
    copy_image_ID("pp3", "pp3c", 0);
    IDc = image_ID("pp3c");
    for(ii=0; ii<size*size; ii++)
    {
        tmp1 = data.image[IDc].array.F[ii];
        tmp1 -= 2.0*PI/3.0;
        if(tmp1<-PI)
            tmp1 += 2.0*PI;
        data.image[IDc].array.F[ii] = tmp1;
    }

    piaaconfdirection = 1;
    coronagraphs_PIAA_apodize_beam("pa3","pp3","pa4","pp4");
    coronagraphs_PIAA_apodize_beam("pa3","pp3b","pa4b","pp4b");
    coronagraphs_PIAA_apodize_beam("pa3","pp3c","pa4c","pp4c");
    piaaconfdirection = 0;

    ID = image_ID("pp4");
    IDb = image_ID("pp4b");
    IDc = image_ID("pp4c");
    for(ii=0; ii<size*size; ii++)
    {
        tmp1 = data.image[ID].array.F[ii];
        tmp2 = data.image[IDb].array.F[ii];
        tmp3 = data.image[IDc].array.F[ii];

        if((fabs(tmp2-(2.0*PI/3.0-PI))<0.3*PI)&&(fabs(tmp3+(2.0*PI/3.0-PI))<0.3*PI))
        {
            tmp1 = PI+ 0.5*(tmp2+tmp3);
            if(tmp1>PI)
                tmp1 -= 2.0*PI;
            data.image[ID].array.F[ii] = tmp1;
        }
    }

    delete_image_ID("pp3b");
    delete_image_ID("pp3c");
    delete_image_ID("pa4b");
    delete_image_ID("pp4b");
    delete_image_ID("pa4c");
    delete_image_ID("pp4c");

    if(CHEAT==1)
    {
        printf("------------ MASK SUBTRACTION ----------\n");
        ID1 = image_ID("pa1");
        ID2 = image_ID("pp1");
        ID3 = image_ID("pa4");
        ID4 = image_ID("pp4");

        for(ii=10; ii<size-10; ii++)
            for(jj=10; jj<size-10; jj++)
            {
                index1 = jj*size+ii;
                index2 = (size-jj)*size+(size-ii);
                re1 = data.image[ID1].array.F[index2]*cos(data.image[ID2].array.F[index2]);
                im1 = data.image[ID1].array.F[index2]*sin(data.image[ID2].array.F[index2]);
                re2 = data.image[ID3].array.F[index1]*cos(data.image[ID4].array.F[index1]);
                im2 = data.image[ID3].array.F[index1]*sin(data.image[ID4].array.F[index1]);

                re = re1*size*size-re2;
                im = im1*size*size-im2;

                tmp1 = 1.0*ii-size/2;
                tmp2 = 1.0*jj-size/2;
                r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
                if(r>PIAAPPMASKRAD)
                {
                    re = 0.0;
                    im = 0.0;
                }

                data.image[ID3].array.F[index1] = sqrt(re*re+im*im);
                data.image[ID4].array.F[index1] = atan2(im,re);
            }
    }

    mk_complex_from_amph("pa4","pp4","pc4", 0);

    /*  save_fl_fits("pa4","!pa4");*/

    permut("pc4");
    do2dfft("pc4","fc4");
    permut("fc4");
    mk_amph_from_complex("fc4","fa4","fp4", 0);

    arith_image_mult("fa4","fa4",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pa2");
    delete_image_ID("pa3");
    delete_image_ID("pa4");
    delete_image_ID("pp1");
    delete_image_ID("pp2");
    delete_image_ID("pp3");
    delete_image_ID("pp4");
    delete_image_ID("pc3");
    delete_image_ID("pc4");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    return(0);
}


int coronagraph_simul_STRIPC(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID,IDc;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    double rad1,rad2,rad3,rad4;
    long index1,index2,index3,index4;
    double a1,p1,a2,p2,a3,p3,a4,p4;
    double r1,r2,i1,i2,r3,i3,r4,i4;
    long OFFSET;
    long totcnt;
    /*  long IDpi;*/

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_disk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");

    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }
    IDc = create_2DCimage_ID("pc2",size,size);

    OFFSET = (long) (trad_pix*STRIPCOFFSET);

    /*  IDpi = create_2Dimage_ID("pupimstripc",size,size);*/

    totcnt = 0;
    for(ii=0; ii<trad_pix; ii++)
    {
        for(jj=0; jj<trad_pix; jj++)
        {
            rad1 = sqrt(ii*ii+jj*jj);
            rad2 = sqrt((2.0*OFFSET-ii)*(2.0*OFFSET-ii)+(2.0*OFFSET-jj)*(2.0*OFFSET-jj));
            rad3 = sqrt((2.0*OFFSET-ii)*(2.0*OFFSET-ii)+jj*jj);
            rad4 = sqrt(ii*ii+(2.0*OFFSET-jj)*(2.0*OFFSET-jj));
            if((rad1<(trad_pix))&&(rad2<(trad_pix))&&(rad3<(trad_pix))&&(rad4<(trad_pix))&&(ii<2.0*OFFSET)&&(jj<2.0*OFFSET))
            {
                totcnt++;
                index1 = (size/2+jj)*size+(size/2+ii);
                index2 = (size/2-2*OFFSET+jj)*size+(size/2+ii);
                index3 = (size/2-2*OFFSET+jj)*size+(size/2-2*OFFSET+ii);
                index4 = (size/2+jj)*size+(size/2-2*OFFSET+ii);

                /*		data.image[IDpi].array.F[index1] += 1.0;
                data.image[IDpi].array.F[index2] += 1.0;
                data.image[IDpi].array.F[index3] += 1.0;
                data.image[IDpi].array.F[index4] += 1.0;*/

                a1 = data.image[IDa1].array.F[index1];
                p1 = data.image[IDp1].array.F[index1];
                a2 = data.image[IDa1].array.F[index2];
                p2 = data.image[IDp1].array.F[index2];
                a3 = data.image[IDa1].array.F[index3];
                p3 = data.image[IDp1].array.F[index3];
                a4 = data.image[IDa1].array.F[index4];
                p4 = data.image[IDp1].array.F[index4];
                r1 = a1*cos(p1);
                i1 = a1*sin(p1);
                r2 = a2*cos(p2);
                i2 = a2*sin(p2);
                r3 = a1*cos(p3);
                i3 = a1*sin(p3);
                r4 = a2*cos(p4);
                i4 = a2*sin(p4);
                data.image[IDc].array.CF[index1].re = 0.25*((r1-r4)-(r2-r3));
                data.image[IDc].array.CF[index1].im = 0.25*((i1-i4)-(i2-i3));

                data.image[IDc].array.CF[index2].re = 0.25*((r2-r3)-(r1-r4));
                data.image[IDc].array.CF[index2].im = 0.25*((i2-i3)-(i1-i4));

                data.image[IDc].array.CF[index3].re = 0.25*((r3-r2)-(r4-r1));
                data.image[IDc].array.CF[index3].im = 0.25*((i3-i2)-(i4-i1));

                data.image[IDc].array.CF[index4].re = 0.25*((r4-r1)-(r3-r2));
                data.image[IDc].array.CF[index4].im = 0.25*((i4-i1)-(i3-i2));
            }
        }
    }

    /*  save_fl_fits("pupimstripc","!pupimstripc");
        delete_image_ID("pupimstripc");*/

    /*  printf("T = %f\n",totcnt/(PI/4*trad_pix*trad_pix));*/

    mk_amph_from_complex("pc2","pa2","pp2", 0);
    /* save_fl_fits("pa2","!pa2");
       save_fl_fits("pp2","!pp2");*/
    delete_image_ID("pa2");
    delete_image_ID("pp2");

    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    mk_amph_from_complex("fc2","fa2","fp2", 0);
    arith_image_mult("fa2","fa2",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fa2");
    delete_image_ID("fp2");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2;

    return(0);
}

int coronagraph_simul_SIMXY(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID,IDc;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    double rad1;
    long index1,index2,index3,index4;
    double a1,p1,a2,p2,a3,p3,a4,p4;
    double r1,r2,i1,i2,r3,i3,r4,i4;

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_disk("pa1",size,size,0.5*size,0.5*size,trad_pix);
    total = arith_image_total("pa1");

    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }
    IDc = create_2DCimage_ID("pc2",size,size);

    for(ii=0; ii<trad_pix; ii++)
    {
        for(jj=0; jj<trad_pix; jj++)
        {
            rad1 = sqrt(ii*ii+jj*jj);
            if(rad1<trad_pix)
            {
                index1 = (size/2+jj)*size+(size/2+ii);
                index2 = (size/2-jj)*size+(size/2+ii);
                index3 = (size/2-jj)*size+(size/2-ii);
                index4 = (size/2+jj)*size+(size/2-ii);

                a1 = data.image[IDa1].array.F[index1];
                p1 = data.image[IDp1].array.F[index1];
                a2 = data.image[IDa1].array.F[index2];
                p2 = data.image[IDp1].array.F[index2];
                a3 = data.image[IDa1].array.F[index3];
                p3 = data.image[IDp1].array.F[index3];
                a4 = data.image[IDa1].array.F[index4];
                p4 = data.image[IDp1].array.F[index4];
                r1 = a1*cos(p1);
                i1 = a1*sin(p1);
                r2 = a2*cos(p2);
                i2 = a2*sin(p2);
                r3 = a1*cos(p3);
                i3 = a1*sin(p3);
                r4 = a2*cos(p4);
                i4 = a2*sin(p4);
                data.image[IDc].array.CF[index1].re = 0.25*((r1-r4)-(r2-r3));
                data.image[IDc].array.CF[index1].im = 0.25*((i1-i4)-(i2-i3));
                data.image[IDc].array.CF[index2].re = 0.25*((r2-r3)-(r1-r4));
                data.image[IDc].array.CF[index2].im = 0.25*((i2-i3)-(i1-i4));
                data.image[IDc].array.CF[index3].re = 0.25*((r3-r2)-(r4-r1));
                data.image[IDc].array.CF[index3].im = 0.25*((i3-i2)-(i4-i1));
                data.image[IDc].array.CF[index4].re = 0.25*((r4-r1)-(r3-r2));
                data.image[IDc].array.CF[index4].im = 0.25*((i4-i1)-(i3-i2));
            }
        }
    }

    mk_amph_from_complex("pc2","pa2","pp2", 0);
    /*save_fl_fits("pa2","!pa2");
      save_fl_fits("pp2","!pp2");*/
    delete_image_ID("pa2");
    delete_image_ID("pp2");

    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");
    mk_amph_from_complex("fc2","fa2","fp2", 0);
    arith_image_mult("fa2","fa2",psfname);
    ID = image_ID(psfname);
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fa2");
    delete_image_ID("fp2");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2;

    return(0);
}


int coronagraph_simul_AIC_PIAAC(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long IDa1,IDp1;
    long ID;
    double trad_pix;
    long ii,jj;
    double total;
    long size2=size*size;
    double r,tmp1,tmp2;
    long IDaicfield,IDaicc,IDaicd;
    double r11,r12,i11,i12;
    long IDpiaac;

    if(initPIAA==0)
    {
        coronagraph_init_PIAA();
        initPIAA = 1;
    }

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    IDa1 = make_subpixdisk("pa1",size,size,0.5*size,0.5*size,trad_pix*PIAAOVERSIZE);
    total = arith_image_total("pa1")/PIAAOVERSIZE/PIAAOVERSIZE;

    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);
    ID = image_ID("corphase");
    if(ID!=-1)
    {
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }


    piaaconfdirection = 0;
    coronagraphs_PIAA_apodize_beam("pa1","pp1","pa2","pp2");

    mk_complex_from_amph("pa2","pp2","pc2", 0);
    permut("pc2");
    do2dfft("pc2","fc2");
    permut("fc2");

    IDaicfield = create_2DCimage_ID("fc2_aic",size,size);
    ID = image_ID("fc2");
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)*CORONAGRAPHS_PIXSCALE;

            if(r<PIAAAIC_FIELDLIMIT)
            {
                data.image[IDaicfield].array.CF[jj*size+ii].re = data.image[ID].array.CF[jj*size+ii].re;
                data.image[IDaicfield].array.CF[jj*size+ii].im = data.image[ID].array.CF[jj*size+ii].im;

                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }

    /* Process PIAAC part */

    permut("fc2");
    do2dfft("fc2","pc3");
    permut("pc3");
    ID = image_ID("pc3");

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            tmp1 = 1.0*ii-size/2;
            tmp2 = 1.0*jj-size/2;
            r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
            if(r>PIAAPPMASKRAD)
            {
                data.image[ID].array.CF[jj*size+ii].re = 0.0;
                data.image[ID].array.CF[jj*size+ii].im = 0.0;
            }
        }

    mk_amph_from_complex("pc3","pa3","pp3", 0);
    piaaconfdirection = 1;
    coronagraphs_PIAA_apodize_beam("pa3","pp3","pa4","pp4");
    piaaconfdirection = 0;
    mk_complex_from_amph("pa4","pp4","pc4", 0);

    permut("pc4");
    do2dfft("pc4","fc4");
    permut("fc4");
    mk_amph_from_complex("fc4","fa4","fp4", 0);

    arith_image_mult("fa4","fa4","psfPIAAC");
    ID = image_ID("psfPIAAC");
    delete_image_ID("pc2");
    delete_image_ID("fc2");
    delete_image_ID("fc4");
    delete_image_ID("fa4");
    delete_image_ID("fp4");
    delete_image_ID("pa1");
    delete_image_ID("pa2");
    delete_image_ID("pa3");
    delete_image_ID("pa4");
    delete_image_ID("pp1");
    delete_image_ID("pp2");
    delete_image_ID("pp3");
    delete_image_ID("pp4");
    delete_image_ID("pc3");
    delete_image_ID("pc4");
    for(ii=0; ii<size2; ii++)
        data.image[ID].array.F[ii] /= total*size2*size2*size2;

    /* Process AIC part */
    IDaicc = create_2Dimage_ID("psfAICc",size,size);
    IDaicd = create_2Dimage_ID("psfAICd",size,size);
    for(ii=1; ii<size; ii++)
        for(jj=1; jj<size; jj++)
        {
            r11 = data.image[IDaicfield].array.CF[jj*size+ii].re;
            r12 = data.image[IDaicfield].array.CF[(size-jj)*size+(size-ii)].re;
            i11 = data.image[IDaicfield].array.CF[jj*size+ii].im;
            i12 = data.image[IDaicfield].array.CF[(size-jj)*size+(size-ii)].im;
            data.image[IDaicd].array.F[jj*size+ii] = 0.25*((r11-r12)*(r11-r12)+(i11-i12)*(i11-i12));
            data.image[IDaicc].array.F[jj*size+ii] = 0.25*((r11+r12)*(r11+r12)+(i11+i12)*(i11+i12));
            data.image[IDaicd].array.F[jj*size+ii] /= total*size2;
            data.image[IDaicc].array.F[jj*size+ii] /= total*size2;
        }

    delete_image_ID("fc2_aic");

    /* combine the 3 images */
    IDpiaac = image_ID("psfPIAAC");

    ID = create_2Dimage_ID(psfname,size,size);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size/2; jj++)
            data.image[ID].array.F[jj*size+ii] = data.image[IDpiaac].array.F[(jj+size/4)*size+ii];

    for(ii=0; ii<size/2; ii++)
        for(jj=size/2; jj<size; jj++)
            data.image[ID].array.F[jj*size+ii] = data.image[IDaicd].array.F[(jj-size/4)*size+(ii+size/4)];

    for(ii=size/2; ii<size; ii++)
        for(jj=size/2; jj<size; jj++)
            data.image[ID].array.F[jj*size+ii] = data.image[IDaicc].array.F[(jj-size/4)*size+(ii-size/4)];

    delete_image_ID("psfPIAAC");
    delete_image_ID("psfAICc");
    delete_image_ID("psfAICd");

    return(0);
}


int coronagraph_RRPM_optimize()
{
    double best_RRPM_RADIUS = RRPM_RADIUS;
    double best_RRPM_P2 = RRPM_P2;
    double best_RRPM_P3 = RRPM_P3;
    double best_RRPM_P4 = RRPM_P4;
    double best_RRPM_P5 = RRPM_P5;
    double best_value;
    double value;
    char command[500];

    int result = 0;

    sprintf(command,"rm %s/rrpmref_fm.ref %s/rrpmref_pm.ref",CORONAGRAPHSDATALOCAL,CORONAGRAPHSDATALOCAL);
    result = system(command);
    delete_image_ID("rrpmref_fm.ref");
    delete_image_ID("rrpmref_pm.ref");
    coronagraph_simul_RRPM(0.0,0.0,"psf");
    best_value = arith_image_total("psf");
    delete_image_ID("psf");


    while(0==0)
    {
        RRPM_RADIUS += 0.01*pow((0.5-ran1()),3.0);
        RRPM_P2 += 0.0005*pow((0.5-ran1()),3.0);
        RRPM_P3 += 0.0001*pow((0.5-ran1()),3.0);
        RRPM_P4 += 0.0001*pow((0.5-ran1()),3.0);
        RRPM_P5 += 0.0001*pow((0.5-ran1()),3.0);

        sprintf(command,"rm %s/rrpmref_fm.ref %s/rrpmref_pm.ref",CORONAGRAPHSDATALOCAL,CORONAGRAPHSDATALOCAL);
        result = system(command);
        delete_image_ID("rrpmref_fm.ref");
        delete_image_ID("rrpmref_pm.ref");
        coronagraph_simul_RRPM(0.0,0.0,"psf");
        value = arith_image_total("psf");
        delete_image_ID("psf");

        printf("[%g %g %g %g %g] %g     [%g %g %g %g %g] %g\n",RRPM_RADIUS,RRPM_P2,RRPM_P3,RRPM_P4,RRPM_P5,value,best_RRPM_RADIUS,best_RRPM_P2,best_RRPM_P3,best_RRPM_P4,best_RRPM_P5,best_value);
        if(value<best_value)
        {
            best_RRPM_RADIUS = RRPM_RADIUS;
            best_RRPM_P2 = RRPM_P2;
            best_RRPM_P3 = RRPM_P3;
            best_RRPM_P4 = RRPM_P4;
            best_RRPM_P5 = RRPM_P5;
            best_value = value;
        }
        else
        {
            RRPM_RADIUS = best_RRPM_RADIUS;
            RRPM_P2 = best_RRPM_P2;
            RRPM_P3 = best_RRPM_P3;
            RRPM_P4 = best_RRPM_P4;
            RRPM_P5 = best_RRPM_P5;
        }

    }

    return(0);
}


int coronagraph_simul_MULTISTEP_APLC(double xld, double yld, const char *psfname)
{
    FILE *fp;
    int result = 0;
    long size = CORONAGRAPHS_ARRAYSIZE;
    char fname[200];
    char fname1[200];
    double trad_pix;
    long IDa1,IDp1,ID1;
    long IDapo = -1;
    long ii,jj;
    long ID;
    long IDfpm;
    long IDpupa0;
    double total0;
    long size2 = size*size;
    long iter;
    long IDout;
    long IDb, IDc, ID2, ID3, ID4;
    double tmp1, tmp2, tmp3, tmp4;
    int CHEAT = 0;
    long index1,index2;
    double re1,im1,re2,im2,re,im,r;
    long ii1,jj1,ii2,jj2;
    double tot1, tot2;
    double v1, v2;
    int FPMZONES;
    double MASKCAMULT;
    long kx, ky;
    int ok1;
    int WriteFiles = 1;
    struct stat info;
    long IDlm;

    long IDfpm_a, IDfpm_p;

    long iisize, jjsize;
    long iico, jjfpmrad;

    long n;



    if((ID=variable_ID("DFTZFACTOR"))!=-1)
        DFTZFACTOR = data.variable[ID].value.f;


    if(useDFT==0)
        DFTZFACTOR = 1;
    else
        printf("DFTZFACTOR = %f\n", DFTZFACTOR);


    ID = variable_ID("WRITEFILES");
    if(ID != -1)
        WriteFiles = 1;

    printf("SIZE = %ld\n",size);
    fflush(stdout);





    //
    // READ APODIZATION PARAMETERS FOR PURELY CENTRALLY OBSCURED CIRCULAR PUPIL
    //

    IDprol_init = image_ID("pinit");
    if(IDprol_init == -1)
    {
        sprintf(fname, "%s/APLCapo_init.fits.gz", CORONAGRAPHSDATADIR);
        IDprol_init = load_fits(fname, "pinit", 1);
    }

    IDprol_ffrac = image_ID("pffrac");
    if(IDprol_ffrac == -1)
    {
        sprintf(fname, "%s/APLCapo_ffrac.fits.gz", CORONAGRAPHSDATADIR);
        IDprol_ffrac = load_fits(fname, "pffrac", 1);
    }

    IDprol_transm = image_ID("ptransm");
    if(IDprol_transm == -1)
    {
        sprintf(fname, "%s/APLCapo_transm.fits.gz", CORONAGRAPHSDATADIR);
        IDprol_transm = load_fits(fname,"ptransm", 1);
    }

    IDprol_peak = image_ID("ppeak");
    if(IDprol_peak == -1)
    {
        sprintf(fname, "%s/APLCapo_peak.fits.gz", CORONAGRAPHSDATADIR);
        IDprol_peak = load_fits(fname, "ppeak", 1);
    }



    IDprol_fitapo_a = image_ID("fitapoa");
    if(IDprol_fitapo_a == -1)
    {
        sprintf(fname, "%s/APLCapo_fitapo_a.fits.gz", CORONAGRAPHSDATADIR);
        IDprol_fitapo_a = load_fits(fname, "fitapoa", 1);
    }

    IDprol_fitapo_b = image_ID("fitapob");
    if(IDprol_fitapo_b == -1)
    {
        sprintf(fname, "%s/APLCapo_fitapo_b.fits.gz", CORONAGRAPHSDATADIR);
        IDprol_fitapo_b = load_fits(fname, "fitapob", 1);
    }

    IDprol_fitapo_c = image_ID("fitapoc");
    if(IDprol_fitapo_c == -1)
    {
        sprintf(fname, "%s/APLCapo_fitapo_c.fits.gz", CORONAGRAPHSDATADIR);
        IDprol_fitapo_c = load_fits(fname, "fitapoc", 1);
    }




    IDprol_fitfit = image_ID("fitfit");
    if(IDprol_fitfit == -1)
    {
        sprintf(fname, "%s/APLCapo_fitfit.fits.gz", CORONAGRAPHSDATADIR);
        IDprol_fitfit = load_fits(fname, "fitfit", 1);
    }

    iisize = data.image[IDprol_init].md[0].size[0];
    jjsize = data.image[IDprol_init].md[0].size[1];


    iico = (long) ((APLC_CentOBS1-APLCapo_CO_START)/APLCapo_CO_STEP + 0.1);
    jjfpmrad = (long) ((APLC_FPMASKsize-APLCapo_FPMRAD_START)/APLCapo_FPMRAD_STEP + 0.1);


    printf("%ld %ld\n", iico, jjfpmrad);

    if((iico<0)||(iico>iisize-1))
    {
        printf("ERROR: iico out of range : %ld / %ld\n", iico, iisize);
        exit(0);
    }

    if((jjfpmrad<0)||(jjfpmrad>jjsize-1))
    {
        printf("ERROR: iico out of range : %ld / %ld\n", jjfpmrad, jjsize);
        exit(0);
    }




    if(data.image[IDprol_init].array.F[jjfpmrad*iisize+iico] > 0.1)
    {
        tmp4 = data.image[IDprol_peak].array.F[jjfpmrad*iisize+iico];
        MASKCAMULT = -(1.0-tmp4)/tmp4;
        printf("Focal plane mask CA = %f (peak = %f)\n", MASKCAMULT, tmp4);

        if(fitapoINIT==0)
        {
            fitapoN = data.image[IDprol_fitapo_a].md[0].size[2];
            fitapo_a = (double*) malloc(sizeof(double)*fitapoN);
            fitapo_b = (double*) malloc(sizeof(double)*fitapoN);
            fitapo_c = (double*) malloc(sizeof(double)*fitapoN);
            fitapo_c1 = (double*) malloc(sizeof(double)*fitapoN);
            fitapoINIT = 1;
        }
        for(n=0; n<fitapoN; n++)
        {
            fitapo_a[n] = data.image[IDprol_fitapo_a].array.F[n*iisize*jjsize+jjfpmrad*iisize+iico];
            fitapo_b[n] = data.image[IDprol_fitapo_b].array.F[n*iisize*jjsize+jjfpmrad*iisize+iico];
            fitapo_c[n] = data.image[IDprol_fitapo_c].array.F[n*iisize*jjsize+jjfpmrad*iisize+iico];
            fitapo_c1[n] = fitapo_c[n];
            if(fitapo_c1[n]<fitapo_minc)
                fitapo_c1[n] = fitapo_minc;
        }
        n = 0;
        printf("f(x) = %15.10f * exp(%15.10f * x**%15.10f)", fitapo_a[n], fitapo_b[n], fitapo_c1[n]);
        for(n=1; n<fitapoN; n++)
            printf("+ %15.10f * exp(%15.10f * x**%15.10f)", fitapo_a[n], fitapo_b[n], fitapo_c1[n]);
        printf("\n");


    }
    else
    {
        printf("NO DATA FOR THIS POINT\n");
        exit(0);
    }



/*
    
    if(APLC_CentOBS1>0.0001)
      sprintf(fname,"%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, APLC_FPMASKsize, APLC_CentOBS1, size);
    else
      sprintf(fname,"%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz", CORONAGRAPHSDATALOCAL, APLC_FPMASKsize, size);

    sprintf(fname1, "aplcapo");
    printf("Searching for Apodization file \"%s\" ... ", fname1);
    IDapo = image_ID(fname1);
    if(IDapo==-1)
      printf("not found\n");
    else
      printf("found\n");


    ID = image_ID("apomask");

    if(ID==-1)
      {
        if(APLC_CentOBS1>0.0001)
    sprintf(PIAAAPO_FNAME, "APLCapo_%.3f.%.3f.%ld.ref.prof", APLC_FPMASKsize, APLC_CentOBS1, size);
        else
    sprintf(PIAAAPO_FNAME, "APLCapo_%.3f.%ld.ref.prof", APLC_FPMASKsize, size);

        if(APLC_CentOBS1>0.0001)
    sprintf(PIAAAPODIZE_2DAPOFNAME, "APLCapo_%.3f.%.3f.%ld.ref.gz", APLC_FPMASKsize, APLC_CentOBS1, size);
        else
    sprintf(PIAAAPODIZE_2DAPOFNAME, "APLCapo_%.3f.%ld.ref.gz", APLC_FPMASKsize, size);

        if(IDapo==-1)
    {
      if(file_exists(fname)==1)
        {
          printf("LOADING APODIZATION FILE \"%s\"\n",fname);
          IDapo = load_fits(fname, fname1);
        }
      else
        {
          if(APLC_CentOBS1<0.01)
    	{
    	  printf("CREATE APODIZATION\n");
    	  coronagraph_make_2Dprolate_DFT(APLC_FPMASKsize/CORONAGRAPHS_PIXSCALE, 0.0, fname1);
    	  save_fl_fits(fname1, fname);
    	  IDapo = image_ID(fname1);
    	}
          else
    	{
    	  printf("APODIZATION \"%s\" MISSING - CANNOT CONTINUE\n", fname);
    	  exit(0);
    	}
        }
    }


        sprintf(fname, "APLCapo.%.3f.%.3f.info", APLC_FPMASKsize, APLC_CentOBS1);
        if(stat(fname, &info)!=0)
    {
      if(APLC_CentOBS1>0.0001)
        sprintf(fname,"%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, APLC_FPMASKsize, APLC_CentOBS1, size);
      else
        sprintf(fname,"%s/APLCapo/APLCapo_%.3f.%ld.ref.info", CORONAGRAPHSDATALOCAL, APLC_FPMASKsize, size);
    }


        if((fp = fopen(fname,"r"))==NULL)
    {
      printf("ERROR: cannot open file \"%s\"\n",fname);
      exit(0);
    }
        result = fscanf(fp,"%lf %lf %lf %lf\n", &tmp1, &tmp2, &tmp3, &tmp4);
        fclose(fp);
        if((NB_APLC_STEP==1)&&(APLC_PMASK==1))
    MASKCAMULT = -(1.0-tmp4)/tmp4;
        else
    MASKCAMULT = 0.0;
        if(MASKCAMULT<-1.0)
    MASKCAMULT = -1.0;

        printf("Focal plane mask CA = %f (%f)\n", MASKCAMULT, tmp4);
      }
    */









    //
    // FOCAL PLANE MASK
    //

    if(FPMASKSIZE_ERROR==0)
    {
        //      printf("--------------- NOMINAL FP SIZE ------------\n");
        //  sprintf(fname,"%s/APLCfpmask_%f.%ld.ref", CORONAGRAPHSDATALOCAL, APLC_FPMASKsize, size);
        sprintf(fname1, "aplcfpm");
        IDfpm = image_ID(fname1);

        if(IDfpm==-1) // if focal plane mask not already in memory, create it
        {
            //	  if(file_exists(fname)==1)
            //  IDfpm = load_fits(fname, fname1);
            // else
            // {
            printf("MASK RADIUS = %f pix\n",DFTZFACTOR*APLC_FPMASKsize/CORONAGRAPHS_PIXSCALE );
            IDfpm = make_subpixdisk(fname1, size, size, size/2, size/2, DFTZFACTOR*APLC_FPMASKsize/CORONAGRAPHS_PIXSCALE);

            for(ii=0; ii<size2; ii++)
                data.image[IDfpm].array.F[ii] = 1.0-(1.0-MASKCAMULT*(1.0+FPMASK_transm_error))*data.image[IDfpm].array.F[ii];

            save_fl_fits(fname1, "!FPmask.fits");
            IDfpm = image_ID(fname1);
            //}
        }
    }
    else
    {
        sprintf(fname1,"aplcfpm");
        IDfpm = image_ID(fname1);
        if(IDfpm==-1)
        {
            FPMZONES = 1;
            printf("-------------- FPMASKSIZE_ERROR = 1, FPMASK_FACTOR = %f, FPMZONES = %d ----------- (%f %f)\n", FPMASK_FACTOR, FPMZONES, FPMASK_transm_error, FPMASK_size_error);
            if(FPMZONES==1) // if focal plane mask is a single zone
            {
                IDfpm = make_subpixdisk("aplcfpm", size, size, size/2, size/2, DFTZFACTOR*APLC_FPMASKsize*FPMASK_FACTOR/CORONAGRAPHS_PIXSCALE*(1.0+FPMASK_size_error));
                //    for(ii=0;ii<size2;ii++)
                //data.image[IDfpm].array.F[ii] = 1.0-(1.0-MASKCAMULT*(1.0+FPMASK_transm_error))*data.image[IDfpm].array.F[ii];

                save_fl_fits("aplcfpm", "!fpm.test.fits");
                //	      sprintf(fname, "!%s/APLCfpmask_%f.%ld.ref", CORONAGRAPHSDATALOCAL, APLC_FPMASKsize, size);
                //	      save_fl_fits("aplcfpm", fname);
                IDfpm = image_ID("aplcfpm");
            }
            else // 2 zones mask
            {
                // here, FPMASK_FACTOR1 < FPMASK_FACTOR2
                printf("--------------- 2 zones mask  -------------\n");
                IDfpm = make_subpixdisk("aplcfpm", size, size, size/2, size/2, DFTZFACTOR*APLC_FPMASKsize*FPMASK_FACTOR1/CORONAGRAPHS_PIXSCALE);
                ID2 = make_subpixdisk("aplcfpm2", size, size, size/2, size/2, DFTZFACTOR*APLC_FPMASKsize*FPMASK_FACTOR2/CORONAGRAPHS_PIXSCALE);
                for(ii=0; ii<size2; ii++)
                {
                    v1 = data.image[IDfpm].array.F[ii];
                    v2 = data.image[ID2].array.F[ii];
                    data.image[IDfpm].array.F[ii] = 1.0-v1*(1.0-FPM_TRANSM1)-(v2-v1)*(1.0-FPM_TRANSM2);
                }
                delete_image_ID("aplcfpm2");
                save_fl_fits("aplcfpm","!aplcfpm");
            }
        }
    }




    // Make input pupil
    printf("MAKE INPUT PUPIL\n");

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    printf("RADIUS = %f\n", trad_pix);
    IDa1 = make_disk("pa1", size, size, size/2, size/2, trad_pix+2.0);

    //  save_fl_fits("pa1", "!pa110.fits");
    //  list_image_ID();

    IDpupa0 = make_disk("pa0", size, size, size/2, size/2, trad_pix);

    //  total0 = arith_image_total("pa1");

    IDp1 = make_slopexy("pp1", size, size, PI*xld/trad_pix, PI*yld/trad_pix);

    ID = image_ID("corpha");
    if(ID!=-1)
    {
        printf("Applying pre phase errors\n");
        for(ii=0; ii<size2; ii++)
            data.image[IDp1].array.F[ii] += data.image[ID].array.F[ii];
    }

    ID = image_ID("coramp");
    if(ID!=-1)
    {
        printf("Applying pre amplitude errors\n");
        for(ii=0; ii<size2; ii++)
            data.image[IDa1].array.F[ii] *= data.image[ID].array.F[ii];
    }

    total0 = 0.0;
    for(ii=0; ii<size2; ii++)
        total0 += data.image[IDa1].array.F[ii]*data.image[IDa1].array.F[ii];



    if(APLC_PIAA==0)
    {
        ID = image_ID("apomask");
        if(ID!=-1)
        {
            printf("Applying amplitude apodization mask\n");
            for(ii=0; ii<size2; ii++)
                data.image[IDa1].array.F[ii] *= data.image[ID].array.F[ii];
        }
        else // default apodization
        {
			if(IDapo!=-1)
				for(ii=0; ii<size2; ii++)
					data.image[IDa1].array.F[ii] *= data.image[IDapo].array.F[ii];
        }

        mk_complex_from_amph("pa1","pp1","pc1", 0);
        delete_image_ID("pa1");
        delete_image_ID("pp1");


        // AUTO DEFINE LYOT MASK
        ID = image_ID("coramp");
        if(ID!=-1)
        {
            for(ii=3; ii<size-3; ii++)
                for(jj=3; jj<size-3; jj++)
                {
                    ok1 = 1;
                    kx = 0;
                    ky = 0;
                    if(data.image[ID].array.F[(jj+ky)*size+ii+kx]<0.001)
                        ok1 = 0;
                    if(ok1==0)
                        data.image[IDpupa0].array.F[jj*size+ii] = 0.0;
                }
        }
    }
    else
    {
        PIAAAPO_NBPOINTS = (long) (CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0);
        if(initPIAA==0)
        {
            // initialize PIAA using PIAAAPO_FNAME apodization profile
            // this creates the relationship between input and output pupil radii
            coronagraph_init_PIAA();
            initPIAA = 1;
        }
        //    save_fl_fits("pa1","!pa111.fits");

        coronagraphs_PIAA_apodize_beam("pa1", "pp1", "pa1a", "pp1a");

        //create_2Dimage_ID("pa1a", size, size);
        //create_2Dimage_ID("pp1a", size, size);
        //      save_fl_fits("pa1a", "!pa112.fits");

        // AUTO DEFINE LYOT MASK
        ID1 = image_ID("pa1a");
        if(WriteFiles==1)
        {
            save_fl_fits("pa1","!pa1.fits");
            save_fl_fits("pp1","!pp1.fits");

            save_fl_fits("pa1a","!pa1a.fits");
            save_fl_fits("pp1a","!pp1a.fits");
        }


        /*      ID = image_ID("coramp");
        if(ID!=-1)
        {
          for(ii=3;ii<size-3;ii++)
            for(jj=3;jj<size-3;jj++)
              {
        	ok1 = 1;
        	kx = 0;
        	ky = 0;
        	if(data.image[ID].array.F[(jj+ky)*size+ii+kx]<0.001)
        	  ok1 = 0;
        	if(ok1==0)
        	  data.image[IDpupa0].array.F[jj*size+ii] = 0.0;
              }
        }
             */

        delete_image_ID("pa1");
        delete_image_ID("pp1");

        ID1 = image_ID("pa1a");
        ID = image_ID("postPIAAmask");
        if(ID!=-1)
        {
            printf("Applying postPIAA amplitude mask\n");
            for(ii=0; ii<size2; ii++)
                data.image[ID1].array.F[ii] *= data.image[ID].array.F[ii];
            if(WriteFiles==1)
            {
                save_fl_fits("pa1a","!pa1a_post.fits");
                save_fl_fits("pp1a","!pp1a_post.fits");
            }
        }

        mk_complex_from_amph("pa1a","pp1a","pc1", 0);
        if(WriteFiles==1)
        {
            save_fl_fits("pa1a","!pa1a");
            save_fl_fits("pp1a","!pp1a");
        }
        delete_image_ID("pa1a");
        delete_image_ID("pp1a");
    }



    for(iter=0; iter<NB_APLC_STEP; iter++)
    {
        //
        // pupil -> pupil propagation with focal plane mask
        //


        if(useDFT==0) // regular FFT propagation
        {
            IDfpm = image_ID("aplcfpm");

            //  for(ii=0;ii<size2;ii++)
            // data.image[IDfpm].array.F[ii] = 1.0-(1.0-MASKCAMULT*(1.0+FPMASK_transm_error))*data.image[IDfpm].array.F[ii];

            permut("pc1");
            do2dfft("pc1","fc1");
            delete_image_ID("pc1");
            permut("fc1");
            ID = image_ID("fc1");
            tot1 = 0.0;
            tot2 = 0.0;
            for(ii=0; ii<size2; ii++)
            {
                re = data.image[ID].array.CF[ii].re;
                im = data.image[ID].array.CF[ii].im;
                tot1 += re*re+im*im;
                data.image[ID].array.CF[ii].re *= data.image[IDfpm].array.F[ii];
                data.image[ID].array.CF[ii].im *= data.image[IDfpm].array.F[ii];
                re = data.image[ID].array.CF[ii].re;
                im = data.image[ID].array.CF[ii].im;
                tot2 += re*re+im*im;
            }


            mk_amph_from_complex("fc1","fc1_a","fc1_p", 0);

            if(WriteFiles==1)
            {
                save_fl_fits("aplcfpm","!aplcfpm.fits");
                save_fl_fits("fc1_a","!fc1_a.fits");
                save_fl_fits("fc1_p","!fc1_p.fits");
            }

            delete_image_ID("fc1_a");
            delete_image_ID("fc1_p");

            //      printf("Focal plane mask %ld : throughput = %g\n",iter,tot2/tot1);
            permut("fc1");
            do2dfft("fc1","pc1");
            delete_image_ID("fc1");
            permut("pc1");

            if(WriteFiles==1)
            {
                mk_amph_from_complex("pc1","pc1a","pc1p", 0);
                save_fl_fits("pc1a","!pc1a.fits");
                save_fl_fits("pc1p","!pc1p.fits");
                delete_image_ID("pc1a");
                delete_image_ID("pc1p");
            }

            ID = image_ID("pc1");


            // APPLY LYOT MASK
            IDlm = image_ID("LyotMask");
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    if((ii<5)||(jj<5)||(ii>size-5)||(jj>size-5))
                    {
                        data.image[ID].array.CF[jj*size+ii].re = 0.0;
                        data.image[ID].array.CF[jj*size+ii].im = 0.0;
                    }
                    else
                    {
                        data.image[ID].array.CF[jj*size+ii].re *= data.image[IDlm].array.F[(size-jj)*size+(size-ii)]/size2;
                        data.image[ID].array.CF[jj*size+ii].im *= data.image[IDlm].array.F[(size-jj)*size+(size-ii)]/size2;
                        //			data.image[ID].array.CF[jj*size+ii].re /= size2;
                        //			data.image[ID].array.CF[jj*size+ii].im /= size2;
                    }
                }
            if(WriteFiles==1)
            {
                //	      save_fl_fits("pa0","!pa0_mask.fits");
                mk_amph_from_complex("pc1","pc1am","pc1pm", 0);
                save_fl_fits("pc1am","!pc1am.fits");
                save_fl_fits("pc1pm","!pc1pm.fits");
                delete_image_ID("pc1am");
                delete_image_ID("pc1pm");
            }
        }
        else // DFT propagation
        {
            printf("STARTING DFT  [%f %f]\n", MASKCAMULT, -(1.0-MASKCAMULT));
            fflush(stdout);

            if(WriteFiles==1)
                save_fl_fits("aplcfpm", "!TEST_aplcfpm.fits");
            ID = image_ID("aplcfpm");
            ID1 = create_2Dimage_ID("aplcfpm1", size, size);
            for(ii=0; ii<size2; ii++)
                data.image[ID1].array.F[ii] = 1.0-data.image[ID].array.F[ii];
            if(WriteFiles==1)
                save_fl_fits("aplcfpm1", "!TEST_aplcfpm1.fits");

            IDfpm_a = image_ID("aplcfpm_a");
            IDfpm_p = image_ID("aplcfpm_p");
            if((IDfpm_a!=-1)&&(IDfpm_p!=-1))
            {
                printf("============= USING CUSTOM COMPLEX FOCAL PLANE MASK ============\n");
                mk_complex_from_amph("aplcfpm_a", "aplcfpm_p", "aplcfpm_c", 0);
                ID = image_ID("aplcfpm_c");
                for(ii=0; ii<size2; ii++)
                {
                    data.image[ID].array.CF[ii].re = 1.0 - data.image[ID].array.CF[ii].re;
                    data.image[ID].array.CF[ii].im = -data.image[ID].array.CF[ii].im;
                }
                fft_DFTinsertFPM( "pc1", "aplcfpm_c", DFTZFACTOR, "pc2");
                delete_image_ID("aplcfpm_c");
            }
            else
                fft_DFTinsertFPM_re( "pc1", "aplcfpm1", DFTZFACTOR, "pc2");
            //	  list_image_ID();

            if(0) // testing
            {
                mk_amph_from_complex("pc1", "dftia", "dftip", 0);
                save_fl_fits("dftia", "!dftia.fits");
                save_fl_fits("dftip", "!dftip.fits");
                delete_image_ID("dftia");
                delete_image_ID("dftip");
                mk_amph_from_complex("pc2", "dftoa", "dftop", 0);
                save_fl_fits("dftoa", "!dftoa.fits");
                save_fl_fits("dftop", "!dftop.fits");
                delete_image_ID("dftoa");
                delete_image_ID("dftop");
                //mk_amph_from_complex("aplcfpm1", "dftfpma", "dftfpmp", 0);
                save_fl_fits("aplcfpm1", "!dftfpma.fits");
                //	      save_fl_fits("dftfpmp", "!dftfpmp.fits");
                // delete_image_ID("dftfpma");
                //delete_image_ID("dftfpmp");
            }
            delete_image_ID("aplcfpm1");
            //	  exit(0);


            ID1 = image_ID("pc1");
            ID2 = image_ID("pc2");
            for(ii=0; ii<size2; ii++)
            {
                re1 = data.image[ID1].array.CF[ii].re;
                im1 = data.image[ID1].array.CF[ii].im;

                re2 = data.image[ID2].array.CF[ii].re;
                im2 = data.image[ID2].array.CF[ii].im;

                re1 = re1-(1.0-0.0*MASKCAMULT)*re2;
                im1 = im1-(1.0-0.0*MASKCAMULT)*im2;

                data.image[ID1].array.CF[ii].re = re1;
                data.image[ID1].array.CF[ii].im = im1;
            }
            delete_image_ID("pc2");
            if(WriteFiles==1)
            {
                mk_amph_from_complex("pc1","pc1a","pc1p", 0);
                save_fl_fits("pc1a","!pc1a_01.fits");
                save_fl_fits("pc1p","!pc1p_01.fits");
                delete_image_ID("pc1a");
                delete_image_ID("pc1p");
            }
        }
    }


    if(APLC_PIAA==2)  // INVERSE PIAA
    {
        mk_amph_from_complex("pc1","pa3","pp3", 0);

        if(WriteFiles==1)
            save_fl_fits("pa3","!pa3.fits");

        copy_image_ID("pp3", "pp3b", 0);
        IDb = image_ID("pp3b");
        for(ii=0; ii<size*size; ii++)
        {
            tmp1 = data.image[IDb].array.F[ii];
            tmp1 += 2.0*PI/3.0;
            if(tmp1>PI)
                tmp1 -= 2.0*PI;
            data.image[IDb].array.F[ii] = tmp1;
        }
        copy_image_ID("pp3", "pp3c", 0);
        IDc = image_ID("pp3c");
        for(ii=0; ii<size*size; ii++)
        {
            tmp1 = data.image[IDc].array.F[ii];
            tmp1 -= 2.0*PI/3.0;
            if(tmp1<-PI)
                tmp1 += 2.0*PI;
            data.image[IDc].array.F[ii] = tmp1;
        }

        piaaconfdirection = 1;

        if(0) // testing
        {
            execute_arith("zeroim=0*postPIAAmask");
            coronagraphs_PIAA_apodize_beam("postPIAAmask", "postPIAAmask", "prePIAAmask1", "prePIAAmask");
            save_fl_fits("postPIAAmask", "!test_postPIAAmask.fits");
            save_fl_fits("prePIAAmask", "!test_prePIAAmask.fits");
        }

        coronagraphs_PIAA_apodize_beam("pa3", "pp3", "pa4", "pp4");
        coronagraphs_PIAA_apodize_beam("pa3", "pp3b", "pa4b", "pp4b");
        coronagraphs_PIAA_apodize_beam("pa3", "pp3c", "pa4c", "pp4c");
        piaaconfdirection = 0;

        ID = image_ID("pp4");
        IDb = image_ID("pp4b");
        IDc = image_ID("pp4c");
        for(ii=0; ii<size*size; ii++)
        {
            tmp1 = data.image[ID].array.F[ii];
            tmp2 = data.image[IDb].array.F[ii];
            tmp3 = data.image[IDc].array.F[ii];

            if((fabs(tmp2-(2.0*PI/3.0-PI))<0.3*PI)&&(fabs(tmp3+(2.0*PI/3.0-PI))<0.3*PI))
            {
                tmp1 = PI+ 0.5*(tmp2+tmp3);
                if(tmp1>PI)
                    tmp1 -= 2.0*PI;
                data.image[ID].array.F[ii] = tmp1;
            }
        }

        delete_image_ID("pp3b");
        delete_image_ID("pp3c");
        delete_image_ID("pa4b");
        delete_image_ID("pp4b");
        delete_image_ID("pa4c");
        delete_image_ID("pp4c");


        printf(" CENT OBS 0 = %f x %f -------------------\n", APLC_CentOBS0, trad_pix);

        if(CHEAT==1)
        {
            printf("------------CHEAT----------\n");
            ID1 = image_ID("pa1");
            ID2 = image_ID("pp1");
            ID3 = image_ID("pa4");
            ID4 = image_ID("pp4");

            for(ii=10; ii<size-10; ii++)
                for(jj=10; jj<size-10; jj++)
                {
                    index1 = jj*size+ii;
                    index2 = (size-jj)*size+(size-ii);
                    re1 = data.image[ID1].array.F[index2]*cos(data.image[ID2].array.F[index2]);
                    im1 = data.image[ID1].array.F[index2]*sin(data.image[ID2].array.F[index2]);
                    re2 = data.image[ID3].array.F[index1]*cos(data.image[ID4].array.F[index1]);
                    im2 = data.image[ID3].array.F[index1]*sin(data.image[ID4].array.F[index1]);

                    re = re1*size*size-re2;
                    im = im1*size*size-im2;

                    tmp1 = 1.0*ii-size/2;
                    tmp2 = 1.0*jj-size/2;
                    r = sqrt(tmp1*tmp1+tmp2*tmp2)/trad_pix;
                    if((r>PIAAPPMASKRAD)||(r<APLC_CentOBS0))
                    {
                        re = 0.0;
                        im = 0.0;
                    }

                    data.image[ID3].array.F[index1] = sqrt(re*re+im*im);
                    data.image[ID4].array.F[index1] = atan2(im,re);
                }
        }
        else
        {
            ID3 = image_ID("pa4");
            ID4 = image_ID("pp4");
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    index1 = jj*size+ii;
                    re = data.image[ID3].array.F[index1]*cos(data.image[ID4].array.F[index1]);
                    im = data.image[ID3].array.F[index1]*sin(data.image[ID4].array.F[index1]);
                    tmp1 = 1.0*ii-size/2;
                    tmp2 = 1.0*jj-size/2;
                    r = sqrt(tmp1*tmp1+tmp2*tmp2);
                    if((r>PIAAPPMASKRAD*trad_pix)||(r<APLC_CentOBS0*trad_pix+3.0)) // 3 pixel buffer to avoid edge effects
                    {
                        re = 0.0;
                        im = 0.0;
                    }

                    data.image[ID3].array.F[index1] = sqrt(re*re+im*im);
                    data.image[ID4].array.F[index1] = atan2(im,re);
                }
        }

        mk_complex_from_amph("pa4","pp4","pc4", 0);
        if(WriteFiles==1)
        {
            save_fl_fits("pa4","!pa4.fits");
            save_fl_fits("pp4","!pp4.fits");
        }

        permut("pc4");
        do2dfft("pc4","fc4");
        permut("fc4");

        if(APLC_FLIP==1)
        {
            printf("-------- FLIP --------\n");
            ID = image_ID("fc4");
            for(ii=1; ii<size; ii++)
                for(jj=1; jj<size/2+1; jj++)
                {
                    ii1 = ii;
                    jj1 = jj;
                    ii2 = size-ii;
                    jj2 = size-jj;
                    re1 = data.image[ID].array.CF[jj1*size+ii1].re;
                    im1 = data.image[ID].array.CF[jj1*size+ii1].im;
                    re2 = data.image[ID].array.CF[jj2*size+ii2].re;
                    im2 = data.image[ID].array.CF[jj2*size+ii2].im;
                    re = 0.5*(re1+re2);
                    im = 0.5*(im1+im2);
                    data.image[ID].array.CF[jj1*size+ii1].re = re;
                    data.image[ID].array.CF[jj1*size+ii1].im = im;
                    data.image[ID].array.CF[jj2*size+ii2].re = re;
                    data.image[ID].array.CF[jj2*size+ii2].im = im;
                }
        }
        mk_amph_from_complex("fc4", "fa4", "fp4", 0);

        arith_image_mult("fa4", "fa4", psfname);
        ID = image_ID(psfname);
        for(ii=0; ii<size2; ii++)
            data.image[ID].array.F[ii] /= size2*total0;

        delete_image_ID("pa3");
        delete_image_ID("pp3");
        delete_image_ID("pc4");
        delete_image_ID("pa4");
        delete_image_ID("pp4");
        delete_image_ID("fc4");
        if(WriteFiles==1)
        {
            save_fl_fits("fa4","!fa4.fits");
            save_fl_fits("fp4","!fp4.fits");
        }
        delete_image_ID("fa4");
        delete_image_ID("fp4");
    }
    else
    {
        do2dfft("pc1","fc2");
        permut("fc2");

        if(APLC_FLIP==1)
        {
            ID = image_ID("fc2");
            for(ii=1; ii<size; ii++)
                for(jj=1; jj<size/2+1; jj++)
                {
                    ii1 = ii;
                    jj1 = jj;
                    ii2 = size-ii;
                    jj2 = size-jj;
                    re1 = data.image[ID].array.CF[jj1*size+ii1].re;
                    im1 = data.image[ID].array.CF[jj1*size+ii1].im;
                    re2 = data.image[ID].array.CF[jj2*size+ii2].re;
                    im2 = data.image[ID].array.CF[jj2*size+ii2].im;
                    re = 0.5*(re1+re2);
                    im = 0.5*(im1+im2);
                    data.image[ID].array.CF[jj1*size+ii1].re = re;
                    data.image[ID].array.CF[jj1*size+ii1].im = im;
                    data.image[ID].array.CF[jj2*size+ii2].re = re;
                    data.image[ID].array.CF[jj2*size+ii2].im = im;
                }
        }

        mk_amph_from_complex("fc2","fa2","fp2", 0);

        ID = image_ID("fa2");
        delete_image_ID("fp2");
        IDout = create_2Dimage_ID(psfname,size,size);
        for(ii=0; ii<size2; ii++)
            data.image[IDout].array.F[ii] = data.image[ID].array.F[ii]/size2/total0*data.image[ID].array.F[ii];

        delete_image_ID("fc2");
        delete_image_ID("fa2");
    }

    delete_image_ID("pc1");
    delete_image_ID("pa0");

    // list_image_ID();

    return(0);
}





int coronagraph_init_EXTERNAL_OCCULTER(double D, double l, double lambda, long FACTOR)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    double PARAM_a = 12.5;
    double PARAM_b = 12.5;
    long PARAM_n = 6;
    double pixscale;
    double trad_pix;
    long ID;

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    pixscale = D/2.0/trad_pix;
    make_offsetHyperGaussian(size*FACTOR,trad_pix/(D/2.0)*PARAM_a,trad_pix/(D/2.0)*PARAM_b,PARAM_n,"extocc");
    printf("pix:  %g  %g\n",trad_pix/(D/2.0)*PARAM_a,trad_pix/(D/2.0)*PARAM_b);
    save_fl_fits("extocc","!extocc");
    create_2Dimage_ID("zero",size*FACTOR,size*FACTOR);
    mk_complex_from_amph("extocc","zero","extocc_c", 0);
    printf("pixscale = %g    l = %g    lambda = %g\n",pixscale,l,lambda);
    Fresnel_propagate_wavefront("extocc_c","extoccp",pixscale,l,lambda);
    ID = image_ID("extoccp");

    delete_image_ID("zero");
    delete_image_ID("extocc");
    delete_image_ID("extocc_c");

    return(ID);
}

int coronagraph_simul_EXTERNAL_OCCULTER(double xld, double yld, const char *psfname)
{
    long size = CORONAGRAPHS_ARRAYSIZE;
    long size2 = size*size;
    double D = 6.5; /* in m */
    double l = 50000000.0; /* in m */
    double lambda = 0.0000005; /* in m */
    long IDocc,IDocca,IDoccp;
    double shadow_pixscale; /* rad/pixel */
    double pixscale;
    double trad_pix;
    double xrad,yrad;
    long IDa1,IDp1;
    long ii,jj,ii1,jj1;
    long xtransl,ytransl;
    double r;
    long ID,IDout;
    double total0;
    long FACTOR=2;

    IDocc = image_ID("extoccp");
    if(IDocc==-1)
    {
        IDocc = coronagraph_init_EXTERNAL_OCCULTER(D,l,lambda,FACTOR);
    }

    mk_amph_from_complex("extoccp","exta","extp", 0);
    IDocca = image_ID("exta");
    IDoccp = image_ID("extp");
    //  save_fl_fits("exta","!exta");


    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    pixscale = D/2.0/trad_pix; /* m/pixel */
    shadow_pixscale = pixscale/l; /* rad/pixel */

    xrad = xld*lambda/D;
    yrad = yld*lambda/D;

    printf("rad = %g %g   [%g]\n",xrad,yrad,shadow_pixscale);

    xtransl = (long) (xrad/shadow_pixscale);
    ytransl = (long) (yrad/shadow_pixscale);

    printf("transl = %ld %ld\n",xtransl,ytransl);

    IDa1 = create_2Dimage_ID("pa1",size,size);
    IDp1 = make_slopexy("pp1",size,size,PI*xld/trad_pix,PI*yld/trad_pix);

    total0 = 0.0;
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            ii1 = (ii-size/2)+size*FACTOR/2-xtransl;
            jj1 = (jj-size/2)+size*FACTOR/2-ytransl;
            if((ii1>0)&&(ii1<size*FACTOR)&&(jj1>0)&&(jj1<size*FACTOR))
            {
                data.image[IDa1].array.F[jj*size+ii] = data.image[IDocca].array.F[jj1*size*FACTOR+ii1];
                data.image[IDp1].array.F[jj*size+ii] += data.image[IDoccp].array.F[jj1*size*FACTOR+ii1];
            }
            else
            {
                data.image[IDa1].array.F[jj*size+ii] = 1.0;
            }
            r = sqrt((ii-size/2)*(ii-size/2)+(jj-size/2)*(jj-size/2));
            if(r>trad_pix)
                data.image[IDa1].array.F[jj*size+ii] = 0.0;
            else
                total0 += 1.0;
        }

    //  save_fl_fits("pa1","!pa1");
    mk_complex_from_amph("pa1","pp1","pc1", 0);
    do2dfft("pc1","fc2");
    permut("fc2");
    mk_amph_from_complex("fc2","fa2","fp2", 0);

    ID = image_ID("fa2");
    delete_image_ID("fp2");
    IDout = create_2Dimage_ID(psfname,size,size);
    for(ii=0; ii<size2; ii++)
        data.image[IDout].array.F[ii] = data.image[ID].array.F[ii]/size2/total0*data.image[ID].array.F[ii];

    delete_image_ID("exta");
    delete_image_ID("extp");
    delete_image_ID("pa1");
    delete_image_ID("pp1");
    delete_image_ID("pc1");
    delete_image_ID("fc2");
    delete_image_ID("fa2");

    return(0);
}



int coronagraph_simulPSF(double xld, double yld, const char *psfname, long coronagraph_type, const char *options)
{
    FILE *fp;
    int result = 0;
    int str_pos,i;
    char input[200];
    char str1[200];
    int ok1;
    double tmpf0, tmpf1;

    printf("SCALE : %f l/D per pix\n", CORONAGRAPHS_PIXSCALE);
    printf("CORONAGRAPHSDATADIR = %s\n", CORONAGRAPHSDATADIR);
    printf("CORONAGRAPHSDATALOCAL = %s\n", CORONAGRAPHSDATALOCAL);
    printf("CORONAGRAPHS_ARRAYSIZE = %d\n", CORONAGRAPHS_ARRAYSIZE);
    printf("OPTIONS = \"%s\"\n",options);




    /* this function can simulate various coronagraphs
    1 : AIC
    2 : phase mask
    3 : 4 quadrants
    4 : PIAA/CPA hybrid
    5 : PIAAC/CPA hybrid
    6 : 8 order BL, circular
    7 : 8 order BL, linear
    8 : AIC/PIAAC combination
    9 : CPA
    10 : APLC
    11 : APLC, 2 step
    12 : APLC, 3 step
    13 : APLC, 4 step
    14 : APLC, 5 step
    15 : APLC, 0 step (= CPA)
    16 : ODC
    19 : STRIPC
    20 : SUMXY
    */


    if (strstr(options,"-FPMSERR ")!=NULL)
    {
        str_pos=strstr(options,"-FPMSERR ")-options;
        str_pos = str_pos + strlen("-FPMSERR ");
        i=0;
        while((options[i+str_pos]!=' ')&&(options[i+str_pos]!='\n')&&(options[i+str_pos]!='\0'))
        {
            input[i] = options[i+str_pos];
            i++;
        }
        input[i] = '\0';
        FPMASKSIZE_ERROR = 1;
        FPMASK_size_error = atof(input);
        printf("-------------- MASK size error = %f\n",FPMASK_size_error);
    }



    switch (coronagraph_type)
    {
    case -1:
        coronagraph_simul_NOCORO(xld,yld,psfname);
        break;

    case 0:
        coronagraph_simul_DICC(xld,yld,psfname);
        break;

    case 1:
        coronagraph_simul_AIC(xld,yld,psfname);
        break;

    case 2: // Phase mask coronagraph
        coronagraph_simul_RRPM(xld,yld,psfname);
        break;

    case 3:
        coronagraph_simul_4QPM(xld,yld,psfname);
        break;

    case 4:
        coronagraph_simul_PIAA(xld,yld,psfname);
        break;

    case 5:
        PIAAFPMASKRAD = 5.5; // IWA ~ 2.2 l/D
        coronagraph_simul_PIAAC(xld,yld,psfname);
        break;

    case 6:
        BL8MODE = 0;
        coronagraph_simul_BL8(xld,yld,psfname);
        break;

    case 7:
        BL8MODE = 1;
        coronagraph_simul_BL8(xld,yld,psfname);
        break;

    case 8:
        coronagraph_simul_AIC_PIAAC(xld,yld,psfname);
        break;

    case 9:
        coronagraph_simul_CPA(xld,yld,psfname);
        break;

    case 10:
        NB_APLC_STEP = 1;
        APLC_FPMASKsize = 1.8; // 1.8 is the best mask size for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 11:
        NB_APLC_STEP = 2;
        APLC_FPMASKsize = 1.4; // 1.4 is the best mask size for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 12:
        NB_APLC_STEP = 3;
        APLC_FPMASKsize = 1.2; // 1.2 is the best mask size for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 13:
        NB_APLC_STEP = 4;
        APLC_FPMASKsize = 1.0; // 1.0 is the best mask size for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 14:
        NB_APLC_STEP = 5;
        APLC_FPMASKsize = 1.0; // 1.0 is the best mask size for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 15:
        NB_APLC_STEP = 0;
        APLC_FPMASKsize = 4.2; // this is the best mask size for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 16:
        coronagraph_simul_ODC(xld,yld,psfname);
        break;

    case 17:
        SHEAR4_SHEAR = 0.3;
        coronagraph_simul_SHEAR4(xld,yld,psfname);
        break;

    case 18:
        coronagraph_simul_BL4(xld,yld,psfname);
        break;

    case 19:
        STRIPCOFFSET = 0.4;
        coronagraph_simul_STRIPC(xld,yld,psfname);
        break;

    case 20:
        coronagraph_simul_SIMXY(xld,yld,psfname);
        break;

    case 21:
        PIAAFPMASKRAD = 5.0;
        coronagraph_simul_PIAAC(xld,yld,psfname);
        break;

    case 22:
        STRIPCOFFSET = 0.35;
        coronagraph_simul_STRIPC(xld,yld,psfname);
        break;

    case 23:
        STRIPCOFFSET = 0.45;
        coronagraph_simul_STRIPC(xld,yld,psfname);
        break;

    case 24: // SMEX baseline
        PIAAFPMASKRAD = 3.5;
        coronagraph_simul_PIAAC(xld,yld,psfname);
        break;

    case 25:
        OVC_CHARGE = 2;
        coronagraph_simul_OVC(xld,yld,psfname);
        break;

    case 26:
        OVC_CHARGE = 4;
        coronagraph_simul_OVC(xld,yld,psfname);
        break;

    case 27:
        OVC_CHARGE = 5;
        coronagraph_simul_OVC(xld,yld,psfname);
        break;

    case 28:
        OVC_CHARGE = 6;
        coronagraph_simul_OVC(xld,yld,psfname);
        break;

    case 29:
        OVC_CHARGE = 8;
        coronagraph_simul_OVC(xld,yld,psfname);
        break;

    case 30:
        coronagraph_simul_PPA(xld,yld,psfname);
        break;

    case 31:
        APLC_PIAA = 1;
        NB_APLC_STEP = 2;
        APLC_FPMASKsize = 1.6; // for 1e10 contrast with IWA = 1.2 l/D
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 32:
        APLC_PIAA = 1;
        NB_APLC_STEP = 3;
        APLC_FPMASKsize = 1.3; // for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 33:
        APLC_PIAA = 1;
        NB_APLC_STEP = 4;
        APLC_FPMASKsize = 1.2; // for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 34:
        APLC_PIAA = 1;
        NB_APLC_STEP = 2;
        APLC_FPMASKsize = 2.0; // for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 35:
        APLC_PIAA = 1;
        //      APLC_FLIP = 1;
        NB_APLC_STEP = 1;
        APLC_FPMASKsize = 2.4; // for 1e10 contrast
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    /* coronagraphs for 1e6 contrast - ground-based */

    case 36:  /* optimal for 1e6 contrast */
        APLC_PIAA = 0;
        NB_APLC_STEP = 0;
        APLC_FPMASKsize = 2.6;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 37:  /* optimal for 1e6 contrast */
        APLC_PIAA = 0;
        NB_APLC_STEP = 1;
        APLC_FPMASKsize = 1.4;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 38: /* optimal for 1e6 contrast */
        APLC_PIAA = 0;
        NB_APLC_STEP = 2;
        APLC_FPMASKsize = 1.0;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 39:  /* optimal for 1e6 contrast */
        APLC_PIAA = 1;
        NB_APLC_STEP = 0;
        APLC_FPMASKsize = 2.6;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 40:  // optimal for 1e6 contrast baselined for SMEX - but needs to check chromaticity
        APLC_PIAA = 1;
        NB_APLC_STEP = 1;
        APLC_FPMASKsize = 1.4;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 41: /* optimal for 1e6 contrast */
        APLC_PIAA = 1;
        NB_APLC_STEP = 2;
        APLC_FPMASKsize = 1.0;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    /* with central obscuration = 0.275 */

    case 42:  /* optimal for 1e6 contrast */
        PIAACENTOBS = 0.275;
        APLC_PIAA = 1;
        NB_APLC_STEP = 0;
        APLC_FPMASKsize = 2.6;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 43:  /* optimal for 1e6 contrast */
        PIAACENTOBS = 0.275;
        APLC_PIAA = 1;
        NB_APLC_STEP = 1;
        APLC_FPMASKsize = 1.4;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 44: /* optimal for 1e6 contrast */
        PIAACENTOBS = 0.275;
        APLC_PIAA = 1;
        NB_APLC_STEP = 2;
        APLC_FPMASKsize = 1.0;
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 45: // optimal for 1e6 contrast - BASELINE FOR EXCEDE
        APLC_PIAA = 2; // PIAA -> APLC -> inverse PIAA
        NB_APLC_STEP = 1;
        APLC_FPMASKsize = 2.0;
        fprintf(stderr,"MASK SIZE ERROR = %f\n",FPMASK_size_error);
        coronagraph_simul_MULTISTEP_APLC(xld,yld,psfname);
        break;

    case 50: // generalized PIAA APLC, read parameters from file
        APLC_PIAA = 1;
        APLC_PMASK = 1; // if NB_APLC_STEP = 1, put phase on focal plane mask
        APLC_CentOBS0 = 0;
        APLC_CentOBS1 = 0;
        fp = fopen("APLCPIAA.conf", "r");
        if(fp==NULL)
        {
            printERROR(__FILE__,__func__,__LINE__,"Cannot read configuration file \"APLCPIAA.conf\"\n");
            exit(0);
        }
        //     fscanf(fp,"%ld %f %d %d\n", &NB_APLC_STEP, &APLC_FPMASKsize, &APLC_PIAA, &APLC_PMASK);
        result = fscanf(fp,"%ld %s %d %d\n", &NB_APLC_STEP, str1, &APLC_PIAA, &APLC_PMASK);

        //      printf("str1 = \"%s\"\n", str1);
        ok1 = 0;
        for(i=0; i<strlen(str1); i++)
            if(str1[i]=='-')
            {
                ok1 = 1;
                str1[i] = ' ';
            }
        if(ok1==1)
        {
            sscanf(str1, "%lf %lf %lf", &APLC_FPMASKsize, &tmpf0, &tmpf1);
            APLC_CentOBS0 = tmpf0;
            APLC_CentOBS1 = tmpf1;
        }
        else
            sscanf(str1, "%lf", &APLC_FPMASKsize);
        fclose(fp);

        //printf("str1 = \"%s\"\n", str1);

        fprintf(stderr, "NB_APLC_STEP = %ld\n", NB_APLC_STEP); // number of APLC steps (1 for PIAACMC)
        fprintf(stderr, "APLC_FPMASKsize = %f\n", APLC_FPMASKsize); // focal plane mask size (l/D)
        fprintf(stderr, "APLC_CentOBS0 = %f\n", APLC_CentOBS0); // central obstruction in input
        fprintf(stderr, "APLC_CentOBS1 = %f\n", APLC_CentOBS1); // central obstruction in output
        fprintf(stderr, "APLC_PIAA = %d\n", APLC_PIAA); // 2: PIAA->APLC->invPIAA
        fprintf(stderr, "APLC_PMASK = %d\n", APLC_PMASK); // 1 if fpmask has phase

        if((fp=fopen("MASKerr.conf","r"))!=NULL)
        {
            result = fscanf(fp,"%lf %lf\n", &FPMASK_transm_error, &FPMASK_size_error);
            fclose(fp);
            fprintf(stderr,"MASK TRANSM ERROR = %f\n",FPMASK_transm_error);
            fprintf(stderr,"MASK SIZE ERROR = %f\n",FPMASK_size_error);
        }
        else
        {
            printf("No focal plane mask error file - using perfect focal plane mask\n");
            FPMASK_transm_error = 0.0;
            FPMASK_size_error = 0.0;
        }

        coronagraph_simul_MULTISTEP_APLC(xld, yld, psfname);

        break;

    case 60: /* external occulter */
        coronagraph_simul_EXTERNAL_OCCULTER(xld,yld,psfname);
        break;
    }

    return(0);
}


int coronagraph_transm(const char *fname, long coronagraph_type, double logcontrast, const char *options)
{
    FILE *fp;
    double sepstep = 0.05; // 0.05;
    double sepmax = 6.0; //8.0;
    double sep;
    double total,total1,total2,total3,total4;
    long i,j,ii,jj,kk;
    double x,y,r;
    long size = CORONAGRAPHS_ARRAYSIZE;
    long size2 = size*size;
    long ID,ID0;
    double flux0,flux;
    double TTerror = 0.001; /* in l/d */
    double FOCerror = 0.001; /* in l/d */
    double Zerror = 0.0;
    long Zindex = 0;
    int Zerror_included = 0;
    int TTerror_included = 0;
    int FOCerror_included = 0;
    char input[50];
    int str_pos;
    int nbTTpts = 20;
    int nbFOCpts = 10;
    long cnt;
    long ID1;
    double trad_pix;
    double tmp;
    double epsilon = 1.0e-10;
    double peak;
    double *arrayS;
    double *arrayS1;
    double *arrayS2;
    int RADAVER = 0;
    double angle;
    long anglestep;
    long NBanglestep = 16;
    double MaxAngle = 0.5*PI;
    double maxtransm ;
    double iwa,sepold,total2old;
    int iwaOK;

    y = 0.0;
    printf("CORONAGRAPH TYPE %ld    OPTION : \"%s\"\n",coronagraph_type,options);
    printf("LOGCONTRAST = %f\n",logcontrast);
    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;

    if (strstr(options,"-RADAVER")!=NULL)
    {
        RADAVER = 1;
        printf("RADIAL AVERAGE\n");
    }

    if (strstr(options,"-TT ")!=NULL)
    {
        str_pos=strstr(options,"-TT ")-options;
        str_pos = str_pos + strlen("-TT ");
        i=0;
        while((options[i+str_pos]!=' ')&&(options[i+str_pos]!='\n')&&(options[i+str_pos]!='\0'))
        {
            input[i] = options[i+str_pos];
            i++;
        }
        input[i] = '\0';
        TTerror = atof(input);
        if(TTerror>epsilon)
            TTerror_included = 1;
        printf("TT error is %f\n",TTerror);
        /*      if(TTerror<0.02)
        PIAAFPMASKRAD = 4.6+TTerror*70.0;
             else
             PIAAFPMASKRAD = 6.0+TTerror*10.0;*/
    }

    if (strstr(options,"-FPMSERR ")!=NULL)
    {
        str_pos=strstr(options,"-FPMSERR ")-options;
        str_pos = str_pos + strlen("-FPMSERR ");
        i=0;
        while((options[i+str_pos]!=' ')&&(options[i+str_pos]!='\n')&&(options[i+str_pos]!='\0'))
        {
            input[i] = options[i+str_pos];
            i++;
        }
        input[i] = '\0';
        FPMASKSIZE_ERROR = 1;
        FPMASK_FACTOR = atof(input);
        printf("--------------- FPMASK_FACTOR = %f ----------------------------\n",FPMASK_FACTOR);
    }


    if (strstr(options,"-FOC ")!=NULL)
    {
        str_pos=strstr(options,"-FOC ")-options;
        str_pos = str_pos + strlen("-FOC ");
        i=0;
        while((options[i+str_pos]!=' ')&&(options[i+str_pos]!='\n')&&(options[i+str_pos]!='\0'))
        {
            input[i] = options[i+str_pos];
            i++;
        }
        input[i] = '\0';
        FOCerror = atof(input);
        if(FOCerror>epsilon)
            FOCerror_included = 1;
        printf("FOC error is %f  (lambda at edge)\n",FOCerror);
    }

    if (strstr(options,"-Znb ")!=NULL)
    {
        str_pos=strstr(options,"-Znb ")-options;
        str_pos = str_pos + strlen("-Znb ");
        i=0;
        while((options[i+str_pos]!=' ')&&(options[i+str_pos]!='\n')&&(options[i+str_pos]!='\0'))
        {
            input[i] = options[i+str_pos];
            i++;
        }
        input[i] = '\0';
        Zindex = atol(input);
        Zerror_included = 1;
        printf("Z index is %ld\n",Zindex);
    }

    if (strstr(options,"-Zerr ")!=NULL)
    {
        str_pos=strstr(options,"-Zerr ")-options;
        str_pos = str_pos + strlen("-Zerr ");
        i=0;
        while((options[i+str_pos]!=' ')&&(options[i+str_pos]!='\n')&&(options[i+str_pos]!='\0'))
        {
            input[i] = options[i+str_pos];
            i++;
        }
        input[i] = '\0';
        Zerror = atof(input);
        if(Zerror>epsilon)
            Zerror_included = 1;
        printf("Zerr error is %f rad RMS\n",Zerror);
    }


    if((AUTOPIAACMASK==1)&&(coronagraph_type==5))
    {
        PIAAFPMASKRAD = 0.0;
        if((FOCerror_included==1)||(TTerror_included==1))
        {
            ID0 = create_2Dimage_ID("psf0",size,size);
            cnt = 0;
            if(TTerror_included==1)
            {
                for(i=-nbTTpts; i<nbTTpts; i++)
                    for(j=-nbTTpts; j<nbTTpts; j++)
                    {
                        printf("TT point %ld %ld\n",i,j);
                        x = 1.0*i*TTerror/nbTTpts;
                        y = 1.0*j*TTerror/nbTTpts;
                        r=sqrt(x*x+y*y);
                        if(r<TTerror)
                        {
                            coronagraph_simulPSF(x,y,"psf1",4,"");
                            ID1 = image_ID("psf1");
                            for(ii=0; ii<size2; ii++)
                                data.image[ID0].array.F[ii] += data.image[ID1].array.F[ii];
                            delete_image_ID("psf1");
                            printf(".");
                            fflush(stdout);
                            cnt++;
                        }
                    }
                printf("\n");
            }
            else
            {
                ID = create_2Dimage_ID("corphase",size,size);
                for(i=-nbFOCpts; i<nbFOCpts; i++)
                {
                    tmp = 1.0*i*FOCerror/nbFOCpts;
                    for(ii=0; ii<size; ii++)
                        for(jj=0; jj<size; jj++)
                        {
                            x = 1.0*(ii-size/2)/trad_pix;
                            y = 1.0*(jj-size/2)/trad_pix;
                            data.image[ID].array.F[jj*size+ii] = (x*x+y*y)*tmp*2.0*PI;
                        }
                    coronagraph_simulPSF(0.0,0.0,"psf1",4,"");
                    ID1 = image_ID("psf1");
                    for(ii=0; ii<size2; ii++)
                        data.image[ID0].array.F[ii] += data.image[ID1].array.F[ii];
                    delete_image_ID("psf1");
                    printf(".");
                    fflush(stdout);
                    cnt++;
                }
                delete_image_ID("corphase");
                printf("\n");
            }
            for(ii=0; ii<size2; ii++)
                data.image[ID0].array.F[ii] /= cnt;
        }
        else
        {
            coronagraph_simulPSF(0.0,0.0,"psf0",4,"");
            ID0 = image_ID("psf0");
        }
        peak = data.image[ID0].array.F[(size/2)*size+size/2];
        ii = 0;
        printf("%g %g\n",data.image[ID0].array.F[(size/2)*size+ii],peak);
        while(data.image[ID0].array.F[(size/2)*size+ii]<2.0*peak/pow(10.0,logcontrast))
            ii ++;
        PIAAFPMASKRAD = CORONAGRAPHS_PIXSCALE*(size/2-(ii-1));
        printf("PIAAC:  MASK RADIUS = %f lambda/d  (ii = %ld)\n",PIAAFPMASKRAD,ii);
        delete_image_ID("psf0");
    }


    if((FOCerror_included==1)||(TTerror_included==1))
    {
        ID0 = create_2Dimage_ID("psf0",size,size);
        cnt = 0;
        if(TTerror_included==1)
        {
            for(i=-nbTTpts; i<nbTTpts; i++)
                for(j=-nbTTpts; j<nbTTpts; j++)
                {
                    printf("TT point %ld %ld\n",i,j);
                    x = 1.0*i*TTerror/nbTTpts;
                    y = 1.0*j*TTerror/nbTTpts;
                    r=sqrt(x*x+y*y);
                    if(r<TTerror)
                    {
                        coronagraph_simulPSF(x,y,"psf1",coronagraph_type,"");
                        ID1 = image_ID("psf1");
                        for(ii=0; ii<size2; ii++)
                            data.image[ID0].array.F[ii] += data.image[ID1].array.F[ii];
                        delete_image_ID("psf1");
                        printf(".");
                        fflush(stdout);
                        cnt++;
                    }
                }
            printf("\n");
        }
        else
        {
            ID = create_2Dimage_ID("corphase",size,size);
            for(i=-nbFOCpts; i<nbFOCpts; i++)
            {
                tmp = 1.0*i*FOCerror/nbFOCpts;
                for(ii=0; ii<size; ii++)
                    for(jj=0; jj<size; jj++)
                    {
                        x = 1.0*(ii-size/2)/trad_pix;
                        y = 1.0*(jj-size/2)/trad_pix;
                        data.image[ID].array.F[jj*size+ii] = (x*x+y*y)*tmp*2.0*PI;
                    }
                coronagraph_simulPSF(0.0,0.0,"psf1",coronagraph_type,"");
                ID1 = image_ID("psf1");
                for(ii=0; ii<size2; ii++)
                    data.image[ID0].array.F[ii] += data.image[ID1].array.F[ii];
                delete_image_ID("psf1");
                printf(".");
                fflush(stdout);
                cnt++;
            }
            delete_image_ID("corphase");
            printf("\n");
        }
        for(ii=0; ii<size2; ii++)
            data.image[ID0].array.F[ii] /= cnt;
    }
    else
    {
        if(Zerror_included==1)
        {
            mk_zer("corphase",size,Zindex,trad_pix+1);
            ID = image_ID("corphase");
            for(ii=0; ii<size*size; ii++)
                data.image[ID].array.F[ii] *= Zerror;
            coronagraph_simulPSF(0.0,0.0,"psf0",coronagraph_type,"");
            ID0 = image_ID("psf0");
            save_fl_fits("corphase","!corphase");
            delete_image_ID("corphase");
        }
        else
        {
            coronagraph_simulPSF(0.0,0.0,"psf0",coronagraph_type,"");
            ID0 = image_ID("psf0");
        }
    }

    save_fl_fits("psf0","!psf0");
    //    exit(0);


    arrayS = (double*) malloc(sizeof(double)*size*size);
    arrayS1 = (double*) malloc(sizeof(double)*size*size);
    arrayS2 = (double*) malloc(sizeof(double)*size*size);

    if((fp=fopen(fname,"w"))==NULL)
    {
        printf("ERROR: cannot create file \"%s\"\n",fname);
        exit(0);
    }
    sep = 0.0;
    fclose(fp);


    if(coronagraph_type==30)
    {
        MaxAngle = 2.0*PI;
        NBanglestep = 64;
    }
    COROTMP1 = 0.0;
    total2old = 0.0;
    if(RADAVER==0)
    {
        iwaOK = 0;
        sepold = 0.0;
        total2old = 0.0;
        while(sep<sepmax)
        {
            if((coronagraph_type==3)||(coronagraph_type==17)||(coronagraph_type==19)||(coronagraph_type==20)||(coronagraph_type==22)||(coronagraph_type==23)||(coronagraph_type==30))
            {
                if(coronagraph_type==30)
                    coronagraph_simulPSF(-sep/sqrt(2.0),-sep/sqrt(2.0),"psf",coronagraph_type,"");
                else
                    coronagraph_simulPSF(sep/sqrt(2.0),sep/sqrt(2.0),"psf",coronagraph_type,"");
            }
            else
                coronagraph_simulPSF(sep,0.0,"psf",coronagraph_type,"");



            ID = image_ID("psf");

            total = 0.0;
            peak = 0.0;
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    flux0 =  data.image[ID0].array.F[jj*size+ii]*pow(10.0,logcontrast);
                    flux =  data.image[ID].array.F[jj*size+ii];
                    if(flux>peak)
                        peak = flux;
                    if(flux>flux0)
                        total += data.image[ID].array.F[jj*size+ii];
                }
            kk = 0;
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    flux =  data.image[ID].array.F[jj*size+ii];
                    if(flux>0.000001*peak)
                    {
                        arrayS1[kk] = data.image[ID0].array.F[jj*size+ii]*pow(10.0,logcontrast);
                        arrayS2[kk] = data.image[ID].array.F[jj*size+ii];
                        arrayS[kk] = arrayS1[kk]/arrayS2[kk];
                        kk ++;
                    }
                }
            quick_sort3(arrayS,arrayS1,arrayS2,kk);
            total1 = 0.0;
            total2 = 0.0;
            ii = 0;
            total1 += arrayS1[ii];
            total2 += arrayS2[ii];
            ii++;
            while((total1<total2)&&(ii<kk))
            {
                total1 += arrayS1[ii];
                total2 += arrayS2[ii];
                ii++;
            }

            COROTMP1 += sqrt(total2);
            if((fp=fopen(fname,"a"))==NULL)
            {
                printf("ERROR: cannot create file \"%s\"\n",fname);
                exit(0);
            }
            fprintf(fp,"%g %g %g\n",sep,total2,total);
            fclose(fp);
            printf("%g %g %g\n",sep,total2,total);
            delete_image_ID("psf");

            //IWA measurement
            if((total2>0.5)&&(iwaOK==0))
            {
                iwa = sep-(sep-sepold)*(total2-0.5)/(total2-total2old);
                iwaOK = 1;
                fp = fopen(fname,"a");
                fprintf(fp,"#IWA50 %f %f %f %f\n",TTerror,logcontrast,APLC_FPMASKsize,iwa);
                fclose(fp);
                //	      sep += 100.0;
            }

            total2old = total2;
            sepold = sep;

            sep += sepstep;
            sep *= 1.00;

            // to speed things up
            // if(total2<0.01)
            // sep += sepstep;
            //if(total2>0.9)
            // sep += sepstep;
            //	  if(sep>2.0)
            //  sepstep = 0.1;
        }
    }
    else
    {
        while(sep<sepmax)
        {
            maxtransm = 0.0;

            total3 = 0.0;
            total4 = 0.0;
            for(anglestep=0; anglestep<NBanglestep; anglestep++)
            {
                angle = MaxAngle*anglestep/NBanglestep;
                coronagraph_simulPSF(sep*cos(angle),sep*sin(angle),"psf",coronagraph_type,"");
                ID = image_ID("psf");

                total = 0.0;
                peak = 0.0;
                for(ii=0; ii<size; ii++)
                    for(jj=0; jj<size; jj++)
                    {
                        flux0 =  data.image[ID0].array.F[jj*size+ii]*pow(10.0,logcontrast);
                        flux =  data.image[ID].array.F[jj*size+ii];
                        if(flux>peak)
                            peak = flux;
                        if(flux>flux0)
                            total += data.image[ID].array.F[jj*size+ii];
                    }
                kk = 0;
                for(ii=0; ii<size; ii++)
                    for(jj=0; jj<size; jj++)
                    {
                        flux =  data.image[ID].array.F[jj*size+ii];
                        if(flux>0.000001*peak)
                        {
                            arrayS1[kk] = data.image[ID0].array.F[jj*size+ii]*pow(10.0,logcontrast);
                            arrayS2[kk] = data.image[ID].array.F[jj*size+ii];
                            arrayS[kk] = arrayS1[kk]/arrayS2[kk];
                            kk ++;
                        }
                    }
                quick_sort3(arrayS,arrayS1,arrayS2,kk);
                total1 = 0.0;
                total2 = 0.0;
                ii = 0;
                total1 += arrayS1[ii];
                total2 += arrayS2[ii];
                ii++;
                while((total1<total2)&&(ii<kk))
                {
                    total1 += arrayS1[ii];
                    total2 += arrayS2[ii];
                    ii++;
                }
                total3 += total2;
                total4 += total;
                /*      fprintf(fp,"%g %g %g\n",sep,total2,total);*/
                printf("%g (%ld/%ld) %g %g\n",sep,anglestep,NBanglestep,total2,total);
                delete_image_ID("psf");
                if(total2>maxtransm)
                    maxtransm = total2;
            }
            if((fp=fopen(fname,"a"))==NULL)
            {
                printf("ERROR: cannot create file \"%s\"\n",fname);
                exit(0);
            }
            fprintf(fp,"%g %g %g %g\n",sep,total3/NBanglestep,total4/NBanglestep,maxtransm);
            fclose(fp);
            printf("%g %g %g\n",sep,total3/NBanglestep,total4/NBanglestep);
            sep += sepstep;
            sep *= 1.00;

            //if(sep>2.0)
            // sepstep = 0.1;
        }
    }

    delete_image_ID("psf0");


    free(arrayS);
    free(arrayS1);
    free(arrayS2);

    return(0);
}


int coronagraph_userfunc()
{
    FILE *fp;
    long size = CORONAGRAPHS_ARRAYSIZE;
    double v1,x,y,r;
    long ii,jj;
    double trad_pix;
    long ID;

    double *zarray;
    long *zindex;
    long NBzern = 1;
    double ratio = 10.0;
    long k;

    fprintf(stderr,"CORONAGRAPH USER FUNCTION\n");

    zarray = (double*) malloc(sizeof(double)*NBzern);
    zindex = (long*) malloc(sizeof(long)*NBzern);
    zindex[0] = 12; //24;
    zarray[0] = 1.0e-5;
    coronagraph_PIAAperturbation(zarray, zindex, NBzern, ratio);
    free(zarray);
    free(zindex);
    exit(0);

    trad_pix = CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/2.0;
    ID = make_subpixdisk("coramp",size,size,0.5*size,0.5*size,trad_pix);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = 1.0*(ii-size/2);
            y = 1.0*(jj-size/2);
            r = sqrt(x*x+y*y);
            data.image[ID].array.F[jj*size+ii] *= fabs(x)/trad_pix;
        }
    save_fl_fits("coramp","!coramp");


    ID = make_subpixdisk("corpha",size,size,0.5*size,0.5*size,trad_pix);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = 1.0*(ii-size/2);
            y = 1.0*(jj-size/2);
            r = sqrt(x*x+y*y);
            if(x>0)
                data.image[ID].array.F[jj*size+ii] = 0.0;
            else
                data.image[ID].array.F[jj*size+ii] = PI;
        }
    save_fl_fits("corpha","!corpha");

    coronagraph_simulPSF(0.0,0.0,"psf",25,"");
    printf("Total = %g\n",arith_image_total("psf"));
    save_fl_fits("psf","!psf");

    while(0)
    {
        fprintf(stderr,"TESTING NEW CONFIGURATION\n");
        FPM_TRANSM1 = 0.0;
        FPM_TRANSM2 = 1.0; //ran1();
        FPMASK_FACTOR1 = 1.0+0.01*(0.5-ran1());
        FPMASK_FACTOR2 = 1.0+0.1*ran1();

        coronagraph_transm("transm.dat",35,10.0,"-FPMSERR 1.0");
        delete_image_ID("aplcfpm");
        v1 = COROTMP1;

        fp = fopen("result3.log","a");
        fprintf(fp,"%f %f\n",FPMASK_FACTOR1,v1);
        fclose(fp);

    }
    return(0);
}



int coronagraph_compute_limitcoeff()
{
    FILE *fp;
    long Csize=200;
    long i,k;
    double fact;
    double x,y,r;
    long ii,jj;
    long cnt;
    long NBcoeff=50;
    double *coeffmax;
    double *coeff;
    double *bestcoeff;
    double *xpow;
    double *xval;
    double *pval;
    double value;
    long NBpts = 100; /* number of evaluation points */
    double xmax = 5.0; /* maximum separation in l/d */
    double maxpval;
    int OK;
    double *pvalmax;
    long i0;
    double amp;
    double xeval = 0.5;
    double tmp;

    coeffmax = (double*) malloc(sizeof(double)*NBcoeff);
    coeff = (double*) malloc(sizeof(double)*NBcoeff);
    bestcoeff = (double*) malloc(sizeof(double)*NBcoeff);

    xval = (double*) malloc(sizeof(double)*NBpts);
    pval = (double*) malloc(sizeof(double)*NBpts);
    xpow = (double*) malloc(sizeof(double)*NBcoeff*NBpts);

    pvalmax = (double*) malloc(sizeof(double)*NBpts);

    /* initialization */

    fact = 1.0;
    for(k=0; k<NBcoeff; k++)
    {
        if(k!=0)
            fact *= k;
        value = 0.0;
        cnt = 0;
        for(ii=0; ii<Csize; ii++)
            for(jj=0; jj<Csize; jj++)
            {
                x = 1.0*ii/Csize;
                y = 1.0*jj/Csize;
                r = sqrt(x*x+y*y);
                if(r<1.0)
                {
                    value += pow(x,2*k);
                    cnt++;
                }
            }
        value /= cnt;
        coeffmax[k] = sqrt(value)*pow(PI,k)/fact;
        /*      printf("Coeffmax %ld = %g\n",k,coeffmax[k]);*/
    }



    for(i=0; i<NBpts; i++)
    {
        xval[i] = 1.0*i*xmax/NBpts;
        pvalmax[i] = 0.0;
    }

    for(i=0; i<NBpts; i++)
        for(k=0; k<NBcoeff; k++)
        {
            xpow[i*NBcoeff+k] = coeffmax[k]*pow(xval[i],k);
        }

    for(i=0; i<NBpts; i++)
    {
        tmp = 0.0;
        for(k=1; k<NBcoeff; k++)
            tmp += xpow[i*NBcoeff+k];
        printf("%f %g\n",xval[i],tmp);
    }
    exit(0);

    for(k=0; k<NBcoeff; k++)
        bestcoeff[k] = 0.0;


    /* initialize coeffs */
    while(1==1)
    {
        amp = pow(ran1(),5.0);
        for(k=0; k<NBcoeff; k++)
            coeff[k] = bestcoeff[k]+amp*(ran1()-0.5);


        for(k=0; k<NBcoeff; k++)
        {
            if(coeff[k]>1.0)
                coeff[k] = 1.0;
            if(coeff[k]<-1.0)
                coeff[k] = -1.0;
        }

        coeff[0] = 0.0;


        for(i=0; i<NBpts; i++)
        {
            pval[i] = 0.0;
            for(k=0; k<NBcoeff; k++)
                pval[i] += xpow[i*NBcoeff+k]*coeff[k];
        }
        maxpval = 0.0;
        for(i=0; i<NBpts; i++)
            if(fabs(pval[i])>fabs(maxpval))
                maxpval = pval[i];
        for(i=0; i<NBpts; i++)
            pval[i] /= maxpval;

        OK = 0;
        /*      for(i=0;i<NBpts;i++)
        {
          if(pval[i]>pvalmax[i])
            {
              pvalmax[i] = pval[i];
              OK = 1;
            }
            }*/

        i0 = (long) (xeval/xmax*NBpts);
        if(pval[i0]>pvalmax[i0])
        {
            for(i=0; i<NBpts; i++)
                pvalmax[i] = pval[i];

            for(k=0; k<NBcoeff; k++)
                bestcoeff[k] = coeff[k];

            OK = 1;
        }

        if(OK==1)
        {
            printf("%g\n",pvalmax[i0]*pvalmax[i0]);
            fflush(stdout);
            fp=fopen("transm.txt","w");
            for(i=0; i<NBpts; i++)
                fprintf(fp,"%f %f\n",xval[i],pvalmax[i]*pvalmax[i]);
            fclose(fp);

            fp=fopen("bestcoeff.txt","w");
            for(k=0; k<NBcoeff; k++)
                fprintf(fp,"%ld %f\n",k,fabs(coeff[k]));
            fclose(fp);
        }
        /*      else
        {
          printf(".");
          fflush(stdout);
          }*/
    }

    free(coeffmax);
    free(coeff);
    free(xpow);
    free(bestcoeff);
    free(xval);
    free(pval);
    free(pvalmax);

    return(0);
}




int CORONAGRAPHS_scanPIAACMC_centObs_perf( double obs0input )
{
    FILE *fpconf;
    FILE *fp;
    FILE *fpresult;
    int result = 0;

    char command[500];

    double obs0;
    double obs1; // central obstruction in input and output
    double fpmaskrad; // focal plane mask radius

    long kmax = 18;
    double *separray;
    double *transmarray;
    long k;

    char fname[500];
    char fname1[500];
    char fnamec[500];
    char cmdline[500];
    char fnameprol[500];
    char fnameprolinfo[500];
    char profname[500];
    char transmname[500];
    char fnamenew[500];

    struct stat buffer;

    long size = 256;
    long ii, jj;
    double x, y, r;
    int CentralObstructionFlag = 0;

    double pixscale = 10.25; // pix per l/D
    double peak;

    double starprof_step = 0.1; // l/D
    double starprof_maxrad = 10.0; // l/D
    long starprofNBstep;
    double *starprofrad;
    double *starprofarray;
    double *starprofarraycnt;
    long i;

    double star_radius = 0.01; // l/D
    long ID;

    int COMPUTEDATA = 1;


    float *obs1array = NULL;
    long obs1arraysize = 1;
    long obs1_index = 0;
    float obs1min = 0.0;
    float obs1max = 0.4;
    float obs1step = 0.01;

    float *fpmaskradarray = NULL;
    long fpmaskradarraysize = 1;
    long fpmaskrad_index = 0;
    float fpmaskradmin = 0.8;
    float fpmaskradmax = 3.0;
    float fpmaskradstep = 0.1;

    float tmpf0, tmpf1, tmpf2, tmpf3;
    long IDtransm, IDiwa, IDcontrastmax, IDfpmasktransm;
    float alpha, iwa;
    float tmax = 1.0;

    double eval_contrast_start = 2.5;
    double eval_contrast_end = 2.8;
    double maxval;
    double fpmasktransm;

    long OBS1INDEXSTART = 0;

    int psfsok = 1;
    int skipPSFcomp = 0;

    int ret;


    printf("RUNNING PIAACMC DESIGN SCAN FOR CENTRAL OBSTRUCTION %f\n", obs0input  );
    fflush(stdout);



    if(1) // FAST MODE, 10x faster
    {
        obs1step = 0.05;
        fpmaskradstep = 0.1;
    }


    starprofNBstep = (long) (starprof_maxrad/starprof_step);
    starprofrad =  (double*) malloc(sizeof(double)*starprofNBstep);
    starprofarray = (double*) malloc(sizeof(double)*starprofNBstep);
    starprofarraycnt = (double*) malloc(sizeof(double)*starprofNBstep);

    transmarray = (double*) malloc(sizeof(double)*kmax);

    separray = (double*) malloc(sizeof(double)*kmax);
    separray[0] = 0.0;
    separray[1] = 0.005;
    separray[2] = 0.01;
    separray[3] = 0.02;
    separray[4] = 0.1;
    separray[5] = 0.2;
    separray[6] = 0.5;
    separray[7] = 1.0;
    separray[8] = 1.5;
    separray[9] = 2.0;
    separray[10] = 2.5;
    separray[11] = 3.0;
    separray[12] = 3.5;
    separray[13] = 4.0;
    separray[14] = 5.0;
    separray[15] = 6.0;
    separray[16] = 7.0;
    separray[17] = 8.0;

    // precompute arrays size
    obs1arraysize = (long) ((obs1max-obs1min)/obs1step);
    obs1array = (float*) malloc(sizeof(float)*obs1arraysize);
    obs1_index = 0;
    for(obs1_index=0; obs1_index<obs1arraysize; obs1_index++)
        obs1array[obs1_index] = obs1min + obs1step*obs1_index;

    fpmaskradarraysize = (long) ((fpmaskradmax-fpmaskradmin)/fpmaskradstep);
    fpmaskradarray = (float*) malloc(sizeof(float)*fpmaskradarraysize);
    fpmaskrad_index = 0;
    for(fpmaskrad_index=0; fpmaskrad_index<fpmaskradarraysize; fpmaskrad_index++)
        fpmaskradarray[fpmaskrad_index] = fpmaskradmin + fpmaskradstep*fpmaskrad_index;




    if((fp = fopen("obs1indexstart.txt", "r"))!=NULL)
    {
        ret = fscanf(fp, "%ld", &OBS1INDEXSTART);
        fclose(fp);
    }




    if(obs0input<0.0)
    {
        COMPUTEDATA = 1;
        obs0 = -obs0input;
    }
    else
        obs0 = obs0input;



    if(0) // rename
    {
        for(fpmaskrad_index=0; fpmaskrad_index<fpmaskradarraysize; fpmaskrad_index++)
            for(obs1_index=OBS1INDEXSTART; obs1_index<obs1arraysize; obs1_index++)
            {
                obs1 = obs1array[obs1_index];
                fpmaskrad = fpmaskradarray[fpmaskrad_index];

                for(k=0; k< kmax; k++)
                {
                    printf("=============================================================================\n");
                    sprintf(command, "mv ./PSF/PSF_%4.2f_%5.3f_%5.3f_%5.3f.fits  ./PSF/PSF_%5.3f_%5.3f_%5.3f_%5.3f.fits\n", fpmaskrad, obs0, obs1, separray[k], fpmaskrad, obs0, obs1, separray[k]);
                    printf("%s", command);
                    ret = system(command);

                    sprintf(command, "mv ./PSF/PSFc_%4.2f_%5.3f_%5.3f_%5.3f.fits  ./PSF/PSFc_%5.3f_%5.3f_%5.3f_%5.3f.fits\n", fpmaskrad, obs0, obs1, separray[k], fpmaskrad, obs0, obs1, separray[k]);
                    printf("%s", command);
                    ret = system(command);
                }

                sprintf(command, "mv ./PROF/star_%4.2f_%5.3f_%5.3f.prof ./PROF/star_%5.3f_%5.3f_%5.3f.prof\n", fpmaskrad, obs0, obs1, fpmaskrad, obs0, obs1);
                ret = system(command);


                //	    sprintf(command, "mv ./TRANSM/transm_%4.2f_%5.3f_%5.3f.prof ./TRANSM/transm_%5.3f_%5.3f_%5.3f.prof\n",  fpmaskrad, obs0, obs1, fpmaskrad, obs0, obs1);
                //printf(command);
                //system(command);
            }

        exit(0);
    }




    if(COMPUTEDATA==1)
    {
        for(fpmaskrad_index=0; fpmaskrad_index<fpmaskradarraysize; fpmaskrad_index++)
        {
            //      for(fpmaskrad_index=fpmaskradarraysize-1; fpmaskrad_index>0; fpmaskrad_index--)


            for(obs1_index=OBS1INDEXSTART; obs1_index<obs1arraysize; obs1_index++)
            {
                obs1 = obs1array[obs1_index];
                if(obs1 > 0.001)
                {
                    CentralObstructionFlag = 1;
                    //  printf("Central obstruction = %f\n", obs1);
                }
                fpmaskrad = fpmaskradarray[fpmaskrad_index];


                if((fpmaskrad>0.8)&&(fpmaskrad<2.0)&&(obs1>0.0)&&(obs1<0.4)) // put limits here if required
                {
                    // write configuration file
                    fpconf = fopen("APLCPIAA.conf","w");
                    fprintf(fpconf,"1 %5.3f-%5.3f-%5.3f 2 1\n", fpmaskrad, obs0, obs1);
                    printf("1 %5.3f-%5.3f-%5.3f 2 1\n", fpmaskrad, obs0, obs1);
                    fclose(fpconf);

                    if(CentralObstructionFlag == 0)
                        sprintf(fnameprol, "%s/APLCapo/raw/APLCapo_%5.3f.4096.ref.gz", CORONAGRAPHSDATALOCAL, fpmaskrad);
                    else
                        sprintf(fnameprol, "%s/APLCapo/raw/APLCapo_%5.3f.%5.3f.4096.ref.gz", CORONAGRAPHSDATALOCAL, fpmaskrad, obs1);

                    sprintf(fnamenew,  "%s/APLCapo/APLCapo_%5.3f.%5.3f.4096.ref.new1", CORONAGRAPHSDATALOCAL, fpmaskrad, obs1);

                    if(stat(fnameprol,&buffer)==0) // if the prolate file exists
                    {
                        printf("FILE %s EXISTS... MOVING AHEAD\n", fnameprol);
                        psfsok = 1;
                        for(k=0; k< kmax; k++)
                        {
                            if(stat("STOP",&buffer)==0)
                            {
                                result = system("rm STOP");
                                exit(0);
                            }
                            sprintf(fname, "./PSF/PSF_%5.3f_%5.3f_%5.3f_%5.3f.fits", fpmaskrad, obs0, obs1, separray[k]);
                            sprintf(fnamec, "./PSF/PSFc_%5.3f_%5.3f_%5.3f_%5.3f.fits", fpmaskrad, obs0, obs1, separray[k]);

                            if(stat(fname,&buffer)!=0) // if PSF file does not exist, create it
                            {
                                if(skipPSFcomp==1)
                                    psfsok = 0;
                                else
                                {
                                    printf("============= CREATING %s ===========\n", fname);
                                    fpconf = fopen("APLCPIAA.conf","w");
                                    fprintf(fpconf,"1 %5.3f-%5.3f-%5.3f 2 1\n", fpmaskrad, obs0, obs1);
                                    printf("1 %5.3f-%5.3f-%5.3f 2 1\n", fpmaskrad, obs0, obs1);
                                    fclose(fpconf);

                                    sprintf(cmdline, "./mkpsf %f %f", separray[k], 0.0);
                                    result = system(cmdline);
                                    sprintf(cmdline, "mv psf.fits %s", fname);
                                    result = system(cmdline);
                                }
                            }

                            if(psfsok==1)
                            {
                                sprintf(fname1, "psf%02ld", k);
                                load_fits(fname, fname1, 1);
                            }
                        }


                        if(psfsok==1)
                        {
                            for(i=0; i<starprofNBstep; i++)
                            {
                                starprofarraycnt[i] = 0.0;
                                starprofarray[i] = 0.0;
                            }

                            k = 2; // 0.01 l/D
                            sprintf(fname1, "psf%02ld", k);
                            ID = image_ID(fname1);
                            for(ii=0; ii<size; ii++)
                                for(jj=0; jj<size; jj++)
                                {
                                    x = 1.0*ii-size/2;
                                    y = 1.0*jj-size/2;
                                    x /= pixscale;
                                    y /= pixscale;
                                    r = sqrt(x*x+y*y);
                                    i = (long) (r/starprof_step+0.5);
                                    if(i<starprofNBstep)
                                    {
                                        starprofarray[i] += data.image[ID].array.F[jj*size+ii];
                                        starprofarraycnt[i] += 1.0;
                                    }
                                }



                            sprintf(fname1, "psf%02d", (int) (kmax-1));
                            peak = arith_image_max(fname1);
                            printf("PEAK = %f\n", peak);

                            sprintf(profname, "./PROF/star_%5.3f_%5.3f_%5.3f.prof",  fpmaskrad, obs0, obs1);
                            fp = fopen(profname, "w");
                            for(i=0; i<starprofNBstep; i++)
                            {
                                starprofarray[i] /= (starprofarraycnt[i]+0.000001);
                                starprofarray[i] *= 0.5 * pow((star_radius/separray[k]),2.0)/peak;
                                fprintf(fp, "%5.3f %g\n", starprof_step*i, starprofarray[i] );
                            }
                            fclose(fp);


                            // compute transmission curve if it does not exist
                            sprintf(transmname, "./TRANSM/transm_%5.3f_%5.3f_%5.3f.prof",  fpmaskrad, obs0, obs1);
                            if(stat(transmname, &buffer)!=0) // if does not exist
                            {
                                tmax = 0.0;
                                for(k=0; k< kmax; k++)
                                {
                                    sprintf(fname1, "psf%02ld", k);
                                    transmarray[k] = arith_image_total(fname1)/(1.0-obs0*obs0);
                                    if(transmarray[k]>tmax)
                                        tmax = transmarray[k];
                                }

                                for(k=0; k< kmax; k++)
                                    transmarray[k] /= tmax;

                                fp = fopen(transmname, "w");
                                for(k=0; k< kmax; k++)
                                {
                                    sprintf(fname1, "psf%02ld", k);
                                    fprintf(fp,"%5.3f %15.13f\n", separray[k], transmarray[k]);
                                }
                                fclose(fp);
                            }
                            for(k=0; k< kmax; k++)
                            {
                                sprintf(fname1, "psf%02ld", k);
                                delete_image_ID(fname1);
                            }
                        }
                    }
                    else
                        printf("--------- FILE %s MISSING -> NO PROCESSING ----\n", fnameprol);
                }
            }
        }
    }


    // ANALYSIS
    //



    // create transmission cube and IWA image
    //
    fpresult = fopen("result.log", "w");
    IDtransm = create_3Dimage_ID("PIAACMCtransm", obs1arraysize, fpmaskradarraysize, kmax);
    IDiwa = create_2Dimage_ID("PIAACMCiwa", obs1arraysize, fpmaskradarraysize);
    IDcontrastmax = create_2Dimage_ID("PIAACMCcontrastmax", obs1arraysize, fpmaskradarraysize);
    IDfpmasktransm = create_2Dimage_ID("PIAACMCfpmasktransm", obs1arraysize, fpmaskradarraysize);

    for(obs1_index=0; obs1_index<obs1arraysize; obs1_index++)
    {
        obs1 = obs1array[obs1_index];
        if(obs1 > 0.001)
        {
            CentralObstructionFlag = 1;
            printf("Central obstruction = %f\n", obs1);
        }
        for(fpmaskrad_index=0; fpmaskrad_index<fpmaskradarraysize; fpmaskrad_index++)
        {
            //for(fpmaskrad_index=fpmaskradarraysize-1; fpmaskrad_index>0; fpmaskrad_index--)

            fpmaskrad = fpmaskradarray[fpmaskrad_index];

            // load info file for apodization
            if(CentralObstructionFlag == 0)
                sprintf(fnameprolinfo, "%s/APLCapo/APLCapo_%5.3f.4096.ref.info", CORONAGRAPHSDATALOCAL, fpmaskrad);
            else
                sprintf(fnameprolinfo, "%s/APLCapo/APLCapo_%5.3f.%5.3f.4096.ref.info", CORONAGRAPHSDATALOCAL, fpmaskrad, obs1);
            if(stat(fnameprolinfo, &buffer)==0) // if file exists
            {
                fp = fopen(fnameprolinfo, "r");
                ret = fscanf(fp, "%f %f %f %f\n", &tmpf0, &tmpf1, &tmpf2, &tmpf3);
                fclose(fp);
                fpmasktransm =  -(1.0-tmpf3)/tmpf3;
                data.image[IDfpmasktransm].array.F[fpmaskrad_index*obs1arraysize + obs1_index] = fpmasktransm;
            }

            // load and read radial profile for star
            maxval = 0.0;
            sprintf(profname, "./PROF/star_%5.3f_%5.3f_%5.3f.prof",  fpmaskrad, obs0, obs1);
            if(stat(profname, &buffer)==0) // if file exists
            {
                printf("Reading file %s\n", profname);
                fp = fopen(profname, "r");
                for(i=0; i<starprofNBstep; i++)
                {
                    ret = fscanf(fp, "%f %f\n", &tmpf0, &tmpf1 );
                    //		  printf(" -- %5.3f %g\n", tmpf0, tmpf1 );
                    if((tmpf0>eval_contrast_start)&&(tmpf0<eval_contrast_end))
                        if(tmpf1>maxval)
                            maxval = tmpf1;
                }
                fclose(fp);
            }
            else
            {

            }


            printf("== %s ==   %f  %f  %f          %e   [%ld %ld]\n", profname, fpmaskrad, obs0, obs1, maxval, fpmaskrad_index, obs1_index);
            //if(obs1>0.181)
            //exit(0);

            data.image[IDcontrastmax].array.F[fpmaskrad_index*obs1arraysize + obs1_index] = maxval;


            // load and read transmission
            sprintf(transmname, "./TRANSM/transm_%5.3f_%5.3f_%5.3f.prof",  fpmaskrad, obs0, obs1);
            if(stat(transmname, &buffer)==0) // if file exists
            {
                printf("Reading file %s\n", transmname);
                fp = fopen(transmname, "r");
                for(k=0; k< kmax; k++)
                {
                    result = fscanf(fp,"%f %f\n", &tmpf0, &tmpf1);
                    data.image[IDtransm].array.F[k*fpmaskradarraysize*obs1arraysize + fpmaskrad_index*obs1arraysize + obs1_index] = tmpf1;
                    transmarray[k] = tmpf1;
                }
                fclose(fp);

                k = 1;
                tmpf0 = transmarray[k-1];
                tmpf1 = transmarray[k];
                while((tmpf1<0.5)&&(k<kmax-1))
                {
                    tmpf0 = transmarray[k-1];
                    tmpf1 = transmarray[k];
                    k++;
                }
                k--;
                alpha = (0.5-transmarray[k-1])/(transmarray[k]-transmarray[k-1]);
                iwa =  separray[k-1] + alpha*(separray[k]-separray[k-1]);
                data.image[IDiwa].array.F[fpmaskrad_index*obs1arraysize + obs1_index] = iwa;
                printf("%ld %f %f   %ld %f %f  -> [%f %f] %f %f\n", k-1, separray[k-1], tmpf0, k, separray[k], tmpf1, (0.5-transmarray[k-1]), (transmarray[k]-transmarray[k-1]), alpha, iwa);
            }
            else
            {
                //    printf("Skipping file %s\n", transmname);
                iwa = 0.0;
            }

            if(iwa>0.5) // point is good
                fprintf(fpresult, "%5.3f %5.3f %g %g %g\n", obs1, fpmaskrad, iwa, maxval, fpmasktransm);

        }
    }
    fclose(fpresult);

    sprintf(fname, "!PIAACMCtransm_%5.3f.fits", obs0input);
    save_fl_fits("PIAACMCtransm", fname);

    sprintf(fname, "!PIAACMCiwa_%5.3f.fits", obs0input);
    save_fl_fits("PIAACMCiwa", fname);

    sprintf(fname, "!PIAACMCcontrastmax_%5.3f.fits", obs0input);
    save_fl_fits("PIAACMCcontrastmax", fname);

    sprintf(fname, "!PIAACMCfpmasktransm_%5.3f.fits", obs0input);
    save_fl_fits("PIAACMCfpmasktransm", fname);


    free(separray);
    free(transmarray);
    free(starprofarray);
    free(starprofrad);

    free(obs1array);
    free(fpmaskradarray);

    return(0);
}

