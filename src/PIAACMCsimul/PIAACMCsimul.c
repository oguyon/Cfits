#include <fitsio.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <assert.h>


#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "info/info.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate/WFpropagate.h"
#include "statistic/statistic.h"
#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "image_filter/image_filter.h"
#include "image_basic/image_basic.h"
#include "coronagraphs/coronagraphs.h"
#include "PIAACMCsimul/PIAACMCsimul.h"
#include "OptSystProp/OptSystProp.h"


# ifdef HAVE_LIBGOMP
#include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
#endif



/**
 * @file PIAACMCsimul.c
 * @author Olivier Guyon
 */



int WRITE_OK = 1;

extern DATA data;

#define SBUFFERSIZE 2000
 
///  Current configuration directory
char piaacmcconfdir[300];

OPTSYST *optsyst;
int optsystinit = 0;
long IDx, IDy, IDr, IDPA;
long ID_CPAfreq;

double FPMSCALEFACTOR = 0.9; // undersize mask in array to avoid edge clipping


// this makes 20% bandwidth, from 0.55/1.1 to 0.55*1.1
double LAMBDASTART = 0.5e-6;
double LAMBDAEND = 0.605e-6;
#define NBLAMBDA 5

OPTPIAACMCDESIGN *piaacmc;


int FORCE_CREATE_Cmodes = 0;
int CREATE_Cmodes = 0;
int FORCE_CREATE_Fmodes = 0;
int CREATE_Fmodes = 0;

int FORCE_CREATE_fpmzmap = 0;
int CREATE_fpmzmap = 0;
int FORCE_CREATE_fpmzt = 0;
int CREATE_fpmzt = 0;

int FORCE_CREATE_fpmza = 0;
int CREATE_fpmza;

int FORCE_MAKE_PIAA0shape = 0;
int MAKE_PIAA0shape = 0;
int FORCE_MAKE_PIAA1shape = 0;
int MAKE_PIAA1shape = 0;

int focmMode = -1; // if != -1, compute only impulse response to corresponding zone
int PIAACMC_FPMsectors = 0; // 1 if focal plane mask should have sectors


// declared here for speed
double evalval;
long evali;
long evalk, evalki, evalki1, evalmz, evalii, evalii1, evalii2, evalkv;
double evalcosp, evalsinp, evalre, evalim, evalre1, evalim1, evalpha;
double evalv1;
double PIAACMCSIMUL_VAL;
double PIAACMCSIMUL_VAL0;
double PIAACMCSIMUL_VALREF;

// for minimization
double *fpmresp_array;
double *zonez_array;
double *zonez0_array;
double *zonez1_array;
double *zonezbest_array;
double *dphadz_array;
double *outtmp_array;
long NBoptVar;
static long LOOPCNT = 0;
long vsize;
double cval0;

double CnormFactor = 1.0; // for contrast normalization
double THICKRANGE = 2.0e-6;

int computePSF_FAST_FPMresp = 0;
int computePSF_ResolvedTarget = 0; // source size = 1e-{0.1*computePSF_ResolvedTarget}
int computePSF_ResolvedTarget_mode = 0; // 0: source is simulated as 3 points, 1: source is simulated as 6 points
int PIAACMC_FPM_FASTDERIVATIVES = 0;

 
long NBsubPix = 8;

double SCORINGTOTAL = 1.0;
double MODampl = 1.0e-6;

int SCORINGMASKTYPE = 0;

int PIAACMC_save = 1;



float PIAACMC_MASKRADLD = 0.0; // not initialized yet
float PIAACMC_MASKregcoeff = 1.0;
int PIAACMC_fpmtype = 0; // 0 for idealized PIAACMC focal plane mask, 1 for physical focal plane mask
long PIAACMC_FPMresp_mp;
long PIAACMC_FPMresp_thread;


long PIAACMC_MAXRINGOPTNB = 100; // maximum number of rings to optimize, from inside out
long PIAACMC_RINGOPTNB;

int PIAACMC_CIRC = 0; // 1 if PIAA optics must be circular symmetric



// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string, not existing image
// 4: existing image
// 5: string 
//



int PIAACMCsimul_rings2sectors_cli()
{
    if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,3)==0)
    {
        PIAACMCsimul_rings2sectors(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string);
        return 0;
    }
    else
        return 1;
}




int PIAACMCsimul_run_cli()
{

    if(CLI_checkarg(1,3)+CLI_checkarg(2,2)==0)
    {
        PIAACMCsimul_run(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl);
        return 0;
    }
    else
        return 1;
}




int PIAACMC_FPMresp_rmzones_cli()
{
    if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,2)==0)
    {
        PIAACMC_FPMresp_rmzones(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl);
        return 0;
    }
    else
        return 1;
}




int PIAACMC_FPMresp_resample_cli()
{
    if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,2)+CLI_checkarg(4,2)==0)
    {
        PIAACMC_FPMresp_resample(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl);
        return 0;
    }
    else
        return 1;
}





/**
 * Initializes module
 */
int init_PIAACMCsimul()
{
    strcpy(data.module[data.NBmodule].name, __FILE__);
    strcpy(data.module[data.NBmodule].info, "PIAACMC system simulation");
    data.NBmodule++;

   strcpy(data.cmd[data.NBcmd].key,"piaacmcsimring2sect");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = PIAACMCsimul_rings2sectors_cli;
    strcpy(data.cmd[data.NBcmd].info,"turn ring fpm design into sectors");
    strcpy(data.cmd[data.NBcmd].syntax,"<input ring fpm> <zone-ring table> <output sector fpm>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmcsimring2sect");
    strcpy(data.cmd[data.NBcmd].Ccall,"long PIAACMCsimul_rings2sectors(char *IDin_name, char *sectfname, char *IDout_name)");
    data.NBcmd++;


    strcpy(data.cmd[data.NBcmd].key,"piaacmcsimrun");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = PIAACMCsimul_run_cli;
    strcpy(data.cmd[data.NBcmd].info,"Simulate PIAACMC");
    strcpy(data.cmd[data.NBcmd].syntax,"<configuration index [string]> <mode[int]>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmcsimrun");
    strcpy(data.cmd[data.NBcmd].Ccall,"int PIAACMCsimul_run(char *confindex, long mode)");
    data.NBcmd++;


    strcpy(data.cmd[data.NBcmd].key,"piaacmsimfpmresprm");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = PIAACMC_FPMresp_rmzones_cli;
    strcpy(data.cmd[data.NBcmd].info,"remove zones in FPM resp matrix");
    strcpy(data.cmd[data.NBcmd].syntax,"<input FPMresp> <output FPMresp> <NBzone removed>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmsimfpmresprm FPMresp FPMrespout 125");
    strcpy(data.cmd[data.NBcmd].Ccall,"long PIAACMC_FPMresp_rmzones(char *FPMresp_in_name, char *FPMresp_out_name, long NBzones)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"piaacmsimfpmresprs");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = PIAACMC_FPMresp_resample_cli;
    strcpy(data.cmd[data.NBcmd].info,"resample FPM resp matrix");
    strcpy(data.cmd[data.NBcmd].syntax,"<input FPMresp> <output FPMresp> <NBlambda> <EvalPts step>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmsimfpmresprs FPMresp FPMrespout 10 2");
    strcpy(data.cmd[data.NBcmd].Ccall,"long PIAACMC_FPMresp_resample(char *FPMresp_in_name, char *FPMresp_out_name, long NBlambda, long PTstep)");
    data.NBcmd++;



    // add atexit functions here
    atexit(PIAACMCsimul_free);

    return 0;

}

/**
 * Frees memory for module
 */
void PIAACMCsimul_free( void )
{
    if(optsystinit ==1)
    {
        free(optsyst);
    }
}




















// ************************************************************************************
//
//           FOCAL PLANE MASK
//
// ************************************************************************************

/**
 * @param[out]  IDname  Name of output image
 */

long PIAACMCsimul_mkFPM_zonemap(char *IDname)
{
    FILE *fp;
    char fname[500];
    long NBzones;
    long ID;
    double x, y, r, PA;
    long ii, jj;
    long zi;
    long *sizearray;

    long ring;
    long *nbsector;
    long *nbsectorcumul;
    double PAf;
    double eps = 1.0e-6;

    unsigned short int zoneindex;
    long cnt, cnt1;
    long nbzonescc;

    double hexstep = 1.0;
    double hexgap = -0.0001;
    double hexsteppix;
    long ii1, jj1, ii1max, jj1max;
    double hx, hy;
    double hex_x[10000];
    double hex_y[10000];
    long hex_ring[10000];
    long hex_number[10000];
    long hcnt;
    long hindex;
    long hindexMax;
    long ID1;
	
	

    sizearray = (long*) malloc(sizeof(long)*2);
    sizearray[0] = piaacmc[0].fpmarraysize;
    sizearray[1] = piaacmc[0].fpmarraysize;
    ID = create_image_ID(IDname, 2, sizearray, USHORT, 0, 0);
    free(sizearray);

    nbsector = (long*) malloc(sizeof(long)*piaacmc[0].NBrings);
    nbsectorcumul = (long*) malloc(sizeof(long)*piaacmc[0].NBrings);


    switch (PIAACMC_FPMsectors) {

    case 0: // rings
        nbsectorcumul[0] = 1;
        for(ring=1; ring<piaacmc[0].NBrings; ring++)
            nbsectorcumul[ring] = nbsectorcumul[ring-1] + 1;
        break;

    case 1: // rings broken in sectors
        sprintf(fname, "%s/fpmsectors%d_%03ld.txt", piaacmcconfdir, PIAACMC_FPMsectors, piaacmc[0].NBrings);
        fp = fopen(fname, "w");
        NBzones = 0;
        cnt = 0;
        fprintf(fp, "0 0\n");
        nbsector[0] = 1;
        nbsectorcumul[0] = 1;
        for(ring=1; ring<piaacmc[0].NBrings; ring++)
        {

            nbsector[ring] = 2*(ring+1);
            nbsectorcumul[ring] = nbsectorcumul[ring-1]+nbsector[ring];
            for(cnt1=0; cnt1<nbsector[ring]; cnt1++)
            {
                cnt++;
                fprintf(fp, "%ld %ld\n", cnt, ring);
            }
        }

        fclose(fp);
        for(ring=0; ring<piaacmc[0].NBrings; ring++)
            printf("ring %ld : %ld %ld\n", ring, nbsector[ring], nbsectorcumul[ring]);
        break;

    case 2: // rings of hexagons
         sprintf(fname, "%s/fpmsectors%d_%03ld.txt", piaacmcconfdir, PIAACMC_FPMsectors, piaacmc[0].NBrings);
        fp = fopen(fname, "w");
        nbsector[0] = 1;
        nbsectorcumul[0] = 1;
       for(ring=1; ring<piaacmc[0].NBrings; ring++)
        {
            nbsector[ring] = 0;
            nbsectorcumul[ring] = 0;
        }
        hindex = 0;
        // hegagon side = s = ring unit
        ii1max = (long) (piaacmc[0].NBrings/3+2);
        jj1max = (long) (piaacmc[0].NBrings/sqrt(3.0)+2);
        for(ii1=-ii1max; ii1<ii1max; ii1++)
            for(jj1=-jj1max; jj1<jj1max; jj1++)
            {
                hx = hexstep*ii1*3;
                hy = hexstep*sqrt(3.0)*jj1;
                ring = (long) sqrt(hx*hx+hy*hy);
                if(ring<piaacmc[0].NBrings)
                {
                    nbsector[ring] ++;
                    hex_x[hindex] = hx;
                    hex_y[hindex] = hy;
                    hex_ring[hindex] = ring;
                    hindex++;
                }


                hx += hexstep*1.5;
                hy += hexstep*sqrt(3.0)/2.0;
                ring = (long) sqrt(hx*hx+hy*hy);
                if(ring<piaacmc[0].NBrings)
                {
                    nbsector[ring] ++;
                    hex_x[hindex] = hx;
                    hex_y[hindex] = hy;
                    hex_ring[hindex] = ring;
                    hindex++;
                }
            }
        hindexMax = hindex;
        fprintf(fp, "0 0\n");
        hcnt = 1;
        for(ring=1; ring<piaacmc[0].NBrings; ring++)
        {
            for(hindex=0; hindex<hindexMax; hindex++)
                if(hex_ring[hindex]==ring)
                {
                    hex_number[hindex] = hcnt;
                    fprintf(fp, "%ld %ld\n", hcnt, ring);
                    hcnt++;
                }
            if(ring>0)
                nbsectorcumul[ring] = nbsectorcumul[ring-1]+nbsector[ring];
        }
        fclose(fp);
         for(ring=0; ring<piaacmc[0].NBrings; ring++)
            printf("ring %ld : %ld %ld\n", ring, nbsector[ring], nbsectorcumul[ring]);
        break;

    default:
        printf("ERROR: FPMsector mode (%d) not recognized\n", PIAACMC_FPMsectors);
        exit(0);
        break;
    }



    if(piaacmc[0].NBringCentCone>0)
        nbzonescc = nbsectorcumul[piaacmc[0].NBringCentCone-1];
    else
        nbzonescc = 0;



    //  ID = create_2Dimage_ID(IDname, piaacmc[0].fpmarraysize, piaacmc[0].fpmarraysize);


    if((PIAACMC_FPMsectors==0)||(PIAACMC_FPMsectors==1))
    {
        for(ii=0; ii<piaacmc[0].fpmarraysize; ii++)
            for(jj=0; jj<piaacmc[0].fpmarraysize; jj++)
            {
                x = (2.0*ii-1.0*piaacmc[0].fpmarraysize)/piaacmc[0].fpmarraysize / FPMSCALEFACTOR;
                y = (2.0*jj-1.0*piaacmc[0].fpmarraysize)/piaacmc[0].fpmarraysize / FPMSCALEFACTOR;
                r = sqrt(x*x+y*y);
                PA = atan2(y,x);
                PAf = 0.5*((PA/M_PI)+1.0);
                if(PAf<eps)
                    PAf = eps;
                if(PAf>1.0-eps)
                    PAf = 1.0-eps;

                zi = (long) ceil((1.0-r)*piaacmc[0].NBrings);


                if(zi<0.1)
                    zi = 0;
                if(zi>piaacmc[0].NBrings)
                    zi = piaacmc[0].NBrings;

                ring = piaacmc[0].NBrings-zi; // 0 for inner disk, increases outward
                if(zi==0)
                    ring = -1;


                if(PIAACMC_FPMsectors==0)
                {
                    zoneindex = (unsigned short int) (piaacmc[0].NBrings-zi+1);
                    if(zi==0)
                        zoneindex = 0;
                }
                else
                {
                    if(ring==-1)
                        zoneindex = 0;
                    else
                    {
                        if(ring==0) // inner disk
                            zoneindex = 1;
                        else
                        {
                            zoneindex = (unsigned short int) nbsectorcumul[ring-1]+1;
                            zoneindex += (unsigned short int) (PAf*nbsector[ring]);
                        }

                        if(piaacmc[0].NBrings>1)
                        {
                            if(ring<piaacmc[0].NBringCentCone)
                                zoneindex = 0;
                            else
                                zoneindex -= nbzonescc;
                        }
                    }
                }

                data.image[ID].array.U[jj*piaacmc[0].fpmarraysize+ii] = zoneindex;

            }

        if(PIAACMC_FPMsectors==0)
        {
            if(piaacmc[0].NBrings>1)
                piaacmc[0].focmNBzone = piaacmc[0].NBrings - nbzonescc;
            else
                piaacmc[0].focmNBzone = 1;
        }
        else
        {
            if(piaacmc[0].NBrings>1)
                piaacmc[0].focmNBzone = nbsectorcumul[piaacmc[0].NBrings-1] - nbzonescc;
            else
                piaacmc[0].focmNBzone = 1;
        }
    }


    if(PIAACMC_FPMsectors==2)
        {
            hexsteppix = 0.5*piaacmc[0].fpmarraysize/piaacmc[0].NBrings * FPMSCALEFACTOR;
            for(hindex=0;hindex<hindexMax;hindex++)
            {
                ID1 = make_hexagon("_TMPhex", piaacmc[0].fpmarraysize, piaacmc[0].fpmarraysize, 0.5*piaacmc[0].fpmarraysize + hex_x[hindex]*hexsteppix, 0.5*piaacmc[0].fpmarraysize + hex_y[hindex]*hexsteppix, hexsteppix*(1.0-hexgap)*(sqrt(3.0)/2.0));
                
                for(ii=0; ii<piaacmc[0].fpmarraysize*piaacmc[0].fpmarraysize; ii++)                   
                        {
                            if(data.image[ID1].array.F[ii] > 0.5)
                                data.image[ID].array.U[ii] = (unsigned int) hex_number[hindex]+1;
                        }
                delete_image_ID("_TMPhex");
            }
            
            if(piaacmc[0].NBrings>1)
                piaacmc[0].focmNBzone = nbsectorcumul[piaacmc[0].NBrings-1];
            else
                piaacmc[0].focmNBzone = 1;
        }


    printf("[%d] piaacmc[0].focmNBzone  =  %ld %ld    %ld %ld   ->  %ld   (%ld)\n", PIAACMC_FPMsectors, piaacmc[0].NBrings, nbsectorcumul[piaacmc[0].NBrings-1], piaacmc[0].NBringCentCone, nbzonescc, piaacmc[0].focmNBzone, piaacmc[0].NBrings);


   if(PIAACMC_FPMsectors!=0)
    {
      //  printf("Saving %s ....\n", IDname);
      //  save_fits(IDname, "!__test_zonemap_00.fits"); //TEST
       // sleep(100000);
    }

    free(nbsector);
    free(nbsectorcumul);

    return ID;
}





/// @param[in] IDin_name	input image: circular mask design
/// @param[in] sectfname	text file specifying which zones belong to which rings
/// @param[out] IDout_name	output sector mask design
long PIAACMCsimul_rings2sectors(char *IDin_name, char *sectfname, char *IDout_name)
{
    long IDin, IDout;
    FILE *fp;
    long nbring, nbzone;
    long tmpl1, tmpl2;
    long zone;
    long arrayring[5000];

    IDin = image_ID(IDin_name);
    nbring = data.image[IDin].md[0].size[0];

    nbzone = 0;
    nbring = 0;
    fp = fopen(sectfname,"r");
    while(fscanf(fp,"%ld %ld\n", &tmpl1, &tmpl2)==2)
    {
        arrayring[tmpl1] = tmpl2;
        if(tmpl2>nbring)
            nbring = tmpl2;
        if(tmpl1>nbzone)
            nbzone = tmpl1;
    }
    fclose(fp);
    nbring++;
    nbzone++;

    IDout = create_2Dimage_ID_double(IDout_name, nbzone, 1);
    for(zone=0; zone<nbzone; zone++)
        data.image[IDout].array.D[zone] = data.image[IDin].array.D[arrayring[zone]];

    printf("%ld zones in %ld rings\n", nbzone, nbring);

    return(IDout);
}





///
/// @param[in]  IDzonemap_name	zones
/// @param[in]  ID_name
/// @param[in]  mode       if mode = -1, make whole 1-fpm, if mode = zone, make only 1 zone with CA = (1.0, 0.0)
/// @param[in]  saveMask   1 if mask saved to file system
//
// makes 1-fpm CA
// zone numbering starts here from 1 (zone 1 = outermost ring)
//
long PIAACMCsimul_mkFocalPlaneMask(char *IDzonemap_name, char *ID_name, int mode, int saveMask)
{
    long ID, IDm;
    long IDz;
    long IDsag, IDzone;
    long size;
    long nblambda;
    long k;
    long ii, jj;
    double x, y, r; // in meter
    long ii1, jj1;
    double fpscale; // [m/pix]
    int zi;
    double t, a, amp;
    long size2;
    double pha, re, im;
    long iii, jjj;
    double retmp, imtmp, ttmp, zonetmp;

    int CentCone = 0;
    int OuterCone = 0;
	char fname[1000];


    size = optsyst[0].size;
    size2 = size*size;
    nblambda = optsyst[0].nblambda;


    IDz = image_ID(IDzonemap_name);
    ID = create_3DCimage_ID(ID_name, size, size, nblambda);
    IDsag = create_3Dimage_ID("fpmsag", size, size, nblambda);
    IDzone = create_3Dimage_ID("fpmzone", size, size, nblambda);

/*    printf("Saving image %s\n", IDzonemap_name);
    save_fits(IDzonemap_name, "!__test_zonemap_02.fits");
    sleep(10);//TEST
  */
    
    if(piaacmc[0].NBrings>2)
    {
        CentCone = 1;
        OuterCone = 1;
    }

    printf("===================== Make focal plane mask  %s %s %d\n", IDzonemap_name, ID_name, mode);


    // CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/CORONAGRAPHS_ARRAYSIZE

    /*  printf("%ld %ld\n", piaacmc[0].zoneaID, piaacmc[0].zonezID);
    for(k=0;k<data.image[piaacmc[0].zonezID].md[0].size[0];k++)
    {
    t = data.image[piaacmc[0].zonezID].array.D[k]; // thickness
    a = data.image[piaacmc[0].zoneaID].array.D[k]; // amplitude transmission
    amp = a;
    pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, t, optsyst[0].lambdaarray[0]);

    re = amp*cos(pha);
    im = amp*sin(pha);

    printf("ZONE %2ld  %12g %12g   %12g %12g    %12g %12g\n", k, a, t, amp, pha, 1.0-re, -im);
    }
    */

    fpscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size/piaacmc[0].fpzfactor*optsyst[0].lambdaarray[0]*piaacmc[0].Fratio;
    printf("piaacmc[0].fpmRad = %g m    fpscale[0] = %g    mode = %d\n", piaacmc[0].fpmRad, fpscale, mode);

    


# ifdef HAVE_LIBGOMP
    #pragma omp parallel default(shared) private(ii, jj, x, y, r, retmp, imtmp, iii, jjj, ii1, jj1, zi, t, a, fpscale, amp, pha, ttmp, zonetmp)
    {
        #pragma omp for
# endif
        for(k=0; k<nblambda; k++)
        {
            fpscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size/piaacmc[0].fpzfactor*optsyst[0].lambdaarray[k]*piaacmc[0].Fratio;
            printf("LAMBDA %3ld / %3ld = %10.5g m    SCALE = %10.5g m/pix   size=%4ld  rad=%g\n", k, nblambda, optsyst[0].lambdaarray[k], fpscale, size, piaacmc[0].fpmRad);
            printf("Zone 0 amplitude [%ld]: %lf\n", piaacmc[0].zoneaID, data.image[piaacmc[0].zoneaID].array.D[0]);
            printf("Zone 0 thickness: %g\n", data.image[piaacmc[0].zonezID].array.D[0]);
            printf("piaacmc[0].fpmRad = %g m\n", piaacmc[0].fpmRad);
            printf("piaacmc[0].fpmCentConeRad = %g m\n", piaacmc[0].fpmCentConeRad);
            printf("piaacmc[0].fpmOuterConeRad [%d] = %g m\n",  OuterCone, piaacmc[0].fpmOuterConeRad);

            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    //printf("[ %4ld %4ld ] ", ii, jj);

                    x = (1.0*ii-size/2)*fpscale; // [m]
                    y = (1.0*jj-size/2)*fpscale; // [m]
                    r = sqrt(x*x+y*y); // [m]

                    t = 0.0;
                    a = 1.0;

                    // outer part
                    if(OuterCone==1)
                    {
                        t = 0.0; //piaacmc[0].fpmOuterConeZ;
                        a = 1.0;
                    }

                    if((r>0.9*piaacmc[0].fpmRad)&&(r<piaacmc[0].fpmOuterConeRad)) // outer cone
                    {
                        t = piaacmc[0].fpmOuterConeZ*(piaacmc[0].fpmOuterConeRad-r)/(piaacmc[0].fpmOuterConeRad-piaacmc[0].fpmRad);
                        // (r-piaacmc[0].fpmRad)/(piaacmc[0].fpmOuterConeRad-piaacmc[0].fpmRad);
                        a = 1.0;
                    }

                    data.image[IDzone].array.F[k*size2+jj*size+ii] = 0;

                    if(r<1.1*piaacmc[0].fpmOuterConeRad) //     piaacmc[0].fpmRad)
                    {
                        //       printf("pix outercone ...");
                        //      fflush(stdout);

                        retmp = 0.0;
                        imtmp = 0.0;
                        ttmp = 0.0;
                        zonetmp = 0.0;
                        for(iii=0; iii<NBsubPix; iii++)
                        {
                            for(jjj=0; jjj<NBsubPix; jjj++)
                            {
                                // physical coordinates on mask
                                // x and y in [m]
                                x = (1.0*ii - size/2 + 1.0*(0.5+iii)/NBsubPix-0.5)*fpscale;
                                y = (1.0*jj - size/2 + 1.0*(0.5+jjj)/NBsubPix-0.5)*fpscale;
                                r = sqrt(x*x+y*y); // [m]

                                zi = 0; // default
                                if((r>0.9*piaacmc[0].fpmRad)&&(r<piaacmc[0].fpmOuterConeRad)) // outer cone
                                {
                                    t = piaacmc[0].fpmOuterConeZ*(piaacmc[0].fpmOuterConeRad-r)/(piaacmc[0].fpmOuterConeRad-piaacmc[0].fpmRad);
                                  //  piaacmc[0].fpmOuterConeZ*(r-piaacmc[0].fpmRad)/(piaacmc[0].fpmOuterConeRad-piaacmc[0].fpmRad);
                                    a = 1.0;
                                }

                                ii1 = (long) ( (0.5 + 0.5*x/piaacmc[0].fpmRad*FPMSCALEFACTOR)*piaacmc[0].fpmarraysize + 0.5);
                                jj1 = (long) ( (0.5 + 0.5*y/piaacmc[0].fpmRad*FPMSCALEFACTOR)*piaacmc[0].fpmarraysize + 0.5);
                                if((ii1>-1)&&(ii1<piaacmc[0].fpmarraysize)&&(jj1>-1)&&(jj1<piaacmc[0].fpmarraysize))
                                {
                                    if(CentCone==1)
                                    {
                                        // central cone
                                        if((r<0.99*piaacmc[0].fpmRad)&&(r<piaacmc[0].fpmCentConeRad))
                                        {
                                            t = piaacmc[0].fpmCentConeZ + (r/piaacmc[0].fpmCentConeRad)*(0.5*(piaacmc[0].fpmminsag+piaacmc[0].fpmmaxsag)-piaacmc[0].fpmCentConeZ);
                                            // piaacmc[0].fpmCentConeZ*(piaacmc[0].fpmCentConeRad-r)/(piaacmc[0].fpmCentConeRad); //piaacmc[0].fpmCentConeZ
                                            a = 1.0;
                                        }
                                    }

                                    // Zone number
                                    zi = (long) (data.image[IDz].array.U[jj1*piaacmc[0].fpmarraysize+ii1]);
                                    if(zi-1>data.image[piaacmc[0].zonezID].md[0].size[0]-1)
                                    {
                                        printf("ERROR: Zone %d does not exist (image %s has size %ld %ld)   pix %ld %ld   %ld\n", zi, data.image[piaacmc[0].zonezID].md[0].name, data.image[piaacmc[0].zonezID].md[0].size[0], data.image[piaacmc[0].zonezID].md[0].size[1], jj1, jj1, piaacmc[0].fpmarraysize);
                                        exit(0);
                                    }
                                    if(zi>0)
                                    {
                                        t = data.image[piaacmc[0].zonezID].array.D[zi-1]; // thickness
                                        a = data.image[piaacmc[0].zoneaID].array.D[zi-1]; // amplitude transmission
                                    }
                                }


                                if(mode == -1)   // make 1-fpm
                                {
                                    // if(zi>0.1)
                                    //{
                                    amp = a;
                                    pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, t, optsyst[0].lambdaarray[k]);

                                    retmp += 1.0-amp*cos(pha);
                                    imtmp += amp*sin(pha);
                                    //								data.image[ID].array.CF[k*size2+jj*size+ii].re = 1.0-re;
                                    //								data.image[ID].array.CF[k*size2+jj*size+ii].im = -im;

                                }
                                else // impulse response from single zone
                                {
                                    if(mode == zi)
                                    {
                                        amp = 1.0;
                                        pha = 0.0; //OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, t, optsyst[0].lambdaarray[k]);
                                        retmp += amp*cos(pha);
                                        imtmp += amp*sin(pha);
                                        //								data.image[ID].array.CF[k*size2+jj*size+ii].re = 1.0;
                                        //								data.image[ID].array.CF[k*size2+jj*size+ii].im = 0.0;
                                    }
                                }

                                //bad location
                                ttmp += t;
                                zonetmp += 1.0*zi;
                            }
                        }


                        data.image[ID].array.CF[k*size2+jj*size+ii].re = retmp/(NBsubPix*NBsubPix);
                        data.image[ID].array.CF[k*size2+jj*size+ii].im = imtmp/(NBsubPix*NBsubPix);
                        data.image[IDsag].array.F[k*size2+jj*size+ii] = ttmp/(NBsubPix*NBsubPix);
                        data.image[IDzone].array.F[k*size2+jj*size+ii] = zonetmp/(NBsubPix*NBsubPix);
                    }
                    else
                    {
                        if(mode == -1)   // make 1-fpm
                        {
                            amp = a;
                            pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, t, optsyst[0].lambdaarray[k]);

                            data.image[ID].array.CF[k*size2+jj*size+ii].re =  1.0-amp*cos(pha);
                            data.image[ID].array.CF[k*size2+jj*size+ii].im = amp*sin(pha);
                        }
                        else
                        {
                            data.image[ID].array.CF[k*size2+jj*size+ii].re = 0.0;
                            data.image[ID].array.CF[k*size2+jj*size+ii].im = 0.0;
                        }
                        data.image[IDsag].array.F[k*size2+jj*size+ii] = t;
                        data.image[IDzone].array.F[k*size2+jj*size+ii] = 0.0;
                    }
                }
        }
# ifdef HAVE_LIBGOMP
    }
# endif




   if(saveMask==1)
    {
        /* save mask sag */
        save_fits("fpmsag", "!tmp_fpmsag.fits");
        
        sprintf(fname, "!%s/fpm_sagmap_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
		save_fits("fpmsag", fname);


        /* save zones */
        save_fits("fpmzone", "!tmp_fpmzone.fits");


        /* save mask transmission */
        IDm = create_3DCimage_ID("fpmCA", size, size, nblambda);
        for(k=0; k<nblambda; k++)
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    data.image[IDm].array.CF[k*size2+jj*size+ii].re = 1.0-data.image[ID].array.CF[k*size2+jj*size+ii].re;
                    data.image[IDm].array.CF[k*size2+jj*size+ii].im = data.image[ID].array.CF[k*size2+jj*size+ii].im;
                    // [re,im] = 1-fpm  -> fpm = [(1.0-re), im]
                }


        mk_amph_from_complex("fpmCA", "tfpma", "tfpmp", 0);
        delete_image_ID("fpmCA");
        save_fits("tfpma", "!tmp_fpmCA_ampl.fits");
        save_fits("tfpmp", "!tmp_fpmCA_pha.fits");
        delete_image_ID("tfpma");
        delete_image_ID("tfpmp");
    }

    delete_image_ID("fpmsag");

    return(ID);
}























///
/// initializes the optsyst structure to simulate PIAACMC system
/// Fills in an OPTSYST global optsyst (see OptSysProp.h) which describes
/// the optical system as a series of planes based on the input design structure
///

// TTxld and TTyld are tip/tilt x-y coordinates specifying the location of the source relative
// to the optical axis in units of lambda/D
// index allows multiple configurations, but it's always 0.  Nonzero values are untested
void PIAACMCsimul_init( OPTPIAACMCDESIGN *design, long index, double TTxld, double TTyld )
{
    FILE *fp;
    FILE *fpri;
    long k, i;
    long size;
    double x, y, PA;
    long ii, jj;
    long nblambda;
    long size2;
    double beamradpix;
    long kx, ky, kxy;
    long IDpiaaz0, IDpiaaz1;
    long surf;
    long IDa;
    char fname_pupa0[500];
    long ID;
    long elem;
    char fname[500];
    long IDv;
    long IDopderr;
    long iDM; // DM index

    int savefpm;

    //    double ri, ri0, sag2opd_coeff;
    //    long IDpiaar0zsag, IDpiaar1zsag;
    //    int mkpiaar0zsag, mkpiaar1zsag;
    //    double sag2opd_coeff0;
    int IDpiaam0z, IDpiaam1z;

    int ret;
    char command[1000];

    assert(index == 0);   // test that index is always 0

    optsyst[0].nblambda = design[index].nblambda;
    nblambda = optsyst[0].nblambda;

    if(PIAACMC_save==1)
    {
        sprintf(fname, "%s/lambdalist.txt", piaacmcconfdir);
        fp = fopen(fname, "w");
    }


    printf("lambda = %g\n", design[index].lambda);
    printf("LAMBDASTART = %g\n", LAMBDASTART);
    printf("LAMBDAEND = %g\n", LAMBDAEND);

    // sets up the wavelengths over specifed bandwidth
    for(k=0; k<optsyst[0].nblambda; k++)
    {
        optsyst[0].lambdaarray[k] = LAMBDASTART + (0.5+k)*(LAMBDAEND-LAMBDASTART)/optsyst[0].nblambda;
        if(PIAACMC_save==1)
            fprintf(fp, "%02ld %20g\n", k, optsyst[0].lambdaarray[k]);
    }
    if(PIAACMC_save==1)
        fclose(fp);

    // the physical pupil size in meters
    optsyst[0].beamrad = design[index].beamrad;
    // the number of pixels in each side of the square complex amplitude arrays at each plane
    optsyst[0].size = design[index].size;
    size = optsyst[0].size;
    size2 = size*size; // area
    optsyst[0].pixscale = design[index].pixscale;
    optsyst[0].DFTgridpad = 0; // 0 for full DFT sampling, >0 for faster execution
    // beam radius in pixels
    beamradpix = optsyst[0].beamrad/optsyst[0].pixscale;




    // printf("BEAM RADIUS = %f / %f  = %f pix,   piaacmc[0].beamrad = %f\n", optsyst[0].beamrad, optsyst[0].pixscale, beamradpix, piaacmc[0].beamrad );
    // sleep(10);

    // parameter that determines sampling of DFTs for progation onto the FPM
    // 0 => full sampling
    // 1 => every other pixel in each dimension
    // 2 => every third pixel in each dimension
    // etc: n => every (n+1)th pixel in each dimension
    // allows subsampling to speed up the DFT computation
    if((IDv=variable_ID("PIAACMC_dftgrid"))!=-1)
        optsyst[0].DFTgridpad = (long) (data.variable[IDv].value.f+0.001);


    // define optical elements and locations
    // have at least two aspheric mirrors in addition to the DMs
    // tyically design[index].nbDM = 0
    optsyst[0].NB_asphsurfm = 2+design[index].nbDM;
    // no aspheric lenses
    optsyst[0].NB_asphsurfr = 0;

    optsyst[0].NBelem = 100; // to be updated later

    if(PIAACMC_save==1)
    {
        sprintf(fname, "%s/conjugations.txt", piaacmcconfdir);
        fp = fopen(fname, "w");
    }





    elem = 0;
    // ------------------- elem 0: input pupil -----------------------
    sprintf(optsyst[0].name[elem], "input pupil");
    optsyst[0].elemtype[elem] = 1; // pupil mask
    // input pupil from file - will always exist
    sprintf(fname_pupa0, "%s/pupa0_%ld.fits", piaacmcconfdir, size);

    if(file_exists(fname_pupa0)==1)
        load_fits(fname_pupa0, "pupa0", 1);

    IDa = image_ID("pupa0");
    if(IDa==-1) // if pupil does not exist, use circular one (this occurs in initial design steps)
    { 
        printf("CREATING INPUT PUPIL\n");
        if(IDa!=-1)
            delete_image_ID("pupa0");
        IDa = create_3Dimage_ID("pupa0", size, size, nblambda);

        ID = image_ID("telpup");
        if(ID==-1)
            if(file_exists("telpup.fits")==1)
                ID = load_fits("telpup.fits", "telpup", 1);


        if(ID==-1)
        {
            for(k=0; k<nblambda; k++)
                for(ii=0; ii<size2; ii++)
                {
                    if((data.image[IDr].array.F[ii]>design[index].centObs0)&&(data.image[IDr].array.F[ii]<1.0))
                        data.image[IDa].array.F[k*size2+ii] = 1.0;
                    else
                        data.image[IDa].array.F[k*size2+ii] = 0.0;
                }
        }
        else
            for(k=0; k<nblambda; k++)
                for(ii=0; ii<size2; ii++)
                {
                    if(data.image[ID].array.F[ii]>0.5)
                        data.image[IDa].array.F[k*size2+ii] = 1.0;
                    else
                        data.image[IDa].array.F[k*size2+ii] = 0.0;
                }

        sprintf(fname_pupa0, "!%s/pupa0_%ld.fits", piaacmcconfdir, size);
        save_fl_fits("pupa0", fname_pupa0);
    }
    optsyst[0].elemarrayindex[elem] = IDa;
    optsyst[0].elemZpos[elem] = 0.0; // pupil is at z = 0
    if(PIAACMC_save==1)
        fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
    elem++;







    // set up the shape of a mirror to insert tip/tilt and optical error
    
    // initialize this mirror by setting pointing (simulated as mirror shape), defining off-axis source
    ID = create_2Dimage_ID("TTm", size, size);

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = (1.0*ii-0.5*size)/beamradpix; // position of pixel (ii,jj) in pupil radius units
            y = (1.0*jj-0.5*size)/beamradpix;
            // set the mirror shape as a linear tilt of pixel position reflecting the tilt
            data.image[ID].array.F[jj*size+ii] = 0.25*(TTxld*x+TTyld*y)*(LAMBDAEND+LAMBDASTART)*0.5; // xld -> half-OPD
        }


    // add OPD error on TTM if it exists
    IDopderr = image_ID("opderr");
    if(IDopderr != -1)
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                x = (1.0*ii-0.5*size)/beamradpix;
                y = (1.0*jj-0.5*size)/beamradpix;
                // add the error shape to the mirror shape
                data.image[ID].array.F[jj*size+ii] += data.image[IDopderr].array.F[jj*size+ii]*0.5;
            }
    }

    // sprintf(fname, "!%s/TTm.fits", piaacmcconfdir);
    // save_fits("TTm", fname);

    // finish the definition of the TT mirror specifying various properties
    sprintf(optsyst[0].name[elem], "TT mirror");
    optsyst[0].elemtype[elem] = 3; // reflective mirror
    optsyst[0].elemarrayindex[elem] = 0; // not used because this field is only relevant for
                                        // DM or aspheric mirrors, which this mirror is not
                                        // this mirror is "flat" except for possible injected OPD error
    optsyst[0].ASPHSURFMarray[0].surfID = ID; // store array ID
    optsyst[0].elemZpos[elem] = 0.0; // put it at the entrance pupil
    if(PIAACMC_save==1)
        fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
    //        fprintf(fp,"%02ld  %f    Fold mirror used to induce pointing offsets\n", elem, optsyst[0].elemZpos[elem]);
    elem++;




    // set up the deformable mirrors
    // tyically design[index].nbDM = 0 so we will skip this
    for(iDM=0; iDM<design[index].nbDM; iDM++)
    {
        // ----------------- DM (s) -----------------------------------------------
        sprintf(optsyst[0].name[elem], "DM %ld", iDM);
        optsyst[0].elemtype[elem] = 3; // reflective element
        optsyst[0].elemarrayindex[elem] = 3+iDM; // index
        optsyst[0].ASPHSURFMarray[optsyst[0].elemarrayindex[elem]].surfID = design[index].ID_DM[iDM];
        optsyst[0].elemZpos[elem] = design[index].DMpos[iDM];
        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //          fprintf(fp,"%02ld  %f    DM %ld\n", elem, optsyst[0].elemZpos[elem], iDM);
        elem++;
    }




    // shape/sag for the first aspheric mirror
    IDpiaam0z = image_ID("piaam0z");  // nominal sag (mirror equivalent)
    // shape/sag for the second aspheric mirror
    IDpiaam1z = image_ID("piaam1z");  //







    // ------------------- [OPTIONAL] pre-apodizer  -----------------------
    // typically not present, skipping
    ID = image_ID("prePIAA0mask");
    if(ID==-1)
        ID = load_fits("prePIAA0mask.fits", "prePIAA0mask", 1);
    if(ID!=-1)
    {
        // tell the design that this element exists (was found on disk)
        design[index].prePIAA0mask = 1;
        sprintf(optsyst[0].name[elem], "pupil plane apodizer");
        optsyst[0].elemtype[elem] = 1; // opaque mask
        optsyst[0].elemarrayindex[elem] = ID;
        optsyst[0].elemZpos[elem] = design[index].prePIAA0maskpos;

        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        elem++;
    }
    else
        design[index].prePIAA0mask = 0;

    // PIAAmode = 1 => this is a PIAA system
    // PIAAmode = 0 => this is not a PIAA system
    if(piaacmc[0].PIAAmode == 1)
    {
        // ------------------- elem 2:  PIAA M/L 0  -----------------------
        // (M/L is "mirror or lens" - in our case it's a mirror)
        sprintf(optsyst[0].name[elem], "PIAA optics 0");
        optsyst[0].elemtype[elem] = 3; // reflective PIAA M0
        // this is an aspheric mirror, so we need actual sag shapes
        optsyst[0].elemarrayindex[elem] = 1; // index = 1 implied aspheric
        optsyst[0].elemZpos[elem] = design[index].PIAA0pos;  // location of this element relative to pupil

        printf("============ (2) PIAA0pos = %f ==================\n", optsyst[0].elemZpos[elem]);
//        sleep(5);
        // set up the element properties
        if(design[index].PIAAmaterial_code == 0) // mirror
            // specify the sag array, put in data global by routine named something like createPIAA_mirror_shapes
            optsyst[0].ASPHSURFMarray[optsyst[0].elemarrayindex[elem]].surfID = IDpiaam0z;
        else // lens
        {
            optsyst[0].elemtype[elem] = 4;
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].surfID = image_ID("piaar0zsag"); //IDpiaar0zsag;
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].mat0 = 100;
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].mat1 = design[0].PIAAmaterial_code; // vacuum
        }
        // make sure the above did something
        if(optsyst[0].ASPHSURFMarray[1].surfID==-1)
        {
            printf("ERROR: surface 0 not identified\n");
            list_image_ID();
            exit(0);
        }
        // print this element to tracking file if desired
        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        elem++;
    }





    // ------------------- [OPTIONAL] opaque mask after last elem -----------------------
    // get opaque mask from the file, with a standard filename for the first mask
    // we don't have one in the nominal design
    ID = load_fits("postPIAA0mask.fits", "postPIAA0mask", 1);
    if(ID!=-1)
    {
        // tell the design that this element exists (was found on disk)
        design[index].postPIAA0mask = 1;
        sprintf(optsyst[0].name[elem], "opaque mask after PIAA element 0");
        optsyst[0].elemtype[elem] = 1; // opaque mask
        optsyst[0].elemarrayindex[elem] = ID;
        optsyst[0].elemZpos[elem] = design[index].postPIAA0maskpos; // get position from design input

        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        elem++;
    }
    else
        design[index].postPIAA0mask = 0;



    // if we're a PIAA
    if(piaacmc[0].PIAAmode == 1)
    {
        // add one more mirror and mask
        // ------------------- elem 3: reflective PIAA M1  -----------------------
        sprintf(optsyst[0].name[elem], "PIAA optics 1");
        optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
        optsyst[0].elemarrayindex[elem] = 2;
        optsyst[0].elemZpos[elem] = design[index].PIAA0pos + design[index].PIAAsep;

        if(design[index].PIAAmaterial_code == 0) // mirror
            optsyst[0].ASPHSURFMarray[2].surfID = IDpiaam1z;
        else // lens
        {
            optsyst[0].elemtype[elem] = 4;
            optsyst[0].ASPHSURFRarray[2].surfID = image_ID("piaar1zsag"); //IDpiaar0zsag;
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].mat0 = 100; // vacuum
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].mat1 = design[0].PIAAmaterial_code;
        }


        if(optsyst[0].ASPHSURFMarray[2].surfID==-1)
        {
            printf("ERROR: surface 1 not identified\n");
            list_image_ID();
            exit(0);
        }

        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //       fprintf(fp,"%02ld  %f    PIAAM1\n", elem, optsyst[0].elemZpos[elem]);
        elem++;





        // ------------------- elem 4 opaque mask at reflective PIAA M1  -----------------------
        sprintf(optsyst[0].name[elem], "opaque mask at PIAA elem 1");
        optsyst[0].elemtype[elem] = 1; // opaque mask
        ID = load_fits("piaa1mask.fits", "piaa1mask", 1);
        if(ID==-1)
            ID = make_disk("piaa1mask", size, size, 0.5*size, 0.5*size, design[index].r1lim*beamradpix);
        optsyst[0].elemarrayindex[elem] = ID;
        optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1];

        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //        fprintf(fp,"%02ld  %f    PIAAM1 edge opaque mask\n", elem, optsyst[0].elemZpos[elem]);
        elem++;
    }



    // --------------------  elem 5: focal plane mask ------------------------
    if((IDv=variable_ID("PIAACMC_NOFPM"))==-1)
    {
        sprintf(optsyst[0].name[elem], "post focal plane mask pupil");
        optsyst[0].elemtype[elem] = 5; // focal plane mask
        optsyst[0].elemarrayindex[elem] = 0;

        printf("=========== MAKE FOCAL PLANE MASK ===========\n");
        //sleep(5);


        savefpm = 0;
        if((IDv=variable_ID("PIAACMC_SAVE_fpm"))!=-1)
            savefpm = (int) (data.variable[IDv].value.f+0.001);

        // make the focal plane mask here
        optsyst[0].FOCMASKarray[0].fpmID = PIAACMCsimul_mkFocalPlaneMask("fpmzmap", "piaacmcfpm", focmMode, savefpm); // if -1, this is 1-fpm; otherwise, this is impulse response from single zone


        // zfactor is the zoom factor for the DFT, driving sample resolution in the z direction
        // to allow faster DFTs.  Similar to DFTgridpad.
        optsyst[0].FOCMASKarray[0].zfactor = design[index].fpzfactor;
        // set the position of the pupil from which the DFT propagates to the FPM
        // NOT the position along the beam of the FPM.  That's OK and intended.
        // For this element, this defines the conjugation of the pupil from which we are computing the DFT
        optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1]; // plane from which FT is done
        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //      fprintf(fp,"%02ld  %f    post-focal plane mask pupil\n", elem, optsyst[0].elemZpos[elem]);
        elem++;
    }



    // if we're a PIAA
    if(piaacmc[0].PIAAmode == 1)
    {
        // if the the inverse PIAA is prior to the Lyot stop
        // (invPIAAmode = 0 has no inverse PIAA, = 1 has Lyot stop prior to inverse PIAA)
        // there is no inverse PIAA in the WFIRST design
        if(design[index].invPIAAmode == 2) // inv PIAA -> Lyot stops
        {
            // --------------------  elem 8: inv PIAA1 ------------------------
            sprintf(optsyst[0].name[elem], "invPIAA optics 1");

            if(design[index].PIAAmaterial_code == 0) // mirror
                optsyst[0].elemtype[elem] = 3; // reflective PIAA M/L 1
            else
                optsyst[0].elemtype[elem] = 4; // refractive PIAA M/L 1

            optsyst[0].elemarrayindex[elem] = 2;
            // put an element at z=0 in conjugation space (conjugate to the pupil)
            optsyst[0].elemZpos[elem] = 0.0;
            if(PIAACMC_save==1)
                fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
            //          fprintf(fp,"%02ld  %f    invPIAA1\n", elem, optsyst[0].elemZpos[elem]);
            elem++;

            // --------------------  elem 9: inv PIAA0 ------------------------
            sprintf(optsyst[0].name[elem], "invPIAA optics 0");

            if(design[index].PIAAmaterial_code == 0) //  mirror
                optsyst[0].elemtype[elem] = 3; // reflective PIAA M/L 0
            else
                optsyst[0].elemtype[elem] = 4; // refractive PIAA M/L 0

            optsyst[0].elemarrayindex[elem] = 1;
            // previous element + PIAAsep
            optsyst[0].elemZpos[elem] = design[index].PIAAsep;
            if(PIAACMC_save==1)
                fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
            //         fprintf(fp,"%02ld  %f    invPIAA0\n", elem, optsyst[0].elemZpos[elem]);
            elem++;
        }
    }


    // --------------------  Lyot masks  ------------------------
    // add Lyot masks as specified in the design
    for(i=0; i<design[index].NBLyotStop; i++)
    {
        sprintf(optsyst[0].name[elem], "Lyot mask %ld", i);
        optsyst[0].elemtype[elem] = 1; // Lyot mask
        optsyst[0].elemarrayindex[elem] = design[index].IDLyotStop[i];
        printf("elem %ld  Lyot mask %ld : %ld\n", elem, i, design[index].IDLyotStop[i]);
        optsyst[0].elemZpos[elem] =  design[index].LyotStop_zpos[i];
        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //          fprintf(fp,"%02ld  %f  Lyot Stop %ld\n", elem, optsyst[0].elemZpos[elem], i);
        elem++;
    }

    // if we're a PIAA
    if(piaacmc[0].PIAAmode == 1)
    {
        // add stops for the inverse PIAA
        // not in WFIRST design, skipping...
        if(design[index].invPIAAmode == 1) // Lyot masks -> inv PIAA
        {
            // --------------------  elem 8: inv PIAA1 ------------------------
            sprintf(optsyst[0].name[elem], "invPIAA optics 1");
            if(design[index].PIAAmaterial_code == 0) // mirror
                optsyst[0].elemtype[elem] = 3; // reflective PIAA M/L 1
            else
                optsyst[0].elemtype[elem] = 4; // refractive PIAA M/L 1
            optsyst[0].elemarrayindex[elem] = 2;
            optsyst[0].elemZpos[elem] = 0.0;
            if(PIAACMC_save==1)
                fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
            //           fprintf(fp,"%02ld  %f    invPIAA1\n", elem, optsyst[0].elemZpos[elem]);
            elem++;

            // --------------------  elem 9: inv PIAA0 ------------------------
            sprintf(optsyst[0].name[elem], "invPIAA optics 0");
            if(design[index].PIAAmaterial_code == 0) //  mirror
                optsyst[0].elemtype[elem] = 3; // reflective PIAA M/L 0
            else
                optsyst[0].elemtype[elem] = 4; // refractive PIAA M/L 0
            optsyst[0].elemarrayindex[elem] = 1;
            optsyst[0].elemZpos[elem] = design[index].PIAAsep;
            if(PIAACMC_save==1)
                fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
            //           fprintf(fp,"%02ld  %f    invPIAA0\n", elem, optsyst[0].elemZpos[elem]);
            elem++;
        }
    }


    // if we're a PIAA
    if(piaacmc[0].PIAAmode == 1)
    {
        // --------------------  elem 9: back end mask  ------------------------
        // not in WFIRST design, skipping, but it looks very straightforward

        sprintf(optsyst[0].name[elem], "back end pupil stop  (rad = %f)", design[index].pupoutmaskrad);

        optsyst[0].elemtype[elem] = 1;
        ID = make_disk("pupoutmask", size, size, 0.5*size, 0.5*size, design[index].pupoutmaskrad*design[index].beamrad/design[index].pixscale);
        optsyst[0].elemarrayindex[elem] = ID;
        optsyst[0].elemZpos[elem] =  optsyst[0].elemZpos[elem-1];
        if(PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //     fprintf(fp,"%02ld  %f   back end mask\n", elem, optsyst[0].elemZpos[elem]);
        elem++;
    }



    if(PIAACMC_save==1)
        fclose(fp);

    optsyst[0].NBelem = elem;
    optsyst[0].endmode = 0;

    optsystinit = 1;
}













//
// RADIAL PIAACMC SYSTEM DESIGN (geometrical optics)
//



//
// load and fit radial apodization profile
// modal basis is mk(r) : cos(r*k*M_PI/1.3)
//
int PIAACMCsimul_load2DRadialApodization(char *IDapo_name, float beamradpix, char *IDapofit_name)
{
    long NBpts;
    long IDm;
    long sizem;
    long kmax = 10;
    long ID, IDmask, IDin;
    long ii, jj;
    long offset;
    long sizein;
    float eps = 1.0e-4;
    char fname[500];
    int ret;
    char command[1000];
    int debug = 0;


    sizem = (long) (beamradpix*2);

    // CREATE MODES IF THEY DO NOT EXIST
    if((IDm=image_ID("APOmodesCos"))==-1)
    {
        IDm = linopt_imtools_makeCosRadModes("APOmodesCos", sizem, kmax, ApoFitCosFact*beamradpix, 1.0);
        sprintf(fname, "!%s/APOmodesCos.fits", piaacmcconfdir);
        save_fits("APOmodesCos", fname);
    }

    // CREATE MASK AND CROP INPUT
    IDmask = create_2Dimage_ID("fitmaskapo", sizem, sizem);

    IDin = image_ID(IDapo_name);
    sizein = data.image[IDin].md[0].size[0];
    ID = create_2Dimage_ID("_apoincrop", sizem, sizem);
    offset = (sizein-sizem)/2;
    for(ii=0; ii<sizem; ii++)
        for(jj=0; jj<sizem; jj++)
        {
            data.image[ID].array.F[jj*sizem+ii] = data.image[IDin].array.F[(jj+offset)*sizein+(ii+offset)];
            if((data.image[ID].array.F[jj*sizem+ii]>eps)&&(ii%1==0)&&(jj%1==0))
                data.image[IDmask].array.F[jj*sizem+ii] = 1.0;
        }

    if(debug==1)
    {
        sprintf(fname, "!%s/_apoincrop.fits", piaacmcconfdir);
        save_fits("_apoincrop", fname);

        sprintf(fname, "!%s/fitmaskapo.fits", piaacmcconfdir);
        save_fits("fitmaskapo", fname);
    }

    linopt_imtools_image_fitModes("_apoincrop", "APOmodesCos", "fitmaskapo", 1.0e-8, IDapofit_name, 0);
    sprintf(command, "mv %s/eigenv.dat %s/eigenv_APOmodesCos.dat", piaacmcconfdir, piaacmcconfdir);
    ret = system(command);

    if(debug==1) // test fit quality
    {
        linopt_imtools_image_construct("APOmodesCos", IDapofit_name, "testapofitsol");

        sprintf(fname, "!%s/testapofitsol.fits", piaacmcconfdir);
        save_fits("testapofitsol", fname);

        arith_image_sub("_apoincrop", "testapofitsol", "apofitres");
        arith_image_mult("apofitres", "fitmaskapo", "apofitresm");

        sprintf(fname, "!%s/apofitres.fits", piaacmcconfdir);
        save_fits("apofitres", fname);

        sprintf(fname, "!%s/apofitresm.fits", piaacmcconfdir);
        save_fits("apofitresm", fname);

        // linopt_imtools_image_fitModes("apofitres", "APOmodesCos", "fitmaskapo", 1.0e-5, "test2c", 0);
        info_image_stats("apofitresm", "");
    }

    delete_image_ID("_apoincrop");
    delete_image_ID("fitmaskapo");


    return 0;
}











/**
 * computes radial PIAA optics sag
 *
 * this function only works for circular PIAA
 * uses radial PIAACMC design to initialize PIAA optics shapes and focal plane mask
 */

int PIAACMCsimul_init_geomPIAA_rad(char *IDapofit_name)
{
    long i, ii, k;
    double *pup0;
    double *pup1;
    double *flux0cumul;
    double *flux1cumul;

    long IDcoeff;
    long nbcoeff;
    double r;
    FILE *fp;
    double total;

    // to convert r ro r1 (assymptotic outer radius on pup1)
    double coeffa = 3.0; // convergence rate from linear to assymptotic value
    double coeffa1 = 0.5; // convergence radius limit (added to 1.0)

    double r0, r1;
    double FLUX0_in = 0.0; // inside central obstruction
    double FLUX0_out = 0.0; // outside beam edge
    double FLUX1_in = 0.0; // inside central obstruction
    double FLUX1_out = 0.0; // outside beam edge
    double normcoeff;

    // inner profile adjustment
    double a, b, verr, bstep, x, eps1;
    double t0, t0cnt, value;
    int dir, odir;
    long NBistep;
    long iioffset;

    double fluxdens, F0, F1, F2;
    double dr0, dr1, ndr0, ndr1;

    long NBpoints;
    double *piaar00;
    double *piaar11;
    double *piaar01;
    double *piaar10;
    long cnt;
    double tmp;
    double epsilon = 0.000000000001;

    double *piaaM0z;
    double *piaaM1z;
    double r0c, r1c, dx, dy, dist, y3, r0n, slope, dz;

    char fname[500];

    pup0 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    pup1 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    flux0cumul = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    flux1cumul = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);


    // STEP 1: CREATE AMPLITUDE AND CUMULATIVE INTENSITY PROFILES


    // CREATE OUTPUT AMPLITUDE APODIZATION PROFILE AND ITS CUMUL

    IDcoeff = image_ID(IDapofit_name);
    nbcoeff = data.image[IDcoeff].md[0].size[0];
    printf("%ld coefficients\n", nbcoeff);

    total = 0.0;
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        pup1[ii] = 0.0;
        r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
        if(r<1.0)
            r1 = r;
        else
            r1 = 1.0 + (r-1) / pow((1.0 + pow(1.0/coeffa1 * (r-1),coeffa)), 1.0/coeffa);

        for(k=0; k<nbcoeff; k++)
            pup1[ii] += data.image[IDcoeff].array.F[k]*cos(r1*k*M_PI/ApoFitCosFact);
        if(r<piaacmc[0].centObs1)
            FLUX1_in += pup1[ii]*pup1[ii]*r;
        if(r>1.0)
            FLUX1_out += pup1[ii]*pup1[ii]*r;
        total += pup1[ii]*pup1[ii]*r;
        flux1cumul[ii] = total;
    }


    normcoeff = 1.0/(total-FLUX1_in-FLUX1_out);

    FLUX1_in *= normcoeff;
    FLUX1_out *= normcoeff;
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
        flux1cumul[ii] *= normcoeff;
    }


    printf("outer fluxes 1: %lf %lf\n", FLUX1_in, FLUX1_out);




    // CREATE FLUX0

    total = 0.0;
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
        pup0[ii] = 1.0;

        if(r<piaacmc[0].centObs0)
            FLUX0_in += pup0[ii]*pup0[ii]*r;
        if(r>1.0)
            FLUX0_out += pup0[ii]*pup0[ii]*r;

        total += pup0[ii]*pup0[ii]*r;
        flux0cumul[ii] = total;
    }
    normcoeff = 1.0/(total-FLUX0_in-FLUX0_out);

    FLUX0_in *= normcoeff;
    FLUX0_out *= normcoeff;

    printf("outer fluxes 0: %lf (%lf)    %lf\n", FLUX0_in, FLUX1_in, FLUX0_out);


    //
    // Compute inner pseudo profile
    //
    b = 0.5;
    bstep = 0.1;
    verr = 1.0;
    NBistep = piaacmc[0].centObs0*piaacmc[0].NBradpts/piaacmc[0].r0lim;
    //  innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);
	dir = 1.0; // initial direction
    while(fabs(verr)>1.0e-9)
    {
        t0 = 0.0;
        t0cnt = 0.0;
        for(ii=0; ii<NBistep; ii++)
        {
            x = 1.0*ii/NBistep;
            r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
            a = 0.5;
            eps1 = 1e-8;
            if(x<eps1)
                x = eps1;
            if(x>1.0-eps1)
                x = 1.0-eps1;
            pup0[ii] =  b + (1.0-b)*(0.5+atan(-a/b/(x*x) + a/pow(x-1.0,2))/M_PI);
            t0 += r*pup0[ii]* pup0[ii];
            flux0cumul[ii] = t0;
        }

        verr = t0*normcoeff - FLUX1_in;

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
    printf("TOTAL = %f -> %g (%g %g)\n", b, t0*normcoeff, bstep, verr);


    // outer region
    b = 0.5;
    bstep = 0.1;
    verr = 1.0;
    NBistep = piaacmc[0].NBradpts*(piaacmc[0].r0lim-1.0)/piaacmc[0].r0lim;
    //  innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);
    iioffset = (long) (1.0*piaacmc[0].NBradpts/piaacmc[0].r0lim);
    NBistep = piaacmc[0].NBradpts-iioffset;
    while(fabs(verr)>1.0e-9)
    {
        t0 = 0.0;
        t0cnt = 0.0;
        for(ii=0; ii<NBistep; ii++)
        {
            x = 1.0-1.0*ii/NBistep;
            r = 1.0+1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
            a = 0.5;
            eps1 = 1e-8;
            if(x<eps1)
                x = eps1;
            if(x>1.0-eps1)
                x = 1.0-eps1;
            pup0[ii+iioffset] =  b + (1.0-b)*(0.5+atan(-a/b/(x*x) + a/pow(x-1.0,2))/M_PI);
            t0 += r*pup0[ii+iioffset]* pup0[ii+iioffset];
            flux0cumul[ii+iioffset] = t0;
        }

        verr = t0*normcoeff - FLUX1_out;

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
    printf("TOTAL = %f -> %g (%g %g)\n", b, t0*normcoeff, bstep, verr);



    total = 0.0;
    FLUX0_in = 0.0;
    FLUX0_out = 0.0;
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
        if(r<piaacmc[0].centObs0)
            FLUX0_in += pup0[ii]*pup0[ii]*r;
        if(r>1.0)
            FLUX0_out += pup0[ii]*pup0[ii]*r;

        total += pup0[ii]*pup0[ii]*r;
        flux0cumul[ii] = total;
        flux0cumul[ii] *= normcoeff;;
    }
    FLUX0_in *= normcoeff;
    FLUX0_out *= normcoeff;

    printf("outer fluxes 0: %lf (%lf)    %lf\n", FLUX0_in, FLUX1_in, FLUX0_out);




    sprintf(fname, "%s/pup01.prof", piaacmcconfdir);
    fp = fopen(fname, "w");
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        r0 = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
        r1 = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
        fprintf(fp, "%f %f %g %g %g %g\n", r0, r1, pup0[ii], pup1[ii], flux0cumul[ii], flux1cumul[ii]);
    }
    fclose(fp);




    // STEP 2: COMPUTE r0 - r1 CORRESPONDANCE

    piaar00 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r0 as a function of r0 index
    piaar11 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r1 as a function of r1 index
    piaar10 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r1 as a function of r0 index
    piaar01 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r0 as a function of r1 index

    /* computing r0 and r1 */
    /* r0 and r1 are dimensionless */

    /* first, r0 is evenly distributed on the first optic */
    for(i=0; i<piaacmc[0].NBradpts; i++)
    {
        piaar00[i] = piaacmc[0].r0lim*i/piaacmc[0].NBradpts;
        piaar11[i] = piaacmc[0].r1lim*i/piaacmc[0].NBradpts;
    }

    i=0;
    ii=0;
    cnt = 0;
    piaar00[0] = 0.0;
    piaar10[0] = 0.0;
    //  fp = fopen("test0.txt", "w");
    for(i=1; i<piaacmc[0].NBradpts; i++)
    {
        F0 = flux0cumul[i];
        while((flux1cumul[ii]<flux0cumul[i])&&(ii<piaacmc[0].NBradpts))
            ii++;
        F1 = flux1cumul[ii-1];
        F2 = flux1cumul[ii];

        /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
        if(fabs(F2-F1)>0.0000000001)
            fluxdens = (F2-F1)/(2.0*ii-1.0);
        else
            fluxdens = 0.0000000001;
        x = sqrt((F0-F1)/fluxdens+(1.0*ii*ii-2.0*ii+1.0)) + 1.0 - 1.0*ii;

        piaar10[i] = piaacmc[0].r1lim*(1.0*ii-1.0+x)/piaacmc[0].NBradpts;
        //  fprintf(fp, "%lf %lf %lf\n", piaar00[i], piaar10[i], F0);
    }
    //  fclose(fp);



    i=0;
    ii=0;
    cnt = 0;
    piaar01[0] = 0.0;
    piaar11[0] = 0.0;
    //  fp = fopen("test1.txt", "w");
    for(i=1; i<piaacmc[0].NBradpts; i++)
    {
        F0 = flux1cumul[i];
        while((flux0cumul[ii]<flux1cumul[i])&&(ii<piaacmc[0].NBradpts))
            ii++;
        F1 = flux0cumul[ii-1];
        F2 = flux0cumul[ii];

        /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
        if(fabs(F2-F1)>0.0000000001)
            fluxdens = (F2-F1)/(2.0*ii-1.0);
        else
            fluxdens = 0.0000000001;
        x = sqrt((F0-F1)/fluxdens+(1.0*ii*ii-2.0*ii+1.0)) + 1.0 - 1.0*ii;

        piaar01[i] = piaacmc[0].r0lim*(1.0*ii-1.0+x)/piaacmc[0].NBradpts;
        //  fprintf(fp, "%lf %lf %lf\n", piaar11[i], piaar01[i], F0);
    }
    //  fclose(fp);





    printf("======== Compute PIAA optics shapes ============\n");
    piaaM0z = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    piaaM1z = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);

    piaaM0z[0] = 0.0;
    piaaM1z[0] = piaacmc[0].PIAAsep;


    for(i=0; i<piaacmc[0].NBradpts-1; i++)
    {
        r0c = piaar00[i];
        r1c = piaar10[i];
        dx = (r0c-r1c)*piaacmc[0].beamrad;
        dz = piaaM1z[i]-piaaM0z[i];
   //     dist = sqrt(dx*dx+dz*dz);
        dist = dz * sqrt(1. + (dx/dz)*(dx/dz)); // preserve sign of dz
        y3 = dist - dz;
        if(fabs(dx)>0.000000001)
            slope = y3/dx;
        else
            slope = 0.0;
        r0n = piaacmc[0].r0lim*(i+1)/piaacmc[0].NBradpts;
        piaaM0z[i+1] = piaaM0z[i] + slope*(r0n-r0c)*piaacmc[0].beamrad;

        if(fabs(dx)>0.000000001)
            slope = y3/dx;
        else
            slope = 0.0;
        piaaM1z[i+1] = piaaM1z[i] + slope*(piaar10[i+1]-r1c)*piaacmc[0].beamrad;
    }

    sprintf(fname, "%s/PIAA_Mshapes.txt", piaacmcconfdir);
    fp = fopen(fname, "w");
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
        fprintf(fp, "%18.16f %18.16f %18.16f %18.16f\n", piaar00[ii]*piaacmc[0].beamrad, piaaM0z[ii], piaar10[ii]*piaacmc[0].beamrad, piaaM1z[ii]);
    fclose(fp);




    free(flux0cumul);
    free(flux1cumul);
    free(pup0);
    free(pup1);

    free(piaar00);
    free(piaar10);
    free(piaar01);
    free(piaar11);


    return(0);
}
















//
// make PIAA OPD screens from radial sag profile
//
int PIAACMCsimul_mkPIAAMshapes_from_RadSag(char *fname, char *ID_PIAAM0_name, char *ID_PIAAM1_name)
{
    FILE *fp;
    long size;
    long ii, jj;
    long ID_PIAAM0, ID_PIAAM1;

    long k;

    double x, y, r, r1;

    double *r0array;
    double *z0array;
    double *r1array;
    double *z1array;

    double alpha;
    double r00, r01;
    double val;

    double beamradpix;
    int ret;


    size = piaacmc[0].size;
    beamradpix = piaacmc[0].beamrad/piaacmc[0].pixscale;
    printf("SIZE = %ld, beamrad = %f pix, sep = %f m\n", size, beamradpix, piaacmc[0].PIAAsep);
    fflush(stdout);



    r0array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    z0array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    r1array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    z1array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);

    fp = fopen(fname, "r");
    for(k=0; k<piaacmc[0].NBradpts; k++)
        ret = fscanf(fp,"%lf %lf %lf %lf\n", &r0array[k], &z0array[k], &r1array[k], &z1array[k]);
    fclose(fp);

    //  for(k=0;k<nbpt;k++)
    //  printf("%ld %.8lf %.8lf %.8lf %.8lf\n", k, r0array[k], z0array[k], r1array[k], z1array[k]);


    for(k=0; k<piaacmc[0].NBradpts; k++)
        z1array[k] -= piaacmc[0].PIAAsep;




    ID_PIAAM0 = create_2Dimage_ID(ID_PIAAM0_name, size, size);
    ID_PIAAM1 = create_2Dimage_ID(ID_PIAAM1_name, size, size);

    printf("\n\n");

# ifdef HAVE_LIBGOMP
    #pragma omp parallel default(shared) private(ii, jj, x, y, r, k, r00, r01, alpha, val)
    {
# endif


# ifdef HAVE_LIBGOMP
        #pragma omp for
# endif
        for(ii=0; ii<size; ii++)
        {
            //      printf("\r %ld / %ld     ", ii, size);
            //fflush(stdout);


            for(jj=0; jj<size; jj++)
            {
                x = (1.0*ii-0.5*size)/beamradpix;
                y = (1.0*jj-0.5*size)/beamradpix;
                r = sqrt(x*x+y*y)*piaacmc[0].beamrad;

                if(r<piaacmc[0].r0lim*piaacmc[0].beamrad)
                {
                    k = 1;
                    while((r0array[k]<r)&&(k<piaacmc[0].NBradpts-2))
                        k++;
                    r00 = r0array[k-1];
                    r01 = r0array[k];
                    alpha = (r-r00)/(r01-r00);
                    if(alpha>1.0)
                        alpha = 1.0;
                    val = (1.0-alpha)*z0array[k-1] + alpha*z0array[k];
                    data.image[ID_PIAAM0].array.F[jj*size+ii] = val;
                }
                else
                    data.image[ID_PIAAM0].array.F[jj*size+ii] = 0.0;

                if(r<piaacmc[0].r1lim*piaacmc[0].beamrad)
                {
                    k = 1;
                    while((r1array[k]<r)&&(k<piaacmc[0].NBradpts-2))
                        k++;
                    r00 = r1array[k-1];
                    r01 = r1array[k];
                    alpha = (r-r00)/(r01-r00);
                    if(alpha>1.0)
                        alpha = 1.0;
                    val = (1.0-alpha)*z1array[k-1] + alpha*z1array[k];
                    data.image[ID_PIAAM1].array.F[jj*size+ii] = -val;//-piaacmc[0].PIAAsep);
                }
                else
                    data.image[ID_PIAAM1].array.F[jj*size+ii] = 0.0;
            }
        }
# ifdef HAVE_LIBGOMP
    }
# endif



    printf("\n\n");

    free(r0array);
    free(z0array);
    free(r1array);
    free(z1array);

    return 0;
}






//
// Detailed simulation of PIAACMC
//
//
//
//


// transmits between rin and rout
long PIAAsimul_mkSimpleLyotStop(char *ID_name, float rin, float rout)
{
    long size;
    long size2;
    long ii, k;
    long ID;
    float r;


    size = piaacmc[0].size;
    size2 = size*size;

    ID = create_3Dimage_ID(ID_name, size, size, piaacmc[0].nblambda);
    for(k=0; k<piaacmc[0].nblambda; k++)
        for(ii=0; ii<size2; ii++)
        {
            if((data.image[IDr].array.F[ii]<rout)&&(data.image[IDr].array.F[ii]>rin))
                data.image[ID].array.F[k*size2+ii] = 1.0;
            else
                data.image[ID].array.F[k*size2+ii] = 0.0;
        }

    return(ID);
}





/**
 * @brief Creates/initializes piaacmcconf structure and directory
 *
 * @param[in] piaacmctype  Type of system: 0=idealized mask, 1=physical mask
 * @param[in] fpmradld     Focal plane mask nominal radius
 * @param[in] centobs0     Input central obstruction
 * @param[in] centobs1     Output central obstruction
 * @param[in] WFSmode      Number of DMs (0: no WFC)
 * @param[in] load         if 1, attempt to load configuration from file
 * 
 * piaacmctype:
 * - 0: if configuration does not exist, create Monochromatic idealized PIAACMC, otherwise, read configuration
 * - 1: physical mask
 *
 * 
 */

int PIAAsimul_initpiaacmcconf(long piaacmctype, double fpmradld, double centobs0, double centobs1, int WFCmode, int load)
{
    FILE *fp;
    float beamradpix;
    long NBpiaacmcdesign = 1;
    long ii, jj, k, i;
    double x, y;
    long size, size0;
    long Cmsize;
    long Fmsize;
    long ID, ID0, ID1;
    long size2;
    int r;
    char command[1000];
    long IDv1, IDv2;
    char fname[500];
    char name[500];
    float tmpf;
    double pha0, t0, t;
    int loaded = 0; // has the configuration been loaded ?

    int saveconf = 0; // if 1, save conf at end of this function
    long IDv;
    int ret;
    double tmplf;

    int IDlscumul;
    long iDM; // DM index
    long *sizearray;


    long IDapo;
    long xsize = 0;
    long ysize = 0;
    long IDapo_PIAA, IDapo_CPA;
    double coeff;




    if(piaacmc == NULL)
    {
        piaacmc = (OPTPIAACMCDESIGN*) malloc(sizeof(OPTPIAACMCDESIGN)*NBpiaacmcdesign);


        // Default Values for PIAACMC (will adopt them unless configuration file exists)
        piaacmc[0].nblambda = 8;

        //piaacmc[0].nblambda = NBLAMBDA;

        // high resolution
        //piaacmc[0].size = 4096;
        //piaacmc[0].pixscale = 0.000055;

        // mid resolution
        piaacmc[0].size = 2048;
        piaacmc[0].pixscale = 0.000055;

        // low resolution
        //piaacmc[0].size = 1024;
        //piaacmc[0].pixscale = 0.00011;


        // very low resolution
        //    piaacmc[0].size = 512;
        // piaacmc[0].pixscale = 0.00022;


        piaacmc[0].beamrad = 0.01; // beam physical radius
        piaacmc[0].PIAA0pos = 1.0; // piaa 0 position [m]
        piaacmc[0].PIAAsep = 1.00; // [m]
        piaacmc[0].fpzfactor = 16.0;
        piaacmc[0].Fratio = 80.0; // default
        strcpy(piaacmc[0].PIAAmaterial_name, "Mirror");  // mirrors
        piaacmc[0].prePIAA0mask = 0;
        piaacmc[0].prePIAA0maskpos = 0.0;
        piaacmc[0].postPIAA0mask = 0;
        piaacmc[0].postPIAA0maskpos = 0.0;
        piaacmc[0].piaaNBCmodesmax =  40;
        piaacmc[0].piaaCPAmax = 10.0;

        piaacmc[0].centObs0 = centobs0; // input central obstruction
        piaacmc[0].centObs1 = centobs1; // output central obstruction
        piaacmc[0].NBradpts = 50000;
        piaacmc[0].r0lim = 1.15; //1425; // outer radius after extrapolation, piaa optics 0
        piaacmc[0].r1lim = 1.5; // outer radius after extrapolation, piaa optics 1


        /// Wavefront control
        piaacmc[0].nbDM = WFCmode; // number of deformable mirrors (10 max)
        for(iDM=0; iDM<piaacmc[0].nbDM; iDM++)
        {
            piaacmc[0].DMpos[iDM] = 0.0 + 0.6*iDM/(0.01+piaacmc[0].nbDM-1.0); // DM conjugation in collimated space
            piaacmc[0].ID_DM[iDM] = -1;  // DM image identifier - to be updated later
        }


        piaacmc[0].NBLyotStop = 2;
        for(i=0; i<10; i++)
        {
            piaacmc[0].LyotStop_zpos[i] = 0.0;
            piaacmc[0].IDLyotStop[i] = -1;
        }

        piaacmc[0].fpmaskradld = fpmradld; // to compute prolate spheroidal function
        piaacmc[0].fpmarraysize = 2048;

        piaacmc[0].fpmRad = 100.0e-6; // focal plane radius [m]
        piaacmc[0].NBrings = 4; // number of rings in focal plane mask
        piaacmc[0].fpmminsag = -1.0e-5;
        piaacmc[0].fpmmaxsag = 1.0e-5;
        piaacmc[0].fpmsagreg_coeff = 1.0;
        piaacmc[0].fpmsagreg_alpha = 1.0;
        piaacmc[0].NBringCentCone = 0; // central cone
        piaacmc[0].fpmCentConeZ = piaacmc[0].fpmmaxsag;
        piaacmc[0].fpmOuterConeZ = piaacmc[0].fpmmaxsag;
        piaacmc[0].fpmOuterConeRadld = 80.0;
        piaacmc[0].fpmmaterial_code = 0;  // 0: mirror
        piaacmc[0].fpmaskamptransm = 1.0;


        if((fp = fopen("conf/conf_peakPSF.txt", "r"))!=NULL)
        {
            ret = fscanf(fp, "%f", &piaacmc[0].peakPSF);
            fclose(fp);
        }
        else
            piaacmc[0].peakPSF = -1.0;



        piaacmc[0].PIAAmode = 1;
        if((IDv=variable_ID("PIAACMC_PIAAmode"))!=-1)
            piaacmc[0].PIAAmode = (int) (data.variable[IDv].value.f+0.01);

        piaacmc[0].PIAAcoeff = 1.0;
        if((IDv=variable_ID("PIAACMC_PIAAcoeff"))!=-1)
            piaacmc[0].PIAAcoeff = data.variable[IDv].value.f;



        if(piaacmc[0].PIAAmode == 0)
        {
            piaacmc[0].invPIAAmode = 0;
        }
        else
        {
            piaacmc[0].invPIAAmode = 1;
            if((IDv=variable_ID("PIAACMC_invPIAAmode"))!=-1)
                piaacmc[0].invPIAAmode = (long) (data.variable[IDv].value.f+0.001);
        }



        
                sprintf(fname, "%s/conf_fpmmaterial_name.txt", piaacmcconfdir );
                if( (fp = fopen(fname, "r")) != NULL)
                {
                    ret = fscanf(fp, "%s", name);
                    strcpy(piaacmc[0].fpmmaterial_name, name);
                    printf("Reading %s   piaacmc[0].fpmmaterial_name : %s\n", fname, piaacmc[0].fpmmaterial_name);
                    fclose(fp);
                }
                else
                {
                    sprintf(piaacmc[0].fpmmaterial_name, "Mirror");
                    sprintf(fname, "%s/conf_fpmmaterial_name.txt", piaacmcconfdir);
                    printf("Writing %s   piaacmc[0].fpmmaterial_name : %s\n", fname, piaacmc[0].fpmmaterial_name);
                    if((fp=fopen(fname,"w"))!=NULL)
                    {
                        fprintf(fp, "%s\n", piaacmc[0].fpmmaterial_name);
                        fclose(fp);
                    }
                    else
                    {
                        printf("ERROR: cannot create file \"%s\"\n", fname);
                        exit(0);
                    }
                }

                printf("piaacmc[0].fpmmaterial_name : %s\n", piaacmc[0].fpmmaterial_name);
                piaacmc[0].fpmmaterial_code = OPTICSMATERIALS_code(piaacmc[0].fpmmaterial_name);

                sprintf(fname, "%s/conf_fpmmaterial_code.txt", piaacmcconfdir);
                if((fp=fopen(fname,"w"))!=NULL)
                {
                    fprintf(fp, "%d\n", piaacmc[0].fpmmaterial_code);
                    fclose(fp);
                }
                else
                {
                    printf("ERROR: cannot create file \"%s\"\n", fname);
                    exit(0);
                }
        



        if((IDv=variable_ID("PIAACMC_beamrad"))!=-1)
            piaacmc[0].beamrad = data.variable[IDv].value.f; // beam physical radius

        if((IDv=variable_ID("PIAACMC_Fratio"))!=-1)
            piaacmc[0].Fratio = data.variable[IDv].value.f; // Focal ratio
        if((IDv=variable_ID("PIAACMC_r0lim"))!=-1)
            piaacmc[0].r0lim = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_r1lim"))!=-1)
            piaacmc[0].r1lim = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_PIAAsep"))!=-1)
            piaacmc[0].PIAAsep = data.variable[IDv].value.f; // piaa separation
        if((IDv=variable_ID("PIAACMC_PIAA0pos"))!=-1)
            piaacmc[0].PIAA0pos = data.variable[IDv].value.f; // piaa elem 0 position



        if((IDv=variable_ID("PIAACMC_prePIAA0maskpos"))!=-1)
            piaacmc[0].prePIAA0maskpos = data.variable[IDv].value.f; // pre piaa elem 0 mask position
        if((IDv=variable_ID("PIAACMC_postPIAA0maskpos"))!=-1)
            piaacmc[0].postPIAA0maskpos = data.variable[IDv].value.f; // post piaa elem 0 mask position



        piaacmc[0].LyotZmin = -3.0;
        if((IDv=variable_ID("PIAACMC_LyotZmin"))!=-1)
            piaacmc[0].LyotZmin = data.variable[IDv].value.f;
        piaacmc[0].LyotZmax = 3.0;
        if((IDv=variable_ID("PIAACMC_LyotZmax"))!=-1)
            piaacmc[0].LyotZmax = data.variable[IDv].value.f;

        piaacmc[0].pupoutmaskrad = 0.95;
        if((IDv=variable_ID("PIAACMC_pupoutmaskrad"))!=-1)
            piaacmc[0].pupoutmaskrad = data.variable[IDv].value.f;




        if((IDv=variable_ID("PIAACMC_piaaNBCmodesmax"))!=-1)
            piaacmc[0].piaaNBCmodesmax = (long) (data.variable[IDv].value.f +0.01); // max number of Cosine terms
        if((IDv=variable_ID("PIAACMC_piaaCPAmax"))!=-1)
            piaacmc[0].piaaCPAmax = data.variable[IDv].value.f; // max CPA for PIAA shapes tuning


        piaacmc[0].NBLyotStop = 1;
        if(piaacmc[0].PIAAmode == 0)
        {
            piaacmc[0].NBLyotStop = 1;
        }
        else
        {
            if((IDv=variable_ID("PIAACMC_nblstop"))!=-1)
                piaacmc[0].NBLyotStop = (long) data.variable[IDv].value.f+0.01;
        }




        if((IDv=variable_ID("PIAACMC_lambda"))!=-1)
            piaacmc[0].lambda = 1.0e-9*data.variable[IDv].value.f; // central wavelength [m]
        //             printf("lambda = %g\n", piaacmc[0].lambda);

        if((IDv=variable_ID("PIAACMC_lambdaB"))!=-1)
            piaacmc[0].lambdaB = data.variable[IDv].value.f; // spectral bandwidth [%]

        LAMBDASTART = piaacmc[0].lambda * (1.0 - 0.005*piaacmc[0].lambdaB);
        LAMBDAEND = piaacmc[0].lambda * (1.0 + 0.005*piaacmc[0].lambdaB);


        if((IDv=variable_ID("PIAACMC_nblambda"))!=-1)
            piaacmc[0].nblambda = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_NBrings"))!=-1)
            piaacmc[0].NBrings = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_fpmminsag"))!=-1)
            piaacmc[0].fpmminsag = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_fpmmaxsag"))!=-1)
            piaacmc[0].fpmmaxsag = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_fpmsagreg_coeff"))!=-1)
            piaacmc[0].fpmsagreg_coeff = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_fpmsagreg_alpha"))!=-1)
            piaacmc[0].fpmsagreg_alpha = data.variable[IDv].value.f;
         
         
         
            
        if((IDv=variable_ID("PIAACMC_NBringCentCone"))!=-1)
            piaacmc[0].NBringCentCone = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_fpmCentConeZ"))!=-1)
            piaacmc[0].fpmCentConeZ = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_fpmOuterConeZ"))!=-1)
            piaacmc[0].fpmOuterConeZ = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_fpmOuterConeRadld"))!=-1)
            piaacmc[0].fpmOuterConeRadld = data.variable[IDv].value.f;
        piaacmc[0].fpmOuterConeRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio*piaacmc[0].fpmOuterConeRadld;  // [l/D] radius

        if((IDv=variable_ID("PIAACMC_size"))!=-1)
            piaacmc[0].size = (long) (data.variable[IDv].value.f+0.01);

        if((IDv=variable_ID("PIAACMC_pixscale"))!=-1)
            piaacmc[0].pixscale = data.variable[IDv].value.f;


        if(piaacmctype==0) // idealized focal plane mask
        {
            FORCE_CREATE_fpmzt = 1; // force making the focal plane mask
            piaacmc[0].NBrings = 1;
            piaacmc[0].NBringCentCone = 0;
            piaacmc[0].fpmOuterConeZ = 0.0;
            piaacmc[0].fpmminsag = -1e-5;
            piaacmc[0].fpmmaxsag = 1e-5;
            piaacmc[0].fpmsagreg_coeff = 1.0;
            piaacmc[0].fpmsagreg_alpha = 1.0;
            piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio*fpmradld;  // [l/D] radius
            printf("Idealized focal plane mask  radius = %f l/D  = %g m    [lambda = %g - %g]\n", fpmradld, piaacmc[0].fpmRad, LAMBDASTART, LAMBDAEND);
        }
        else
        {
            if(PIAACMC_MASKRADLD<0.2) // not initialized
                PIAACMC_MASKRADLD = 0.1*( (long) (10.0*1.2*fpmradld)); // 1.2x nominal radius, rounded to nearest 0.1 l/D

            piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio * PIAACMC_MASKRADLD;
            printf("Physical focal plane mask - rad = %f l/D -> %g    [lambda = %g - %g]\n", PIAACMC_MASKRADLD, piaacmc[0].fpmRad, LAMBDASTART, LAMBDAEND);
        }

        piaacmc[0].CmodesID = -1; // Cosine radial mode
        piaacmc[0].FmodesID = -1; // Fourier 2D modes
        piaacmc[0].piaa0CmodesID = -1;
        piaacmc[0].piaa0FmodesID = -1;
        piaacmc[0].piaa1CmodesID = -1;
        piaacmc[0].piaa1FmodesID = -1;
        piaacmc[0].zonezID = -1;  // focm zone material thickness, double precision image
        piaacmc[0].zoneaID = -1;  // focm zone amplitude transmission, double precision image
    }

    piaacmc[0].fpmCentConeRad = piaacmc[0].fpmRad*piaacmc[0].NBringCentCone/piaacmc[0].NBrings;

    printf("fpmRad = %g m\n", 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio*fpmradld);
    printf("factor = %f   (%ld / %ld)\n", 1.0*piaacmc[0].NBringCentCone/piaacmc[0].NBrings, piaacmc[0].NBringCentCone, piaacmc[0].NBrings);
    printf("fpmCentConeRad =  %g\n", (0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio*fpmradld)*piaacmc[0].NBringCentCone/piaacmc[0].NBrings);




    if(load==1)
    {
        printf("Loading PIAACMC configuration\n");
        fflush(stdout);
        sprintf(command, "mkdir -p %s", piaacmcconfdir);
        r = system(command);
        loaded = PIAAsimul_loadpiaacmcconf(piaacmcconfdir);
        if(loaded==0)
        {
            printf("Saving default configuration\n");
            fflush(stdout);
            saveconf = 1;
        }
    }





    sprintf(fname, "%s/conf_PIAAmaterial_name.txt", piaacmcconfdir );
    if( (fp = fopen(fname, "r")) != NULL)
    {
        ret = fscanf(fp, "%s", name);
        strcpy(piaacmc[0].PIAAmaterial_name, name);
        fclose(fp);
    }
    else
    {
        sprintf(piaacmc[0].PIAAmaterial_name, "Mirror");
        sprintf(fname, "%s/conf_PIAAmaterial_name.txt", piaacmcconfdir);
        if((fp=fopen(fname,"w"))!=NULL)
        {
            fprintf(fp, "%s\n", piaacmc[0].PIAAmaterial_name);
            fclose(fp);
        }
        else
        {
            printf("ERROR: cannot create file \"%s\"\n", fname);
            exit(0);
        }
    }

    printf("piaacmc[0].PIAAmaterial_name = %s\n", piaacmc[0].PIAAmaterial_name);
    piaacmc[0].PIAAmaterial_code = OPTICSMATERIALS_code(piaacmc[0].PIAAmaterial_name);

    sprintf(fname, "%s/conf_PIAAmaterial_code.txt", piaacmcconfdir);
    if((fp=fopen(fname,"w"))!=NULL)
    {
        fprintf(fp, "%d\n", piaacmc[0].PIAAmaterial_code);
        fclose(fp);
    }
    else
    {
        printf("ERROR: cannot create file \"%s\"\n", fname);
        exit(0);
    }



    printf("lambda = %g\n", piaacmc[0].lambda);
    printf("LAMBDASTART = %g\n", LAMBDASTART);
    printf("LAMBDAEND = %g\n", LAMBDAEND);


    for(k=0; k<piaacmc[0].nblambda; k++)
        piaacmc[0].lambdaarray[k] = LAMBDASTART + (0.5+k)*(LAMBDAEND-LAMBDASTART)/piaacmc[0].nblambda;





    // create modes for aspheric optical surfaces description
    beamradpix = piaacmc[0].beamrad/piaacmc[0].pixscale;
    size = piaacmc[0].size;

    printf("BEAM RADIUS :  %f / %f =  %f pix, size = %ld\n", piaacmc[0].beamrad, piaacmc[0].pixscale, beamradpix, size);
    fflush(stdout);

    // x, y, r and PA coordinates in beam (for convenience & speed)
    IDx = create_2Dimage_ID("xcoord", size, size);
    IDy = create_2Dimage_ID("ycoord", size, size);
    IDr = create_2Dimage_ID("rcoord", size, size);
    IDPA = create_2Dimage_ID("PAcoord", size, size);
    printf("pre-computing x, y, r, and PA\n");
    fflush(stdout);
    list_image_ID();

    for(ii=0; ii<size; ii++)
    {
        for(jj=0; jj<size; jj++)
        {
            x = (1.0*ii-0.5*size)/beamradpix;
            y = (1.0*jj-0.5*size)/beamradpix;
            data.image[IDx].array.F[jj*size+ii] = x;
            data.image[IDy].array.F[jj*size+ii] = y;
            data.image[IDr].array.F[jj*size+ii] = sqrt(x*x+y*y);
            data.image[IDPA].array.F[jj*size+ii] = atan2(y,x);
        }
    }

    // ==================== CREATE DMs ===============
    printf("%d DM(s)\n", piaacmc[0].nbDM);
    fflush(stdout);
    for(iDM=0; iDM<piaacmc[0].nbDM; iDM++)
    {
        printf("DM # %ld\n",iDM);
        sprintf(fname, "wfcDM%ld", iDM);
        piaacmc[0].ID_DM[iDM] = image_ID(fname);  // DM image identifier - to be updated later
        printf("ID = %ld", piaacmc[0].ID_DM[iDM]);

        if(piaacmc[0].ID_DM[iDM] == -1)
        {
            read_sharedmem_image(fname);
            piaacmc[0].ID_DM[iDM] = image_ID(fname);
        }

        if(piaacmc[0].ID_DM[iDM] == -1)
        {
            sizearray = (long*) malloc(sizeof(long)*2);
            sizearray[0] = size;
            sizearray[1] = size;
            piaacmc[0].ID_DM[iDM] = create_image_ID(fname, 2, sizearray, FLOAT, 1, 0);
            free(sizearray);
        }
    }

    // ==================== CREATE MODES USED TO FIT AND DESCRIBE PIAA SHAPES ===============
    printf("Creating / loading Cmodes and Fmodes ...\n");
    fflush(stdout);


    CREATE_Cmodes = 0;
    //   sprintf(fname, "%s/Cmodes.fits", piaacmcconfdir);
    sprintf(fname, "Cmodes_%ld.fits", piaacmc[0].size);
    if(FORCE_CREATE_Cmodes==0)
    {
        piaacmc[0].CmodesID = image_ID("Cmodes");
        if(piaacmc[0].CmodesID==-1)
            piaacmc[0].CmodesID = load_fits(fname, "Cmodes", 0);
        if(piaacmc[0].CmodesID==-1)
            CREATE_Cmodes = 1;
    }
    else
        CREATE_Cmodes = 1;
    if(CREATE_Cmodes == 1)
    {
        if(piaacmc[0].CmodesID!=-1)
            delete_image_ID("Cmodes");
        Cmsize = (long) (beamradpix*4);
        if(Cmsize>size)
			Cmsize = size;
        printf("beamradpix = %f -> Cmsize = %ld\n", beamradpix, Cmsize);
        // make sure Cmsize if even
        if (Cmsize%2 == 1)
            Cmsize++;
        piaacmc[0].Cmsize = Cmsize;
        linopt_imtools_makeCosRadModes("Cmodes", Cmsize, piaacmc[0].piaaNBCmodesmax, ApoFitCosFact*beamradpix, 2.0);
        piaacmc[0].CmodesID = image_ID("Cmodes");
        save_fits("Cmodes", fname);
    }
    piaacmc[0].NBCmodes = data.image[piaacmc[0].CmodesID].md[0].size[2];
    piaacmc[0].Cmsize = data.image[piaacmc[0].CmodesID].md[0].size[0];



    CREATE_Fmodes = 0;
    //    sprintf(fname, "%s/Fmodes.fits", piaacmcconfdir);
    sprintf(fname, "Fmodes_%ld.fits", piaacmc[0].size);
    if(FORCE_CREATE_Fmodes == 0)
    {
        piaacmc[0].FmodesID = image_ID("Fmodes");
        if(piaacmc[0].FmodesID==-1)
            piaacmc[0].FmodesID = load_fits(fname, "Fmodes", 0);
        if(piaacmc[0].FmodesID==-1)
            CREATE_Fmodes = 1;
    }
    else
        CREATE_Fmodes = 1;
    if(CREATE_Fmodes == 1)
    {
        Fmsize = (long) (beamradpix*4);
         if(Fmsize>size)
			Fmsize = size;
        piaacmc[0].Fmsize = Fmsize;
        linopt_imtools_makeCPAmodes("Fmodes",  Fmsize, piaacmc[0].piaaCPAmax, 0.8, beamradpix, 2.0, 1);
        piaacmc[0].FmodesID = image_ID("Fmodes");
        save_fits("Fmodes", fname);
		save_fits("cpamodesfreq", "!cpamodesfreq.fits");
        sprintf(command, "mv ModesExpr_CPA.txt %s/", piaacmcconfdir);
        
        r = system(command);
    }
    piaacmc[0].NBFmodes = data.image[piaacmc[0].FmodesID].md[0].size[2];
    piaacmc[0].Fmsize = data.image[piaacmc[0].FmodesID].md[0].size[0];

    printf("DONE Creating / loading Cmodes and Fmodes\n");
    fflush(stdout);









    // =================== IMPORT / CREATE PIAA SHAPES =====================
    sprintf(command, "mkdir -p %s/piaaref/", piaacmcconfdir);
    ret = system(command);

    if(piaacmc[0].PIAAmode == 1)
    {
        piaacmc[0].piaa0CmodesID = image_ID("piaa0Cmodescoeff");
        piaacmc[0].piaa0FmodesID = image_ID("piaa0Fmodescoeff");
        piaacmc[0].piaa1CmodesID = image_ID("piaa1Cmodescoeff");
        piaacmc[0].piaa1FmodesID = image_ID("piaa1Fmodescoeff");

        sprintf(command, "mkdir -p %s/piaaref/", piaacmcconfdir);
        ret = system(command);

        if((piaacmc[0].piaa0CmodesID==-1)||( piaacmc[0].piaa0FmodesID==-1)||(piaacmc[0].piaa1CmodesID==-1)||( piaacmc[0].piaa1FmodesID==-1))
        {
            sprintf(fname, "%s/piaaref/piaa0Cmodes.fits", piaacmcconfdir);
            piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff", 1);

            sprintf(fname, "%s/piaaref/piaa0Fmodes.fits", piaacmcconfdir);
            piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff", 1);

            sprintf(fname, "%s/piaaref/piaa1Cmodes.fits", piaacmcconfdir);
            piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff", 1);

            sprintf(fname, "%s/piaaref/piaa1Fmodes.fits", piaacmcconfdir);
            piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff", 1);

            sprintf(fname, "%s/piaaref/APLCmaskCtransm.txt", piaacmcconfdir);
            fp = fopen(fname, "r");
            if(fp!=NULL)
            {
                ret = fscanf(fp, "%f", &tmpf);
                piaacmc[0].fpmaskamptransm = tmpf;
                fclose(fp);
            }
        }
    }

    if((piaacmc[0].piaa0CmodesID==-1)||( piaacmc[0].piaa0FmodesID==-1)||(piaacmc[0].piaa1CmodesID==-1)||( piaacmc[0].piaa1FmodesID==-1))
    {
        sprintf(fname, "%s/apo2Drad.fits", piaacmcconfdir);
        if(load_fits(fname, "apo2Drad", 1)==-1)  // CREATE APODIZATION
        {
            sprintf(command, "cp %s/piaaref/apo2Drad.fits %s/apo2Drad.fits", piaacmcconfdir, piaacmcconfdir);
            ret = system(command);

            sprintf(fname, "%s/apo2Drad.fits", piaacmcconfdir);
            if(load_fits(fname, "apo2Drad", 1)==-1)
            {

                printf("Creating 2D apodization for idealized circular monochromatic PIAACMC\n");
                fflush(stdout);

                // first iteration: half size image, 2x zoom, without pupil mask
                IDv1 = create_variable_ID("DFTZFACTOR", 2);
                IDv2 = create_variable_ID("PNBITER", 15);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix*0.5, piaacmc[0].centObs1, "apotmp1", size/2, "NULLim");

                // expand solution to full size
                basic_resizeim("apotmp1", "apostart", size, size);
                delete_image_ID("apotmp1");

                // full size, 4x zoom
                IDv1 = create_variable_ID("DFTZFACTOR", 4);
                IDv2 = create_variable_ID("PNBITER", 5);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size, "pupmaskim");

                // full size, 8x zoom
                chname_image_ID("apo", "apostart");
                IDv1 = create_variable_ID("DFTZFACTOR", 8);
                IDv2 = create_variable_ID("PNBITER", 5);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size, "pupmaskim");


                // full size, 16x zoom
                chname_image_ID("apo", "apostart");
                IDv1 = create_variable_ID("DFTZFACTOR", 16);
                IDv2 = create_variable_ID("PNBITER", 10);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size, "pupmaskim");


                chname_image_ID("apo", "apo2Drad");
                sprintf(fname, "!%s/apo2Drad.fits", piaacmcconfdir);
                save_fits("apo2Drad", fname);

                if(piaacmc[0].PIAAmode == 1)
                {
                    sprintf(fname, "!%s/piaaref/apo2Drad.fits", piaacmcconfdir);
                    save_fits("apo2Drad", fname);
                }



                if((piaacmctype==0)&&(loaded==0)) // idealized focal plane mask
                {
                    piaacmc[0].fpmaskamptransm =  -data.variable[variable_ID("APLCmaskCtransm")].value.f;
                    printf("FOCAL PLANE MASK TRANSM = %f\n", piaacmc[0].fpmaskamptransm);
                    printf("Saving default configuration\n");
                    fflush(stdout);
                    saveconf = 1;

                    sprintf(fname, "%s/piaaref/APLCmaskCtransm.txt", piaacmcconfdir);
                    fp = fopen(fname, "w");
                    fprintf(fp, "%.20f\n", piaacmc[0].fpmaskamptransm);
                    fclose(fp);
                }

            }


        }


        // split apodization in conventional pupil apodizer (apoCPA) and PIAA apodization (apo2Drad_PIAA)
        IDapo = image_ID("apo2Drad");
        xsize = data.image[IDapo].md[0].size[0];
        ysize = data.image[IDapo].md[0].size[1];
        IDapo_PIAA = create_2Dimage_ID("apo2Drad_PIAA", xsize, ysize);
        IDapo_CPA = create_2Dimage_ID("apo2Drad_CPA", xsize, ysize);

        if(piaacmc[0].PIAAmode==0)
        {
            for(ii=0; ii<xsize*ysize; ii++)
            {
                data.image[IDapo_PIAA].array.F[ii] = 1.0;
                data.image[IDapo_CPA].array.F[ii] = data.image[IDapo].array.F[ii];
            }
        }
        else
        {
            for(ii=0; ii<xsize*ysize; ii++)
            {
                coeff = piaacmc[0].PIAAcoeff; // fraction of apodization done by PIAA - between 0 and 1
                data.image[IDapo_PIAA].array.F[ii] = pow(data.image[IDapo].array.F[ii], coeff);
                data.image[IDapo_CPA].array.F[ii] = pow(data.image[IDapo].array.F[ii], 1.0-coeff);
            }
        }

        copy_image_ID("apo2Drad_CPA", "prePIAA0mask", 0);
        save_fits("prePIAA0mask", "!test_prePIAA0mask.fits");


        // load PIAA apodization profile and fit it a series of cosines
        PIAACMCsimul_load2DRadialApodization("apo2Drad_PIAA", beamradpix, "outApofit");

        // compute radial PIAA sag -> <piaacmcconfdir>/PIAA_Mshapes.txt
        PIAACMCsimul_init_geomPIAA_rad("outApofit");


        // make 2D sag maps
        sprintf(fname, "%s/PIAA_Mshapes.txt", piaacmcconfdir);
        PIAACMCsimul_mkPIAAMshapes_from_RadSag(fname, "piaam0z", "piaam1z");


        if(PIAACMC_save==1)
        {
            sprintf(fname, "!%s/piaam0z.fits", piaacmcconfdir);
            save_fits("piaam0z", fname);

            sprintf(fname, "!%s/piaam1z.fits", piaacmcconfdir);
            save_fits("piaam1z", fname);
        }

        // crop piaam0z and piaam1z to Cmodes size
        ID0 = image_ID("Cmodes");
        size0 = data.image[ID0].md[0].size[0];
        ID1 = create_2Dimage_ID("piaa0zcrop", size0, size0);
        ID = image_ID("piaam0z");
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID1].array.F[jj*size0+ii] = data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)];

        ID1 = create_2Dimage_ID("piaa1zcrop", size0, size0);
        ID = image_ID("piaam1z");
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID1].array.F[jj*size0+ii] = data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)];

        make_disk("maskd", size0, size0, 0.5*size0, 0.5*size0, beamradpix);
        make_2Dgridpix("gridpix", size0, size0, 1, 1, 0, 0);
        arith_image_mult("maskd", "gridpix", "maskfit");

        //sprintf(fname, "!%s/maskfit.fits", piaacmcconfdir);
        //save_fits("maskfit", fname);

        printf("--------- FITTING COSINE MODES ---------\n");
        fflush(stdout);

        linopt_imtools_image_fitModes("piaa0zcrop", "Cmodes", "maskfit", 1.0e-6, "piaa0Cmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa0Cmodes.dat", piaacmcconfdir, piaacmcconfdir);
        ret = system(command);


        //      sprintf(fname, "!%s/piaa0Cmodescoeff.fits", piaacmcconfdir);
        //      save_fits("piaa0Cmodescoeff", fname);

        linopt_imtools_image_fitModes("piaa1zcrop", "Cmodes", "maskfit", 1.0e-6, "piaa1Cmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa1Cmodes.dat", piaacmcconfdir, piaacmcconfdir);
        ret = system(command);

        //      sprintf(fname, "!%s/piaa1Cmodescoeff.fits", piaacmcconfdir);
        //      save_fits("piaa1Cmodescoeff", fname);

        linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");

        if(PIAACMC_save==1)
        {
            sprintf(fname, "!%s/piaa0Cz.fits", piaacmcconfdir);
            save_fits("piaa0Cz", fname);
        }

        linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");

        if(PIAACMC_save==1)
        {
            sprintf(fname, "!%s/piaa1Cz.fits", piaacmcconfdir);
            save_fits("piaa1Cz", fname);
        }


        ID0 = image_ID("piaa0Cz");
        size0 = data.image[ID0].md[0].size[0];
        ID1 = image_ID("piaam0z");
        ID = create_2Dimage_ID("piaa0Cres", size0, size0);
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[jj*size0+ii] = data.image[ID1].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)]-data.image[ID0].array.F[jj*size0+ii];
        if(PIAACMC_save==1)
        {
            sprintf(fname, "!%s/piaa0Cres.fits", piaacmcconfdir);
            save_fits("piaa0Cres", fname);
        }




        ID0 = image_ID("piaa1Cz");
        size0 = data.image[ID0].md[0].size[0];
        ID1 = image_ID("piaam1z");
        ID = create_2Dimage_ID("piaa1Cres", size0, size0);
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[jj*size0+ii] = data.image[ID1].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)]-data.image[ID0].array.F[jj*size0+ii];

        if(PIAACMC_save==1)
        {   sprintf(fname, "!%s/piaa1Cres.fits", piaacmcconfdir);
            save_fits("piaa1Cres", fname);
        }



        printf("--------- FITTING FOURIER MODES ---------\n");
        fflush(stdout);

        linopt_imtools_image_fitModes("piaa0Cres", "Fmodes", "maskfit", 0.01, "piaa0Fmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa0Fmodes.dat", piaacmcconfdir, piaacmcconfdir);
        ret = system(command);

        linopt_imtools_image_fitModes("piaa1Cres", "Fmodes", "maskfit", 0.01, "piaa1Fmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa1Fmodes.dat", piaacmcconfdir, piaacmcconfdir);
        ret = system(command);


        // save_fits("piaa1Fmodescoeff", "!piaa1Fmodescoeff.fits");

        //linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
        //   save_fits("piaa0Fz", "!piaa0Fz.fits");
        //arith_image_sub("piaa0Cres", "piaa0Fz", "piaa0CFres");
        //save_fits("piaa0CFres", "!piaa0CFres.fits");
        delete_image_ID("piaa0zcrop");

        //linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
        //save_fits("piaa1Fz", "!piaa1Fz.fits");
        //arith_image_sub("piaa1Cres", "piaa1Fz", "piaa1CFres");
        //save_fits("piaa1CFres", "!piaa1CFres.fits");
        delete_image_ID("piaa1zcrop");

        delete_image_ID("maskfit");



        piaacmc[0].piaa0CmodesID = image_ID("piaa0Cmodescoeff");
        piaacmc[0].piaa0FmodesID = image_ID("piaa0Fmodescoeff");
        piaacmc[0].piaa1CmodesID = image_ID("piaa1Cmodescoeff");
        piaacmc[0].piaa1FmodesID = image_ID("piaa1Fmodescoeff");



        sprintf(fname, "!%s/piaaref/piaa0Cmodes.fits", piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa0CmodesID].name, fname);

        sprintf(fname, "!%s/piaaref/piaa0Fmodes.fits", piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa0FmodesID].name, fname);

        sprintf(fname, "!%s/piaaref/piaa1Cmodes.fits", piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa1CmodesID].name, fname);

        sprintf(fname, "!%s/piaaref/piaa1Fmodes.fits", piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa1FmodesID].name, fname);

        sprintf(command, "cp %s/piaaref/* %s/", piaacmcconfdir, piaacmcconfdir);
        ret = system(command);
    }


    // ============ MAKE FOCAL PLANE MASK ===============

    /*    CREATE_fpmzmap = 0;
        if(FORCE_CREATE_fpmzmap == 0)
        {
            if(image_ID("fpmzmap")==-1)
                {
                    sprintf(fname, "%s/fpmzmap%d_%03ld_%03ld.fits", piaacmcconfdir, PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone);
                    load_fits(fname, "fpmzmap", 1);
                    if(image_ID("fpmzmap")==-1)
                        CREATE_fpmzmap = 1;
                }
        }
        else
            CREATE_fpmzmap = 1;

        if(CREATE_fpmzmap == 1)
        {
            if(image_ID("fpmzmap")!=-1)
                delete_image_ID("fpmzmap");
            PIAACMCsimul_mkFPM_zonemap("fpmzmap");
            sprintf(fname, "!%s/fpmzmap%d_%03ld_%03ld.fits", piaacmcconfdir, PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone);
            save_fits("fpmzmap", fname);
        }
    */

    if(image_ID("fpmzmap")==-1)
    {
        printf("Make zonemap ...\n");
        fflush(stdout);
        PIAACMCsimul_mkFPM_zonemap("fpmzmap");
    }
    else
    {
        printf("zonemap already exists\n");
        fflush(stdout);
    }

    //    sprintf(fname, "!%s/fpmzmap.fits", piaacmcconfdir);
    //    save_fits("fpmzmap", fname);



    /* sprintf(fname, "!%s/fpmzmap.fits", piaacmcconfdir);
     save_fits("fpmzmap", fname);
     exit(0);*/


    // zones thickness

    CREATE_fpmzt = 0;
    if(FORCE_CREATE_fpmzt == 0)
    {
        piaacmc[0].zonezID = image_ID("fpmzt");
        if(piaacmc[0].zonezID == -1)
        {            
            sprintf(fname, "%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

            printf("LOADING FILE NAME : \"%s\"  -  %ld %d \n", fname, piaacmctype, loaded);

            piaacmc[0].zonezID = load_fits(fname, "fpmzt", 1);
            if(piaacmc[0].zonezID == -1)
                CREATE_fpmzt = 1;
        }
    }
    else
        CREATE_fpmzt = 1;



    if(CREATE_fpmzt == 1)
    {
        printf("Creating fpmzt, saving as fpm_zonez.fits - %ld %d\n", piaacmctype, loaded);
        piaacmc[0].zonezID = image_ID("fpmzt");
        if(piaacmc[0].zonezID!=-1)
            delete_image_ID("fpmzt");

        piaacmc[0].zonezID = create_2Dimage_ID_double("fpmzt", piaacmc[0].focmNBzone, 1);
        t = 1.0e-9;

        if(piaacmctype==0) // idealized focal plane mask
        {
            printf("IDEALIZED FOCAL PLANE MASK\n");
            fflush(stdout);
            //          exit(0);

            // measure dpha/dt
            t0 = 1.0e-8;
            pha0 = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, t0, 0.5*(LAMBDASTART+LAMBDAEND));
            // set t to get PI phase
            t = (M_PI/pha0)*t0;
            printf("t = %g m (%lf %g) -> %g %g\n", t, pha0, t0, OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, t, 0.5*(LAMBDASTART+LAMBDAEND)), 0.5*(LAMBDASTART+LAMBDAEND));

            printf(" -- lambda = %g\n", piaacmc[0].lambda);
            printf(" -- lambdaB = %g\n", piaacmc[0].lambdaB);
            printf(" -- LAMBDASTART = %g\n", LAMBDASTART);
            printf(" -- LAMBDAEND = %g\n", LAMBDAEND);
        }
        else
        {
            printf("CREATING EXAMPLE FOCAL PLANE MASK  %ld %d\n", piaacmctype, loaded);
            fflush(stdout);
            //                exit(0);
        }

        for(ii=0; ii<piaacmc[0].focmNBzone; ii++)
            data.image[piaacmc[0].zonezID].array.D[ii] = t;

		sprintf(fname, "!%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

        printf("Writing %s\n", fname);
        save_fits("fpmzt", fname);
    }


    // zones transmission amplitude

    printf("CREATE_fpmza = %d\n", CREATE_fpmza);
    if(FORCE_CREATE_fpmza == 0)
    {
        piaacmc[0].zoneaID = image_ID("fpmza");
        if(piaacmc[0].zoneaID == -1)
        {
            sprintf(fname, "%s/fpm_zonea_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

            printf("LOADING FILE NAME : \"%s\"\n", fname);
            piaacmc[0].zoneaID = load_fits(fname, "fpmza", 1);

            if(piaacmc[0].zoneaID == -1)
                CREATE_fpmza = 1;
        }
    }
    else
        CREATE_fpmza = 1;



    if(CREATE_fpmza == 1)
    {
        if(piaacmc[0].zoneaID != -1)
            delete_image_ID("fpmza");
        piaacmc[0].zoneaID = create_2Dimage_ID_double("fpmza", piaacmc[0].focmNBzone, 1);

        if(PIAACMC_MASKRADLD>0.2) // physical mask
        {
            printf("PHYSICAL MASK ... %ld zones\n", piaacmc[0].focmNBzone);
            for(ii=0; ii<piaacmc[0].focmNBzone; ii++)
                data.image[piaacmc[0].zoneaID].array.D[ii] = 1.0;
        }
        else // idealized mask
        {
            printf("IDEALIZED MASK ... %ld zones\n", piaacmc[0].focmNBzone);
            for(ii=0; ii<piaacmc[0].focmNBzone; ii++)
                data.image[piaacmc[0].zoneaID].array.D[ii] = piaacmc[0].fpmaskamptransm;
        }


        sprintf(fname, "!%s/fpm_zonea_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

        printf("Writing %s\n", fname);
        save_fits("fpmza", fname);
    }

    //   printf("%d piaacmc[0].fpmaskamptransm = %f       %lf\n", CREATE_fpmza, piaacmc[0].fpmaskamptransm, data.image[piaacmc[0].zoneaID].array.D[0]);
    //   sleep(10);


    // ============= MAKE LYOT STOPS =======================
    printf("LOADING/CREATING LYOT MASK  - %ld masks  (PIAAmode = %d, %ld x %ld)\n", piaacmc[0].NBLyotStop, piaacmc[0].PIAAmode, xsize, ysize);
    list_image_ID();
    size2 = size*size;

	
    if(piaacmc[0].PIAAmode == 1)
    {
        for(i=0; i<piaacmc[0].NBLyotStop; i++)
        {
            printf("LYOT MASK %ld\n", i);
            fflush(stdout);

            sprintf(fname, "%s/LyotStop%ld.fits", piaacmcconfdir, i);
            sprintf(name, "lyotstop%ld", i);

            piaacmc[0].IDLyotStop[i] = image_ID(name);
            if(piaacmc[0].IDLyotStop[i]==-1)
            {
                sprintf(fname, "!%s/LyotStop%ld.fits", piaacmcconfdir, i);
                if ((i == 0) && (piaacmc[0].NBLyotStop > 1))
                    piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, -0.01, 0.98);
                else if ((i == 1) && (piaacmc[0].NBLyotStop > 2))
                    piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, piaacmc[0].centObs1+0.02, 1.2);
                else
                    piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, piaacmc[0].centObs1+0.02, 0.98);

                save_fl_fits(name, fname);
            }
        }
    }
    else
    {
        i = 0;
        sprintf(fname, "!%s/LyotStop%ld.fits", piaacmcconfdir, i);
        sprintf(name, "lyotstop%ld", i);
        
        piaacmc[0].IDLyotStop[i] = image_ID(name);
        if(piaacmc[0].IDLyotStop[i]==-1)
			{
				piaacmc[0].IDLyotStop[i] = create_2Dimage_ID(name, xsize, ysize);
				ID = image_ID("pupmaskim");
				for(ii=0; ii<xsize*ysize; ii++)
					if(data.image[ID].array.F[ii] < 0.999)
						data.image[piaacmc[0].IDLyotStop[i]].array.F[ii] = 0.0;
					else
						data.image[piaacmc[0].IDLyotStop[i]].array.F[ii] = 1.0;
				save_fl_fits(name, fname);
			}
    }

    if(saveconf==1)
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);



    return(0);
}











///
/// mirror sag shapes:  piaam0z, piaam1z
///
/// piaa0Cmodescoeff -> piaa0Cz
/// piaa0Fmodescoeff -> piaa0Fz
/// piaa0Cz + piaa0Fz -> piaam0z
///

int PIAACMCsimul_makePIAAshapes(OPTPIAACMCDESIGN *design, long index)
{
    long ID, ID0, ID1;
    long size, size0, size1;
    long ii, jj;
    char fname[500];
    long IDpiaam0z, IDpiaam1z;
    long IDpiaar0zsag, IDpiaar1zsag;
    double ri, ri0, sag2opd_coeff;
    int mkpiaar0zsag, mkpiaar1zsag;
    double sag2opd_coeff0;
    long k;
    FILE *fpri;

    size = piaacmc[0].size;

    // ============ construct PIAA shapes from fitting coefficients ==================



	if(piaacmc[0].PIAAmode==1)
	{



    MAKE_PIAA0shape = 0;
    if(FORCE_MAKE_PIAA0shape == 0)
    {
        ID = image_ID("piaam0z");
        if(ID==-1)
            MAKE_PIAA0shape = 1;
    }
    else
        MAKE_PIAA0shape = 1;






    MAKE_PIAA1shape = 0;
    if(FORCE_MAKE_PIAA1shape == 0)
    {
        ID = image_ID("piaam1z");
        if(ID==-1)
            MAKE_PIAA1shape = 1;
    }
    else
        MAKE_PIAA1shape = 1;






    if(MAKE_PIAA0shape == 1)
    {
        // assemble piaa0z and piaa1z images
        ID0 = linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");
        ID1 = linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
        ID = image_ID("piaam0z");
        if(ID==-1)
            ID = create_2Dimage_ID("piaam0z", size, size);



        size0 = data.image[ID0].md[0].size[0];
        size1 = data.image[ID1].md[0].size[0];
        for(ii=0; ii<size*size; ii++)
            data.image[ID].array.F[ii] = 0.0;

        printf("========================== STEP 01a  %ld %ld %ld\n", ID, ID0, ID1);
		fflush(stdout);
		

        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)] += data.image[ID0].array.F[jj*size0+ii];
        for(ii=0; ii<size1; ii++)
            for(jj=0; jj<size1; jj++)
                data.image[ID].array.F[(jj+(size-size1)/2)*size+(ii+(size-size1)/2)] += data.image[ID1].array.F[jj*size1+ii];

        if(PIAACMC_save==1)
        {   sprintf(fname, "!%s/piaa0Cz.fits", piaacmcconfdir);
            save_fits("piaa0Cz", fname);

            sprintf(fname, "!%s/piaa0Fz.fits", piaacmcconfdir);
            save_fits("piaa0Fz", fname);

            sprintf(fname, "!%s/piaam0z.fits", piaacmcconfdir);
            save_fits("piaam0z", fname);
        }
        delete_image_ID("piaa0Cz");
        delete_image_ID("piaa0Fz");

        IDpiaam0z = ID;

        // make lense shapes if applicable
        if(design[index].PIAAmaterial_code != 0) // refractive PIAA
        {
            // if piaar0zsag does not exist or is wrong size, create it
            IDpiaar0zsag = image_ID("piaar0zsag");
            if(IDpiaar0zsag == -1)
                mkpiaar0zsag = 1;
            else
            {
                if((data.image[IDpiaar0zsag].md[0].size[0] != size)||(data.image[IDpiaar0zsag].md[0].size[1] != size)) //||(data.image[IDpiaar0zsag].md[0].size[2] != design[index].nblambda))
                {
                    delete_image_ID("piaar0zsag");
                    mkpiaar0zsag = 1;
                }
            }
            if(mkpiaar0zsag == 1)
                IDpiaar0zsag = create_2Dimage_ID("piaar0zsag", size, size); //, design[index].nblambda);



            sprintf(fname, "%s/ri_array.txt", piaacmcconfdir);
            if( (fpri=fopen(fname, "w")) == NULL)
            {
                printf("ERROR: cannot open file \"%s\"\n", fname);
                exit(0);
            }
            ri0 = OPTICSMATERIALS_n(design[index].PIAAmaterial_code, design[index].lambda); // refractive index at central lambda
            sag2opd_coeff0 = (ri0-1.0)/2.0;
            for(k=0; k<design[index].nblambda; k++)
            {
                // sag to OPD coeff
                ri = OPTICSMATERIALS_n(design[index].PIAAmaterial_code, piaacmc[0].lambdaarray[k]); // refractive index
                sag2opd_coeff = (ri-1.0)/2.0;
                fprintf(fpri, "%g %.16f %.16f %.16f %.16f\n", piaacmc[0].lambdaarray[k], ri, ri0, sag2opd_coeff, sag2opd_coeff/sag2opd_coeff0);
                //                for(ii=0; ii<size*size; ii++)
                //                  data.image[IDpiaar0zsag].array.F[k*size*size+ii] = data.image[IDpiaam0z].array.F[ii] * sag2opd_coeff/sag2opd_coeff0; //sag2opd_coeff * data.image[IDpiaam0z].array.F[ii] / sag2opd_coeff0;
            }
            fclose(fpri);

            for(ii=0; ii<size*size; ii++)
                data.image[IDpiaar0zsag].array.F[ii] = data.image[IDpiaam0z].array.F[ii] / sag2opd_coeff0;

            sprintf(fname, "!%s/piaar0zsag.fits", piaacmcconfdir);
            if(PIAACMC_save==1)
                save_fl_fits("piaar0zsag", fname);
            printf("Saved piaar0zsag to %s\n", fname);
        }

    }

	

    if(MAKE_PIAA1shape == 1)
    {
        ID0 = linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");
        ID1 = linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
        ID = image_ID("piaam1z");
        if(ID==-1)
            ID = create_2Dimage_ID("piaam1z", size, size);
        for(ii=0; ii<size*size; ii++)
            data.image[ID].array.F[ii] = 0.0;
        size0 = data.image[ID0].md[0].size[0];
        size1 = data.image[ID1].md[0].size[0];
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)] += data.image[ID0].array.F[jj*size0+ii];
        for(ii=0; ii<size1; ii++)
            for(jj=0; jj<size1; jj++)
                data.image[ID].array.F[(jj+(size-size1)/2)*size+(ii+(size-size1)/2)] += data.image[ID1].array.F[jj*size1+ii];

        if(PIAACMC_save==1)
        {   sprintf(fname, "!%s/piaa1Cz.fits", piaacmcconfdir);
            save_fits("piaa1Cz", fname);

            sprintf(fname, "!%s/piaa1Fz.fits", piaacmcconfdir);
            save_fits("piaa1Fz", fname);

            sprintf(fname, "!%s/piaam1z.fits", piaacmcconfdir);
            save_fits("piaam1z", fname);
        }
        delete_image_ID("piaa1Cz");
        delete_image_ID("piaa1Fz");

        IDpiaam1z = ID;


        // make lense shapes if applicable
        if(design[index].PIAAmaterial_code != 0) // refractive PIAA
        {

            // if piaar1zsag does not exist or is wrong size, create it
            IDpiaar1zsag = image_ID("piaar0zsag");
            if(IDpiaar1zsag == -1)
                mkpiaar1zsag = 1;
            else
            {
                if((data.image[IDpiaar1zsag].md[0].size[0] != size)||(data.image[IDpiaar1zsag].md[0].size[1] != size)) //||(data.image[IDpiaar1zsag].md[0].size[2] != design[index].nblambda))
                {
                    delete_image_ID("piaar1zsag");
                    mkpiaar1zsag = 1;
                }
            }
            if(mkpiaar1zsag == 1)
                IDpiaar1zsag = create_2Dimage_ID("piaar1zsag", size, size); //, design[index].nblambda);



            sprintf(fname, "%s/ri_array.txt", piaacmcconfdir);
            if( (fpri=fopen(fname, "w")) == NULL)
            {
                printf("ERROR: cannot open file \"%s\"\n", fname);
                exit(0);
            }
            ri0 = OPTICSMATERIALS_n(design[index].PIAAmaterial_code, design[index].lambda); // refractive index at central lambda
            sag2opd_coeff0 = (ri0-1.0)/2.0;
            for(k=0; k<piaacmc[0].nblambda; k++)
            {
                // sag to OPD coeff
                ri = OPTICSMATERIALS_n(design[index].PIAAmaterial_code, piaacmc[0].lambdaarray[k]); // refractive index
                sag2opd_coeff = (ri-1.0)/2.0;
                fprintf(fpri, "%g %.16f %.16f %.16f %.16f\n", piaacmc[0].lambdaarray[k], ri, ri0, sag2opd_coeff, sag2opd_coeff/sag2opd_coeff0);
                // for(ii=0; ii<size*size; ii++)
                //    data.image[IDpiaar1zsag].array.F[k*size*size+ii] = sag2opd_coeff * data.image[IDpiaam1z].array.F[ii] / sag2opd_coeff0;
            }
            fclose(fpri);

            for(ii=0; ii<size*size; ii++)
                data.image[IDpiaar1zsag].array.F[ii] = data.image[IDpiaam1z].array.F[ii] / sag2opd_coeff0;

            sprintf(fname, "!%s/piaar1zsag.fits", piaacmcconfdir);
            if(PIAACMC_save==1)   save_fl_fits("piaar1zsag", fname);
            printf("Saved piaar1zsag to %s\n", fname);

        }
    }


	}

    return 0;
}






///
/// returns average contrast in evaluation zone
///
/// source size = 1e-{sourcesize*0.1}, except if sourcesize = 0 (point source)
/// sourcesize is a 2-digit number ( 10 = 0.1 l/D, 20 = 0.01 l/D etc..)
///
/// extmode = 0 : 3 point sources, 120 apart on circle radius = source size
/// extmode = 1 : 6 point sources. 3 as above on circle radius 1/sqrt(2.5) + 3 on outer circle, radius 2/sqrt(2.5), 120 apart, clockled 60 deg off inner points
///
/// if opderrcube exists, include each slice as a WF mode  
///
/// PSF is held in shared memory by default
///

double PIAACMCsimul_computePSF(float xld, float yld, long startelem, long endelem, int savepsf, int sourcesize, int extmode, int outsave)
{
    FILE *fp;
    FILE *fpflux;
    double x, y;
    long IDa, IDp;
    long size;
    long nblambda;
    long size2;
    long ii, jj, k;
    long IDpiaa1z, IDpiaa2z;
    long elem;
    long kl;

    char fname_piaa1z[500];
    char fname_piaa2z[500];
    char fname_pupa0[500];
    char fname_pupp0[500];
    char fname[500];

    long ID;
    long index;

    double proplim = 1.0e-4;
    double total;

    long size0, size1;
    long Cmsize, Fmsize;


    // how to measure quality
    float focscale; // l/D per pix
    float scoringIWA = 1.5; 
    float scoringOWA = 20.0;
    float scoringOWAhr = 8.0;
    float scoringIWAx = -20.5;
    long IDsm;
    float r;

    double value;
    double avContrast;
    double peakcontrast;
    double tmpv;
    double value1;


    double dld;
    long nbelem;
    long ID1;
    long offset;
    long offset1;
    double normcoeff = 1.0;
    double pha;

    int ret;
    char command[1000];
    double rad1, rad2;

    float val;
    long IDv;

	long nbOPDerr = 0;
	long OPDmode;
	long IDopderrC = -1;
	long IDopderr;
	
	long imindex = 0;
	long NBimindex = 0;
	char imname[200];
	

    // size of one side of each image array
    size = piaacmc[0].size;
    // number of pixels in each image
    size2 = size*size;

	//printf("Loading (optional) OPDerr file\n");
	//fflush(stdout);
	// load an error if it exists
	IDopderrC = image_ID("OPDerrC");
	if(IDopderrC == -1)
		IDopderrC = load_fits("OPDerrC.fits", "OPDerrC", 0);

	if(IDopderrC != -1)
		{
			nbOPDerr = data.image[IDopderrC].md[0].size[2];  // number of error arrays
			//printf("INCLUDING %ld OPD ERROR MODES\n", nbOPDerr);
			//fflush(stdout);
		}
	else
		{
			//printf("NO OPD ERROR MODES\n");
			//fflush(stdout);
			nbOPDerr = 0;
		}
    // focal plane plate scale in lambda/D per pixel
    focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;


    // CREATE SCORING MASK IF IT DOES NOT EXIST
    // which is the array of evaluation zones on the focal plane
    if((IDsm=image_ID("scoringmask"))==-1)
    {
		printf("CREATING SCORING MASK\n");
        printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);
		fflush(stdout);
        IDsm = create_2Dimage_ID("scoringmask", size, size);

        if(SCORINGMASKTYPE==0) // high density, wide
        {
            // draw an array of points, clip to desired subregions
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {   // create regular array and the raduis of each point
                    x = (1.0*ii-0.5*size)*focscale;
                    y = (1.0*jj-0.5*size)*focscale;
                    r = sqrt(x*x+y*y);

                    // clip the regular array to the desired annulus inside high-resolution
                    // region defined by scoringOWAhr
                    // use every other point
                    // and clip to an x > scoringIWAx part of the annulus if desired
                    if((r>scoringIWA)&&(r<scoringOWAhr)&&(x>scoringIWAx)&&((ii+jj)%2==0))
                        data.image[IDsm].array.F[jj*size+ii] = 1.0;
                    // pick every other row and column between scoringOWAhr and scoringOWA
                    if((r>scoringOWAhr)&&(r<scoringOWA)&&(x>scoringIWAx)&&(ii%2==0)&&(jj%2==0))
                        data.image[IDsm].array.F[jj*size+ii] = 1.0;
                    // draw a single radial line of points out to IWA = 70 (every other point)
                    if((x>scoringOWA)&&(fabs(y)<scoringIWA*0.25)&&(r<70.0)&&((ii+jj)%2==0)) // single line
                        data.image[IDsm].array.F[jj*size+ii] = 1.0;
                }
        }
        else // focused on central pixels, fewer pixels for faster convergence - used for FPMresp based optimization
        {
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    x = (1.0*ii-0.5*size)*focscale;
                    y = (1.0*jj-0.5*size)*focscale;
                    r = sqrt(x*x+y*y);
                    // clip from scoringIWA to scoringOWAhr only using every other column and row
                    if((r>scoringIWA)&&(r<scoringOWAhr)&&(x>scoringIWAx)&&(ii%2==0)&&(jj%2==0))
                        data.image[IDsm].array.F[jj*size+ii] = 1.0;
                }
        }
        if(PIAACMC_save==1)
        {
            // save a disgnostic image
            sprintf(fname, "!%s/scoringmask%d.fits", piaacmcconfdir, SCORINGMASKTYPE);
            save_fits("scoringmask", fname);
        }
        // a pixtable is a list of non-zero pixels with their coordinates
        linopt_imtools_mask_to_pixtable("scoringmask", "pixindex", "pixmult");
        // sums the image, giving the total number of pixels in the scoring mask
        SCORINGTOTAL = arith_image_total("scoringmask");

        //exit(0);

    }






    if(computePSF_FAST_FPMresp==1) // only possible if mode 11 has already been executed
    {
        // compute the PSF as the complex amplitude for the evaluation points on the focal plane
        // for a given FPM zone thickness based on the FPMresp array computed in mode 11
        value1 = PIAACMCsimul_achromFPMsol_eval(fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], optsyst[0].nblambda);
        
        // PSF result is stored in outtmp_array

        value = 0.0;
        peakcontrast = 0.0;

        ID = image_ID("imvect"); // use imvect if it exists
        if(ID==-1)
            ID = create_2Dimage_ID("imvect", vsize*optsyst[0].nblambda, 1);
        // write the result into imvect (= ID)
        for(ii=0; ii<vsize*optsyst[0].nblambda; ii++) // for each wavelength
        {
            data.image[ID].array.F[ii] = outtmp_array[ii];
            // square to give intensity
            tmpv = outtmp_array[ii]*outtmp_array[ii];
            // total intensity = sum(intensity_i) = sum(Re_i^2 + Im_i^2)
            value += tmpv;
        }
        // here value is the total flux in the output vector
        // store in global as total flux
        PIAACMCSIMUL_VAL0 = value;

        // morally: value -> average value per area normalized to flux
        value = value/size/size/optsyst[0].flux[0]; // flux[0] is proportional to the number of lambda channels, so this normalization makes value independant of number of spectral channels
        // here value is the total light (averaged across spectral channels) in the measurement points, normalized to the input flux
        // actual average contrast, averaging over # of pixels and physical area
        avContrast = value/(SCORINGTOTAL*focscale*focscale);


        //        printf("Peak constrast (rough estimate)= %g\n", peakcontrast/size/size/optsyst[0].flux[0]/focscale/focscale*3.0);
        //        printf("value1 = %g\n", value1);
        //		printf("Total light in scoring field = %g  -> Average contrast = %g   (%g)\n", value, value/(arith_image_total("scoringmask")*focscale*focscale), value1/CnormFactor/optsyst[0].nblambda);

    }
    else // we need to create the PSF from scratch
    {
        // if non-zero we are going to approximate the PSF for an extended source as
        // a collection of point sources
        // sourcesize determines the separation of the point sources
        if(sourcesize!=0) // sourcesize > 0 only if in linear optimization (step >= 100)
        {
            printf("COMPUTING RESOLVED SOURCE PSF / ADDING OPD MODES\n");	
            fflush(stdout);

            // dld is the radius of the circle containing the point sources in lamba/D
            dld = 1.0/pow(10.0, 0.1*sourcesize); // nominal pointing offset [l/D]
            // extmode controls how many point sources we use to model the extended source
            // extmode == 0 => three point sources
            // extmode == 1 => six point sources
            if (extmode==1)
            {
                // if I have six sources, put them in two rings of three each
                rad1 = dld/sqrt(2.5);
                rad2 = 2.0*dld/sqrt(2.5);
            }
            else
            {
                // if I have three sources keep them at the same radii
                rad1 = dld;
                rad2 = dld;
            }
            
            // we will collect propagation results in a set of vectors called imvectp<1,2,...>
            // which will ultimately be collected into a single long imvect

            // image index, counts the number of PSFs we make, one for each point source
			imindex = 0;
            // name of the image vector to hold the propagation result
			sprintf(imname, "imvectp%02ld", imindex);
            // xld and yld are the input positions of the input source in lamba/D
            // initialize first point source, which sets optsyst
            PIAACMCsimul_init(piaacmc, 0, xld+rad1, yld);	
			PIAACMCsimul_makePIAAshapes(piaacmc, 0);
			
            // propagate it (optsyst is a global), output in psfc0 (complex amlitude)
            // and psfi0 (intensity)
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0);
            // convert the image to the vector
            linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
            // save the intensity of the first point
            copy_image_ID("psfi0", "psfi0ext", 0);
			sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcconfdir, imindex);
			save_fits("psfi0", fname);
            
            // do the same for the second point
			imindex++;
			sprintf(imname, "imvectp%02ld", imindex);
            pha = 2.0*M_PI/3.0; // 1/3 around the circle
            PIAACMCsimul_init(piaacmc, 0, xld+rad1*cos(pha), yld+rad1*sin(pha)); 
            PIAACMCsimul_makePIAAshapes(piaacmc, 0);
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0);
            linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
            // add the intensity to build up PSF for extended source
            arith_image_add_inplace("psfi0ext", "psfi0");
			sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcconfdir, imindex);
			save_fits("psfi0", fname);
			
            // do the same for the third point
			imindex++;
			sprintf(imname, "imvectp%02ld", imindex);
            pha = 4.0*M_PI/3.0; // 2/3 around the circle
            PIAACMCsimul_init(piaacmc, 0, xld+rad1*cos(pha), yld+rad1*sin(pha));
            PIAACMCsimul_makePIAAshapes(piaacmc, 0);
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0);
            linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
            // add the intensity to build up PSF for extended source
            arith_image_add_inplace("psfi0ext", "psfi0");
       		sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcconfdir, imindex);
			save_fits("psfi0", fname);
			
            if (extmode==1)
            { // keep going for the other three points if desired, on the outer radius
				imindex++;
				sprintf(imname, "imvectp%02ld", imindex);
                pha = M_PI/3.0;
                PIAACMCsimul_init(piaacmc, 0, xld+rad2*cos(pha), yld+rad2*sin(pha));
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0);
                linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
                arith_image_add_inplace("psfi0ext","psfi0");
				sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcconfdir, imindex);
				save_fits("psfi0", fname);
			

				imindex++;
				sprintf(imname, "imvectp%02ld", imindex);
                pha = 2.0*M_PI/3.0 + M_PI/3.0;
                PIAACMCsimul_init(piaacmc, 0, xld+rad2*cos(pha), yld+rad2*sin(pha));
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0);
                linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
                arith_image_add_inplace("psfi0ext", "psfi0");
				sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcconfdir, imindex);
				save_fits("psfi0", fname);

				imindex++;
				sprintf(imname, "imvectp%02ld", imindex);
                pha = 4.0*M_PI/3.0 + M_PI/3.0;
                PIAACMCsimul_init(piaacmc, 0, xld+rad2*cos(pha), yld+rad2*sin(pha));
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0);
                linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
                arith_image_add_inplace("psfi0ext", "psfi0");
				sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcconfdir, imindex);
				save_fits("psfi0", fname);
			
                // but multiply by 0.5 'cause we have twice as many points
				//    arith_image_cstmult_inplace("psfi0ext", 0.5);
            }
			
			printf("Adding optional OPD error modes (%ld modes)\n", nbOPDerr);
			fflush(stdout);
			// add error modes if any
            // add new evaluation points for the error to imvect so we minimize the error
            // as well as the non-error
			for(OPDmode=0; OPDmode < nbOPDerr; OPDmode++)
			{
				imindex++;
				sprintf(imname, "imvectp%02ld", imindex);

				IDopderr = create_2Dimage_ID("opderr", size, size);
                // "opderr" is a standard name read by PIAACMCsimul_init
				for(ii=0;ii<size*size;ii++)
					data.image[IDopderr].array.F[ii] = data.image[IDopderrC].array.F[size*size*OPDmode + ii];
				PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0); // add error to the data
				PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0);
                linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
                arith_image_add_inplace("psfi0ext","psfi0");
	       	//	sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcconfdir, imindex);
			//	save_fits("psfi0", fname);
				delete_image_ID("opderr");
			}
			
            // now average over all the PSFs we've created to simulate this extended source
			NBimindex = imindex;
			arith_image_cstmult_inplace("psfi0ext", 1.0/NBimindex);


            if(outsave==1)
            {
                sprintf(fname, "!%s/psfi0_exsrc%3d_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, sourcesize, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);    
                
                save_fits("psfi0ext", fname);
            }

            // get the number of elements in a single-PSF vector
            ID = image_ID("imvectp00");
            nbelem = data.image[ID].md[0].nelement;
            // make big vector to collect the complex amplitudes of all the above PSFs
            ID = image_ID("imvect");
            if(ID!=-1)
                delete_image_ID("imvect");

            offset = nbelem/piaacmc[0].nblambda; // number of pixels per lambda x2 (re, im)
            printf("offset = %ld\n", offset);

            // number of pixels per lambda x 2 times the number of PSFs
            offset1 = NBimindex*offset;
            normcoeff = 1.0/sqrt(NBimindex);

            // make an imvect for each lambda
            ID = create_2Dimage_ID("imvect", offset1, piaacmc[0].nblambda);


            // fill in with each imvectp* created above
			for(imindex=0; imindex<NBimindex; imindex++)
			{
				sprintf(imname, "imvectp%02ld", imindex);
				ID1 = image_ID(imname);
				for(kl=0; kl<piaacmc[0].nblambda; kl++)
					for(ii=0; ii<offset; ii++)
						data.image[ID].array.F[kl*offset1 + imindex*offset + ii] = data.image[ID1].array.F[kl*offset+ii]*normcoeff;
				delete_image_ID(imname);
			}
/*
            ID1 = image_ID("imvectp0");
            for(kl=0; kl<piaacmc[0].nblambda; kl++)
                for(ii=0; ii<offset; ii++)
                    data.image[ID].array.F[kl*offset1+ii] = data.image[ID1].array.F[kl*offset+ii]*normcoeff;
            delete_image_ID("imvectp0");

            ID1 = image_ID("imvectp1");
            for(kl=0; kl<piaacmc[0].nblambda; kl++)
                for(ii=0; ii<offset; ii++)
                    data.image[ID].array.F[kl*offset1+offset+ii] = data.image[ID1].array.F[kl*offset+ii]*normcoeff;
            delete_image_ID("imvectp1");

            ID1 = image_ID("imvectp2");
            for(kl=0; kl<piaacmc[0].nblambda; kl++)
                for(ii=0; ii<offset; ii++)
                    data.image[ID].array.F[kl*offset1+2*offset+ii] = data.image[ID1].array.F[kl*offset+ii]*normcoeff;
            delete_image_ID("imvectp2");


            if(extmode==1)
            {
                ID1 = image_ID("imvectp3");
                for(kl=0; kl<piaacmc[0].nblambda; kl++)
                    for(ii=0; ii<offset; ii++)
                        data.image[ID].array.F[kl*offset1+3*offset+ii] = data.image[ID1].array.F[kl*offset+ii]*normcoeff;
                delete_image_ID("imvectp3");

                ID1 = image_ID("imvectp4");
                for(kl=0; kl<piaacmc[0].nblambda; kl++)
                    for(ii=0; ii<offset; ii++)
                        data.image[ID].array.F[kl*offset1+4*offset+ii] = data.image[ID1].array.F[kl*offset+ii]*normcoeff;
                delete_image_ID("imvectp4");

                ID1 = image_ID("imvectp5");
                for(kl=0; kl<piaacmc[0].nblambda; kl++)
                    for(ii=0; ii<offset; ii++)
                        data.image[ID].array.F[kl*offset1+5*offset+ii] = data.image[ID1].array.F[kl*offset+ii]*normcoeff;
                delete_image_ID("imvectp5");
            }        
*/



            //linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", "imvect");
            
            // measure the contrast for all aimplitudes in imvect
            value = 0.0;
            peakcontrast = 0.0;
            ID = image_ID("imvect");
            for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
            {
                tmpv = data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
                value += tmpv;
                if(tmpv>peakcontrast)
                    peakcontrast = tmpv;
            }

            for(elem=0; elem<optsyst[0].NBelem; elem++)
                printf("    FLUX %3ld   %12.4lf %8.6lf\n", elem, optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0]);

            //            sprintf(fname,"%s/flux.txt", piaacmcconfdir);
            
            if(outsave==1)
            {
                sprintf(fname, "!%s/flux_exsrc%3d_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcconfdir, sourcesize, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda); 
                
                fpflux = fopen(fname, "w");
                for(elem=0; elem<optsyst[0].NBelem; elem++)
                    fprintf(fpflux, "%18.16lf %18.16lf  %d\n", optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0], optsyst[0].nblambda);
                fprintf(fpflux, "W0  %d\n", optsyst[0].nblambda);
                fclose(fpflux);
            }



            // compute average contrast
            value = value/size/size/optsyst[0].flux[0];
            avContrast = value/(SCORINGTOTAL*focscale*focscale);

            //           CnormFactor = size*size*optsyst[0].flux[0]*arith_image_total("scoringmask")*focscale*focscale; // /optsyst[0].nblambda;
            CnormFactor = focscale*focscale*size*size*optsyst[0].flux[0]/optsyst[0].nblambda;

            if(WRITE_OK==1)
            {
                sprintf(fname, "%s/CnormFactor.txt", piaacmcconfdir);
                fp = fopen(fname, "w");
                fprintf(fp, "%g\n", CnormFactor);
                fprintf(fp, "0      %g %ld %g %d\n", focscale, size, optsyst[0].flux[0], optsyst[0].nblambda);
                fclose(fp);
            }

            printf("COMPUTING RESOLVED SOURCE PSF\n");
            printf("SCORINGTOTAL = %f  %f\n", SCORINGTOTAL, arith_image_total("scoringmask"));
            printf("Peak constrast (rough estimate)= %g\n", peakcontrast/size/size/optsyst[0].flux[0]/focscale/focscale/normcoeff/normcoeff);
            avContrast = value/(arith_image_total("scoringmask")*focscale*focscale);
            printf("Total light in scoring field = %g, peak PSF = %g   -> Average contrast = %g\n", value, piaacmc[0].peakPSF, avContrast);
            

            if(outsave==1)
            {
                sprintf(fname, "%s/contrast_exsrc%3d_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcconfdir, sourcesize, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);    
                
                fp = fopen(fname, "w");
                fprintf(fp, "%g", avContrast);
                fclose(fp);
            }
        }
        else // called for step 0 through 15.  Does not use OPDerr
        { // compute the PSF for a single point source at offset xld, yld
            printf("COMPUTING UNRESOLVED SOURCE PSF [%f x %f]\n", xld, yld);


            // ========== initializes optical system to piaacmc design ===========
            // xld and yld are the input positions of the input source in lamba/D
            // initialize first point source, which sets optsyst
            PIAACMCsimul_init(piaacmc, 0, xld, yld);
            PIAACMCsimul_makePIAAshapes(piaacmc, 0);



            // ============ perform propagations ================
            // propagate it (optsyst is a global), output in psfc0 (complex amlitude)
            // and psfi0 (intensity)
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0);

            if(outsave==1)
            {
              sprintf(fname, "!%s/psfi0_ptsrc_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);    
                
               save_fits("psfi0", fname);
            }


 //           list_image_ID();
            // linearize the result into imvect
            linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", "imvect");
           
           
            save_fits("imvect", "!test_imvect.fits");
            // extract amplitude and phase for diagnostics
            mk_amph_from_complex("psfc0", "psfc0a", "psfc0p", 0);
            save_fits("psfc0a", "!test_psfc0a.fits");
            list_image_ID();
            delete_image_ID("psfc0a");
            delete_image_ID("psfc0p");
            printf("saved -> test_psfc0a.fits\n");
            fflush(stdout);
            // compute average contrast
            value = 0.0;
            peakcontrast = 0.0;
            ID = image_ID("imvect");
            for(ii=0; ii<data.image[ID].md[0].nelement; ii+=2)
            {
                // intensity as Re^2 + Im^2
                tmpv = data.image[ID].array.F[ii]*data.image[ID].array.F[ii] + data.image[ID].array.F[ii+1]*data.image[ID].array.F[ii+1];
                value += tmpv;
                if(tmpv>peakcontrast)
                    peakcontrast = tmpv;
            }
            // report the contrast
            for(elem=0; elem<optsyst[0].NBelem; elem++)
                printf("    FLUX %3ld   %12.4lf %8.6lf\n", elem, optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0]);
//            value = value/size/size/optsyst[0].flux[0];

            if(outsave==1)
            {
                sprintf(fname, "%s/flux_ptsrc_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda); 

                fpflux = fopen(fname, "w");
                for(elem=0; elem<optsyst[0].NBelem; elem++)
                    fprintf(fpflux, "%18.16lf %18.16lf  %d\n", optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0], optsyst[0].nblambda);
                fprintf(fpflux, "W1\n");
                fclose(fpflux);
        
                sprintf(command, "cp %s %s/flux.txt", fname, piaacmcconfdir);
                ret = system(command);

           }

  
                
            //         CnormFactor = size*size*optsyst[0].flux[0]*arith_image_total("scoringmask")*focscale*focscale; // /optsyst[0].nblambda;
            CnormFactor = focscale*focscale*size*size*optsyst[0].flux[0]/optsyst[0].nblambda;
            sprintf(fname, "%s/CnormFactor.txt", piaacmcconfdir);
            fp = fopen(fname, "w");
            fprintf(fp, "%g\n", CnormFactor);
            fprintf(fp, "1      %g %ld %g %d\n", focscale, size, optsyst[0].flux[0], optsyst[0].nblambda);
            fclose(fp);
            
            // here we're essentially done!


            printf("COMPUTING UNRESOLVED SOURCE PSF -*- [%f x %f]\n", xld, yld);
            printf("SCORINGTOTAL = %f  %f\n", SCORINGTOTAL, arith_image_total("scoringmask"));
   
            if((IDv=variable_ID("PIAACMC_NOFPM"))!=-1)
                {
                    ID = image_ID("psfc0");
                    piaacmc[0].peakPSF = 0.0;
                    for(ii=0;ii<size*size;ii++)
                        {
                            val = data.image[ID].array.CF[ii].re*data.image[ID].array.CF[ii].re + data.image[ID].array.CF[ii].im*data.image[ID].array.CF[ii].im;
                            if(val>piaacmc[0].peakPSF)
                                piaacmc[0].peakPSF = val;
                        }

                    fp = fopen("conf/conf_peakPSF.txt", "w");
                    fprintf(fp, "%g\n", piaacmc[0].peakPSF);
                    fclose(fp);
                }
            if(piaacmc[0].peakPSF<1.0)
                {
                    printf("Peak constrast (rough estimate)= %g -> %g\n", peakcontrast, peakcontrast/(optsyst[0].flux[0]*optsyst[0].flux[0]));
                    printf("optsyst[0].flux[0]  = %g\n", optsyst[0].flux[0]);
                    printf("SCORINGMASKTYPE = %d\n", SCORINGMASKTYPE);
                    avContrast = value/(optsyst[0].flux[0]*optsyst[0].flux[0])/SCORINGTOTAL;
                    printf("[0] Total light in scoring field = %g, peak PSF = %g, SCOTINGTOTAL = %g   -> Average contrast = %g\n", value, piaacmc[0].peakPSF,  SCORINGTOTAL, avContrast);
                }
            else
                {
                    printf("Peak constrast = %g -> %g\n", peakcontrast, peakcontrast/piaacmc[0].peakPSF);
                    printf("SCORINGMASKTYPE = %d\n", SCORINGMASKTYPE);
                    avContrast = value/piaacmc[0].peakPSF/SCORINGTOTAL;
                    printf("[1] Total light in scoring field = %g, peak PSF = %g, SCOTINGTOTAL = %g  -> Average contrast = %g\n", value, piaacmc[0].peakPSF, SCORINGTOTAL, avContrast);
                    
                    if((fp=fopen("PSFcontrastval.txt","w"))!=NULL)
                    {
                        fprintf(fp, "%g %g\n", peakcontrast/piaacmc[0].peakPSF, value/piaacmc[0].peakPSF/SCORINGTOTAL);
                        fclose(fp);
                    }
                }

            if(outsave==1)
            {    
                sprintf(fname, "%s/contrast_ptsrc_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);    
                printf("saving contrast value [%g] -> %s\n", avContrast, fname);
                
                fp = fopen(fname, "w");
                fprintf(fp, "%g", avContrast);
                fclose(fp);
            }
        }
    }

    return(avContrast);
}











int PIAAsimul_savepiaacmcconf(char *dname)
{
    char command[1000];
    int r;
    FILE *fp;
    char fname[500];
    long i;

    
    sprintf(command, "mkdir -p %s", dname);
    r = system(command);

    sprintf(fname,"%s/piaacmcparams.conf", dname);
    fp = fopen(fname, "w");


    fprintf(fp, "%10ld   NBradpts\n", piaacmc[0].NBradpts);

    for(i=0; i<10; i++)
    {
        if(i<piaacmc[0].NBLyotStop)
        {
            sprintf(fname, "!%s/LyotStop%ld.fits", dname, i);
            if(piaacmc[0].IDLyotStop[i]!=-1)
                save_fits(data.image[piaacmc[0].IDLyotStop[i]].name, fname);
            fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
        }
        else
            fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
    }

    sprintf(fname, "!%s/piaa0Cmodes.fits", dname);
    if(piaacmc[0].piaa0CmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa0CmodesID].name, fname);

    sprintf(fname, "!%s/piaa0Fmodes.fits", dname);
    if(piaacmc[0].piaa0FmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa0FmodesID].name, fname);

    sprintf(fname, "!%s/piaa1Cmodes.fits", dname);
    if(piaacmc[0].piaa1CmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa1CmodesID].name, fname);

    sprintf(fname, "!%s/piaa1Fmodes.fits", dname);
    if(piaacmc[0].piaa1FmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa1FmodesID].name, fname);


    sprintf(fname, "!%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

    
    if(piaacmc[0].zonezID!=-1)
        save_fits(data.image[piaacmc[0].zonezID].name, fname);

    fprintf(fp, "%10.6f    fpmaskamptransm\n", piaacmc[0].fpmaskamptransm);

    sprintf(fname, "!%s/fpm_zonea_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

    
    if(piaacmc[0].zoneaID!=-1)
        save_fits(data.image[piaacmc[0].zoneaID].name, fname);

    fclose(fp);

    

    return(0);
}



int PIAAsimul_loadpiaacmcconf(char *dname)
{
    char command[1000];
    int r;
    FILE *fp;
    char fname[500];
    char imname[500];
    long i;

    int tmpi;
    long tmpl;
    float tmpf;
    double tmplf;



    sprintf(fname,"%s/piaacmcparams.conf", dname);
    printf("%s\n", fname);


    fp = fopen(fname, "r");
    if(fp==NULL)
    {
        printf("Configuration file \"%s\" does not exist (yet), using previously set configuration\n", fname);
        fflush(stdout);
        r = 0;
    }
    else
    {

        r = fscanf(fp, "%ld   NBradpts\n", &tmpl);
        piaacmc[0].NBradpts = tmpl;

         for(i=0; i<10; i++)
        {
            if(i<piaacmc[0].NBLyotStop)
            {
                sprintf(fname, "%s/LyotStop%ld.fits", dname, i);
                sprintf(imname, "lyotstop%ld", i);
                printf("Loading \"%s\" as \"%s\"\n", fname, imname);
                piaacmc[0].IDLyotStop[i] = load_fits(fname, imname, 1);
                r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                piaacmc[0].LyotStop_zpos[i] = tmplf;
            }
            else
            {
                r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                piaacmc[0].LyotStop_zpos[i] = tmplf;
            }
        printf("LYOT STOP %ld POS : %lf\n", i, tmplf);
            
        }


        sprintf(fname, "%s/piaa0Cmodes.fits", dname);
        piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff", 1);
        if(piaacmc[0].piaa0CmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa0Cmodes.fits", dname);
            piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff", 1);
        }


        sprintf(fname, "%s/piaa0Fmodes.fits", dname);
        piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff", 1);
        if(piaacmc[0].piaa0FmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa0Fmodes.fits", dname);
            piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff", 1);
        }

        sprintf(fname, "%s/piaa1Cmodes.fits", dname);
        piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff", 1);
        if(piaacmc[0].piaa1CmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa1Cmodes.fits", dname);
            piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff", 1);
        }

        sprintf(fname, "%s/piaa1Fmodes.fits", dname);
        piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff", 1);
        if(piaacmc[0].piaa1FmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa1Fmodes.fits", dname);
            piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff", 1);
        }


        r = fscanf(fp, "%f   fpmaskamptransm\n",    &tmpf);
        piaacmc[0].fpmaskamptransm = tmpf;

        r = 1;

        fclose(fp);
    }



    return(r);
}


/// Make Lyot stop geometry
/// param[in] IDincoh_name   Incoherent Lyot pupil intensity response to off-axis sources
/// parampin] IDmc_name      Intensity Lyot pupil image for on-axis source
///
/// explores two thresholding methods applied together :
/// (1) keeps pixels for which offaxisLight / onaxisLight > rsl
/// (2) keeps pixels for which onaxisLight < v0
/// selects the mask that achieves the strongest on-axis rejection while satifying the throughput constraint

long PIAACMCsimul_mkLyotMask(char *IDincoh_name, char *IDmc_name, char *IDzone_name, double throughput, char *IDout_name)
{
    long ID, ID1;
    long IDmc, IDincoh, IDzone;
    double val, val1, v, v0, bestval, v_best, rsl_best;
    double rsl, rsl0;
    long iter, NBiter;
    long ii;
    long xsize = 0;
    long ysize = 0;
    long IDout;
    float sigma = 4.0;
    int filter_size = 10;



    NBiter = 100;

    sigma = 0.01*piaacmc[0].beamrad/piaacmc[0].pixscale;
    filter_size = (long) (sigma*2.0);

    printf("IDincoh_name : %s   %ld\n", IDincoh_name, image_ID(IDincoh_name));
    printf("IDmc_name    : %s   %ld\n", IDmc_name, image_ID(IDmc_name));
    printf("IDzone_name  : %s   %ld\n", IDzone_name, image_ID(IDzone_name));

    //IDincoh = gauss_filter(IDincoh_name, "incohg", sigma, filter_size);
    IDincoh = image_ID(IDincoh_name);

    IDmc = image_ID(IDmc_name);
    //	IDmc = gauss_filter(IDmc_name, "mcg", sigma, filter_size);

    IDzone = image_ID(IDzone_name);
    xsize = data.image[IDmc].md[0].size[0];
    ysize = data.image[IDmc].md[0].size[1];

    IDout = create_2Dimage_ID(IDout_name, xsize, ysize);

    // normalize both images to 1.0
    val = 0.0;
    for(ii=0; ii<xsize*ysize; ii++)
        val += data.image[IDmc].array.F[ii];
    for(ii=0; ii<xsize*ysize; ii++)
        data.image[IDmc].array.F[ii] /= val;

    val = 0.0;
    for(ii=0; ii<xsize*ysize; ii++)
        val += data.image[IDincoh].array.F[ii];
    for(ii=0; ii<xsize*ysize; ii++)
        data.image[IDincoh].array.F[ii] /= val;



    // estimate iteratively rsl0, the threshold in offaxis/onaxis starlight that will achieve the goal throughput
    rsl = 1.0;
    for(iter=0; iter<NBiter; iter++)
    {
        val = 0.0;
        val1 = 0.0;

        for(ii=0; ii<xsize*ysize; ii++)
        {
            if((data.image[IDzone].array.F[ii]>-1)&&(data.image[IDincoh].array.F[ii]/data.image[IDmc].array.F[ii]>rsl))
            {
                val += data.image[IDincoh].array.F[ii];
                val1 += data.image[IDmc].array.F[ii];
                data.image[IDout].array.F[ii] = 1.0;
            }
            else
                data.image[IDout].array.F[ii] = 0.0;
        }
        printf("rsl = %f  ->  %f %f   (%f)\n", rsl, val, val1, throughput);
        if(val>throughput) // too much light came through
            rsl *= 1.1;
        else
            rsl *= 0.9;
    }
    rsl0 = rsl;

    // v0 = img_percentile("mcg", 0.99);
    v0 = img_percentile(IDmc_name, 0.99);
    printf("v0 = %lf\n", v0);


    

    bestval = 1.0; // to be minized: total starlight transmitted
    for(rsl=0.0*rsl0; rsl< 2.0*rsl0; rsl+=0.02*rsl0)
        for(v=0.00000001*v0; v<50.0*v0; v*=1.2)
        {
            val = 0.0;
            val1 = 0.0;

            for(ii=0; ii<xsize*ysize; ii++)
            {
                if((data.image[IDzone].array.F[ii]>-1)&&(data.image[IDincoh].array.F[ii]/data.image[IDmc].array.F[ii]>rsl)&&(data.image[IDmc].array.F[ii]<v))
                {
                    val += data.image[IDincoh].array.F[ii];
                    val1 += data.image[IDmc].array.F[ii];
                }
            }

            if(val>throughput)
            {
                if(val1<bestval)
                {
                    bestval = val1;
                    rsl_best = rsl;
                    v_best = v;
                    printf("BEST SOLUTION: %.12lf / %.12lf    %.12lf / %.12lf  -> %.12lf  %.12lf\n", rsl_best, rsl0, v_best, v0, val, bestval);
                }

            }
        }

    for(ii=0; ii<xsize*ysize; ii++)
    {
        if((data.image[IDzone].array.F[ii]>-1)&&(data.image[IDincoh].array.F[ii]/data.image[IDmc].array.F[ii]>rsl_best)&&(data.image[IDmc].array.F[ii]<v_best))
            data.image[IDout].array.F[ii] = 1.0;
        else
            data.image[IDout].array.F[ii] = 0.0;
    }

    if(0)
    {
        ID1 = create_2Dimage_ID("postLMim", xsize, ysize);
        for(ii=0; ii<xsize*ysize; ii++)
            data.image[ID1].array.F[ii] = data.image[IDmc].array.F[ii]*data.image[IDout].array.F[ii];
        save_fits("postLMim", "!postLMim.fits");
        delete_image_ID("postLMim");
    }


    return(IDout);
}




///
/// propagate complex amplitude image into intensity map cube
///
long PIAACMCsimul_CA2propCubeInt(char *IDamp_name, char *IDpha_name, float zmin, float zmax, long NBz, char *IDout_name)
{
    long IDout;
    long l;
    long IDa, IDp;
    long xsize = 0;
    long ysize = 0;
    long ii;
    long nblambda, k;
    float *zarray;
    float zprop;
    long IDre, IDim;
    double amp, pha;
    long IDreg, IDimg;
    double re, im;
    long IDintg, IDintgg;


    printf("PIAACMCsimul_CA2propCubeInt   : %s %s %f %f %ld %s\n", IDamp_name, IDpha_name, zmin, zmax, NBz, IDout_name);
    fflush(stdout);
    
    list_image_ID();
    
    IDa = image_ID(IDamp_name);
    xsize = data.image[IDa].md[0].size[0];
    ysize = data.image[IDa].md[0].size[1];

    IDre = create_2Dimage_ID("retmpim", xsize, ysize);
    IDim = create_2Dimage_ID("imtmpim", xsize, ysize);

    if(data.image[IDa].md[0].naxis==3)
    nblambda = data.image[IDa].md[0].size[2];
    else
    nblambda = 1;

    IDout = create_3Dimage_ID(IDout_name, xsize, ysize, NBz);
    IDintg = create_2Dimage_ID("tmpintg", xsize, ysize);

    
   // initialize zarray
    zarray = (float*) malloc(sizeof(float)*NBz);
    for(l=0; l<NBz; l++)
        zarray[l] = zmin + (zmax-zmin)*l/(NBz-1);


    
    for(l=0; l<NBz; l++)
    {
        printf("l = %ld/%ld\n", l, NBz);
        fflush(stdout);
        
        zprop = zarray[l];
        OptSystProp_propagateCube(optsyst, 0, IDamp_name, IDpha_name, "_tmppropamp", "_tmpproppha", zprop, 0);

        IDa = image_ID("_tmppropamp");
        IDp = image_ID("_tmpproppha");


        // write intensity
        for(k=0; k<nblambda; k++)
            for(ii=0; ii<xsize*ysize; ii++)
                data.image[IDout].array.F[l*xsize*ysize+ii] += data.image[IDa].array.F[ii]*data.image[IDa].array.F[k*xsize*ysize+ii]; 
        


        delete_image_ID("_tmppropamp");
        delete_image_ID("_tmpproppha");
    }
    
    free(zarray);
    delete_image_ID("retmpim");
    delete_image_ID("imtmpim");
    
    return IDout;
}



/// make Lyot stops using off-axis light minimums
/// finds minumum flux level in intensity data cube

long PIAACMCsimul_optimizeLyotStop_OAmin(char *IDincohc_name)
{
    long IDincohc;
    
    long IDindex;
    long IDminflux;
    long ii, jj, kk;
    long xsize = 0;
    long ysize = 0;
    long xysize = 0;
    double minv;
    long minindex;
    double tmpv;
    long NBz;
    
    
    IDincohc = image_ID(IDincohc_name);

    xsize = data.image[IDincohc].md[0].size[0];
    ysize = data.image[IDincohc].md[0].size[1];
    xysize = xsize*ysize;
    NBz = data.image[IDincohc].md[0].size[2];
    
    IDindex = create_2Dimage_ID("oals_index", xsize, ysize);
    IDminflux = create_2Dimage_ID("oals_val", xsize, ysize);
    
    for(ii=0;ii<xsize;ii++)
        for(jj=0;jj<ysize;jj++)
            {
                minv = data.image[IDincohc].array.F[jj*xsize+ii];
                minindex = 0;

                for(kk=1;kk<NBz;kk++)
                {
                    tmpv = data.image[IDincohc].array.F[kk*xysize+jj*xsize+ii];
                    if(tmpv<minv)
                        {
                            minv = tmpv;
                            minindex = kk;
                        }
                }
                data.image[IDminflux].array.F[jj*xsize+ii] = minv;
                data.image[IDindex].array.F[jj*xsize+ii] = minindex;
            }    
    
    save_fits("oals_index", "!test_oals_index.fits");
    save_fits("oals_val", "!test_oals_val.fits");
    
    return 0;
}



///
/// Lyot stops positions from zmin to zmax relative to current, working back (light goes from 0 to zmax)
/// @param[in] FluxTOT  total flux in current plane
/// @param[in] FluxLim   max flux allowed from star
///
double PIAACMCsimul_optimizeLyotStop(char *IDamp_name, char *IDpha_name, char *IDincohc_name, float zmin, float zmax, double throughput, long NBz, long NBmasks)
{
    // initial guess places Lyot stops regularly from zmin to zmax
    // light propagates from zmin to zmax
    // we start with a single mask in zmax, and work back
    //
    double ratio = 1.0; // output metric... not used yet, currently a placeholder

    long ID, IDa, IDp;
    long nblambda; // number of wavelengths, read from input cube
    float *zarray;
    long l;
    double zprop;

    char nameamp[500];
    char namepha[500];
    char nameint[500];
    char fname[500];
    char fname1[500];
    long xsize = 0;
    long ysize = 0;
    long ii, jj, k, m;

    float *rinarray;
    float *routarray;
    float dr = 0.02;
    double tot;


    double *totarray;
    double *tot2array;
    long IDzone;
    double x, y, r;

    double zbest, valbest, val;
    long lbest;

    long IDincohc, IDint, IDmc, IDmc1;
    float rsl;
    long iter;
    long NBiter = 100;
    double val1;
    long IDm;
    char name[500];

    long IDre, IDim, IDreg, IDimg;
    double amp, pha, re, im;
    float sigma, sigma1;
    int filter_size;

    long IDlscumul;
    long IDintg, IDintgg;

    float *rarray;


    FILE *fp;
    double alpha = 1.01; // norm alpha used to identify best plane


    sigma = 0.01*piaacmc[0].beamrad/piaacmc[0].pixscale;
    filter_size = (long) (sigma*2.0);

    zarray = (float*) malloc(sizeof(float)*NBz);

    rinarray = (float*) malloc(sizeof(float)*NBmasks);
    routarray = (float*) malloc(sizeof(float)*NBmasks);

    totarray = (double*) malloc(sizeof(double)*NBmasks*NBz);
    tot2array = (double*) malloc(sizeof(double)*NBmasks*NBz);


    routarray[0] = 1.0;
    rinarray[0] = 1.0 - 1.0/NBmasks;
    for(m=1; m<NBmasks; m++)
    {
        routarray[m] = rinarray[m-1];
        rinarray[m] = routarray[m] - (1.0-piaacmc[0].centObs1)/NBmasks;
    }
    rinarray[NBmasks-1] = 0.0;

    for(m=0; m<NBmasks; m++)
        printf("annulus %ld : %f - %f\n", m, routarray[m], rinarray[m]);


    IDa = image_ID(IDamp_name);
    IDp = image_ID(IDpha_name);
    xsize = data.image[IDa].md[0].size[0];
    ysize = data.image[IDa].md[0].size[1];

    
    rarray = (float*) malloc(sizeof(float)*xsize*ysize);
    IDincohc = image_ID(IDincohc_name);



    if(data.image[IDa].md[0].naxis==3)
        nblambda = data.image[IDa].md[0].size[2];
    else
        nblambda = 1;

    IDzone = create_2Dimage_ID("LMzonemap", xsize, ysize);
    for(ii=0; ii<xsize; ii++)
        for(jj=0; jj<ysize; jj++)
        {
            data.image[IDzone].array.F[jj*xsize+ii] = -2;
            x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
            y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
            r = sqrt(x*x+y*y);
            rarray[jj*xsize+ii] = r;
            for(m=0; m<NBmasks; m++)
                if((r>rinarray[m]-0.0001)&&(r<routarray[m]+0.0001))
                    data.image[IDzone].array.F[jj*xsize+ii] = m;
        }
    if(PIAACMC_save==1)
    {
        sprintf(fname, "!%s/LMzonemap.fits", piaacmcconfdir);
        save_fits("LMzonemap", fname);
    }
    // initialize zarray
    for(l=0; l<NBz; l++)
        zarray[l] = zmin + (zmax-zmin)*l/(NBz-1);

    //  save_fits(nameint, fname);


    IDre = create_2Dimage_ID("retmpim", xsize, ysize);
    IDim = create_2Dimage_ID("imtmpim", xsize, ysize);

    ID = create_3Dimage_ID("LMintC", xsize, ysize, NBz);
    IDintg = create_2Dimage_ID("tmpintg", xsize, ysize);

    for(l=0; l<NBz; l++)
    {
        sprintf(nameamp, "LMPamp%02ld", l);
        sprintf(namepha, "LMPpha%02ld", l);
        zprop = zarray[l];
        OptSystProp_propagateCube(optsyst, 0, IDamp_name, IDpha_name, nameamp, namepha, zprop, 0);

     
        IDa = image_ID(nameamp);
        IDp = image_ID(namepha);
     
        for(k=0; k<nblambda; k++)
        {
        for(ii=0; ii<xsize*ysize; ii++)
            {
                amp = data.image[IDa].array.F[k*xsize*ysize+ii];
                pha = data.image[IDp].array.F[k*xsize*ysize+ii];
                data.image[IDre].array.F[ii] = amp*cos(pha);
                data.image[IDim].array.F[ii] = amp*sin(pha);
            }
            IDreg = gauss_filter("retmpim", "retmpimg", sigma, filter_size);
            IDimg = gauss_filter("imtmpim", "imtmpimg", sigma, filter_size);

            for(ii=0; ii<xsize*ysize; ii++)
            {
                re = data.image[IDreg].array.F[ii];
                im = data.image[IDimg].array.F[ii];
                data.image[IDintg].array.F[ii] = re*re + im*im;
            }
            IDintgg = gauss_filter("tmpintg", "tmpintgg", 2.0*sigma, filter_size);


            for(ii=0; ii<xsize*ysize; ii++)
                data.image[ID].array.F[l*xsize*ysize+ii] += data.image[IDintgg].array.F[ii];
                //data.image[IDa].array.F[ii]*data.image[IDa].array.F[ii]; //data.image[IDintgg].array.F[ii];

            delete_image_ID("retmpimg");
            delete_image_ID("imtmpimg");
            delete_image_ID("tmpintgg");
        }

    /*    for(ii=0; ii<xsize*ysize; ii++)
        {
            m = (long) (data.image[IDzone].array.F[ii]+0.1);

    
            if((m>-1)&&(m<NBmasks)&&(rarray[ii]<1.0)&&(rarray[ii]>0.9*piaacmc[0].centObs1))
            {
                totarray[l*NBmasks+m] += data.image[ID].array.F[l*xsize*ysize+ii];
                tot2array[l*NBmasks+m] += pow(data.image[ID].array.F[l*xsize*ysize+ii], alpha);
            }
        }*/
        delete_image_ID(nameamp);
        delete_image_ID(namepha);
    }
    delete_image_ID("retmpim");
   delete_image_ID("imtmpim");


    sprintf(fname,  "!%s/LMintC.fits", piaacmcconfdir);
    save_fits("LMintC", fname);




    for(l=0; l<NBz; l++)
        for(m=0;m<NBmasks;m++)
            {
                totarray[l*NBmasks+m] = 0.0;
                tot2array[l*NBmasks+m] = 0.0;
            }

    for(l=0; l<NBz; l++)
    {
        for(ii=0; ii<xsize*ysize; ii++)
        {
            m = (long) (data.image[IDzone].array.F[ii]+0.1);
    
            if((m>-1)&&(m<NBmasks)&&(rarray[ii]<1.0)&&(rarray[ii]>0.9*piaacmc[0].centObs1))
            {
                totarray[l*NBmasks+m] += data.image[IDincohc].array.F[l*xsize*ysize+ii];
                tot2array[l*NBmasks+m] += pow(data.image[IDincohc].array.F[l*xsize*ysize+ii], alpha);
            }
        }
    }

    IDmc = create_2Dimage_ID("Lcomb", xsize, ysize);
    IDmc1 = create_2Dimage_ID("LcombOA", xsize, ysize);
    for(ii=0; ii<xsize*ysize; ii++)
        data.image[IDmc1].array.F[ii] = 0.0;

    sprintf(fname1, "%s/LyotMasks_zpos.txt", piaacmcconfdir);
    fp = fopen(fname1, "w");
    IDint = image_ID("LMintC");
    for(m=0; m<NBmasks; m++)
    {
        valbest = 0.0;
        lbest = 0;
        zbest = 0.0;
        for(l=0; l<NBz; l++)
        {
            val =  tot2array[l*NBmasks+m]/pow(totarray[l*NBmasks+m], alpha);
            printf("MASK %ld   z(%ld)= %f  ->  %g   ( %g %g) \n", m, l, zarray[l], val, tot2array[l*NBmasks+m], totarray[l*NBmasks+m]);
            if(val>valbest)
            {
                valbest = val;
                zbest = zarray[l];
                lbest = l;
            }
        }
        printf(" ==========  MASK %ld   BEST CONJUGATION : %ld %f (%g)\n", m, lbest, zbest, valbest);
        piaacmc[0].LyotStop_zpos[m] = zbest; // relative to starting plane
        fprintf(fp, "%02ld %f\n", lbest, zbest);

        for(ii=0; ii<xsize*ysize; ii++)
            if(m==data.image[IDzone].array.F[ii])
                {
                    data.image[IDmc].array.F[ii] = data.image[IDint].array.F[lbest*xsize*ysize+ii];
                    data.image[IDmc1].array.F[ii] = data.image[IDincohc].array.F[lbest*xsize*ysize+ii];
                }
    }
    fclose(fp);

    if(PIAACMC_save==1)
    {
        sprintf(fname, "!%s/Lcomb.fits", piaacmcconfdir);
        save_fits("Lcomb", fname);
        sprintf(fname, "!%s/LcombOA.fits", piaacmcconfdir);
        save_fits("LcombOA", fname);
    }
    
    ID = PIAACMCsimul_mkLyotMask("LcombOA", "Lcomb", "LMzonemap", throughput, "LMask");
    
    if(PIAACMC_save==1)
    {
        sprintf(fname, "!%s/LMask.fits", piaacmcconfdir);
        save_fits("LMask", fname);
    }
    delete_image_ID("Lcomb");

    IDlscumul = create_2Dimage_ID("LMcumul", xsize, ysize);
    for(ii=0; ii<xsize*ysize; ii++)
        data.image[IDlscumul].array.F[ii] = 1.0;

    for(m=0; m<NBmasks; m++)
    {
        sprintf(name, "optLM%02ld", m);
        IDm = create_2Dimage_ID(name, xsize, ysize);
        for(ii=0; ii<xsize; ii++)
            for(jj=0; jj<ysize; jj++)
            {
                x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
                y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
                r = sqrt(x*x+y*y);

                if((r>rinarray[m]-dr)&&(r<routarray[m]+dr))
                    data.image[IDm].array.F[jj*xsize+ii] = data.image[ID].array.F[jj*xsize+ii];
                else
                    data.image[IDm].array.F[jj*xsize+ii] = 1.0;
            }
        if(m==0)
            for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
                {
                    x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
                    y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
                    r = sqrt(x*x+y*y);
                    if(r>1.0)
                        data.image[IDm].array.F[jj*xsize+ii] = 0.0;
                }


        for(ii=0; ii<xsize*ysize; ii++)
        {
            data.image[IDm].array.F[ii] *= data.image[IDlscumul].array.F[ii];
            data.image[IDlscumul].array.F[ii] = data.image[IDm].array.F[ii];
        }


        sprintf(fname, "!%s/optLM%02ld.fits", piaacmcconfdir, m);
        save_fits(name, fname);
    }


    delete_image_ID("LMcumul");

    free(totarray);
    free(tot2array);
    free(rinarray);
    free(routarray);
    free(zarray);
    free(rarray);

    delete_image_ID("LMzonemap");

    return(ratio);
}




///
/// solves for focal plane mask solution using pre-computed zone responses
///
/// @param[in] fpmresp_array   Mask zones responses, double array
/// @param[in] zonez_array     zone thicknesses, double array
/// @param[in] dphadz_array    for each lambda, pha = thickness x dphadt_array[lambdaindex]
/// @param[out] outtmp_array    output temp array
///
/// written to be fast, no checking of array sizes
/// all arrays pre-allocated outside this function
///
double PIAACMCsimul_achromFPMsol_eval(double *fpmresp_array, double *zonez_array, double *dphadz_array, double *outtmp_array, long vsize, long nbz, long nbl)
{

    // axis 0: eval pts (ii)   size = data.image[IDfpmresp].md[0].size[0] -> vsize
    // axis 1: zones (mz)      size = data.image[piaacmc[0].zonezID].md[0].size[0]+1 = nbz+1
    // axis 3: lambda (k)      size = piaacmc[0].nblambda -> nbl
    //
    // indexing :  k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*vsize + mz*vsize + ii




    for(evalk=0; evalk<nbl; evalk++)
    {
        evalki = evalk*(nbz+1)*vsize;

        if(optsyst[0].FOCMASKarray[0].mode == 1) // include outer zone
        {
            // outer zone
            for(evalii=0; evalii<vsize/2; evalii++)
            {
                outtmp_array[evalk*vsize+2*evalii] = fpmresp_array[evalk*(nbz+1)*vsize+2*evalii];
                outtmp_array[evalk*vsize+2*evalii+1] = fpmresp_array[evalk*(nbz+1)*vsize+2*evalii+1];
            }


            for(evalmz=0; evalmz<nbz; evalmz++)
            {
                evalpha = -zonez_array[evalmz]*dphadz_array[evalk];
                evalcosp = cos(evalpha);
                evalsinp = sin(evalpha);
                evalki1 = evalki + (evalmz+1)*vsize;
                evalkv = evalk*vsize;

                for(evalii=0; evalii<vsize/2; evalii++)
                {
                    evalii1 = 2*evalii;
                    evalii2 = 2*evalii+1;
                    evalre = fpmresp_array[evalki1 + evalii1];
                    evalim = fpmresp_array[evalki1 + evalii2];
                    evalre1 = evalre*evalcosp - evalim*evalsinp;
                    evalim1 = evalre*evalsinp + evalim*evalcosp;
                    outtmp_array[evalkv + evalii1] += evalre1;
                    outtmp_array[evalkv + evalii2] += evalim1;
                }
            }

        }
        else  // single zone impulse
        {
            evalmz = focmMode-1;
            evalpha = zonez_array[evalmz]*dphadz_array[evalk];
            evalcosp = 1.0; //cos(evalpha);
            evalsinp = 0.0; //sin(evalpha);
            evalki1 = evalki + (evalmz+1)*vsize;
            evalkv = evalk*vsize;
            for(evalii=0; evalii<vsize/2; evalii++)
            {
                evalii1 = 2*evalii;
                evalii2 = 2*evalii+1;
                evalre = fpmresp_array[evalki1 + evalii1];
                evalim = fpmresp_array[evalki1 + evalii2];
                evalre1 = evalre*evalcosp - evalim*evalsinp;
                evalim1 = evalre*evalsinp + evalim*evalcosp;
                outtmp_array[evalkv + evalii1] = evalre1;
                outtmp_array[evalkv + evalii2] = evalim1;
            }

        }
    }


    //	for(evalmz=0; evalmz<nbz; evalmz++)
    //	outtmp_array[nbl*vsize + evalmz] = PIAACMC_MASKregcoeff*zonez_array[evalmz]*sqrt(vsize*nbl/nbz);


    evalval = 0.0;
    for(evalii=0; evalii<vsize*nbl; evalii++)
    {
        evalv1 = outtmp_array[evalii];
        evalval += evalv1*evalv1;
    }
    //  evalval /= vsize*nbl;

    // note that evalval is prop to bumber of spectral channels x number of evaluation pixels 
    return evalval;
}




double PIAACMCsimul_achromFPMsol_eval_zonezderivative(long zone, double *fpmresp_array, double *zonez_array, double *dphadz_array, double *outtmp_array, long vsize, long nbz, long nbl)
{

    // axis 0: eval pts (ii)   size = data.image[IDfpmresp].md[0].size[0] -> vsize
    // axis 1: zones (mz)      size = data.image[piaacmc[0].zonezID].md[0].size[0]+1 = nbz+1
    // axis 3: lambda (k)      size = piaacmc[0].nblambda -> nbl
    //
    // indexing :  k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*vsize + mz*vsize + ii



    for(evalk=0; evalk<nbl; evalk++) // lambda loop
    {
        evalki = evalk*(nbz+1)*vsize;


        // outer zone
        //  for(evalii=0; evalii<vsize; evalii++)
        //    outtmp_array[evalk*vsize+evalii] = 0.0; //fpmresp_array[evalk*(nbz+1)*vsize+evalii];


        evalmz = zone;
        
        // compute derivative as
        // dphadz is a function of wavelength
        // zonez is the current zone thickness
        // -zonez*dphadz sets the phase gives the resulting phase
        // Re = Re * -dphadz*sin(-zonez*dphadz) - Im * dphadz*cos(-zonez*dphadz)
        // Im = Re * dphadz*cos(-zonez*dphadz) + Im * -dphadz*sin(-zonez*dphadz)

        evalpha = -zonez_array[evalmz]*dphadz_array[evalk];
        // !!! note that cos is sin and sin is cos !!!
        // this implements a 90 degree pre-rotation so that this is a
        // derivative
        evalcosp = sin(evalpha)*dphadz_array[evalk]; //cos(evalpha);
        evalsinp = -cos(evalpha)*dphadz_array[evalk]; //sin(evalpha);
        evalki1 = evalki + (evalmz+1)*vsize;
        evalkv = evalk*vsize;

        for(evalii=0; evalii<vsize/2; evalii++)
        {
            evalii1 = 2*evalii;
            evalii2 = 2*evalii+1;
            evalre = fpmresp_array[evalki1 + evalii1];
            evalim = fpmresp_array[evalki1 + evalii2];
            evalre1 = evalre*evalcosp - evalim*evalsinp;
            evalim1 = evalre*evalsinp + evalim*evalcosp;
            outtmp_array[evalkv + evalii1] = evalre1;
            outtmp_array[evalkv + evalii2] = evalim1;
        }
    }

    return 0.0;
}











double f_evalmask (const gsl_vector *v, void *params)
{
    double *p = (double *)params;
    double value;
    long k;

    for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
        zonez_array[k] = gsl_vector_get(v, k);


    value = PIAACMCsimul_achromFPMsol_eval(fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
    value /= CnormFactor*SCORINGTOTAL*piaacmc[0].nblambda;

    if(LOOPCNT==0)
        cval0 = value;
    LOOPCNT++;

    return (value);

}






// remove outer zones to FPMresp
long PIAACMC_FPMresp_rmzones(char *FPMresp_in_name, char *FPMresp_out_name, long NBzones)
{
    long ID, IDout;
    long ii, jj, kk;
    long xsize = 0;
    long ysize = 0;
    long zsize = 0;
    long ysize1 = 0;

    ID = image_ID(FPMresp_in_name);
    xsize = data.image[ID].md[0].size[0];
    ysize = data.image[ID].md[0].size[1];
    zsize = data.image[ID].md[0].size[2];


    ysize1 = data.image[ID].md[0].size[1]-NBzones;

    IDout = create_3Dimage_ID_double(FPMresp_out_name, xsize, ysize1, zsize);

    for(ii=0; ii<xsize; ii++)
        for(kk=0; kk<zsize; kk++)
        {
            for(jj=0; jj<ysize1; jj++)
                data.image[IDout].array.D[kk*xsize*ysize1 + jj*xsize + ii] = data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
            for(jj=ysize1; jj<ysize; jj++)
                data.image[IDout].array.D[kk*xsize*ysize1 + ii] += data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
        }

    return(IDout);
}



// compress FPMresp to smaller number of lambda and points
long PIAACMC_FPMresp_resample(char *FPMresp_in_name, char *FPMresp_out_name, long NBlambda, long PTstep)
{
    long ID = -1;
    long IDout = -1;
    long ii, jj, kk, ii1, kk1;
    long xsize = 0;
    long ysize = 0;
    long zsize = 0;
    long xsize1 = 0;
    long zsize1 = 0;
    double alpha;
    long kk1xi;
    double kk1x;
    double re, im, re0, im0, re1, im1;


    ID = image_ID(FPMresp_in_name);
    xsize = data.image[ID].md[0].size[0];
    ysize = data.image[ID].md[0].size[1];
    zsize = data.image[ID].md[0].size[2];


    xsize1 = (long) (xsize/PTstep);
    zsize1 = NBlambda;


    IDout = create_3Dimage_ID_double(FPMresp_out_name, xsize1, ysize, zsize1);
    /*	for(kk1=0;kk1<zsize1;kk1++)
    		for(ii=0;ii<xsize;ii++)
    			data.image[IDout].array.D[kk1*xsize1*ysize + ii] = 1;
    	*/



    for(kk1=0; kk1<zsize1; kk1++)
    {
        kk1x = 1.0*kk1/(zsize1-1);
        kk1xi = (long) (kk1x*(zsize-1));
        alpha = kk1x*(zsize-1)-kk1xi;

        if(kk1xi==zsize-1)
        {
            kk1xi = zsize-2;
            alpha = 1.0;
        }

        printf("lambda index %ld  (%lf)   :  %ld %lf (%lf)\n", kk1, kk1x, kk1xi, alpha, (kk1xi+alpha)/(zsize-1));
        ii1 = 0;



        for(ii=0; ii<xsize; ii+=2*PTstep)
        {

            for(jj=0; jj<ysize; jj++)
            {
                re0 = data.image[ID].array.D[kk1xi*xsize*ysize + jj*xsize + ii];
                im0 = data.image[ID].array.D[kk1xi*xsize*ysize + jj*xsize + ii+1];

                re1 = data.image[ID].array.D[(kk1xi+1)*xsize*ysize + jj*xsize + ii];
                im1 = data.image[ID].array.D[(kk1xi+1)*xsize*ysize + jj*xsize + ii+1];

                re = (1.0-alpha)*re0 + alpha*re1;
                im = (1.0-alpha)*im0 + alpha*im1;


                data.image[IDout].array.D[kk1*xsize1*ysize + jj*xsize1 + ii1] = re; //*sqrt(PTstep);
                data.image[IDout].array.D[kk1*xsize1*ysize + jj*xsize1 + ii1+1] = im; //*sqrt(PTstep);
            }

            ii1+=2;
        }
    }


    /*		for(kk=0;kk<zsize;kk++)
    			{
    				for(jj=0;jj<ysize1;jj++)
    					data.image[IDout].array.D[kk*xsize*ysize1 + jj*xsize + ii] = data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
    				for(jj=ysize1;jj<ysize;jj++)
    					data.image[IDout].array.D[kk*xsize*ysize1 + ii] += data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
    			}
    */


    return(IDout);
}






/**
 *
 * @brief Main simulation routine
 *
 * @param[in] confindex  PIAACMC configuration index pointing to the input/output directory number
 * @param[in] mode       Type of operation to be performed
 *
 */

int PIAACMCsimul_exec(char *confindex, long mode)
{
    long NBparam;
    FILE *fp;
    FILE *fpt;
    FILE *fptest;
    char fptestname[1000];
    char command[1000];
    char dirname[500];
    int paramtype[10000]; // FLOAT or DOUBLE
    double *paramval[10000]; // array of pointers, double
    float *paramvalf[10000]; // array of pointers, float
    double paramrefval[10000];

    double paramdelta[10000];
    double parammaxstep[10000]; // maximum single iteration step
    double parammin[10000]; // minimum value
    double parammax[10000]; // maximum value

    double paramdeltaval[10000];

    double valref, valref1, valbest, val0;
    double parambest[10000]; // for scanning
    double paramref[10000];
    long i, ii, jj, kk;
    long IDv, ID, ID1, IDref, IDa, IDcomb;
    long IDmodes, IDmodes2D;
    long xsize = 0;
    long ysize = 0;
    long zsize = 0;
    long k;

    long iter;
    long NBiter = 1000;

    long IDfpmresp, IDref1;
    double t, a, dpha, amp;
    int zi;

    char fname[800];
    char fnamecomb[500];
    char fname1[500];
    char fname2[500];
    char fnamet[500];
    char fnametmp[500];
    char fnametransm[500];
    char fnamelog[500];
    long IDm, ID1D, ID1Dref;
    long size1Dvec, size1Dvec0;

    char fnamea[500];
    char fnamep[500];
    long elem0;


    // OPTIMIZATION PARAMETERS
    int REGPIAASHAPES = 0;

    float piaa0C_regcoeff = 0.0e-7; // regularization coeff
    float piaa1C_regcoeff = 0.0e-7; // regularization coeff
    float piaa0C_regcoeff_alpha = 1.0; // regularization coeff power
    float piaa1C_regcoeff_alpha = 1.0; // regularization coeff power

    float piaa0F_regcoeff = 0.0e-7; // regularization coeff
    float piaa1F_regcoeff = 0.0e-7; // regularization coeff
    float piaa0F_regcoeff_alpha = 1.0; // regularization coeff power
    float piaa1F_regcoeff_alpha = 1.0; // regularization coeff power


    int REGFPMSAG = 0; // regularization for FPM sag
 //   float fpmsagreg_coeff = 1.0e-8;
//    float fpmsagreg_coeff_alpha = 1.0;




    int r;

    double val, v1, v2, pha, cosp, sinp, re, im, re1, im1;


    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    double range, stepsize;
    int loopOK;
    long ls, ls1, NBpropstep;
    double lstransm;
    long mz;

    int LINOPT = 0; // 1 if perform linear optimization
    long ii1, ii2, ki, kv, ki1;

    double scangain;
    double scanstepgain = 0.001;
    int linscanOK;
    double valold, oldval;
    double bestgain;
    double linoptgainarray[100];
    double linoptvalarray[100];
    int linoptlimflagarray[100];
    long NBlinoptgain;
    long kmax;


    double valtest;
    int ret;
    int kmaxC, kmaxF;


    long st;
    long NBoptVar;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function fpmeval_func;
    int status;

    // simulated annealing loop(s)
    long NBITER_SA = 1000000;
    long ITERMAX;
    long iter1;
    double SAcoeff;
    double amp1;
    double val1;
    int SA_MV = 0;
    double tmp1;
    double bestvalue;
    FILE *fpbest;

    // downhill simplex
    long NBITER_DS = 1000;
    long iterMax;
    int OK;
    int KEEP;
    double KEEPlimit = 1.0;
    double KEEPlimit1 = 1.0;
    double eps1;
    double mmsize;

    long elem;


    double tmplf1, tmplf2;
    float *peakarray;
    float avpeak;

    double *eval_contrastCurve;
    double *eval_contrastCurve_cnt;
    double eval_sepstepld = 0.2; // in l/D
    double eval_sepmaxld = 20.0; // in l/D
    long eval_sepNBpt;
    double focscale;
    double xc, yc, rc;
    long ri;
    long IDps;
    double xld, stepld;
    double ldoffset;

    // spreading computations over multiple processes for resp matrix
    long mzoffset = 0;
    long mzstep = 1;
    long thr;
    long tmpl1;
    int tmpd1;

    int iterOK;


    char stopfile[500];
    struct stat buf;

    double aveC;
    long aveCcnt;

    double alphareg;
    double bestval = 0.0;
    long IDoptvec = -1;

    long index;

    long tmpNBrings;
    int tmpnblambda;
    float dx, dy;

    double zmin, zmax, sigma;
    long IDc;
    double oaoffset; // off-axis offset for Lyot stop optimization
    long k1, kr;
    long NBincpt, NBincpt1, NBkr;
    long cnt;

    int PIAACMC_WFCmode = 0;

    double xpos, ypos, fval;
    int initscene;
    long IDscene;

    double acoeff0, acoeff1, acoeff2;
    float scangainfact;
    int alphascaninit = 0;


    long *sizearray;
    long IDstatus;

    double valContrast;
    double tmp;
    int initbestval = 0;


	long size;
	long IDopderrC, nbOPDerr, OPDmode, IDopderr;

	int tmpn1;
	char userinputstr[100];
	

    // Create status shared variable
    // this allows realtime monitoring of the code by other processes
    // sets status at different points in the code
    IDstatus = image_ID("stat_PIAACMCsimulexec");
    if(IDstatus == -1)
        IDstatus = read_sharedmem_image("stat_PIAACMCsimulexec");

    sizearray = (long*) malloc(sizeof(long)*2);
    sizearray[0] = 1;
    sizearray[1] = 1;
    IDstatus = create_image_ID("stat_PIAACMCsimulexec", 2, sizearray, USHORT, 1, 0);
    free(sizearray);



    piaacmc = NULL; // set the pointer to the piaacmc structure to null

    // if the optical system pointer is empty, create an empty version
    if(optsyst==NULL)
    {
        optsyst = (OPTSYST*) malloc(sizeof(OPTSYST));
        optsyst[0].SAVE = 1;
    }

    for(elem=0; elem<100; elem++)
        optsyst[0].keepMem[elem] = 0; // flag that says save this element for reuse


    // set the result directories
    sprintf(piaacmcconfdir, "%s", confindex);
    sprintf(data.SAVEDIR, "%s", piaacmcconfdir);



    optsyst[0].SAVE = PIAACMC_save;

    // get variables from command line, possibly sets globals
    if((IDv=variable_ID("PIAACMC_centobs0"))!=-1)
        centobs0 = data.variable[IDv].value.f;
    if((IDv=variable_ID("PIAACMC_centobs1"))!=-1)
        centobs1 = data.variable[IDv].value.f;
    if((IDv=variable_ID("PIAACMC_fpmradld"))!=-1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    // start a log of code mode entry/exit times
    sprintf(command, "echo \"%03ld     $(date)\" >> ./log/PIAACMC_mode_log.txt", mode);
    ret = system(command);
    printf("command = %s\n", command);


    // set the name of the stopfile
    sprintf(stopfile, "%s/stopmode%ld.txt", piaacmcconfdir, mode);

    switch (mode) {
    case 0 :  // Run existing config for on-axis point source. If new, create centrally obscured idealized PIAACMC
        // compatible with wavefront control
        printf("=================================== mode 000 ===================================\n");
        // Either load a set of point sources from "scene.txt" or use a single on-axis point source,
        // and create the image for these sources by computing and adding their PSFs
            
        // load some more cli variables
        PIAACMC_fpmtype = 0; // idealized (default)
        if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
            PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);
        printf("PIAACMC_fpmtype = %d\n", PIAACMC_fpmtype);

        PIAACMC_WFCmode = 0; // number of DMs
        if((IDv=variable_ID("PIAACMC_WFCmode"))!=-1)
            PIAACMC_WFCmode = (int) (data.variable[IDv].value.f + 0.1);
        printf("PIAACMC_WFCmode = %d\n", PIAACMC_WFCmode);

        // force creation of the FPM zone amplitudes by called functions
        FORCE_CREATE_fpmza = 1;

        // main initialization function to set up the piaacmc structure
        PIAAsimul_initpiaacmcconf(PIAACMC_fpmtype, fpmradld, centobs0, centobs1, PIAACMC_WFCmode, 1);
        // make the mirror or lenses shapes (only mirrors for WFIRST piaacmc)
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm normalization for efficiency


        // if file "scene.txt" exists, compute series of PSFs and sum
        fp = fopen("SCENE.txt", "r");
        if(fp!=NULL)
        {
            initscene = 0;
            // for each source in the scene, read position and flux
            while(fscanf(fp, "%lf %lf %lf\n", &xpos, &ypos, &fval) == 3)
            {
                printf("COMPUTING PSF AT POSITION %lf %lf, flux  = %g\n", xpos, ypos, fval);
                // make the actual PSF
                PIAACMCsimul_computePSF(xpos, ypos, 0, optsyst[0].NBelem, 1, 0, 0, 1);
                // get the image "psfi0" index, which was created in PIAACMCsimul_computePSF
                ID = image_ID("psfi0");
                // get image size.  3rd dimension is wavelength
                xsize = data.image[ID].md[0].size[0];
                ysize = data.image[ID].md[0].size[1];
                zsize = data.image[ID].md[0].size[2];

                if(initscene==0)
                {
                    initscene = 1;
                    // create 3D image to sum the PSFs into
                    IDscene = create_3Dimage_ID("scene", xsize, ysize, zsize);
                }
                ID = image_ID("psfi0");
                // sum the current PSF into the image: summed image is IDscene, source is ID
                for(ii=0; ii<xsize*ysize*zsize; ii++)
                    data.image[IDscene].array.F[ii] += fval*data.image[ID].array.F[ii];


            }
            fclose(fp);
            // we're done!  Save it, overwriting previous scene.fits file
            save_fits("scene", "!scene.fits");
        }
        else // scene.txt does not exist, just do an on-axis source
        {
            valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 1);
            printf("valref = %g\n", valref);
        }

        printf("EXEC CASE 0 COMPLETED\n");
        fflush(stdout);

        break;












    case 1 : // optimize Lyot stop positions
        // Lyot stop positions are encoded as piaacmc[0].LyotStop_zpos
        // there can be multiple LyotStop_zpos
        // Vary these zpos, looking for the best contrast returned by PIAACMCsimul_computePSF
        // Search is performed by iterative refined marching
        printf("=================================== mode 001 ===================================\n");

        // init as in mode 0
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        // initialization
		PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0); // necessary to initialize optical design
        // set initial lyot stop marching range (current position +- range)
        if((IDv=variable_ID("PIAACMC_lsoptrange"))!=-1)
            range = data.variable[IDv].value.f; // from cli
        else
            range = 3.0; // default, in meters
        stepsize = range/3.0; // initial march stepsize
        // store initial Lyot stop positions
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++) // NBLyotStop = length(LyotStop_zpos)
            paramref[ls] = piaacmc[0].LyotStop_zpos[ls];
        NBiter = 4; // number of iterations

        // start up a log
        sprintf(fnamelog, "%s/result_LMpos.log", piaacmcconfdir);
        fp = fopen(fnamelog, "w");
        fclose(fp);



        // pick another initial march stepsize
        stepsize = range/5.0;
        // start the iterative refined march
        for(iter=0; iter<NBiter; iter++)
        {
            // for each Lyot stop, find its best position
            for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            {
                // start position for march.  paramref is current best value
                piaacmc[0].LyotStop_zpos[ls] = paramref[ls]-range;
                // current best position
                parambest[ls] = piaacmc[0].LyotStop_zpos[ls];

                loopOK = 1;
                valbest = 1.0;
                
                // march to the other other end of range
                while(piaacmc[0].LyotStop_zpos[ls]<paramref[ls]+range)
                {
                    elem0 = 6; // elem0 is the starting point of the PSF propagation.  This is a staring default

                    // look for the element called "Lyot mask 0" as the actual starting point
					elem0 = -1;
					printf("Number of elements = %ld\n", optsyst[0].NBelem);
					assert ( optsyst[0].NBelem > 0 );
					for(elem=0; elem<optsyst[0].NBelem; elem++)
					{
						printf("elem %ld :  %s\n", elem, optsyst[0].name[elem]);
                        if(strcmp("Lyot mask 0", optsyst[0].name[elem])==0)
                            elem0 = elem;
					}
					assert( elem0 != -1); // throw a message if this was not found


                    optsyst[0].keepMem[elem0] = 1; // save this element and reuse

                    // compute the PSF for this Lyot stop position, returning contrast in the evaluation zone
                    val = PIAACMCsimul_computePSF(0.0, 0.0, elem0, optsyst[0].NBelem, 0, 0, 0, 0);

                    // if this is the best contrast for this stop, save it for it and the position of this stop
                    if(val<valbest)
                    {
                        parambest[ls] = piaacmc[0].LyotStop_zpos[ls];
                        valbest = val;
                    }

                    // say what's happening
                    fp = fopen(fnamelog, "a");
                    for(ls1=0; ls1<piaacmc[0].NBLyotStop; ls1++)
                        fprintf(fp," %lf", piaacmc[0].LyotStop_zpos[ls1]);
                    fprintf(fp, " %g\n", val);
                    fclose(fp);

                    // march along by the step size
                    piaacmc[0].LyotStop_zpos[ls] += stepsize;
                }
                printf("BEST SOLUTION :  ");
                paramref[ls] = parambest[ls]; // update best position for this stop
                piaacmc[0].LyotStop_zpos[ls] = paramref[ls]; // store in case this is last iteration
                printf(" %lf", parambest[ls]);
                printf(" %g\n", valbest);
            }

            fp = fopen(fnamelog, "a");
            fprintf(fp, "\n");
            fclose(fp);

            // reduce the range and stepsize, refining the march
            range *= 0.3;
            stepsize = range/3.0;
        }
        // store all best positions  Done!!
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            piaacmc[0].LyotStop_zpos[ls] = parambest[ls];
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);
        break;






    case 2 : // optimize focal plane mask transmission for monochromatic idealized PIAACMC
        // for monochromatic, idealized PIAACMC, find the scalar transimssion of the uniform focal plane mask
        // that provides best contrast in the evaluation zone
        // very similar to the Lyot stop search in mode 1: iterative refined marching, changing the
        // the transmission value piaacmc[0].fpmaskamptransm, which
        // is between 0 and 1 (with Olivier's sign convention)
        // uses single on-axis light source
        printf("=================================== mode 002 ===================================\n");

        // init as in mode 0
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        // initialization, see mode 1
        range = 0.3;
        stepsize = range/3.0;
        paramref[0] = piaacmc[0].fpmaskamptransm; // initialized in PIAAsimul_initpiaacmcconf
        NBiter = 6;

        sprintf(fnamelog, "%s/result_fpmt.log", piaacmcconfdir);
        fp = fopen(fnamelog, "w");
        fclose(fp);

        for(iter=0; iter<NBiter; iter++)
        {
            // starting point of march
            piaacmc[0].fpmaskamptransm = paramref[0]-range;
            // store current value as best
            parambest[0] = piaacmc[0].fpmaskamptransm;

            loopOK = 1;
            valbest = 1.0;

            // while within the march range
            while(loopOK==1)
            {
                FORCE_CREATE_fpmza = 1; // forces creation of new focal plane mask in the next two routines
                PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0);
                // compute on-axis PSF of all optical elements returning contrast in evaluation zone
                // ************************* need to do all optsyst[0].NBelem?
                val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);

                if(val<valbest)
                {
                    // we have a better contrast!  Store it
                    parambest[0] = piaacmc[0].fpmaskamptransm;
                    valbest = val;
                }

                fp = fopen(fnamelog, "a");
                fprintf(fp," %lf", piaacmc[0].fpmaskamptransm);
                fprintf(fp, " %g  %ld %g %g\n", val, iter, range, stepsize);
                fclose(fp);

                ls = 0; // probably a copy and paste from mode 1 ************************************
                // keep marching along
                piaacmc[0].fpmaskamptransm += stepsize;
                // if we've reached the end of the range stop the loop
                if(piaacmc[0].fpmaskamptransm>paramref[0]+range+0.001*stepsize)
                    loopOK = 0;
            }

            printf("BEST SOLUTION :  ");
            // store best solution
            paramref[0] = parambest[0];
            printf(" %lf", parambest[0]);

            printf(" %g\n", valbest);


            fp = fopen(fnamelog, "a");
            fprintf(fp, "\n");
            fclose(fp);
            // refine range and stepsize
            range *= 0.3;
            stepsize = range/3.0;
        }
        // save final result
        piaacmc[0].fpmaskamptransm = parambest[0];
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0); // why? **************************
        PIAAsimul_savepiaacmcconf(piaacmcconfdir); // save final result to disk
        FORCE_CREATE_fpmza = 0; // turning off to be good citizens
        break;






    case 3 : // calibrate, no focal plane mask (not currently used)
        // compute PSF and contrast with no focal plane mask with the current design
        // provides the denominator for the contrast estimate
        // saved by PIAACMCsimul_computePSF as fits file "psfi0"
        printf("=================================== mode 003 ===================================\n");

        // init as in mode 0
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        paramref[0] = piaacmc[0].fpmaskamptransm;

        piaacmc[0].fpmaskamptransm = -1.0;  // Remove focal plane mask
        FORCE_CREATE_fpmza = 1;
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0);
        // compute the PSF for an on-axis source, all optical elements
        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);

        // restore original configuration
        piaacmc[0].fpmaskamptransm = paramref[0];
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0);
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);
        FORCE_CREATE_fpmza = 0;

        break;





    case 4 : // optimize PIAA optics shapes, cosine modes only (not currently used, replaced by mode 40. skipping)
        printf("=================================== mode 004 ===================================\n");

        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        LINOPT = 1; // perform linear optimization
        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            NBiter = (long) data.variable[IDv].value.f+0.01;
        else
            NBiter = 1000;

        kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
            kmax = (long) data.variable[IDv].value.f+0.01;

        if(kmax>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
            kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];

        NBparam = 0;
        for(k=0; k<kmax; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-5;
            parammax[NBparam] = 1.0e-5;
            NBparam++;
        }

        for(k=0; k<kmax; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-5;
            parammax[NBparam] = 1.0e-5;
            NBparam++;
        }
        FORCE_MAKE_PIAA0shape = 1;
        FORCE_MAKE_PIAA1shape = 1;
        break;








    case 5 : // optimize Lyot stops shapes and positions
        printf("=================================== mode 005 ===================================\n");
        // init as in mode 0
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        // load cli variables as appropriate
        // # of propagation steps along the beam
        NBpropstep = 150;
        if((IDv=variable_ID("PIAACMC_nbpropstep"))!=-1)
            NBpropstep = (long) data.variable[IDv].value.f+0.01;
        // desired Lyot stop transmission
        lstransm = 0.85;
        if((IDv=variable_ID("PIAACMC_lstransm"))!=-1)
            lstransm = (double) data.variable[IDv].value.f;
        printf("lstransm  = %f\n", lstransm);

        /// identify post focal plane pupil plane (first pupil after focal plane mask)
        // provides reference complex amplitude plane for downstream analysis
        PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);
        printf("=========== %ld elements ======================================================\n", optsyst[0].NBelem);
        // find the ID of the "post focal plane mask pupil" element
        for(elem=0; elem<optsyst[0].NBelem; elem++)
        {
            if(strcmp("post focal plane mask pupil", optsyst[0].name[elem])==0)
            {
                elem0 = elem;
                printf("post focal plane mask pupil = %ld\n", elem);
            }
            else
                printf("elem %ld : %s\n", elem, optsyst[0].name[elem]);
        }
        optsyst[0].keepMem[elem0] = 1; // keep it for future use

        oaoffset = 20.0; // off axis amplitude
        // compute the reference on-axis PSF
        PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);

        // filenames of the complex amplitude and phase in the post FPM pupil plane indexed by elem0
        sprintf(fnamea, "WFamp0_%03ld", elem0);
        sprintf(fnamep, "WFpha0_%03ld", elem0);

        printf("elem0 = %ld\n", elem0);

        // args for the PIAACMCsimul_CA2propCubeInt function
       
        // set range of propagation
        zmin = piaacmc[0].LyotZmin;
        zmax = piaacmc[0].LyotZmax;
            
        // args that determine "extended" off-axis source
        // number of off-axis sources on each circle
        NBincpt = 15;
        // number of circle radii
        NBkr = 5;

        // propagate complex amplitude in a range from zmin to zmax, where 0 is elem0
        // computes the diffracted light from the on-axis source
        ID1 = PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iproptmp");
        // complex amplitude at elem0, only used to determine image size
        IDa = image_ID(fnamea);

        xsize = data.image[IDa].md[0].size[0];
        ysize = data.image[IDa].md[0].size[1];

        // OAincohc is the summed light "all" from off-axis sources in the pupil,
        // including the on-axis source(!),
        // giving intensity contribution of all off-axis sources
        // in order to preserve the intensity of the off-axis in the design.
        // load OAincohc if exist, maybe we've been here before
        sprintf(fname, "%s/OAincohc.fits", piaacmcconfdir);
        IDc = load_fits(fname, "OAincohc", 1);


        if(IDc==-1) // OAincohc does not exist so we have to make it
        {
            // create image to receive sum
            IDc = create_3Dimage_ID("OAincohc", xsize, ysize, NBpropstep);

            cnt = 0; // initialize counter so later we can normalize by number of sources
            // loop over radii
            for(kr=0; kr<NBkr; kr++)
            {
                // loop over points at current radius
                for(k1=0; k1<NBincpt; k1++)
                {
                    // compute PSF for a point at this angle with scaled offset
                    // PIAACMCsimul_computePSF changes fnamea and fnamep (in call to OptSystProp_run)!
                    PIAACMCsimul_computePSF(oaoffset*(1.0+kr)/NBkr*cos(2.0*M_PI*k1/NBincpt), oaoffset*(1.0+kr)/NBkr*sin(2.0*M_PI*k1/NBincpt), 0, optsyst[0].NBelem, 0, 0, 0, 0);
                    // propagate that elem0 from zmin to zmax with new PSF
                    ID1 = PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iproptmp");
                    for(ii=0; ii<xsize*ysize; ii++) // ii is indexing x-y plane
                        for(k=0; k<NBpropstep; k++) // k is indexing z-direction. adding to IDc which looks right
                            data.image[IDc].array.F[k*xsize*ysize+ii] += data.image[ID1].array.F[k*xsize*ysize+ii];
                    delete_image_ID("iproptmp");
                    cnt ++;
                }
            }
            // scale by the number of sources to give average
            for(ii=0; ii<xsize*ysize; ii++)
                for(k=0; k<NBpropstep; k++)
                    data.image[IDc].array.F[k*xsize*ysize+ii] /= cnt;


            sprintf(fname, "!%s/OAincohc.fits", piaacmcconfdir);
            // save the final result
            save_fits("OAincohc", fname);
        }
        // compute on-axis PSF to define light to reject
        PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        // propagate it into the optical system, with result in image named "iprop00"
        PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iprop00");
        //  save_fits("iprop00", "!test_iprop00.fits");

        // make image that has the min along z of OAincohc at each x,y
        PIAACMCsimul_optimizeLyotStop_OAmin("OAincohc");
        // do the actual Lyot stop shape and location optimization, producing optimal Lyot stops in optLM*.fits
        // and position relative to elem0 in piaacmc[0].LyotStop_zpos
        PIAACMCsimul_optimizeLyotStop(fnamea, fnamep, "OAincohc", zmin, zmax, lstransm, NBpropstep, piaacmc[0].NBLyotStop);

        sprintf(fptestname, "conj_test.txt");
        fptest = fopen(fptestname, "w");
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            fprintf(fptest, "%ld  %f  %f     %f  %f\n", ls, zmin, zmax, piaacmc[0].LyotStop_zpos[ls], optsyst[0].elemZpos[elem0]);
        fclose(fptest);

        // convert Lyot stop position from relative to elem0 to absolute
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            piaacmc[0].LyotStop_zpos[ls] += optsyst[0].elemZpos[elem0];

        // and we're done!  save.
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);
        // copy to the final Lyot stop file for this mode
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
        {
            sprintf(command, "cp ./%s/optLM%02ld.fits ./%s/LyotStop%ld.fits", piaacmcconfdir, ls, piaacmcconfdir, ls);
            r = system(command);
        }

        break;




 




    case 11 : // setup multizone ring mask and Compute polychromatic response to zones, store result in FPMresp
        // here we compute how the light propagates from each individual mask zone to the focal plane
        // (where each mask zone is completely tranparent)
        printf("=================================== mode 011 ===================================\n");

        printf("STEP01\n");

        // get cli variables
        if((IDv=variable_ID("PIAACMC_nblambda"))!=-1)
            tmpnblambda = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_NBrings"))!=-1)
            tmpNBrings = data.variable[IDv].value.f;


        // initialize
        PIAAsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1);


        printf("piaacmcconfdir     : %s\n", piaacmcconfdir);
        fflush(stdout);
        printf("SCORINGMASKTYPE    : %d\n", SCORINGMASKTYPE);
        fflush(stdout);
        printf("PIAACMC_FPMsectors : %d\n", PIAACMC_FPMsectors);
        fflush(stdout);
        printf("lamda              : %ld nm\n", (long) (1.0e9*piaacmc[0].lambda + 0.1));
        fflush(stdout);
        printf("lamdaB             : %ld \n", (long) (1.0*piaacmc[0].lambdaB + 0.1));
        fflush(stdout);
        printf("piaacmc[0].NBrings : %ld\n", piaacmc[0].NBrings);
        fflush(stdout);
        printf("mask rad           : %ld\n", (long) (100.0*PIAACMC_MASKRADLD+0.1));
        fflush(stdout);
        printf("computePSF_ResolvedTarget : %d\n", computePSF_ResolvedTarget);
        fflush(stdout);
        printf("computePSF_ResolvedTarget_mode : %d\n", computePSF_ResolvedTarget_mode);
        fflush(stdout);
        printf("piaacmc[0].fpmmaterial_name : %s\n", piaacmc[0].fpmmaterial_name);
        fflush(stdout);
        printf("piaacmc[0].nblambda         : %d\n", piaacmc[0].nblambda);
        fflush(stdout);

        // set output filename of the combined focal plane mask response file
       // sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ccnbr%03ld_ccz%06ld_ocr%04ld_ocz%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBringCentCone, (long) (1.0e9*piaacmc[0].fpmCentConeZ+0.1), (long) (100.0*piaacmc[0].fpmOuterConeRadld+0.1), (long) (1.0e9*piaacmc[0].fpmOuterConeZ+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

		sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
  
        printf("fname = %s\n", fname);
        fflush(stdout);


        // get the combined focal plane mask response
        ID = load_fits(fname, "FPMresp", 1);
        // if it did not exist, create it
        if(ID==-1)
        {

            // get the number of tmux threads from cli
            PIAACMC_FPMresp_mp = 1; // 1: all computations on a single thread
            if((IDv=variable_ID("PIAACMC_FPMresp_mp"))!=-1) // multi threaded
                PIAACMC_FPMresp_mp = (long) data.variable[IDv].value.f+0.01;
            printf("PIAACMC_FPMresp_mp = %ld\n", PIAACMC_FPMresp_mp);

            printf("------------------------------------- STEP02\n");
            printf("piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
            fflush(stdout);

            // get our tmux thread number in [0 PIAACMC_FPMresp_mp]
            // where the master thread has PIAACMC_FPMresp_thread == PIAACMC_FPMresp_mp
            PIAACMC_FPMresp_thread = 0;
            if((IDv=variable_ID("PIAACMC_FPMresp_thread"))!=-1) // multi threaded
                PIAACMC_FPMresp_thread = (long) data.variable[IDv].value.f+0.01;
            printf("PIAACMC_FPMresp_thread = %ld\n", PIAACMC_FPMresp_thread);


            
            index = 0;
            if((PIAACMC_FPMresp_mp==1)||(PIAACMC_FPMresp_thread>PIAACMC_FPMresp_mp-1))  // main or combine process
                                            // why not test PIAACMC_FPMresp_thread == PIAACMC_FPMresp_mp?
            {
                // we're the parent set up the FPM zone map
                FORCE_CREATE_fpmzmap = 1;
                FORCE_CREATE_fpmzt = 1;
                FORCE_CREATE_fpmza = 1;
                PIAAsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1);
            }
            else
            {
                printf("NO initOK file created\n");
                // we're a child tmux thread, do not set up the FPM zone map, get it from parent via file
                FORCE_CREATE_fpmzmap = 0;
                FORCE_CREATE_fpmzt = 0;
                FORCE_CREATE_fpmza = 0;
                PIAAsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1);
            }



            // if we're the parent load
            if((PIAACMC_FPMresp_mp==1)||(PIAACMC_FPMresp_thread>PIAACMC_FPMresp_mp-1))
            {
//                sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ccnbr%03ld_ccz%06ld_ocr%04ld_ocz%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBringCentCone, (long) (1.0e9*piaacmc[0].fpmCentConeZ+0.1), (long) (100.0*piaacmc[0].fpmOuterConeRadld+0.1), (long) (1.0e9*piaacmc[0].fpmOuterConeZ+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

                sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

                sprintf(fname1, "!%s.tmp", fname);
                sprintf(fname2, "!%s", fname);
                mzoffset = 0;
                mzstep = 1;

                ID = load_fits(fname, "FPMresp", 1);  // this will always fail in the current state (see line 6606) ************************
                IDcomb = ID;
            }
            else // we're a child tmux thread.
            {
                // combined FPMresp file
//                sprintf(fnamecomb, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ccnbr%03ld_ccz%06ld_ocr%04ld_ocz%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBringCentCone, (long) (1.0e9*piaacmc[0].fpmCentConeZ+0.1), (long) (100.0*piaacmc[0].fpmOuterConeRadld+0.1), (long) (1.0e9*piaacmc[0].fpmOuterConeZ+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
               sprintf(fnamecomb, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);


                // partial FPMresp file
//                sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ccnbr%03ld_ccz%06ld_ocr%04ld_ocz%06ld_ssr%02d_ssm%d_%s_wb%02d_mp%02ld_thread%02ld.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBringCentCone, (long) (1.0e9*piaacmc[0].fpmCentConeZ+0.1), (long) (100.0*piaacmc[0].fpmOuterConeRadld+0.1), (long) (1.0e9*piaacmc[0].fpmOuterConeZ+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda, PIAACMC_FPMresp_mp, PIAACMC_FPMresp_thread);
               sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d_mp%02ld_thread%02ld.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda, PIAACMC_FPMresp_mp, PIAACMC_FPMresp_thread);

                // stash the filename of the partial file for later
                sprintf(fname1, "!%s.tmp", fname);
                sprintf(fname2, "%s", fname);
                // set region of the partial file that this child computes
                mzoffset = PIAACMC_FPMresp_thread;
                mzstep = PIAACMC_FPMresp_mp;

                ID = load_fits(fname, "FPMresp", 1); // may exist from a previous execution with restart
                IDcomb = load_fits(fnamecomb, "FPMresp", 1); // will always fail
            }
            // at this point IDcomb==-1, and in the parent ID==-1 always, and in the child ID==-1 if this is not a restart
            // actually create the FPMresp file either as a part by a child or combined by the parent
            if((IDcomb==-1)&&(ID==-1))  // this will always fire for the parent thread,
                                            // and will always fire for children in a fresh run
            {
                //                printf("--------------------------------------------------------STEP 0005 File \"%s\" does not exist: creating\n", fname);
                //   printf("piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
                //  sleep(3);

                // usual initialzation
                optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                //piaacmc[0].fpmaskamptransm = 1.0;
                // set the physical size of the FPM as mean(lambda/D)*mask radius in units of lambda/D
                piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio * PIAACMC_MASKRADLD; // PIAACMC_MASKRADLD l/D radius at central lambda
                // initialize the optical system
                PIAAsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 0);

                //     printf("-------------------------- STEP 0005a  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
                //            sleep(3);
                
                // computes or loads the piaa optics from the piaacmc structure
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                //   printf("-------------------------- STEP 0005b  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
                //          sleep(3);
                
                // initialize the optical system to be on axis
                PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);
                // printf("-------------------------- STEP 0005c  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
                //        sleep(3);

                // make the shapes again why?  *****************************************
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);

                /*               printf("piaacmc[0].NBrings =                             %ld\n", piaacmc[0].NBrings);
                                printf("piaacmc[0].focmNBzone =                          %ld\n", piaacmc[0].focmNBzone);
                                printf("piaacmc[0].nblambda =                            %d\n", piaacmc[0].nblambda);

                                fflush(stdout);
                sleep(5);*/
                
                
                // focmMode controls which part of the FPM is propagated
                // if focmMode is a legal zone index, that zone index is propagated
                // here set focmMode beyond a legal zone index, so all zones are transparent and all
                // light including that which misses the FPM is propagated.
                // Later, we will subtract off the actual zone contributions, which will leave only
                // the light that misses the FPM.
                focmMode = data.image[piaacmc[0].zonezID].md[0].size[0]+10;  // response for no focal plane mask
                optsyst[0].FOCMASKarray[0].mode = 1; // 1-fpm
                // compute the on-axis PSF to see what on-axis light goes around the FPM, return contrast in evaluation zone
                val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, 0);
                printf("val = %g\n", val);
                ID = image_ID("imvect");

                // FPMresp geometry:
                // first dimension (size[0]) is twice the number of evaluation points in the focal plane, giving Re and Im
                //      of the field at that evaluation point
                // second dimension (size[1]) is nbzones+1 zone indices, where nbzones is the number of mask zones (hexagons)
                // +1 because the first zone index stores the response to light that misses the FPM
                // WARNING: FPMresp size[1] is nbzones+1, as first vector stored is the response for light outside the mask
                // third dimension (size[2]) is wavelength


                // axis 0: eval pts (ii) - size = data.image[ID].md[0].size[0]
                // axis 1: zones (mz) - size = data.image[piaacmc[0].zonezID].md[0].size[0]+1
                // axis 3: lambda (k) - size = piaacmc[0].nblambda

                // allocate the combined FPMresp 3D array
                // ID is the "imvect" array created by PIAACMCsimul_computePSF and contains the pixels in the
                // evaluation set as a 1D vector (0th index) per wavelength
                IDfpmresp = create_3Dimage_ID_double("FPMresp", data.image[ID].md[0].size[0], piaacmc[0].focmNBzone+1, piaacmc[0].nblambda);
                //     list_image_ID();
                //    sleep(100);

                // light outside mask
                for(k=0; k<piaacmc[0].nblambda; k++) // loop over wavelengths
                    for(ii=0; ii<data.image[ID].md[0].size[0]; ii++) // loop over evaluation points
                        // set the 0th zone to be the light from the above on-axis PSF computation with a
                        // black FPM on the evaluation pixels in "imvect"
                        data.image[IDfpmresp].array.D[k*(piaacmc[0].focmNBzone+1)*data.image[ID].md[0].size[0] + ii] = data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];


                if(PIAACMC_FPMresp_thread>PIAACMC_FPMresp_mp-1) // if we're the parent, combine files
                {
                    if((IDv=variable_ID("PID"))!=-1)
                        index = (long) data.variable[IDv].value.f+0.01;
                    // this file is looked for in the bash script, which waits for this
                    // file to spawn the tmux child processes
                    sprintf(command, "touch initOK_%ld", index);
                    printf("EXECUTING : %s\n", command);
                    r = system(command);
                    
                    // now the tmux children have been kicked off by the bash script and are
                    // computing their partial FPMresp files

                    // begin combining the partial files into the final FPMresp file when ready
                    printf("COMBINING FILES\n");
                    fflush(stdout);
                    sprintf(fnamet, "%s/FPMthreadstatus.txt", piaacmcconfdir);
                    fpt = fopen(fnamet, "w");
                    fclose(fpt);

                    // we now wait for the children to complete their partial FPM resp files, then combine them
                    for(thr=0; thr<PIAACMC_FPMresp_mp; thr++)
                    {
                        printf("thr = %ld\n", thr);
                        fflush(stdout);
                        
                        // each child creates an FPMresp...thread*.fits.tmp file, which is moved to
                        // FPMresp...thread*.fits (set as fname in the next sprintf) when the child is done
                        // signaling to the parent process that this part is ready to ingest.

                        ID1 = -1;
                        while(ID1 == -1) // wait for the partial FPMresp file from each child
                        {
                            // name of final child partial FPMresp file
                            //sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ccnbr%03ld_ccz%06ld_ocr%04ld_ocz%06ld_ssr%02d_ssm%d_%s_wb%02d_mp%02ld_thread%02ld.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBringCentCone, (long) (1.0e9*piaacmc[0].fpmCentConeZ+0.1), (long) (100.0*piaacmc[0].fpmOuterConeRadld+0.1), (long) (1.0e9*piaacmc[0].fpmOuterConeZ+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda, PIAACMC_FPMresp_mp, thr);
                            sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d_mp%02ld_thread%02ld.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda, PIAACMC_FPMresp_mp, thr);
                             
                            printf("Waiting for file \"%s\" ...\n", fname);
                            fflush(stdout);
                            // update thread status file
                            fpt = fopen(fnamet,"a");
                            fprintf(fpt, "Process %ld (thread %ld) --- Waiting for file \"%s\" ...\n", (long) getpid(), thr, fname);
                            fclose(fpt);
                            sleep(1.0);

                            // safely remove image with this name
                            delete_image_ID("tmpFPMresp");
                            list_image_ID();
                            // try to load the final child's partial FPMresp file
                            ID1 = load_fits(fname, "tmpFPMresp", 1);
                            list_image_ID();
                        }
                        // we found this child's partial FPMresp file!
                        // now insert it into our combined FPMresp file

                        fpt = fopen(fnamet,"a");
                        fprintf(fpt, "READING %s\n", fname);
                        fclose(fpt);

                        /*     list_image_ID();

                             printf("piaacmc[0].NBrings =                             %ld\n", piaacmc[0].NBrings);
                             printf("piaacmc[0].focmNBzone =                          %ld\n", piaacmc[0].focmNBzone);
                             printf("piaacmc[0].nblambda =                            %d\n", piaacmc[0].nblambda);
                             printf("data.image[ID].md[0].size[0] =                   %ld\n", data.image[ID].md[0].size[0]);
                             printf("data.image[piaacmc[0].zonezID].md[0].size[0] =   %ld\n", data.image[piaacmc[0].zonezID].md[0].size[0]);
                             fflush(stdout);
                             sleep(100);
                        */
                        
                        assert( ID1 != -1 );
                        if(ID1!=-1) // this should always be true
                        {
                            mzstep = PIAACMC_FPMresp_mp; // total number of tmux threads
                            mzoffset = thr; // the thread number of the child that just delivered its result
                            // insert the partial result of child thr into the combined FPMresp array
                            // be sure to skip the first line 'cause we already set it to be the light the went around the FPM
                            for(mz=1+mzoffset; mz<piaacmc[0].focmNBzone+1; mz+=mzstep) // loop over zone, do every PIAACMC_FPMresp_mp line
                            {

                                printf("mz = %ld    %ld %ld\n", mz, IDfpmresp, ID1);
                                fflush(stdout);
                                for(k=0; k<piaacmc[0].nblambda; k++) // for each wavelenth
                                    for(ii=0; ii<data.image[ID].md[0].size[0]; ii++) // for each evaluation point
                                    {
                                        // index of this evaluation point and wavelength and zone
                                        // tmpl1 = k*(nzones+1)*nEvaluationPoints) + zoneIndex*nEvaluationPoints + evaluationPoint
                                        tmpl1 = k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + mz*data.image[ID].md[0].size[0] + ii;
                                        // set the combined array value from the partial file (both are same shape and size, of course)
                                        data.image[IDfpmresp].array.D[tmpl1] = data.image[ID1].array.D[tmpl1];
                                        // subtract the current zone value from the first zone line, which contained all light
                                        // (with no mask).  Eventually this will contain only light that misses the FPM.
                                        data.image[IDfpmresp].array.D[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + ii] -= data.image[ID1].array.D[tmpl1];
                                    }
                            }
                            delete_image_ID("tmpFPMresp"); // we're dont with the partial array, so delete it
                        }
                   
                    }
                    // write out the current state of the combined FPMresp file
//                    sprintf(fname, "!%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ccnbr%03ld_ccz%06ld_ocr%04ld_ocz%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBringCentCone, (long) (1.0e9*piaacmc[0].fpmCentConeZ+0.1), (long) (100.0*piaacmc[0].fpmOuterConeRadld+0.1), (long) (1.0e9*piaacmc[0].fpmOuterConeZ+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
                   sprintf(fname, "!%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);


                    //                   sprintf(fname, "!%s/FPMresp%d_%02d_%d_%d_%02ld_%03ld_%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, PIAACMC_FPMsectors, (long) (10.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBrings, piaacmc[0].nblambda);
                    save_fits("FPMresp", fname);
                    // remove the child's .tmp file just in case (it should no longer exist 'cause we renamed, not copied, the .tmp file)
                    sprintf(command, "rm %s/FPMresp*.fits.tmp", piaacmcconfdir);
                    r = system(command);
                }
                else // we're a tmux child, so compute the response for our portion
                {
                    // name of the child partial FPMresp file (to become .tmp)
                 //   sprintf(fname,"!%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ccnbr%03ld_ccz%06ld_ocr%04ld_ocz%06ld_ssr%02d_ssm%d_%s_wb%02d_mp%02ld_thread%02ld.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBringCentCone, (long) (1.0e9*piaacmc[0].fpmCentConeZ+0.1), (long) (100.0*piaacmc[0].fpmOuterConeRadld+0.1), (long) (1.0e9*piaacmc[0].fpmOuterConeZ+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda, PIAACMC_FPMresp_mp, PIAACMC_FPMresp_thread);
                   sprintf(fname,"!%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d_mp%02ld_thread%02ld.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda, PIAACMC_FPMresp_mp, PIAACMC_FPMresp_thread);


                    // diagnostic file to make sure the child is working with the right zones
                    sprintf(fnametmp, "!%s/fpmzmap_thread%02ld.fits", piaacmcconfdir, PIAACMC_FPMresp_thread);
                    save_fits("fpmzmap", fnametmp);

                    printf("Making component %ld / %ld\n", PIAACMC_FPMresp_thread, PIAACMC_FPMresp_mp);
                    fflush(stdout);
                    WRITE_OK = 0;
                    // for each FPM zone, compute the response
                    // skip the first one 'cause it is not computed by the children
                    for(mz=1+mzoffset; mz<piaacmc[0].focmNBzone+1; mz+=mzstep)
                    {
                        focmMode = mz;  // focmMode can be a zone index, in which case operations are on that zone
                        optsyst[0].FOCMASKarray[0].mode = 0; // direct focal plane mask response

                        // default the reference to the 4th element
                        // but look for the element called "opaque mask at PIAA elem 1"
                        elem0 = 4;
                        for(elem=0; elem<optsyst[0].NBelem; elem++)
                        {
                            if(strcmp("opaque mask at PIAA elem 1", optsyst[0].name[elem])==0)
                            {
                                elem0 = elem;
                                printf("opaque mask at PIAA elem 1 = %ld\n", elem);
                            }
                             // raise an alarm if "opaque mask at PIAA elem 1" is not found *************************************
                        }

                        optsyst[0].keepMem[elem0] = 1; // keep it in memory

                        printf("piaacmc[0].NBrings =                             %ld\n", piaacmc[0].NBrings);
                        printf("piaacmc[0].focmNBzone =                          %ld\n", piaacmc[0].focmNBzone);
                        printf("piaacmc[0].nblambda =                            %d\n", piaacmc[0].nblambda);
                        printf("data.image[ID].md[0].size[0] =                   %ld\n", data.image[ID].md[0].size[0]);
                        printf("data.image[piaacmc[0].zonezID].md[0].size[0] =   %ld\n", data.image[piaacmc[0].zonezID].md[0].size[0]);
                        fflush(stdout);
                        //  sleep(100);

                        // compute the on-axis PSF
                        val = PIAACMCsimul_computePSF(0.0, 0.0, elem0, optsyst[0].NBelem, 0, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, 0);
                        // The PSF result for the evaluation points is put in array "imvect" which previously was
                        // assigned to another PSF result.
                        // Should put another ID = image_ID("imvect") here *************************************

                        // set the response of this zone from the PSF result
                        for(k=0; k<piaacmc[0].nblambda; k++) // loop over wavelength
                            for(ii=0; ii<data.image[ID].md[0].size[0]; ii++) // loop over evaluation points
                            {
                                // see previous example for explanation of indexing
                                // save response, which is just the value of the on-axis PSF at each evaluation point
                                data.image[IDfpmresp].array.D[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + mz*data.image[ID].md[0].size[0] + ii] = data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];
                                if(PIAACMC_FPMresp_mp==1) // if we're single threaded (no children)
                                    // subtract the current zone value from the first zone line, which contained all light
                                    // (with no mask).  Eventually this will contain only light that misses the FPM.
                                    data.image[IDfpmresp].array.D[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + ii] -= data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];
                            }


                        printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp"), fname1);
                        fflush(stdout);
                        // fname1 is the .tmp name
                        save_fits("FPMresp", fname1);
                        printf("Done \n");
                        fflush(stdout);
                    }
                    
                    //*************************** make up our minds about single threading, in the meantime say we don't support it

                    // partial file complete!  move it to the final file name so parent can see it
                    printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp"), fname2);
                    fflush(stdout);
                    // fname2 is the final name
                    save_fits_atomic("FPMresp", fname2);
                    printf("Done \n");
                    fflush(stdout);
                    WRITE_OK = 0;

                 //   sprintf(command, "mv %s %s", fname1, fname2);
                 //   r = system(command);
                }
            }
            else
                printf("File \"%s\" or \"%s\" exists\n", fname, fnamecomb);
        }
        focmMode = -1;
        break;




    case 13 : // optimize focal plane mask zones only
        // uses "fast" mode:
        // after mode 11, we can use the (complex) light propagated from each zone to compute the impact of
        // any thickness (sag) of that zone: the zone thickness induces a phase rotation for that zone,
        // which is applied to the unobstructed light from that zone as a complex rotation.
        //
        // The search is via steepest descent from random starting points
        //
        // this mode only sets up the optimization that actually happens after exiting the switch statement
        // if LINOPT = 1 (as does mode 40)
        printf("=================================== mode 013 ===================================\n");


        // set the name of the stopfile
        sprintf(stopfile, "%s/stoploop13.txt", piaacmcconfdir);



        // get cli variables
        // FPM sag regularization control flag, do if == 1
        REGFPMSAG = 1; // default
        if((IDv=variable_ID("REGFPMSAG"))!=-1)
            REGFPMSAG = (long) data.variable[IDv].value.f+0.01;

     
  

        // set current state for statistical tracking
        data.image[IDstatus].array.U[0] = 0;

        // usual initialization
        PIAAsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);

        PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);

        // set current state for statistical tracking
        data.image[IDstatus].array.U[0] = 1;

        // tracking diagnostic, giving the total flux in each plane
        sprintf(fname,"%s/flux.txt", piaacmcconfdir);
        fp = fopen(fname, "r");
        for(elem=0; elem<optsyst[0].NBelem; elem++)
        {
            ret = fscanf(fp, "%lf %lf  %d\n", &tmplf1, &tmplf2, &tmpd1);
            optsyst[0].flux[elem] = tmplf1/tmpd1*optsyst[0].nblambda;   // scale flux to current number of lambda
        }
        fclose(fp);


        LINOPT = 1; // perform linear optimization after the switch exits
        // get the number of iterations
        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            NBiter = (long) data.variable[IDv].value.f+0.01;
        else
            NBiter = 50; // default number of iterations

        // set current state for statistical tracking
        data.image[IDstatus].array.U[0] = 2;

        // get the FPMresp array computed in mode 11
//        sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ccnbr%03ld_ccz%06ld_ocr%04ld_ocz%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBringCentCone, (long) (1.0e9*piaacmc[0].fpmCentConeZ+0.1), (long) (100.0*piaacmc[0].fpmOuterConeRadld+0.1), (long) (1.0e9*piaacmc[0].fpmOuterConeZ+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
        sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);



        IDfpmresp = load_fits(fname, "FPMresp", 1);

        vsize = data.image[IDfpmresp].md[0].size[0]; // number of eval pts x2
        // make an array that holds the resulting light for evaluation point given the FPM solution, for each wavelenth
        ID = create_2Dimage_ID("imvect1", vsize, piaacmc[0].nblambda);

        // allocate arrays for fast routine
        // define convenient array variables
        fpmresp_array = data.image[IDfpmresp].array.D;
        zonez_array = data.image[piaacmc[0].zonezID].array.D;
        // allocate derivative of phase against thickness array
        dphadz_array = (double*) malloc(sizeof(double)*piaacmc[0].nblambda);
        // compute this derivative
        for(k=0; k<piaacmc[0].nblambda; k++)
        {
            // OPTICSMATERIALS_pha_lambda computes change in phase per unit thickness at specified wavelength
            // second arg is thickness, so 1.0 meters determines result per meter thickness
            dphadz_array[k] = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, 1.0, optsyst[0].lambdaarray[k]);
            printf("%ld  %g %g\n", k, optsyst[0].lambdaarray[k], dphadz_array[k]);
        }
        outtmp_array = (double*) malloc(sizeof(double)*(vsize*piaacmc[0].nblambda+data.image[piaacmc[0].zonezID].md[0].size[0]));

        // do the fast optimization using the results of mode 11
        computePSF_FAST_FPMresp = 1;

        // set current state for statistical tracking
        data.image[IDstatus].array.U[0] = 3;

        // read the contrast normalization factor into CnormFactor
        sprintf(fname, "%s/CnormFactor.txt", piaacmcconfdir);
        fp = fopen(fname, "r");
        ret = fscanf(fp, "%lf", &CnormFactor);
        fclose(fp);
        // for each zone, add a random offset in range +- MODampl
        // this randomizes the starting point for each zone
        // data.image[piaacmc[0].zonezID].array.D[k] is set in PIAACMCsimul_run()
        for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
            data.image[piaacmc[0].zonezID].array.D[k] += MODampl*(1.0-2.0*ran1());

        // set up optimization parameters for each zone
        // uses abstract specification of optimization parameters called paramval, ...
        NBparam = 0;
        for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++)
        {
            // parameter type
            paramtype[NBparam] = DOUBLE;
            // value: sag of each zone
            paramval[NBparam] = &data.image[piaacmc[0].zonezID].array.D[mz];
            // derivative step size
            paramdelta[NBparam] = 3.0e-9;
            // max parameter step size
            parammaxstep[NBparam] = 1.0e-6;
            // max and min allowed values for parameter
            parammin[NBparam] = piaacmc[0].fpmminsag;
            parammax[NBparam] = piaacmc[0].fpmmaxsag;
            // move on to next parameter
            NBparam++;
        }
        PIAACMC_FPM_FASTDERIVATIVES = 1; // for fast execution using analytic derivatives

        // set current state for statistical tracking
        data.image[IDstatus].array.U[0] = 4;
        // Now on to the actual optimization, after exit from the switch statement
        // I hope you have a lot of time...
        break;








    case 40 : // optimize PIAA optics shapes (and focal plane mask transmission for idealized PIAACMC)
        printf("=================================== mode 040 ===================================\n");
        //		FORCE_CREATE_fpmza = 1;
        PIAACMC_fpmtype = 0; // idealized (default)
        if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
            PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);
        printf("PIAACMC_fpmtype = %d\n", PIAACMC_fpmtype);


        FORCE_CREATE_fpmza = 1;
        PIAAsimul_initpiaacmcconf(PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);
        //       printf("data.image[piaacmc[0].zoneaID].array.D[0] = %lf\n", data.image[piaacmc[0].zoneaID].array.D[0]);
        //        sleep(10);



        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm


        if(0) // TEST
        {
            valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 1);
            printf("valref = %g\n", valref);
            printf("EXEC CASE 40 COMPLETED\n");
        }
        else
        {

            LINOPT = 1; // perform linear optimization
            if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
                NBiter = (long) data.variable[IDv].value.f+0.01;
            else
                NBiter = 1000;


            kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
            if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
                kmaxC = (long) data.variable[IDv].value.f+0.01;
            if(kmaxC>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
                kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];


            kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];
            if((IDv=variable_ID("PIAACMC_maxoptFterm"))!=-1)
                kmaxF = (long) data.variable[IDv].value.f+0.01;
            if(kmaxF>data.image[piaacmc[0].piaa0FmodesID].md[0].size[0])
                kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];




            // PIAA shapes regularization

            REGPIAASHAPES = 0; // default
            if((IDv=variable_ID("REGPIAASHAPES"))!=-1)
                REGPIAASHAPES = (long) data.variable[IDv].value.f+0.01;

            piaa0C_regcoeff = 0.0e-7; // regularization coeff
            piaa1C_regcoeff = 0.0e-7; // regularization coeff
            if((IDv=variable_ID("REGPIAA_C_COEFF"))!=-1)
            {
                piaa0C_regcoeff = data.variable[IDv].value.f;
                piaa1C_regcoeff = data.variable[IDv].value.f;
            }

            piaa0C_regcoeff_alpha = 1.0; // regularization coeff power
            piaa1C_regcoeff_alpha = 1.0; // regularization coeff power
            if((IDv=variable_ID("REGPIAA_C_ALPHA"))!=-1)
            {
                piaa0C_regcoeff_alpha = data.variable[IDv].value.f;
                piaa1C_regcoeff_alpha = data.variable[IDv].value.f;
            }

            piaa0F_regcoeff = 0.0e-7; // regularization coeff
            piaa1F_regcoeff = 0.0e-7; // regularization coeff
            if((IDv=variable_ID("REGPIAA_F_COEFF"))!=-1)
            {
                piaa0F_regcoeff = data.variable[IDv].value.f;
                piaa1F_regcoeff = data.variable[IDv].value.f;
            }

            piaa0F_regcoeff_alpha = 1.0; // regularization coeff power
            piaa1F_regcoeff_alpha = 1.0; // regularization coeff power
            if((IDv=variable_ID("REGPIAA_F_ALPHA"))!=-1)
            {
                piaa0F_regcoeff_alpha = data.variable[IDv].value.f;
                piaa1F_regcoeff_alpha = data.variable[IDv].value.f;
            }

            if(REGPIAASHAPES==1)
            {
				printf("loading CPA modes frequ\n");
				fflush(stdout);
                ID_CPAfreq = image_ID("cpamodesfreq");
                if(ID_CPAfreq == -1)
                    ID_CPAfreq = load_fits("cpamodesfreq.fits", "cpamodesfreq", 2);
            }
			

            // FPM SAG regularization

            REGFPMSAG = 1; // default
            if((IDv=variable_ID("REGFPMSAG"))!=-1)
                REGFPMSAG = (long) data.variable[IDv].value.f+0.01;




            NBparam = 0;

            if(PIAACMC_fpmtype==0) // ideal mask
            {
                paramtype[NBparam] = DOUBLE;
                paramval[NBparam] = &data.image[piaacmc[0].zoneaID].array.D[0];
                paramdelta[NBparam] = 1.0e-3;
                parammaxstep[NBparam] = 2.0e-1;
                parammin[NBparam] = -1.0;
                parammax[NBparam] = 1.0;
                NBparam++;
            }
            else // real physical mask
            {
                if(variable_ID("PIAACMC_mzOPT")!=-1) // optimize zones
                {
                    for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++)
                    {
                        paramtype[NBparam] = DOUBLE;
                        paramval[NBparam] = &data.image[piaacmc[0].zonezID].array.D[mz];
                        paramdelta[NBparam] = 1.0e-9;
                        parammaxstep[NBparam] = 5.0e-8;
                        parammin[NBparam] = -2.0e-6;
                        parammax[NBparam] = 2.0e-6;
                        NBparam++;
                    }
                }
            }


            for(k=0; k<kmaxC; k++)
            {
                paramtype[NBparam] = FLOAT;
                paramvalf[NBparam] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
                paramdelta[NBparam] = 1.0e-10;
                parammaxstep[NBparam] = 1.0e-7;
                parammin[NBparam] = -1.0e-3;
                parammax[NBparam] = 1.0e-3;
                NBparam++;
            }

            for(k=0; k<kmaxC; k++)
            {
                paramtype[NBparam] = FLOAT;
                paramvalf[NBparam] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
                paramdelta[NBparam] = 1.0e-10;
                parammaxstep[NBparam] = 1.0e-7;
                parammin[NBparam] = -1.0e-3;
                parammax[NBparam] = 1.0e-3;
                NBparam++;
            }

            for(k=0; k<kmaxF; k++)
            {
                paramtype[NBparam] = FLOAT;
                paramvalf[NBparam] = &data.image[piaacmc[0].piaa0FmodesID].array.F[k];
                paramdelta[NBparam] = 1.0e-10;
                parammaxstep[NBparam] = 1.0e-7;
                parammin[NBparam] = -1.0e-3;
                parammax[NBparam] = 1.0e-3;
                NBparam++;
            }

            for(k=0; k<kmaxF; k++)
            {
                paramtype[NBparam] = FLOAT;
                paramvalf[NBparam] = &data.image[piaacmc[0].piaa1FmodesID].array.F[k];
                paramdelta[NBparam] = 1.0e-10;
                parammaxstep[NBparam] = 1.0e-7;
                parammin[NBparam] = -1.0e-3;
                parammax[NBparam] = 1.0e-3;
                NBparam++;
            }

            FORCE_MAKE_PIAA0shape = 1;
            FORCE_MAKE_PIAA1shape = 1;
        }
        break;








    case 100 : // evaluate current design: polychromatic contrast, pointing sensitivity
        printf("=================================== mode 100 ===================================\n");
		
		

		// measure sensitivity to errors
	//	printf("Loading (optional) OPDerr file\n");
	//	fflush(stdout);
		// load an error if it exists
		IDopderrC = image_ID("OPDerrC");
		if(IDopderrC == -1)
			IDopderrC = load_fits("OPDerrC.fits", "OPDerrC", 0);
	

		if(IDopderrC != -1)
		{
			nbOPDerr = data.image[IDopderrC].md[0].size[2];  // number of error arrays
			//printf("INCLUDING %ld OPD ERROR MODES\n", nbOPDerr);
			//fflush(stdout);
		}
		else
		{
			//printf("NO OPD ERROR MODES\n");
			//fflush(stdout);
			nbOPDerr = 0;
		}

	
  
		printf("Will add optional OPD error modes (%ld modes)\n", nbOPDerr);
		fflush(stdout);
		


        PIAACMC_fpmtype = 0; // idealized (default)
        if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
            PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);


	

        FORCE_CREATE_fpmza = 1;
        PIAAsimul_initpiaacmcconf(PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);





        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        //        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0);
        //       printf("valref = %g\n", valref);

        //exit(0);




        ldoffset = 0.01; // default
        if((IDv=variable_ID("PIAACMC_ldoffset"))!=-1)
            ldoffset = data.variable[IDv].value.f;

        printf("ldoffset = %f\n", ldoffset);


        valref = PIAACMCsimul_computePSF(5.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_x50_y00.fits", piaacmcconfdir);
        save_fits("psfi0", fname);
        //load_fits(fname, "psfi");




        ID = image_ID("psfi0");
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];
        zsize = data.image[ID].md[0].size[2];
        peakarray = (float*) malloc(sizeof(float)*zsize);
        for(kk=0; kk<zsize; kk++)
        {
            peakarray[kk] = 0.0;
            for(ii=0; ii<xsize*ysize; ii++)
            {
                val = data.image[ID].array.F[kk*xsize*ysize+ii];
                if(val>peakarray[kk])
                    peakarray[kk] = val;
            }
        }
        avpeak = 0.0;
        for(kk=0; kk<zsize; kk++)
        {
            printf("peak %02ld  %10lf\n", kk, peakarray[kk]);
            avpeak += peakarray[kk];
        }
        avpeak /= zsize;
        free(peakarray);
        delete_image_ID("psfi0");

        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 1);
        sprintf(fname,"!%s/psfi0_x00_y00.fits", piaacmcconfdir);
        save_fits("psfi0", fname);
        //load_fits(fname, "psfi");


        ID = image_ID("psfi0");




        /// compute contrast curve
        focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;
        printf("focscale = %f\n", focscale);
        eval_sepNBpt = (long) (eval_sepmaxld/eval_sepstepld);
        eval_contrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
        eval_contrastCurve_cnt = (double*) malloc(sizeof(double)*eval_sepNBpt);
        for(ri=0; ri<eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] = 0.0;
            eval_contrastCurve_cnt[ri] = 0.0;
        }

        aveC = 0.0;
        aveCcnt = 0;

        for(kk=0; kk<zsize; kk++)
            for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
                {
                    xc = 1.0*ii-0.5*xsize;
                    yc = 1.0*jj-0.5*ysize;
                    xc *= focscale;
                    yc *= focscale;
                    rc = sqrt(xc*xc+yc*yc);
                    ri = (long) (rc/eval_sepstepld-0.5);
                    if(ri<0)
                        ri = 0;
                    if(ri<eval_sepNBpt)
                    {
                        eval_contrastCurve[ri] += data.image[ID].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        eval_contrastCurve_cnt[ri] += 1.0;
                    }
                    if((rc>2.0)&&(rc<6.0))
                    {
                        aveC += data.image[ID].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        aveCcnt++;
                    }
                }



        sprintf(fname, "%s/ContrastCurve_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d_tt000.txt", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

        fp = fopen(fname, "w");
        for(ri=0; ri<eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
            fprintf(fp, "%10f %10g %10g\n", eval_sepstepld*ri, eval_contrastCurve[ri], eval_contrastCurve_cnt[ri]);
        }
        fclose(fp);
        free(eval_contrastCurve);
        free(eval_contrastCurve_cnt);

       sprintf(fname, "%s/ContrastVal_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d_tt000.txt", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);


        fp = fopen(fname, "w");
        fprintf(fp, "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %02d %d\n", valref, aveC/aveCcnt, piaacmc[0].fpmaskradld, PIAACMC_MASKRADLD, piaacmc[0].centObs0, piaacmc[0].centObs1, (long) (piaacmc[0].lambda*1e9), (long) (piaacmc[0].lambdaB+0.1), PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda, (long) 0, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode);
        fclose(fp);




        // measure pointing sensitivity
        IDps = create_3Dimage_ID("starim", piaacmc[0].size, piaacmc[0].size, zsize);

        valref = 0.25*PIAACMCsimul_computePSF(ldoffset, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_p0.fits", piaacmcconfdir);
        save_fits("psfi0", fname);
        ID = image_ID("psfi0");
        for(ii=0; ii<xsize*ysize*zsize; ii++)
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];

        valref += 0.25*PIAACMCsimul_computePSF(-ldoffset, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_m0.fits", piaacmcconfdir);
        save_fits("psfi0", fname);
        ID = image_ID("psfi0");
        for(ii=0; ii<xsize*ysize*zsize; ii++)
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];

        valref += 0.25*PIAACMCsimul_computePSF(0.0, ldoffset, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_0p.fits", piaacmcconfdir);
        save_fits("psfi0", fname);
        ID = image_ID("psfi0");
        for(ii=0; ii<xsize*ysize*zsize; ii++)
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];

        valref += 0.25*PIAACMCsimul_computePSF(0.0, -ldoffset, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_0m.fits", piaacmcconfdir);
        save_fits("psfi0", fname);
        ID = image_ID("psfi0");
        for(ii=0; ii<xsize*ysize*zsize; ii++)
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];

        for(ii=0; ii<xsize*ysize*zsize; ii++)
            data.image[IDps].array.F[ii] /= 4.0;

        sprintf(fname, "!%s/psfi0_starim.fits", piaacmcconfdir);
        save_fits("starim", fname);

	
	
	 // measure sensitivity to errors

//	printf("Adding optional OPD error modes (%ld modes)\n", nbOPDerr);
//	fflush(stdout);
	// add error modes if any
	for(OPDmode=0; OPDmode < nbOPDerr; OPDmode++)
		{
			size = data.image[IDopderrC].md[0].size[0];
			IDopderr = create_2Dimage_ID("opderr", size, size);
            // "opderr" is a standard name read by PIAACMCsimul_init
			for(ii=0;ii<size*size;ii++)
				data.image[IDopderr].array.F[ii] = data.image[IDopderrC].array.F[size*size*OPDmode + ii];			
			PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0); // add error to the data
			PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
			sprintf(fname, "!%s/psfi0_opderr%02ld.fits", piaacmcconfdir, OPDmode);
			save_fits("psfi0", fname);
			delete_image_ID("opderr");
		}

        sprintf(fname, "!%s/psfi0_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

        save_fits("starim", fname);



        /// compute contrast curve
        /// measure average contrast value, 2-6 lambda/D
        focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;
        printf("focscale = %f\n", focscale);
        eval_sepNBpt = (long) (eval_sepmaxld/eval_sepstepld);
        eval_contrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
        eval_contrastCurve_cnt = (double*) malloc(sizeof(double)*eval_sepNBpt);
        for(ri=0; ri<eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] = 0.0;
            eval_contrastCurve_cnt[ri] = 0.0;
        }

        aveC = 0.0;
        aveCcnt = 0;

        for(kk=0; kk<zsize; kk++)
            for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
                {
                    xc = 1.0*ii-0.5*xsize;
                    yc = 1.0*jj-0.5*ysize;
                    xc *= focscale;
                    yc *= focscale;
                    rc = sqrt(xc*xc+yc*yc);
                    ri = (long) (rc/eval_sepstepld-0.5);
                    if(ri<0)
                        ri = 0;
                    if(ri<eval_sepNBpt)
                    {
                        eval_contrastCurve[ri] += data.image[IDps].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        eval_contrastCurve_cnt[ri] += 1.0;
                    }
                    if((rc>2.0)&&(rc<6.0))
                    {
                        aveC += data.image[IDps].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        aveCcnt++;
                    }
                }

        sprintf(fname, "%s/ContrastCurve_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);


        fp = fopen(fname, "w");
        for(ri=0; ri<eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
            fprintf(fp, "%10f %10g %10g\n", eval_sepstepld*ri, eval_contrastCurve[ri], eval_contrastCurve_cnt[ri]);
        }
        fclose(fp);

        free(eval_contrastCurve);
        free(eval_contrastCurve_cnt);

        sprintf(fname, "%s/ContrastVal_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

        fp = fopen(fname, "w");
        fprintf(fp, "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %7.3f %02d %d\n", valref, aveC/aveCcnt, piaacmc[0].fpmaskradld, PIAACMC_MASKRADLD, piaacmc[0].centObs0, piaacmc[0].centObs1, (long) (piaacmc[0].lambda*1e9), (long) (piaacmc[0].lambdaB+0.1), PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda, (long) (1000.0*ldoffset), piaacmc[0].fpmsagreg_coeff, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode);
        fclose(fp);

        break;



    case 101 : // transmission as a function of angular separation
        printf("=================================== mode 101 ===================================\n");

        printf("101: transm as a function of angular separation  ldoffset = %f\n", ldoffset);

        PIAACMC_fpmtype = 0; // idealized (default)
        if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
            PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);

        FORCE_CREATE_fpmza = 1;
        PIAAsimul_initpiaacmcconf(PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);

        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm


        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 0);
        sprintf(fname,"!%s/psfi0test_x00_y00.fits", piaacmcconfdir);
        save_fits("psfi0", fname);


       sprintf(fnametransm, "%s/transmCurve_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);


        fpt = fopen(fnametransm, "w");
        fclose(fpt);

        stepld = 0.001;
        for(xld=0.0; xld<10.0; xld+=stepld)
        {
            valref = PIAACMCsimul_computePSF(xld, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
            ID = image_ID("psfi0");
            //            sprintf(fname, "!psfi0transm_%04.1f.fits", xld);
            //           save_fits("psfi0", fname);
            printf("ID = %ld\n", ID);
            xsize = data.image[ID].md[0].size[0];
            ysize = data.image[ID].md[0].size[1];
            zsize = data.image[ID].md[0].size[2];
            printf("image size = %ld %ld %ld\n", xsize, ysize, zsize);
            val = 0.0;
            for(kk=0; kk<zsize; kk++)
            {
                for(ii=0; ii<xsize; ii++)
                    for(jj=0; jj<ysize; jj++)
                    {
                        dx = 1.0*ii-0.5*xsize;
                        dy = 1.0*jj-0.5*ysize;
                        if((dx*dx+dy*dy)<30.0*30.0)
                            val += data.image[ID].array.F[kk*xsize*ysize+jj*ysize+ii];
                    }
            }
            val /= zsize;

            fpt = fopen(fnametransm, "a");
            fprintf(fpt, "%10f %.18f\n", xld, val);
            fclose(fpt);
            delete_image_ID("psfi0");

            stepld = 0.001;
            stepld += 0.1*xld;
            if(stepld>0.2)
                stepld = 0.2;
        }
        break;



    case 300 : // import FPM configuration setting from parent directory
        printf("=================================== mode 300 ===================================\n");

        /*  sprintf(command, "cp conf_MASKRADLD.txt %s/", piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_FPMsectors.txt %s/", piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_NBrings.txt %s/", piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_nblambda.txt %s/", piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_resolved.txt %s/", piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_extmode.txt %s/", piaacmcconfdir);
          r = system(command);*/
        break;

    case 301 : // remove configuration settings
        printf("=================================== mode 301 ===================================\n");

        /*   sprintf(command, "mv %s/conf_MASKRADLD.txt %s/saveconf/conf_MASKRADLD.txt", piaacmcconfdir, piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_FPMsectors.txt %s/saveconf/conf_FPMsectors.txt", piaacmcconfdir, piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_NBrings.txt %s/saveconf/conf_NBrings.txt", piaacmcconfdir, piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_nblambda.txt %s/saveconf/conf_nblambda.txt", piaacmcconfdir, piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_resolved.txt %s/saveconf/conf_resolved.txt", piaacmcconfdir, piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_extmode.txt %s/saveconf/conf_extmode.txt", piaacmcconfdir, piaacmcconfdir);
           r = system(command);*/
        break;

    case 302 : // restore configuration settings
        printf("=================================== mode 302 ===================================\n");

        sprintf(command, "cp %s/saveconf/conf_*.txt %s/", piaacmcconfdir, piaacmcconfdir);
        r = system(command);
        break;






    default :
        printERROR(__FILE__,__func__,__LINE__, "mode not recognized");
        break;
    }



































    // linear optimization set up in modes 13 and 40
    //
    // the output parameters start as the evaluation zone values in "imvect"
    // If we are regularizing, we supplement the output parameters by adding
    // penalty terms as additional output parameters
    if(LINOPT == 1) // linear optimization
    {
        // for state tracking and statistics
        data.image[IDstatus].array.U[0] = 5;

        // Compute Reference on-axis performance contrast (valref)
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, 0);


        // save the current configuration to the _linopt directory
        sprintf(dirname, "%s_linopt", piaacmcconfdir);
        PIAAsimul_savepiaacmcconf(dirname);

        // import configuration from _linopt directory
        sprintf(command, "rsync -au --progress %s/* ./%s/", dirname, piaacmcconfdir);
        r = system(command);

        // "cp -n" will only copy the file if destination does not exist
        // this will ensure that the .ref.fits files are the PIAA shapes before any linear optimization
        // these will be the "reference" shapes used to regularize PIAA shapes: the deviation from this reference will be kept small by the regularization

        sprintf(command, "cp -n %s/piaa0Cmodes.fits %s/piaa0Cmodes.ref.fits", dirname, dirname);
        r = system(command);
        sprintf(command, "cp -n %s/piaa0Fmodes.fits %s/piaa0Fmodes.ref.fits", dirname, dirname);
        r = system(command);

        sprintf(command, "cp -n %s/piaa1Cmodes.fits %s/piaa1Cmodes.ref.fits", dirname, dirname);
        r = system(command);
        sprintf(command, "cp -n %s/piaa1Fmodes.fits %s/piaa1Fmodes.ref.fits", dirname, dirname);
        r = system(command);

        sprintf(command, "cp -n %s/piaacmcparams.conf %s/piaacmcparams.ref.conf", dirname, dirname);
        r = system(command);


        // load PIAA reference shapes (used for regularization)
        sprintf(fname, "%s/piaa0Cmodes.ref.fits", dirname);
        load_fits(fname, "piaa0Cmref", 1);

        sprintf(fname, "%s/piaa1Cmodes.ref.fits", dirname);
        load_fits(fname, "piaa1Cmref", 1);

        sprintf(fname, "%s/piaa0Fmodes.ref.fits", dirname);
        load_fits(fname, "piaa0Fmref", 1);

        sprintf(fname, "%s/piaa1Fmodes.ref.fits", dirname);
        load_fits(fname, "piaa1Fmref", 1);

        // we have now saved the starting point of the optimization for future comparison
        // in the <piaacmcconfdir>_linopt directory

        // here we compute regularization value of piaashapes and store it in the val0 variable
        // if regularization is of PIAA shapes is ON, then val0 will be computed and added to the overal performance metric valref
        val0 = 1.0;

        // regularize the piaashapes via a penalty added to the reference contrast valref
        // The optimization minimizes the summed contrast + val0 + val1.
        // Regularization is via adding a constant val0 + val1 to the contrast we're minimizing
        // note that here we're setting output parameters.
        if(REGPIAASHAPES==1)
        {
            // first we compute the starting regularization constant
            val0 = 0.0;
            // index of the PIAA element 0 (first mirror) shapes via cosine modes
            ID = piaacmc[0].piaa0CmodesID;
            // index of PIAA shapes reference image
            IDref = image_ID("piaa0Cmref");
            if(IDref==-1)
            {   // error message if we get here?  ***************************************
                // if the reference image doesn't exist, create it
                IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                // initialize to zero shape
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            // This section of code does not actually fill in the regularization terms in the output vector
            // filling in is done later.  Here we are only computing the initial reference scalar objective value
            
            // for each cosine mode set the optimization parameter = cosine mode modified by regularization
            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
            {
                // compute square of C*(deviation from reference)*(mode index)^(alpha)
                // so higher-index zones (higher spatial frequency) are more
                // heavily penalized
                tmp = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa0C_regcoeff_alpha);
                val0 += tmp*tmp;
            }


            // do the same for PIAA element 1
            ID = piaacmc[0].piaa1CmodesID;
            IDref = image_ID("piaa1Cmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
            {
                tmp = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa1C_regcoeff_alpha);
                val0 += tmp*tmp;
            }

            // get spatial frequency of each mode in cycles/aperture
            ID_CPAfreq = image_ID("cpamodesfreq");

            // do the same for PIAA element 0 and 1 for the Fourier modes
            // this time use the actual spatial frequency rather than mode index as proxy for frequency
            ID = piaacmc[0].piaa0FmodesID;
            IDref = image_ID("piaa0Fmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
            {
                tmp = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                val0 += tmp*tmp;
            }

            ID = piaacmc[0].piaa1FmodesID;
            IDref = image_ID("piaa1Fmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
            {
                tmp = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                val0 += tmp*tmp;
            }

            printf("VALREF = %g + %g -> %g\n", valref, val0, valref+val0);
            valref += val0;
        }


        // val1 is the regularization value for the focal plane mask sag values
        // same as above, but not dependent on position
        val1 = 1.0;
        if(REGFPMSAG == 1)
        {
            ID = piaacmc[0].zonezID;
            val1 = 0.0;
            for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
            {
                // compute the square of (sag/coeff)^alpha
                tmp = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                val1 += tmp*tmp;
            }
            valref += val1;
        }
        // At this point all we've done is compute the overall performance metric including
        // regularization in valref.


        // for state tracking and statistics
        data.image[IDstatus].array.U[0] = 6;
        printf("================================ Reference = %g\n", valref);

        // copy imvect to vecDHref "vector dark hole reference"
        // vecDHref is the dark hole complex amplitude state at the beginning of the linear optimization
        // (this dark hole is nominally a full annulus from 1.5 to ~8 lambda/D, created at the top
        // of PIAACMCsimul_computePSF with size controlled by scoringIWA and scoringOWA
        // the corresponding performance metric is valref
        chname_image_ID("imvect", "vecDHref"); // note: imvect was computed by PIAACMCsimul_computePSF called ~150 lines above
        ID = image_ID("vecDHref"); // ID changed identity
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];



        // for state tracking and statistics
        data.image[IDstatus].array.U[0] = 7;

        // now we will just determine the size of the size of the
        // optimization vectors that we will actually fill in later
        // save vecDHref initial state as a reference
        sprintf(fname, "!%s/vecDMref.fits", piaacmcconfdir);
        save_fits("vecDHref", fname);
        // get and copy size of vecDHref, 'cause we're manipulating size1Dvec
        size1Dvec = data.image[ID].md[0].nelement;
        size1Dvec0 = size1Dvec;

        // PIAA shapes regularization
        // if regularization is turned on, the size of the evaluation vector is increased to include PIAA shape coefficients, in addition to complex amplitude in focal plane
        // the optimization code will then simultaneously minimize the light in the focal plane AND the PIAA shape deviations from the nominal shape
        if(REGPIAASHAPES==1)
        {
            // there are 4 groups of PIAA shape parameters: 2 sets of cosine modes and 2 sets of Fourier modes
            size1Dvec += data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa1CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa1FmodesID].md[0].size[0];
        }

        // The same approach is used for regularization of the focal plane mask sag values
        // the sag values are appended to the evaluation vector
        if(REGFPMSAG==1)
        {
            size1Dvec += data.image[piaacmc[0].zonezID].md[0].size[0];
        }


        // re-package vector into 1D array and add regularization terms
        // the resulting vector is stored in vecDHref1D
        // we also create a mask image which is intended to mask out some of the pixels from the evaluation
        // the mask is currently not used, so we will write 1.0 in all of its pixels
        // DHmask and vecDHref1D contain all the optimization parameters as set above
        IDm = create_2Dimage_ID("DHmask", size1Dvec, 1); // "ID of evaluation mode mask"
        ID1Dref = create_2Dimage_ID("vecDHref1D", size1Dvec, 1); // "ID of 1D dark zone reference"

        // we first write 1.0 into the focal plane complex amplitudes in the vector
        ID = image_ID("vecDHref");
        for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
        {
            // imbed vecDHref into the evaluation zone part of of the full parameter vector vecDHref1D
            data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
            // sets the evaluation zone part of of the full parameter vector vecDHref1D to 1
            // 1 means the evaluation zone is on
            data.image[IDm].array.F[ii] = 1.0;
        }
        // !!!!!! WARNING !!!!!!!
        // the state of ii at this point drives the code below and will evolve until the comment
        // that says we're done with ii

        // for state tracking and statistics
        data.image[IDstatus].array.U[0] = 8;

        // Now actually fill in the regularized output vector.
        // If we are not regularizing, the output evaluation zone values are filled in by the
        // PSF calls in the optimization loop
        // and then append regularization vectors in the main evaluation vector
        // Note: index ii is incremented as we add "sub-vectors" into the main evaluation vector
        if(REGPIAASHAPES == 1)
        {
            // initialize by filling in the regularization terms of the output,
            // !! starting at the current value of ii !!
            // This means we're actually writing the output vector regularization terms.
            // Otherwise this is the same as the if(REGPIAASHAPES == 1) code block above

            ID = piaacmc[0].piaa0CmodesID;
            IDref = image_ID("piaa0Cmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa0C_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }

            ID = piaacmc[0].piaa1CmodesID;
            IDref = image_ID("piaa1Cmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa1C_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }

            ID_CPAfreq = image_ID("cpamodesfreq");
            ID = piaacmc[0].piaa0FmodesID;
            IDref = image_ID("piaa0Fmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }

            ID = piaacmc[0].piaa1FmodesID;
            IDref = image_ID("piaa1Fmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }
        }




        // same for the sags, starting at the current value of ii
        if(REGFPMSAG == 1)
        {
            ID = piaacmc[0].zonezID;
            for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = pow( data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }
        }
        // !!!!!!
        // we're done with ii


        delete_image_ID("vecDHref"); // vecDHref has beem embedded into vecDHref1D

        // at this point, we have completed the initialization, and the optimization loop starts


        initbestval = 0;
        // file that will track optimization loop progress
        sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
        fp = fopen(fname, "w");
        fclose(fp);

        //	list_image_ID();
        //
        // LINEAR OPTIMIZATION AROUND CURRENT POINT
        //

        iterOK=1;
        iter = 0;
        oldval = 1.0;
        data.image[IDstatus].array.U[0] = 9;

        // while # of iterations < NBiter
        //  and the ojective changes by more than 2% after the second iteration
        //  and something about NBlinoptgain ???????
        while(iterOK==1)//        for(iter=0; iter<NBiter; iter++)
        {
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 10;
            printf("Iteration %ld/%ld\n", iter, NBiter);
            fflush(stdout);
            // array for collecting dark hole mode derivatives
            // stores derivative of output vector against input parameters
            IDmodes = create_3Dimage_ID("DHmodes", size1Dvec, 1, NBparam);
            // 2D array for diagnostic display
            IDmodes2D = create_2Dimage_ID("DHmodes2D", size1Dvec, NBparam); //TEST

            // get ready to update optimization tracking file
            sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
            fp = fopen(fname, "a");
            fprintf(fp, "### PIAACMC_FPM_FASTDERIVATIVES = %d\n", PIAACMC_FPM_FASTDERIVATIVES);
            fclose(fp);
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 11;

            // compute local derivatives of output vector against input focal plane mask zones (sags)
            if(PIAACMC_FPM_FASTDERIVATIVES == 1) // TO BE USED ONLY FOR FOCAL PLANE MASK OPTIMIZATION
            { // this only happens in mode 13
                // the fast derivative mode only works for focal plane mask optimization, for which derivatives against sag values can be comptuted by simple rotation of pre-computed vectors from mode 11
                // for state tracking and statistics
                data.image[IDstatus].array.U[0] = 12;
                optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                //				ID = create_2Dimage_ID("DHmodes2Dtest", size1Dvec, NBparam);
                for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++) // loop over mask zones
                {
                    // actually compute the derivative
                    // fpmresp_array is results from mode 11
                    // from mode 13 above:
                    //      fpmresp_array = data.image[IDfpmresp].array.D;
                    //      zonez_array = data.image[piaacmc[0].zonezID].array.D;
                    // dphadz_array was computed in mode 13 shortly afterwards
                    // outtmp_array is output
                    PIAACMCsimul_achromFPMsol_eval_zonezderivative(mz, fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
                    for(ii=0; ii<size1Dvec0; ii++)
                        data.image[IDmodes].array.F[mz*size1Dvec+ii] = outtmp_array[ii]*paramdelta[mz];
                }
                // derivatives of regularization values against sag values can also be computed analytically without requiring diffraction propagation
                if(REGFPMSAG == 1)
                {
                    ID = piaacmc[0].zonezID;
                    // following should be derivative of (sag/coeff)^alpha
                    // w.r.t. sag
                    for(mz=0; mz < data.image[ID].md[0].size[0]; mz++)
                        data.image[IDmodes].array.F[mz*size1Dvec + (size1Dvec0+mz)] = (piaacmc[0].fpmsagreg_alpha/piaacmc[0].fpmsagreg_coeff) * pow( data.image[ID].array.D[mz]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha-1.0)*paramdelta[mz];
                }


                // TEST diagnostic
                memcpy(data.image[IDmodes2D].array.F, data.image[IDmodes].array.F, sizeof(float)*size1Dvec*NBparam);
                save_fl_fits("DHmodes2D", "!test_DHmodes2D.fits");

                // for state tracking and statistics
                data.image[IDstatus].array.U[0] = 13;
            }
            else // ONLY FOR PIAA SHAPES OPTIMIZATION
            {
                // derivatives against PIAA shapes must be computed numerically
                
                // for state tracking and statistics
                data.image[IDstatus].array.U[0] = 14;
                for(i=0; i<NBparam; i++)
                {
                    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                    // get delta on-axis PSF response to the change in paramdelta
                    // to later give derivative w.r.t. paramdelta
                    if(paramtype[i]==FLOAT)
                    {
                        *(paramvalf[i]) += (float) paramdelta[i];
                        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, 0);
                    }
                    else
                    {
                        *(paramval[i]) += paramdelta[i];
                        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, 0);
                    }

                    //      sprintf(fname,"!%s/imvect_%02ld.fits", piaacmcconfdir, i);
                    //       save_fits("imvect", fname);
                    ID = image_ID("imvect");


                    sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
                    fp = fopen(fname, "a");
                    fprintf(fp, "# %5ld/%5ld %5ld/%5ld %20.15g %20.15g %20.15g      %20.15g\n", iter, NBiter, i, NBparam, paramdelta[i], val, valref, bestval);
                    fclose(fp);

                    // re-package vector into 1D array and add regularization terms
					// evaluation vector is "imvect1D", ID = ID1D
                    // similar to vecDHref before
                    ID1D = create_2Dimage_ID("imvect1D", size1Dvec, 1);
                    // fill in the evaluation point portion
                    for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
                        data.image[ID1D].array.F[ii] = data.image[ID].array.F[ii];

                    if(REGPIAASHAPES==1)
                    { // fill in the shape regularization value's response to paramdelta using the
                        // same formulas as before with the delta PSF as input
                        ID = piaacmc[0].piaa0CmodesID;
                        IDref = image_ID("piaa0Cmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        {
                            data.image[ID1D].array.F[ii] = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj,piaa0C_regcoeff_alpha);
                            ii++;
                        }

                        ID = piaacmc[0].piaa1CmodesID;
                        IDref = image_ID("piaa1Cmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                        {
                            data.image[ID1D].array.F[ii] = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj,piaa1C_regcoeff_alpha);
                            ii++;
                        }

                        ID_CPAfreq = image_ID("cpamodesfreq");

                        ID = piaacmc[0].piaa0FmodesID;
                        IDref = image_ID("piaa0Fmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                        {
                            data.image[ID1D].array.F[ii] = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                            ii++;
                        }

                        ID = piaacmc[0].piaa1FmodesID;
                        IDref = image_ID("piaa1Fmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                        {
                            data.image[ID1D].array.F[ii] = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                            ii++;
                        }
                    }


                    delete_image_ID("imvect"); // has been imbedded into imvect1D

					// restore original state (return to original staring point)
                    if(paramtype[i]==FLOAT)
                        *(paramvalf[i]) -= (float) paramdelta[i];
                    else
                        *(paramval[i]) -= paramdelta[i];


                    // compute actual derivative as first difference from reference
                    // this is the starting derivative
                    for(ii=0; ii<data.image[ID1D].md[0].nelement; ii++)
                        data.image[IDmodes].array.F[i*data.image[ID1D].md[0].nelement+ii] = (data.image[ID1D].array.F[ii] - data.image[ID1Dref].array.F[ii]);


                    //    printf("%3ld %g %g\n", i, val, valref);

                    // create diagnostic image
                    ID = create_2Dimage_ID("DHmodes2D", size1Dvec, NBparam);
                    for(ii=0; ii<data.image[IDmodes].md[0].nelement; ii++)
                        data.image[ID].array.F[ii] = data.image[IDmodes].array.F[ii];

                    sprintf(fname, "!%s/DMmodes.fits", piaacmcconfdir);
                    save_fits("DHmodes2D", fname);

                    delete_image_ID("DHmodes2D");
                }
                // for state tracking and statistics
                data.image[IDstatus].array.U[0] = 15;
            }


            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 16;
            // print the results to file for human tracking
            sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
            fp = fopen(fname, "a");
            fprintf(fp, "### scanning gain \n");
            fprintf(fp, "### <alphareg>  <gain>  <contrast>\n");
            fclose(fp);
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 17;

            // first three arguments are names of the input arrays
            // vecDHref1D is the input data
            // DHmodes is the basis of modes to expand vecDHref1D into
            // DHmask weights elements of vecDHref1D (nominally all weights = 1)
            // 4th arg is the pseudoinverse eigenvalue (via eigenvalue decomposition)
            // this decomposes vecDHref1D into the DHmodes
            // 5th arg is the output: optcoeff0*DHmodes = vecDHref1D
            // computed via pseudoinverse of DHmodes
            // This decomposition facilitates the cancellation of vecDHref1D by
            // searching in the DHmodes basis
            //
            // use three cutoff values to give three options for future evaluation
            // smallest cutoff values produce the largest changes (are least well conditioned)
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.1, "optcoeff0", 0);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.01, "optcoeff1", 0);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.001, "optcoeff2", 0);
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 18;

            // initialize zero "optimal" vector optvec giving direction to move in the search for the min
            arith_image_cstmult("optcoeff0", 0.0, "optvec"); // create optimal vector
            IDoptvec = image_ID("optvec");
            // initialize the objective value
            initbestval = 0;
            bestval = valref;

            alphareg = 1.0; // has no effect (see next loop)

            // say something here ???
            NBlinoptgain = 0;


            scangainfact = 1.2;
            alphascaninit = 0;
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 19;
            // alphareg controls linear combinations of the directions optcoeff0,1,2
            // alphareg = 0 => moving along optcoeff0
            // alphareg = 0.5 => moving along optcoeff1
            // alphareg = 1 => moving along optcoeff2
            // with interpolation
            // look in 5 steps if alphareg += 0.2
            for(alphareg=0.0; alphareg<1.01; alphareg += 0.2)
            {
                // produce piecewise linear interoplation coefficients
                acoeff0 = 1.0 - 2.0*alphareg; // ranges from -1 to 1
                if(acoeff0<0.0)
                    acoeff0 = 0.0; // clip to 0 to 1, = 0 if alphareg > 0.5

                acoeff1 = 1.0 - fabs(2.0*(alphareg-0.5)); // two lines: from 0,1 to .5,0 to 1,1

                acoeff2 = 2.0*alphareg - 1.0; // ranges from -1 to 1
                if(acoeff2<0.0)
                    acoeff2 = 0.0; // clip to 0 to 1, = 0 if alphareg < 0.5
                
                // sum of acoeff0,1,2 = 1 at all values of alphareg.

                // for state tracking and statistics
                data.image[IDstatus].array.U[0] = 20;
                
                // optcoeff0m = acoeff0*optcoeff0 etc.
                arith_image_cstmult("optcoeff0", acoeff0, "optcoeff0m");
                arith_image_cstmult("optcoeff1", acoeff1, "optcoeff1m");
                arith_image_cstmult("optcoeff2", acoeff2, "optcoeff2m");

                // diagnostic
                save_fl_fits("optcoeff0", "!optcoeff0.fits");//TEST
                save_fl_fits("optcoeff1", "!optcoeff1.fits");
                save_fl_fits("optcoeff2", "!optcoeff2.fits");

                // optcoeff01m = acoeff0*optcoeff0 + acoeff1*optcoeff1
                arith_image_add("optcoeff0m", "optcoeff1m", "optcoeff01m");
                // optcoeff = acoeff0*optcoeff0 + acoeff1*optcoeff1 + acoeff2*optcoeff2
                arith_image_add("optcoeff01m", "optcoeff2m", "optcoeff");
                // optcoeff now has our search direction
                delete_image_ID("optcoeff0m");
                delete_image_ID("optcoeff1m");
                delete_image_ID("optcoeff2m");
                delete_image_ID("optcoeff01m");

                ID = image_ID("optcoeff");
                // for state tracking and statistics
                data.image[IDstatus].array.U[0] = 21;

                // do linear scan along the direction optcoeff from current parameter location
                linscanOK = 1;
                // size of step in this direction
                scangain = 0.0; //scanstepgain; overriden below
                val = 100000000000.0; // initialize minimization objective to a big number
                bestgain = 0.0;
                // iteration counter of steps in the current direction optcoeff
                k = 0;

                // if(alphascaninit==1)
                //		scangain = bestgain/scangainfact/scangainfact/scangainfact/scangainfact/scangainfact;
                //	alphascaninit = 1;

                // scangain is our location along the direction optcoeff
                scangain = 0.001;
                //              scanstepgain = 0.000001; // TEST
                //              scangainfact = 1.00001; // TEST

                // while objective value < previous value and we've taken no more than than 90 steps
                while(linscanOK==1)
                {
                    // for state tracking and statistics
                    data.image[IDstatus].array.U[0] = 22;

                    // compute offsets
                    ID = image_ID("optcoeff"); // direction vector
                    linoptlimflagarray[k] = 0;
                    // step each parameter by optcoeff
                    for(i=0; i<NBparam; i++) // looping over parameters
                    {
                        // compute step delta for this parameter
                        // image[ID] = optcoeff is a derivative w.r.t. paramdelta, so has dimension
                        // (parameter dimension)/(paramdelta dimension) so we have to mulitply by
                        // paramdelta to put our step in physical parameter units
                        // negative because we want to cancel the value from the delta PSF
                        paramdeltaval[i] = -scangain * data.image[ID].array.F[i] * paramdelta[i];
                        if(paramdeltaval[i]<-parammaxstep[i]) // if the step is too large in the negative direction
                        {
                            printf("MIN LIMIT [%3ld   %20g]   %20g -> ", i, paramdelta[i], paramdeltaval[i]); //TEST
                            paramdeltaval[i] = -parammaxstep[i]; // set it to the negative largest allowed step
                            printf(" %20g\n", paramdeltaval[i]); //TEST
                            linoptlimflagarray[k] = 1;
                        }
                        if(paramdeltaval[i]>parammaxstep[i])// if the step is too large in the positive direction
                        {
                            printf("MAX LIMIT [%3ld   %20g]   %20g -> ", i, paramdelta[i], paramdeltaval[i]); //TEST
                            paramdeltaval[i] = parammaxstep[i]; // set it to the positive largest allowed step
                            printf(" %20g\n", paramdeltaval[i]); //TEST
                            linoptlimflagarray[k] = 1;
                        }

                        // apply offsets to the global data object via the pointers paramvalf, which
                        // point into the (hopefully) coorect locations of each parameter's value in the
                        // data object
                        if(paramtype[i]==FLOAT)
                        {
                            if(  *(paramvalf[i]) + (float) paramdeltaval[i]  > parammax[i] )
                                // if we're about to step too far, set the step to the
                                // limit - current parameter value, so we step to the limit
                                paramdeltaval[i] = parammax[i] - *(paramvalf[i]);

                            if(  *(paramvalf[i]) + (float) paramdeltaval[i]  < parammin[i] )
                                // if we're about to step too far, set the step to the
                                // limit - current parameter value, so we step to the limit
                                paramdeltaval[i] = parammin[i] - *(paramvalf[i]);
                            // take the actual step (paramvalf is a 1D array of pointers)
                            *(paramvalf[i]) += (float) paramdeltaval[i];
                        }
                        else // same for the double case
                        {
                            if(  *(paramval[i]) + paramdeltaval[i]  > parammax[i] )
                                paramdeltaval[i] = parammax[i] - *(paramval[i]);

                            if(  *(paramval[i]) + paramdeltaval[i]  < parammin[i] )
                                paramdeltaval[i] = parammin[i] - *(paramval[i]);

                            *(paramval[i]) += paramdeltaval[i];
                        }
                    }
                    // store the current objective value for later comparison
                    valold = val;
                    // for state tracking and statistics
                    data.image[IDstatus].array.U[0] = 23;

					// compute new state and compute assossiated evaluation metric
                    // using the modified global data object
                    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, 0);
                    valContrast = val; // contrast component of the evaluation metric
                    // we've now only done the light portion
                    
                    // add regularization component of the evaluation metrix
                    // first, compute and add PIAA shape regularization value (val0) if applicable
                    // this is the same as the previous computation of val0 around line 7627
                    val0 = 1.0;
                    if(REGPIAASHAPES==1)
                    {
                        val0 = 0.0;
                        ID = piaacmc[0].piaa0CmodesID;
                        IDref = image_ID("piaa0Cmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        {
                            tmp = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa0C_regcoeff_alpha);
                            val0 += tmp*tmp;
                        }

                        ID = piaacmc[0].piaa1CmodesID;
                        IDref = image_ID("piaa1Cmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                        {
                            tmp = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa1C_regcoeff_alpha);
                            val0 += tmp*tmp;
                        }


                        ID_CPAfreq = image_ID("cpamodesfreq");

                        ID = piaacmc[0].piaa0FmodesID;
                        IDref = image_ID("piaa0Fmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                        {
                            tmp = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                            val0 += tmp*tmp;
                        }

                        ID = piaacmc[0].piaa1FmodesID;
                        IDref = image_ID("piaa1Fmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                        {
                            tmp = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                            val0 += tmp*tmp;
                        }
                        val += val0;
                    }

					// add sag regularization (val1) if applicable, as before
                    val1 = 1.0;
                    if(REGFPMSAG == 1)
                    {
                        ID = piaacmc[0].zonezID;
                        val1 = 0.0;
                        for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
                        {
                            tmp = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                            val1 += tmp*tmp;
                        }
                        val += val1;
                    }
                    // val is now our complete objective!! Yay!!

                    
                    // for state tracking and statistics
                    data.image[IDstatus].array.U[0] = 24;
                    // print it for monitoring
                    sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
                    fp = fopen(fname, "a");
                    // printf the first part of the line reporting current values
                    fprintf(fp, "##  %5.3f   %20lf           %20g    (reg = %12g %12g   contrast = %20g)       [%d] [%ld]", alphareg, scangain, val, val0, val1, valContrast, linoptlimflagarray[k], NBparam);
 
                    // now add text indicating status and complete line
                    // and store all parameters for the current best solution
                    if((val<bestval)||(initbestval==0))
                    {
                        for(i=0; i<NBparam; i++)
                            if(paramtype[i]==FLOAT)
                                data.image[IDoptvec].array.F[i] = *(paramvalf[i]);
                            else
                                data.image[IDoptvec].array.F[i] = (float) *(paramval[i]);
                        bestval = val;
                        if(initbestval == 0)
                            fprintf(fp, " ===== START POINT =====\n");
                        else
                            fprintf(fp, "  -> BEST VECTOR =======\n");
                        bestgain = scangain;
                        initbestval = 1;
                    }
                    else
                    {
                        fprintf(fp, " bestval = %12g\n", bestval);
                    }
                    fclose(fp);

                    // remove offsets returning the global data object to its original state
                    for(i=0; i<NBparam; i++)
                    {
                        if(paramtype[i]==FLOAT)
                            *(paramvalf[i]) -= (float) paramdeltaval[i];
                        else
                            *(paramval[i]) -= paramdeltaval[i];
                    }
                    // for state tracking and statistics
                    data.image[IDstatus].array.U[0] = 25;

                    // store the current position and value
                    linoptgainarray[k] = scangain;
                    linoptvalarray[k] = val;
                    k++; // next step

                    // test to see if we're no longer getting better
                    if(val<valold)
                    {
                        linscanOK = 1; // if we're getting better keep going
                        // bestgain = scangain;
                        // scangain += scanstepgain;
                    }
                    else // otherwise stop stepping
                        linscanOK = 0;



                    if(k>90)  // stop if we've taken too many steps
                        linscanOK = 0;

                    // increment our location along the line
                    scangain += scanstepgain; // scanstepgain is an initilizaed function local (currently 0.001)
                    scangain *= scangainfact; // causes later steps to be larger
                                                // (implicit scangainfact^n for the nth step)
                }
                // NBlinoptgain is counting the largest number of steps needed in this inner loop
                // stepping in the current direction.  When this is small (< 3) we declare victory
                // and stop the outer linear optimization
                if(k>NBlinoptgain)
                    NBlinoptgain = k;

                delete_image_ID("optcoeff"); // delete the current direction
                // for state tracking and statistics
                data.image[IDstatus].array.U[0] = 26;
            }
            // best solution after this linear linescan is stored in IDoptvec
            delete_image_ID("optcoeff0");
            delete_image_ID("optcoeff1");
            delete_image_ID("DHmodes");

            // we've now found the minimum using the three directions from the
            // alternative decompositions of the parameter space with the DHmodes basis
            // (with different conditioning).
            // Now we check the result by recomputing from scratch the objective with the current
            // (hopefully) optimal parameters.  Belts and suspenders
            
            // At this point the global data object has been restored to its original
            // (non-optimal) state.
            
            // (re-)compute best solution identified in previous linescan
			// update state to best solution (IDoptvec), setting the global data object
            // to the optimal state
            for(i=0; i<NBparam; i++)
            {
                if(paramtype[i]==FLOAT)
                    *(paramvalf[i]) = data.image[IDoptvec].array.F[i];
                else
                    *(paramval[i]) = (double) data.image[IDoptvec].array.F[i];
            }
            valold = val;
			// compute contrast metric component -> val using the data object in the latest optimal state
            val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, 0);

			// add PIAA shape regularization component (val0) if applicable
            // same val0 and val1 computations are before, using the latest optimal state
            val0 = 1.0;
            if(REGPIAASHAPES==1)
            {
                val0 = 0.0;
                ID = piaacmc[0].piaa0CmodesID;
                IDref = image_ID("piaa0Cmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                {
                    tmp = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa0C_regcoeff_alpha);
                    val0 += tmp*tmp;
                }

                ID = piaacmc[0].piaa1CmodesID;
                IDref = image_ID("piaa1Cmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                {
                    tmp = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa1C_regcoeff_alpha);
                    val0 += tmp*tmp;
                }


                ID_CPAfreq = image_ID("cpamodesfreq");

                ID = piaacmc[0].piaa0FmodesID;
                IDref = image_ID("piaa0Fmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                {
                    tmp = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                    val0 += tmp*tmp;
                }

                ID = piaacmc[0].piaa1FmodesID;
                IDref = image_ID("piaa1Fmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                {
                    tmp = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                    val0 += tmp*tmp;
                }
                val += val0;
            }
            
            // add sag regularization component if applicable
            val1 = 1.0;
            if(REGFPMSAG == 1)
            {
                ID = piaacmc[0].zonezID;
                val1 = 0.0;
                for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
                {
                    tmp = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                    val1 += tmp*tmp;
                }
                val += val1;
            }
            
            // now val is the objective including any desired regularization terms using the
            // latest optimal solution


            printf("gain: %lf -> val = %20g\n", bestgain, val);
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 27;



			// update reference state evaluation vector of optimal results including evaluation zone values
            // and desired regularization terms
            // which sets up starting the next iteration at the best solution
            
            ID1Dref = image_ID("vecDHref1D");
            ID = image_ID("imvect");
            // first fill in evaluation zone (complex) values
            for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
                data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 28;

            // now fill in the regularization terms if desired
            // same code as before.
            if(REGPIAASHAPES==1)
            {
                ID = piaacmc[0].piaa0CmodesID;
                IDref = image_ID("piaa0Cmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }

                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj,piaa0C_regcoeff_alpha);
                    ii++;
                }

                ID = piaacmc[0].piaa1CmodesID;
                IDref = image_ID("piaa1Cmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj,piaa1C_regcoeff_alpha);
                    ii++;
                }


                ID_CPAfreq = image_ID("cpamodesfreq");

                ID = piaacmc[0].piaa0FmodesID;
                IDref = image_ID("piaa0Fmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) *pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                    ii++;
                }

                ID = piaacmc[0].piaa1FmodesID;
                IDref = image_ID("piaa1Fmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) *pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                    ii++;
                }
            }


            if(REGFPMSAG == 1)
            {
                ID = piaacmc[0].zonezID;
                for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                    ii++;
                }
            }

            delete_image_ID("imvect");
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 29;


            

       




            // print out current best value for tracking
            sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
            fp = fopen(fname, "a");
            if(fp==NULL)
            {
                printf("ERROR: cannot open file \"%s\"\n", fname);
                exit(0);
            }
            fprintf(fp, "-> %5ld    %20g <- %20g \n", iter, val, valref);
            printf("%5ld %20g %20g \n", iter, val, valref);
            fflush(stdout);
            fclose(fp);

            // save current best value and reference value in globals
            PIAACMCSIMUL_VAL = val;
            PIAACMCSIMUL_VALREF = valref;



            // Nominally if we're in this linear
            // optimization PIAACMC_fpmtype = 1, so the next line is not executed
            if(PIAACMC_fpmtype==0) // in the idealized PIAACMC case
                piaacmc[0].fpmaskamptransm = data.image[piaacmc[0].zoneaID].array.D[0]; // required to ensure that the new optimal focal plane mask transmission is written to disk

            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 30;


            // tracking diagnostics giving behavior of the modes by iteration
            sprintf(dirname, "%s_linopt", piaacmcconfdir);
            PIAAsimul_savepiaacmcconf(dirname); // staging area
            sprintf(command, "rsync -au --progress %s/* ./%s/", dirname, piaacmcconfdir);
            r = system(command);

            sprintf(command, "cp %s/piaa0Cmodes.fits %s/piaa0Cmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);
            sprintf(command, "cp %s/piaa0Fmodes.fits %s/piaa0Fmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);

            sprintf(command, "cp %s/piaa1Cmodes.fits %s/piaa1Cmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);
            sprintf(command, "cp %s/piaa1Fmodes.fits %s/piaa1Fmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);

            sprintf(command, "cp %s/piaacmcparams.conf %s/piaacmcparams.%04ld.conf", dirname, dirname, iter);
            r = system(command);

			if(file_exists(stopfile)==1)
				iterOK = 0;

            // Figure out if current loop should continue optimization
            // if optimization ends, then iterOK set to 0
            // if we've reached the allowed number of iterations
            if(iter==NBiter)
                iterOK = 0;
            if(iter>2)
            {
                if(val>0.98*oldval) // if after second iteration and we'be improved by less than 10%
                    iterOK = 0;
            }
            
            if(NBlinoptgain<3) // if we've stopped moving much
                iterOK = 0;


            // set up for next iteration
            oldval = val;
            iter++;

            printf("END OF LOOP ITERATION\n");
            fflush(stdout);
            // for state tracking and statistics
            data.image[IDstatus].array.U[0] = 31;
        }
        printf(" ============ END OF OPTIMIZATION LOOP ======= \n");
        // for state tracking and statistics
        data.image[IDstatus].array.U[0] = 32;
    } // end of if (LINOPT==1): done with the linear optimization





    //  PIAAsimul_savepiaacmcconf("piaacmc0");
    //  PIAAsimul_loadpiaacmcconf("piaacmc0");
    // PIAAsimul_savepiaacmcconf("piaacmc1");
    //exit(0);

    return 0;
}





/// @param[in] confindex	configuration index (sets name of directory for results)
/// @param[in] mode			operation to be executed
/*
    entry point for PIAACMCsimul from the cli
*/
int PIAACMCsimul_run(char *confindex, long mode)
{
    long i;
    FILE *fp;
    char fname[500];
    char fname1[500];
    char fnamebestval[500];
    double bestval = 1.0;
    int ret;
    char command[1000];
    long k;
    int fOK = 0;
    int bOK = 0;
    int zeroST=0;
    int r;
    int loopOK;
    char stopfile[500];
    long IDv;
    long IDbestsoltmp, IDbestsol;
    char fnamebestsol[500];
    int loopin = 0;
    struct timeval start, end;
    long secs_used,micros_used;

    double searchtime = 3600.0*10.0; // [second] default 10 hours


    IDbestsol = -1; // data array index of current best solution

    
    // read various cli variables, possibly setting globals
    if((IDv=variable_ID("PIAACMC_MASKRADLD"))!=-1)
        PIAACMC_MASKRADLD = data.variable[IDv].value.f;


    if((IDv=variable_ID("PIAACMC_FPMsectors"))!=-1)
        PIAACMC_FPMsectors = (long) data.variable[IDv].value.f+0.01;
    printf("PIAACMC_FPMsectors = %d\n", PIAACMC_FPMsectors);

    if((IDv=variable_ID("SCORINGMASKTYPE"))!=-1)
        SCORINGMASKTYPE = (long) data.variable[IDv].value.f+0.01;
    printf("SCORINGMASKTYPE = %d\n", SCORINGMASKTYPE);

    if((IDv=variable_ID("PIAACMC_save"))!=-1)
        PIAACMC_save = (long) data.variable[IDv].value.f+0.01;
    printf("PIAACMC_save = %d\n", PIAACMC_save);



    if((IDv=variable_ID("PIAACMC_resolved"))!=-1)
        computePSF_ResolvedTarget = (long) (data.variable[IDv].value.f+0.01);
    if((IDv=variable_ID("PIAACMC_extmode"))!=-1)
        computePSF_ResolvedTarget_mode = (long) (data.variable[IDv].value.f+0.01);


    printf("mode = %ld\n", mode);


    // mode 13: optimize focal plane mask zones only, setting the sag values for each mask zone
    // This outer loop is to choose more different starting points for the exec loop
    if(mode==13) // loop to keep looking for optimal solution
    {
        sprintf(fname, "searchtime.txt");
        fp = fopen(fname,"r");
        if(fp!=NULL)
        {
            r = fscanf(fp, "%lf\n", &searchtime);
            fclose(fp);
        }


        if(searchtime<0.1)
            loopOK = 0;
        else
            loopOK = 1;


        gettimeofday(&start, NULL);
        i = 0;

        // while not exceed searchtime or no stop file
        while((loopOK==1)&&(i<1000000))
        {
            // read in the real searchtime nominally set by the bash script
            sprintf(fname, "searchtime.txt");
            fp = fopen(fname,"r");
            if(fp!=NULL)
            {
                r = fscanf(fp, "%lf\n", &searchtime);
                fclose(fp);
            }


            loopin = 1; // loop has been initialized
            if((i<1))
                MODampl = 0.0; // MODampl is a global
            else
                MODampl = 1.0e-6*pow(ran1(),4.0); // pick amplitude for random optimization starting point




            // after the first iteration, half the time set zeroST=0
            // 1/4th the time when a best solution exists set zeroST = 1
            // 1/4th the time set zeroST = 2
            // zeroST is an information-only flag that does not control activity: it reflects
            //  the settings in each conditional
            // zeroST = 0 => starting point uses random distribution of FP mask zone sags
            // zeroST = 1 => starting point is a mask that has no sag: blank focal plane mask
            // zeroST = 2 => starting point is best solution found so far
            //
            // data.image[piaacmc[0].zonezID].array is ???? *******************************
            if((i>1)&&(ran1()>0.5))
            {
                if((ran1()>0.5)&&(IDbestsol!=-1))
                {
                    zeroST = 2; // starting point = optimal solution
                    // copy the best solution to the current zoneID of array of sags
                    for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                        data.image[piaacmc[0].zonezID].array.D[k] = data.image[IDbestsol].array.D[k];
                }
                else
                {
                    zeroST = 1; // starting point = 0
                    // zero out the current zoneID of array of sags
                    for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                        data.image[piaacmc[0].zonezID].array.D[k] = 0.0;
                }
            }
            else
                zeroST = 0;


            // zeroST = 3 => starting point is best solution found so far.  Same as zeroST=2
            // this flags that it's this value 'cause it's third iteration
            if(i==3)
            {
                zeroST = 3;
                for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                    data.image[piaacmc[0].zonezID].array.D[k] = data.image[IDbestsol].array.D[k];
                MODampl = 0.0;
            }

            // actually do the optmization
            PIAACMCsimul_exec(confindex, mode);
            bOK = 0; // initialize have better value flag for printing "best" in a nice place

            // if there is no best _solution_, load the current solution
            if(IDbestsol==-1)
            {
               sprintf(fnamebestsol, "%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.best.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

                printf("LOADING \"%s\"...\n", fnamebestsol);
                fflush(stdout);
                IDbestsol = load_fits(fnamebestsol, "fpmbestsol", 0);
            }

            // set the name of the stopfile
            sprintf(stopfile, "%s/stoploop13.txt", piaacmcconfdir);

            // on first iteration load the best _value_ if it exists
            if(i==0)
            {
               sprintf(fnamebestval, "%s/mode13_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.bestval.txt", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

                printf("READING FILE \"%s\"\n", fnamebestval);
                fflush(stdout);
                fp = fopen(fnamebestval, "r");
                if(fp != NULL)
                {
                    r = fscanf(fp, "%lf", &bestval); // this loads only the first value on the line
                    fclose(fp);
                }
            }


            printf("\n\n\n\n======= val = %g [%g]\n", PIAACMCSIMUL_VAL, bestval);
            fflush(stdout);

            if(PIAACMCSIMUL_VAL<bestval) // PIAACMCSIMUL_VAL was set in PIAACMCsimul_exec()
            {
                // we have a better solution!
                bOK = 1;
                bestval = PIAACMCSIMUL_VAL; // record it
                printf("============================================================   SAVING BEST MASK SOLUTION -> fpm_zonez.best.fits\n");
                fflush(stdout);
                
                // if a previous best solution has not been identified with an index, set its index
                // by loading the current best solution.  This probably never happens
                if(IDbestsol==-1)
                {
                   sprintf(fnamebestsol, "%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.best.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

                    IDbestsol = load_fits(fnamebestsol, "fpmbestsol", 0);
                }
                else // otherwise load the temporary best solution.  This is probably what always happens
                {
                    IDbestsoltmp = load_fits(fnamebestsol, "fpmbestsoltmp", 0);
                    for(k=0; k<data.image[IDbestsol].md[0].size[0]; k++)
                        data.image[IDbestsol].array.D[k] = data.image[IDbestsoltmp].array.D[k];
                    delete_image_ID("fpmbestsoltmp");
                }
                
                // fname1 is the name of the current solution, which is now the best solution
               sprintf(fname1, "%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

                // fnamebestsol is the name of the stored best solution, should always be the same
                // as the name in line 8599 (if(IDbestsol==-1)...)
               sprintf(fnamebestsol, "%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.best.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

                // copy the current solution to the best solution
                sprintf(command, "cp %s %s", fname1, fnamebestsol);
                ret = system(command);

                // write new best value in file
                fp = fopen(fnamebestval, "w");
                fprintf(fp, "%30g %d %04ld %02ld %03ld %04ld %03ld %02d %d %s %02d\n", bestval, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, piaacmc[0].focmNBzone, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
                fclose(fp);

                // advertise the existence of new best solution via file signaling.  Currently no listeners?
                sprintf(command, "touch %s/newbestsol.txt", piaacmcconfdir);
                r = system(command);
            }

            // add current solution (possibly not best) to the mode13...opt.txt file
            sprintf(fname, "%s/mode13_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.opt.txt", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

            // first time through, open mode13...opt.txt for additional writing
            // for additional writing.  possibly redundant on next line.
            if(fOK==0)
            {
                fp = fopen(fname, "a");
                fclose(fp);
                fOK = 1;
            }


            // open mode13...opt.txt for adding and write current value
            fp = fopen(fname, "a");
            fprintf(fp,"%10ld %20.5g   %16.5g -> %16.5g   (%16.5g) %d  [%12g %2d %12g %12g  %12g]", i, MODampl, PIAACMCSIMUL_VALREF, PIAACMCSIMUL_VAL, bestval, zeroST, CnormFactor, piaacmc[0].nblambda, optsyst[0].flux[0], SCORINGTOTAL, PIAACMCSIMUL_VAL0);
            if(bOK==1) // mark it as best if it is
                fprintf(fp, " BEST\n");
            else
                fprintf(fp, "\n");
            fclose(fp);

            // if(PIAACMCSIMUL_VAL>PIAACMCSIMUL_VALREF)
            //	exit(0);

            i++; // increment iteration counter (!!)

            // stop iterations if stopfile exists
            if(file_exists(stopfile)==1)
            {
                printf("FILE \"%s\" found\n", stopfile);
                loopOK = 0;
                sprintf(command, "rm %s", stopfile);
                ret = system(command);
            }
            else
                printf("File \"%s\" not found\n", stopfile);


            gettimeofday(&end, NULL);

            printf("start: %ld secs, %ld usecs\n", (long) start.tv_sec, (long) start.tv_usec);
            printf("end: %ld secs, %ld usecs\n", (long) end.tv_sec, (long) end.tv_usec);

            secs_used=(end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
            micros_used= ((secs_used*1000000) + end.tv_usec) - (start.tv_usec);

            fp=fopen("timeused.txt", "w");
            fprintf(fp, "%12.3f    %12.3f\n", 1.0e-6*micros_used, searchtime);
            fclose(fp);

            // check to see if time has run out
            if(micros_used > 1000000.0*searchtime) // searchtime is in seconds
                loopOK = 0; // stop loop flag
        }



        // initialize loop.  loopin is always set to 1 above.
        if(loopin == 1)
        {
            printf("piaacmcconfdir              : %s\n", piaacmcconfdir);
            fflush(stdout);

            printf("computePSF_ResolvedTarget   : %d\n", computePSF_ResolvedTarget);
            fflush(stdout);


            printf("computePSF_ResolvedTarget_mode   : %d\n", computePSF_ResolvedTarget_mode);
            fflush(stdout);

            printf("PIAACMC_FPMsectors   : %d\n", PIAACMC_FPMsectors);
            fflush(stdout);

            printf("(long) (10.0*PIAACMC_MASKRADLD+0.1)   : %ld\n", (long) (10.0*PIAACMC_MASKRADLD+0.1));
            fflush(stdout);

            printf("piaacmc[0].NBrings   : %ld\n", piaacmc[0].NBrings);
            fflush(stdout);

            printf("piaacmc[0].nblambda   : %d\n", piaacmc[0].nblambda);
            fflush(stdout);

            // copy current solution to best solution ************************** why?
            sprintf(fname1, "%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

            sprintf(fnamebestsol, "%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.best.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

            sprintf(command, "cp %s %s", fnamebestsol, fname1);
            printf("COMMAND: %s\n", command);
            fflush(stdout);

            ret = system(command);
        }
    }
    else
        PIAACMCsimul_exec(confindex, mode);


    free(piaacmc);

    return 0;
}










