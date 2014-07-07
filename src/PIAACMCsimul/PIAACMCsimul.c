#include <fitsio.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>



#include "CLIcore.h"
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




extern DATA data;

#define SBUFFERSIZE 2000

///  Current configuration directory
char piaacmcconfdir[200];

OPTSYST *optsyst;
int optsystinit = 0;
long IDx, IDy, IDr, IDPA;

// this makes 20% bandwidth, from 0.55/1.1 to 0.55*1.1
double LAMBDASTART = 0.5e-6;
double LAMBDAEND = 0.605e-6;
#define NBLAMBDA 5

MIRRORPIAACMCDESIGN *piaacmc;
OPTSYST *optsyst;


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

int invPIAAmode = 1; // 0: no inv PIAA, 1: inv PIAA after Lyot stops, 2: inv PIAA before Lyot stops



// declared here for speed
double evalval;
long evali;
long evalk, evalki, evalki1, evalmz, evalii, evalii1, evalii2, evalkv;
double evalcosp, evalsinp, evalre, evalim, evalre1, evalim1, evalpha;
double evalv1;
double PIAACMCSIMUL_VAL;
double PIAACMCSIMUL_VALREF;

// for minimization
float *fpmresp_array;
double *zonez_array;
double *zonez0_array;
double *zonez1_array;
double *zonezbest_array;
double *dphadz_array;
float *outtmp_array;
long NBoptVar;
static long LOOPCNT = 0;
long vsize;
double cval0;

double CnormFactor = 1.0; // for contrast normalization
double THICKRANGE = 2.0e-6;

int computePSF_FAST_FPMresp = 0;
int computePSF_ResolvedTarget = 0;
int PIAACMC_FPM_FASTDERIVATIVES = 0;


long NBsubPix = 8;

double SCORINGTOTAL = 1.0;
double MODampl = 1.0e-6;


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//

int PIAACMCsimul_run_cli()
{

    if(CLI_checkarg(1,2)+CLI_checkarg(2,2)==0)
    {
        PIAACMCsimul_run(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl);
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



    strcpy(data.cmd[data.NBcmd].key,"piaacmcsimrun");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = PIAACMCsimul_run_cli;
    strcpy(data.cmd[data.NBcmd].info,"Simulate PIAACMC");
    strcpy(data.cmd[data.NBcmd].syntax,"<configuration index [long]>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmcsimrun");
    strcpy(data.cmd[data.NBcmd].Ccall,"int PIAACMCsimul_run(int confindex)");
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
    long NBzones;
    long ID;
    double x, y, r, PA;
    long ii, jj;
    long zi;
    long *sizearray;

	int sectors = 1; // build sectors as well as rings
	long ring;
	long *nbsector;
	long *nbsectorcumul;
	double PAf;
	double eps = 1.0e-6;
	
    sizearray = (long*) malloc(sizeof(long)*2);
    sizearray[0] = piaacmc[0].fpmarraysize;
    sizearray[1] = piaacmc[0].fpmarraysize;
    ID = create_image_ID(IDname, 2, sizearray, USHORT, 0, 0);
    free(sizearray);

	nbsector = (long*) malloc(sizeof(long)*piaacmc[0].NBrings);
	nbsectorcumul = (long*) malloc(sizeof(long)*piaacmc[0].NBrings);

	if(sectors==1)
	{
		NBzones = 0;
		nbsector[0] = 1;
		nbsectorcumul[0] = 1;
		for(ring=1;ring<piaacmc[0].NBrings;ring++)
		{
			nbsector[ring] = 2*(ring+1);
			nbsectorcumul[ring] = nbsectorcumul[ring-1]+nbsector[ring];
		}
	}

	for(ring=0;ring<piaacmc[0].NBrings;ring++)
	printf("ring %ld : %ld %ld\n", ring, nbsector[ring], nbsectorcumul[ring]);

    //  ID = create_2Dimage_ID(IDname, piaacmc[0].fpmarraysize, piaacmc[0].fpmarraysize);
    for(ii=0; ii<piaacmc[0].fpmarraysize; ii++)
        for(jj=0; jj<piaacmc[0].fpmarraysize; jj++)
        {
            x = (2.0*ii-1.0*piaacmc[0].fpmarraysize)/piaacmc[0].fpmarraysize;
            y = (2.0*jj-1.0*piaacmc[0].fpmarraysize)/piaacmc[0].fpmarraysize;
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
			
			
			if(sectors==0)
				data.image[ID].array.U[jj*piaacmc[0].fpmarraysize+ii] = (unsigned short int) zi;
			else
			{
				if(ring==-1)
					data.image[ID].array.U[jj*piaacmc[0].fpmarraysize+ii] = 0;
				else
				{
					if(ring==0) // inner disk
						data.image[ID].array.U[jj*piaacmc[0].fpmarraysize+ii] = 1;
					else
						{
							data.image[ID].array.U[jj*piaacmc[0].fpmarraysize+ii] = (unsigned short int) nbsectorcumul[ring-1]+1;
							data.image[ID].array.U[jj*piaacmc[0].fpmarraysize+ii] += (unsigned short int) (PAf*nbsector[ring]);
						}
				}
			}	
        }
        
	if(sectors==0)
		piaacmc[0].focmNBzone = piaacmc[0].NBrings;
	else
		piaacmc[0].focmNBzone = nbsectorcumul[piaacmc[0].NBrings-1];
		
	free(nbsector);
	free(nbsectorcumul);

    return ID;
}




//
// makes 1-fpm CA
// if mode = -1, make whole 1-fpm
// if mode = zone, make only 1 zone with CA = (1.0, 0.0)
// zone numbering starts here from 1 (zone 1 = outermost ring)
//
long PIAACMCsimul_mkFocalPlaneMask(char *IDzonemap_name, char *ID_name, int mode)
{
    long ID;
    long IDz;
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
	double retmp, imtmp;

    size = optsyst[0].size;
    size2 = size*size;
    nblambda = optsyst[0].nblambda;


    IDz = image_ID(IDzonemap_name);
    ID = create_3DCimage_ID(ID_name, size, size, nblambda);

    // CORONAGRAPHS_TDIAM/CORONAGRAPHS_PSCALE/CORONAGRAPHS_ARRAYSIZE

    /*  printf("%ld %ld\n", piaacmc[0].zoneaID, piaacmc[0].zonezID);
    for(k=0;k<data.image[piaacmc[0].zonezID].md[0].size[0];k++)
    {
    t = data.image[piaacmc[0].zonezID].array.D[k]; // thickness
    a = data.image[piaacmc[0].zoneaID].array.D[k]; // amplitude transmission
    amp = a;
    pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[0]);

    re = amp*cos(pha);
    im = amp*sin(pha);

    printf("ZONE %2ld  %12g %12g   %12g %12g    %12g %12g\n", k, a, t, amp, pha, 1.0-re, -im);
    }
    */


# ifdef HAVE_LIBGOMP
    #pragma omp parallel default(shared) private(ii, jj, x, y, r, retmp, imtmp, iii, jjj, ii1, jj1, zi, t, a, fpscale, amp, pha)
    {
        #pragma omp for
# endif
    for(k=0; k<nblambda; k++)
    {
        fpscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size/piaacmc[0].fpzfactor*optsyst[0].lambdaarray[k]*piaacmc[0].Fratio;
        printf("LAMBDA = %10.5g m    SCALE = %10.5g m/pix   size=%4ld  rad=%g\n", optsyst[0].lambdaarray[k], fpscale, size, piaacmc[0].fpmRad);


        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                x = (1.0*ii-size/2)*fpscale; // [m]
                y = (1.0*jj-size/2)*fpscale; // [m]
                r = sqrt(x*x+y*y); // [m]
	
				retmp = 0.0;
				imtmp = 0.0;
				if(r<1.1*piaacmc[0].fpmRad)
				for(iii=0;iii<NBsubPix;iii++)
					for(jjj=0;jjj<NBsubPix;jjj++)
					{
						// x and y in [m]
						x = (1.0*ii-size/2+1.0*(0.5+iii)/NBsubPix-0.5)*fpscale;
						y = (1.0*jj-size/2+1.0*(0.5+jjj)/NBsubPix-0.5)*fpscale;


						ii1 = (long) ( (0.5 + 0.5*x/piaacmc[0].fpmRad)*piaacmc[0].fpmarraysize + 0.5);
						jj1 = (long) ( (0.5 + 0.5*y/piaacmc[0].fpmRad)*piaacmc[0].fpmarraysize + 0.5);
						if((ii1>-1)&&(ii1<piaacmc[0].fpmarraysize)&&(jj1>-1)&&(jj1<piaacmc[0].fpmarraysize))
						{
							zi = data.image[IDz].array.U[jj1*piaacmc[0].fpmarraysize+ii1];
							t = data.image[piaacmc[0].zonezID].array.D[zi-1]; // thickness
							a = data.image[piaacmc[0].zoneaID].array.D[zi-1]; // amplitude transmission
						}
						else
						{
							zi = 0;
							t = 0.0;
							a = 1.0;
						}

						if(mode == -1)
						{
							if(zi>0.1)
							{
								amp = a;
								pha = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[k]);
								//pha = 0.0; //TEST
								retmp += 1.0-amp*cos(pha);
								imtmp += amp*sin(pha);
//								data.image[ID].array.CF[k*size2+jj*size+ii].re = 1.0-re;
//								data.image[ID].array.CF[k*size2+jj*size+ii].im = -im;
							}
							else
							{
//								data.image[ID].array.CF[k*size2+jj*size+ii].re = 0.0;
//								data.image[ID].array.CF[k*size2+jj*size+ii].im = 0.0;
							}
						}
						else // impusle response from single zone
						{
							if(mode == zi)
							{
								amp = 1.0;
								pha = 0.0; //OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t, optsyst[0].lambdaarray[k]);
								retmp += amp*cos(pha);
								imtmp += amp*sin(pha);
//								data.image[ID].array.CF[k*size2+jj*size+ii].re = 1.0;
//								data.image[ID].array.CF[k*size2+jj*size+ii].im = 0.0;
							}
							else
							{
//								data.image[ID].array.CF[k*size2+jj*size+ii].re = 0.0;
//								data.image[ID].array.CF[k*size2+jj*size+ii].im = 0.0;
							}
						}

					}
			data.image[ID].array.CF[k*size2+jj*size+ii].re = retmp/(NBsubPix*NBsubPix);
			data.image[ID].array.CF[k*size2+jj*size+ii].im = imtmp/(NBsubPix*NBsubPix);
			}
	}
# ifdef HAVE_LIBGOMP
    }
# endif

	

    return(ID);
}


















/// initializes the optsyst structure to simulate reflective PIAACMC system

void PIAACMCsimul_init( MIRRORPIAACMCDESIGN *design, long index, double TTxld, double TTyld )
{
    FILE *fp;
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
    char fname_pupa0[200];
    long ID;
    long elem;
    char fname[200];
    long IDv;

    optsyst[0].nblambda = design[index].nblambda;
    nblambda = optsyst[0].nblambda;

    sprintf(fname, "%s/lambdalist.txt", piaacmcconfdir);
    fp = fopen(fname, "w");
    for(k=0; k<optsyst[0].nblambda; k++)
    {
        optsyst[0].lambdaarray[k] = LAMBDASTART + (0.5+k)*(LAMBDAEND-LAMBDASTART)/optsyst[0].nblambda;
        fprintf(fp, "%02ld %20g\n", k, optsyst[0].lambdaarray[k]);
    }
    fclose(fp);

    optsyst[0].beamrad = design[index].beamrad; // 8mm
    optsyst[0].size = design[index].size;
    size = optsyst[0].size;
    size2 = size*size;
    optsyst[0].pixscale = design[index].pixscale;
    optsyst[0].DFTgridpad = 0; // 0 for full DFT sampling, >0 for faster execution

    beamradpix = optsyst[0].beamrad/optsyst[0].pixscale;
    printf("BEAM RADIUS = %f pix\n", beamradpix);


    if((IDv=variable_ID("PIAACMC_invPIAAmode"))!=-1)
        invPIAAmode = (long) (data.variable[IDv].value+0.001);

    if((IDv=variable_ID("PIAACMC_dftgrid"))!=-1)
        optsyst[0].DFTgridpad = (long) (data.variable[IDv].value+0.001);


    // define optical elements and locations
    optsyst[0].NB_DM = 0;
    optsyst[0].NB_asphsurfm = 2;
    optsyst[0].NB_asphsurfr = 0;

    optsyst[0].NBelem = 100; // to be updated later

    sprintf(fname, "%s/conjugations.txt", piaacmcconfdir);
    fp = fopen(fname, "w");

    elem = 0;
    // ------------------- elem 0: input pupil -----------------------

    optsyst[0].elemtype[elem] = 1; // pupil mask
    // input pupil
    sprintf(fname_pupa0, "%s/pupa0_%ld.fits", piaacmcconfdir, size);

    if(file_exists(fname_pupa0)==1)
        load_fits(fname_pupa0, "pupa0");

    IDa = image_ID("pupa0");
    if(IDa==-1)
    {
        printf("CREATING INPUT PUPIL\n");
        if(IDa!=-1)
            delete_image_ID("pupa0");
        IDa = create_3Dimage_ID("pupa0", size, size, nblambda);

        ID = image_ID("telpup");
        if(ID==-1)
            if(file_exists("telpup.fits")==1)
                ID = load_fits("telpup.fits", "telpup");


        if(ID==-1)
        {
            for(k=0; k<nblambda; k++)
                for(ii=0; ii<size2; ii++)
                {
                    if((data.image[IDr].array.F[ii]>piaacmc[0].centObs0)&&(data.image[IDr].array.F[ii]<1.0))
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
    optsyst[0].elemZpos[elem] = 0.0;
    fprintf(fp,"%02ld  %f    Input pupil\n", elem, optsyst[0].elemZpos[elem]);
    elem++;





    // pointing (simulated as mirror)
    ID = create_2Dimage_ID("TTm", size, size);

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = (1.0*ii-0.5*size)/beamradpix;
            y = (1.0*jj-0.5*size)/beamradpix;
            data.image[ID].array.F[jj*size+ii] = 0.25*(TTxld*x+TTyld*y)*(LAMBDAEND+LAMBDASTART)*0.5;
        }

    sprintf(fname, "!%s/TTm.fits", piaacmcconfdir);
    save_fits("TTm", fname);

    optsyst[0].elemtype[elem] = 3; // reflective mirror
    optsyst[0].elemarrayindex[elem] = 0; // index
    optsyst[0].ASPHSURFMarray[0].surfID = ID;
    optsyst[0].elemZpos[elem] = 0.0;
    fprintf(fp,"%02ld  %f    Fold mirror used to induce pointing offsets\n", elem, optsyst[0].elemZpos[elem]);
    elem++;


    IDpiaaz0 = image_ID("piaa0z");
    IDpiaaz1 = image_ID("piaa1z");


    // ------------------- elem 2: reflective PIAA M0  -----------------------
    optsyst[0].elemtype[elem] = 3; // reflective PIAA M0
    optsyst[0].elemarrayindex[elem] = 1; // index
    optsyst[0].ASPHSURFMarray[1].surfID = IDpiaaz0;
    optsyst[0].elemZpos[elem] = 0.0;
   fprintf(fp,"%02ld  %f    PIAAM0\n", elem, optsyst[0].elemZpos[elem]);
    elem++;

    // ------------------- elem 3: reflective PIAA M1  -----------------------
    optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
    optsyst[0].elemarrayindex[elem] = 2;
    optsyst[0].ASPHSURFMarray[2].surfID = IDpiaaz1;
    optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1]+design[index].piaasep;
   fprintf(fp,"%02ld  %f    PIAAM1\n", elem, optsyst[0].elemZpos[elem]);
    elem++;


    // ------------------- elem 4 opaque mask at reflective PIAA M1  -----------------------
    optsyst[0].elemtype[elem] = 1; // opaque mask
    ID = make_disk("piaam1mask", size, size, 0.5*size, 0.5*size, design[index].r1lim*beamradpix);
    optsyst[0].elemarrayindex[elem] = ID;
    optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1];
    sprintf(fname, "!%s/piaam1mask.fits", piaacmcconfdir);
    save_fits("piaam1mask", fname);
    fprintf(fp,"%02ld  %f    PIAAM1 edge opaque mask\n", elem, optsyst[0].elemZpos[elem]);
    elem++;


    // --------------------  elem 5: focal plane mask ------------------------
    optsyst[0].elemtype[elem] = 5; // focal plane mask
    optsyst[0].elemarrayindex[elem] = 0;

    optsyst[0].FOCMASKarray[0].fpmID = PIAACMCsimul_mkFocalPlaneMask("fpmzmap", "piaacmcfpm", focmMode); // if -1, this is 1-fpm; otherwise, this is impulse response from single zone

    if(0)// testing
    {
        sprintf(fname, "!focma_%d.fits", focmMode);
        mk_amph_from_complex("piaacmcfpm", "fpma", "fpmp");
        save_fits("fpma", fname);
        save_fits("fpmp", "!fpmp.fits");
        delete_image_ID("fpma");
        delete_image_ID("fpmp");
    }



    optsyst[0]. FOCMASKarray[0].zfactor = design[index].fpzfactor;
    optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1]; // plane from which FT is done
    fprintf(fp,"%02ld  %f    post-focal plane mask pupil\n", elem, optsyst[0].elemZpos[elem]);
    elem++;






    if(invPIAAmode == 2) // inv PIAA -> Lyot stops
    {
        // --------------------  elem 8: inv PIAA1 ------------------------
        optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
        optsyst[0].elemarrayindex[elem] = 2;
        optsyst[0].ASPHSURFMarray[2].surfID = IDpiaaz1;
        optsyst[0].elemZpos[elem] = 0.0;
       fprintf(fp,"%02ld  %f    invPIAA1\n", elem, optsyst[0].elemZpos[elem]);
        elem++;

        // --------------------  elem 9: inv PIAA0 ------------------------
        optsyst[0].elemtype[elem] = 3; // reflective PIAA M0
        optsyst[0].elemarrayindex[elem] = 1;
        optsyst[0].ASPHSURFMarray[1].surfID = IDpiaaz0;
        optsyst[0].elemZpos[elem] = design[index].piaasep;
       fprintf(fp,"%02ld  %f    invPIAA0\n", elem, optsyst[0].elemZpos[elem]);
        elem++;
    }



    // --------------------  Lyot masks  ------------------------
    for(i=0; i<design[index].NBLyotStop; i++)
    {
        optsyst[0].elemtype[elem] = 1; // Lyot mask
        optsyst[0].elemarrayindex[elem] = design[index].IDLyotStop[i];
        printf("elem %ld  Lyot mask %ld : %ld\n", elem, i, design[index].IDLyotStop[i]);
        optsyst[0].elemZpos[elem] =  design[index].LyotStop_zpos[i];
       fprintf(fp,"%02ld  %f  Lyot Stop %ld\n", elem, optsyst[0].elemZpos[elem], i);
        elem++;
    }



    if(invPIAAmode == 1) // Lyot masks -> inv PIAA
    {
        // --------------------  elem 8: inv PIAA1 ------------------------
        optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
        optsyst[0].elemarrayindex[elem] = 2;
        optsyst[0].ASPHSURFMarray[2].surfID = IDpiaaz1;
        optsyst[0].elemZpos[elem] = 0.0;
        fprintf(fp,"%02ld  %f    invPIAA1\n", elem, optsyst[0].elemZpos[elem]);
        elem++;

        // --------------------  elem 9: inv PIAA0 ------------------------
        optsyst[0].elemtype[elem] = 3; // reflective PIAA M0
        optsyst[0].elemarrayindex[elem] = 1;
        optsyst[0].ASPHSURFMarray[1].surfID = IDpiaaz0;
        optsyst[0].elemZpos[elem] = design[index].piaasep;
        fprintf(fp,"%02ld  %f    invPIAA0\n", elem, optsyst[0].elemZpos[elem]);
        elem++;
    }



    // --------------------  elem 9: back end mask  ------------------------

    optsyst[0].elemtype[elem] = 1; // Lyot mask
    ID = make_disk("outmask", size, size, 0.5*size, 0.5*size, 0.92*design[index].beamrad/design[index].pixscale);
    optsyst[0].elemarrayindex[elem] = ID;
    optsyst[0].elemZpos[elem] =  design[index].piaasep;
    fprintf(fp,"%02ld  %f   back end mask\n", elem, optsyst[0].elemZpos[elem]);
    elem++;

    fclose(fp);

    optsyst[0].NBelem = elem;


    optsystinit = 1;
}











//
// RADIAL PIAACMC SYSTEM DESIGN (geometrical optics)
//



//
// load and fit radial apodization profile
// modal basis is mk(r) : cos(r*k*M_PI/1.3)
//
int PIAACMCsimul_load2DRadialApodization(char *IDapo_name, float beamradpix, float centralObs, char *IDapofit_name)
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
    char fname[200];
    int ret;
    char command[500];

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


    sprintf(fname, "!%s/_apoincrop.fits", piaacmcconfdir);
    save_fits("_apoincrop", fname);

    sprintf(fname, "!%s/fitmaskapo.fits", piaacmcconfdir);
    save_fits("fitmaskapo", fname);


    linopt_imtools_image_fitModes("_apoincrop", "APOmodesCos", "fitmaskapo", 1.0e-8, IDapofit_name, 0);
    sprintf(command, "mv %s/eigenv.dat %s/eigenv_APOmodesCos.dat", piaacmcconfdir, piaacmcconfdir);
    ret = system(command);

    if(0) // test fit quality
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
    double *piaar0;
    double *piaar1;

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

    char fname[200];

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





    printf("======== Compute PIAA mirror shapes ============\n");
    piaaM0z = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    piaaM1z = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);

    piaaM0z[0] = 0.0;
    piaaM1z[0] = piaacmc[0].piaasep;


    for(i=0; i<piaacmc[0].NBradpts-1; i++)
    {
        r0c = piaar00[i];
        r1c = piaar10[i];
        dx = (r0c-r1c)*piaacmc[0].beamrad;
        dz = piaaM1z[i]-piaaM0z[i];
        dist = sqrt(dx*dx+dz*dz);
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

    free(piaar0);
    free(piaar1);
    free(piaar00);
    free(piaar10);
    free(piaar01);
    free(piaar11);


    return(0);
}


















//
// make PIAA shapes from radial sag profile
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
    printf("SIZE = %ld, beamrad = %f pix, sep = %f m\n", size, beamradpix, piaacmc[0].piaasep);
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
        z1array[k] -= piaacmc[0].piaasep;




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
                    data.image[ID_PIAAM1].array.F[jj*size+ii] = -val;//-piaacmc[0].piaasep);
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
 * @param[in] piaacmctype  Type of system
 * @param[in] fpmradld     Focal plane mask nominal radius
 * @param[in] centobs0     Input central obstruction
 * @param[in] centobs1     Output central obstruction
 * @param[in] load         if 1, attempt to load configuration from file
 *
 * piaacmctype:
 * - 0: if configuration does not exist, create Monochromatic idealized PIAACMC, otherwise, read configuration
 * - 1: change array size to variable PIAACMC_size
 *
 */

int PIAAsimul_initpiaacmcconf(long piaacmctype, double fpmradld, double centobs0, double centobs1, int load)
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
    char command[500];
    long IDv1, IDv2;
    char fname[200];
    char name[200];
	float tmpf;
    double pha0, t0, t;
    int loaded;

    int saveconf = 0; // if 1, save conf at end of this function
    long IDv;
    int ret;

    int IDlscumul;

    if(piaacmc == NULL)
    {
        piaacmc = (MIRRORPIAACMCDESIGN*) malloc(sizeof(MIRRORPIAACMCDESIGN)*NBpiaacmcdesign);


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


        piaacmc[0].beamrad = 0.022; // 44 mm diam
        piaacmc[0].piaasep = 2.302606; // [m]
        piaacmc[0].fpzfactor = 16.0;
        piaacmc[0].Fratio = 80.0;

        piaacmc[0].centObs0 = centobs0; // input central obstruction
        piaacmc[0].centObs1 = centobs1; // output central obstruction
        piaacmc[0].NBradpts = 50000;
        piaacmc[0].r0lim = 1.15; //1425; // outer radius after extrapolation, piaa mirror 0
        piaacmc[0].r1lim = 1.5; // outer radius after extrapolation, piaa mirror 1

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
        piaacmc[0].fpmmaterial = 0;  // 0: mirror, 4: PMMA
        piaacmc[0].fpmaskamptransm = 1.0;


        if(piaacmctype==0) // idealized focal plane mask
        {
            piaacmc[0].NBrings = 1;
            piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio;
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

    if((IDv=variable_ID("PIAACMC_nblambda"))!=-1)
        piaacmc[0].nblambda = data.variable[IDv].value;

    if((IDv=variable_ID("PIAACMC_NBrings"))!=-1)
        piaacmc[0].NBrings = data.variable[IDv].value;

    if((IDv=variable_ID("PIAACMC_size"))!=-1)
        piaacmc[0].size = (long) (data.variable[IDv].value+0.01);

    if((IDv=variable_ID("PIAACMC_pixscale"))!=-1)
        piaacmc[0].pixscale = data.variable[IDv].value;


    // create modes for aspheric optical surfaces description
    beamradpix = piaacmc[0].beamrad/piaacmc[0].pixscale;
    size = piaacmc[0].size;

    printf("BEAM RADIUS:  %f / %f =  %f pix\n", piaacmc[0].beamrad, piaacmc[0].pixscale, beamradpix);

    // x, y, r and PA coordinates in beam (for convenience & speed)
    IDx = create_2Dimage_ID("xcoord", size, size);
    IDy = create_2Dimage_ID("ycoord", size, size);
    IDr = create_2Dimage_ID("rcoord", size, size);
    IDPA = create_2Dimage_ID("PAcoord", size, size);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = (1.0*ii-0.5*size)/beamradpix;
            y = (1.0*jj-0.5*size)/beamradpix;
            data.image[IDx].array.F[jj*size+ii] = x;
            data.image[IDy].array.F[jj*size+ii] = y;
            data.image[IDr].array.F[jj*size+ii] = sqrt(x*x+y*y);
            data.image[IDPA].array.F[jj*size+ii] = atan2(y,x);
        }



    // ==================== CREATE MODES USED TO FIT AND DESCRIBE PIAA SHAPES ===============
    CREATE_Cmodes = 0;
    //   sprintf(fname, "%s/Cmodes.fits", piaacmcconfdir);
    sprintf(fname, "Cmodes.fits");
    if(FORCE_CREATE_Cmodes==0)
    {
        piaacmc[0].CmodesID = image_ID("Cmodes");
        if(piaacmc[0].CmodesID==-1)
            piaacmc[0].CmodesID = load_fits(fname, "Cmodes");
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
        piaacmc[0].Cmsize = Cmsize;
        linopt_imtools_makeCosRadModes("Cmodes", Cmsize, 40, ApoFitCosFact*beamradpix, 2.0);
        piaacmc[0].CmodesID = image_ID("Cmodes");
        save_fits("Cmodes", fname);
    }
    piaacmc[0].NBCmodes = data.image[piaacmc[0].CmodesID].md[0].size[2];
    piaacmc[0].Cmsize = data.image[piaacmc[0].CmodesID].md[0].size[0];

    CREATE_Fmodes = 0;
    //    sprintf(fname, "%s/Fmodes.fits", piaacmcconfdir);
    sprintf(fname, "Fmodes.fits");
    if(FORCE_CREATE_Fmodes == 0)
    {
        piaacmc[0].FmodesID = image_ID("Fmodes");
        if(piaacmc[0].FmodesID==-1)
            piaacmc[0].FmodesID = load_fits(fname, "Fmodes");
        if(piaacmc[0].FmodesID==-1)
            CREATE_Fmodes = 1;
    }
    else
        CREATE_Fmodes = 1;
    if(CREATE_Fmodes == 1)
    {
        Fmsize = (long) (beamradpix*4);
        piaacmc[0].Fmsize = Fmsize;
        linopt_imtools_makeCPAmodes("Fmodes",  Fmsize, 10.0, 0.8, beamradpix, 2.0, 1);
        piaacmc[0].FmodesID = image_ID("Fmodes");
        save_fits("Fmodes", fname);
        sprintf(command, "mv ModesExpr_CPA.txt %s/", piaacmcconfdir);
        r = system(command);
    }
    piaacmc[0].NBFmodes = data.image[piaacmc[0].FmodesID].md[0].size[2];
    piaacmc[0].Fmsize = data.image[piaacmc[0].FmodesID].md[0].size[0];




    // =================== IMPORT / CREATE PIAA SHAPES =====================

    piaacmc[0].piaa0CmodesID = image_ID("piaa0Cmodescoeff");
    piaacmc[0].piaa0FmodesID = image_ID("piaa0Fmodescoeff");
    piaacmc[0].piaa1CmodesID = image_ID("piaa1Cmodescoeff");
    piaacmc[0].piaa1FmodesID = image_ID("piaa1Fmodescoeff");

    sprintf(command, "mkdir -p %s/piaaref/", piaacmcconfdir);
    ret = system(command);

	if((piaacmc[0].piaa0CmodesID==-1)||( piaacmc[0].piaa0FmodesID==-1)||(piaacmc[0].piaa1CmodesID==-1)||( piaacmc[0].piaa1FmodesID==-1))
		{
	     sprintf(fname, "%s/piaaref/piaa0Cmodes.fits", piaacmcconfdir);
		piaacmc[0].piaa0CmodesID= load_fits(fname, "piaa0Cmodescoeff");

        sprintf(fname, "%s/piaaref/piaa0Fmodes.fits", piaacmcconfdir);
        piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff");

        sprintf(fname, "%s/piaaref/piaa1Cmodes.fits", piaacmcconfdir);
		piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff");

        sprintf(fname, "%s/piaaref/piaa1Fmodes.fits", piaacmcconfdir);
        piaacmc[0].piaa1FmodesID =load_fits(fname, "piaa1Fmodescoeff");
        
        sprintf(fname, "%s/piaaref/APLCmaskCtransm.txt", piaacmcconfdir);
        fp = fopen(fname, "r");
        if(fp!=NULL)
        {
			ret = fscanf(fp, "%f", &tmpf);
            piaacmc[0].fpmaskamptransm = tmpf;
            fclose(fp);
		}
		}


    if((piaacmc[0].piaa0CmodesID==-1)||( piaacmc[0].piaa0FmodesID==-1)||(piaacmc[0].piaa1CmodesID==-1)||( piaacmc[0].piaa1FmodesID==-1))
    {
        sprintf(fname, "%s/apo2Drad.fits", piaacmcconfdir);
        if(load_fits(fname, "apo2Drad")==-1)
        {
            sprintf(command, "cp %s/piaaref/apo2Drad.fits %s/apo2Drad.fits", piaacmcconfdir, piaacmcconfdir);
            ret = system(command);

            sprintf(fname, "%s/apo2Drad.fits", piaacmcconfdir);
            if(load_fits(fname, "apo2Drad")==-1)
            {

                printf("Creating 2D apodization for idealized circular monochromatic PIAACMC\n");
                fflush(stdout);

                // first iteration: half size image, 2x zoom
                IDv1 = create_variable_ID("DFTZFACTOR", 2);
                IDv2 = create_variable_ID("PNBITER", 15);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix*0.5, piaacmc[0].centObs1, "apotmp1", size/2);

                // expand solution to full size
                basic_resizeim("apotmp1", "apostart", size, size);
                delete_image_ID("apotmp1");

                // full size, 4x zoom
                IDv1 = create_variable_ID("DFTZFACTOR", 4);
                IDv2 = create_variable_ID("PNBITER", 5);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size);

                // full size, 8x zoom
				chname_image_ID("apo", "apostart");
                IDv1 = create_variable_ID("DFTZFACTOR", 8);
				IDv2 = create_variable_ID("PNBITER", 5);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size);


                // full size, 16x zoom
				chname_image_ID("apo", "apostart");
                IDv1 = create_variable_ID("DFTZFACTOR", 16);
                IDv2 = create_variable_ID("PNBITER", 10);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size);

                //  sprintf(command, "mv _DFT* %s/", piaacmcconfdir);
                //   r = system(command);
                //  sprintf(command, "mv APLCapo* %s/", piaacmcconfdir);
                //  r = system(command);
                // sprintf(command, "mv FPmask.tmp.fits %s/", piaacmcconfdir);
                // r = system(command);

                chname_image_ID("apo", "apo2Drad");
				
				sprintf(fname, "!%s/apo2Drad.fits", piaacmcconfdir);
                save_fits("apo2Drad", fname);
 
				sprintf(fname, "!%s/piaaref/apo2Drad.fits", piaacmcconfdir);
                save_fits("apo2Drad", fname);



                if((piaacmctype==0)&&(loaded==0)) // idealized focal plane mask
                {
                    piaacmc[0].fpmaskamptransm =  -data.variable[variable_ID("APLCmaskCtransm")].value;
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

        // load apodization profile and fit it a series of cosines
        PIAACMCsimul_load2DRadialApodization("apo2Drad", beamradpix, 0.3, "outApofit");

        // compute radial PIAA sag
        PIAACMCsimul_init_geomPIAA_rad("outApofit");

        // make 2D sag maps
        sprintf(fname, "%s/PIAA_Mshapes.txt", piaacmcconfdir);
        PIAACMCsimul_mkPIAAMshapes_from_RadSag(fname, "piaam0z", "piaam1z");

        sprintf(fname, "!%s/piaam0z.fits", piaacmcconfdir);
        save_fits("piaam0z", fname);

        sprintf(fname, "!%s/piaam1z.fits", piaacmcconfdir);
        save_fits("piaam1z", fname);



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

        sprintf(fname, "!%s/maskfit.fits", piaacmcconfdir);
        save_fits("maskfit", fname);

        printf("--------- FITTING COSINE MODES ---------\n");
        fflush(stdout);

        linopt_imtools_image_fitModes("piaa0zcrop", "Cmodes", "maskfit", 1.0e-6, "piaa0Cmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa0Cmodes.dat", piaacmcconfdir, piaacmcconfdir);
        ret = system(command);


        sprintf(fname, "!%s/piaa0Cmodescoeff.fits", piaacmcconfdir);
        save_fits("piaa0Cmodescoeff", fname);

        linopt_imtools_image_fitModes("piaa1zcrop", "Cmodes", "maskfit", 1.0e-6, "piaa1Cmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa1Cmodes.dat", piaacmcconfdir, piaacmcconfdir);
        ret = system(command);

        sprintf(fname, "!%s/piaa1Cmodescoeff.fits", piaacmcconfdir);
        save_fits("piaa1Cmodescoeff", fname);

        linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");

        sprintf(fname, "!%s/piaa0Cz.fits", piaacmcconfdir);
        save_fits("piaa0Cz", fname);

        linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");

        sprintf(fname, "!%s/piaa1Cz.fits", piaacmcconfdir);
        save_fits("piaa1Cz", fname);

        ID0 = image_ID("piaa0Cz");
        size0 = data.image[ID0].md[0].size[0];
        ID1 = image_ID("piaam0z");
        ID = create_2Dimage_ID("piaa0Cres", size0, size0);
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[jj*size0+ii] = data.image[ID1].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)]-data.image[ID0].array.F[jj*size0+ii];
        sprintf(fname, "!%s/piaa0Cres.fits", piaacmcconfdir);
        save_fits("piaa0Cres", fname);

        ID0 = image_ID("piaa1Cz");
        size0 = data.image[ID0].md[0].size[0];
        ID1 = image_ID("piaam1z");
        ID = create_2Dimage_ID("piaa1Cres", size0, size0);
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[jj*size0+ii] = data.image[ID1].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)]-data.image[ID0].array.F[jj*size0+ii];

        sprintf(fname, "!%s/piaa1Cres.fits", piaacmcconfdir);
        save_fits("piaa1Cres", fname);




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
        save_fits(data.image[piaacmc[0].piaa0CmodesID].md[0].name, fname);

        sprintf(fname, "!%s/piaaref/piaa0Fmodes.fits", piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa0FmodesID].md[0].name, fname);

        sprintf(fname, "!%s/piaaref/piaa1Cmodes.fits", piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa1CmodesID].md[0].name, fname);

        sprintf(fname, "!%s/piaaref/piaa1Fmodes.fits", piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa1FmodesID].md[0].name, fname);

        sprintf(command, "cp %s/piaaref/* %s/", piaacmcconfdir, piaacmcconfdir);
        ret = system(command);
    }




    // ============ MAKE FOCAL PLANE MASK ===============

    CREATE_fpmzmap = 0;
    if(FORCE_CREATE_fpmzmap == 0)
    {
        if(image_ID("fpmzmap")==-1)
            CREATE_fpmzmap = 1;
    }
    else
        CREATE_fpmzmap = 1;
    if(CREATE_fpmzmap == 1)
    {
        if(image_ID("fpmzmap")!=-1)
            delete_image_ID("fpmzmap");
        PIAACMCsimul_mkFPM_zonemap("fpmzmap");
        sprintf(fname, "!%s/fpmzmap.fits", piaacmcconfdir);
        save_fits("fpmzmap", fname);
    }

    CREATE_fpmzt = 0;
    if(FORCE_CREATE_fpmzt == 0)
    {
        piaacmc[0].zonezID = image_ID("fpmzt");
        if(piaacmc[0].zonezID == -1)
        {
            sprintf(fname, "%s/fpm_zonez.fits", piaacmcconfdir);
            piaacmc[0].zonezID = load_fits(fname, "fpmzt");
            if(piaacmc[0].zonezID == -1)
                CREATE_fpmzt = 1;
        }
    }
    else
        CREATE_fpmzt = 1;
    if(CREATE_fpmzt == 1)
    {
        printf("Creating fpmzt, saving as fpm_zonez.fits\n");
        piaacmc[0].zonezID = image_ID("fpmzt");
        if(piaacmc[0].zonezID!=-1)
            delete_image_ID("fpmzt");

        piaacmc[0].zonezID = create_2Dimagedouble_ID("fpmzt", piaacmc[0].focmNBzone, 1);
        t = 1.0e-9;
        if((piaacmctype==0)&&(loaded==0)) // idealized focal plane mask
        {
            t0 = 1.0e-8;
            pha0 = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, t0, 0.5*(LAMBDASTART+LAMBDAEND));
            t = (M_PI/pha0)*t0;
            printf("t = %g m (%lf %g)\n", t, pha0, t0);
        }
        for(ii=0; ii<piaacmc[0].focmNBzone; ii++)
            data.image[piaacmc[0].zonezID].array.D[ii] = t;
        sprintf(fname, "!%s/fpm_zonez.fits", piaacmcconfdir);
        save_fits("fpmzt", fname);
    }

    if(FORCE_CREATE_fpmza == 0)
    {
        piaacmc[0].zoneaID = image_ID("fpmza");
        if(piaacmc[0].zoneaID == -1)
        {
            sprintf(fname, "%s/fpm_zonea.fits", piaacmcconfdir);
            load_fits(fname, "fpmza");

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
        piaacmc[0].zoneaID = create_2Dimagedouble_ID("fpmza", piaacmc[0].focmNBzone, 1);
        for(ii=0; ii<piaacmc[0].focmNBzone; ii++)
            data.image[piaacmc[0].zoneaID].array.D[ii] = piaacmc[0].fpmaskamptransm;
        sprintf(fname, "!%s/fpm_zonea.fits", piaacmcconfdir);
        save_fits("fpmza", fname);
    }



    // ============= MAKE LYOT STOPS =======================
    printf("LOADING/CREATING LYOT MASK  - %ld masks\n", piaacmc[0].NBLyotStop);
    size2 = size*size;

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
            switch (i) {
            case 0 :
                piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, -0.01, 0.98);
                break;
            case 1 :
                piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, piaacmc[0].centObs1+0.02, 1.2);
                break;
            default :
                piaacmc[0].IDLyotStop[i] = PIAAsimul_mkSimpleLyotStop(name, piaacmc[0].centObs1+0.02, 0.98);
                break;
            }
     /*       ID = image_ID(name);
            for(ii=0; ii<size2; ii++)
            {
                data.image[ID].array.F[ii] *= data.image[IDlscumul].array.F[ii];
                data.image[IDlscumul].array.F[ii] = data.image[ID].array.F[ii];
            }*/
            save_fl_fits(name, fname);
        }
    }

    if(saveconf==1)
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);

    return(0);
}






int PIAACMCsimul_makePIAAshapes()
{
    long ID, ID0, ID1;
    long size, size0, size1;
    long ii, jj;
    char fname[200];

    size = piaacmc[0].size;

    // ============ construct PIAA shapes from fitting coefficients ==================

    MAKE_PIAA0shape = 0;
    if(FORCE_MAKE_PIAA0shape == 0)
    {
        ID = image_ID("piaa0z");
        if(ID==-1)
            MAKE_PIAA0shape = 1;
    }
    else
        MAKE_PIAA0shape = 1;

    if(MAKE_PIAA0shape == 1)
    {
        // assemble piaa0z and piaa1z images
        ID0 = linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");
        ID1 = linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
        ID = image_ID("piaa0z");
        if(ID==-1)
            ID = create_2Dimage_ID("piaa0z", size, size);

        printf("========================== STEP 01\n");
        fflush(stdout);


        size0 = data.image[ID0].md[0].size[0];
        size1 = data.image[ID1].md[0].size[0];
        for(ii=0; ii<size*size; ii++)
            data.image[ID].array.F[ii] = 0.0;

        printf("========================== STEP 01a  %ld %ld %ld\n", ID, ID0, ID1);
        printf(" %ld %ld\n", size0, size1);
        fflush(stdout);

        list_image_ID();

        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)] += data.image[ID0].array.F[jj*size0+ii];
        for(ii=0; ii<size1; ii++)
            for(jj=0; jj<size1; jj++)
                data.image[ID].array.F[(jj+(size-size1)/2)*size+(ii+(size-size1)/2)] += data.image[ID1].array.F[jj*size1+ii];

        sprintf(fname, "!%s/piaa0Cz.fits", piaacmcconfdir);
        save_fits("piaa0Cz", fname);

        sprintf(fname, "!%s/piaa0Fz.fits", piaacmcconfdir);
        save_fits("piaa0Fz", fname);

        sprintf(fname, "!%s/piaa0z.fits", piaacmcconfdir);
        save_fits("piaa0z", fname);

        delete_image_ID("piaa0Cz");
        delete_image_ID("piaa0Fz");
    }


    MAKE_PIAA1shape = 0;
    if(FORCE_MAKE_PIAA1shape == 0)
    {
        ID = image_ID("piaa1z");
        if(ID==-1)
            MAKE_PIAA1shape = 1;
    }
    else
        MAKE_PIAA1shape = 1;

    if(MAKE_PIAA1shape == 1)
    {
        ID0 = linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");
        ID1 = linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
        ID = image_ID("piaa1z");
        if(ID==-1)
            ID = create_2Dimage_ID("piaa1z", size, size);
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

        sprintf(fname, "!%s/piaa1Cz.fits", piaacmcconfdir);
        save_fits("piaa1Cz", fname);

        sprintf(fname, "!%s/piaaF1z.fits", piaacmcconfdir);
        save_fits("piaa1Fz", fname);

        sprintf(fname, "!%s/piaa1z.fits", piaacmcconfdir);
        save_fits("piaa1z", fname);

        delete_image_ID("piaa1Cz");
        delete_image_ID("piaa1Fz");
    }

    return 0;
}



///
/// returns average contrast in evaluation zone
///

double PIAACMCsimul_computePSF(float xld, float yld, long startelem, long endelem)
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

    char fname_piaa1z[200];
    char fname_piaa2z[200];
    char fname_pupa0[200];
    char fname_pupp0[200];
    char fname[200];

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

    size = piaacmc[0].size;
    size2 = size*size;


    focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;


    // CREATE SCORING MASK IF IT DOES NOT EXIST
    if((IDsm=image_ID("scoringmask"))==-1)
    {
        printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);
        IDsm = create_2Dimage_ID("scoringmask", size, size);
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                x = (1.0*ii-0.5*size)*focscale;
                y = (1.0*jj-0.5*size)*focscale;
                r = sqrt(x*x+y*y);
                
                if((r>scoringIWA)&&(r<scoringOWAhr)&&(x>scoringIWAx))
                    data.image[IDsm].array.F[jj*size+ii] = 1.0;

                 if((r>scoringIWA)&&(r<scoringOWA)&&(x>scoringIWAx)&&(ii%2==0)&&(jj%2==0))
                    data.image[IDsm].array.F[jj*size+ii] = 1.0;
                  
                if((x>scoringIWA)&&(fabs(y)<scoringIWA*0.5)&&(r<50.0)) // single line
                    data.image[IDsm].array.F[jj*size+ii] = 1.0;
            }

        sprintf(fname, "!%s/scoringmask.fits", piaacmcconfdir);
        save_fits("scoringmask", fname);

        linopt_imtools_mask_to_pixtable("scoringmask", "pixindex", "pixmult");

        SCORINGTOTAL = arith_image_total("scoringmask");



    }


    if(computePSF_ResolvedTarget==1)
    {
		dld = 0.01; // pointing offset [l/D]

        PIAACMCsimul_init(piaacmc, 0, xld-dld, yld-dld);
        PIAACMCsimul_makePIAAshapes();
        OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir);
        linopt_imtools_Image_to_vec("psfc", "pixindex", "pixmult", "imvectmm");

        PIAACMCsimul_init(piaacmc, 0, xld+dld, yld-dld);
        PIAACMCsimul_makePIAAshapes();
        OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir);
        linopt_imtools_Image_to_vec("psfc", "pixindex", "pixmult", "imvectpm");

        PIAACMCsimul_init(piaacmc, 0, xld-dld, yld+dld);
        PIAACMCsimul_makePIAAshapes();
        OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir);
        linopt_imtools_Image_to_vec("psfc", "pixindex", "pixmult", "imvectmp");

        PIAACMCsimul_init(piaacmc, 0, xld+dld, yld+dld);
        PIAACMCsimul_makePIAAshapes();
        OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir);
        linopt_imtools_Image_to_vec("psfc", "pixindex", "pixmult", "imvectpp");

        ID = image_ID("imvectpp");
		nbelem = data.image[ID].md[0].nelement;

        ID = image_ID("imvect");
		if(ID!=-1)
		delete_image_ID("imvect");
		
		ID = create_2Dimage_ID("imvect", nbelem*4, 1);
		
		ID1 = image_ID("imvectpp");
		for(ii=0; ii<nbelem; ii++)
			data.image[ID].array.F[ii] = data.image[ID1].array.F[ii]/2.0;

		ID1 = image_ID("imvectpm");
		for(ii=0; ii<nbelem; ii++)
			data.image[ID].array.F[nbelem+ii] = data.image[ID1].array.F[ii]/2.0;

		ID1 = image_ID("imvectmp");
		for(ii=0; ii<nbelem; ii++)
			data.image[ID].array.F[2*nbelem+ii] = data.image[ID1].array.F[ii]/2.0;
	
 		ID1 = image_ID("imvectmm");
		for(ii=0; ii<nbelem; ii++)
			data.image[ID].array.F[3*nbelem+ii] = data.image[ID1].array.F[ii]/2.0;

		delete_image_ID("imvectpp");
		delete_image_ID("imvectpm");
		delete_image_ID("imvectmp");
		delete_image_ID("imvectmm");
		
		  value = 0.0;
            ID = image_ID("imvect");
            for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
            {
                tmpv = data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
                value += tmpv;
            }
		           value = value/size/size/optsyst[0].flux[0];
          avContrast = value/(arith_image_total("scoringmask")*focscale*focscale);


   }
    else
    {
        if(computePSF_FAST_FPMresp==0)
        {
            // ========== initializes optical system to piaacmc design ===========
            PIAACMCsimul_init(piaacmc, 0, xld, yld);
            PIAACMCsimul_makePIAAshapes();

            // ============ perform propagations ================
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir);

            linopt_imtools_Image_to_vec("psfc", "pixindex", "pixmult", "imvect");
            //save_fits("imvect", "!imvect.fits");


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
            value = value/size/size/optsyst[0].flux[0];

            sprintf(fname,"%s/flux.txt", piaacmcconfdir);
            fpflux = fopen(fname, "w");
            for(elem=0; elem<optsyst[0].NBelem; elem++)
                fprintf(fpflux, "%18.16lf %18.16lf\n", optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0]);
            fclose(fpflux);


            avContrast = value/(arith_image_total("scoringmask")*focscale*focscale);

            CnormFactor = size*size*optsyst[0].flux[0]*arith_image_total("scoringmask")*focscale*focscale; // /optsyst[0].nblambda;
            CnormFactor = arith_image_total("scoringmask")*focscale*focscale*size*size*optsyst[0].flux[0]/optsyst[0].nblambda;
            sprintf(fname, "%s/CnormFactor.txt", piaacmcconfdir);
            fp = fopen(fname, "w");
            fprintf(fp, "%g\n", CnormFactor);
            fclose(fp);


            printf("Peak constrast (rough estimate)= %g\n", peakcontrast/size/size/optsyst[0].flux[0]/focscale/focscale*3.0);
            printf("Total light in scoring field = %g  -> Average contrast = %g\n", value, value/(arith_image_total("scoringmask")*focscale*focscale));
        }
        else
        {
            value1 = PIAACMCsimul_achromFPMsol_eval(fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], optsyst[0].nblambda);
            value = 0.0;
            peakcontrast = 0.0;

            ID = image_ID("imvect");
            if(ID==-1)
                ID = create_2Dimage_ID("imvect", vsize*optsyst[0].nblambda, 1);
            for(ii=0; ii<vsize*optsyst[0].nblambda; ii++)
            {
                data.image[ID].array.F[ii] = outtmp_array[ii];
                tmpv = outtmp_array[ii]*outtmp_array[ii];
                value += tmpv;
                //           if(tmpv>peakcontrast)
                //             peakcontrast = tmpv;
            }
            value = value/size/size/optsyst[0].flux[0];
            avContrast = value/(SCORINGTOTAL*focscale*focscale);
            //        printf("*********************************************************************************\n");
            //        printf("Peak constrast (rough estimate)= %g\n", peakcontrast/size/size/optsyst[0].flux[0]/focscale/focscale*3.0);
            //        printf("value1 = %g\n", value1);
            //		printf("Total light in scoring field = %g  -> Average contrast = %g   (%g)\n", value, value/(arith_image_total("scoringmask")*focscale*focscale), value1/CnormFactor/optsyst[0].nblambda);
        }
    }

    return(avContrast);
}







int PIAAsimul_savepiaacmcconf(char *dname)
{
    char command[200];
    int r;
    FILE *fp;
    char fname[200];
    long i;


    sprintf(command, "mkdir -p %s", dname);
    r = system(command);

    sprintf(fname,"%s/piaacmcparams.conf", dname);
    fp = fopen(fname, "w");


    fprintf(fp, "%10.6f   beamrad\n", piaacmc[0].beamrad);
    fprintf(fp, "%10ld    size\n", piaacmc[0].size);
    fprintf(fp, "%10.6g   pixscale\n", piaacmc[0].pixscale);
    fprintf(fp, "%10.6f   piaasep\n", piaacmc[0].piaasep);

    fprintf(fp, "%10.6f   centObs0\n", piaacmc[0].centObs0);
    fprintf(fp, "%10.6f   centObs1\n", piaacmc[0].centObs1);
    fprintf(fp, "%10ld   NBradpts\n", piaacmc[0].NBradpts);
    fprintf(fp, "%10.6f   r0lim\n", piaacmc[0].r0lim);
    fprintf(fp, "%10.6f   r1lim\n", piaacmc[0].r1lim);




    fprintf(fp, "%10ld    NBLyotStop\n", piaacmc[0].NBLyotStop);

    printf("%10ld    NBLyotStop\n", piaacmc[0].NBLyotStop);
    fflush(stdout);


    for(i=0; i<10; i++)
    {
        if(i<piaacmc[0].NBLyotStop)
        {
            sprintf(fname, "!%s/LyotStop%ld.fits", dname, i);
            if(piaacmc[0].IDLyotStop[i]!=-1)
                save_fits(data.image[piaacmc[0].IDLyotStop[i]].md[0].name, fname);
            fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
        }
        else
            fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
    }

    sprintf(fname, "!%s/piaa0Cmodes.fits", dname);
    if(piaacmc[0].piaa0CmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa0CmodesID].md[0].name, fname);

    sprintf(fname, "!%s/piaa0Fmodes.fits", dname);
    if(piaacmc[0].piaa0FmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa0FmodesID].md[0].name, fname);

    sprintf(fname, "!%s/piaa1Cmodes.fits", dname);
    if(piaacmc[0].piaa1CmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa1CmodesID].md[0].name, fname);

    sprintf(fname, "!%s/piaa1Fmodes.fits", dname);
    if(piaacmc[0].piaa1FmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa1FmodesID].md[0].name, fname);


    fprintf(fp, "%10.6f    fpmaskradld\n", piaacmc[0].fpmaskradld);
    fprintf(fp, "%10ld    focmNBzone\n", piaacmc[0].focmNBzone);
    fprintf(fp, "%10.6f   Fratio\n", piaacmc[0].Fratio);

    sprintf(fname, "!%s/fpm_zonez.fits", dname);
    if(piaacmc[0].zonezID!=-1)
        save_fits(data.image[piaacmc[0].zonezID].md[0].name, fname);

    fprintf(fp, "%10.6f    fpmaskamptransm\n", piaacmc[0].fpmaskamptransm);

    sprintf(fname, "!%s/fpm_zonea.fits", dname);
    if(piaacmc[0].zoneaID!=-1)
        save_fits(data.image[piaacmc[0].zoneaID].md[0].name, fname);



    fprintf(fp, "%10.6f   fpzfactor\n", piaacmc[0].fpzfactor);
    fprintf(fp, "%10.6g   fpmRad\n", piaacmc[0].fpmRad);
    fprintf(fp, "%10ld    NBrings\n", piaacmc[0].NBrings);
    fprintf(fp, "%10ld    fpmarraysize \n", piaacmc[0].fpmarraysize);
    fprintf(fp, "%10d    fpmmaterial\n", piaacmc[0].fpmmaterial);


    fclose(fp);

    return(0);
}



int PIAAsimul_loadpiaacmcconf(char *dname)
{
    char command[200];
    int r;
    FILE *fp;
    char fname[200];
    char imname[200];
    long i;

    int tmpi;
    long tmpl;
    float tmpf;
    double tmplf;



    sprintf(fname,"%s/piaacmcparams.conf", dname);

    fp = fopen(fname, "r");
    if(fp==NULL)
    {
        printf("Configuration file \"%s\" does not exist (yet), using previously set configuration\n", fname);
        fflush(stdout);
        r = 0;
    }
    else
    {
        r = fscanf(fp, "%f   beamrad\n", &tmpf);
        piaacmc[0].beamrad = tmpf;

        r = fscanf(fp, "%ld    size\n", &tmpl);
        piaacmc[0].size = tmpl;

        r = fscanf(fp, "%g   pixscale\n", &tmpf);
        piaacmc[0].pixscale = tmpf;

        r = fscanf(fp, "%f   piaasep\n", &tmpf);
        piaacmc[0].piaasep = tmpf;

        r = fscanf(fp, "%lf   centObs0\n", &tmplf);
        piaacmc[0].centObs0 = tmplf;

        r = fscanf(fp, "%lf   centObs1\n", &tmplf);
        piaacmc[0].centObs1 = tmplf;

        r = fscanf(fp, "%ld   NBradpts\n", &tmpl);
        piaacmc[0].NBradpts = tmpl;

        r = fscanf(fp, "%lf   r0lim\n", &tmplf);
        piaacmc[0].r0lim = tmplf;

        r = fscanf(fp, "%lf   r1lim\n", &tmplf);
        piaacmc[0].r1lim = tmplf;




        r = fscanf(fp, "%10ld    NBLyotStop\n", &piaacmc[0].NBLyotStop);
        for(i=0; i<10; i++)
        {
            if(i<piaacmc[0].NBLyotStop)
            {
                sprintf(fname, "%s/LyotStop%ld.fits", dname, i);
                sprintf(imname, "lyotstop%ld", i);
                piaacmc[0].IDLyotStop[i] = load_fits(fname, imname);
                r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                piaacmc[0].LyotStop_zpos[i] = tmplf;
            }
            else
            {
                r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                piaacmc[0].LyotStop_zpos[i] = tmplf;
                printf("LYOT STOP %ld POS : %lf\n", i, tmplf);
            }
        }


        sprintf(fname, "%s/piaa0Cmodes.fits", dname);
        piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff");
		if(piaacmc[0].piaa0CmodesID==-1)
		{
			sprintf(fname, "%s/piaaref/piaa0Cmodes.fits", dname);
			piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff");
		}			


        sprintf(fname, "%s/piaa0Fmodes.fits", dname);
        piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff");
		if(piaacmc[0].piaa0FmodesID==-1)
		{
			sprintf(fname, "%s/piaaref/piaa0Fmodes.fits", dname);
			piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff");
		}			

        sprintf(fname, "%s/piaa1Cmodes.fits", dname);
        piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff");
		if(piaacmc[0].piaa1CmodesID==-1)
		{
			sprintf(fname, "%s/piaaref/piaa1Cmodes.fits", dname);
			piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff");
		}			

        sprintf(fname, "%s/piaa1Fmodes.fits", dname);
        piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff");
		if(piaacmc[0].piaa1FmodesID==-1)
		{
			sprintf(fname, "%s/piaaref/piaa1Fmodes.fits", dname);
			piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff");
		}			




        r = fscanf(fp, "%f    fpmaskradld\n", &tmpf);
        piaacmc[0].fpmaskradld = tmpf;

        r = fscanf(fp, "%ld    focmNBzone\n",  &tmpl);
        piaacmc[0].focmNBzone = tmpl;

        r = fscanf(fp, "%f   Fratio\n",      &tmpf);
        piaacmc[0].Fratio = tmpf;

        sprintf(fname, "%s/fpm_zonez.fits", dname);
        delete_image_ID("fpmzt");
        piaacmc[0].zonezID = load_fits(fname, "fpmzt");

        r = fscanf(fp, "%f   fpmaskamptransm\n",    &tmpf);
        piaacmc[0].fpmaskamptransm = tmpf;

        sprintf(fname, "%s/fpm_zonea.fits", dname);
        delete_image_ID("fpmza");
        piaacmc[0].zoneaID = load_fits(fname, "fpmza");

        r = fscanf(fp, "%f   fpzfactor\n",   &tmpf);
        piaacmc[0].fpzfactor = tmpf;

        r = fscanf(fp, "%f   fpmRad\n",      &tmpf);
        piaacmc[0].fpmRad = tmpf;

        r = fscanf(fp, "%ld    NBrings\n",     &tmpl);
        piaacmc[0].NBrings = tmpl;

        r = fscanf(fp, "%ld    fpmarraysize \n", &tmpl);
        piaacmc[0].fpmarraysize = tmpl;

        r = fscanf(fp, "%d    fpmmaterial\n",  &tmpi);
        piaacmc[0].fpmmaterial = tmpi;

        r = 1;

        fclose(fp);
    }

    return(r);
}


/// Make Lyot stop geometry
/// param[in] IDincoh_name   Incoherent Lyot pupil intensity response to off-axis sources
/// parampin] IDmc_name      Intensity Lyot pupil image for on-axis source
//

long PIAACMCsimul_mkLyotMask(char *IDincoh_name, char *IDmc_name, char *IDzone_name, double throughput, char *IDout_name)
{
    long ID;
    long IDmc, IDincoh, IDzone;
    double val, val1, v, v0, bestval, v_best, rsl_best;
    double rsl, rsl0;
    long iter, NBiter;
    long ii;
    long xsize, ysize;
    long IDout;
    float sigma = 4.0;
    int filter_size = 10;



    NBiter = 100;

    sigma = 0.01*piaacmc[0].beamrad/piaacmc[0].pixscale;
    filter_size = (long) (sigma*2.0);

    printf("IDincoh_name : %s   %ld\n", IDincoh_name, image_ID(IDincoh_name));
    printf("IDmc_name    : %s   %ld\n", IDmc_name, image_ID(IDmc_name));
    printf("IDzone_name  : %s   %ld\n", IDzone_name, image_ID(IDzone_name));

    IDincoh = gauss_filter(IDincoh_name, "incohg", sigma, filter_size);
    //IDincoh = image_ID(IDincoh_name);
 
//    IDmc = image_ID(IDmc_name);
	IDmc = gauss_filter(IDmc_name, "mcg", sigma, filter_size);

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

    v0 = img_percentile("mcg", 0.99);
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
					else
					{
			//			printf("               %.12lf / %.12lf    %.12lf / %.12lf  -> %.12lf  %.12lf\n", rsl, rsl0, v, v0, val, val1);
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
    delete_image_ID("incohg");
    delete_image_ID("mcg");

    return(IDout);
}





//
/// Lyot stops positions from zmin to zmax relative to current, working back (light goes from 0 to zmax)
/// @param[in] FluxTOT  total flux in current plane
/// @param[in] FluxLim   max flux allowed from star
//
double PIAACMCsimul_optimizeLyotStop(char *IDamp_name, char *IDpha_name, char *IDincoh_name, float zmin, float zmax, double throughput, long NBz, long NBmasks)
{
    // initial guess places Lyot stops regularly from zmin to zmax
    // light propagates from zmin to zmax
    // we start with a single mask in zmax, and work back
    //

    double ratio;

    long ID, IDa, IDp;
    long nblambda; // number of wavelengths, read from input cube
    float *zarray;
    long l;
    double zprop;

    char nameamp[200];
    char namepha[200];
    char nameint[200];
    char fname[200];
    char fname1[200];
    long xsize, ysize;
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

    long IDincoh, IDint, IDmc;
    float rsl;
    long iter;
    long NBiter = 100;
    double val1;
    long IDm;
    char name[200];

	long IDre, IDim, IDreg, IDimg;
	double amp, pha, re, im;
    float sigma;
    int filter_size;

long IDlscumul;

    FILE *fp;
    double alpha = 1.05; // norm alpha used to identify best plane


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
    IDincoh = image_ID(IDincoh_name);

    xsize = data.image[IDa].md[0].size[0];
    ysize = data.image[IDa].md[0].size[1];

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
            for(m=0; m<NBmasks; m++)
                if((r>rinarray[m]-0.0001)&&(r<routarray[m]+0.0001))
                    data.image[IDzone].array.F[jj*xsize+ii] = m;
        }

    sprintf(fname, "!%s/LMzonemap.fits", piaacmcconfdir);
    save_fits("LMzonemap", fname);

    // initialize zarray
    for(l=0; l<NBz; l++)
        zarray[l] = zmin + (zmax-zmin)*l/(NBz-1);

    //  save_fits(nameint, fname);


	IDre = create_2Dimage_ID("retmpim", xsize, ysize);
	IDim = create_2Dimage_ID("imtmpim", xsize, ysize);

    ID = create_3Dimage_ID("LMintC", xsize, ysize, NBz);
    for(l=0; l<NBz; l++)
    {
        sprintf(nameamp, "LMPamp%02ld", l);
        sprintf(namepha, "LMPpha%02ld", l);
        zprop = zarray[l];
        OptSystProp_propagateCube(optsyst, 0, IDamp_name, IDpha_name, nameamp, namepha, zprop);

        // collapse broadband intensities
        //  sprintf(nameint, "LMPint%02ld", l);
        IDa = image_ID(nameamp);
        IDp = image_ID(namepha);
        //ID = create_2Dimage_ID(nameint, xsize, ysize);


        for(m=0; m<NBmasks; m++)
        {
            tot2array[l*NBmasks+m] = 0.0;
            totarray[l*NBmasks+m] = 0.0;
        }

		// convolve in complex amplitude
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
				data.image[ID].array.F[l*xsize*ysize+ii] += re*re+im*im;
			}
		}
		delete_image_ID("retmpimg");
		delete_image_ID("imtmpimg");
		 

        for(ii=0; ii<xsize*ysize; ii++)
        {
            m = (long) (data.image[IDzone].array.F[ii]+0.1);

  //          for(k=0; k<nblambda; k++)
    //        {
				// convolve in complex amplitude
				
				
//                data.image[ID].array.F[l*xsize*ysize+ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii];
                if((m>-1)&&(m<NBmasks))
                {
                    totarray[l*NBmasks+m] += data.image[ID].array.F[l*xsize*ysize+ii];
                    tot2array[l*NBmasks+m] += pow(data.image[ID].array.F[l*xsize*ysize+ii], alpha);
                }
      //      }
        }
        //sprintf(fname, "!LMPint%02ld.fits", l);
        //      save_fits(nameint, fname);

        delete_image_ID(nameamp);
        delete_image_ID(namepha);
    }
	delete_image_ID("retmpim");
	delete_image_ID("imtmpim");


    sprintf(fname,  "!%s/LMintC.fits", piaacmcconfdir);
   // save_fits("LMintC", fname);

    IDmc = create_2Dimage_ID("Lcomb", xsize, ysize);

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
            printf("MASK %ld   z= %f  ->  %g   ( %g %g) \n", m, zarray[l], val, tot2array[l*NBmasks+m], totarray[l*NBmasks+m]);
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
        //     sprintf(nameint, "LMPint%02ld", lbest);
        //      IDint = image_ID(nameint);

        for(ii=0; ii<xsize*ysize; ii++)
            if(m==data.image[IDzone].array.F[ii])
                data.image[IDmc].array.F[ii] = data.image[IDint].array.F[lbest*xsize*ysize+ii];
    }
    fclose(fp);

    sprintf(fname, "!%s/Lcomb.fits", piaacmcconfdir);
    save_fits("Lcomb", fname);

    ID = PIAACMCsimul_mkLyotMask(IDincoh_name, "Lcomb", "LMzonemap", throughput, "LMask");
    sprintf(fname, "!%s/LMask.fits", piaacmcconfdir);
    save_fits("LMask", fname);

    delete_image_ID("Lcomb");

	IDlscumul = create_2Dimage_ID("LMcumul", xsize, ysize);
      	for(ii=0;ii<xsize*ysize;ii++)
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


        	for(ii=0;ii<xsize*ysize;ii++)
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

    delete_image_ID("LMzonemap");

    return(ratio);
}


///
/// solves for focal plane mask solution using pre-computed zone responses
///
/// @param[in] fpmresp_array   Mask zones responses, float array
/// @param[in] zonez_array     zone thicknesses, double array
/// @param[in] dphadz_array    for each lambda, pha = thickness x dphadt_array[lambdaindex]
/// @param[out] outtmp_array    output temp array
///
/// written to be fast, no checking of array sizes
/// all arrays pre-allocated outside this function
///
double PIAACMCsimul_achromFPMsol_eval(float *fpmresp_array, double *zonez_array, double *dphadz_array, float *outtmp_array, long vsize, long nbz, long nbl)
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
			for(evalii=0; evalii<vsize; evalii++)
				outtmp_array[evalk*vsize+evalii] = fpmresp_array[evalk*(nbz+1)*vsize+evalii];


 
      			for(evalmz=0; evalmz<nbz; evalmz++)
			{
				evalpha = zonez_array[evalmz]*dphadz_array[evalk];
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

    evalval = 0.0;
    for(evalii=0; evalii<vsize*nbl; evalii++)
    {
        evalv1 = outtmp_array[evalii];
        evalval += evalv1*evalv1;
    }
    //  evalval /= vsize*nbl;

    return evalval;
}




double PIAACMCsimul_achromFPMsol_eval_zonezderivative(long zone, float *fpmresp_array, double *zonez_array, double *dphadz_array, float *outtmp_array, long vsize, long nbz, long nbl)
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

        evalpha = zonez_array[evalmz]*dphadz_array[evalk];
        evalcosp = -sin(evalpha)*dphadz_array[evalk]; //cos(evalpha);
        evalsinp = cos(evalpha)*dphadz_array[evalk]; //sin(evalpha);
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
    value /= CnormFactor*piaacmc[0].nblambda;

    if(LOOPCNT==0)
        cval0 = value;
    LOOPCNT++;

    return (value);

}




/**
 *
 * @brief Main simulation routine
 *
 * @param[in] confindex  PIAACMC configuration index pointing to the input/output directory number
 * @param[in] mode       Type of operation to be performed
 *
 * mode values:
 * - 0 : compute on-axis propagation for specified configuration. If configuration does not exist, create idealized monochromatic PIAACMC (intended to create a new index) and compute on-axis propagation
 * - 1 : optimize Lyot stop(s) locations
 * - 2 : optimize focal plane mask transmission for idealized monochromatic PIAACMC
 * - 3 : run without focal plane mask (for testing and calibration)
 * - 4 : linear optimization around current design, free parameters = PIAA optics shapes (runs for infinite number of iterations, track progress by looking at val.opt file)
* - 5 : optimize Lyot stops shapes and positions
* - 6 : test off-axis performance (to be written)
* - 10 : setup polychromatic optimization
* - 11 : Compute polychromatic response to zones, store result in FPMresp
* - 12 : search for best mask solution using FPMresp, random search
* - 40 : optimize PIAA optics shapes and focal plane mask transmission (idealized PIAACMC)
 */

int PIAACMCsimul_exec(long confindex, long mode)
{
    long NBparam;
    FILE *fp;
    char command[500];
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

    double valref, valbest, val0;
    double parambest[10000]; // for scanning
    double paramref[10000];
    long i, ii, jj, kk;
    long IDv, ID, IDref, IDa;
    long IDmodes;
    long xsize, ysize, zsize;
    long k;

    long iter;
    long NBiter = 1000;

    long IDfpmresp, IDref1;
    double t, a, dpha, amp;
    int zi;

    char fname[200];
    char fnamelog[200];
    long IDm, ID1D, ID1Dref;
    long size1Dvec;


    // OPTIMIZATION PARAMETERS
    int REGPIAASHAPES = 0;
    float piaa0C_regcoeff = 0.0e-7; // regularization coeff
    float piaa1C_regcoeff = 0.0e-7; // regularization coeff

    float piaa0C_regcoeff_alpha = 1.0; // regularization coeff power
    float piaa1C_regcoeff_alpha = 1.0; // regularization coeff power
    int r;

    double val, v1, pha, cosp, sinp, re, im, re1, im1;


    double fpmradld = 0.9;
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    double range, stepsize;
    int loopOK;
    long NBls, ls, ls1, NBpropstep;
    double lstransm;
    long mz;

    int LINOPT = 0; // 1 if perform linear optimization
    long ii1, ii2, ki, kv, ki1;

    double scangain;
    double scanstepgain = 0.01;
    int linscanOK;
    double valold;
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

    piaacmc = NULL;

    if(optsyst==NULL)
        optsyst = (OPTSYST*) malloc(sizeof(OPTSYST));
	for(elem=0;elem<100;elem++)
		optsyst[0].keepMem[i] = 0;



    sprintf(piaacmcconfdir, "piaacmcconf%03ld", confindex);
	sprintf(data.SAVEDIR, "%s", piaacmcconfdir);

    switch (mode) {

    case 0 :  // Run existing config for on-axis point source. If new, create centrally obscured idealized PIAACMC
        if((IDv=variable_ID("PIAACMC_centobs0"))!=-1)
            centobs0 = data.variable[IDv].value;
        if((IDv=variable_ID("PIAACMC_centobs1"))!=-1)
            centobs1 = data.variable[IDv].value;
        if((IDv=variable_ID("PIAACMC_fpmradld"))!=-1)
            fpmradld = data.variable[IDv].value;
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        PIAACMCsimul_makePIAAshapes();
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
        printf("valref = %g\n", valref);       
        break;

	case 100 : // evaluate current design: polychromatic contrast, pointing sensitivity
		PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        PIAACMCsimul_makePIAAshapes();
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
   
    /*    valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
		sprintf(fname,"!%s/psfi_00_00.fits", piaacmcconfdir);
        save_fits("psfi", fname);
        
        valref = PIAACMCsimul_computePSF(1.0, 0.0, 0, optsyst[0].NBelem);
 		sprintf(fname,"!%s/psfi_10_00.fits", piaacmcconfdir);
        save_fits("psfi", fname);
        valref = PIAACMCsimul_computePSF(2.0, 0.0, 0, optsyst[0].NBelem);
 		sprintf(fname,"!%s/psfi_20_00.fits", piaacmcconfdir);
        save_fits("psfi", fname);
        valref = PIAACMCsimul_computePSF(3.0, 0.0, 0, optsyst[0].NBelem);
 		sprintf(fname,"!%s/psfi_30_00.fits", piaacmcconfdir);
        save_fits("psfi", fname);
        valref = PIAACMCsimul_computePSF(4.0, 0.0, 0, optsyst[0].NBelem);
 		sprintf(fname,"!%s/psfi_40_00.fits", piaacmcconfdir);
        save_fits("psfi", fname);*/
        
      valref = PIAACMCsimul_computePSF(5.0, 0.0, 0, optsyst[0].NBelem);
 		sprintf(fname,"!%s/psfi_50_00.fits", piaacmcconfdir);
        save_fits("psfi", fname);
//		load_fits(fname, "psfi");

		
		ID = image_ID("psfi");
		xsize = data.image[ID].md[0].size[0];
		ysize = data.image[ID].md[0].size[1];
		zsize = data.image[ID].md[0].size[2];
		peakarray = (float*) malloc(sizeof(float)*zsize);
		for(kk=0;kk<zsize;kk++)
		{
			peakarray[kk] = 0.0;
			for(ii=0;ii<xsize*ysize;ii++)
			{
				val = data.image[ID].array.F[kk*xsize*ysize+ii];
				if(val>peakarray[kk])
					peakarray[kk] = val;
			}
		}
		avpeak = 0.0;
		for(kk=0;kk<zsize;kk++)
		{
			printf("peak %02ld  %10lf\n", kk, peakarray[kk]);
			avpeak += peakarray[kk];
		}
		avpeak /= zsize;
		free(peakarray);
		delete_image_ID("psfi");
 
        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
		sprintf(fname,"!%s/psfi_00_00.fits", piaacmcconfdir);
        save_fits("psfi", fname);
		//load_fits(fname, "psfi");

		ID = image_ID("psfi");
		

		/// compute contrast curve
		focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;
		printf("focscale = %f\n", focscale);
		eval_sepNBpt = (long) (eval_sepmaxld/eval_sepstepld);
		eval_contrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
		eval_contrastCurve_cnt = (double*) malloc(sizeof(double)*eval_sepNBpt);
		for(ri=0;ri<eval_sepNBpt;ri++)
		{
			eval_contrastCurve[ri] = 0.0;
			eval_contrastCurve_cnt[ri] = 0.0;
		}


		for(kk=0;kk<zsize;kk++)
			for(ii=0;ii<xsize;ii++)
				for(jj=0;jj<ysize;jj++)
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
				}
		sprintf(fname, "%s/ContrastCurve.txt", piaacmcconfdir);
		fp = fopen(fname, "w");
		for(ri=0;ri<eval_sepNBpt;ri++)
		{
			eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
			fprintf(fp, "%10f %10g %10g\n", eval_sepstepld*ri, eval_contrastCurve[ri], eval_contrastCurve_cnt[ri]);
		}
		fclose(fp);
		free(eval_contrastCurve);
		free(eval_contrastCurve_cnt);


		// pointing sensitivity
		IDps = create_3Dimage_ID("starim", piaacmc[0].size, piaacmc[0].size, zsize);

		valref = PIAACMCsimul_computePSF(0.01, 0.0, 0, optsyst[0].NBelem);
		sprintf(fname,"!%s/psfi_p0.fits", piaacmcconfdir);
        save_fits("psfi", fname);
		ID = image_ID("psfi");
		for(ii=0;ii<xsize*ysize*zsize;ii++)
			data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];

		valref = PIAACMCsimul_computePSF(-0.01, 0.0, 0, optsyst[0].NBelem);
		sprintf(fname,"!%s/psfi_m0.fits", piaacmcconfdir);
        save_fits("psfi", fname);
		ID = image_ID("psfi");
		for(ii=0;ii<xsize*ysize*zsize;ii++)
			data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];

		valref = PIAACMCsimul_computePSF(0.0, 0.01, 0, optsyst[0].NBelem);
		sprintf(fname,"!%s/psfi_0p.fits", piaacmcconfdir);
        save_fits("psfi", fname);
		ID = image_ID("psfi");
		for(ii=0;ii<xsize*ysize*zsize;ii++)
			data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];

		valref = PIAACMCsimul_computePSF(0.00, -0.01, 0, optsyst[0].NBelem);
		sprintf(fname,"!%s/psfi_0m.fits", piaacmcconfdir);
        save_fits("psfi", fname);
		ID = image_ID("psfi");
		for(ii=0;ii<xsize*ysize*zsize;ii++)
			data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];

		for(ii=0;ii<xsize*ysize*zsize;ii++)
			data.image[IDps].array.F[ii] /= 4.0;

		sprintf(fname,"!%s/psfi_starim.fits", piaacmcconfdir);
        save_fits("starim", fname);

		/// compute contrast curve
		focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;
		printf("focscale = %f\n", focscale);
		eval_sepNBpt = (long) (eval_sepmaxld/eval_sepstepld);
		eval_contrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
		eval_contrastCurve_cnt = (double*) malloc(sizeof(double)*eval_sepNBpt);
		for(ri=0;ri<eval_sepNBpt;ri++)
		{
			eval_contrastCurve[ri] = 0.0;
			eval_contrastCurve_cnt[ri] = 0.0;
		}


		for(kk=0;kk<zsize;kk++)
			for(ii=0;ii<xsize;ii++)
				for(jj=0;jj<ysize;jj++)
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
				}
		sprintf(fname, "%s/ContrastCurve_ps.txt", piaacmcconfdir);
		fp = fopen(fname, "w");
		for(ri=0;ri<eval_sepNBpt;ri++)
		{
			eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
			fprintf(fp, "%10f %10g %10g\n", eval_sepstepld*ri, eval_contrastCurve[ri], eval_contrastCurve_cnt[ri]);
		}
		fclose(fp);
		free(eval_contrastCurve);
		free(eval_contrastCurve_cnt);		
	break;


    case 1 : // optimize Lyot stop positions
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        PIAACMCsimul_makePIAAshapes();
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        // initialization
        if((IDv=variable_ID("PIAACMC_lsoptrange"))!=-1)
            range = data.variable[IDv].value;
        else
            range = 2.0;
        stepsize = range/3.0;
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            paramref[ls] = piaacmc[0].LyotStop_zpos[ls];
        NBiter = 4;


        sprintf(fnamelog, "%s/result_LMpos.log", piaacmcconfdir);
		fp = fopen(fnamelog, "w");
		fclose(fp);




		stepsize = range/5.0;
		for(iter=0; iter<NBiter; iter++)
		{
			for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
			{
				piaacmc[0].LyotStop_zpos[ls] = paramref[ls]-range;
				parambest[ls] = piaacmc[0].LyotStop_zpos[ls];

				loopOK = 1;
				valbest = 1.0;

				while(piaacmc[0].LyotStop_zpos[ls]<paramref[ls]+range)
				{
					if(invPIAAmode == 2) // invPIAA -> Lyot masks
					{
						optsyst[0].keepMem[8] = 1;
						val = PIAACMCsimul_computePSF(0.0, 0.0, 8, optsyst[0].NBelem);
					}
					else // Lyot masks -> invPIAA
					{
						optsyst[0].keepMem[6] = 1;
						val = PIAACMCsimul_computePSF(0.0, 0.0, 6, optsyst[0].NBelem);
					}

					if(val<valbest)
					{
						parambest[ls] = piaacmc[0].LyotStop_zpos[ls];
						valbest = val;
					}

					fp = fopen(fnamelog, "a");
					for(ls1=0; ls1<piaacmc[0].NBLyotStop; ls1++)
						fprintf(fp," %lf", piaacmc[0].LyotStop_zpos[ls1]);
					fprintf(fp, " %g\n", val);
					fclose(fp);

					piaacmc[0].LyotStop_zpos[ls] += stepsize;
				}
				printf("BEST SOLUTION :  ");
				paramref[ls] = parambest[ls];
				piaacmc[0].LyotStop_zpos[ls] = paramref[ls];
				printf(" %lf", parambest[ls]);
				printf(" %g\n", valbest);
			}

			fp = fopen(fnamelog, "a");
			fprintf(fp, "\n");
			fclose(fp);

			range *= 0.3;
			stepsize = range/3.0;
		}
		for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            piaacmc[0].LyotStop_zpos[ls] = parambest[ls];
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);
        break;



    case 2 : // optimize focal plane mask transmission for monochromatic idealized PIAACMC
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        PIAACMCsimul_makePIAAshapes();
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        // initialization
        range = 0.2;
        stepsize = range/3.0;
        paramref[0] = piaacmc[0].fpmaskamptransm;
        NBiter = 5;

        sprintf(fnamelog, "%s/result_fpmt.log", piaacmcconfdir);
        fp = fopen(fnamelog, "w");
        fclose(fp);

        for(iter=0; iter<NBiter; iter++)
        {
            piaacmc[0].fpmaskamptransm = paramref[0]-range;
            parambest[0] = piaacmc[0].fpmaskamptransm;

            loopOK = 1;
            valbest = 1.0;

            while(loopOK==1)
            {
                FORCE_CREATE_fpmza = 1;
                PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
                val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);

                if(val<valbest)
                {
                    parambest[0] = piaacmc[0].fpmaskamptransm;
                    valbest = val;
                }

                fp = fopen(fnamelog, "a");
                fprintf(fp," %lf", piaacmc[0].fpmaskamptransm);
                fprintf(fp, " %g  %ld %g %g\n", val, iter, range, stepsize);
                fclose(fp);

                ls = 0;
                piaacmc[0].fpmaskamptransm += stepsize;
                if(piaacmc[0].fpmaskamptransm>paramref[0]+range+0.001*stepsize)
                    loopOK = 0;
            }

            printf("BEST SOLUTION :  ");

            paramref[0] = parambest[0];
            printf(" %lf", parambest[0]);

            printf(" %g\n", valbest);


            fp = fopen(fnamelog, "a");
            fprintf(fp, "\n");
            fclose(fp);

            range *= 0.3;
            stepsize = range/3.0;
        }

        piaacmc[0].fpmaskamptransm = parambest[0];
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);
        FORCE_CREATE_fpmza = 0;
        break;


    case 3 : // calibrate
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        PIAACMCsimul_makePIAAshapes();
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        paramref[0] = piaacmc[0].fpmaskamptransm;

        piaacmc[0].fpmaskamptransm = -1.0;
        FORCE_CREATE_fpmza = 1;
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);

        // restore original configuration
        piaacmc[0].fpmaskamptransm = paramref[0];
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);
        FORCE_CREATE_fpmza = 0;

        break;

    case 4 : // optimize PIAA optics shapes
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        LINOPT = 1; // perform linear optimization
        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            NBiter = (long) data.variable[IDv].value+0.01;
        else
            NBiter = 1000;

        kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
            kmax = (long) data.variable[IDv].value+0.01;

        if(kmax>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
            kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];

        NBparam = 0;
        for(k=0; k<kmax; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }

        for(k=0; k<kmax; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }
        FORCE_MAKE_PIAA0shape = 1;
        FORCE_MAKE_PIAA1shape = 1;
        break;
        

        case 5 : // optimize Lyot stops shapes and positions
		PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
		PIAACMCsimul_makePIAAshapes();
		optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

		NBls = 3;
		if((IDv=variable_ID("PIAACMC_nblstop"))!=-1)
			NBls = (long) data.variable[IDv].value+0.01;

		NBpropstep = 100;
		if((IDv=variable_ID("PIAACMC_nbpropstep"))!=-1)
			NBpropstep = (long) data.variable[IDv].value+0.01;
		
		lstransm = 0.85;
		if((IDv=variable_ID("PIAACMC_lstransm"))!=-1)
			lstransm = (double) data.variable[IDv].value;
		printf("lstransm  = %f\n", lstransm);

		if(invPIAAmode==2)
			optsyst[0].keepMem[7] = 1;
		else
			optsyst[0].keepMem[5] = 1;
		PIAACMCsimul_computePSF(3.0, 0.0, 0, optsyst[0].NBelem);
        
		if(invPIAAmode==2)
            IDa = image_ID("WFamp_007");
        else
            IDa = image_ID("WFamp_005");

        xsize = data.image[IDa].md[0].size[0];
        ysize = data.image[IDa].md[0].size[1];

        ID = create_2Dimage_ID("OAincoh", xsize, ysize);
        for(ii=0; ii<xsize*ysize; ii++)
            for(k=0; k<optsyst[0].nblambda; k++)
                data.image[ID].array.F[ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii]/4;

        PIAACMCsimul_computePSF(-3.0, 0.0, 0, optsyst[0].NBelem);
        if(invPIAAmode==2)
            IDa = image_ID("WFamp_007");
        else
            IDa = image_ID("WFamp_005");
        for(ii=0; ii<xsize*ysize; ii++)
            for(k=0; k<optsyst[0].nblambda; k++)
                data.image[ID].array.F[ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii]/4;

        PIAACMCsimul_computePSF(0.0, 3.0, 0, optsyst[0].NBelem);
        if(invPIAAmode==2)
            IDa = image_ID("WFamp_007");
        else
            IDa = image_ID("WFamp_005");
        for(ii=0; ii<xsize*ysize; ii++)
            for(k=0; k<optsyst[0].nblambda; k++)
                data.image[ID].array.F[ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii]/4;

        PIAACMCsimul_computePSF(0.0, -3.0, 0, optsyst[0].NBelem);
        if(invPIAAmode==2)
            IDa = image_ID("WFamp_007");
        else
            IDa = image_ID("WFamp_005");
        for(ii=0; ii<xsize*ysize; ii++)
            for(k=0; k<optsyst[0].nblambda; k++)
                data.image[ID].array.F[ii] += data.image[IDa].array.F[k*xsize*ysize+ii]*data.image[IDa].array.F[k*xsize*ysize+ii]/4;

        sprintf(fname, "!%s/OAincoh.fits", piaacmcconfdir);
        save_fits("OAincoh", fname);

        PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
        if(invPIAAmode==2)
            PIAACMCsimul_optimizeLyotStop("WFamp_007", "WFpha_007", "OAincoh", -3.0, 0.0, lstransm, NBpropstep, NBls);
        else
            PIAACMCsimul_optimizeLyotStop("WFamp_005", "WFpha_005", "OAincoh", -3.0, 0.0, lstransm, NBpropstep, NBls);
        delete_image_ID("OAincoh");

        piaacmc[0].NBLyotStop = NBls;

        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            piaacmc[0].LyotStop_zpos[ls] += optsyst[0].elemZpos[5];

        PIAAsimul_savepiaacmcconf(piaacmcconfdir);
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
        {
            sprintf(command, "cp ./%s/optLM%02ld.fits ./%s/LyotStop%ld.fits", piaacmcconfdir, ls, piaacmcconfdir, ls);
            r = system(command);
        }
        break;


    case 6: // test off-axis performance
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        PIAACMCsimul_makePIAAshapes();
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
        PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
        break;


    case 10 : // setup multizone ring mask
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
        piaacmc[0].fpmaskamptransm = 1.0;
		piaacmc[0].NBrings = 16;
        piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio * 2.0; // 2 l/D radius at central lambda
        FORCE_CREATE_fpmzmap = 1;
        FORCE_CREATE_fpmza = 1;
        FORCE_CREATE_fpmzt = 1;
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0);
        PIAAsimul_savepiaacmcconf(piaacmcconfdir);
        break;

    case 11 : // Compute polychromatic response to zones, store result in FPMresp
		FORCE_CREATE_fpmzmap = 1;
		FORCE_CREATE_fpmzt = 1;
		FORCE_CREATE_fpmza = 1;
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        PIAACMCsimul_makePIAAshapes();		
		focmMode = data.image[piaacmc[0].zonezID].md[0].size[0]+10;  // response for no focal plane mask
        optsyst[0].FOCMASKarray[0].mode = 1; // 1-fpm
        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
        printf("val = %g\n", val);
        ID = image_ID("imvect");
        
		// WARNING: FPMresp size[1] is nbzones+1, as fist vector stored is the response for light outside the mask


        // axis 0: eval pts (ii) - size = data.image[ID].md[0].size[0]
        // axis 1: zones (mz) - size = data.image[piaacmc[0].zonezID].md[0].size[0]+1
        // axis 3: lambda (k) - size = piaacmc[0].nblambda
        //
        // indexing :  k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + mz*data.image[ID].md[0].size[0] + ii
        IDfpmresp = image_ID("FPMresp");
        if(IDfpmresp==-1)
            IDfpmresp = create_3Dimage_ID("FPMresp", data.image[ID].md[0].size[0], piaacmc[0].focmNBzone+1, piaacmc[0].nblambda);

		// light outside mask
        for(k=0; k<piaacmc[0].nblambda; k++)
            for(ii=0; ii<data.image[ID].md[0].size[0]; ii++)
                data.image[IDfpmresp].array.F[k*(piaacmc[0].focmNBzone+1)*data.image[ID].md[0].size[0] + ii] = data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];

        sprintf(fname, "!%s/FPMresp%d_%02ld_%02d.fits", piaacmcconfdir, computePSF_ResolvedTarget, piaacmc[0].focmNBzone, piaacmc[0].nblambda);
        save_fits("FPMresp", fname);

        for(mz=1; mz<piaacmc[0].focmNBzone+1; mz++)
        {
            focmMode = mz;
            optsyst[0].FOCMASKarray[0].mode = 0; // direct focal plane mask response
            optsyst[0].keepMem[4] = 1;
            val = PIAACMCsimul_computePSF(0.0, 0.0, 4, optsyst[0].NBelem);

            for(k=0; k<piaacmc[0].nblambda; k++)
                for(ii=0; ii<data.image[ID].md[0].size[0]; ii++)
                {
                    data.image[IDfpmresp].array.F[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + mz*data.image[ID].md[0].size[0] + ii] = data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];
                    data.image[IDfpmresp].array.F[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + ii] -= data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];
                }

            sprintf(fname, "!%s/FPMresp%d_%02ld_%02d.fits", piaacmcconfdir, computePSF_ResolvedTarget, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
            save_fits("FPMresp", fname);

            sprintf(fname, "!%s/psfi_z%02ld.fits", piaacmcconfdir, mz);
            save_fits("psfi", fname);
        }

        /*	 sprintf(fname, "!piaacmcfpm_z%02ld.fits", mz);
        printf("SAVING NOW --------------- %ld ---------------->  %s\n", image_ID("piaacmcfpm"), fname);
        mk_amph_from_complex("piaacmcfpm", "fpma", "fpmp");
        save_fits("fpma", fname);
        delete_image_ID("fpma");
        delete_image_ID("fpmp");
        */


        focmMode = -1;
        break;


    case 12 : // search for best mask solution using FPMresp
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
		PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);
        PIAACMCsimul_makePIAAshapes();
        optsyst[0].FOCMASKarray[0].mode = 1; // 1-fpm

		sprintf(fname,"%s/flux.txt", piaacmcconfdir);
		fp = fopen(fname, "r");
        for(elem=0; elem<optsyst[0].NBelem; elem++)
			{
				ret = fscanf(fp, "%lf %lf\n", &tmplf1, &tmplf2);
				optsyst[0].flux[elem] = tmplf1;
			}
		fclose(fp);

		computePSF_FAST_FPMresp = 1;
        sprintf(fname, "%s/CnormFactor.txt", piaacmcconfdir);
        fp = fopen(fname, "r");
        ret = fscanf(fp, "%lf", &CnormFactor);
        fclose(fp);

  /*      val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
        valbest = val;
        val0 = val;

        sprintf(fname, "!%s/imvect.fits", piaacmcconfdir);
        save_fits("imvect", fname);


      mk_amph_from_complex("piaacmcfpm", "fpma", "fpmp");
        sprintf(fname, "!%s/piaacmcfpma.fits", piaacmcconfdir);
        save_fits("fpma", fname);
        sprintf(fname, "!%s/piaacmcfpmp.fits", piaacmcconfdir);
        save_fits("fpmp", fname);
        delete_image_ID("fpma");
        delete_image_ID("fpmp");


        ID = image_ID("imvect");
*/


        sprintf(fname, "%s/FPMresp%d_%02ld_%02d.fits", piaacmcconfdir, computePSF_ResolvedTarget, piaacmc[0].focmNBzone, piaacmc[0].nblambda);
        IDfpmresp = load_fits(fname, "FPMresp");

        vsize = data.image[IDfpmresp].md[0].size[0]; // number of eval pts x2
        ID = create_2Dimage_ID("imvect1", vsize, piaacmc[0].nblambda);
        // measured speed:
        //
        // [nblambda=5]
        // 4.35 kHz on single thread (without omp)
        // 9.09 kHz with omp, 8 threads
        // -> better to launch multiple instances
        //
        // [nblambda=8]
        // 2.78 kHz
        // 13.89 kHz with omp, 8 threads
        // -> x5 speedup

        // allocate arrays for fast routine


        fpmresp_array = data.image[IDfpmresp].array.F;
        zonez_array = data.image[piaacmc[0].zonezID].array.D;
        zonez0_array = (double*) malloc(sizeof(double)*data.image[piaacmc[0].zonezID].md[0].size[0]); // reference point
        zonez1_array = (double*) malloc(sizeof(double)*data.image[piaacmc[0].zonezID].md[0].size[0]); // reference point
        zonezbest_array = (double*) malloc(sizeof(double)*data.image[piaacmc[0].zonezID].md[0].size[0]); // best point

        dphadz_array = (double*) malloc(sizeof(double)*piaacmc[0].nblambda);
        for(k=0; k<piaacmc[0].nblambda; k++)
            dphadz_array[k] = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, 1.0, optsyst[0].lambdaarray[k]);
        outtmp_array = (float*) malloc(sizeof(float)*vsize*piaacmc[0].nblambda);

        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
 
        printf("Preparing optimization ... \n");
        fflush(stdout);
        NBoptVar = data.image[piaacmc[0].zonezID].md[0].size[0];
        x = gsl_vector_alloc (NBoptVar);
        ss = gsl_vector_alloc (NBoptVar);
        fpmeval_func.n = NBoptVar;  /* number of function components */
        fpmeval_func.f = &f_evalmask;
        fpmeval_func.params = (void *) NULL;
        s = gsl_multimin_fminimizer_alloc (T, NBoptVar);

        // we assume the starting point is the best point
        for(mz=0; mz<NBoptVar; mz++)
            zonezbest_array[mz] = zonez_array[mz];

		if((IDv=variable_ID("PIAACMC_nbiterSA"))!=-1)
			NBITER_SA = (long) data.variable[IDv].value+0.01;


        ITERMAX = NBITER_SA; // SIMULATED ANNEALING - STARTING FROM ZERO
        printf("----------- STARTING SIMULATED ANNEALING ------------- [%ld]\n",ITERMAX);
        SAcoeff = 1.0e-6;
        //  SAcoeff = 1.0e-2*pow(cos(0.00005*iter1),8.0);
        for(iter1=0; iter1<ITERMAX; iter1++)
        {
            if((iter1==0)||(ran1()<0.0001)) // Start at best point and go back to it every once in a while
            {
                printf("[%05ld] Starting at best point ...\n", iter1);
                for(mz=0; mz<NBoptVar; mz++)
                {
                    zonez0_array[mz] = zonezbest_array[mz];
                    zonez_array[mz] = zonez0_array[mz];
                }
                val0 = PIAACMCsimul_achromFPMsol_eval(fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
                val0 /= CnormFactor*piaacmc[0].nblambda;
                if(iter1==0)
                    bestvalue = val0;
                printf("%10ld  val0 = %g  (%g)\n", iter1, val0, bestvalue);
            }

            amp = ran1()*pow(ran1(),3.0)*THICKRANGE;
            for(mz=0; mz<NBoptVar; mz++)
            {
                amp1 = pow(ran1(), 4.0);
                zonez1_array[mz] = zonez0_array[mz] + (1.0-2.0*ran1())*amp*amp1;
                zonez_array[mz] = zonez1_array[mz];
            }
            val1 = PIAACMCsimul_achromFPMsol_eval(fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
            val1 /= CnormFactor*piaacmc[0].nblambda;
            //			printf("%10ld  val1 = %g   (%g)\n", iter1, val1, bestvalue);


            if(val1<val0) // new point is best
                SA_MV = 1;
            else
            {
                tmp1 = exp(-(val1-val0)/SAcoeff); // between 0 and 1, close to 1 if val1~val0
                if(tmp1>ran1()) // if tmp1 close to 1, move
                {
                    // printf("[M] ");
                    SA_MV = 1;
                    SAcoeff *= 0.95; // decrease temperature
                }
                else
                {
                    //      printf("[S] ");
                    SA_MV = 0;
                    SAcoeff *= 1.01;
                }
            }

            if(SA_MV==1)
            {
                for(mz=0; mz<NBoptVar; mz++)
                    zonez0_array[mz] = zonez1_array[mz];
                val0 = val1;
            }

            if(val1<bestvalue)
            {
                bestvalue = val1;
                printf("[%05ld] -> %8g ---------------------------------------\n", iter1, bestvalue);
                sprintf(fname,"%s/bestmask.txt", piaacmcconfdir);
                fpbest = fopen(fname,"w");
                fprintf(fpbest,"%.20g", bestvalue);
                for(k=0; k<NBoptVar; k++)
                {
                    zonezbest_array[k] = zonez_array[k];
                    fprintf(fpbest," %g", zonezbest_array[k]);
                }
                fprintf(fpbest,"\n");
                fclose(fpbest);
                sprintf(fname, "!%s/fpm_zonez.fits", piaacmcconfdir);
                save_fits(data.image[piaacmc[0].zonezID].md[0].name, fname);
            }
        }





		if((IDv=variable_ID("PIAACMC_nbiterDS"))!=-1)
			NBITER_DS = (long) data.variable[IDv].value+0.01;


        ITERMAX = NBITER_DS;
        iterMax = 50000;
        printf("----------- STARTING DOWNHILL SIMPLEX ------------- [%ld]\n",ITERMAX);
		sprintf(fname, "%s/maskres,txt", piaacmcconfdir);
        fp = fopen(fname, "w");
        OK = 0; // did it improve ?
        KEEPlimit = bestvalue;
        KEEPlimit1 = 3.0*bestvalue;
        KEEP = 0;
        eps1 = 0.99999999;
        for(iter1=0; iter1<ITERMAX; iter1++) // DOWNHILL SIMPLEX METHOD
        {
            printf("%05ld ", iter1);

            printf("[KL %e %e]  ", KEEPlimit, KEEPlimit1);


            if((iter1==0)||(OK==1)) // set starting point = best point if 1st iteration or loop is still making progress
            {
                printf(" NEW ");
                for(mz=0; mz<NBoptVar; mz++)
                    gsl_vector_set (x, mz, zonezbest_array[mz]);
            }
            else
            {
                if(KEEP==0)
                {
                    printf("     ");
                    for(mz=0; mz<NBoptVar; mz++)
                        gsl_vector_set (x, mz, pow(ran1(),4.0)*(1.0-2.0*ran1())*THICKRANGE);
                }
                else
                {
                    printf("KEEP ");
                    KEEPlimit = (0.95*KEEPlimit + 0.05*bestvalue)*0.7 + 0.3*s->fval;
                    for(mz=0; mz<NBoptVar; mz++)
                        gsl_vector_set (x, mz, zonez_array[mz]);
                }
            }
            //      printf("%ld INIT: %e  ", iter1, f_evalmask (x, (void*) NULL));
            //  exit(0);

            /* Set initial step sizes to 1e-8 */
            gsl_vector_set_all (ss, 1.0e-9);

            /* Initialize method and iterate */
            iter = 0;

            gsl_multimin_fminimizer_set (s, &fpmeval_func, x, ss);
            LOOPCNT = 0;
            OK = 0;
            do
            {
                iter++;
                status = gsl_multimin_fminimizer_iterate(s);
                if (status)
                    break;

                mmsize = gsl_multimin_fminimizer_size (s);
                status = gsl_multimin_test_size (mmsize, 1e-11);


                //	  printf ("............[%05ld] %e ->  %e  %e  (%e)\n", iter, cval0, s->fval, size, bestvalue);



                if ((status == GSL_SUCCESS)||(s->fval < bestvalue * eps1))
                {
                    if(status == GSL_SUCCESS)
                    {
                        printf (" %e ->  %e [%5ld]  (%e)", cval0, s->fval, iter, bestvalue);
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
                    fprintf(fp, "%e %ld ", s->fval, iter);
                    for(mz=0; mz<NBoptVar; mz++)
                        fprintf(fp, "%e ", zonez_array[mz]);
                    fprintf(fp, "\n");

                    if(s->fval < bestvalue * eps1) //if(s->f<bestvalue)
                    {
                        OK = 1;
                        sprintf(fname,"%s/bestmask_DS.txt.CONF", piaacmcconfdir);
                        fpbest = fopen(fname, "w");
                        bestvalue = s->fval;
                        fprintf(fpbest,"%.20g", bestvalue);
                        for(mz=0; mz<NBoptVar; mz++)
                        {
                            zonezbest_array[mz] = zonez_array[mz];
                            fprintf(fpbest," %g", zonezbest_array[mz]);
                        }
                        fprintf(fpbest, "\n");

                        /*		  printf("%.20g", bestvalue);
                        for(k=0;k<NBzones;k++)
                          {
                            Zthickbest[k] = Zthick[k];
                            printf(" %g", Zthickbest[k]);
                          }
                          printf("\n");*/
                        fclose(fpbest);
						sprintf(fname, "!%s/fpm_zonez.fits", piaacmcconfdir);
						save_fits(data.image[piaacmc[0].zonezID].md[0].name, fname);
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
                    }
                }
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
        fclose(fp);
	
		

        exit(0);

        val = 1.0;
        for(i=0; i<10000000; i++)
        {
            for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                data.image[piaacmc[0].zonezID].array.D[k] = 1.0e-6*(1.0-2.0*ran1());


            val = PIAACMCsimul_achromFPMsol_eval(fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
            val /= CnormFactor*piaacmc[0].nblambda;

            if(val<valbest)
            {
                printf("%10ld  best value = %20g  (%20g)\n", i, val, val0);
                valbest = val;
                sprintf(fname, "!%s/fpm_zonez.fits", piaacmcconfdir);
                save_fits(data.image[piaacmc[0].zonezID].md[0].name, fname);
            }
        }

        sprintf(fname, "!%s/imvect1.fits", piaacmcconfdir);
        save_fits("imvect1", fname);
        break;




    case 13 : // optimize focal plane mask zones
		PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);
        
		sprintf(fname,"%s/flux.txt", piaacmcconfdir);
		fp = fopen(fname, "r");
        for(elem=0; elem<optsyst[0].NBelem; elem++)
			{
				ret = fscanf(fp, "%lf %lf\n", &tmplf1, &tmplf2);
				optsyst[0].flux[elem] = tmplf1;
			}
		fclose(fp);

        LINOPT = 1; // perform linear optimization
        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            NBiter = (long) data.variable[IDv].value+0.01;
        else
            NBiter = 20;


        sprintf(fname, "%s/FPMresp%d_%02ld_%02d.fits", piaacmcconfdir, computePSF_ResolvedTarget, piaacmc[0].focmNBzone, piaacmc[0].nblambda);
        IDfpmresp = load_fits(fname, "FPMresp");

        vsize = data.image[IDfpmresp].md[0].size[0]; // number of eval pts x2
        ID = create_2Dimage_ID("imvect1", vsize, piaacmc[0].nblambda);
    
        // allocate arrays for fast routine
        fpmresp_array = data.image[IDfpmresp].array.F;
        zonez_array = data.image[piaacmc[0].zonezID].array.D;
        dphadz_array = (double*) malloc(sizeof(double)*piaacmc[0].nblambda);
        for(k=0; k<piaacmc[0].nblambda; k++)
            {
				dphadz_array[k] = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial, 1.0, optsyst[0].lambdaarray[k]);
				printf("%ld  %g %g\n", k, optsyst[0].lambdaarray[k], dphadz_array[k]);
			}
        outtmp_array = (float*) malloc(sizeof(float)*vsize*piaacmc[0].nblambda);

		computePSF_FAST_FPMresp = 1;
        sprintf(fname, "%s/CnormFactor.txt", piaacmcconfdir);
        fp = fopen(fname, "r");
        ret = fscanf(fp, "%lf", &CnormFactor);
        fclose(fp);
        
		 for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
			data.image[piaacmc[0].zonezID].array.D[k] = MODampl*(1.0-2.0*ran1());
				
		NBparam = 0;
		for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++)
        {
            paramtype[NBparam] = DOUBLE;
            paramval[NBparam] = &data.image[piaacmc[0].zonezID].array.D[mz];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 2.0e-7;
            parammin[NBparam] = -1.0e-6;
            parammax[NBparam] = 1.0e-6;
            NBparam++;
        }
        break;
  



    case 40 : // optimize PIAA optics shapes and focal plane mask transmission (idealized PIAACMC)
//		FORCE_CREATE_fpmza = 1;

		if((IDv=variable_ID("PIAACMC_resolved"))!=-1)
            computePSF_ResolvedTarget = 1;

        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        LINOPT = 1; // perform linear optimization
        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            NBiter = (long) data.variable[IDv].value+0.01;
        else
            NBiter = 1000;

        kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
            kmaxC = (long) data.variable[IDv].value+0.01;

        if(kmaxC>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
            kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];


        kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptFterm"))!=-1)
            kmaxF = (long) data.variable[IDv].value+0.01;

        if(kmaxF>data.image[piaacmc[0].piaa0FmodesID].md[0].size[0])
            kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];



        NBparam = 0;

        paramtype[NBparam] = DOUBLE;
        paramval[NBparam] = &data.image[piaacmc[0].zoneaID].array.D[0];
//        piaacmc[0].fpmaskamptransm;
        paramdelta[NBparam] = 1.0e-2;
        parammaxstep[NBparam] = 1.0e-1;
        parammin[NBparam] = -1.0;
        parammax[NBparam] = 1.0;
        NBparam++;


        for(k=0; k<kmaxC; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }

        for(k=0; k<kmaxC; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }

        for(k=0; k<kmaxF; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa0FmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }

        for(k=0; k<kmaxF; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa1FmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }

        FORCE_MAKE_PIAA0shape = 1;
        FORCE_MAKE_PIAA1shape = 1;

        break;



 


   case 41 : // optimize PIAA optics shapes and focal plane mask zones (polychromatic)
        PIAAsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 1);
        LINOPT = 1; // perform linear optimization
        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            NBiter = (long) data.variable[IDv].value+0.01;
        else
            NBiter = 1000;

        kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
            kmaxC = (long) data.variable[IDv].value+0.01;

        if(kmaxC>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
            kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];


        kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptFterm"))!=-1)
            kmaxF = (long) data.variable[IDv].value+0.01;

        if(kmaxF>data.image[piaacmc[0].piaa0FmodesID].md[0].size[0])
            kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];



        NBparam = 0;

        paramtype[NBparam] = DOUBLE;
        paramval[NBparam] = &piaacmc[0].fpmaskamptransm;
        paramdelta[NBparam] = 1.0e-6;
        parammaxstep[NBparam] = 1.0e-3;
        parammin[NBparam] = -1.0;
        parammax[NBparam] = 1.0;
        NBparam++;


        for(k=0; k<kmaxC; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }

        for(k=0; k<kmaxC; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }

        for(k=0; k<kmaxF; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa0FmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }

        for(k=0; k<kmaxF; k++)
        {
            paramtype[NBparam] = FLOAT;
            paramvalf[NBparam] = &data.image[piaacmc[0].piaa1FmodesID].array.F[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-8;
            parammin[NBparam] = -1.0e-7;
            parammax[NBparam] = 1.0e-7;
            NBparam++;
        }



        FORCE_MAKE_PIAA0shape = 1;
        FORCE_MAKE_PIAA1shape = 1;
        break;



    default :
        printERROR(__FILE__,__func__,__LINE__, "mode not recognized");
        break;
    }

































    if(LINOPT == 1) // linear optimization
    {

        // Compute Reference
        PIAACMCsimul_makePIAAshapes();
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);

        printf("Reference = %g\n", valref);
        chname_image_ID("imvect", "vecDHref");
        ID = image_ID("vecDHref");
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];

        sprintf(fname, "!%s/vecDMref.fits", piaacmcconfdir);
        save_fits("vecDHref", fname);

        size1Dvec = data.image[ID].md[0].nelement;
        if(REGPIAASHAPES==1)
        {
            size1Dvec += data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa1CmodesID].md[0].size[0];
        }


        // re-package vector into 1D array and add regularization terms
        IDm = create_2Dimage_ID("DHmask", size1Dvec, 1);
        ID1Dref = create_2Dimage_ID("vecDHref1D", size1Dvec, 1);


        ID = image_ID("vecDHref");
        for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
        {
            data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
            data.image[IDm].array.F[ii] = 1.0;
        }
        if(REGPIAASHAPES == 1)
        {
            ID = piaacmc[0].piaa0CmodesID;
            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa0C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa0C_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }

            ID = piaacmc[0].piaa1CmodesID;
            for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa1C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa1C_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }
        }
        delete_image_ID("vecDHref");


	


        sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
        fp = fopen(fname, "w");
        fclose(fp);

list_image_ID();
        //
        // LINEAR OPTIMIZATION AROUND CURRENT POINT
        //

        for(iter=0; iter<NBiter; iter++)
        {
			printf("Iteration %ld/%ld\n", iter, NBiter);
            IDmodes = create_3Dimage_ID("DHmodes", size1Dvec, 1, NBparam);



			// compute local derivatives
			if(PIAACMC_FPM_FASTDERIVATIVES == 1)
			{
				ID = create_2Dimage_ID("DHmodes2Dtest", size1Dvec, NBparam);
	for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++)
		{
			PIAACMCsimul_achromFPMsol_eval_zonezderivative(mz, fpmresp_array, zonez_array, dphadz_array, outtmp_array, vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
			for(ii=0; ii<data.image[ID1D].md[0].nelement; ii++)
                    data.image[ID].array.F[mz*data.image[ID1D].md[0].nelement+ii] = outtmp_array[ii]*paramdelta[mz];
		}
	            sprintf(fname, "!%s/DHmodes_test.fits", piaacmcconfdir);
                save_fits("DHmodes2Dtest", fname);

                delete_image_ID("DHmodes2Dtest");
			}
			else
			{
            for(i=0; i<NBparam; i++)
            {
                optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                if(paramtype[i]==FLOAT)
                {
                    *(paramvalf[i]) += (float) paramdelta[i];
                    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
                }
                else
                {
                    *(paramval[i]) += paramdelta[i];
                    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
                }

          //      sprintf(fname,"!%s/imvect_%02ld.fits", piaacmcconfdir, i);
         //       save_fits("imvect", fname);
                ID = image_ID("imvect");


                sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
                fp = fopen(fname, "a");
                fprintf(fp, "# %5ld/%5ld %5ld/%5ld %20g %20g \n", iter, NBiter, i, NBparam, val, valref);
                fclose(fp);



                // re-package vector into 1D array and add regularization terms
                ID1D = create_2Dimage_ID("imvect1D", size1Dvec, 1);

                for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
                    data.image[ID1D].array.F[ii] = data.image[ID].array.F[ii];

                if(REGPIAASHAPES==1)
                {
                    ID = piaacmc[0].piaa0CmodesID;
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    {
                        data.image[ID1D].array.F[ii] = piaa0C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa0C_regcoeff_alpha);
                        ii++;
                    }

                    ID = piaacmc[0].piaa1CmodesID;
                    for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                    {
                        data.image[ID1D].array.F[ii] = piaa1C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa1C_regcoeff_alpha);
                        ii++;
                    }
                }
                delete_image_ID("imvect");



                if(paramtype[i]==FLOAT)
                    *(paramvalf[i]) -= (float) paramdelta[i];
                else
                    *(paramval[i]) -= paramdelta[i];



                for(ii=0; ii<data.image[ID1D].md[0].nelement; ii++)
                    data.image[IDmodes].array.F[i*data.image[ID1D].md[0].nelement+ii] = (data.image[ID1D].array.F[ii] - data.image[ID1Dref].array.F[ii]);


            //    printf("%3ld %g %g\n", i, val, valref);


                ID = create_2Dimage_ID("DHmodes2D", size1Dvec, NBparam);
                for(ii=0; ii<data.image[IDmodes].md[0].nelement; ii++)
                    data.image[ID].array.F[ii] = data.image[IDmodes].array.F[ii];

                sprintf(fname, "!%s/DMmodes.fits", piaacmcconfdir);
                save_fits("DHmodes2D", fname);

                delete_image_ID("DHmodes2D");
            }
			}
			

            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 1.0e-5, "optcoeff", 0);
           // sprintf(command, "mv eigenv.dat %s/eigenv_DH.dat", piaacmcconfdir);
            //ret = system(command);


            // delete_image_ID("vecDHref1D");

            ID = image_ID("optcoeff");

            // do linear scan
            linscanOK = 1;
            scangain = 0.0; //scanstepgain;
            val = 100000000000.0; // big number
            bestgain = 0.0;
            k = 0;
            while(linscanOK==1)
            {
                linoptlimflagarray[k] = 0;
                for(i=0; i<NBparam; i++)
                {
                    paramdeltaval[i] = -scangain*data.image[ID].array.F[i]*paramdelta[i];
                    if(paramdeltaval[i]<-parammaxstep[i])
                    {
                        paramdeltaval[i] = -parammaxstep[i];
                        linoptlimflagarray[k] = 1;
                    }
                    if(paramdeltaval[i]>parammaxstep[i])
                    {
                        paramdeltaval[i] = parammaxstep[i];
                        linoptlimflagarray[k] = 1;
                    }
                    if(paramtype[i]==FLOAT)
                        *(paramvalf[i]) += (float) paramdeltaval[i];
                    else
                        *(paramval[i]) += paramdeltaval[i];

                    if(paramtype[i]==FLOAT)
						printf("PARAM %ld  %g\n", i, *(paramvalf[i]));
					else
						printf("PARAM %ld  %g\n", i, *(paramval[i]));
                }
                valold = val;

				


                val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);

                sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
                fp = fopen(fname, "a");
                fprintf(fp, "##  %20lf %20g            %20g   %20g\n", scangain, val, data.image[piaacmc[0].zoneaID].array.D[0], *(paramval[0]));
                fclose(fp);


                for(i=0; i<NBparam; i++)
                {
                    if(paramtype[i]==FLOAT)
                        *(paramvalf[i]) -= (float) paramdeltaval[i];
                    else
                        *(paramval[i]) -= paramdeltaval[i];
                }

                linoptgainarray[k] = scangain;
                linoptvalarray[k] = val;
                k++;


                if(val<valold)
                {
                    linscanOK = 1;
                    bestgain = scangain;
                    scangain += scanstepgain;
                }
                else
                    linscanOK = 0;

                if(k>90)
                    linscanOK = 0;
                scangain += scanstepgain;
            }
            NBlinoptgain = k;
            for(k=0; k<NBlinoptgain; k++)
                printf("%2ld %10f %20g %d\n", k, linoptgainarray[k], linoptvalarray[k], linoptlimflagarray[k]);


            for(i=0; i<NBparam; i++)
            {
                paramdeltaval[i] = -bestgain*data.image[ID].array.F[i]*paramdelta[i];
                if(paramdeltaval[i]<-parammaxstep[i])
                    paramdeltaval[i] = -parammaxstep[i];
                if(paramdeltaval[i]>parammaxstep[i])
                    paramdeltaval[i] = parammaxstep[i];

                if(paramtype[i]==FLOAT)
                    *(paramvalf[i]) += (float) paramdeltaval[i];
                else
                    *(paramval[i]) += paramdeltaval[i];
            }
            valold = val;

            val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem);
            printf("gain: %lf -> val = %20g\n", bestgain, val);



            r = system("cp psfi.fits psfi_ref.fits");
            printf(" %g -> %g\n", valref, val);

            delete_image_ID("DHmodes");


            ID1Dref = image_ID("vecDHref1D"); //create_2Dimage_ID("vecDHref1D", size1Dvec, 1);
            ID = image_ID("imvect");
            for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
                data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];


            if(REGPIAASHAPES==1)
            {
                ID = piaacmc[0].piaa0CmodesID;
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa0C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa0C_regcoeff_alpha);
                    ii++;
                }

                ID = piaacmc[0].piaa1CmodesID;
                for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa1C_regcoeff*data.image[ID].array.F[jj]*pow(1.0*jj,piaa1C_regcoeff_alpha);
                    ii++;
                }
            }
            delete_image_ID("imvect");

            printf("Writing results\n");
            fflush(stdout);

            ID1Dref = image_ID("vecDHref1D");

            ID = image_ID("optcoeff");
			
			sprintf(fname, "%s/param.opt", piaacmcconfdir);
            fp = fopen(fname, "w");
            for(i=0; i<NBparam; i++)
            {
                if(paramtype[i]==FLOAT)
                    fprintf(fp, "%5ld %20g %20g %20g %20g\n", i, *(paramvalf[i]), data.image[ID].array.F[i], data.image[ID].array.F[i]*paramdelta[i], paramdeltaval[i]);
                else
                    fprintf(fp, "%5ld %20g %20g %20g %20g\n", i, *(paramval[i]), data.image[ID].array.F[i], data.image[ID].array.F[i]*paramdelta[i], paramdeltaval[i]);

            }
            fclose(fp);
            delete_image_ID("optcoeff");

            sprintf(fname, "%s/linoptval.txt", piaacmcconfdir);
            fp = fopen(fname, "a");
            fprintf(fp, "%5ld %20g %20g                %20g  %20g\n", iter, val, valref, data.image[piaacmc[0].zoneaID].array.D[0], *(paramval[0]));
            fclose(fp);
			
			PIAACMCSIMUL_VAL = val;
			PIAACMCSIMUL_VALREF = valref;


			sprintf(dirname, "%s_linopt", piaacmcconfdir);
			PIAAsimul_savepiaacmcconf(dirname); // staging area
			sprintf(command, "rsync -au --progress %s/* ./%s/", dirname, piaacmcconfdir);
			r = system(command);
        }
    }





    //  PIAAsimul_savepiaacmcconf("piaacmc0");
    //  PIAAsimul_loadpiaacmcconf("piaacmc0");
    // PIAAsimul_savepiaacmcconf("piaacmc1");
    //exit(0);



    if(0) // Lyot mask #0 position
    {
        paramtype[NBparam] = DOUBLE;
        paramval[NBparam] = &piaacmc[0].LyotStop_zpos[0];
        paramdelta[NBparam] = 0.05;
        parammaxstep[NBparam] = 0.05;
        parammin[NBparam] = 0.0;
        parammax[NBparam] = 2.5;
        NBparam++;
    }

    if(0) // Lyot mask #1 position
    {
        paramtype[NBparam] = DOUBLE;
        paramval[NBparam] = &piaacmc[0].LyotStop_zpos[1];
        paramdelta[NBparam] = 0.05;
        parammaxstep[NBparam] = 0.05;
        parammin[NBparam] = -0.5;
        parammax[NBparam] = 0.5;
        NBparam++;
    }

    if(0) // Focal plane mask radius
    {
        paramtype[NBparam] = DOUBLE;
        paramval[NBparam] = &piaacmc[0].fpmRad;
        paramdelta[NBparam] = 1.0e-6;
        parammaxstep[NBparam] = 5.0e-6;
        parammin[NBparam] = 1.0e-6;
        parammax[NBparam] = 1.0e-4;
        NBparam++;
    }

    if(0) // Focal plane material thickness
    {
        for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
        {
            paramtype[NBparam] = DOUBLE;
            paramval[NBparam] = &data.image[piaacmc[0].zonezID].array.D[k];
            paramdelta[NBparam] = 1.0e-9;
            parammaxstep[NBparam] = 1.0e-7;
            parammin[NBparam] = -1.0e-5;
            parammax[NBparam] = 1.0e-5;
            NBparam++;
        }
    }

    if(0) // Focal plane material transmission
    {
        for(k=0; k<data.image[piaacmc[0].zoneaID].md[0].size[0]; k++)
        {
            paramtype[NBparam] = DOUBLE;
            paramval[NBparam] = &data.image[piaacmc[0].zoneaID].array.D[k];
            paramdelta[NBparam] = 1.0e-4;
            parammaxstep[NBparam] = 5.0e-2;
            parammin[NBparam] = 1.0e-5;
            parammax[NBparam] = 0.99;
            NBparam++;
        }
    }




 




    free(piaacmc);

    return 0;
}




int PIAACMCsimul_run(long confindex, long mode)
{
long i;
FILE *fp;
char fname[200];
double bestval = 1.0;
int ret;
char command[500];


	if(mode==13)
	{
		sprintf(fname,"mode13.opt.txt");
		fp = fopen(fname, "w");
		fclose(fp);
		for(i=0;i<100000;i++)
		{
			if(i<3)
				MODampl = 0.0;
			else
				MODampl = 1.0e-6*ran1()*ran1()*ran1();
			PIAACMCsimul_exec(confindex, mode);
			if(PIAACMCSIMUL_VAL<bestval)
				{
					bestval = PIAACMCSIMUL_VAL;
					sprintf(command, "cp %s/fpm_zonez.fits %s/fpm_zonez_best.fits", piaacmcconfdir, piaacmcconfdir);
					ret = system(command);
				}
		fp = fopen(fname, "a");
		fprintf(fp,"%20.5g %20.5g %20.5g %20.5g\n", MODampl, PIAACMCSIMUL_VALREF, PIAACMCSIMUL_VAL, bestval);
		fclose(fp);
	
		sprintf(command, "cp %s/fpm_zonez_best.fits %s/fpm_zonez.fits", piaacmcconfdir, piaacmcconfdir);
		ret = system(command);		
		}
	}
	else
		PIAACMCsimul_exec(confindex, mode);
		
	return 0;
	}


