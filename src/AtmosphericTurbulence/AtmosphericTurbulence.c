#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <ctype.h>
#include <math.h>


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
    return 0;
}
#else
#include <time.h>
#endif



#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"

#include "fft/fft.h"
#include "image_basic/image_basic.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "image_filter/image_filter.h"
#include "info/info.h"
#include "image_gen/image_gen.h"
#include "statistic/statistic.h"

#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "psf/psf.h"
#include "WFpropagate/WFpropagate.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "AtmosphereModel/AtmosphereModel.h"


#ifdef _OPENMP
#include <omp.h>
#endif



#define SWAP(x,y)  temp=(x);x=(y);y=temp;

#define PI 3.14159265358979323846264338328


extern DATA data;

char CONFFILE[200] = "WFsim.conf";
 
//int TimeDayOfYear;
//float TimeLocalSolarTime;


float SiteLat;
float SiteLong;
float SiteAlt;
//float CO2_ppm;

//float SiteH2OMethod;
//float SiteTPW;
//float SiteRH;
//float SitePWSH;
//float alpha1H2O;






// constants
double C_me = 9.10938291e-31; // electon mass [kg]
double C_e0 = 8.854187817620e-12; // Vacuum permittivity [F.m-1]
double C_Na = 6.0221413e23; // Avogadro number
double C_e = 1.60217657e-19; // electron charge [C]
double C_ls = 2.686777447e25; // Loschmidt constant


double rhocoeff = 1.0;








// CONFIGURATION

// ------------ TURBULENCE AND ATMOSPHERE PARAMETERS -------------------------------------------
float CONF_LAMBDA;
float CONF_SEEING;
char CONF_TURBULENCE_PROF_FILE[200];
float CONF_ZANGLE;
float CONF_SOURCE_Xpos;
float CONF_SOURCE_Ypos;


// ------------ LOW_FIDELITY OUTPUT AT REF LAMBDA ------------------------------------------------

int CONF_WFOUTPUT = 1;
char CONF_WF_FILE_PREFIX[200];
int CONF_SHM_OUTPUT = 0;


// ------------- HIGH FIDELITY OUTPUT WAVELENGTH ---------------------------------

int CONF_MAKE_SWAVEFRONT = 0;
float CONF_SLAMBDA;
int CONF_SWF_WRITE2DISK = 0;
char CONF_SWF_FILE_PREFIX[200];

int CONF_SHM_SOUTPUT = 0;
char CONF_SHM_SPREFIX[100];
int CONF_SHM_SOUTPUTM = 0; // 1: output in [meter]


// ------------ OUTPUT PARAMETERS --------------------------------------------
long CONF_WFsize;
float CONF_PUPIL_SCALE;


// ------------- TIMING ---------------------------------------
int CONF_ATMWF_REALTIME;
double CONF_ATMWF_REALTIMEFACTOR = 1.0;
float CONF_WFTIME_STEP;
float CONF_TIME_SPAN;
long CONF_NB_TSPAN;
long CONF_SIMTDELAY = 0;
int CONF_WAITFORSEM = 0;
char CONF_WAITSEMIMNAME[100];


// ------------ COMPUTATION PARAMETERS, MODES --------------------------------
int CONF_SKIP_EXISTING;
long CONF_WF_RAW_SIZE;
long CONF_MASTER_SIZE;


// ------------ WAVEFRONT AMPLITUDE -------------------------------------------
int CONF_FRESNEL_PROPAGATION;
int CONF_WAVEFRONT_AMPLITUDE;
float CONF_FRESNEL_PROPAGATION_BIN;













// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


int make_AtmosphericTurbulence_wavefront_series_cli()
{
  make_AtmosphericTurbulence_wavefront_series();

  return 0;
}


int make_AtmosphericTurbulence_vonKarmanWind_cli()
{
    if(CLI_checkarg(1,2)+CLI_checkarg(2,1)+CLI_checkarg(3,1)+CLI_checkarg(4,1)+CLI_checkarg(5,2)+CLI_checkarg(6,3)==0)
        make_AtmosphericTurbulence_vonKarmanWind(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numf, data.cmdargtoken[5].val.numl, data.cmdargtoken[6].val.string);
    else
        return 1;
}   


int AtmosphericTurbulence_mkmastert_cli()
{
  
  if(CLI_checkarg(1,3)+CLI_checkarg(2,3)+CLI_checkarg(3,2)+CLI_checkarg(4,1)+CLI_checkarg(5,1)==0)
    {
      make_master_turbulence_screen(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl, 1.0*data.cmdargtoken[4].val.numf, 1.0*data.cmdargtoken[5].val.numf);
    }
  else
    return 1;
}

int AtmosphericTurbulence_makeHV_CN2prof_cli()
{
  if(CLI_checkarg(1,1)+CLI_checkarg(2,1)+CLI_checkarg(3,1)+CLI_checkarg(4,2)+CLI_checkarg(5,3)==0)
    {
      AtmosphericTurbulence_makeHV_CN2prof(data.cmdargtoken[1].val.numf, data.cmdargtoken[2].val.numf, data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numl, data.cmdargtoken[5].val.string);
    }
  else
    return 1;
}



//int AtmosphericTurbulence_makeHV_CN2prof(double wspeed, double r0, double sitealt, char *outfile)


int init_AtmosphericTurbulence()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "Atmospheric Turbulence");
  data.NBmodule++;


  strcpy(data.cmd[data.NBcmd].key,"mkwfs");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = make_AtmosphericTurbulence_wavefront_series_cli;
  strcpy(data.cmd[data.NBcmd].info,"make wavefront series");
  strcpy(data.cmd[data.NBcmd].syntax,"no arguments");
  strcpy(data.cmd[data.NBcmd].example,"mkwfs");
  strcpy(data.cmd[data.NBcmd].Ccall,"int make_AtmosphericTurbulence_wavefront_series()");
  data.NBcmd++;


    strcpy(data.cmd[data.NBcmd].key,"mkvonKarmanWind");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = make_AtmosphericTurbulence_vonKarmanWind_cli;
  strcpy(data.cmd[data.NBcmd].info,"make vonKarman wind model");
  strcpy(data.cmd[data.NBcmd].syntax,"<pixsize> <pixscale [m/pix]> <sigma windspeed [m/s]> <scale [m]> <size [long]> <output name>");
  strcpy(data.cmd[data.NBcmd].example,"mkvonKarmanWind 8192 0.1 20.0 50.0 512 vKmodel");
  strcpy(data.cmd[data.NBcmd].Ccall,"long make_AtmosphericTurbulence_vonKarmanWind(long vKsize, float pixscale, float sigmawind, float Lwind, long size, char *IDout_name)");
  data.NBcmd++;

 


  strcpy(data.cmd[data.NBcmd].key,"mkmastert");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AtmosphericTurbulence_mkmastert_cli;
  strcpy(data.cmd[data.NBcmd].info,"make 2 master phase screens");
  strcpy(data.cmd[data.NBcmd].syntax,"<screen0> <screen1> <size> <outerscale> <innerscale>");
  strcpy(data.cmd[data.NBcmd].example,"mkmastert scr0 scr1 2048 50.0 2.0");
  strcpy(data.cmd[data.NBcmd].Ccall,"int make_master_turbulence_screen(char *ID_name1, char *ID_name2, long size, float outercale, float innercale)");
  data.NBcmd++;


  strcpy(data.cmd[data.NBcmd].key,"mkHVturbprof");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = AtmosphericTurbulence_makeHV_CN2prof_cli;
  strcpy(data.cmd[data.NBcmd].info,"make Hufnager-Valley turbulence profile");
  strcpy(data.cmd[data.NBcmd].syntax,"<high wind speed [m/s]> <r0 [m]> <site alt [m]> <NBlayers> <output file>");
  strcpy(data.cmd[data.NBcmd].example,"mkHVturbprof 21.0 0.15 4200 100 turbHV.prof");
  strcpy(data.cmd[data.NBcmd].Ccall,"int AtmosphericTurbulence_makeHV_CN2prof(double wspeed, double r0, double sitealt, long NBlayer, char *outfile)");
  data.NBcmd++;



  return 0;
}


















int AtmosphericTurbulence_change_configuration_file(char *fname)
{
  sprintf(CONFFILE, "%s", fname);
  
  return(0);
}






//
// pixscale [m/pix]
// sigmawind [m/s]
// Lwind [m]
//
long make_AtmosphericTurbulence_vonKarmanWind(long vKsize, float pixscale, float sigmawind, float Lwind, long size, char *IDout_name)
{
    long ID, IDc;
    long ii, jj;
    double dx, dy, r;
    double Lwindpix;
    double rms = 0.0;
    
    double sigmau, sigmav, sigmaw;
    
    IDc = create_3Dimage_ID(IDout_name, vKsize, 1, 3);
    sigmau = sigmawind;
    sigmav = sigmawind;
    sigmaw = sigmawind;
    
    // longitudinal (u)
    make_rnd("tmppha0",vKsize, 1,"");
    arith_image_cstmult("tmppha0", 2.0*PI,"tmppha");
    delete_image_ID("tmppha0");
  
    printf("vK wind outer scale = %f m\n", Lwind);
    printf("pixscale            = %f m\n", pixscale);
    printf("Image size          = %f m\n", vKsize*pixscale);

  
    ID = create_2Dimage_ID("tmpamp0", vKsize, 1);
    for(ii=0; ii<vKsize; ii++)
        {
            dx = 1.0*ii-vKsize/2;
            r = sqrt(dx*dx); // period = (size*pixscale)/r
            // spatial frequency = 2 PI / period
            data.image[ID].array.F[ii] = sqrt( 1.0/pow(1.0 + pow(1.339*2.0*M_PI*r/(vKsize*pixscale)*Lwind,2.0), 5.0/6.0) );
        }
    
    
    make_rnd("tmpg", vKsize, 1,"-gauss");
    arith_image_mult("tmpg", "tmpamp0", "tmpamp");
    save_fits("tmpamp0","!vKwind_tmpamp0.fits");
    delete_image_ID("tmpamp0");
    delete_image_ID("tmpg");
    arith_set_pixel("tmpamp", 0.0, vKsize/2, 0); 
    mk_complex_from_amph("tmpamp", "tmppha", "tmpc", 0);
    delete_image_ID("tmpamp");
    delete_image_ID("tmppha");
    permut("tmpc");
    do2dfft("tmpc","tmpcf");
    delete_image_ID("tmpc");
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf");
    delete_image_ID("tmpo2");
    ID = image_ID("tmpo1");
    rms = 0.0;
    for(ii=0; ii<vKsize; ii++)
        rms += data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
    rms = sqrt(rms/vKsize);
    
    for(ii=0; ii<vKsize; ii++)
        data.image[IDc].array.F[ii] = data.image[ID].array.F[ii]/rms*sigmau;
    delete_image_ID("tmpo1");

  
     // tangential (v)
    make_rnd("tmppha0", vKsize, 1, "");
    arith_image_cstmult("tmppha0", 2.0*PI,"tmppha");
    delete_image_ID("tmppha0");
  
    ID = create_2Dimage_ID("tmpamp0", vKsize, 1);
    for(ii=0; ii<vKsize; ii++)
        {
            dx = 1.0*ii-vKsize/2;
            r = sqrt(dx*dx); 
            data.image[ID].array.F[ii] = sqrt( (1.0 + 8.0/3.0*pow(2.678*2.0*M_PI*r/(vKsize*pixscale)*Lwind,2.0)) /pow(1.0 + pow(2.678*2.0*M_PI*r/(vKsize*pixscale)*Lwind,2.0), 11.0/6.0) );
        }
    make_rnd("tmpg", vKsize, 1,"-gauss");
    arith_image_mult("tmpg", "tmpamp0", "tmpamp");
    delete_image_ID("tmpamp0");
    delete_image_ID("tmpg");
    arith_set_pixel("tmpamp", 0.0, vKsize/2, 0);
    mk_complex_from_amph("tmpamp", "tmppha", "tmpc", 0);
    delete_image_ID("tmpamp");
    delete_image_ID("tmppha");
    permut("tmpc");
    do2dfft("tmpc","tmpcf");
    delete_image_ID("tmpc");
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf");
    delete_image_ID("tmpo2");
    ID = image_ID("tmpo1");
    rms = 0.0;
    for(ii=0; ii<vKsize; ii++)
        rms += data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
    rms = sqrt(rms/(vKsize));
    
    for(ii=0; ii<vKsize; ii++)
        data.image[IDc].array.F[vKsize+ii] = data.image[ID].array.F[ii]/rms*sigmav;
    delete_image_ID("tmpo1");
      
      
      
      
    // vertical (w)
    make_rnd("tmppha0", vKsize, 1, "");
    arith_image_cstmult("tmppha0", 2.0*PI,"tmppha");
    delete_image_ID("tmppha0");
  
    ID = create_2Dimage_ID("tmpamp0", vKsize, 1);
    for(ii=0; ii<vKsize; ii++)
        {
            dx = 1.0*ii-size/2;
            r = sqrt(dx*dx);
            data.image[ID].array.F[ii] = sqrt( (1.0 + 8.0/3.0*pow(2.678*2.0*M_PI*r/(vKsize*pixscale)*Lwind,2.0)) /pow(1.0 + pow(2.678*2.0*M_PI*r/(vKsize*pixscale)*Lwind,2.0), 11.0/6.0) );
        }
    make_rnd("tmpg", vKsize, 1,"-gauss");
    arith_image_mult("tmpg", "tmpamp0", "tmpamp");
    delete_image_ID("tmpamp0");
    delete_image_ID("tmpg");
    arith_set_pixel("tmpamp", 0.0, vKsize/2, 0); 
    mk_complex_from_amph("tmpamp", "tmppha", "tmpc", 0);
    delete_image_ID("tmpamp");
    delete_image_ID("tmppha");
    permut("tmpc");
    do2dfft("tmpc","tmpcf");
    delete_image_ID("tmpc");
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf");
    delete_image_ID("tmpo2");
    ID = image_ID("tmpo1");
    rms = 0.0;
    for(ii=0; ii<vKsize; ii++)
        rms += data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
    rms = sqrt(rms/vKsize);
    
    for(ii=0; ii<size; ii++)
        data.image[IDc].array.F[vKsize*2+ii] = data.image[ID].array.F[ii]/rms*sigmav;
    delete_image_ID("tmpo1");

    
    
    
    return(IDc);
}




//
// innerscale and outerscale in pixel
// von Karman spectrum
//
int make_master_turbulence_screen(char *ID_name1, char *ID_name2, long size, float outerscale, float innerscale)
{
    long ID,ii,jj;
    float value,C1,C2;
    long cnt;
    long Dlim = 3;
    long IDv;

    int OUTERSCALE_MODE = 1; // 1 if outer scale
    double OUTERscale_f0;
    double INNERscale_f0;
    double dx, dy, r;
    double rlim = 0.0;
    int RLIMMODE = 0;
    double iscoeff;

    /*  IDv = variable_ID("OUTERSCALE");
      if(IDv!=-1)
        {
          outerscale = data.variable[IDv].value.f;
          printf("Outer scale = %f pix\n", outerscale);
        }
     */

    IDv = variable_ID("RLIM");
    if(IDv!=-1)
    {
        RLIMMODE = 1;
        rlim = data.variable[IDv].value.f;
        printf("R limit = %f pix\n",rlim);
    }

    OUTERscale_f0 = 1.0*size/outerscale; // [1/pix] in F plane
    INNERscale_f0 = (5.92/(2.0*M_PI))*size/innerscale;

    make_rnd("tmppha",size,size,"");
    arith_image_cstmult("tmppha", 2.0*PI,"tmppha1");
    delete_image_ID("tmppha");
    //  make_dist("tmpd",size,size,size/2,size/2);
    ID = create_2Dimage_ID("tmpd",size,size);
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            dx = 1.0*ii-size/2;
            dy = 1.0*jj-size/2;

            if(RLIMMODE==1)
            {
                r = sqrt(dx*dx + dy*dy);
                if(r<rlim)
                    data.image[ID].array.F[jj*size+ii] = 0.0;
                else
                    data.image[ID].array.F[jj*size+ii] = sqrt(dx*dx + dy*dy + OUTERscale_f0*OUTERscale_f0);
            }
            else
                data.image[ID].array.F[jj*size+ii] = sqrt(dx*dx + dy*dy + OUTERscale_f0*OUTERscale_f0);
        }
    //  data.image[ID].array.F[size/2*size+size/2+10] = 1.0;

    // period [pix] = size/sqrt(dx*dx+dy*dy)
    // f [1/pix] = sqrt(dx*dx+dy*dy)/size
    // f [1/pix] * size = sqrt(dx*dx+dy*dy)

    make_rnd("tmpg",size,size,"-gauss");
    ID = image_ID("tmpg");
    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            dx = 1.0*ii-size/2;
            dy = 1.0*jj-size/2;
            iscoeff = exp(-(dx*dx+dy*dy)/INNERscale_f0/INNERscale_f0);
            data.image[ID].array.F[jj*size+ii] *= sqrt(iscoeff); // power -> amplitude : sqrt 
        }

    arith_image_cstpow("tmpd", 11.0/6.0, "tmpd1");
    delete_image_ID("tmpd");
    arith_image_div("tmpg", "tmpd1", "tmpamp");
    delete_image_ID("tmpg");
    delete_image_ID("tmpd1");
    arith_set_pixel("tmpamp", 0.0, size/2, size/2);
    mk_complex_from_amph("tmpamp", "tmppha1", "tmpc", 0);
    delete_image_ID("tmpamp");
    delete_image_ID("tmppha1");
    permut("tmpc");
    do2dfft("tmpc","tmpcf");
    delete_image_ID("tmpc");
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf");

    /* compute the scaling factor in the power law of the structure function */
    fft_structure_function("tmpo1", "strf");
    ID = image_ID("strf");
    value = 0.0;
    cnt = 0;
    for(ii = 1; ii<Dlim; ii++)
        for(jj = 1; jj<Dlim; jj++)
        {
            value += log10(data.image[ID].array.F[jj*size+ii])-5.0/3.0*log10(sqrt(ii*ii+jj*jj));
            cnt++;
        }
    // save_fl_fits("strf","!strf.fits");
    delete_image_ID("strf");
    C1 = pow(10.0,value/cnt);

    fft_structure_function("tmpo2", "strf");
    ID=image_ID("strf");
    value = 0.0;
    cnt = 0;
    for(ii=1; ii<Dlim; ii++)
        for(jj=1; jj<Dlim; jj++)
        {
            value += log10(data.image[ID].array.F[jj*size+ii])-5.0/3.0*log10(sqrt(ii*ii+jj*jj));
            cnt++;
        }
    delete_image_ID("strf");
    C2 = pow(10.0,value/cnt);

    printf("%f %f\n", C1, C2);

    arith_image_cstmult("tmpo1",1.0/sqrt(C1),ID_name1);
    arith_image_cstmult("tmpo2",1.0/sqrt(C2),ID_name2);
    delete_image_ID("tmpo1");
    delete_image_ID("tmpo2");

    return(0);
}





int make_master_turbulence_screen_pow(char *ID_name1, char *ID_name2, long size, float power)
{
    long ID,ii,jj;
    float value,C1,C2;
    long cnt;
    long Dlim = 3;

    make_rnd("tmppha",size,size,"");
    arith_image_cstmult("tmppha",2.0*PI,"tmppha1");
    delete_image_ID("tmppha");
    make_dist("tmpd",size,size,size/2,size/2);
    make_rnd("tmpg",size,size,"-gauss");

    arith_image_cstpow("tmpd",power,"tmpd1");
    delete_image_ID("tmpd");
    arith_image_div("tmpg","tmpd1","tmpamp");
    delete_image_ID("tmpg");
    delete_image_ID("tmpd1");
    arith_set_pixel("tmpamp",0.0,size/2,size/2);
    mk_complex_from_amph("tmpamp","tmppha1","tmpc", 0);
    delete_image_ID("tmpamp");
    delete_image_ID("tmppha1");
    permut("tmpc");
    do2dfft("tmpc","tmpcf");
    delete_image_ID("tmpc");
    mk_reim_from_complex("tmpcf","tmpo1","tmpo2", 0);
    delete_image_ID("tmpcf");

    /* compute the scaling factor in the power law of the structure function */
    fft_structure_function("tmpo1","strf");
    ID=image_ID("strf");
    value = 0.0;
    cnt = 0;
    for(ii=1; ii<Dlim; ii++)
        for(jj=1; jj<Dlim; jj++)
        {
            value += log10(data.image[ID].array.F[jj*size+ii])-power*log10(sqrt(ii*ii+jj*jj));
            /*	printf("%ld %ld %f\n",ii,jj,log10(data.image[ID].array.F[jj*size+ii])-5.0/3.0*log10(sqrt(ii*ii+jj*jj)));*/
            cnt++;
        }
    delete_image_ID("strf");
    C1=pow(10.0,value/cnt);

    fft_structure_function("tmpo2","strf");
    ID=image_ID("strf");
    value = 0.0;
    cnt = 0;
    for(ii=1; ii<Dlim; ii++)
        for(jj=1; jj<Dlim; jj++)
        {
            value += log10(data.image[ID].array.F[jj*size+ii])-power*log10(sqrt(ii*ii+jj*jj));
            cnt++;
        }
    delete_image_ID("strf");
    C2=pow(10.0,value/cnt);
    /*  printf("%f %f\n",C1,C2);*/
    arith_image_cstmult("tmpo1",1.0/sqrt(C1),ID_name1);
    arith_image_cstmult("tmpo2",1.0/sqrt(C2),ID_name2);
    delete_image_ID("tmpo1");
    delete_image_ID("tmpo2");

    return(0);
}






int contract_wavefront_cube(char *ina_file, char *inp_file, char *outa_file, char *outp_file, int factor)
{
    /* contracts the wavefront series by a factor of 2^factor */
    long IDamp,IDpha,IDoutamp,IDoutpha;
    long ii,jj,kk;
    long i,j;
    long naxes[3];
    long naxes_out[3];
    float re, im, amp, pha;
    float P;
    long LARGE = 10000;
    float pharef, ampref;
    int pfactor;

    pfactor=1;
    for(i=0; i<factor; i++)
        pfactor *= 2;

    load_fits(inp_file, "tmpwfp", 1);
    IDpha=image_ID("tmpwfp");
    load_fits(ina_file, "tmpwfa", 1);
    IDamp=image_ID("tmpwfa");
    naxes[0] = data.image[IDpha].md[0].size[0];
    naxes[1] = data.image[IDpha].md[0].size[1];
    naxes[2] = data.image[IDpha].md[0].size[2];
    naxes_out[0] = data.image[IDpha].md[0].size[0]/pfactor;
    naxes_out[1] = data.image[IDpha].md[0].size[1]/pfactor;
    naxes_out[2] = data.image[IDpha].md[0].size[2];
    IDoutpha = create_3Dimage_ID("tmpwfop",naxes_out[0],naxes_out[1],naxes_out[2]);
    IDoutamp = create_3Dimage_ID("tmpwfoa",naxes_out[0],naxes_out[1],naxes_out[2]);

    ii=0;
    jj=0;
    kk=0;
    amp = 0.0;
    pha = 0.0;

    //  # ifdef _OPENMP
    //  #pragma omp parallel
    //  {
    //  # endif

    //  # ifdef _OPENMP
    //  #pragma omp for private(kk,ii,jj,i,j,amp,pha,ampref,pharef,re,im,P) collapse(3)
    //  # endif
    for(kk=0; kk<naxes[2]; kk++)
    {
        for(ii=0; ii<naxes[0]/pfactor; ii++)
            for(jj=0; jj<naxes[1]/pfactor; jj++)
            {
                re=0.0;
                im=0.0;
                pharef = 0.0;
                ampref = 0.0;
                for(i=0; i<pfactor; i++)
                    for(j=0; j<pfactor; j++)
                    {
                        amp = data.image[IDamp].array.F[kk*naxes[0]*naxes[1]+(pfactor*jj+j)*naxes[0]+pfactor*ii+i];
                        pha = data.image[IDpha].array.F[kk*naxes[0]*naxes[1]+(pfactor*jj+j)*naxes[0]+pfactor*ii+i];
                        pharef += data.image[IDamp].array.F[kk*naxes[0]*naxes[1]+(pfactor*jj+j)*naxes[0]+pfactor*ii+i]*data.image[IDpha].array.F[kk*naxes[0]*naxes[1]+(pfactor*jj+j)*naxes[0]+pfactor*ii+i];
                        ampref += data.image[IDamp].array.F[kk*naxes[0]*naxes[1]+(pfactor*jj+j)*naxes[0]+pfactor*ii+i];
                        re += amp*cos(pha);
                        im += amp*sin(pha);
                    }
                amp = sqrt(re*re+im*im);
                pha = atan2(im,re);
                pharef /= ampref;
                P = 2.0*PI*( ((long) (0.5+1.0*LARGE+(pharef-pha)/2.0/PI)) - LARGE);
                if(ampref<0.01)
                    P = 0.0;
                data.image[IDoutpha].array.F[kk*naxes_out[0]*naxes_out[1]+jj*naxes_out[0]+ii] = pha+P;
                data.image[IDoutamp].array.F[kk*naxes_out[0]*naxes_out[1]+jj*naxes_out[0]+ii] = amp/pfactor/pfactor;
            }
    }

    //  # ifdef _OPENMP
    //  }
    //  # endif


    save_fl_fits("tmpwfop",outp_file);
    save_fl_fits("tmpwfoa",outa_file);

    delete_image_ID("tmpwfa");
    delete_image_ID("tmpwfp");
    delete_image_ID("tmpwfoa");
    delete_image_ID("tmpwfop");


    return(0);
}


int contract_wavefront_cube_phaseonly(char *inp_file, char *outp_file, int factor)
{
    /* contracts the wavefront series by a factor of 2^factor */
    long IDpha,IDoutpha;
    long ii,jj,kk;
    long i,j;
    long naxes[3];
    long naxes_out[3];
    float re,im,amp,pha;
    float P;
    long LARGE = 10000;
    float pharef,ampref;
    int pfactor;
    long l1,l2,l3,l4;

    //  printf("CONTRACT, FACTOR = %d\n",factor);
    //fflush(stdout);

    pfactor=1;
    for(i=0; i<factor; i++)
        pfactor *= 2;

    load_fits(inp_file, "tmpwfp", 1);
    IDpha=image_ID("tmpwfp");
    naxes[0] = data.image[IDpha].md[0].size[0];
    naxes[1] = data.image[IDpha].md[0].size[1];
    naxes[2] = data.image[IDpha].md[0].size[2];
    naxes_out[0] = data.image[IDpha].md[0].size[0]/pfactor;
    naxes_out[1] = data.image[IDpha].md[0].size[1]/pfactor;
    naxes_out[2] = data.image[IDpha].md[0].size[2];
    IDoutpha = create_3Dimage_ID("tmpwfop",naxes_out[0],naxes_out[1],naxes_out[2]);

    ii=0;
    jj=0;
    kk=0;
    amp = 0.0;
    pha = 0.0;
    ampref = 1.0*pfactor*pfactor;

# ifdef _OPENMP
    #pragma omp parallel for shared(naxes,naxes_out,data,pfactor,IDoutpha,IDpha,ampref) private(kk,ii,jj,i,j,amp,pha,pharef,re,im,P,l1,l2,l3,l4)
# endif
    for(kk=0; kk<naxes[2]; kk++)
    {
        l1 = kk*naxes_out[0]*naxes_out[1];
        l2 = kk*naxes[0]*naxes[1];
        for(ii=0; ii<naxes[0]/pfactor; ii++)
            for(jj=0; jj<naxes[1]/pfactor; jj++)
            {
                re=0.0;
                im=0.0;
                pharef = 0.0;
                l4=l2+pfactor*ii;
                for(j=0; j<pfactor; j++)
                {
                    l3=l4+(pfactor*jj+j)*naxes[0];
                    for(i=0; i<pfactor; i++)
                    {
                        pha = data.image[IDpha].array.F[l3+i];
                        pharef += data.image[IDpha].array.F[l3+i];
                        re += cos(pha);
                        im += sin(pha);
                    }
                }
                /*	    amp = sqrt(re*re+im*im);*/
                pha = atan2(im,re);
                pharef /= ampref;
                P = 2.0*PI*( ((long) (0.5+1.0*LARGE+(pharef-pha)/2.0/PI)) - LARGE);
                if(ampref<0.01)
                    P = 0.0;
                data.image[IDoutpha].array.F[l1+jj*naxes_out[0]+ii] = pha+P;
            }

    }

    save_fl_fits("tmpwfop",outp_file);

    delete_image_ID("tmpwfp");
    delete_image_ID("tmpwfop");


    return(0);
}


// P [Pa]
// T [C]
// RH [%]

// ref: Techniques and Topics in Flow Measurement
// By Frank E. Jones
//
// coeffs An, Bn, Cn adopted from Nitrogen
//
// Compressibility factors:
// Davies 1991
//
double Z_Air(double P, double T, double RH)
{
    double Z, Z0, rho, rho0;
    double tc; // C
    double rhocoeff1;
    double An = 4.446e-6; 
    double Bn = 6.4e-13;
    double Cn = -1.07e-16;
    double P0 = 101325.0;
    
    double TK;
    double f; // enhancement factor
    double Psv; // water vapor saturation
    double xv;
    
    double A = 1.2378847e-5;
    double B = -1.9121316e-2;
    double C = 33.93711047;
    double D = -6.3431645e3;
    
    double alpha = 1.00062;
    double beta = 3.14e-8;
    double gamma = 5.6e-7;
    
    double a0 = 1.58123e-6;
    double a1 = -2.9331e-8;
    double a2 = 1.1043e-10;
    double b0 = 5.707e-6;
    double b1 = -2.051e-8;
    double c0 = 1.9898e-4;
    double c1 = -2.376e-6;
    double d = 1.83e-11;
    double e = -0.765e-8;
    
    double RH0 = 0.0;
    double T0, TK0, f0, Psv0, xv0;
    
    TK = T + 273.15;

    Z = 1.00001 - 5.8057e-9*P + 2.6402e-16*P*P - 3.3297e-7*T + 1.2420e-10*P*T - 2.0158e-18*P*P*T + 2.4925e-9*T*T - 6.2873e-13*P*T*T + 5.4174e-21*P*P*T*T - 3.5e-7*RH - 5.0e-9*RH*RH;
    
    Psv = exp(A*TK*TK + B*TK + C + D/TK);
    printf("Water vapor saturation pressure Psv0 = %f Pa\n", Psv);
    f = alpha + beta*P + gamma*T*T;
    printf("enhancement factor f = %f\n", f);
    xv = RH*0.01 * f * Psv/P;
    Z = 1.0 - P/TK * (a0 + a1*T + a2*T*T + (c0+b1*T)*xv + (c0+c1*T)*xv*xv) + P*P/TK/TK*(d+e*xv*xv);
    
    Z0 = 1.00001 - 5.8057e-9*P0 + 2.6402e-16*P0*P0;
    
    T0 = 0.0;
    TK0 = 273.15;
    Psv0 = exp(A*TK0*TK0 + B*TK0 + C + D/TK0);
    printf("Water vapor saturation pressure Psv0 = %f Pa\n", Psv0);
    f0 = alpha + beta*P0 + gamma*T0*T0;
    printf("enhancement factor f0 = %f\n", f0);
    xv0 = RH0*0.01 * f0 * Psv0/P0;
    Z0 = 1.0 - P0/TK0 * (a0 + a1*T0 + a2*T0*T0 + (c0+b1*T0)*xv0 + (c0+c1*T0)*xv0*xv0) + P0*P0/TK0/TK0*(d+e*xv0*xv0);
    
    
    rhocoeff1 = 101325.0/P;
    rhocoeff1 *= (T+273.15)/273.15;  // particle count ratio : more particles at lower temp
    rhocoeff1 *= Z/Z0;
    
    rho0 = C_ls/C_Na; // [mol.m^-3]
    
    rhocoeff = rhocoeff1;
    
    printf("Z = %.8f   Z0 = %.8f\n", Z, Z0);
    
    rho = rho0/rhocoeff1;
    rhocoeff = rhocoeff1; // * (1.0 + Bn/An*rho0 + Cn/An*rho0*rho0) / (1.0 + Bn/An*rho + Cn/An*rho*rho) / (1.0 + Bn/An*rho + Cn/An*rho*rho);
                
    return(Z);
}




double Z_N2(double P, double T)
{
    double Z, Z0, rho, rho0;
    double tc; // C
    double rhocoeff1;
    double An = 4.446e-6; // Achterman et al. 1986
    double Bn = 6.4e-13;
    double Cn = -1.07e-16;
    
    tc = T-273.15;

    // J. D. Wright, “Gas property equations for the NIST fluid flow group gas flow measurement calibration services” (National Institute of Standards and Technology, 2004), http://www.cstl.nist.gov/div836/836.01/PDFs/2004/Gas_Properties.pdf.)
    Z = 1.0 - 101325.0 * (P/101325.0) * (0.449805 - 0.01177*tc + 0.00006*tc*tc)*1e-8;
    Z0 = 1.0 - 101325.0 * (0.449805)*1e-8;
    
    rhocoeff1 = 101325.0/P;
    rhocoeff1 *= T/273.15;  // particle count ratio : more particles at lower temp
    rhocoeff1 *= Z/Z0;
    rho0 = C_ls/C_Na; // [mol.m^-3]
    rho = rho0/rhocoeff1;
    rhocoeff = rhocoeff1 * (1.0 + Bn/An*rho0 + Cn/An*rho0*rho0) / (1.0 + Bn/An*rho + Cn/An*rho*rho) / (1.0 + Bn/An*rho + Cn/An*rho*rho);
                
    return(Z);
}






int AtmosphericTurbulence_ReadConf()
{
    char KEYWORD[200];
    char CONTENT[200];
    
    
   // ------------ TURBULENCE AND ATMOSPHERE PARAMETERS -------------------------------------------
    
    strcpy(KEYWORD,"TURBULENCE_REF_WAVEL");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_LAMBDA = atof(CONTENT)*0.000001;
    
    strcpy(KEYWORD,"TURBULENCE_SEEING");
    read_config_parameter(CONFFILE,KEYWORD,CONTENT);
    CONF_SEEING = atof(CONTENT);

    strcpy(KEYWORD,"TURBULENCE_PROF_FILE");
    read_config_parameter(CONFFILE, KEYWORD, CONF_TURBULENCE_PROF_FILE);
 
    strcpy(KEYWORD,"ZENITH_ANGLE");
    read_config_parameter(CONFFILE,KEYWORD,CONTENT);
    CONF_ZANGLE = atof(CONTENT);

    strcpy(KEYWORD,"SOURCE_XPOS");
    read_config_parameter(CONFFILE,KEYWORD,CONTENT);
    CONF_SOURCE_Xpos = atof(CONTENT);

    strcpy(KEYWORD,"SOURCE_YPOS");
    read_config_parameter(CONFFILE,KEYWORD,CONTENT);
    CONF_SOURCE_Ypos = atof(CONTENT);

 
    
    // ------------ LOW_FIDELITY OUTPUT AT REF LAMBDA ------------------------------------------------
    
    strcpy(KEYWORD,"WFOUTPUT");
    if(read_config_parameter_exists(CONFFILE, KEYWORD)==1)
    {
        read_config_parameter(CONFFILE,KEYWORD, CONTENT);
        CONF_WFOUTPUT = atoi(CONTENT);
    }
    
    strcpy(KEYWORD,"WF_FILE_PREFIX");
    read_config_parameter(CONFFILE,KEYWORD, CONF_WF_FILE_PREFIX);

    strcpy(KEYWORD,"SHM_OUTPUT");
    if(read_config_parameter_exists(CONFFILE, KEYWORD)==1)
    {
        read_config_parameter(CONFFILE, KEYWORD, CONTENT);
        CONF_SHM_OUTPUT = atoi(CONTENT);
    }

    
    // ------------- HIGH FIDELITY OUTPUT WAVELENGTH ---------------------------------
    
    strcpy(KEYWORD,"MAKE_SWAVEFRONT");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_MAKE_SWAVEFRONT = atoi(CONTENT);

    strcpy(KEYWORD,"SLAMBDA");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_SLAMBDA = atof(CONTENT)*0.000001;

    strcpy(KEYWORD,"SWF_WRITE2DISK");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_SWF_WRITE2DISK = atoi(CONTENT);

    strcpy(KEYWORD,"SWF_FILE_PREFIX");
    read_config_parameter(CONFFILE, KEYWORD, CONF_SWF_FILE_PREFIX);

    strcpy(KEYWORD,"SHM_SOUTPUT");
    if(read_config_parameter_exists(CONFFILE, KEYWORD)==1)
    {
        read_config_parameter(CONFFILE, KEYWORD, CONTENT);
        CONF_SHM_SOUTPUT = atoi(CONTENT);
    }
    strcpy(KEYWORD,"SHM_SPREFIX");
    read_config_parameter(CONFFILE,KEYWORD, CONF_SHM_SPREFIX);

    strcpy(KEYWORD,"SHM_SOUTPUTM");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_SHM_SOUTPUTM = atoi(CONTENT);


    // ------------ OUTPUT PARAMETERS --------------------------------------------
    
    strcpy(KEYWORD,"PUPIL_SCALE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_PUPIL_SCALE = atof(CONTENT);
    
    strcpy(KEYWORD,"WFsize");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WFsize = atol(CONTENT);

 
    
    // ------------- TIMING ---------------------------------------
    
    strcpy(KEYWORD,"REALTIME");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_ATMWF_REALTIME = atoi(CONTENT);
    
    strcpy(KEYWORD,"REALTIMEFACTOR");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_ATMWF_REALTIMEFACTOR = atof(CONTENT);

    strcpy(KEYWORD, "WFTIME_STEP");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WFTIME_STEP = atof(CONTENT);

    strcpy(KEYWORD, "TIME_SPAN");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_TIME_SPAN = atof(CONTENT);

    strcpy(KEYWORD, "NB_TSPAN");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_NB_TSPAN = atol(CONTENT);

    strcpy(KEYWORD, "SIMTDELAY");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_SIMTDELAY = atol(CONTENT);

    
    strcpy(KEYWORD, "WAITFORSEM");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WAITFORSEM = atoi(CONTENT);
    
    strcpy(KEYWORD, "WAITSEMIMNAME");
    read_config_parameter(CONFFILE, KEYWORD, CONF_WAITSEMIMNAME);
    
    
    
    // ------------ COMPUTATION PARAMETERS, MODES --------------------------------
    
    
    strcpy(KEYWORD,"SKIP_EXISTING");
    if(read_config_parameter_exists(CONFFILE,KEYWORD)==1)
    {
        read_config_parameter(CONFFILE,KEYWORD,CONTENT);
        CONF_SKIP_EXISTING = 1;
    }
    else
        CONF_SKIP_EXISTING = 0;

    strcpy(KEYWORD,"WF_RAW_SIZE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WF_RAW_SIZE = atol(CONTENT);
 
    strcpy(KEYWORD,"MASTER_SIZE");
    read_config_parameter(CONFFILE, KEYWORD,CONTENT);
    CONF_MASTER_SIZE = atol(CONTENT);



    //  ------------ WAVEFRONT AMPLITUDE -------------------------------------------

    
    strcpy(KEYWORD,"WAVEFRONT_AMPLITUDE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WAVEFRONT_AMPLITUDE = atoi(CONTENT);

    strcpy(KEYWORD,"FRESNEL_PROPAGATION");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_FRESNEL_PROPAGATION = atoi(CONTENT);

    strcpy(KEYWORD,"FRESNEL_PROPAGATION_BIN");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_FRESNEL_PROPAGATION_BIN = atof(CONTENT);



    // ------------ POSTPROCESSING --------------
 

    strcpy(KEYWORD,"PUPIL_AMPL_FILE");
    if(read_config_parameter_exists(CONFFILE,KEYWORD)==1)
    {
        read_config_parameter(CONFFILE, KEYWORD, CONTENT);
        load_fits(CONTENT, "ST_pa", 1);
    }

    strcpy(KEYWORD,"PUPIL_PHA_FILE");
    if(read_config_parameter_exists(CONFFILE,KEYWORD)==1)
    {
        read_config_parameter(CONFFILE,KEYWORD,CONTENT);
        load_fits(CONTENT, "ST_pp", 1);
    }


    
    return 0;
}





int make_AtmosphericTurbulence_wavefront_series()
{
    //  fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status;
    long  fpixel = 1, naxis = 3, nelement;

    long naxes_MASTER[2];
    /*  long ID;*/
    long ID1,ID2;
    long tspan;
    char fname_a[200];
    char fname1[200];
    char fname1_a[200];
    char fname2[200];
    char tmpafname[200];
    char tmppfname[200];

    // float SODIUM_ALT;
    float layer_scale = 1.0;

    long master;
    long NBMASTERS;
    long *ID_TM;
    long *ID_TML;
    long IDout_array_pha;
    long IDout_array_amp;
    complex_float *array;

    // phase only
    long ID_array1;
    long ID_sarray1;
    long ID_carray1;

    // phase + amplitude
    long ID_array2;
    long ID_sarray2;
    long ID_carray2;

    float re,im;


    long IDout_sarray_pha;
    long IDout_sarray_amp;
    complex_float *sarray;
    float Scoeff;
    float Nlambda,Nslambda, l;


    // cone effect wavefront
    /*
        int make_cwavefront = 0;
        char CWF_FILE_PREFIX[100];
        complex_float *carray;
    */

    long NBLAYERS; /* number of layers */
    long layer;
    FILE *fp;
    char command[200];
    char line[2000];
    char word[200];
    char fname[200];
    char name[200];
    float *LAYER_ALT;
    float *LAYER_CN2;
    float *LAYER_SPD;
    float *LAYER_DIR;
    float *LAYER_OUTERSCALE;
    float *LAYER_INNERSCALE;
    float *LAYER_SIGMAWSPEED;
    float *LAYER_LWIND;
    int stop;

    long *naxes;
    long *naxesout;
    double *xpos;
    double *ypos;
    double *xpos0;
    double *ypos0;
    long *xposfcnt; // counter to keep trak of modulo
    long *yposfcnt;
    float *vxpix;
    float *vypix;
    double vpix,PA;
    long vindex;
    long frame;
    long NBFRAMES;
    float fl1, fl2, fl3, fl4, fl5, fl6, fl7, fl8;
    long ii, jj, iim, jjm, ii1, jj1;
    float value;
    float coeff=0.0;

    float *alt_bin_sep;
    long NB_alt_bin_sep;
    float minaltd;
    float *SLAYER_ALT;
    float *SLAYER_CN2;
    long NBSLAYERS;
    int OK;
    long i,j,index;
    long *super_layer_index;

    int contraction_factor;
    int pfactor;
    //   float OuterScale;
    float CN2total;
    float tmp,tmpf,h,Rindex,Roffset,RoffsetS;

    long xref,yref,xrefm,yrefm,iimax,jjmax,y1;
    long ID;
    long start_tspan=0;

    double tot,p1,p2,tmp1;
    long dpix,cnt,k,r0cnt;
    double r0,r0tot;

    long IDshmpha, IDshmamp;
    long IDshmspha, IDshmsamp;

    //
    // US standard atmosphere density, normalized to sea level, 2km step
    // log10 of density, 2km step
    //
    double StdAtmDens[44] = {0.000000, -0.085107, -0.174693, -0.268526, -0.367316, -0.471661, -0.594121, -0.730392, -0.866722, -1.003203, -1.139185, -1.278509, -1.416593, -1.553349, -1.688809, -1.823082, -1.956197, -2.093072, -2.227379, -2.358485, -2.486619, -2.611739, -2.734220, -2.854125, -2.968550, -3.076566, -3.182071, -3.282703, -3.385361, -3.490222, -3.597335, -3.706660, -3.818623, -3.933104, -4.050311, -4.170053, -4.293230, -4.421899, -4.552842, -4.686219, -4.822140, -4.960707, -5.101812, -5.245839};
    double dens_N2, dens_O2, dens_Ar, dens_H2O, dens_CO2, dens_Ne, dens_He, dens_CH4, dens_Kr, dens_H2, dens_O3, dens_N, dens_O, dens_H; // [cm-3]
    double lambda;
    double xN2, xO2, xAr, xH2O, xCO2, xNe, xHe, xCH4, xKr, xH2, xO3, xN, xO, xH;
    double P, T, TC, Pw, CO2ppm, denstot, xtot;
    double LoschmidtConstant =  2.6867805e25;
    double logD;
    double x;
    long i0,i1;
    double iimf, jjmf, iifrac, jjfrac, value_re, value_im;
    long iim1, jjm1;
    double pha;

    int r;
    double Temp;



    // timing
    struct timespec tnow;
    double tnowdouble;
    int kw;


    // optimal phase unwrapping
    double peakpha_re, peakpha_im, peakpha;
    long IDpeakpha_re, IDpeakpha_im, IDpeakpha;
    long IDpeakpha_re_bin, IDpeakpha_im_bin, IDpeakpha_bin, IDpeakpha_bin_ch;
    long xsizepeakpha;
    long ysizepeakpha;
    long ii2, jj2;
    long index1, index2;
    double pv1, pv2;
    long chiter = 0;
    long chcnt = 0;

    double plim = 0.0;
    double pcoeff2 = 0.15;
    long chcnt0cnt = 0;


    long WSPEEDsize;
    double WSPEEDpixscale;
    long vKwindsize;


    naxes = (long*) malloc(sizeof(long)*3);
    naxesout = (long*) malloc(sizeof(long)*3);



    ID1=-1;
    ID2=-1;
    ID_sarray1=-1;
    sarray = NULL;
    //  ID_carray1=-1;
    //  carray = NULL;






    printf("Making the wavefront series...\n");
    fflush(stdout);

    AtmosphericTurbulence_ReadConf();
    naxesout[0] = CONF_WFsize;
    naxesout[1] = CONF_WFsize;
    naxes[0] = CONF_WF_RAW_SIZE;
    naxes[1] = CONF_WF_RAW_SIZE;

 


    /* strcpy(KEYWORD,"TIME_DAY_OF_YEAR");
     read_config_parameter(CONFFILE,KEYWORD,CONTENT);
     TimeDayOfYear = atoi(CONTENT);

     strcpy(KEYWORD,"TIME_LOCAL_SOLAR_TIME");
     read_config_parameter(CONFFILE,KEYWORD,CONTENT);
     TimeLocalSolarTime= atoi(CONTENT);




     strcpy(KEYWORD,"SITE_LATITUDE");
     read_config_parameter(CONFFILE,KEYWORD,CONTENT);
     SiteLat = atof(CONTENT);

     strcpy(KEYWORD,"SITE_LONGITUDE");
     read_config_parameter(CONFFILE,KEYWORD,CONTENT);
     SiteLong = atof(CONTENT);


      strcpy(KEYWORD,"SITE_ALT");
     read_config_parameter(CONFFILE,KEYWORD,CONTENT);
     SiteAlt = atof(CONTENT);


     strcpy(KEYWORD,"CO2_PPM");
     read_config_parameter(CONFFILE,KEYWORD,CONTENT);
     CO2_ppm = atof(CONTENT);
    */

    // water

    /*
       strcpy(KEYWORD,"SITE_H2O_METHOD");
       read_config_parameter(CONFFILE,KEYWORD,CONTENT);
       SiteH2OMethod = atoi(CONTENT);

    strcpy(KEYWORD,"SITE_TPW");
       read_config_parameter(CONFFILE,KEYWORD,CONTENT);
       SiteTPW = atof(CONTENT);

    strcpy(KEYWORD,"SITE_RH");
       read_config_parameter(CONFFILE,KEYWORD,CONTENT);
       SiteRH = atof(CONTENT);

       strcpy(KEYWORD,"SITE_PW_SCALEH");
       read_config_parameter(CONFFILE,KEYWORD,CONTENT);
       SitePWSH = atof(CONTENT);
    */







    /*
        strcpy(KEYWORD,"MAKE_CWAVEFRONT");
        if(read_config_parameter_exists(CONFFILE,KEYWORD)==1)
        {
            read_config_parameter(CONFFILE,KEYWORD,CONTENT);
            make_cwavefront = atoi(CONTENT);
        }
        else
            make_cwavefront = 0;


        strcpy(KEYWORD,"CWF_FILE_PREFIX");
        if(read_config_parameter_exists(CONFFILE,KEYWORD)==1)
            read_config_parameter(CONFFILE,KEYWORD,CWF_FILE_PREFIX);
        else
            sprintf(CWF_FILE_PREFIX,"cwf_");
    */



    printf("Building reference atmosphere model ...\n");
    // create atm.txt file with concentrations as function of altitude
    // load RIA (Refractive index and Absorption) files if available 
    AtmosphereModel_Create_from_CONF(CONFFILE);


    // SOME TESTING
  //  AtmosphereModel_RefractionPath(1.5, CONF_ZANGLE, 0); //79.999/180.0*M_PI);//CONF_ZANGLE);
 

 /*   fp = fopen("Rprof.txt", "w");
    for(h=0; h<100000.0; h+=10.0)
        fprintf(fp, "%8g %.16f %.16f\n", h, AtmosphereModel_stdAtmModel_N(h, 0.6e-6), AtmosphereModel_stdAtmModel_N(h, 1.6e-6));
    fclose(fp);
*/

    // TESTING VALUES IN CIDDOR 1996
    P = 102993.0;
    Pw = 641.0;
    CO2ppm = 450.0;
    TC = 19.173; 
    
    T = TC + 273.15;
    
    // dry air (Harisson 1965)
    dens_N2 = 780840e-6 * LoschmidtConstant*1e-6;
    dens_O2 = 209844e-6 * LoschmidtConstant*1e-6; // assuming CO2 -> O2
    dens_Ar = 9340e-6 * LoschmidtConstant*1e-6;
    dens_Ne = 18.18e-6 * LoschmidtConstant*1e-6;
    dens_He = 5.24e-6 * LoschmidtConstant*1e-6;
    dens_CH4 = 1.774e-6 * LoschmidtConstant*1e-6;
    dens_Kr = 1.14e-6 * LoschmidtConstant*1e-6;
    dens_H2 = 0.56e-6 * LoschmidtConstant*1e-6;
    dens_O3 = 0.0;
    dens_N = 0.0;
    dens_O = 0.0;
    dens_H = 0.0;
    dens_O2 -= CO2ppm*1e-6;
    denstot = dens_N2 + dens_O2 + dens_Ar + dens_Ne + dens_He + dens_CH4 + dens_Kr + dens_H2 + dens_O3 + dens_N + dens_O + dens_H;
    dens_H2O = denstot * (Pw/P) / (1.0 - (Pw/P) - 1e-6*CO2ppm);
    dens_CO2 = denstot * CO2ppm*1e-6 / (1.0 - (Pw/P) - 1e-6*CO2ppm);
    
    denstot = dens_N2 + dens_O2 + dens_Ar + dens_Ne + dens_He + dens_CH4 + dens_Kr + dens_H2 + dens_O3 + dens_N + dens_O + dens_H + dens_H2O + dens_CO2;
    xN2 = dens_N2/denstot;
    xO2 = dens_O2/denstot;
    xAr = dens_Ar/denstot;
    xNe = dens_Ne/denstot;
    xHe = dens_H2/denstot;
    xCH4 = dens_CH4/denstot;
    xKr = dens_Kr/denstot;
    xH2 = dens_H2/denstot;
    xO3 = dens_O3/denstot;
    xN = dens_N/denstot;
    xO = dens_O/denstot;
    xH = dens_H/denstot;
    
    xCO2 = dens_CO2/denstot;
    xH2O = dens_H2O/denstot;
    
    xtot = xN2 + xO2 + xAr + xNe + xHe + xCH4 + xKr + xH2 + xO3 + xN + xO + xH + xCO2 + xH2O;
    xN2 /= xtot;
    xO2 /= xtot;
    xAr /= xtot;
    xNe /= xtot;
    xHe /= xtot;
    xCH4 /= xtot;
    xKr /= xtot;
    xH2 /= xtot;
    xO3 /= xtot;
    xN /= xtot;
    xO /= xtot;
    xH /= xtot;
    xCO2 /= xtot;
    xH2O /= xtot;
 
    printf("CO2 ppm = %.3f\n", xCO2*1e6);
    printf("H2O ppm = %.3f\n", xH2O*1e6);
    Z_Air(P, TC, 27.0); // /Z_N2(101325.0, 273.15);
    //LL *= rhocoeff;
    printf("rhocoeff = %f %f\n", 1.0/rhocoeff, P/101325.0 * (273.15/T));
    
    denstot = P/101325.0 * (273.15/T) * LoschmidtConstant*1e-6;
    denstot = 1.0/rhocoeff * LoschmidtConstant*1e-6;
    
    dens_N2 = xN2 * denstot;
    dens_O2 = xO2 * denstot;
    dens_Ar = xAr * denstot;
    dens_Ne = xNe * denstot;
    dens_He = xHe * denstot;
    dens_CH4 = xCH4 * denstot;
    dens_Kr = xKr * denstot;
    dens_H2 = xH2 * denstot;
    dens_O3 = xO3 * denstot;
    dens_N = xN * denstot;
    dens_O = xO * denstot;
    dens_H = xH * denstot;
    dens_CO2 = xCO2 * denstot;
    dens_H2O = xH2O * denstot;
  
    lambda = 0.633e-6;
    printf("n = %.10g\n", 1.0+AirMixture_N(lambda, dens_N2, dens_O2, dens_Ar, dens_H2O, dens_CO2, dens_Ne, dens_He, dens_CH4, dens_Kr, dens_H2, dens_O3, dens_N, dens_O, dens_H));

  


    fp = fopen("Rlambda.txt", "w");
    for(l=1e-6; l<5e-6; l*=1.0+1e-5)
        fprintf(fp, "%.16f %.16f %.16f %.16f\n", l, AtmosphereModel_stdAtmModel_N(10, l), AtmosphereModel_stdAtmModel_N(1000, l), AtmosphereModel_stdAtmModel_N(4200, l));
    fclose(fp);

    

   fp = fopen("AtmRefrac.txt", "w");
    for(l=0.4e-6; l<2.0e-6; l+=0.01e-6)
        fprintf(fp, "%.16f %.16f\n", l, asin(sin(CONF_ZANGLE)/(1.0+AtmosphereModel_stdAtmModel_N(SiteAlt, l))));
    fclose(fp);

    printf("CONF_ZANGLE = %f  alt = %f\n", CONF_ZANGLE, SiteAlt);

    

   // for(Temp=173.0; Temp<373.0; Temp+=5.0)
     //   printf("T= %lf K    Ps(H2O) [Pa] = %g\n", Temp, AtmosphereModel_H2O_Saturation(Temp));



    //    Scoeff = LAMBDA/SLAMBDA;
    //   Nlambda = 0.0000834213+0.0240603/(130.0-1.0/pow(CONF_LAMBDA*1000000.0,2.0))+0.00015997/(38.9-1.0/pow(CONF_LAMBDA*1000000.0,2.0));
    //   Nslambda = 0.0000834213+0.0240603/(130.0-1.0/pow(CONF_SLAMBDA*1000000.0,2.0))+0.00015997/(38.9-1.0/pow(CONF_SLAMBDA*1000000.0,2.0));

    //printf("method 1 : %f %f\n", Nlambda, Nslambda);

    Nlambda = AtmosphereModel_stdAtmModel_N(0.0, CONF_LAMBDA);
    Nslambda = AtmosphereModel_stdAtmModel_N(0.0, CONF_SLAMBDA);
    Scoeff =  CONF_LAMBDA/CONF_SLAMBDA * Nslambda/Nlambda; // multiplicative coefficient to go from reference lambda phase to science lambda phase


    //  printf("Scoeff is %f (%f)\n",Scoeff,Nslambda/Nlambda);
    //   fflush(stdout);


    // printf("Zenith angle = %f rad\n", CONF_ZANGLE);



    //  pfactor = naxes[0]/CONF_WFsize;
    pfactor = 1;
    contraction_factor = 4;
    if(naxes[0]/CONF_WFsize==1)
        contraction_factor = 0;
    if(naxes[0]/CONF_WFsize==2)
        contraction_factor = 1;
    if(naxes[0]/CONF_WFsize==4)
        contraction_factor = 2;
    if(naxes[0]/CONF_WFsize==8)
        contraction_factor = 3;

    if(contraction_factor==4)
    {
        printf("ERROR: unknown contraction factor\n");
        fflush(stdout);
        exit(0);
    }

    /*  ID=image_ID("ST_pa");*/



    NBFRAMES = (long) (CONF_TIME_SPAN/CONF_WFTIME_STEP);
    naxes[2] = NBFRAMES;
    naxesout[2] = NBFRAMES;

    nelement = naxes[0] * naxes[1] * naxes[2];
    printf("Allocating memory...\n");
    fflush(stdout);

    // OUTPUT ARRAYS
    IDout_array_pha = create_3Dimage_ID("outarraypha",naxesout[0],naxesout[1],naxesout[2]);
    IDout_array_amp = create_3Dimage_ID("outarrayamp",naxesout[0],naxesout[1],naxesout[2]);
    IDout_sarray_pha = create_3Dimage_ID("outsarraypha",naxesout[0],naxesout[1],naxesout[2]);
    IDout_sarray_amp = create_3Dimage_ID("outsarrayamp",naxesout[0],naxesout[1],naxesout[2]);

    //  array_pha = (float*) malloc(naxesout[2]*naxesout[1]*naxesout[0]*sizeof(float));
    //  array_amp = (float*) malloc(naxesout[2]*naxesout[1]*naxesout[0]*sizeof(float));
    //  sarray_pha = (float*) malloc(naxesout[2]*naxesout[1]*naxesout[0]*sizeof(float));
    //  sarray_amp = (float*) malloc(naxesout[2]*naxesout[1]*naxesout[0]*sizeof(float));
    //  printf("Memory allocated\n");
    //  fflush(stdout);
    for(ii=0; ii<naxesout[0]*naxesout[1]*naxesout[2]; ii++)
    {
        data.image[IDout_array_amp].array.F[ii] = 1.0;
        data.image[IDout_array_pha].array.F[ii] = 0.0;

        data.image[IDout_sarray_amp].array.F[ii] = 1.0;
        data.image[IDout_sarray_pha].array.F[ii] = 0.0;
    }



    naxes_MASTER[0] = CONF_MASTER_SIZE;
    naxes_MASTER[1] = CONF_MASTER_SIZE;

    master=0;
    stop=1;
    while(stop)
    {
        sprintf(fname,"t%03ld_%ld",master, CONF_MASTER_SIZE);
        if(!file_exists(fname))
        {
            stop=0;
        }
        else
            master++;
    }

    NBMASTERS = master;
    printf("%ld turbulence master files\n",NBMASTERS);
    fflush(stdout);

    
    
    if((fp=fopen(CONF_TURBULENCE_PROF_FILE,"r"))==NULL)
    {
        printf("Cannot open turbulence profile file \"%s\"\n", CONF_TURBULENCE_PROF_FILE);
        exit(1);
    }
    NBLAYERS=0;
    while(fgets(line,2000,fp)!=NULL)
    {
        sscanf(line,"%s",word);
        if(isdigit(word[0]))
        {
            NBLAYERS+=1;
        }
    }
    fclose(fp);

    LAYER_ALT = (float*) malloc(NBLAYERS*sizeof(float));
    LAYER_CN2 = (float*) malloc(NBLAYERS*sizeof(float));
    LAYER_SPD = (float*) malloc(NBLAYERS*sizeof(float));
    LAYER_DIR = (float*) malloc(NBLAYERS*sizeof(float));
    LAYER_OUTERSCALE = (float*) malloc(NBLAYERS*sizeof(float));
    LAYER_INNERSCALE = (float*) malloc(NBLAYERS*sizeof(float));
    LAYER_SIGMAWSPEED = (float*) malloc(NBLAYERS*sizeof(float));    
    LAYER_LWIND = (float*) malloc(NBLAYERS*sizeof(float));
    
    if((fp=fopen(CONF_TURBULENCE_PROF_FILE, "r"))==NULL)
    {
        printf("Cannot open turbulence profile file \"%s\"\n", CONF_TURBULENCE_PROF_FILE);
        exit(1);
    }
    layer=0;
    while(fgets(line,2000,fp)!=NULL)
    {
        sscanf(line,"%s",word);
        if(isdigit(word[0]))
        {
            sscanf(line,"%f %f %f %f %f %f %f %f", &fl1, &fl2, &fl3, &fl4, &fl5, &fl6, &fl7, &fl8);
            if(fl1>SiteAlt-0.1)
            {
                LAYER_ALT[layer] = fl1;
                LAYER_CN2[layer] = fl2;
                LAYER_SPD[layer] = fl3;
                LAYER_DIR[layer] = fl4;
                LAYER_OUTERSCALE[layer] = fl5;
                LAYER_INNERSCALE[layer] = fl6;
                LAYER_SIGMAWSPEED[layer] = fl7;
                LAYER_LWIND[layer] = fl8;
                layer+=1;
            }
        }
    }
    fclose(fp);

    /* CN2 normalisation for 1024x1024 -> 256x256*/
    /* S<0.7 : x=(S-0.1)/0.06*2
    S>0.7 : x=(S-0.38)/0.033*2 */
    /*  if(seeing<0.7)
    CN2total = (seeing-0.1)/0.03;
    else
    CN2total = (seeing-0.38)/0.0165;*/
    /* S  = sqrt(CN2/16.666)*0.443 = 0.1085 sqrt(CN2) */

    CN2total = 1.0; //84.926*seeing*seeing;

    tmp = 0;
    for(layer=0; layer<NBLAYERS; layer++)
        tmp += LAYER_CN2[layer];
    for(layer=0; layer<NBLAYERS; layer++)
        LAYER_CN2[layer] *= CN2total/tmp;



    for(layer=0; layer<NBLAYERS; layer++)
    {
        printf("Turbulence layer %ld : alt = %f m   CN2 = %f   V = %f m/s   Angle = %f rad   outerscale = %f m    innerscale = %f m  sigmaWind = %f m/s  Lwind = %f m\n", layer, LAYER_ALT[layer], LAYER_CN2[layer], LAYER_SPD[layer], LAYER_DIR[layer], LAYER_OUTERSCALE[layer], LAYER_INNERSCALE[layer], LAYER_SIGMAWSPEED[layer], LAYER_LWIND[layer]);
    }


    


    SLAYER_ALT = (float*) malloc(NBLAYERS*sizeof(float));
    SLAYER_CN2 = (float*) malloc(NBLAYERS*sizeof(float));
    NBSLAYERS = NBLAYERS;
    for(layer=0; layer<NBSLAYERS; layer++)
    {
        SLAYER_ALT[layer] = LAYER_ALT[layer];
        SLAYER_CN2[layer] = LAYER_CN2[layer];
    }


    /// temporary arrays for phase unwrapping
    IDpeakpha_re = create_2Dimage_ID("peakphare", naxesout[0], naxesout[1]);
    IDpeakpha_im = create_2Dimage_ID("peakphaim", naxesout[0], naxesout[1]);
    IDpeakpha = create_2Dimage_ID("peakpha", naxesout[0], naxesout[1]);

    xsizepeakpha = (long) (naxesout[0]/20);
    ysizepeakpha = (long) (naxesout[1]/20);
    IDpeakpha_re_bin = create_2Dimage_ID("peakphare_bin", xsizepeakpha, ysizepeakpha);
    IDpeakpha_im_bin = create_2Dimage_ID("peakphaim_bin", xsizepeakpha, ysizepeakpha);
    IDpeakpha_bin = create_2Dimage_ID("peakpha_bin", xsizepeakpha, ysizepeakpha);
    IDpeakpha_bin_ch = create_2Dimage_ID("peakpha_bin_ch", xsizepeakpha, ysizepeakpha);



    OK=0;
    while(OK==0)
    {
        printf("--------------------\n");
        for(i=0; i<NBSLAYERS; i++)
            printf("Super layer %ld/%ld  alt: %f  CN2: %f\n", i, NBSLAYERS,SLAYER_ALT[i], SLAYER_CN2[i]);

        /* look for minimum altitude difference */
        minaltd = LAYER_ALT[NBLAYERS-1];
        index = 0;
        for(i=0; i<NBSLAYERS-1; i++)
        {
            value = SLAYER_ALT[i+1]-SLAYER_ALT[i];
            if(value<minaltd)
            {
                minaltd=value;
                index = i;
            }
        }
        if((minaltd>CONF_FRESNEL_PROPAGATION_BIN)||(NBSLAYERS==1))
        {
            OK=1;
        }
        else
        {
            /* group SLAYERs i and i+1 */
            printf("Group slayers %ld and %ld\n",index,index+1);
            SLAYER_ALT[index] = (SLAYER_CN2[index]*SLAYER_ALT[index]+SLAYER_CN2[index+1]*SLAYER_ALT[index+1])/(SLAYER_CN2[index]+SLAYER_CN2[index+1]);
            SLAYER_CN2[index] = SLAYER_CN2[index] + SLAYER_CN2[index+1];
            for(i=index+1; i<NBSLAYERS-1; i++)
            {
                SLAYER_ALT[i] = SLAYER_ALT[i+1];
                SLAYER_CN2[i] = SLAYER_CN2[i+1];
            }
            NBSLAYERS-=1;
        }
    }
    for(i=0; i<NBSLAYERS; i++)
    {
        printf("Super layer %ld  alt: %f  CN2: %g\n",i,SLAYER_ALT[i],SLAYER_CN2[i]);
    }

    NB_alt_bin_sep = NBSLAYERS-1;
    alt_bin_sep = (float*) malloc(NB_alt_bin_sep*sizeof(float));
    for(i=0; i<NB_alt_bin_sep; i++)
    {
        alt_bin_sep[i] = 0.5*(SLAYER_ALT[i]+SLAYER_ALT[i+1]);
        printf("Altitude bin separation %ld is %f\n",i,alt_bin_sep[i]);
    }

    free(SLAYER_CN2);

    super_layer_index = (long*) malloc(NBLAYERS*sizeof(long));
    for(layer=0; layer<NBLAYERS; layer++)
    {
        index=0;
        for(i=0; i<NB_alt_bin_sep; i++)
            if((alt_bin_sep[i]<LAYER_ALT[layer])&&(alt_bin_sep[i+1]>LAYER_ALT[layer]))
                index=i+1;
        if(LAYER_ALT[layer]>alt_bin_sep[NB_alt_bin_sep-1])
            index = NB_alt_bin_sep;
        super_layer_index[layer]=index;
        printf("Layer %ld belongs to superlayer %ld/%ld\n",layer,super_layer_index[layer],NBSLAYERS);
    }


    frame=0;
    xpos = (double*) malloc(sizeof(double)*NBLAYERS);
    ypos = (double*) malloc(sizeof(double)*NBLAYERS);
    xpos0 = (double*) malloc(sizeof(double)*NBLAYERS);
    ypos0 = (double*) malloc(sizeof(double)*NBLAYERS);
    xposfcnt = (long*) malloc(sizeof(long)*NBLAYERS);
    yposfcnt = (long*) malloc(sizeof(long)*NBLAYERS);

    vxpix = (float*) malloc(sizeof(float)*NBLAYERS);
    vypix = (float*) malloc(sizeof(float)*NBLAYERS);

    h = 0.0;
    printf("\n\n");
    printf("Refractivity = %g -> %g     %g -> %g\n", CONF_LAMBDA, AtmosphereModel_stdAtmModel_N(h, CONF_LAMBDA), CONF_SLAMBDA, AtmosphereModel_stdAtmModel_N(h, CONF_SLAMBDA));
    

    printf("Computing refraction and position offset for CONF_ZANGLE = %f\n", CONF_ZANGLE);
    for(layer=0; layer<NBLAYERS; layer++)
    {
        xposfcnt[layer] = 0;
        yposfcnt[layer] = 0;

        xpos[layer] = 0.05*CONF_MASTER_SIZE;
        ypos[layer] = 0.05*CONF_MASTER_SIZE;

        // layer shift due to atmospheric refraction
        // computes here what the offset is at a the reference wavelength
        Roffset = 0.0;
        for(h=SiteAlt; h<LAYER_ALT[layer]; h += 1.0)
        {
            Rindex = 1.0 + AtmosphereModel_stdAtmModel_N(h, CONF_LAMBDA);
            tmpf = sin(CONF_ZANGLE)/sqrt(Rindex*Rindex-sin(CONF_ZANGLE)*sin(CONF_ZANGLE));
            tmpf -= sin(CONF_ZANGLE)/sqrt(1.0-sin(CONF_ZANGLE)*sin(CONF_ZANGLE));
            Roffset += tmpf;
 //           printf("h = %12f   Rindex = %12f   Roffset = %12f\n", h, Rindex, Roffset);
        }

        // we compute here the offset at the science wavelength
        RoffsetS = 0.0;
        for(h=SiteAlt; h<LAYER_ALT[layer]; h += 1.0)
        {
            Rindex = 1.0 + AtmosphereModel_stdAtmModel_N(h, CONF_SLAMBDA);
            tmpf = sin(CONF_ZANGLE)/sqrt(Rindex*Rindex-sin(CONF_ZANGLE)*sin(CONF_ZANGLE));
            tmpf -= sin(CONF_ZANGLE)/sqrt(1.0-sin(CONF_ZANGLE)*sin(CONF_ZANGLE));
            RoffsetS += tmpf;
        }

        // the refractive offset is the difference between the reference and science wavelength offsets
        ypos[layer] += (RoffsetS-Roffset)/(CONF_PUPIL_SCALE/pfactor);

        // add here offset due to source position
        // SOURCE_Xpos and SOURCE_Ypos are in radian
        xpos[layer] += CONF_SOURCE_Xpos*LAYER_ALT[layer]/(CONF_PUPIL_SCALE/pfactor);
        ypos[layer] += CONF_SOURCE_Ypos*LAYER_ALT[layer]/(CONF_PUPIL_SCALE/pfactor);

        xpos0[layer] = CONF_SOURCE_Xpos*LAYER_ALT[layer]/(CONF_PUPIL_SCALE/pfactor); // for realtime mode
        ypos0[layer] = CONF_SOURCE_Ypos*LAYER_ALT[layer]/(CONF_PUPIL_SCALE/pfactor);

        vxpix[layer] = LAYER_SPD[layer]*cos(LAYER_DIR[layer])/(CONF_PUPIL_SCALE/pfactor); /* pixel coordinate speed, in pixel per sec, x axis */
        vypix[layer] = LAYER_SPD[layer]*sin(LAYER_DIR[layer])/(CONF_PUPIL_SCALE/pfactor); /* pixel coordinate speed, in pixel per sec, y axis */

        printf("------ layer %5ld, SPEED = %12f x %12f pix/step, offset = %12f m  [ %12f m  %12f m ] ----------\n",layer, vxpix[layer], vypix[layer], RoffsetS-Roffset, RoffsetS, Roffset);
    }



    printf("NBMASTERS = %ld\n",NBMASTERS);
    if(NBMASTERS<NBLAYERS)
        NBMASTERS = NBLAYERS;
    ID_TM = (long*) malloc(sizeof(long)*NBMASTERS);
    for(i=0; i<NBMASTERS; i++)
    {
        sprintf(fname,"t%03ld_%ld", i, CONF_MASTER_SIZE);
        sprintf(fname1,"TM%ld",i);
        if(load_fits(fname, fname1, 1)==-1)
        {
            sprintf(fname,"t%03ld_%ld", i, CONF_MASTER_SIZE);
            sprintf(fname1,"TM%ld", i);
            printf("CREATING %s   (%f - %f)\n", fname, LAYER_OUTERSCALE[i]/CONF_PUPIL_SCALE, LAYER_INNERSCALE[i]/CONF_PUPIL_SCALE);
            make_master_turbulence_screen(fname1, "tursctmp", CONF_MASTER_SIZE, LAYER_OUTERSCALE[i]/CONF_PUPIL_SCALE, LAYER_INNERSCALE[i]/CONF_PUPIL_SCALE);
            save_fl_fits(fname1, fname);
            delete_image_ID("tursctmp");
        }
        ID_TM[i] = image_ID(fname1);
    }
    ID_TML = (long*) malloc(sizeof(long)*NBLAYERS);
    j=0;
    for(i=0; i<NBLAYERS; i++)
    {
        ID_TML[i] = ID_TM[j];
        if(j==NBMASTERS)
        {
            printf("ERROR: number of master turbulence phase screens (%ld) is too small\n",NBMASTERS);
            exit(0);
        }
        j++;
    }







    // Measure r0 (pix) for each master turbulence screen
    dpix = 50;
    r0tot = 0.0;
    r0cnt = 0;
    for(k=0; k<NBMASTERS; k++)
    {
        cnt = 0;
        tot = 0.0;
        for(ii=0; ii<CONF_MASTER_SIZE-dpix; ii++)
            for(jj=0; jj<CONF_MASTER_SIZE; jj++)
            {
                p1 = data.image[ID_TM[k]].array.F[jj*CONF_MASTER_SIZE+ii];
                p2 = data.image[ID_TM[k]].array.F[jj*CONF_MASTER_SIZE+ii+dpix];
                tot += (p1-p2)*(p1-p2);
                cnt++;
            }
        r0 = 1.0*dpix*pow((tot/cnt)/6.88,-3.0/5.0);
        printf("TURBULENCE MASTER %ld    r0 = %g pix\n",k,r0);
        r0tot += r0;
        r0cnt++;
    }
    r0 = r0tot/r0cnt;
    printf("r0 = %g pix -> %g pix\n", r0, CONF_LAMBDA/(CONF_SEEING/3600.0/180.0*PI)/CONF_PUPIL_SCALE*pfactor);

    // renormalize turbulence screens such that a single screen has the right r0
    tmp1 = pow(r0/(CONF_LAMBDA/(CONF_SEEING/3600.0/180.0*PI)/CONF_PUPIL_SCALE*pfactor),5.0/6.0);
    for(k=0; k<NBMASTERS; k++)
        for(ii=0; ii<CONF_MASTER_SIZE*CONF_MASTER_SIZE; ii++)
            data.image[ID_TM[k]].array.F[ii] *= tmp1;

    r0tot = 0.0;
    r0cnt = 0;
    for(k=0; k<NBMASTERS; k++)
    {
        cnt = 0;
        tot = 0.0;
        for(ii=0; ii<CONF_MASTER_SIZE-dpix; ii++)
            for(jj=0; jj<CONF_MASTER_SIZE; jj++)
            {
                p1 = data.image[ID_TM[k]].array.F[jj*CONF_MASTER_SIZE+ii];
                p2 = data.image[ID_TM[k]].array.F[jj*CONF_MASTER_SIZE+ii+dpix];
                tot += (p1-p2)*(p1-p2);
                cnt++;
            }
        r0 = 1.0*dpix*pow((tot/cnt)/6.88,-3.0/5.0);
        printf("TURBULENCE MASTER %ld    r0 = %g pix\n",k,r0);
        r0tot += r0;
        r0cnt++;
    }
    r0 = r0tot/r0cnt;
    printf("r0 = %g pix\n",r0);


    // target seeing = seeing [arcsec]
    // ref lambda = CONF_LAMBDA [m]
    // if single screen:
    // r0[m] = CONF_LAMBDA[m]/seeing[rad]
    // r0[pix] = CONF_LAMBDA*1.0e-6/(seeing/3600.0/180.0*PI)/CONF_PUPIL_SCALE
    // multiply by (r0/r0goal)^6/5


    // each layer coeff mult by sqrt(fracCN2/cos(CONF_ZANGLE)))


    for(i=0; i<NBLAYERS; i++)
    {
        //      coeff = 3.645*183.8115*pow(10.0,-12)/LAMBDA/LAMBDA*CONF_PUPIL_SCALE*PUPIL_SCALE*sqrt(LAYER_CN2[i]);
        //if(pfactor==2)
        //	coeff *= 1.0/(pfactor*pfactor*pfactor*pfactor*pfactor*pfactor);
        for(ii=0; ii<CONF_MASTER_SIZE*CONF_MASTER_SIZE; ii++)
            data.image[ID_TML[i]].array.F[ii] *= sqrt(LAYER_CN2[i]/cos(CONF_ZANGLE));
        printf("Layer %ld, coeff = %g\n",i,sqrt(LAYER_CN2[i]/cos(CONF_ZANGLE)));
    }




    ID_array1 = create_2Dimage_ID("array1",naxes[0],naxes[1]);
    if(CONF_MAKE_SWAVEFRONT==1)
        ID_sarray1 = create_2Dimage_ID("sarray1",naxes[0],naxes[1]);


    /* if(make_cwavefront==1)
         ID_carray1 = create_2Dimage_ID("carray1",naxes[0],naxes[1]);
    */


    if(CONF_WAVEFRONT_AMPLITUDE==1) // includes sub pixel translation
    {
        ID_array2 = create_2DCimage_ID("array2", naxes[0], naxes[1]);
        if(CONF_MAKE_SWAVEFRONT==1)
            ID_sarray2 = create_2DCimage_ID("sarray2", naxes[0], naxes[1]);
        /*     if(make_cwavefront==1)
                 ID_carray2 = create_2DCimage_ID("carray2", naxes[0], naxes[1]);*/
    }


    if((array = (complex_float*) malloc(NBFRAMES*naxes[0]*naxes[1]*sizeof(complex_float)))==NULL)
    {
        printf("Memory allocation error (\"array\" in make_AtmosphericTurbulence_wavefront_series)\n");
        printf("Decrease the size of the wavefront cube\n");
        exit(0);
    }

    if(CONF_MAKE_SWAVEFRONT==1)
    {
        if((sarray = (complex_float*) malloc(NBFRAMES*naxes[0]*naxes[1]*sizeof(complex_float)))==NULL)
        {
            printf("Memory allocation error (\"sarray\" in make_AtmosphericTurbulence_wavefront_series)\n");
            printf("Decrease the size of the wavefront cube\n");
            exit(0);
        }
    }

    /*
        if(make_cwavefront==1)
        {
            if((carray = (complex_float*) malloc(NBFRAMES*naxes[0]*naxes[1]*sizeof(complex_float)))==NULL)
            {
                printf("Memory allocation error (\"carray\" in make_AtmosphericTurbulence_wavefront_series)\n");
                printf("Decrease the size of the wavefront cube\n");
                exit(0);
            }
        }
    */



    if(CONF_SHM_OUTPUT == 1)
    {
        IDshmpha = create_image_ID("atmwfpha", 2, naxesout, FLOAT, 1, 0);
        IDshmamp = create_image_ID("atmwfamp", 2, naxesout, FLOAT, 1, 0);
    }


    if(CONF_SHM_SOUTPUT == 1)
    {
        sprintf(fname, "%spha", CONF_SHM_SPREFIX);
        IDshmspha = create_image_ID(fname, 2, naxesout, FLOAT, 1, 0);
        sprintf(fname, "%samp", CONF_SHM_SPREFIX);
        IDshmsamp = create_image_ID(fname, 2, naxesout, FLOAT, 1, 0);

        kw = 0;
        strcpy(data.image[IDshmspha].kw[kw].name, "TIME");
        data.image[IDshmspha].kw[kw].type = 'D';
        data.image[IDshmspha].kw[kw].value.numf = 0.0;
        strcpy(data.image[IDshmspha].kw[kw].comment, "Physical time [sec]");

        kw = 0;
        strcpy(data.image[IDshmsamp].kw[kw].name, "TIME");
        data.image[IDshmsamp].kw[kw].type = 'D';
        data.image[IDshmsamp].kw[kw].value.numf = 0.0;
        strcpy(data.image[IDshmsamp].kw[kw].comment, "Physical time [sec]");
    }





    printf("SKIP_EXISTING = %d\n", CONF_SKIP_EXISTING);
    if(CONF_SKIP_EXISTING==1)
    {
        start_tspan = 0;
        OK = 1;
        while(OK==1)
        {
            sprintf(fname1,"%s%08ld.pha", CONF_WF_FILE_PREFIX, start_tspan);
            sprintf(fname2,"%s%08ld.pha", CONF_SWF_FILE_PREFIX, start_tspan);
            printf("TESTING FILE %s ... ", fname2);
            if(file_exists(fname2)==1)
            {
                start_tspan ++;
                printf("exists\n");
                OK = 1;
            }
            else
            {
                printf("does not exist\n");
                OK = 0;
            }
        }
    }
    printf("Start TSPAN = %ld\n", start_tspan);


    WSPEEDsize = naxes[0];
    WSPEEDpixscale = CONF_PUPIL_SCALE;
    vKwindsize = 20000.0/WSPEEDpixscale; // 20 km
    
    for(layer=0; layer<NBLAYERS; layer++)
    {
        sprintf(fname, "wspeed_%03ld.fits", layer);
        sprintf(name, "wspeed_%03ld", layer);
        ID = load_fits(fname, name, 1);
        if(ID==-1)
        {
            make_AtmosphericTurbulence_vonKarmanWind(vKwindsize, WSPEEDpixscale, LAYER_SIGMAWSPEED[layer], LAYER_LWIND[layer], WSPEEDsize, name);
            sprintf(fname, "!wspeed_%03ld.fits", layer);
            save_fits(name, fname);
            
            
            // print wind speed as a function of time and position
            sprintf(fname, "wspeed_%03ld.txt", layer);
            ID = image_ID(name);
            if((fp=fopen(fname,"w"))==NULL)
                {
                    printf("ERROR: cannot create file \"%s\"\n", fname);
                    exit(0); 
                }
            for(ii=0;ii<vKwindsize;ii++)
                fprintf(fp, "%6ld   %20f   %20f  %20f  %20f  %20f\n", ii, WSPEEDpixscale*ii, WSPEEDpixscale*ii*LAYER_SPD[layer], data.image[ID].array.F[ii], data.image[ID].array.F[vKwindsize+ii], data.image[ID].array.F[2*vKwindsize+ii]);
            fclose(fp);
        }
    }
    
    
    
    vindex = 0;


    /*
            for(frame=0; frame<NBFRAMES; frame++)
            {
                vindex ++;
                for(layer=NBLAYERS-1; layer!=-1; layer--)
                {
                    // random speed motion
                    vpix = 0.0; //0.01*sin(11.0*vindex*(layer+3))*sqrt(vxpix[layer]*vxpix[layer]+vypix[layer]*vypix[layer]);
                    PA = sin(10.0*vindex*(layer+2));

                    xpos[layer] += vxpix[layer] + vpix*cos(PA);
                    ypos[layer] += vypix[layer] + vpix*sin(PA);

                    while(xpos[layer]<0)
                        xpos[layer] += 1.0*naxes_MASTER[0];
                    while(xpos[layer]>1.0*naxes_MASTER[0])
                        xpos[layer] -= 1.0*naxes_MASTER[0];
                    while(ypos[layer]<0)
                        ypos[layer] += 1.0*naxes_MASTER[1];
                    while(ypos[layer]>1.0*naxes_MASTER[1])
                        ypos[layer] -= 1.0*naxes_MASTER[1];
                }
            }
      */


    printf("WAVEFRONT_AMPLITUDE = %d\n", CONF_WAVEFRONT_AMPLITUDE);
    // exit(0);



    for(tspan=start_tspan; tspan<CONF_NB_TSPAN; tspan++)
    {
        for(frame=0; frame<NBFRAMES; frame++)
        {
            vindex ++;
            if(CONF_MAKE_SWAVEFRONT==1)
            {
                for(ii=0; ii<naxes[0]; ii++)
                    for(jj=0; jj<naxes[1]; jj++)
                    {
                        data.image[ID_array1].array.F[jj*naxes[0]+ii] = 0.0;
                        data.image[ID_sarray1].array.F[jj*naxes[0]+ii] = 0.0;
                    }
                if(CONF_WAVEFRONT_AMPLITUDE==1)
                    for(ii=0; ii<naxes[0]; ii++)
                        for(jj=0; jj<naxes[1]; jj++)
                        {
                            data.image[ID_array2].array.CF[jj*naxes[0]+ii].re = 1.0;
                            data.image[ID_array2].array.CF[jj*naxes[0]+ii].im = 0.0;
                            data.image[ID_sarray2].array.CF[jj*naxes[0]+ii].re = 1.0;
                            data.image[ID_sarray2].array.CF[jj*naxes[0]+ii].im = 0.0;
                        }
            }
            else
            {
                for(ii=0; ii<naxes[0]; ii++)
                    for(jj=0; jj<naxes[1]; jj++)
                        data.image[ID_array1].array.F[jj*naxes[0]+ii] = 0.0;
                if(CONF_WAVEFRONT_AMPLITUDE==1)
                    for(ii=0; ii<naxes[0]; ii++)
                        for(jj=0; jj<naxes[1]; jj++)
                        {
                            data.image[ID_array2].array.CF[jj*naxes[0]+ii].re = 1.0;
                            data.image[ID_array2].array.CF[jj*naxes[0]+ii].im = 0.0;
                        }
            }

            /*       if(make_cwavefront==1)
                   {
                       for(ii=0; ii<naxes[0]; ii++)
                           for(jj=0; jj<naxes[1]; jj++)
                               data.image[ID_carray1].array.F[jj*naxes[0]+ii] = 0.0;
                   }
            */

            usleep(CONF_SIMTDELAY);

            
            if(CONF_WAITFORSEM==1) // wait for semaphore to advance to next WF step
            {
                printf("WAITING for semaphore #0 \"%s\" ...\n", CONF_WAITSEMIMNAME);
                COREMOD_MEMORY_image_set_semwait(CONF_WAITSEMIMNAME, 0);
                printf("Done\n");
            }

            clock_gettime(CLOCK_REALTIME, &tnow);
            tnowdouble = 1.0*tnow.tv_sec + 1.0e-9*tnow.tv_nsec;
            tnowdouble *= CONF_ATMWF_REALTIMEFACTOR;
            for(layer=NBLAYERS-1; layer!=-1; layer--)
            {
                if(CONF_ATMWF_REALTIME==0)
                {
                    tnowdouble = (tspan*NBFRAMES+frame)*CONF_WFTIME_STEP;
                    printf("\rLayer % 2ld/% 2ld, Frame % 4ld/% 4ld, File % 6ld/% 6ld  [TIME = %10.4f s]  ", layer, NBLAYERS, frame, NBFRAMES,tspan, CONF_NB_TSPAN,(tspan*NBFRAMES+frame)*CONF_WFTIME_STEP);
                }
                else
                    printf("\rLayer % 2ld/% 2ld, Frame % 4ld/% 4ld, File % 6ld/% 6ld  [PHYSICAL TIME = %.9lf s]  ",layer,NBLAYERS, frame, NBFRAMES, tspan, CONF_NB_TSPAN, tnowdouble);
                fflush(stdout);

                // recompute Scoeff for this layer
                Nlambda = AtmosphereModel_stdAtmModel_N(LAYER_ALT[layer], CONF_LAMBDA);
                Nslambda = AtmosphereModel_stdAtmModel_N(LAYER_ALT[layer], CONF_SLAMBDA);
                Scoeff =  CONF_LAMBDA/CONF_SLAMBDA * Nslambda/Nlambda; // multiplicative coefficient to go from reference lambda phase to science lambda phase


                if(layer!=NBLAYERS-1)
                {
                    if(super_layer_index[layer+1]!=super_layer_index[layer])
                    {
                        if(CONF_FRESNEL_PROPAGATION==1)
                        {
                            // printf("[%g %g %g]",CONF_PUPIL_SCALE/pfactor,SLAYER_ALT[super_layer_index[layer+1]]-SLAYER_ALT[super_layer_index[layer]], CONF_LAMBDA);
                            Fresnel_propagate_wavefront("array2","array2p", CONF_PUPIL_SCALE/pfactor,(SLAYER_ALT[super_layer_index[layer+1]]-SLAYER_ALT[super_layer_index[layer]])/cos(CONF_ZANGLE), CONF_LAMBDA);
                            delete_image_ID("array2");
                            chname_image_ID("array2p","array2");
                        }
                        ID_array2 = image_ID("array2");
                        if(CONF_MAKE_SWAVEFRONT==1)
                        {
                            if(CONF_FRESNEL_PROPAGATION==1)
                            {
                                Fresnel_propagate_wavefront("sarray2","sarray2p", CONF_PUPIL_SCALE/pfactor, (SLAYER_ALT[super_layer_index[layer+1]]-SLAYER_ALT[super_layer_index[layer]])/cos(CONF_ZANGLE), CONF_SLAMBDA);
                                delete_image_ID("sarray2");
                                chname_image_ID("sarray2p","sarray2");
                            }
                            ID_sarray2 = image_ID("sarray2");
                        }
                    }
                }

                // layer_scale = (SODIUM_ALT-LAYER_ALT[layer])/SODIUM_ALT;

                vpix = 0.0; //0.1*sin(11.0*vindex*(layer+3))*sqrt(vxpix[layer]*vxpix[layer]+vypix[layer]*vypix[layer]);
                PA = sin(10.0*vindex*(layer+2));

                if(CONF_ATMWF_REALTIME==1) // real time
                {
                    xpos[layer] = xpos0[layer] + vxpix[layer]*tnowdouble + 1.0*xposfcnt[layer]*naxes_MASTER[0];
                    ypos[layer] = ypos0[layer] + vypix[layer]*tnowdouble + 1.0*yposfcnt[layer]*naxes_MASTER[0];
                }
                else // non real time
                {
                    
                    xpos[layer] += vxpix[layer]*CONF_WFTIME_STEP;
                    ypos[layer] += vypix[layer]*CONF_WFTIME_STEP;
                }



                xref = (long) (xpos[layer]);
                yref = (long) (ypos[layer]);

                while(xpos[layer]<0)
                {
                    xpos[layer] += 1.0*naxes_MASTER[0];
                    xposfcnt[layer]++;
                    xref = (long) (xpos[layer]);
                }
                while(xpos[layer]>1.0*naxes_MASTER[0])
                {
                    xpos[layer] -= 1.0*naxes_MASTER[0];
                    xposfcnt[layer]--;
                    xref = (long) (xpos[layer]);
                }

                while(ypos[layer]<0)
                {
                    ypos[layer] += 1.0*naxes_MASTER[1];
                    yposfcnt[layer]++;
                    yref = (long) (ypos[layer]);
                }
                while(ypos[layer]>1.0*naxes_MASTER[1])
                {
                    ypos[layer] -= 1.0*naxes_MASTER[1];
                    yposfcnt[layer]--;
                    yref = (long) (ypos[layer]);
                }


                if(xref==naxes_MASTER[0])
                    xref = 0;
                if(yref==naxes_MASTER[1])
                    yref = 0;
                iimax = naxes_MASTER[0]-xref;
                jjmax = naxes_MASTER[1]-yref;
                if(iimax>naxes[0])
                    iimax = naxes[0];
                if(jjmax>naxes[1])
                    jjmax = naxes[1];
                xrefm = xref-naxes_MASTER[0];
                yrefm = yref-naxes_MASTER[1];

                /* make wavefront */

                for(ii=0; ii<naxes[0]; ii++)
                    for(jj=0; jj<naxes[1]; jj++)
                    {
                        iimf = fmod((xpos[layer]+ii),1.0*naxes_MASTER[0]);
                        jjmf = fmod((ypos[layer]+jj),1.0*naxes_MASTER[1]);
                        iim = (long) (iimf);
                        jjm = (long) (jjmf);
                        iifrac = iimf-iim;
                        jjfrac = jjmf-jjm;
                        iim1 = iim+1;
                        jjm1 = jjm+1;
                        if(iim==CONF_MASTER_SIZE)
                            iim = 0;
                        if(jjm==CONF_MASTER_SIZE)
                            jjm = 0;
                        if(iim1>CONF_MASTER_SIZE-1)
                            iim1 -= CONF_MASTER_SIZE;
                        if(jjm1>CONF_MASTER_SIZE-1)
                            jjm1 -= CONF_MASTER_SIZE;

                        //	value = data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim];

                        value = (1.0-iifrac)*(1.0-jjfrac)*data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim];
                        value += (1.0-iifrac)*(jjfrac)*data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim];
                        value += (iifrac)*(jjfrac)*data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim1];
                        value += (iifrac)*(1.0-jjfrac)*data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim1];

                        data.image[ID_array1].array.F[jj*naxes[0]+ii] += value;
                        if(CONF_WAVEFRONT_AMPLITUDE==1)
                        {
                            re = data.image[ID_array2].array.CF[jj*naxes[0]+ii].re;
                            im = data.image[ID_array2].array.CF[jj*naxes[0]+ii].im;
                            data.image[ID_array2].array.CF[jj*naxes[0]+ii].re = re*cos(value)-im*sin(value);
                            data.image[ID_array2].array.CF[jj*naxes[0]+ii].im = re*sin(value)+im*cos(value);
                        }
                    }



                /* make swavefront */
                if(CONF_MAKE_SWAVEFRONT==1)
                {
                    for(ii=0; ii<naxes[0]; ii++)
                        for(jj=0; jj<naxes[1]; jj++)
                        {
                            iimf = fmod((xpos[layer]+ii),1.0*naxes_MASTER[0]);
                            jjmf = fmod((ypos[layer]+jj),1.0*naxes_MASTER[1]);
                            iim = (long) (iimf);
                            jjm = (long) (jjmf);
                            iifrac = iimf-iim;
                            jjfrac = jjmf-jjm;
                            iim1 = iim+1;
                            jjm1 = jjm+1;
                            if(iim==CONF_MASTER_SIZE)
                                iim = 0;
                            if(jjm==CONF_MASTER_SIZE)
                                jjm = 0;
                            if(iim1>CONF_MASTER_SIZE-1)
                                iim1 -= CONF_MASTER_SIZE;
                            if(jjm1>CONF_MASTER_SIZE-1)
                                jjm1 -= CONF_MASTER_SIZE;

                            value = (1.0-iifrac)*(1.0-jjfrac)*data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim];
                            value += (1.0-iifrac)*(jjfrac)*data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim];
                            value += (iifrac)*(jjfrac)*data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim1];
                            value += (iifrac)*(1.0-jjfrac)*data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim1];

                            value *= Scoeff;  // multiplicative coeff to go from ref lambda to science lambda

                            data.image[ID_sarray1].array.F[jj*naxes[0]+ii] += value;

                            if(CONF_WAVEFRONT_AMPLITUDE==1)
                            {
                                re = data.image[ID_sarray2].array.CF[jj*naxes[0]+ii].re;
                                im = data.image[ID_sarray2].array.CF[jj*naxes[0]+ii].im;
                                data.image[ID_sarray2].array.CF[jj*naxes[0]+ii].re = re*cos(value)-im*sin(value);
                                data.image[ID_sarray2].array.CF[jj*naxes[0]+ii].im = re*sin(value)+im*cos(value);
                            }
                        }
                }


                /* cwavefront */
                /*if(make_cwavefront==1)
                {
                    for(ii=0; ii<naxes[0]; ii++)
                        for(jj=0; jj<naxes[1]; jj++)
                        {
                            iim=(long) (fmod((xpos[layer]+naxes[0]/2+(ii-naxes[0]/2)*layer_scale),1.0*naxes_MASTER[0]));
                            jjm=(long) (fmod((ypos[layer]+naxes[1]/2+(jj-naxes[1]/2)*layer_scale),1.0*naxes_MASTER[1]));
                            iim = (long) (iimf);
                            jjm = (long) (jjmf);
                            iifrac = iimf-iim;
                            jjfrac = jjmf-jjm;
                            iim1 = iim+1;
                            jjm1 = jjm+1;
                            if(iim==CONF_MASTER_SIZE)
                                iim = 0;
                            if(jjm==CONF_MASTER_SIZE)
                                jjm = 0;
                            if(iim1>CONF_MASTER_SIZE-1)
                                iim1 -= CONF_MASTER_SIZE;
                            if(jjm1>MASTER_SIZE-1)
                                jjm1 -= CONF_MASTER_SIZE;

                            value = (1.0-iifrac)*(1.0-jjfrac)*data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim];
                            value += (1.0-iifrac)*(jjfrac)*data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim];
                            value += (iifrac)*(jjfrac)*data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim1];
                            value += (iifrac)*(1.0-jjfrac)*data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim1];

                            value *= Scoeff;

                            data.image[ID_carray1].array.F[jj*naxes[0]+ii] += value;

                            if(CONF_WAVEFRONT_AMPLITUDE==1)
                            {
                                re = data.image[ID_carray2].array.CF[jj*naxes[0]+ii].re;
                                im = data.image[ID_carray2].array.CF[jj*naxes[0]+ii].im;
                                data.image[ID_carray2].array.CF[jj*naxes[0]+ii].re = re*cos(value)-im*sin(value);
                                data.image[ID_carray2].array.CF[jj*naxes[0]+ii].im = re*sin(value)+im*cos(value);
                            }
                        }
                }*/

            }



            // REFERENCE LAMBDA
            for(ii=0; ii<naxesout[0]; ii++)
                for(jj=0; jj<naxesout[1]; jj++)
                {
                    ii1 = ii+(naxes[0]-naxesout[0])/2;
                    jj1 = jj+(naxes[1]-naxesout[1])/2;
                    data.image[IDout_array_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii] = data.image[ID_array1].array.F[jj1*naxes[0]+ii1];
                }

            if(CONF_WAVEFRONT_AMPLITUDE==1)
            {
                for(ii=0; ii<naxesout[0]; ii++)
                    for(jj=0; jj<naxesout[1]; jj++)
                    {
                        ii1 = ii+(naxes[0]-naxesout[0])/2;
                        jj1 = jj+(naxes[1]-naxesout[1])/2;
                        array[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].re = data.image[ID_array2].array.CF[jj1*naxes[0]+ii1].re;
                        array[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].im = data.image[ID_array2].array.CF[jj1*naxes[0]+ii1].im;

                        re = array[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].re;
                        im = array[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].im;
                        data.image[IDout_array_amp].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii] = sqrt(re*re+im*im);
                        pha = atan2(im,re);
                        data.image[IDout_array_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii] = pha + 2.0*M_PI*((long) (data.image[IDout_array_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii]/2.0/M_PI+1000.5) - 1000.0);
                    }
            }




            // WRITE CURRENT WF TO SHARED MEMORY

            if(CONF_SHM_OUTPUT == 1)
            {
                if(CONF_WAVEFRONT_AMPLITUDE==0)
                {
                    for(ii=0; ii<naxesout[0]*naxesout[1]; ii++)
                        data.image[IDshmpha].array.F[ii] = data.image[ID_array1].array.F[frame*naxesout[0]*naxesout[1]+ii];
                }
                else
                {
                    for(ii=0; ii<naxesout[0]*naxesout[1]; ii++)
                    {
                        data.image[IDshmpha].array.F[ii] = data.image[IDout_array_pha].array.F[frame*naxesout[0]*naxesout[1]+ii];
                        data.image[IDshmamp].array.F[ii] = data.image[IDout_array_amp].array.F[frame*naxesout[0]*naxesout[1]+ii];
                    }
                }
            }









            // SCIENCE LAMBDA
            if(CONF_MAKE_SWAVEFRONT==1)
            {
                for(ii=0; ii<naxesout[0]; ii++)
                    for(jj=0; jj<naxesout[1]; jj++)
                    {
                        ii1 = ii+(naxes[0]-naxesout[0])/2;
                        jj1 = jj+(naxes[1]-naxesout[1])/2;
                        data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii] = data.image[ID_sarray1].array.F[jj1*naxes[0]+ii1];
                    }
                if(CONF_WAVEFRONT_AMPLITUDE==1)
                {
                    for(ii2=0; ii2<xsizepeakpha*ysizepeakpha; ii2++)
                    {
                        data.image[IDpeakpha_re_bin].array.F[ii2] = 0.0;
                        data.image[IDpeakpha_im_bin].array.F[ii2] = 0.0;
                        data.image[IDpeakpha_bin].array.F[ii2] = 0.0;
                        data.image[IDpeakpha_bin_ch].array.F[ii2] = 0.0;
                    }

                    peakpha_re = 0.0;
                    peakpha_im = 0.0;
                    for(ii=0; ii<naxesout[0]; ii++)
                        for(jj=0; jj<naxesout[1]; jj++)
                        {
                            ii1 = ii+(naxes[0]-naxesout[0])/2;
                            jj1 = jj+(naxes[1]-naxesout[1])/2;
                            sarray[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].re = data.image[ID_sarray2].array.CF[jj1*naxes[0]+ii1].re;
                            sarray[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].im = data.image[ID_sarray2].array.CF[jj1*naxes[0]+ii1].im;
                            re = sarray[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].re;
                            im = sarray[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].im;
                            data.image[IDout_sarray_amp].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii] = sqrt(re*re+im*im);
                            pha = atan2(im,re);
                            data.image[IDpeakpha_re].array.F[jj*naxesout[0]+ii] = cos(data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii]-pha);
                            data.image[IDpeakpha_im].array.F[jj*naxesout[0]+ii] = sin(data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii]-pha);
                            data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii] = pha;
                            ii2 = (long) (1.0*ii/naxesout[0]*xsizepeakpha);
                            jj2 = (long) (1.0*jj/naxesout[1]*ysizepeakpha);

                            if((ii2<xsizepeakpha)&&(jj2<ysizepeakpha))
                            {
                                data.image[IDpeakpha_re_bin].array.F[jj2*xsizepeakpha+ii2] += cos(data.image[ID_sarray1].array.F[jj1*naxes[0]+ii1]-pha);
                                data.image[IDpeakpha_im_bin].array.F[jj2*xsizepeakpha+ii2] += sin(data.image[ID_sarray1].array.F[jj1*naxes[0]+ii1]-pha);
                            }
                            //                                data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii] = pha + 2.0*M_PI*((long) ((data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii]-pha)/2.0/M_PI+1000.5) - 1000.0);
                        }


                    //peakpha = atan2(peakpha_im, peakpha_re);
                    //printf("peak pha = %lf\n", peakpha/2.0/M_PI);


                    peakpha = 0.0;
                    for(ii2=0; ii2<xsizepeakpha*ysizepeakpha; ii2++)
                    {
                        data.image[IDpeakpha_bin].array.F[ii2] = atan2(data.image[IDpeakpha_im_bin].array.F[ii2], data.image[IDpeakpha_re_bin].array.F[ii2]);
                        //	while(data.image[IDpeakpha_bin].array.F[ii2]<0.0)
                        //	data.image[IDpeakpha_bin].array.F[ii2] += 2.0*M_PI;
                    }

                    chcnt = 1;
                    chcnt0cnt = 0;
                    chiter = 0;
                    plim = 1.0;
                    while((plim>0.51)&&(chiter<1000)&&(chcnt0cnt<5))
                    {
                        chiter++;
                        chcnt = 0;
                        for(ii2=0; ii2<xsizepeakpha*ysizepeakpha; ii2++)
                            data.image[IDpeakpha_bin_ch].array.F[ii2] = 0.0;

                        for(ii2=0; ii2<xsizepeakpha-1; ii2++)
                            for(jj2=0; jj2<ysizepeakpha; jj2++)
                            {
                                index1 = jj2*xsizepeakpha+ii2;
                                index2 = jj2*xsizepeakpha+ii2+1;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += 0.2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= 0.2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= 0.2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += 0.2;
                                }
                            }

                        for(ii2=0; ii2<xsizepeakpha; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-1; jj2++)
                            {
                                index1 = jj2*xsizepeakpha+ii2;
                                index2 = (jj2+1)*xsizepeakpha+ii2;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += 0.2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= 0.2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= 0.2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += 0.2;
                                }
                            }


                        for(ii2=0; ii2<xsizepeakpha-1; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-1; jj2++)
                            {
                                index1 = jj2*xsizepeakpha+ii2;
                                index2 = (jj2+1)*xsizepeakpha+ii2+1;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += 0.2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= 0.2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= 0.2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += 0.2;
                                }
                            }


                        for(ii2=0; ii2<xsizepeakpha-1; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-1; jj2++)
                            {
                                index1 = (jj2+1)*xsizepeakpha+ii2;
                                index2 = jj2*xsizepeakpha+ii2+1;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += 0.2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= 0.2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= 0.2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += 0.2;
                                }
                            }


                        pcoeff2;


                        for(ii2=0; ii2<xsizepeakpha-2; ii2++)
                            for(jj2=0; jj2<ysizepeakpha; jj2++)
                            {
                                index1 = jj2*xsizepeakpha+ii2;
                                index2 = jj2*xsizepeakpha+ii2+2;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                }
                            }


                        for(ii2=0; ii2<xsizepeakpha; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-2; jj2++)
                            {
                                index1 = (jj2+2)*xsizepeakpha+ii2;
                                index2 = jj2*xsizepeakpha+ii2;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                }
                            }

                        for(ii2=0; ii2<xsizepeakpha-1; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-2; jj2++)
                            {
                                index1 = (jj2+2)*xsizepeakpha+ii2+1;
                                index2 = jj2*xsizepeakpha+ii2;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                }
                            }

                        for(ii2=0; ii2<xsizepeakpha-1; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-2; jj2++)
                            {
                                index1 = (jj2+2)*xsizepeakpha+ii2;
                                index2 = jj2*xsizepeakpha+ii2+1;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                }
                            }

                        for(ii2=0; ii2<xsizepeakpha-2; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-2; jj2++)
                            {
                                index1 = (jj2+2)*xsizepeakpha+ii2+2;
                                index2 = jj2*xsizepeakpha+ii2;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                }
                            }

                        for(ii2=0; ii2<xsizepeakpha-2; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-2; jj2++)
                            {
                                index1 = (jj2+2)*xsizepeakpha+ii2;
                                index2 = jj2*xsizepeakpha+ii2+2;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                }
                            }

                        for(ii2=0; ii2<xsizepeakpha-2; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-1; jj2++)
                            {
                                index1 = (jj2+1)*xsizepeakpha+ii2;
                                index2 = jj2*xsizepeakpha+ii2+2;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                }
                            }

                        for(ii2=0; ii2<xsizepeakpha-2; ii2++)
                            for(jj2=0; jj2<ysizepeakpha-1; jj2++)
                            {
                                index1 = (jj2+1)*xsizepeakpha+ii2+2;
                                index2 = jj2*xsizepeakpha+ii2;
                                pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                if(pv2>pv1+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                }
                                if(pv1>pv2+M_PI)
                                {
                                    data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                    data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                }
                            }



                        plim = 0.0;
                        for(ii2=0; ii2<xsizepeakpha*ysizepeakpha; ii2++)
                            if(fabs(data.image[IDpeakpha_bin_ch].array.F[ii2])>plim)
                                plim = fabs(data.image[IDpeakpha_bin_ch].array.F[ii2]);
                        plim -= 0.001;

                        if(plim<0.5)
                            plim = 0.5;

                        //	plim = 2.39;

                        //                        save_fits("peakpha_bin", "!peakpha_bin.fits");

                        for(ii2=0; ii2<xsizepeakpha; ii2++)
                            for(jj2=0; jj2<ysizepeakpha; jj2++)
                            {
                                if(data.image[IDpeakpha_bin_ch].array.F[jj2*xsizepeakpha+ii2]>plim)
                                {
                                    if((ii2>1)&&(jj2>1)&&(ii2<xsizepeakpha-2)&&(jj2<ysizepeakpha-2))
                                        chcnt ++;
                                    data.image[IDpeakpha_bin].array.F[jj2*xsizepeakpha+ii2] += 2.0*M_PI;
                                }
                                if(data.image[IDpeakpha_bin_ch].array.F[jj2*xsizepeakpha+ii2]<-plim)
                                {
                                    if((ii2>1)&&(jj2>1)&&(ii2<xsizepeakpha-2)&&(jj2<ysizepeakpha-2))
                                        chcnt ++;
                                    data.image[IDpeakpha_bin].array.F[jj2*xsizepeakpha+ii2] -= 2.0*M_PI;
                                }

                            }
                        //                    printf("chiter = %ld    [%ld]  %f\n", chiter, chcnt, plim);
                        //                      save_fits("peakpha_bin_ch", "!peakpha_bin_ch.fits");

                        if(chcnt==0)
                            chcnt0cnt++;
                        else
                            chcnt0cnt = 0;
                    }


                    /*          list_image_ID();
                              save_fits("peakpha_bin", "!peakpha2_bin.fits");
                              save_fits("peakphare_bin", "!peakpha_re_bin.fits");
                              save_fits("peakphaim_bin", "!peakpha_im_bin.fits");
                    */
                    for(ii=0; ii<naxesout[0]; ii++)
                        for(jj=0; jj<naxesout[1]; jj++)
                        {
                            ii2 = (long) (1.0*ii/naxesout[0]*xsizepeakpha);
                            jj2 = (long) (1.0*jj/naxesout[1]*ysizepeakpha);

                            if((ii2<xsizepeakpha)&&(jj2<ysizepeakpha))
                                peakpha = data.image[IDpeakpha_bin].array.F[jj2*xsizepeakpha+ii2];

                            ii1 = ii+(naxes[0]-naxesout[0])/2;
                            jj1 = jj+(naxes[1]-naxesout[1])/2;

                            pha = data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii];
                            data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii] = pha + 2.0*M_PI*((long) ((data.image[ID_sarray1].array.F[jj1*naxes[0]+ii1]-pha-peakpha)/2.0/M_PI+1000.5) - 1000.0);
                        }
                }
            }






            /*
                        if(make_cwavefront==1)
                        {
                            for(ii=0; ii<naxesout[0]; ii++)
                                for(jj=0; jj<naxesout[1]; jj++)
                                {
                                    ii1 = ii+(naxes[0]-naxesout[0])/2;
                                    jj1 = jj+(naxes[1]-naxesout[1])/2;
                                    carray[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].re = data.image[ID_carray2].array.CF[jj1*naxes[0]+ii1].re;
                                    carray[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].im = data.image[ID_carray2].array.CF[jj1*naxes[0]+ii1].im;
                                }
                        }

            */


            // WRITE CURRENT WF TO SHARED MEMORY
            if((CONF_SHM_SOUTPUT == 1)&&(CONF_MAKE_SWAVEFRONT==1))
            {
                switch (CONF_SHM_SOUTPUTM) {
                case 1 :
                    coeff = CONF_SLAMBDA/2.0/M_PI;
                    break;
                case 2 :
                    coeff = CONF_SLAMBDA/2.0/M_PI*1e6;
                    break;
                default :
                    coeff = 1.0;
                    break;
                }



                if(CONF_WAVEFRONT_AMPLITUDE==0)
                {
                    data.image[IDshmspha].md[0].write = 1;
                    data.image[IDshmspha].kw[0].value.numf = tnowdouble;
                    for(ii=0; ii<naxesout[0]*naxesout[1]; ii++)
                        data.image[IDshmspha].array.F[ii] = data.image[ID_sarray1].array.F[frame*naxesout[0]*naxesout[1]+ii]*coeff;
                    data.image[IDshmspha].md[0].cnt0++;
                    data.image[IDshmspha].md[0].write = 0;
                }
                else
                {
                    data.image[IDshmspha].md[0].write = 1;
                    data.image[IDshmsamp].md[0].write = 1;
                    data.image[IDshmspha].kw[0].value.numf = tnowdouble;
                    data.image[IDshmsamp].kw[0].value.numf = tnowdouble;
                    for(ii=0; ii<naxesout[0]*naxesout[1]; ii++)
                    {
                        data.image[IDshmspha].array.F[ii] = data.image[IDout_sarray_pha].array.F[frame*naxesout[0]*naxesout[1]+ii]*coeff;
                        data.image[IDshmsamp].array.F[ii] = data.image[IDout_sarray_amp].array.F[frame*naxesout[0]*naxesout[1]+ii];
                    }
                    data.image[IDshmspha].md[0].cnt0++;
                    data.image[IDshmsamp].md[0].cnt0++;
                    data.image[IDshmspha].md[0].write = 0;
                    data.image[IDshmsamp].md[0].write = 0;
                }

            }
        }




        if(CONF_WFOUTPUT==1) // WRITE REFERENCE LAMBDA
        {
            sprintf(fname1,"!%s%08ld.pha", CONF_WF_FILE_PREFIX,tspan);
            sprintf(fname2,"!%s%08ld.amp", CONF_WF_FILE_PREFIX,tspan);

            save_fl_fits("outarraypha",fname1);


            if(CONF_WAVEFRONT_AMPLITUDE==1)
                save_fl_fits("outarrayamp",fname2);
        }
        else if (CONF_WFOUTPUT==0)
        {
            // CREATE EMPTY FILES
            sprintf(fname1,"%s%08ld.pha", CONF_WF_FILE_PREFIX,tspan);
            sprintf(fname2,"%s%08ld.amp", CONF_WF_FILE_PREFIX,tspan);
            sprintf(command,"touch %s",fname1);
            r = system(command);
            if(CONF_WAVEFRONT_AMPLITUDE==1)
            {
                sprintf(command,"touch %s",fname2);
                r = system(command);
            }
        }


        if((CONF_MAKE_SWAVEFRONT==1)&&(CONF_SWF_WRITE2DISK==1)) // WRITE SCIENCE LAMBDA
        {
            printf("WRITING SCIENCE WAVEFRONT ...");
            fflush(stdout);
            sprintf(fname1,"!%s%08ld.pha", CONF_SWF_FILE_PREFIX, tspan);
            sprintf(fname2,"!%s%08ld.amp", CONF_SWF_FILE_PREFIX, tspan);


            save_fl_fits("outsarraypha",fname1);

            printf(" - ");
            fflush(stdout);

            if(CONF_WAVEFRONT_AMPLITUDE==1)
                save_fl_fits("outsarrayamp",fname2);

            printf("\n");
            fflush(stdout);
        }

    }

    delete_image_ID("array1");
    free(array);
    if(CONF_MAKE_SWAVEFRONT==1)
    {
        delete_image_ID("sarray1");
        free(sarray);
    }
    /*  if(make_cwavefront==1)
      {
          delete_image_ID("carray1");
          free(carray);
      }*/

    free(SLAYER_ALT);
    free(super_layer_index);
    free(xpos);
    free(ypos);
    free(xpos0);
    free(ypos0);
    free(xposfcnt);
    free(yposfcnt);
    free(vxpix);
    free(vypix);


    free(LAYER_ALT);
    free(LAYER_CN2);
    free(LAYER_SPD);
    free(LAYER_DIR);
    free(LAYER_OUTERSCALE);
    free(LAYER_INNERSCALE);
    free(LAYER_SIGMAWSPEED);
    free(LAYER_LWIND);
    free(ID_TM);
    free(ID_TML);

    free(naxes);
    free(naxesout);

    return(0);
}



















int contract_wavefront_series(char *in_prefix, char *out_prefix, long NB_files)
{
    /* contracts the wavefront series by a factor of 2 */
    char fname[200];
    long IDamp,IDpha,IDoutamp,IDoutpha;
    long ii,jj,kk;
    long i,j;
    long naxes[3];
    long naxes_out[3];
    float re,im,amp,pha;
    long index;
    float P;
    long LARGE = 10000;
    float pharef,ampref;

    for(index=0; index<NB_files; index++)
    {
        printf("INDEX = %ld/%ld\n",index,NB_files);
        sprintf(fname,"%s%08ld.pha",in_prefix,index);
        load_fits(fname, "tmpwfp", 1);
        IDpha=image_ID("tmpwfp");
        sprintf(fname,"%s%08ld.amp",in_prefix,index);
        load_fits(fname, "tmpwfa", 1);
        IDamp=image_ID("tmpwfa");
        naxes[0] = data.image[IDpha].md[0].size[0];
        naxes[1] = data.image[IDpha].md[0].size[1];
        naxes[2] = data.image[IDpha].md[0].size[2];
        naxes_out[0] = data.image[IDpha].md[0].size[0]/2;
        naxes_out[1] = data.image[IDpha].md[0].size[1]/2;
        naxes_out[2] = data.image[IDpha].md[0].size[2];
        IDoutpha = create_3Dimage_ID("tmpwfop",naxes_out[0],naxes_out[1],naxes_out[2]);
        IDoutamp = create_3Dimage_ID("tmpwfoa",naxes_out[0],naxes_out[1],naxes_out[2]);

        ii=0;
        jj=0;
        kk=0;
        amp = 0.0;
        pha = 0.0;
        for(kk=0; kk<naxes[2]; kk++)
        {
            for(ii=0; ii<naxes[0]/2; ii++)
                for(jj=0; jj<naxes[1]/2; jj++)
                {
                    re=0.0;
                    im=0.0;
                    pharef = 0.0;
                    ampref = 0.0;
                    for(i=0; i<2; i++)
                        for(j=0; j<2; j++)
                        {
                            amp = data.image[IDamp].array.F[kk*naxes[0]*naxes[1]+(2*jj+j)*naxes[0]+2*ii+i];
                            pha = data.image[IDpha].array.F[kk*naxes[0]*naxes[1]+(2*jj+j)*naxes[0]+2*ii+i];
                            pharef += data.image[IDamp].array.F[kk*naxes[0]*naxes[1]+(2*jj+j)*naxes[0]+2*ii+i]*data.image[IDpha].array.F[kk*naxes[0]*naxes[1]+(2*jj+j)*naxes[0]+2*ii+i];
                            ampref += data.image[IDamp].array.F[kk*naxes[0]*naxes[1]+(2*jj+j)*naxes[0]+2*ii+i];
                            re += amp*cos(pha);
                            im += amp*sin(pha);
                        }
                    amp = sqrt(re*re+im*im);
                    pha = atan2(im,re);
                    pharef /= ampref;
                    P = 2.0*PI*( ((long) (0.5+1.0*LARGE+(pharef-pha)/2.0/PI)) - LARGE);
                    if(ampref<0.01)
                        P = 0.0;
                    data.image[IDoutpha].array.F[kk*naxes_out[0]*naxes_out[1]+jj*naxes_out[0]+ii] = pha+P;
                    data.image[IDoutamp].array.F[kk*naxes_out[0]*naxes_out[1]+jj*naxes_out[0]+ii] = amp/4.0;
                }
        }
        sprintf(fname,"%s%8ld.pha",out_prefix,index);
        replace_char(fname,' ','0');
        save_fl_fits("tmpwfop",fname);
        sprintf(fname,"%s%8ld.amp",out_prefix,index);
        replace_char(fname,' ','0');
        save_fl_fits("tmpwfoa",fname);

        delete_image_ID("tmpwfa");
        delete_image_ID("tmpwfp");
        delete_image_ID("tmpwfoa");
        delete_image_ID("tmpwfop");
    }

    return(0);
}



//
// analyze WF series: PSF FWHM and aperture photometry
//

int measure_wavefront_series(float factor)
{
    float FOCAL_SCALE;
    double tmp;
    long ID,IDpsf,IDamp,IDpha;
    long IDpupamp;
    long tspan;

    long ID_array1;
    float amp,pha;

    char fnameamp[200];
    char fnamepha[200];

    long naxes[3];
    long frame;
    long NBFRAMES;
    long ii,jj;
    int amplitude_on;

    double puprad = 0.035; // meter
    double pupradpix;
    double psfflux;
    double psfflux1; // within 1 arcsec radius
    double psfflux2; // within 2 arcsec radius
    double psfflux5; // within 5 arcsec radius
    double psfflux10; // within 10 arcsec radius
    double psfflux20; // within 20 arcsec radius

    double dx, dy, r;
    FILE *fpphot;


    AtmosphericTurbulence_ReadConf();

    
    pupradpix = puprad/CONF_PUPIL_SCALE;
    printf("pupradpix = %f m\n",pupradpix);


    FOCAL_SCALE = CONF_LAMBDA/CONF_WFsize/CONF_PUPIL_SCALE/PI*180.0*3600.0; /* in arcsecond per pixel */
    printf("Scale is %f arcsecond per pixel (%ld pixels)\n", FOCAL_SCALE, CONF_WFsize);
    IDpupamp = image_ID("ST_pa");
    naxes[0]=data.image[ID].md[0].size[0];
    naxes[1]=data.image[ID].md[0].size[1];





    NBFRAMES = (long) (CONF_TIME_SPAN/CONF_WFTIME_STEP);
    naxes[2]=NBFRAMES;

    ID_array1 = create_2DCimage_ID("array1", naxes[0], naxes[1]);

    IDpsf = create_2Dimage_ID("PSF", naxes[0], naxes[1]);

    fpphot = fopen("phot.txt","w");
    fclose(fpphot);

    for(tspan=0; tspan<CONF_NB_TSPAN; tspan++)
    {
        printf("%ld/%ld\n", tspan, CONF_NB_TSPAN);
        sprintf(fnamepha,"%s%8ld.pha", CONF_WF_FILE_PREFIX, tspan);
        replace_char(fnamepha,' ','0');
        sprintf(fnameamp,"%s%8ld.amp", CONF_WF_FILE_PREFIX, tspan);
        replace_char(fnameamp,' ','0');
        IDpha = load_fits(fnamepha, "wfpha", 1);
        if(amplitude_on==1)
            IDamp = load_fits(fnameamp, "wfamp", 1);

        for(frame=0; frame<NBFRAMES; frame++)
        {
            psfflux = 0.0;
            psfflux1 = 0.0;
            psfflux2 = 0.0;
            psfflux5 = 0.0;
            psfflux10 = 0.0;
            psfflux20 = 0.0;
            if(amplitude_on==1)
                for(ii=0; ii<naxes[0]; ii++)
                    for(jj=0; jj<naxes[1]; jj++)
                    {
                        amp = data.image[IDamp].array.F[frame*naxes[0]*naxes[1]+jj*naxes[0]+ii]*data.image[IDpupamp].array.F[jj*naxes[0]+ii];
                        pha = factor*data.image[IDpha].array.F[frame*naxes[0]*naxes[1]+jj*naxes[0]+ii];
                        data.image[ID_array1].array.CF[jj*naxes[0]+ii].re = amp*cos(pha);
                        data.image[ID_array1].array.CF[jj*naxes[0]+ii].im = amp*sin(pha);
                        psfflux += amp*amp;
                    }
            else
                for(ii=0; ii<naxes[0]; ii++)
                    for(jj=0; jj<naxes[1]; jj++)
                    {
                        amp = data.image[IDpupamp].array.F[jj*naxes[0]+ii];
                        pha = factor*data.image[IDpha].array.F[frame*naxes[0]*naxes[1]+jj*naxes[0]+ii];
                        data.image[ID_array1].array.CF[jj*naxes[0]+ii].re = amp*cos(pha);
                        data.image[ID_array1].array.CF[jj*naxes[0]+ii].im = amp*sin(pha);
                    }

            do2dfft("array1","im_c");
            permut("im_c");
            ID=image_ID("im_c");
            for(ii=0; ii<naxes[0]; ii++)
                for(jj=0; jj<naxes[1]; jj++)
                {
                    dx = 1.0*ii-naxes[0]/2;
                    dy = 1.0*jj-naxes[1]/2;
                    r = sqrt(dx*dx+dy*dy);
                    tmp = (data.image[ID].array.CF[jj*naxes[0]+ii].re*data.image[ID].array.CF[jj*naxes[0]+ii].re+data.image[ID].array.CF[jj*naxes[0]+ii].im*data.image[ID].array.CF[jj*naxes[0]+ii].im);
                    data.image[IDpsf].array.F[jj*naxes[0]+ii] += tmp;
                    if(r<1.0/FOCAL_SCALE)
                        psfflux1 += tmp;
                    if(r<2.0/FOCAL_SCALE)
                        psfflux2 += tmp;
                    if(r<5.0/FOCAL_SCALE)
                        psfflux5 += tmp;
                    if(r<10.0/FOCAL_SCALE)
                        psfflux10 += tmp;
                    if(r<20.0/FOCAL_SCALE)
                        psfflux20 += tmp;
                }
            delete_image_ID("im_c");
            printf("%.6f %.4f %.4f %.4f %.4f %.4f %.4f\n", (tspan*NBFRAMES+frame)*CONF_WFTIME_STEP, psfflux, psfflux1, psfflux2, psfflux5, psfflux10, psfflux20);
            fpphot = fopen("phot.txt","a");
            fprintf(fpphot,"%.6f %.4f %.4f %.4f %.4f %.4f %.4f\n", (tspan*NBFRAMES+frame)*CONF_WFTIME_STEP, psfflux, psfflux1, psfflux2, psfflux5, psfflux10, psfflux20);
            fclose(fpphot);
        }

        delete_image_ID("wfamp");
        delete_image_ID("wfpha");
    }
    delete_image_ID("array1");
    save_fl_fits("PSF","!PSF.fits");
    tmp = measure_FWHM("PSF",1.0*naxes[0]/2,1.0*naxes[1]/2,1.0,naxes[0]/2);
    printf("FWHM = %f arcseconds (%f pixels)\n",FOCAL_SCALE*tmp,tmp);

    return(0);
}






int measure_wavefront_series_expoframes(float etime, char *outfile)
{
    float FOCAL_SCALE;
    float tmp,tmp1;
    long ID,IDpsf,IDamp,IDpha;
    long tspan;
    FILE *fp;
    char command[200];
    float frac = 0.5;

    long ID_array1;
    float amp,pha;

    char fnameamp[200];
    char fnamepha[200];
    char fname[200];

    long naxes[3];
    long frame;
    long NBFRAMES;
    long ii,jj;
    float etime_left;
    long frame_number;
    double *xcenter;
    double *ycenter;
    long amplitude_on;

    xcenter = (double*) malloc(sizeof(double));
    ycenter = (double*) malloc(sizeof(double));

    printf("Frame exposure time is %f s\n",etime);

 AtmosphericTurbulence_ReadConf();




    FOCAL_SCALE = CONF_LAMBDA/CONF_WFsize/CONF_PUPIL_SCALE/PI*180.0*3600.0; /* in arcsecond per pixel */
    printf("Scale is %f arcsecond per pixel (%ld pixels)\n",FOCAL_SCALE, CONF_WFsize);

    ID=image_ID("ST_pa");
    naxes[0]=data.image[ID].md[0].size[0];
    naxes[1]=data.image[ID].md[0].size[1];
    IDamp = ID;

    NBFRAMES = (long) (CONF_TIME_SPAN/CONF_WFTIME_STEP);
    naxes[2]=NBFRAMES;

    ID_array1 = create_2DCimage_ID("array1", naxes[0], naxes[1]);

    IDpsf = create_2Dimage_ID("PSF",naxes[0],naxes[1]);
    frame_number = 0;
    etime_left = etime;

    sprintf(command,"rm -rf %s",outfile);
    if(system(command)==-1)
    {
        printf("ERROR: system(\"%s\"), %s line %d\n",command,__FILE__,__LINE__);
        exit(0);
    }

    if((fp=fopen(outfile,"w"))==NULL)
    {
        printf("Cannot create file %s\n",outfile);
        exit(0);
    }
    fclose(fp);


    for(tspan=0; tspan<CONF_NB_TSPAN; tspan++)
    {
        printf("%ld/%ld\n", tspan, CONF_NB_TSPAN);
        sprintf(fnamepha,"%s%8ld.pha", CONF_WF_FILE_PREFIX, tspan);
        replace_char(fnamepha,' ','0');
        sprintf(fnameamp,"%s%8ld.amp", CONF_WF_FILE_PREFIX, tspan);
        replace_char(fnameamp,' ','0');
        IDpha = load_fits(fnamepha, "wfpha", 1);
        if(amplitude_on==1)
        {
            printf("reading amp\n");
            fflush(stdout);
            IDamp = load_fits(fnameamp, "wfamp", 1);
        }

        for(frame=0; frame<NBFRAMES; frame++)
        {
            if(amplitude_on==1)
                for(ii=0; ii<naxes[0]; ii++)
                    for(jj=0; jj<naxes[1]; jj++)
                    {
                        amp = data.image[IDamp].array.F[frame*naxes[0]*naxes[1]+jj*naxes[0]+ii];
                        pha = data.image[IDpha].array.F[frame*naxes[0]*naxes[1]+jj*naxes[0]+ii];
                        data.image[ID_array1].array.CF[jj*naxes[0]+ii].re = amp*cos(pha);
                        data.image[ID_array1].array.CF[jj*naxes[0]+ii].im = amp*sin(pha);
                    }
            else
                for(ii=0; ii<naxes[0]; ii++)
                    for(jj=0; jj<naxes[1]; jj++)
                    {
                        amp = data.image[IDamp].array.F[jj*naxes[0]+ii];
                        pha = data.image[IDpha].array.F[frame*naxes[0]*naxes[1]+jj*naxes[0]+ii];
                        data.image[ID_array1].array.CF[jj*naxes[0]+ii].re = amp*cos(pha);
                        data.image[ID_array1].array.CF[jj*naxes[0]+ii].im = amp*sin(pha);
                    }

            do2dfft("array1","im_c");
            permut("im_c");
            ID=image_ID("im_c");
            for(ii=0; ii<naxes[0]; ii++)
                for(jj=0; jj<naxes[1]; jj++)
                {
                    data.image[IDpsf].array.F[jj*naxes[0]+ii] += (data.image[ID].array.CF[jj*naxes[0]+ii].re*data.image[ID].array.CF[jj*naxes[0]+ii].re+data.image[ID].array.CF[jj*naxes[0]+ii].im*data.image[ID].array.CF[jj*naxes[0]+ii].im);
                }
            delete_image_ID("im_c");

            etime_left -= CONF_WFTIME_STEP;
            if(etime_left<0)
            {
                sprintf(fname,"!PSF%ld",frame_number);
                save_fl_fits("PSF",fname);
                gauss_filter("PSF","PSFg",3,10);
                center_PSF("PSFg", xcenter, ycenter, naxes[0]/2);
                tmp = measure_FWHM("PSFg",xcenter[0],ycenter[0],1.0,naxes[0]/2);
                printf("%ld FWHM %f arcseconds (%f pixels)\n",frame_number,FOCAL_SCALE*tmp,tmp);
                tmp1 = measure_enc_NRJ("PSF",xcenter[0],ycenter[0],frac);
                printf("Encircled energy (%f) is %f arcseconds\n",frac,tmp1*FOCAL_SCALE);
                delete_image_ID("PSFg");
                if((fp=fopen(outfile,"a"))==NULL)
                {
                    printf("Cannot open file %s\n",outfile);
                    exit(0);
                }
                fprintf(fp,"%ld %f %f %f %f\n",frame_number,FOCAL_SCALE*tmp,2.0*FOCAL_SCALE*tmp1,xcenter[0],ycenter[0]);
                fclose(fp);

                arith_image_zero("PSF");
                etime_left = etime;
                frame_number++;
                printf("Working on frame %ld\n",frame_number);
            }

        }
        delete_image_ID("wfamp");
        delete_image_ID("wfpha");
    }
    free(xcenter);
    free(ycenter);
    delete_image_ID("array1");
    save_fl_fits("PSF","!PSF");
    tmp = measure_FWHM("PSF",1.0*naxes[0]/2,1.0*naxes[1]/2,1.0,naxes[0]/2);
    printf("FWHM = %f arcseconds (%f pixels)\n",FOCAL_SCALE*tmp,tmp);

    return(0);
}


int frame_select_PSF(char *logfile, long NBfiles, float frac)
{
    /* logfile has the following format:
       <PSF file name> <FWHM> <Enc.ener.0.5> <centerx> <centery>
    */
    /* outputs the following files :
       PSFc:    coadded PSF (no centering/filtering)
       PSFcc:   coadded PSF with centering
       PSFccsf: coadded PSF with centering and selection on FWHM (frac= fraction of the frames kept)
       PSFccse: coadded PSF with centering and selection on Enc.ener.(frac= fraction of the frames kept)
    */
    FILE *fp;
    int OK;
    long i;
    float *FWHM;
    float *ENCE;
    float *xcen;
    float *ycen;
    char fname[200];
    float Xcenter,Ycenter;
    float Xcenter_c,Ycenter_c;
    float Xcenter_cc,Ycenter_cc;
    float Xcenter_ccsf,Ycenter_ccsf;
    float Xcenter_ccse,Ycenter_ccse;
    float limit;
    long cnt;
    float fwhm1,fwhm2,fwhm3,fwhm4;
    float ence1,ence2,ence3,ence4;
    float fs = 0.01128;

    Xcenter_c = 128;
    Ycenter_c = 128;

    FWHM = (float*) malloc(sizeof(float)*NBfiles);
    ENCE = (float*) malloc(sizeof(float)*NBfiles);
    xcen = (float*) malloc(sizeof(float)*NBfiles);
    ycen = (float*) malloc(sizeof(float)*NBfiles);

    /* make PSFc */
    if((fp=fopen(logfile,"r"))==NULL)
    {
        printf("ERROR: cannot open file \"%s\"\n",logfile);
        exit(0);
    }
    for(i=0; i<NBfiles; i++)
    {
        if(fscanf(fp,"%s %f %f %f %f\n",fname,&FWHM[i],&ENCE[i],&xcen[i],&ycen[i])!=5)
        {
            printf("ERROR: fscanf, %s line %d\n",__FILE__,__LINE__);
            exit(0);
        }

        if(i==0)
            load_fits(fname, "PSFc", 1);
        else
        {
            load_fits(fname, "tmppsf", 1);
            execute_arith("PSFc=PSFc+tmppsf");
            delete_image_ID("tmppsf");
        }
    }
    fclose(fp);

    /* make PSFcc */
    Xcenter = 0.0;
    Ycenter = 0.0;
    if((fp=fopen(logfile,"r"))==NULL)
    {
        printf("ERROR: cannot open file \"%s\"\n",logfile);
        exit(0);
    }
    for(i=0; i<NBfiles; i++)
    {
        if(fscanf(fp,"%s %f %f %f %f\n",fname,&FWHM[i],&ENCE[i],&xcen[i],&ycen[i])!=5)
        {
            printf("ERROR: fscanf, %s line %d\n",__FILE__,__LINE__);
            exit(0);
        }

        if(i==0)
        {
            load_fits(fname, "PSFcc", 1);
            Xcenter = xcen[0];
            Ycenter = ycen[0];
        }
        else
        {
            load_fits(fname, "tmppsf" , 1);
            basic_add("PSFcc","tmppsf","nPSFcc",Xcenter-xcen[i],Ycenter-ycen[i]);
            if(Xcenter<xcen[i])
                Xcenter = xcen[i];
            if(Ycenter<ycen[i])
                Ycenter = ycen[i];
            delete_image_ID("PSFcc");
            delete_image_ID("tmppsf");
            chname_image_ID("nPSFcc","PSFcc");
        }
    }
    fclose(fp);
    Xcenter_cc = Xcenter;
    Ycenter_cc = Ycenter;

    /* make PSFccsf */
    quick_sort_float(FWHM, NBfiles);
    limit = FWHM[(long) (frac*NBfiles)];
    Xcenter = 0.0;
    Ycenter = 0.0;
    if((fp=fopen(logfile,"r"))==NULL)
    {
        printf("ERROR: cannot open file \"%s\"\n",logfile);
        exit(0);
    }
    OK = 0;
    cnt = 0;
    for(i=0; i<NBfiles; i++)
    {
        if(fscanf(fp,"%s %f %f %f %f\n",fname,&FWHM[i],&ENCE[i],&xcen[i],&ycen[i])!=5)
        {
            printf("ERROR: fscanf, %s line %d\n",__FILE__,__LINE__);
            exit(0);
        }

        if(FWHM[i]<limit)
        {
            cnt++;
            if(OK==0)
            {
                load_fits(fname, "PSFccsf", 1);
                Xcenter = xcen[i];
                Ycenter = ycen[i];
                OK=1;
            }
            else
            {
                load_fits(fname, "tmppsf", 1);
                basic_add("PSFccsf","tmppsf","nPSFccsf",Xcenter-xcen[i],Ycenter-ycen[i]);
                if(Xcenter<xcen[i])
                    Xcenter = xcen[i];
                if(Ycenter<ycen[i])
                    Ycenter = ycen[i];
                delete_image_ID("PSFccsf");
                delete_image_ID("tmppsf");
                chname_image_ID("nPSFccsf","PSFccsf");
            }
        }
    }
    fclose(fp);
    printf("PSFccsf: %ld frames kept\n",cnt);
    Xcenter_ccsf = Xcenter;
    Ycenter_ccsf = Ycenter;

    /* make PSFccse */
    quick_sort_float(ENCE, NBfiles);
    limit = ENCE[(long) (frac*NBfiles)];
    Xcenter = 0.0;
    Ycenter = 0.0;
    if((fp=fopen(logfile,"r"))==NULL)
    {
        printf("ERROR: cannot open file \"%s\"\n",logfile);
        exit(0);
    }
    OK = 0;
    cnt = 0;
    for(i=0; i<NBfiles; i++)
    {
        if(fscanf(fp,"%s %f %f %f %f\n",fname,&FWHM[i],&ENCE[i],&xcen[i],&ycen[i])!=5)
        {
            printf("ERROR: fscanf, %s line %d\n",__FILE__,__LINE__);
            exit(0);
        }
        if(ENCE[i]<limit)
        {
            cnt++;
            if(OK==0)
            {
                load_fits(fname, "PSFccse", 1);
                Xcenter = xcen[i];
                Ycenter = ycen[i];
                OK=1;
            }
            else
            {
                load_fits(fname, "tmppsf", 1);
                basic_add("PSFccse","tmppsf","nPSFccse",Xcenter-xcen[i],Ycenter-ycen[i]);
                if(Xcenter<xcen[i])
                    Xcenter = xcen[i];
                if(Ycenter<ycen[i])
                    Ycenter = ycen[i];
                delete_image_ID("PSFccse");
                delete_image_ID("tmppsf");
                chname_image_ID("nPSFccse","PSFccse");
            }
        }
    }
    fclose(fp);
    printf("PSFccse : %ld frames kept\n",cnt);
    Xcenter_ccse = Xcenter;
    Ycenter_ccse = Ycenter;


    /* quality evaluation */
    fwhm1 = measure_FWHM("PSFc",Xcenter_c,Ycenter_c,1.0,128);
    fwhm2 = measure_FWHM("PSFcc",Xcenter_cc,Ycenter_cc,1.0,128);
    fwhm3 = measure_FWHM("PSFccsf",Xcenter_ccsf,Ycenter_ccsf,1.0,128);
    fwhm4 = measure_FWHM("PSFccse",Xcenter_ccse,Ycenter_ccse,1.0,128);
    ence1 = measure_enc_NRJ("PSFc",Xcenter_c,Ycenter_c,0.5);
    ence2 = measure_enc_NRJ("PSFcc",Xcenter_cc,Ycenter_cc,0.5);
    ence3 = measure_enc_NRJ("PSFccsf",Xcenter_ccsf,Ycenter_ccsf,0.5);
    ence4 = measure_enc_NRJ("PSFccse",Xcenter_ccse,Ycenter_ccse,0.5);

    printf("PSFc   :  %f %f\n",fwhm1*fs,2.0*ence1*fs);
    printf("PSFcc  :  %f %f\n",fwhm2*fs,2.0*ence2*fs);
    printf("PSFccsf:  %f %f\n",fwhm3*fs,2.0*ence3*fs);
    printf("PSFccse:  %f %f\n",fwhm4*fs,2.0*ence4*fs);

    free(FWHM);
    free(ENCE);
    free(xcen);
    free(ycen);

    return(0);
}




// explore long exposure PSF structure
double AtmosphericTurbulence_makePSF(double Kp, double Ki, double Kd, double Kdgain)
{
    FILE *fp;
    char wf_file_name[200];
    double PupScale = 0.01;



    // SIMULATION PARAMETERS

    long WFSLAMBDA = 600; // [nm]
    double zeroptWFS = 8.354e10; // [ph/s/um/m2] (600 nm)
    //double zeroptWFS = 9.444e9; // [ph/s/um/m2] (1600 nm)


    long SCILAMBDA = 1600; // [nm]
    double zeroptSCI = 9.444e9; // [ph/s/um/m2] (H band)



    double TelDiam = 30.0;
    double etime = 10.0; // end of exposure [s]
    double etimestart = 0.01; // loop closing delay before start of exposure [s]
    double dtime = 0.00001; // internal time step, 0.1 ms
    double WFS_SamplingTime = 0.00025; // WFS sampling time [sec]
    double WFS_Delay = 0.0002; // [sec]  DELAY SHOULD BE SMALLER THAN SAMPLING TIME

    double bandpassWFS = 0.1; // [um]
    double bandpassSCI = 0.1; // [um]
    double throughputWFS = 0.1;
    double throughputSCI = 0.1;
    double magnWFS = 6.0;
    double magnSCI = 6.0;
    double FLUXSCI; // [ph/s]
    double FLUXWFS; // [ph/s]

    double rtime; // running time

    // input parameters
    float FOCAL_SCALE;

    long BINWF = 2; // 2 for 30m telescope
    long WFsize1; // after binning by BINWF

    long size = 512;
    long BINFACTOR = 8; // 8 for 8m telescope
    long sizeb;

    char WFDIRECTORY[200];
    char fnameamp[200];
    char fnamepha[200];
    char fname[200];

    long NBFRAMES;
    long CubeSize;
    long frame0, frame1, frameindex0, frameindex1, cubeindex0, cubeindex1;
    long frame;
    long ii, jj, ii1, jj1;
    long IDac0, IDac1, IDpc0, IDpc1;
    long IDacs0, IDacs1, IDpcs0, IDpcs1;
    char imname[200];
    double framefrac, framef;
    double amp, pha, val0, val1;
    double wfstime, wfstime1;
    long wfscnt;

    int WFSdelayWait;

    long IDtelpup;
    long IDatm_amp, IDatm_opd;
    long IDatm_amp_sci, IDatm_opd_sci;
    long IDwfs_opd, IDwfs_amp, IDwfs_mes_opd, IDwfs_mes_opd_prev;
    long IDwfs_mes_opd_derivative, IDwfs_mes_opd_integral;
    long IDdm_opd, IDdm_opd_tmp;
    long IDsci_amp, IDsci_opd;
    double tot;

    double RMSwf;

    long IDpupa, IDpupp;
    long ID, IDpsfcumul, ID1, IDpsfcumul1, IDre, IDim;
    double peak, re, im, reave, imave, tot0;
    double re_sci, im_sci, errpha, pha_ave, pha_ave_sci;
    double tmpd;
    long i, j;
    long cnt;

    long IDtmp;


    int PIDok = 0;
    double value = 0.0;
    long valuecnt = 0;



    AtmosphericTurbulence_ReadConf();

    FLUXSCI = zeroptSCI*M_PI*TelDiam*TelDiam/4.0*bandpassSCI*throughputSCI/pow(2.511886,magnSCI);
    printf("FLUX SCI = %g ph/s\n", FLUXSCI);
    FLUXWFS = zeroptWFS*M_PI*TelDiam*TelDiam/4.0*bandpassWFS*throughputWFS/pow(2.511886,magnWFS);
    printf("FLUX WFS = %g ph/s\n", FLUXWFS);

    sprintf(WFDIRECTORY,"/media/data1/WFsim/AtmSim_0.01");
    sprintf(CONFFILE, "%s/WF%04ld/AOsim.conf", WFDIRECTORY, WFSLAMBDA);
    AtmosphericTurbulence_change_configuration_file(CONFFILE);

    printf("Frame exposure time is %f s\n", etime);

    
    printf("PUPIL SCALE = %f m\n", CONF_PUPIL_SCALE);
    printf("PUPIL DIAM = %f pix\n", TelDiam/CONF_PUPIL_SCALE);


    WFsize1 = CONF_WFsize/BINWF;

    FOCAL_SCALE = CONF_SLAMBDA/CONF_WFsize/CONF_PUPIL_SCALE/PI*180.0*3600.0/BINWF; /* in arcsecond per pixel */
    printf("Scale is %f arcsecond per pixel (%ld pixels)\n", FOCAL_SCALE, CONF_WFsize);



    CubeSize = (long) (CONF_TIME_SPAN/CONF_WFTIME_STEP-0.001);
    NBFRAMES = (long) (1.0*etime/CONF_WFTIME_STEP);

    printf("CUBE SIZE = %ld\n", CubeSize);


    // MAKE PUPIL MASK
    IDtelpup = make_disk("TelPup", WFsize1, WFsize1, WFsize1/2, WFsize1/2, TelDiam*0.5/CONF_PUPIL_SCALE/BINWF);


    IDatm_opd = create_2Dimage_ID("atmopd", WFsize1, WFsize1); // Atmosphere OPD at WFS lambda
    IDatm_amp = create_2Dimage_ID("atmamp", WFsize1, WFsize1); // Atmosphere amplitude at WFS lambda

    IDatm_opd_sci = create_2Dimage_ID("atmopdsci", WFsize1, WFsize1); // Atmosphere OPD at SCI lambda
    IDatm_amp_sci = create_2Dimage_ID("atmampsci", WFsize1, WFsize1); // Atmosphere amplitude at SCI lambda

    IDwfs_opd = create_2Dimage_ID("wfsopd", WFsize1, WFsize1); // WFS OPD (corrected by DM)
    IDwfs_amp = create_2Dimage_ID("wfsamp", WFsize1, WFsize1); // WFS amplitude (corrected by DM)

    IDdm_opd_tmp = create_2Dimage_ID("dmopdtmp", WFsize1, WFsize1); // DM opd (before being sent to DM)
    IDdm_opd = create_2Dimage_ID("dmopd", WFsize1, WFsize1); // DM opd

    IDwfs_mes_opd = create_2Dimage_ID("wfsmesopd", WFsize1, WFsize1); // measured WFS OPD
    IDwfs_mes_opd_prev = create_2Dimage_ID("wfsmesopdprev", WFsize1, WFsize1); // previous measured WFS OPD
    IDwfs_mes_opd_derivative = create_2Dimage_ID("wfsmesopdder", WFsize1, WFsize1); // derivative of measured WFS OPD
    IDwfs_mes_opd_integral = create_2Dimage_ID("wfsmesopdint", WFsize1, WFsize1); // integral of measured WFS OPD

    IDsci_opd = create_2Dimage_ID("sciopd", WFsize1, WFsize1); // SCI OPD (corrected by DM)
    IDsci_amp = create_2Dimage_ID("sciamp", WFsize1, WFsize1); // SCI amplitude (corrected by DM)

    sizeb = size*BINFACTOR;
    IDpupa = create_2Dimage_ID("pupa", sizeb, sizeb);
    IDpupp = create_2Dimage_ID("pupp", sizeb, sizeb);
    IDpsfcumul = create_2Dimage_ID("PSFcumul", size, size);
    IDpsfcumul1 = create_2Dimage_ID("PSFcumul1", size, size);


    frame = 0;
    frame = 0;
    rtime = 0.0;
    framefrac = 0.0;
    framef = 0.0;
    cubeindex0 = 0;
    cubeindex1 = 0;
    frameindex0 = 0;
    frameindex1 = 1;

    /*  sprintf(fnamepha, "%s/WF0800/WF4096/wf_%08ld.pha", WFDIRECTORY, cubeindex1);
    sprintf(fnameamp, "%s/WF0800/WF4096/wf_%08ld.amp", WFDIRECTORY, cubeindex1);
    ID_WFphaC = load_fits(fnamepha, "WFphaC");
    ID_WFampC = load_fits(fnameamp, "WFampC");
    */

    wfstime = 0.0;
    for(ii=0; ii<WFsize1*WFsize1; ii++)
    {
        data.image[IDwfs_opd].array.F[ii] = 0.0;
        data.image[IDdm_opd].array.F[ii] = 0.0;
        data.image[IDwfs_mes_opd].array.F[ii] = 0.0;
        data.image[IDwfs_mes_opd_prev].array.F[ii] = 0.0;
        data.image[IDwfs_mes_opd_derivative].array.F[ii] = 0.0;
        data.image[IDwfs_mes_opd_integral].array.F[ii] = 0.0;
    }
    wfscnt = 0;
    WFSdelayWait = 1;

    fp = fopen("result.log", "w");
    fclose(fp);

    while (rtime < etime)
    {
        framef = rtime/CONF_WFTIME_STEP;
        frame = (long) framef;
        framefrac = framef-frame;

        // compute current OPD map for WFS
        frame0 = frame;
        frame1 = frame+1;
        frameindex0 = frame0-cubeindex0*CubeSize;
        frameindex1 = frame1-cubeindex1*CubeSize;

        printf("time = %2.10g s  (%04ld %04ld) (%04ld %04ld) %1.6f\n", rtime, frameindex0, cubeindex0, frameindex1, cubeindex1, framefrac);

        while(frameindex0>CubeSize-1)
        {
            cubeindex0 ++;
            frameindex0 -= CubeSize;
        }

        while(frameindex1>CubeSize-1)
        {
            cubeindex1 ++;
            frameindex1 -= CubeSize;
        }


        sprintf(imname, "wfa%08ld",cubeindex0);
        IDac0 = image_ID(imname);
        if(IDac0 == -1)
        {
            sprintf(fname, "%s/WF%04ld/WF4096/wf_%08ld.amp", WFDIRECTORY, WFSLAMBDA, cubeindex0);
            printf("LOADING %s\n", fname);
            IDac0 = load_fits(fname, imname, 1);
        }
        sprintf(imname, "wfa%08ld",cubeindex1);
        IDac1 = image_ID(imname);
        if(IDac1 == -1)
        {
            sprintf(fname, "%s/WF%04ld/WF4096/wf_%08ld.amp", WFDIRECTORY, WFSLAMBDA, cubeindex1);
            printf("LOADING %s\n", fname);
            IDac1 = load_fits(fname, imname, 1);
        }

        sprintf(imname, "wfp%08ld",cubeindex0);
        IDpc0 = image_ID(imname);
        if(IDpc0 == -1)
        {
            sprintf(fname, "%s/WF%04ld/WF4096/wf_%08ld.pha", WFDIRECTORY, WFSLAMBDA, cubeindex0);
            printf("LOADING %s\n", fname);
            IDpc0 = load_fits(fname, imname, 1);
        }
        sprintf(imname, "wfp%08ld",cubeindex1);
        IDpc1 = image_ID(imname);
        if(IDpc1 == -1)
        {
            sprintf(fname, "%s/WF%04ld/WF4096/wf_%08ld.pha", WFDIRECTORY, WFSLAMBDA, cubeindex1);
            printf("LOADING %s\n", fname);
            IDpc1 = load_fits(fname, imname, 1);
        }




        sprintf(imname, "swfa%08ld",cubeindex0);
        IDacs0 = image_ID(imname);
        if(IDacs0 == -1)
        {
            sprintf(fname, "%s/WF%04ld/WF4096/wf_%08ld.amp", WFDIRECTORY, SCILAMBDA, cubeindex0);
            printf("LOADING %s\n", fname);
            IDacs0 = load_fits(fname, imname, 1);
        }
        sprintf(imname, "swfa%08ld",cubeindex1);
        IDacs1 = image_ID(imname);
        if(IDacs1 == -1)
        {
            sprintf(fname, "%s/WF%04ld/WF4096/wf_%08ld.amp", WFDIRECTORY, SCILAMBDA, cubeindex1);
            printf("LOADING %s\n", fname);
            IDacs1 = load_fits(fname, imname, 1);
        }

        sprintf(imname, "swfp%08ld",cubeindex0);
        IDpcs0 = image_ID(imname);
        if(IDpcs0 == -1)
        {
            sprintf(fname, "%s/WF%04ld/WF4096/wf_%08ld.pha", WFDIRECTORY, SCILAMBDA, cubeindex0);
            printf("LOADING %s\n", fname);
            IDpcs0 = load_fits(fname, imname, 1);
        }
        sprintf(imname, "swfp%08ld",cubeindex1);
        IDpcs1 = image_ID(imname);
        if(IDpcs1 == -1)
        {
            sprintf(fname, "%s/WF%04ld/WF4096/wf_%08ld.pha", WFDIRECTORY, SCILAMBDA, cubeindex1);
            printf("LOADING %s\n", fname);
            IDpcs1 = load_fits(fname, imname, 1);
        }



        if(cubeindex0>0)
        {
            sprintf(imname, "wfa%08ld",cubeindex0-1);
            IDtmp = image_ID(imname);
            if(IDtmp!=-1)
                delete_image_ID(imname);

            sprintf(imname, "wfp%08ld",cubeindex0-1);
            IDtmp = image_ID(imname);
            if(IDtmp!=-1)
                delete_image_ID(imname);

            sprintf(imname, "swfa%08ld",cubeindex0-1);
            IDtmp = image_ID(imname);
            if(IDtmp!=-1)
                delete_image_ID(imname);

            sprintf(imname, "swfp%08ld",cubeindex0-1);
            IDtmp = image_ID(imname);
            if(IDtmp!=-1)
                delete_image_ID(imname);
        }


        if(BINWF==1)
        {
            for(ii=0; ii<CONF_WFsize*CONF_WFsize; ii++)
            {
                val0 = data.image[IDac0].array.F[frameindex0*CONF_WFsize*CONF_WFsize+ii];
                val1 = data.image[IDac1].array.F[frameindex1*CONF_WFsize*CONF_WFsize+ii];
                amp = (1.0-framefrac)*val0 + framefrac*val1;

                val0 = data.image[IDpc0].array.F[frameindex0*CONF_WFsize*CONF_WFsize+ii];
                val1 = data.image[IDpc1].array.F[frameindex1*CONF_WFsize*CONF_WFsize+ii];
                pha = (1.0-framefrac)*val0 + framefrac*val1;


                data.image[IDatm_amp].array.F[ii] = amp;
                data.image[IDatm_opd].array.F[ii] = pha/2.0/M_PI*CONF_SLAMBDA;


                val0 = data.image[IDacs0].array.F[frameindex0*CONF_WFsize*CONF_WFsize+ii];
                val1 = data.image[IDacs1].array.F[frameindex1*CONF_WFsize*CONF_WFsize+ii];
                amp = (1.0-framefrac)*val0 + framefrac*val1;

                val0 = data.image[IDpcs0].array.F[frameindex0*CONF_WFsize*CONF_WFsize+ii];
                val1 = data.image[IDpcs1].array.F[frameindex1*CONF_WFsize*CONF_WFsize+ii];
                pha = (1.0-framefrac)*val0 + framefrac*val1;

                data.image[IDatm_amp_sci].array.F[ii] = amp;
                data.image[IDatm_opd_sci].array.F[ii] = pha/2.0/M_PI*(1.0e-9*SCILAMBDA);
            }
        }
        else
        {
            for(ii1=0; ii1<WFsize1; ii1++)
                for(jj1=0; jj1<WFsize1; jj1++)
                {
                    re = 0.0;
                    im = 0.0;
                    re_sci = 0.0;
                    im_sci = 0.0;
                    pha_ave = 0.0;
                    pha_ave_sci = 0.0;

                    for(i=0; i<BINWF; i++)
                        for(j=0; j<BINWF; j++)
                        {
                            ii = ii1*BINWF+i;
                            jj = jj1*BINWF+j;

                            val0 = data.image[IDac0].array.F[frameindex0*CONF_WFsize*CONF_WFsize+jj*CONF_WFsize+ii];
                            val1 = data.image[IDac1].array.F[frameindex1*CONF_WFsize*CONF_WFsize+jj*CONF_WFsize+ii];
                            amp = (1.0-framefrac)*val0 + framefrac*val1;

                            val0 = data.image[IDpc0].array.F[frameindex0*CONF_WFsize*CONF_WFsize+jj*CONF_WFsize+ii];
                            val1 = data.image[IDpc1].array.F[frameindex1*CONF_WFsize*CONF_WFsize+jj*CONF_WFsize+ii];
                            pha = (1.0-framefrac)*val0 + framefrac*val1;

                            re += amp*cos(pha);
                            im += amp*sin(pha);
                            pha_ave += pha;


                            val0 = data.image[IDacs0].array.F[frameindex0*CONF_WFsize*CONF_WFsize+jj*CONF_WFsize+ii];
                            val1 = data.image[IDacs1].array.F[frameindex1*CONF_WFsize*CONF_WFsize+jj*CONF_WFsize+ii];
                            amp = (1.0-framefrac)*val0 + framefrac*val1;

                            val0 = data.image[IDpcs0].array.F[frameindex0*CONF_WFsize*CONF_WFsize+jj*CONF_WFsize+ii];
                            val1 = data.image[IDpcs1].array.F[frameindex1*CONF_WFsize*CONF_WFsize+jj*CONF_WFsize+ii];
                            pha = (1.0-framefrac)*val0 + framefrac*val1;

                            re_sci += amp*cos(pha);
                            im_sci += amp*sin(pha);
                            pha_ave_sci += pha;
                        }
                    re /= (BINWF*BINWF);
                    im /= (BINWF*BINWF);
                    re_sci /= (BINWF*BINWF);
                    im_sci /= (BINWF*BINWF);
                    pha_ave /= (BINWF*BINWF);
                    pha_ave_sci /= (BINWF*BINWF);

                    data.image[IDatm_amp].array.F[jj1*WFsize1+ii1] = sqrt(re*re+im*im);
                    pha = pha_ave;
                    errpha = atan2(im,re)-pha_ave; // close to -2PI, 0, 2PI etc...
                    errpha = modf(errpha/(2.0*M_PI),&tmpd);
                    if(errpha>0.5)
                        errpha -= 1.0;
                    if(errpha<-0.5)
                        errpha += 1.0;
                    pha += errpha*2.0*M_PI;
                    data.image[IDatm_opd].array.F[jj1*WFsize1+ii1] = pha/2.0/M_PI*CONF_SLAMBDA;


                    data.image[IDatm_amp_sci].array.F[jj1*WFsize1+ii1] = sqrt(re_sci*re_sci+im_sci*im_sci);
                    pha = pha_ave_sci;
                    errpha = atan2(im_sci,re_sci)-pha_ave_sci; // close to -2PI, 0, 2PI etc...
                    errpha = modf(errpha/(2.0*M_PI),&tmpd);
                    if(errpha>0.5)
                        errpha -= 1.0;
                    if(errpha<-0.5)
                        errpha += 1.0;
                    pha += errpha*2.0*M_PI;
                    data.image[IDatm_opd_sci].array.F[jj1*WFsize1+ii1] = pha/2.0/M_PI*(1.0e-9*SCILAMBDA);
                }
        }

        //    save_fl_fits("atmopdsci","!test_atmopdsci.fits");
        //save_fl_fits("atmampsci","!test_atmampsci.fits");


        // Apply DM
        for(ii=0; ii<WFsize1*WFsize1; ii++)
        {
            data.image[IDwfs_opd].array.F[ii] = data.image[IDatm_opd].array.F[ii] - data.image[IDdm_opd].array.F[ii];
            data.image[IDsci_opd].array.F[ii] = data.image[IDatm_opd_sci].array.F[ii] - data.image[IDdm_opd].array.F[ii];
        }


        // WFS integration
        for(ii=0; ii<WFsize1*WFsize1; ii++)
            data.image[IDwfs_mes_opd].array.F[ii] += data.image[IDwfs_opd].array.F[ii];
        wfscnt++;



        if(wfstime>WFS_SamplingTime) // WFS measurement
        {

            // ADD WFS NOISE HERE
            //
            //

            if(PIDok == 1)
            {
                for(ii=0; ii<WFsize1*WFsize1; ii++)
                {
                    data.image[IDwfs_mes_opd].array.F[ii] /= wfscnt; // AVERAGE OVER WFS INTEGRATION TIME
                    data.image[IDwfs_mes_opd_derivative].array.F[ii] = (1.0-Kdgain)*data.image[IDwfs_mes_opd_derivative].array.F[ii] + Kdgain*(data.image[IDwfs_mes_opd].array.F[ii]-data.image[IDwfs_mes_opd_prev].array.F[ii])/WFS_SamplingTime;
                    data.image[IDwfs_mes_opd_integral].array.F[ii] += WFS_SamplingTime*data.image[IDwfs_mes_opd].array.F[ii];
                }
            }
            else
            {
                for(ii=0; ii<WFsize1*WFsize1; ii++)
                {
                    data.image[IDwfs_mes_opd].array.F[ii] /= wfscnt; // AVERAGE OVER WFS INTEGRATION TIME
                    data.image[IDwfs_mes_opd_derivative].array.F[ii] = 0.0;
                    data.image[IDwfs_mes_opd_integral].array.F[ii] += WFS_SamplingTime*data.image[IDwfs_mes_opd].array.F[ii];
                }
                PIDok = 1;
            }




            for(ii=0; ii<WFsize1*WFsize1; ii++)
            {
                data.image[IDdm_opd_tmp].array.F[ii] += Kp*data.image[IDwfs_mes_opd].array.F[ii] + Kd*data.image[IDwfs_mes_opd_derivative].array.F[ii]*WFS_SamplingTime;
                data.image[IDwfs_mes_opd].array.F[ii] = 0.0;
            }
            wfscnt = 0.0;
            wfstime = 0.0;
            WFSdelayWait = 1; // start wait
            wfstime1 = 0.0;



            for(ii=0; ii<WFsize1*WFsize1; ii++)
                data.image[IDwfs_mes_opd_prev].array.F[ii] = data.image[IDwfs_mes_opd].array.F[ii]; // PREVIOUS WFS MEASUREMENT
        }

        if((wfstime1>WFS_Delay)&&(WFSdelayWait==1))
        {
            printf("UPDATE DM SHAPE\n");
            WFSdelayWait = 0;
            for(ii=0; ii<WFsize1*WFsize1; ii++)
                data.image[IDdm_opd].array.F[ii] = data.image[IDdm_opd_tmp].array.F[ii];
        }

        // MEASURE WF QUALITY
        val0 = 0.0;
        val1 = 0.0;
        for(ii=0; ii<WFsize1*WFsize1; ii++)
        {
            val0 += data.image[IDwfs_opd].array.F[ii]*data.image[IDtelpup].array.F[ii];
            val1 += data.image[IDtelpup].array.F[ii];
        }
        for(ii=0; ii<WFsize1*WFsize1; ii++) // REMOVE PISTON
            data.image[IDwfs_opd].array.F[ii] -= val0/val1;

        val0 = 0.0;
        val1 = 0.0;
        for(ii=0; ii<WFsize1*WFsize1; ii++) // COMPUTE RMS WF QUALITY
        {
            val0 += data.image[IDsci_opd].array.F[ii]*data.image[IDsci_opd].array.F[ii]*data.image[IDtelpup].array.F[ii];
            val1 += data.image[IDtelpup].array.F[ii];
        }
        RMSwf = sqrt(val0/val1);

        //    tp("0.0");

        printf("WAVEFRONT QUALITY = %g m\n", RMSwf);
        fp = fopen("result.log", "a");
        fprintf(fp,"%10.10g %g\n", rtime, RMSwf);
        fclose(fp);

        if(rtime>etimestart)
        {
            value += RMSwf;
            valuecnt ++;
        }


        // MAKING FOCAL PLANE IMAGE
        //printf("%ld %ld\n", sizeb, WFsize1);
        // list_image_ID();
        //    printf("Image identifiers: %ld %ld %ld %ld %ld\n", IDpupa, IDpupp, IDtelpup, IDatm_amp_sci, IDsci_opd);
        fflush(stdout);

        tot = 0.0;
        for(ii=0; ii<sizeb; ii++)
            for(jj=0; jj<sizeb; jj++)
            {
                ii1 = WFsize1/2-sizeb/2+ii;
                jj1 = WFsize1/2-sizeb/2+jj;
                if((ii1>-1)&&(jj1>-1)&&(ii1<WFsize1)&&(jj1<WFsize1))
                {
                    data.image[IDpupa].array.F[jj*sizeb+ii] = data.image[IDtelpup].array.F[jj1*WFsize1+ii1]*data.image[IDatm_amp_sci].array.F[jj1*WFsize1+ii1];
                    tot += data.image[IDpupa].array.F[jj*sizeb+ii];
                    data.image[IDpupp].array.F[jj*sizeb+ii] = data.image[IDtelpup].array.F[jj1*WFsize1+ii1]*2.0*M_PI*data.image[IDsci_opd].array.F[jj1*WFsize1+ii1]/(1.0e-9*SCILAMBDA);
                }
            }

        //      tp("0.5");

        mk_reim_from_amph("pupa", "pupp", "pupre", "pupim", 0);
        basic_contract("pupre","pupre1",BINFACTOR,BINFACTOR);
        basic_contract("pupim","pupim1",BINFACTOR,BINFACTOR);
        delete_image_ID("pupre");
        delete_image_ID("pupim");
        mk_complex_from_reim("pupre1", "pupim1", "pupc", 0);
        permut("pupc");
        do2dfft("pupc","focc");
        permut("focc");
        delete_image_ID("pupc");
        mk_amph_from_complex("focc", "foca", "focp", 0);
        delete_image_ID("focc");
        delete_image_ID("focp");
        execute_arith("foci=foca*foca");
        delete_image_ID("foca");

        ID = image_ID("foci");
        IDpsfcumul = image_ID("PSFcumul");


        //      tp("1.0");


        // MAKING CORONAGRAPHIC FOCAL PLANE IMAGE
        // USING SIMPLE CORONAGRAPH MODEL REMOVING PERFECTLY MODE 0
        IDre = image_ID("pupre1");
        IDim = image_ID("pupim1");
        peak = 0.0;
        for(ii=0; ii<size*size; ii++)
        {
            re = data.image[IDre].array.F[ii];
            im = data.image[IDim].array.F[ii];
            amp = re*re+im*im;
            if(amp>peak)
                peak = amp;
        }
        reave = 0.0;
        imave = 0.0;
        cnt = 0;
        for(ii=0; ii<size*size; ii++)
        {
            re = data.image[IDre].array.F[ii];
            im = data.image[IDim].array.F[ii];
            reave += re;
            imave += im;
            amp = re*re+im*im;
            if(amp>0.2*peak)
                cnt ++;
        }
        reave /= cnt;
        imave /= cnt;
        for(ii=0; ii<size*size; ii++)
        {
            re = data.image[IDre].array.F[ii];
            im = data.image[IDim].array.F[ii];
            amp = re*re+im*im;
            if(amp>0.2*peak)
            {
                data.image[IDre].array.F[ii] -= reave;
                data.image[IDim].array.F[ii] -= imave;
            }
        }
        mk_complex_from_reim("pupre1", "pupim1", "pupc", 0);
        permut("pupc");
        do2dfft("pupc","focc");
        permut("focc");
        delete_image_ID("pupc");
        mk_amph_from_complex("focc", "foca", "focp", 0);
        delete_image_ID("focc");
        delete_image_ID("focp");
        execute_arith("focic=foca*foca");
        delete_image_ID("foca");
        delete_image_ID("pupre1");
        delete_image_ID("pupim1");


        ID = image_ID("foci");
        ID1 = image_ID("focic");
        peak = 1.0;
        if(rtime>etimestart)
        {
            peak = 0.0;
            for(ii=0; ii<size*size; ii++)
            {
                data.image[IDpsfcumul].array.F[ii] += data.image[ID].array.F[ii];
                data.image[IDpsfcumul1].array.F[ii] += data.image[ID1].array.F[ii];
                if(data.image[IDpsfcumul].array.F[ii]>peak)
                    peak = data.image[IDpsfcumul].array.F[ii];
            }
        }


        // SAVING RESULT
        if(1)
        {
            sprintf(fname, "!psf/psf_%1.10lf.fits", rtime);
            save_fl_fits("foci",fname);
            printf("tot = %g\n", tot);


            if(rtime>etimestart)
            {
                tot = (rtime-etimestart)*FLUXSCI; // total flux in image [ph]
                arith_image_cstmult("PSFcumul",1.0/peak,"PSFcumuln");
                save_fl_fits("PSFcumuln", "!PSFcumul.fits"); // normalized in contrast
                tot0 = arith_image_total("PSFcumuln");
                arith_image_cstmult_inplace("PSFcumuln",tot/tot0); // normalized to tot
                put_poisson_noise("PSFcumuln","PSFcumulnn");
                arith_image_cstmult_inplace("PSFcumulnn",tot0/tot); // re-normalized in contrast
                delete_image_ID("PSFcumuln");
                save_fl_fits("PSFcumulnn", "!PSFcumul_n.fits");

                arith_image_cstmult("PSFcumul1",1.0/peak,"PSFcumul1n");
                save_fl_fits("PSFcumul1n", "!PSFcumul1.fits");
                arith_image_cstmult_inplace("PSFcumul1n",tot/tot0);
                put_poisson_noise("PSFcumul1n","PSFcumul1nn");
                arith_image_cstmult_inplace("PSFcumul1nn",tot0/tot);
                delete_image_ID("PSFcumul1n");
                save_fl_fits("PSFcumul1nn", "!PSFcumul1_n.fits");
            }
        }

        rtime += dtime;
        wfstime += dtime;
        wfstime1 += dtime;

    }


    delete_image_ID("TelPup");
    delete_image_ID("atmopd");
    delete_image_ID("atmamp");
    delete_image_ID("atmopdsci");
    delete_image_ID("atmampsci");
    delete_image_ID("wfsopd");
    delete_image_ID("wfsamp");
    delete_image_ID("dmopdtmp");
    delete_image_ID("dmopd");
    delete_image_ID("wfsmesopd");
    delete_image_ID("wfsmesopdprev");
    delete_image_ID("wfsmesopdder");
    delete_image_ID("wfsmesopdint");
    delete_image_ID("sciopd");
    delete_image_ID("sciamp");
    delete_image_ID("pupa");
    delete_image_ID("pupp");
    delete_image_ID("PSFcumul");
    delete_image_ID("PSFcumul1");
    delete_image_ID("PSFcumulnn");
    delete_image_ID("PSFcumul1nn");
    delete_image_ID("focic");
    delete_image_ID("foci");



    sprintf(imname, "wfa%08ld",cubeindex0);
    IDtmp = image_ID(imname);
    if(IDtmp!=-1)
        delete_image_ID(imname);

    sprintf(imname, "wfp%08ld",cubeindex0);
    IDtmp = image_ID(imname);
    if(IDtmp!=-1)
        delete_image_ID(imname);

    sprintf(imname, "swfa%08ld",cubeindex0);
    IDtmp = image_ID(imname);
    if(IDtmp!=-1)
        delete_image_ID(imname);

    sprintf(imname, "swfp%08ld",cubeindex0);
    IDtmp = image_ID(imname);
    if(IDtmp!=-1)
        delete_image_ID(imname);

    sprintf(imname, "wfa%08ld",cubeindex1);
    IDtmp = image_ID(imname);
    if(IDtmp!=-1)
        delete_image_ID(imname);

    sprintf(imname, "wfp%08ld",cubeindex1);
    IDtmp = image_ID(imname);
    if(IDtmp!=-1)
        delete_image_ID(imname);

    sprintf(imname, "swfa%08ld",cubeindex1);
    IDtmp = image_ID(imname);
    if(IDtmp!=-1)
        delete_image_ID(imname);

    sprintf(imname, "swfp%08ld",cubeindex1);
    IDtmp = image_ID(imname);
    if(IDtmp!=-1)
        delete_image_ID(imname);



    return(value/valuecnt);
}



// custom AO related processing of a series of WF

int AtmosphericTurbulence_WFprocess()
{
    FILE *fp;
    char wf_file_name[200];
    long ID;
    long k, kk;
    double time;
    long size, sizec;
    int OK;
    char command[500];
    double pupsize;


    double Kp, Ki, Kd, Kdgain;
    double val;
    long iter;

    int r;



    // Std loop (no PID)
    Kp = 0.5;
    Ki = 0.0;
    Kd = 0.0;
    val = AtmosphericTurbulence_makePSF(Kp, Ki, Kd, Kdgain);
    exit(0);



    printf("CUSTOM PROCESSING\n");
    for(iter=15; iter<1000; iter++)
    {
        Kp = ran1(); //0.35; // loop gain
        Ki = ran1(); //0.0;
        Kd = ran1(); //0.1;
        Kdgain = ran1(); //0.5;

        Kp = 0.3;
        Ki = 0.0;
        Kd = 0.0;// ran1();
        //Kdgain = 0.5;


        val = AtmosphericTurbulence_makePSF(Kp, Ki, Kd, Kdgain);

        fp = fopen("res.log.txt","a");
        fprintf(fp, "%ld %e %e %e %e %e\n", iter, Kp, Ki, Kd, Kdgain, val);
        fclose(fp);
        sprintf(command, "cp result.log result_%05ld.log", iter);
        r = system(command);

        list_image_ID();
        exit(0);
    }

    OK = 1;
    k = 0;
    while(OK==1)
    {
        sprintf(wf_file_name,"wf550_%08ld.pha",k);
        ID = load_fits(wf_file_name, "tpmwfc", 1);
        if(ID==-1)
        {
            OK = 0;
        }
        else
        {
            size = data.image[ID].md[0].size[0];
            sizec = data.image[ID].md[0].size[2];
            printf("%ld  %ld %ld\n",k,size,sizec);

            for(kk=0; kk<sizec; kk++)
            {

            }

            delete_image_ID("tmpwfc");
        }
        k++;
    }

    return(0);
}


int AtmosphericTurbulence_makeHV_CN2prof(double wspeed, double r0, double sitealt, long NBlayer, char *outfile)
{
    FILE *fp;
    double h;
    double CN2;
    double hstep, hmax;
    double A, A0;
    double CN2sum = 0.0;
    double r0val;
    double lambda = 0.55e-6;
    double Acoeff = 1.0;
    long iter;
    double Astep;
    long k;
    double l0, L0;
    double *layerarray_h;
    double *layerarray_CN2frac;
    double *layerarray_Wspeed;
    double *layerarray_sigmaWindSpeed;
    double *layerarray_Lwind;
    
    hmax = 30000.0;
    hstep = 1.0;
    A0 = 1.7e-14;

    Astep = 1.0;
    for(iter=0;iter<30;iter++)
    {
        A = A0*Acoeff;
        CN2sum = 0.0;
        for(h=sitealt; h<hmax; h+=hstep)
        {
            CN2 = 5.94e-53*pow(wspeed/27.0, 2.0)*pow(h,10.0)*exp(-h/1000.0) + 2.7e-16*exp(-h/1500.0) + A*exp(-(h-sitealt)/100.0);
            CN2sum += CN2/hstep;
        }
        r0val = 1.0/pow( 0.423*pow( 2.0*M_PI/lambda, 2.0)*CN2sum, 3.0/5.0);

      //  printf("Acoeff = %12f -> r0 = %10f m -> seeing = %10f arcsec\n", Acoeff, r0val, (lambda/r0val)/M_PI*180.0*3600.0);

        if(r0val>r0)
            Acoeff *= 1.0+Astep;
        else
            Acoeff /= 1.0+Astep;
        Astep *= 0.8;
    }


    fp = fopen("conf_turb.txt", "w");
    fprintf(fp, "%f\n", (lambda/r0val)/M_PI*180.0*3600.0);
    fclose(fp);

    layerarray_h = (double*) malloc(sizeof(double)*NBlayer);
    layerarray_CN2frac = (double*) malloc(sizeof(double)*NBlayer);
    layerarray_Wspeed = (double*) malloc(sizeof(double)*NBlayer);
    layerarray_sigmaWindSpeed = (double*) malloc(sizeof(double)*NBlayer);
    layerarray_Lwind = (double*) malloc(sizeof(double)*NBlayer);
    
    for(k=0;k<NBlayer;k++)
        {
            layerarray_h[k] = sitealt + pow(1.0*k/(NBlayer-1), 2.0)*(hmax-sitealt);
            layerarray_CN2frac[k] = 0.0;
        }

    hstep = 1.0;
    for(h=sitealt; h<hmax; h+=hstep)
        {
            CN2 = 5.94e-53*pow(wspeed/27.0, 2.0)*pow(h,10.0)*exp(-h/1000.0) + 2.7e-16*exp(-h/1500.0) + A*exp(-(h-sitealt)/100.0);
            k  = (long) (sqrt((h-sitealt)/(hmax-sitealt))*(1.0*NBlayer-1.0)+0.5);
            layerarray_CN2frac[k] += CN2;
        }
    
    fp = fopen(outfile, "w");
    fprintf(fp, "# altitude(m)   relativeCN2     speed(m/s)   direction(rad) outerscale[m] innerscale[m] sigmaWsp[m/s] Lwind[m]\n");
    fprintf(fp, "\n");
    for(k=0;k<NBlayer;k++)
        {
            layerarray_CN2frac[k] /= CN2sum;
            
            l0 = 0.008 + 0.072*pow(layerarray_h[k]/20000.0,1.6);
            
            if(layerarray_h[k]<14000.0)
                L0 = pow(10.0, 2.0 - 0.9*(layerarray_h[k]/14000.0) );
            else
                L0 = pow(10.0, 1.1 + 0.3*(layerarray_h[k]-14000.0)/6000.0);
            
            layerarray_Wspeed[k] = wspeed*(0.3+0.8*sqrt(1.0*k/(1.0+NBlayer)));
            layerarray_sigmaWindSpeed[k] = 0.1*layerarray_Wspeed[k];
            layerarray_Lwind[k] = 500.0;
            
            fprintf(fp, "%12f  %12f  %12f  %12f  %12f  %12f  %12f  %12f\n", layerarray_h[k], layerarray_CN2frac[k], layerarray_Wspeed[k], 2.0*M_PI*k/(1.0+NBlayer), L0, l0, layerarray_sigmaWindSpeed[k], layerarray_Lwind[k]);
        }
    fclose(fp);
    
    free(layerarray_h);
    free(layerarray_CN2frac);
    free(layerarray_Wspeed);
    free(layerarray_sigmaWindSpeed);
    free(layerarray_Lwind);
    
    return(0);
}


