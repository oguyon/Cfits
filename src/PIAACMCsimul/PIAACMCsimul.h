#ifndef _PIAACMCSIMUL_H
#define _PIAACMCSIMUL_H



#define ApoFitCosFact 1.0


//
// *****************************************************************************************************
// -------------------------- structure defining a reflective PIAACMC system ---------------------------
// *****************************************************************************************************



//
// this structure holds parameters to be optimized in the PIAACMC diffractive design
//
typedef struct {

    // ======= SEED RADIAL PIAACMC PARAMETERS ======

    double centObs0;     /**< input central obstruction */
    double centObs1;     /**< output central obstruction */
    double r0lim;        /**< outer radius after extrapolation, piaa mirror 0 */
    double r1lim;        /**< outer radius after extrapolation, piaa mirror 1 */
    long NBradpts;       /**< number of points for common r0, r1, piaa sags 1D table */


    // Wavelength
    int nblambda;
    double lambda; // central wavelength [m]
    double lambdaB; // spectral bandwidth [%]
    double lambdaarray[2000]; // [m]  lambdaarray is also defined in OptSystProp structure


    // ====== Overall OPTICAL Geometry ===============

    float beamrad; // [m]
    long size;
    float pixscale; // [m/pix]

	int PIAAmode; // 0: no PIAA, 1: PIAA

    float PIAA0pos; // conjugation (z) of first PIAA surface [m]
    float PIAAsep;// separation between PIAA surfaces [m]
    int prePIAA0mask; // 1 if mask before PIAA surface 0
    float prePIAA0maskpos; // position of mask before PIAA surface 0 [m]
    int postPIAA0mask; // 1 if mask after PIAA surface 0
    float postPIAA0maskpos; // position of mask after PIAA surface 0 [m]
	float PIAAcoeff; // fraction of apodization done by PIAA
    int invPIAAmode; // 0: no inv PIAA, 1: inv PIAA after Lyot stops, 2: inv PIAA before Lyot stops
    float LyotZmin;
    float LyotZmax;
    float pupoutmaskrad; // output pupil mask radius (scaled to pupil radius)

    // ========== WAVEFRONT CONTROL ==================
    int nbDM; // number of deformable mirrors (10 max)
    double DMpos[10]; // DM conjugation in collimated space
    long ID_DM[10];  // DM image identifier


    // ========= LYOT STOPS ============
    long NBLyotStop;          /**< Number of Lyot stops */
    long IDLyotStop[10];
    double LyotStop_zpos[10];

    // ======= Optics shapes modes ============
    char PIAAmaterial_name[10]; 
    int PIAAmaterial_code; 
    long CmodesID; // Cosine radial mode
    long Cmsize; // cosine modes size
    long NBCmodes;
    long piaaNBCmodesmax; // maximum number of radial cosine modes for PIAA optics

    long FmodesID; // Fourier 2D modes
    long Fmsize;
    long NBFmodes;
    float piaaCPAmax; // maximum spatial frequency (CPA) for PIAA optics

    long piaa0CmodesID;
    long piaa0FmodesID;
    long piaa1CmodesID;
    long piaa1FmodesID;

    // PSF flux calib
    float peakPSF;


    // ========= Focal Plane Mask ============

    double fpmaskradld; // mask radius [l/d] for the idealized PIAACMC starting point
    long focmNBzone; // number of zones
    double Fratio; // beam Fratio at focal plane
    long zonezID;  // focm zone material thickness, double precision image, named fpmzt / fpm_zonez.fits
    double fpmaskamptransm; // mask amplitude transmission (normally 1.0)
    long zoneaID;  // focm zone amplitude transmission, double precision image, named fpmza / fpm_zonea.fits
    double fpzfactor; // focal plane mask DFT zoom factor

    double fpmRad; // outer radius of physical focal plane mask [m]

    long NBrings; // number of rings
    double fpmminsag; // [m]
    double fpmmaxsag; // [m]
    double fpmsagreg_coeff;
    double fpmsagreg_alpha;
    long NBringCentCone; // number of rings that the central cone occupies
    double fpmCentConeRad; // [m]
    double fpmCentConeZ; // peak sag of central cone [m]
    double fpmOuterConeZ; // outer sag offset [m]
    double fpmOuterConeRadld; // outer radius (end of outer cone) [lambda/D] 
    double fpmOuterConeRad; // [m]
    long fpmarraysize;
    char fpmmaterial_name[10]; 
    int fpmmaterial_code; 
 
    // Mask description
    //
    // CENTRAL CONE
    // inner zone (optional) is a cone, covering the central NBringCentCone rings (set to zero for no central cone)
    // The outer edge of the central cone is at sag = 0, and the central peak of the cone is at sag fpmCentConeZ
    //
    // SECTORS
    // The next zone outwards consists of sectors arranged in rings. Each sector has its own sag
    // The zones are between sag fpmminsag and fpmmaxsag
    // 
    // OUTER CONE
    // The outer cone starts at the outer edge of the sectors, where it has sag=0
    // Its outer is at sag fpmOuterZ
    // outside of the outer cone, the sag is constant at fpmOuterZ
    //
 
 
} OPTPIAACMCDESIGN;






/* =============================================================================================== */
/* =============================================================================================== */
/** @name 1. INITIALIZATION, configurations
 *  Allocate memory, import/export configurations
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */

/**
 * @brief Module initialization
 *
 * Registers command line interface (CLI) commands
 * 
 */
int_fast8_t init_PIAACMCsimul();


void PIAACMCsimul_logFunctionCall(char *LogFileName, const char *FunctionName, long line, char *comments);



/**
 * @brief Free PIAACMC memory
 * 
 */
void  PIAACMCsimul_free( void );


/**
 * @brief initializes the optsyst structure to simulate reflective PIAACMC system
 */
void PIAACMCsimul_init( OPTPIAACMCDESIGN *design, long index, double TTxld, double TTyld );


/**
 * @brief initializes configuration
 */
int PIAAsimul_initpiaacmcconf(long piaacmctype, double fpmradld, double centobs0, double centobs1, int WFCmode, int load);

/**
 * @brief Save configuration
 */
int PIAAsimul_savepiaacmcconf(const char *dname);

/**
 * @brief Load configuration
 */
int PIAAsimul_loadpiaacmcconf(const char *dname);


///@}




/* =============================================================================================== */
/* =============================================================================================== */
/** @name 2. Focal plane mask construction 
 *  Define focal plane mask geometry 
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */

long PIAACMCsimul_mkFPM_zonemap(const char *IDname);

long PIAACMCsimul_rings2sectors(const char *IDin_name, const char *sectfname, const char *IDout_name);

long PIAACMCsimul_mkFocalPlaneMask(const char *IDzonemap_name, const char *ID_name,  int mode, int saveMask);

///@}



/* =============================================================================================== */
/* =============================================================================================== */
/** @name  3. PIAA optics  (geometrical optics)
 *  Create PIAA opics according to geometrical optics
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */

uint_fast8_t PIAACMCsimul_load2DRadialApodization(const char *IDapo_name, float beamradpix, const char *IDapofit_name);

int PIAACMCsimul_init_geomPIAA_rad(const char *IDapofit_name);

int PIAACMCsimul_mkPIAAMshapes_from_RadSag(const char *fname, const char *ID_PIAAM0_name, const char *ID_PIAAM1_name);

int PIAACMCsimul_makePIAAshapes(OPTPIAACMCDESIGN *design, long index);

///@}



/* =============================================================================================== */
/* =============================================================================================== */
/** @name  4. Lyot stop(s)
 *  Create, optimize and manage Lyot stop(s)
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */

long PIAAsimul_mkSimpleLyotStop(const char *ID_name, float rin, float rout);

double PIAACMCsimul_optimizeLyotStop(const char *IDamp_name, const char *IDpha_name, const char *IDincoh_name, float zmin, float zmax, double throughput, long NBz, long NBmasks);

long PIAACMCsimul_mkLyotMask(const char *IDincoh_name, const char *IDmc_name, const char *IDzone_name, double throughput, const char *IDout_name);

///@}


/* =============================================================================================== */
/* =============================================================================================== */
/** @name  5. Focal plane mask optimization  
 *  Create, optimize and manage Focal plane solutions
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */

double PIAACMCsimul_achromFPMsol_eval(double *fpmresp_array, double *zonez_array, double *dphadz_array, double *outtmp_array, long vsize, long nbz, long nbl);

double PIAACMCsimul_achromFPMsol_eval_zonezderivative(long zone, double *fpmresp_array, double *zonez_array, double *dphadz_array, double *outtmp_array, long vsize, long nbz, long nbl);

long PIAACMC_FPMresp_rmzones(const char *FPMresp_in_name, const char *FPMresp_out_name, long NBzones);

long PIAACMC_FPMresp_resample(const char *FPMresp_in_name, const char *FPMresp_out_name, long NBlambda, long PTstep);

///@}


/* =============================================================================================== */
/* =============================================================================================== */
/** @name  6. Focal plane processing
 *  Process / resample focal plane solutions
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */


long PIAACMC_FPM_process(const char *FPMsag_name, const char *zonescoord_name, long NBexp, const char *outname);

long PIAACMC_FPMresp_resample(const char *FPMresp_in_name, const char *FPMresp_out_name, long NBlambda, long PTstep);

///@}



/* =============================================================================================== */
/* =============================================================================================== */
/** @name  7. High level routines 
 *  High level optimization and evaluation routines
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */

int PIAACMCsimul_exec(const char *confindex, long mode);

double PIAACMCsimul_computePSF(float xld, float yld, long startelem, long endelem, int savepsf, int sourcesize, int extmode, int outsave);

long PIAACMCsimul_CA2propCubeInt(const char *IDamp_name, const char *IDpha_name, float zmin, float zmax, long NBz, const char *IDout_name);

int PIAACMCsimul_run(const char *confindex, long mode);

///@}



#endif
