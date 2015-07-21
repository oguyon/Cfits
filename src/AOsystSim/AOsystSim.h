#ifndef _AOSYSTSIM_H
#define _AOSYSTSIM_H



int AOsystSim_simpleAOfilter(char *IDin_name, char *IDout_name);

long AOsystSim_mkTelPupDM(char *ID_name, long msize, double xc, double xy, double rin, double rout, double pupPA, double spiderPA, double spideroffset, double spiderthick, double stretchx);

long AOsystSim_fitTelPup(char *ID_name, char *IDtelpup_name);

int init_AOsystSim();

int AOsystSim_run(int syncmode, long delayus);



typedef struct {
    double alpha; // angle
    double alpha_ld;
    double alpha_arcsec;
    double lambda0; // wavelength at seeing measurement [m]
    double lambdai; // imaging wavelength (science) [m]
    double lambdawfs; // WFS wavelength [m]
    
    double CN2layer_h[20]; // [m]
    double CN2layer_coeff[20]; // sum of coeffs = 1
    
    double D; // telescope aperture [m]
    double r0; // Fried parameter at lambda0 [m]
    double windspeed; // [m/s]
    double betapWFS; // WFS efficiency
    double betaaWFS; // WFS efficiency
    double betapWFSsci; // WFS efficiency Sci
    double betaaWFSsci; // WFS efficiency Sci
    double Fwfs; // source flux [ph.s-1.m-2]
    double Fsci;
    
    double framedelay; // delay in unit of frames
    
    double f_wfs; // spatial frequency
    double f_0; // spatial frequency
    double f; // spatial frequency
    double hf; // sine wave component amplitude
    double X;
    double Y;
    double dX;
    double dY;
    double C0;
    double C1;
    double twfs;
    double twfssci;
    double hfc; // corrected sine wave component amplitude
    double twfs_opt; // optimal WFS exposure time - phase
    double twfssci_opt; // optimal WFS exposure time - phase
    double twfs_opt_amp; // optimal WFS exposure time - amplitude
    double hfca; // time lag term
    double hfcb; // WFS photon noise
    double C2; // quadratic sum of C2a and C2b
    double C2_wfs; // quadratic sum of C2a and C2b
    double C3; // corrected amplitude
    double C4; // chromatic phase
    double C5; // chromatic amplitude
    double C6; // refractive index
    double Csum;
    double Csum_detection; // 5 sigma
    
    // near-IR loop:
    double RIC_hfca; // Refractive Index Chromaticity : time lag
    double RIC_hfcb; // Refractive Index Chromaticity : photon noise
    double RIC_hfc;
    double TL_hfca; // time lag correction: time lag
    double TL_hfcb; // time lag correction: photon noise
    double TL_hfc;
    double C7; // post-nearIR correction time lag
    double C8; // post-nearIR correction 
    double C9; // C4 -> 
    double C10; // C5 ->
    double C11; // post-nearIR correction refractive index
    double Csum2; 
    double Csum2ave; 
   
} EXAOSIMCONF;


int AOsystSim_extremeAO_contrast_sim();

#endif
