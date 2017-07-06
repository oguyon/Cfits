/**
 * @file    linARfilterPred.h
 * @brief   Function prototypes for linear autoregressive prediction
 * 
 * Implements Empirical Orthogonal Functions
 * 
 * @author  O. Guyon
 * @date    5 Jul 2017
 *
 * @bug No known bugs. 
 * 
 */

#ifndef _LINARFILTERPRED_H
#define _LINARFILTERPRED_H


int_fast8_t init_linARfilterPred();



/* =============================================================================================== */
/* =============================================================================================== */
/** @name 1. INITIALIZATION, configurations
 *  
 */
///@{                                                                                         
/* =============================================================================================== */
/* =============================================================================================== */

///@}






/* =============================================================================================== */
/* =============================================================================================== */
/** @name 2. I/O TOOLS
 *  
 */
///@{                                                                                         
/* =============================================================================================== */
/* =============================================================================================== */

int NBwords(const char sentence[ ]);

long LINARFILTERPRED_LoadASCIIfiles(double tstart, double dt, long NBpt, long NBfr, const char *IDoutname);

long LINARFILTERPRED_SelectBlock(const char *IDin_name, const char *IDblknb_name, long blkNB, const char *IDout_name);

///@}





/* =============================================================================================== */
/* =============================================================================================== */
/** @name 3. BUILD PREDICTIVE FILTER
 *  
 */
///@{                                                                                         
/* =============================================================================================== */
/* =============================================================================================== */

/** @brief Build predictive filter
 * 
 * 
 * @param[in]  IDin_name     Input telemetry, can be a 2D or 3D image
 * @param[in]  PForder       Filter order: number of time steps in filter 
 * @param[in]  SVDeps        Cutoff limit on singular values
 * @param[in]  RegLambda     Regularization parameter
 * @param[out] IDoutPF_name  Output predictive filter name
 * @param[in]  outMode       Output mode. 0: do not write individual files, 1: write individual files (note: output filter cube is always written)
 * @param[in]  LOOPmode      1 if running in infinite loop waiting for input telemetry
 * @param[in]  LOOPgain      if running in loop, mixing coefficient between previous and current filter
 * 
 * 
 * Optional pixel masks select input and output variables: "inmask" and "outmask"
 * 
 * 
 * 
 * if LOOPmode = 1, operate in a loop, and re-run filter computation everytime IDin_name changes
 * 
 * @note if atmospheric wavefronts, data should be piston-free
 * 
 * @return output filter image index
 */
  
long LINARFILTERPRED_Build_LinPredictor(const char *IDin_name, long PForder, float PFlag, double SVDeps, double RegLambda, const char *IDoutPF_name, int outMode, int LOOPmode, float LOOPgain);

///@}




/* =============================================================================================== */
/* =============================================================================================== */
/** @name 4. APPLY PREDICTIVE FILTER
 *  
 */
///@{                                                                                         
/* =============================================================================================== */
/* =============================================================================================== */


long LINARFILTERPRED_Apply_LinPredictor_RT(const char *IDfilt_name, const char *IDin_name, const char *IDout_name);

long LINARFILTERPRED_Apply_LinPredictor(const char *IDfilt_name, const char *IDin_name, float PFlag, const char *IDout_name);

long LINARFILTERPRED_PF_updatePFmatrix(const char *IDPF_name, const char *IDPFM_name, float alpha);

long LINARFILTERPRED_PF_RealTimeApply(const char *IDmodevalOL_name, long IndexOffset, int semtrig, const char *IDPFM_name, long NBPFstep, const char *IDPFout_name, int nbGPU, long loop, long NBiter, int SAVEMODE, float tlag, long PFindex);

///@}





/* =============================================================================================== */
/* =============================================================================================== */
/** @name 5. MISC TOOLS, DIAGNOSTICS
 *  
 */
///@{                                                                                         
/* =============================================================================================== */
/* =============================================================================================== */

float LINARFILTERPRED_ScanGain(char* IDin_name, float multfact, float framelag);

#endif
