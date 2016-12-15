#ifndef _LINARFILTERPRED_H
#define _LINARFILTERPRED_H


int init_linARfilterPred();


long LINARFILTERPRED_LoadASCIIfiles(double tstart, double dt, long NBpt, long NBfr, char *IDoutname);

long LINARFILTERPRED_SelectBlock(char *IDin_name, char *IDblknb_name, long blkNB, char *IDout_name);

long LINARFILTERPRED_Build_LinPredictor(char *IDin_name, long PForder, float PFlag, double SVDeps, double RegLambda, char *IDoutPF_name, int outMode, int LOOPmode, float LOOPgain);

long LINARFILTERPRED_Apply_LinPredictor_RT(char *IDfilt_name, char *IDin_name, char *IDout_name);

long LINARFILTERPRED_Apply_LinPredictor(char *IDfilt_name, char *IDin_name, float PFlag, char *IDout_name);

float LINARFILTERPRED_ScanGain(char* IDin_name, float multfact, float framelag);

long LINARFILTERPRED_PF_updatePFmatrix(char *IDPF_name, char *IDPFM_name, float alpha);

long LINARFILTERPRED_PF_RealTimeApply(char *IDmodevalOL_name, long IndexOffset, int semtrig, char *IDPFM_name, long NBPFstep, char *IDPFout_name, int nbGPU, long loop);

#endif
