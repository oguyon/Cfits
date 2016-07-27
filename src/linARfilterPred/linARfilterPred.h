#ifndef _LINARFILTERPRED_H
#define _LINARFILTERPRED_H

long LINARFILTERPRED_SelectBlock(char *IDin_name, char *IDblknb_name, long blkNB, char *IDout_name);

long LINARFILTERPRED_Build_LinPredictor(char *IDin_name, long PForder, float PFlag, double SVDeps, double RegLambda, char *IDoutPF_name, int outMode);

long LINARFILTERPRED_Apply_LinPredictor(char *IDfilt_name, char *IDin_name, float PFlag, char *IDout_name);

float LINARFILTERPRED_ScanGain(char* IDin_name, float multfact, float framelag);

#endif
