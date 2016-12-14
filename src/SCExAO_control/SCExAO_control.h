#ifndef _SCEXAOCONTROL_H
#define _SCEXAOCONTROL_H

int init_SCExAO_control();


int SCExAOcontrol_mv_DMstage(long stepXpos, long stepYpos);

long SCExAOcontrol_Average_image(char *imname, long NbAve, char *IDnameout, long semindex);
int SCExAOcontrol_PyramidWFS_AutoAlign_TT(char *WFScam_name);
int SCExAOcontrol_PyramidWFS_AutoAlign_cam(char *WFScam_name);

int SCExAOcontrol_PyramidWFS_Pcenter(char *IDwfsname, float prad, float poffset);
int SCExAOcontrol_Pyramid_flattenRefWF(char *WFScam_name, long zimaxmax, float ampl0);
int SCExAOcontrol_optPSF(char *WFScam_name, long zimaxmax, float alpha);
int SCExAOcontrol_SAPHIRA_cam_process(char *IDinname, char *IDoutname);

long SCExAOcontrol_vib_ComputeCentroid(char *IDin_name, char *IDdark_name, char *IDout_name);
long SCExAOcontrol_vib_mergeData(char *IDacc_name, char *IDttpos_name, char *IDout_name);

#endif
