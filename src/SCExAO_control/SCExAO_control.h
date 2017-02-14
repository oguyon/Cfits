#ifndef _SCEXAOCONTROL_H
#define _SCEXAOCONTROL_H

int_fast8_t init_SCExAO_control();

long SCExAOcontrol_mkSegmentModes(const char *IDdmmap_name, const char *IDout_name);
int SCExAOcontrol_mv_DMstage(long stepXpos, long stepYpos);

long SCExAOcontrol_Average_image(const char *imname, long NbAve, const char *IDnameout, long semindex);
int SCExAOcontrol_PyramidWFS_AutoAlign_TT(const char *WFScam_name, float XposStart, float YposStart);
int SCExAOcontrol_PyramidWFS_AutoAlign_cam(const char *WFScam_name);

int SCExAOcontrol_PyramidWFS_Pcenter(const char *IDwfsname, float prad, float poffset);
int SCExAOcontrol_Pyramid_flattenRefWF(const char *WFScam_name, long zimaxmax, float ampl0);
int SCExAOcontrol_optPSF(const char *WFScam_name, long zimaxmax, float alpha);
int SCExAOcontrol_SAPHIRA_cam_process(const char *IDinname, const char *IDoutname);

long SCExAOcontrol_vib_ComputeCentroid(const char *IDin_name, const char *IDdark_name, const char *IDout_name);
long SCExAOcontrol_vib_mergeData(const char *IDacc_name, const char *IDttpos_name, const char *IDout_name, int mode);

#endif
