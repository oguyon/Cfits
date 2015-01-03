#ifndef _SCEXAOCONTROL_H
#define _SCEXAOCONTROL_H

int init_SCExAO_control();


int SCExAOcontrol_mv_DMstage(long stepXpos, long stepYpos);

long SCExAOcontrol_Average_image(char *imname, long NbAve, char *IDnameout);
int SCExAOcontrol_PyramidWFS_AutoAlign_TT(char *WFScam_name);
int SCExAOcontrol_PyramidWFS_AutoAlign_cam(char *WFScam_name);

int SCExAOcontrol_PyramidWFS_Pcenter(char *IDwfsname, float prad, float poffset);
int SCExAOcontrol_Pyramid_flattenRefWF();

int SCExAOcontrol_SAPHIRA_cam_process(char *IDinname, char *IDoutname);

#endif
