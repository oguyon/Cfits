#ifndef _SCEXAOCONTROL_H
#define _SCEXAOCONTROL_H

int init_SCExAO_control();


int SCExAOcontrol_mv_DMstage(long stepXpos, long stepYpos);

long SCExAOcontrol_Average_image(char *imname, long NbAve, char *IDnameout);
int SCExAOcontrol_PyramidWFS_AutoAlign_TT();
int SCExAOcontrol_PyramidWFS_AutoAlign_cam();

#endif
