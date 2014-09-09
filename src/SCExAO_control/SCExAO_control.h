#ifndef _SCEXAOCONTROL_H
#define _SCEXAOCONTROL_H

int init_SCExAO_control();


int SCExAOcontrol_mv_DMstage(long stepXpos, long stepYpos);

long SCExAOcontrol_TakePyrWFS_image(char *IDname, long NbAve);
int SCExAOcontrol_PyramidWFS_AutoAlign_TT();
int SCExAOcontrol_PyramidWFS_AutoAlign_cam();

#endif
