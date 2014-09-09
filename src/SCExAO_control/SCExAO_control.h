#ifndef _SCEXAOCONTROL_H
#define _SCEXAOCONTROL_H

int init_SCExAO_control();


int SCExAOcontrol_mv_DMstage(long stepXpos, long stepYpos);
int SCExAOcontrol_PyramidWFS_AutoAlign_TT();


#endif
