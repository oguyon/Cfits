#ifndef _FPAOLOOPCONTROL_H
#define _FPAOLOOPCONTROL_H




typedef struct
{
	// DM stream
	long dmRM_ID;
	long dmC_ID;
	long dmxsize;
	long dmysize;
	
	// Focal plane image stream
	long FPim_ID;
	long FPim_dark_ID;
	long fpimxsize;
	long fpimysize;
	
} FPAOLOOPCONTROL_CONF;



int init_FPAOloopControl();





long FPAOloopControl_initMem(long loop, char *IDdmRM_name, char *IDdmC_name, char *IDfpim_name, char *IDfpim_dark_name, int mode);

#endif
