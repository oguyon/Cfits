#ifndef _FPAOLOOPCONTROL_H
#define _FPAOLOOPCONTROL_H




typedef struct
{
    char name[80];
	
	// DM stream
    char dmCname[80];
    char dmRMname[80];
	long dmxsize;
	long dmysize;
	long dmsize;
	
	// Focal plane image stream
	char WFSname[80];
	long sizexWFS;
	long sizeyWFS;
	long sizeWFS;
	
	// timing info
	float loopfrequ; // Hz
	float hardwlatency; // hardware latency between DM command and WFS response [sec] 
	float hardwlatency_frame; // hardware latency between DM command and WFS response 
	
	
	// ============= RESPONSE CALIBRATION ===================
	float fpim_normFlux; // total focal plane flux in the absence of a coronagraph
	float fpim_Xcent;
	float fpim_Ycent;
	

	// ======= LEVEL 1 CALIBRATION ==============
	// Each actuator influence function has the same amplitude, phase is ramp set accordingly to actuator position
	// to be acquired without coronagraph

	
	
	
} FPAOLOOPCONTROL_CONF;



int init_FPAOloopControl();




long FPAOloopControl_InitializeMemory(int mode);
int FPAOloopControl_loadconfigure(long loop, int mode, int level);


// RM Calibration

// level 1
// Each actuator influence function has the same amplitude, phase is ramp set accordingly to actuator position
// to be acquired without coronagraph
long FPAOloopControl_acquireRM_level1(float ampl);



#endif
