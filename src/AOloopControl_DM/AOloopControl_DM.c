#include <fitsio.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <err.h>
#include <fcntl.h>
#include <sched.h>
#include <ncurses.h>
#include <semaphore.h>

#include "CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "AOloopControl_DM/AOloopControl_DM.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"




extern DATA data;

int wcol, wrow; // window size


struct timespec semwaitts;





#define DMSTROKE100 0.7 // um displacement for 100V

long NB_DMindex = 9;

AOLOOPCONTROL_DM_DISPCOMB_CONF *dmdispcombconf; // configuration
int dmdispcomb_loaded = 0;
int SMfd;


AOLOOPCONTROL_DMTURBCONF *dmturbconf; // DM turbulence configuration
int dmturb_loaded = 0;
int SMturbfd;





// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//


int AOloopControl_DM_CombineChannels_cli()
{
    // 1  long DMindex
    // 2  long xsize
    // 3  long ysize
    // 4  int NBchannel
    // 5  int AveMode
    // 6  int dm2dm_mode
    // 7  char *dm2dm_DMmodes
    // 8  char *dm2dm_outdisp
    // 9  int wfsrefmode
    // 10 char *wfsref_WFSRespMat
    // 11 char *wfsref_out
    // 12 int voltmode
    // 13 char *IDvolt_name
    // 14 float maxvolt
    if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,2)+CLI_checkarg(4,2)+CLI_checkarg(5,2)+CLI_checkarg(6,2)+CLI_checkarg(7,3)+CLI_checkarg(8,3)+CLI_checkarg(9,2)+CLI_checkarg(10,3)+CLI_checkarg(11,3)+CLI_checkarg(12,2)+(CLI_checkarg(13,3)*CLI_checkarg(13,4))+CLI_checkarg(14,1)==0)
        AOloopControl_DM_CombineChannels(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl, data.cmdargtoken[5].val.numl, data.cmdargtoken[6].val.numl, data.cmdargtoken[7].val.string, data.cmdargtoken[8].val.string, data.cmdargtoken[9].val.numl, data.cmdargtoken[10].val.string, data.cmdargtoken[11].val.string, data.cmdargtoken[12].val.numl, data.cmdargtoken[13].val.string, data.cmdargtoken[14].val.numf);
    else
        {// DEFAULT: no dm2dm, no wfsref, dmvolt output
            AOloopControl_DM_CombineChannels(00, 50, 50, 8, 1, 0, "dmmodes", "outdisp", 0, "wfsrm", "refout", 1, "dmvolt", 150.0);
        }
        
    return 1;
}



int AOloopControl_DM_chan_setgain_cli()
{
    if(CLI_checkarg(1,2)+CLI_checkarg(2,2)+CLI_checkarg(3,1)==0)
        AOloopControl_DM_chan_setgain(data.cmdargtoken[1].val.numl, data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf);
    else
        return 1;
}


int AOloopControl_DM_dmdispcomboff_cli()
{
        if(CLI_checkarg(1,2)==0)
        AOloopControl_DM_dmdispcomboff(data.cmdargtoken[1].val.numl);
    else
        return 1;
}


int AOloopControl_DM_dmdispcombstatus_cli()
{
    if(CLI_checkarg(1,2)==0)
        AOloopControl_DM_dmdispcombstatus(data.cmdargtoken[1].val.numl);
    else
        return 1;
}

int AOloopControl_DM_dmtrigoff_cli()
{
    if(CLI_checkarg(1,2)==0)
        AOloopControl_DM_dmtrigoff(data.cmdargtoken[1].val.numl);
    else
        return 1;
}


int AOloopControl_DM_dmturb_cli()
{
    if(CLI_checkarg(1,2)==0)
        AOloopControl_DM_dmturb(data.cmdargtoken[1].val.numl);
    else
        return 1;
}

int AOloopControl_DM_dmturboff_cli()
{
    if(CLI_checkarg(1,2)==0)
        AOloopControl_DM_dmturboff(data.cmdargtoken[1].val.numl);
    else
        return 1;
}


int AOloopControl_DM_dmturb_wspeed_cli()
{
    if(CLI_checkarg(1,2)+CLI_checkarg(2,1)==0)
        AOloopControl_DM_dmturb_wspeed(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf);
    else
        return 1;
}

int AOloopControl_DM_dmturb_ampl_cli()
{
    if(CLI_checkarg(1,2)+CLI_checkarg(2,1)==0)
        AOloopControl_DM_dmturb_ampl(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf);
    else
        return 1;
}

int AOloopControl_DM_dmturb_LOcoeff_cli()
{
    if(CLI_checkarg(1,2)+CLI_checkarg(2,1)==0)
        AOloopControl_DM_dmturb_LOcoeff(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf);
    else
        return 1;
}

int AOloopControl_DM_dmturb_tint_cli()
{
    if(CLI_checkarg(1,2)+CLI_checkarg(2,2)==0)
        AOloopControl_DM_dmturb_tint(data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numl);
    else
        return 1;
}





int init_AOloopControl_DM()
{
    strcpy(data.module[data.NBmodule].name, __FILE__);
    strcpy(data.module[data.NBmodule].info, "AO loop Control DM operation");
    data.NBmodule++;


    strcpy(data.cmd[data.NBcmd].key,"aolcontrolDMcomb");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOloopControl_DM_CombineChannels_cli;
    strcpy(data.cmd[data.NBcmd].info,"create and combine DM channels");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)> <xsize> <ysize> <NBchannel> <AveMode (1=if average level removed)> <dm2dm mode> <DMmodes> <outdm stream> <wfsref mode> <WFS resp mat> <wfsref stream> <voltmode (1=dmvolt computed)> <dmvoltname> <maxvolt [V]>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontrolDMcomb 0 50 50 8 0 1 dmmodes outdm 1 wfsrm wfsrefout 1 dmvolt 120.0");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_CombineChannels(long DMindex, long xsize, long ysize, int NBchannel, int AveMode, int dm2dm_mode, char *dm2dm_DMmodes, char *dm2dm_outdisp, int wfsrefmode, char *wfsref_WFSRespMat, char *wfsref_out, int voltmode, char *IDvolt_name, float maxvolt)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"aolcontroldmchgain");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOloopControl_DM_chan_setgain_cli;
    strcpy(data.cmd[data.NBcmd].info,"set gain for DM displacement channel");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)> <chan#> <gain>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmchgain 0 3 0.2");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_chan_setgain(long DMindex, int ch, float gain)");
    data.NBcmd++;



    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmcomboff");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp =  AOloopControl_DM_dmdispcomboff_cli;
    strcpy(data.cmd[data.NBcmd].info,"turn off DM combine");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmcomboff 0");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmdispcomboff(long DMindex)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmcombmon");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp =  AOloopControl_DM_dmdispcombstatus_cli;
    strcpy(data.cmd[data.NBcmd].info,"monitor DM comb program");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmcombmon 0");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmdispcombstatus(long DMindex)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmtrigoff");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp =  AOloopControl_DM_dmtrigoff_cli;
    strcpy(data.cmd[data.NBcmd].info,"turn off DM trigger");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmtrigoff 0");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmtrigoff(long DMindex)");
    data.NBcmd++;



    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmturb");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOloopControl_DM_dmturb_cli;
    strcpy(data.cmd[data.NBcmd].info,"DM turbulence");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontrolDMturb 0");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmturb(long DMindex)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmturboff");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp =  AOloopControl_DM_dmturboff_cli;
    strcpy(data.cmd[data.NBcmd].info,"turn off DM turbulence");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmturboff 0");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmturboff(long DMindex)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmturws");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOloopControl_DM_dmturb_wspeed_cli;
    strcpy(data.cmd[data.NBcmd].info,"set turbulence wind speed");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)> <wind speed [m/s]>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmturws 0 5.2");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmturb_wspeed(long DMindex, double wspeed);");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmturampl");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOloopControl_DM_dmturb_ampl_cli;
    strcpy(data.cmd[data.NBcmd].info,"set turbulence amplitude");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)> <amplitude [um]>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmturampl 0 0.1");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmturb_ampl(long DMindex, double ampl);");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmturlo");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOloopControl_DM_dmturb_LOcoeff_cli;
    strcpy(data.cmd[data.NBcmd].info,"set turbulence low order coefficient");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)> <coeff>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmturlo 0 0.2");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmturb_LOcoeff(long DMindex, double LOcoeff);");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"aoloopcontroldmturtint");
    strcpy(data.cmd[data.NBcmd].module,__FILE__);
    data.cmd[data.NBcmd].fp = AOloopControl_DM_dmturb_tint_cli;
    strcpy(data.cmd[data.NBcmd].info,"set turbulence interval time");
    strcpy(data.cmd[data.NBcmd].syntax,"<DMindex (0-9)> <interval time [us] long>");
    strcpy(data.cmd[data.NBcmd].example,"aoloopcontroldmturtint 0 200");
    strcpy(data.cmd[data.NBcmd].Ccall,"int AOloopControl_DM_dmturb_tint(long DMindex, long tint);");
    data.NBcmd++;





    // add atexit functions here
    atexit((void*) AOloopControl_DM_unloadconf);

    return 0;
}








struct timespec time_diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}



int AOloopControl_DM_disp2V(long DMindex)
{
    long ii;
    float volt;

    data.image[dmdispcombconf[DMindex].IDvolt].md[0].write = 1;
    for(ii=0; ii<dmdispcombconf[DMindex].xysize; ii++)
    {
        volt = 100.0*sqrt(data.image[dmdispcombconf[DMindex].IDdisp].array.F[ii]/DMSTROKE100);
        if(volt>dmdispcombconf[DMindex].MAXVOLT)
            volt = dmdispcombconf[DMindex].MAXVOLT;
        data.image[dmdispcombconf[DMindex].IDvolt].array.U[ii] = (unsigned short int) (volt/300.0*16384.0); //65536.0);
    }
    data.image[dmdispcombconf[DMindex].IDvolt].md[0].write = 0;
    data.image[dmdispcombconf[DMindex].IDvolt].md[0].cnt0++;
    COREMOD_MEMORY_image_set_sempost(data.image[dmdispcombconf[DMindex].IDdisp].name, -1);


    return 0;
}



int AOloopControl_printDMconf()
{
    long DMindex;
    
    printf("DM on   x   y Nbch busy maxvolt  monint stat IDdisp IDvolt voltname\n");
    for(DMindex=0; DMindex<NB_DMindex; DMindex++)
        {
            printf("%ld  %d  %3ld %3ld  %02ld %d     %6.2f  %8ld  %02d  %3ld %3ld  %s\n", DMindex, dmdispcombconf[DMindex].ON, dmdispcombconf[DMindex].xsize, dmdispcombconf[DMindex].ysize, dmdispcombconf[DMindex].NBchannel, dmdispcombconf[DMindex].busy, dmdispcombconf[DMindex].MAXVOLT, dmdispcombconf[DMindex].moninterval, dmdispcombconf[DMindex].status, dmdispcombconf[DMindex].IDdisp, dmdispcombconf[DMindex].IDvolt, dmdispcombconf[DMindex].voltname);
        }
    
    return(0);
}


int AOloopControl_DM_createconf()
{
    int result;
    int ch;
    char fname[200];
    long DMindex;
    char errstr[200];

    sprintf(fname, "/tmp/dmdispcombconf.conf.shm");

    if( dmdispcomb_loaded == 0 )
    {
        printf("Create/read DM configuration, %ld entries\n", NB_DMindex);

        SMfd = open(fname, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
        if (SMfd == -1) {
            sprintf(errstr, "Error opening (O_RDWR | O_CREAT | O_TRUNC) file \"%s\", function AOloopControl_DM_createconf", fname);
            perror(errstr);
            exit(EXIT_FAILURE);
        }

        result = lseek(SMfd, sizeof(AOLOOPCONTROL_DM_DISPCOMB_CONF)*NB_DMindex-1, SEEK_SET);
        if (result == -1) {
            close(SMfd);
            perror("Error calling lseek() to 'stretch' the file");
            exit(EXIT_FAILURE);
        }

        result = write(SMfd, "", 1);
        if (result != 1) {
            close(SMfd);
            perror("Error writing last byte of the file");
            exit(EXIT_FAILURE);
        }

        dmdispcombconf = (AOLOOPCONTROL_DM_DISPCOMB_CONF*)mmap(0, sizeof(AOLOOPCONTROL_DM_DISPCOMB_CONF)*NB_DMindex, PROT_READ | PROT_WRITE, MAP_SHARED, SMfd, 0);
        if (dmdispcombconf == MAP_FAILED) {
            close(SMfd);
            perror("Error mmapping the file");
            exit(EXIT_FAILURE);
        }

        for(DMindex=0; DMindex<NB_DMindex; DMindex++)
        {
            dmdispcombconf[DMindex].ON = 0;
            dmdispcombconf[DMindex].xsize = 0;
            dmdispcombconf[DMindex].ysize = 0;
            dmdispcombconf[DMindex].xysize = 0;
            dmdispcombconf[DMindex].NBchannel = 0;
            dmdispcombconf[DMindex].busy = 0;
            dmdispcombconf[DMindex].MAXVOLT = 150.0;
            dmdispcombconf[DMindex].moninterval = 30000; // 33Hz
            dmdispcombconf[DMindex].status = 0;

            dmdispcombconf[DMindex].IDdisp = -1;
            dmdispcombconf[DMindex].IDvolt = -1;
            sprintf(dmdispcombconf[DMindex].voltname, " ");

            for(ch=0; ch<DM_NUMBER_CHANMAX; ch++)
                {
                    dmdispcombconf[DMindex].dmdispID[ch] = -1;
                    dmdispcombconf[DMindex].dmdispgain[ch] = 1.0;
                }
        }
        dmdispcomb_loaded = 1;

    }
    AOloopControl_printDMconf();
    
    return 0;
}




int AOloopControl_DM_loadconf()
{
    int result;
    char fname[200];
    char errstr[200];

    sprintf(fname, "/tmp/dmdispcombconf.conf.shm");



    if( dmdispcomb_loaded == 0 )
    {
        printf("Create/read DM configuration\n");

        SMfd = open(fname, O_RDWR, (mode_t)0600);
        if (SMfd == -1) {
            AOloopControl_DM_createconf();
        }
        else
        {
            dmdispcombconf = (AOLOOPCONTROL_DM_DISPCOMB_CONF*)mmap(0, sizeof(AOLOOPCONTROL_DM_DISPCOMB_CONF)*NB_DMindex, PROT_READ | PROT_WRITE, MAP_SHARED, SMfd, 0);
            if (dmdispcombconf == MAP_FAILED) {
                close(SMfd);
                printf("Error mmapping the file -> creating it\n");
                AOloopControl_DM_createconf();
                //            exit(EXIT_FAILURE);
            }
        }

        dmdispcomb_loaded = 1;
    }
    AOloopControl_printDMconf();

    return 0;
}





int AOloopControl_DM_unloadconf()
{
    if( dmdispcomb_loaded == 1 )
    {
        if (munmap(dmdispcombconf, sizeof(AOLOOPCONTROL_DM_DISPCOMB_CONF)*NB_DMindex) == -1)
            perror("Error un-mmapping the file");
        close(SMfd);
        dmdispcomb_loaded = 0;
    }
    return 0;
}





//
// DMindex is a unique DM identifier (0-9), so multiple instances can coexist
//
// voltmode = 1 if DM volt computed
//
// AveMode: averaging mode 
//      0: do not appy DC offset command to average, but offset combined average to mid-range, and clip displacement at >0.0
//      1: apply DC offset to remove average
//      2: do not apply DC offset, do not offset sum, do not clip
//
// NOTE: DM displacement is biased to mid displacement
// NOTE: responds immediately to sem[1] in dmdisp
// dmdisp files have 5 semaphores
//
int AOloopControl_DM_CombineChannels(long DMindex, long xsize, long ysize, int NBchannel, int AveMode, int dm2dm_mode, char *dm2dm_DMmodes, char *dm2dm_outdisp, int wfsrefmode, char *wfsref_WFSRespMat, char *wfsref_out, int voltmode, char *IDvolt_name, float maxvolt)
{
    long naxis = 2;
    long *size;
    long ch;
    char name[200];
    long cnt = 0;
    long long cntsumold;
    long long cntsum;
    long ii;
    long IDdisp;
    long IDvolt;
    double ave;
    long ID1;
    int RT_priority = 95; //any number from 0-99
    struct sched_param schedpar;
    int r;
    long sizexy;
    float *dmdispptr;
    float *dmdispptr_array[20];
    long IDdispt;
    char sname[200];
    long nsecwait = 10000; // 10 us
    int vOK;
    float maxmaxvolt = 150.0;
    char errstr[200];
    
    long sizexyDMout;
    long IDtmpoutdm;
    long kk;
    long sizexywfsref;
    long IDtmpoutref;
    
    if(DMindex>NB_DMindex-1)
    {
        printf("ERROR: requested DMindex (%02ld) exceeds maximum number of DMs (%02ld)\n", DMindex, NB_DMindex);
        exit(0);
    }
    
    
    schedpar.sched_priority = RT_priority;
    r = seteuid(euid_called); //This goes up to maximum privileges
    sched_setscheduler(0, SCHED_FIFO, &schedpar); //other option is SCHED_RR, might be faster
    r = seteuid(euid_real);//Go back to normal privileges

    AOloopControl_DM_createconf();
    AOloopControl_DM_loadconf();
    
    dmdispcombconf[DMindex].ON = 1;
    dmdispcombconf[DMindex].xsize = xsize;
    dmdispcombconf[DMindex].ysize = ysize;
    dmdispcombconf[DMindex].xysize = xsize*ysize;
    dmdispcombconf[DMindex].NBchannel = NBchannel;
    dmdispcombconf[DMindex].MAXVOLT = maxvolt;
    sprintf(dmdispcombconf[DMindex].voltname, "%s", IDvolt_name);
    dmdispcombconf[DMindex].status = 0;
    
    dmdispcombconf[DMindex].DClevel = 0.5*(DMSTROKE100*dmdispcombconf[DMindex].MAXVOLT/100.0*dmdispcombconf[DMindex].MAXVOLT/100.0);

    printf("maxvolt = %f\n", maxvolt);


    size = (long*) malloc(sizeof(long)*naxis);
    size[0] = xsize;
    size[1] = ysize;
    sizexy = xsize*ysize;


    dmdispcombconf[DMindex].xsizeout = 0;
    dmdispcombconf[DMindex].ysizeout = 0;

    dmdispcombconf[DMindex].dm2dm_mode = dm2dm_mode;
    

   if(dm2dm_mode == 1) 
   {
        printf("INITIALIZATION AND VERIFICATION FOR dm2dm MODE\n");
        fflush(stdout);

        dmdispcombconf[DMindex].ID_dm2dm_DMmodes = image_ID(dm2dm_DMmodes);
        if(data.image[dmdispcombconf[DMindex].ID_dm2dm_DMmodes].md[0].naxis != 3)
            {
                sprintf(errstr, "image \"%s\" should have naxis = 3", dm2dm_DMmodes);
                printERROR(__FILE__,__func__,__LINE__, errstr);
                exit(0);
            }
        dmdispcombconf[DMindex].xsizeout = data.image[dmdispcombconf[DMindex].ID_dm2dm_DMmodes].md[0].size[0];
        dmdispcombconf[DMindex].ysizeout = data.image[dmdispcombconf[DMindex].ID_dm2dm_DMmodes].md[0].size[1];        
   
        dmdispcombconf[DMindex].ID_dm2dm_outdisp = image_ID(dm2dm_outdisp);
        if(data.image[dmdispcombconf[DMindex].ID_dm2dm_outdisp].md[0].size[0] != dmdispcombconf[DMindex].xsizeout)
            {
                sprintf(errstr, "image \"%s\" should have x axis = %ld", dm2dm_outdisp, dmdispcombconf[DMindex].xsizeout);
                printERROR(__FILE__,__func__,__LINE__, errstr);
                exit(0);
            }
         if(data.image[dmdispcombconf[DMindex].ID_dm2dm_outdisp].md[0].size[1] != dmdispcombconf[DMindex].ysizeout)
            {
                sprintf(errstr, "image \"%s\" should have y axis = %ld", dm2dm_outdisp, dmdispcombconf[DMindex].ysizeout);
                printERROR(__FILE__,__func__,__LINE__, errstr);
                exit(0);
            }
            
        IDtmpoutdm = create_2Dimage_ID("_tmpoutdm", dmdispcombconf[DMindex].xsizeout, dmdispcombconf[DMindex].ysizeout); 
        sizexyDMout = dmdispcombconf[DMindex].xsizeout*dmdispcombconf[DMindex].ysizeout;
        printf("done\n\n");
        fflush(stdout);
   }
   
   
   
    dmdispcombconf[DMindex].wfsrefmode = wfsrefmode;
    if(wfsrefmode == 1) 
    {
        printf("INITIALIZATION AND VERIFICATION FOR wfsref MODE\n");
        fflush(stdout);
    
        dmdispcombconf[DMindex].ID_wfsref_RespMat = image_ID(wfsref_WFSRespMat);
        if(data.image[dmdispcombconf[DMindex].ID_wfsref_RespMat].md[0].naxis != 3)
            {
                sprintf(errstr, "image \"%s\" should have naxis = 3", wfsref_WFSRespMat);
                printERROR(__FILE__,__func__,__LINE__, errstr);
                exit(0);
            }
        dmdispcombconf[DMindex].xsizewfsref = data.image[dmdispcombconf[DMindex].ID_wfsref_RespMat].md[0].size[0];
        dmdispcombconf[DMindex].ysizewfsref = data.image[dmdispcombconf[DMindex].ID_wfsref_RespMat].md[0].size[1];

        dmdispcombconf[DMindex].ID_wfsref_out = image_ID(wfsref_out);
        if(data.image[dmdispcombconf[DMindex].ID_wfsref_out].md[0].size[0] != dmdispcombconf[DMindex].xsizewfsref)
            {
                sprintf(errstr, "image \"%s\" should have x axis = %ld", wfsref_out, dmdispcombconf[DMindex].xsizewfsref);
                printERROR(__FILE__,__func__,__LINE__, errstr);
                exit(0);
            }
        if(data.image[dmdispcombconf[DMindex].ID_wfsref_out].md[0].size[1] != dmdispcombconf[DMindex].ysizewfsref)
            {
                sprintf(errstr, "image \"%s\" should have y axis = %ld", wfsref_out, dmdispcombconf[DMindex].ysizewfsref);
                printERROR(__FILE__,__func__,__LINE__, errstr);
                exit(0);
            }
        IDtmpoutref = create_2Dimage_ID("_tmpoutref", dmdispcombconf[DMindex].xsizewfsref, dmdispcombconf[DMindex].ysizewfsref);
        sizexywfsref = dmdispcombconf[DMindex].xsizewfsref*dmdispcombconf[DMindex].ysizewfsref;

       COREMOD_MEMORY_image_set_createsem(wfsref_out, 5);

        printf("done\n\n");
        fflush(stdout);
    }

    printf("Initialize channels\n");
    printf("Max DM stroke = %f um\n", DMSTROKE100*dmdispcombconf[0].MAXVOLT/100.0*dmdispcombconf[0].MAXVOLT/100.0);
    fflush(stdout);

    for(ch=0; ch<dmdispcombconf[DMindex].NBchannel; ch++)
    {
        sprintf(name, "dm%lddisp%ld", DMindex, ch);
        printf("Channel %ld \n", ch);
        dmdispcombconf[DMindex].dmdispID[ch] = create_image_ID(name, naxis, size, FLOAT, 1, 10);
        COREMOD_MEMORY_image_set_createsem(name, 5);
        dmdispptr_array[ch] = data.image[dmdispcombconf[DMindex].dmdispID[ch]].array.F;
    }


    sprintf(name, "dm%lddisp", DMindex);
    dmdispcombconf[DMindex].IDdisp = create_image_ID(name, naxis, size, FLOAT, 1, 10);
    COREMOD_MEMORY_image_set_createsem(name, 5);
    
    sprintf(name, "dm%lddispt", DMindex);
    IDdispt = create_image_ID(name, naxis, size, FLOAT, 0, 0);
    dmdispptr = data.image[IDdispt].array.F;

    if(voltmode==1)
    {
        IDvolt = image_ID(dmdispcombconf[DMindex].voltname);
        
        vOK = 0;
        if(IDvolt!=-1)
            {
                if((data.image[IDvolt].md[0].atype==USHORT)&&(data.image[IDvolt].md[0].naxis==2)&&(data.image[IDvolt].md[0].size[0]==xsize)&&(data.image[IDvolt].md[0].size[1]==ysize))
                        vOK = 1;
                else
                    delete_image_ID(dmdispcombconf[DMindex].voltname);
            }
        
        printf("vOK = %d\n", vOK);
        if(vOK==0)
        {
            dmdispcombconf[DMindex].IDvolt = create_image_ID(dmdispcombconf[DMindex].voltname, naxis, size, USHORT, 1, 10);
            COREMOD_MEMORY_image_set_createsem(dmdispcombconf[DMindex].voltname, 5);
         }
         else
            dmdispcombconf[DMindex].IDvolt = image_ID(dmdispcombconf[DMindex].voltname);
    }

    cntsumold = 0;

    
    dmdispcombconf[0].status = 1;

    sprintf(name, "dm%lddisp", DMindex);
    COREMOD_MEMORY_image_set_createsem(name, 5);

    if(data.image[dmdispcombconf[DMindex].IDdisp].sem<2)
    {
        printf("ERROR: image %s semaphore %d missing\n", data.image[dmdispcombconf[DMindex].IDdisp].name, 1);
        exit(0);
    }

    dmdispcombconf[DMindex].MAXVOLT = maxvolt;
    if(dmdispcombconf[DMindex].MAXVOLT>maxmaxvolt)
        dmdispcombconf[DMindex].MAXVOLT = maxvolt;
    
    
     AOloopControl_printDMconf();

    
    
    while(dmdispcombconf[DMindex].ON == 1)
    {
        dmdispcombconf[DMindex].status = 2;

        if (clock_gettime(CLOCK_REALTIME, &semwaitts) == -1) {
            perror("clock_gettime");
            exit(EXIT_FAILURE);
        }
        semwaitts.tv_nsec += nsecwait;
        if(semwaitts.tv_nsec >= 1000000000)
            semwaitts.tv_sec = semwaitts.tv_sec + 1;

        sem_timedwait(data.image[dmdispcombconf[DMindex].IDdisp].semptr[1], &semwaitts);

        cntsum = 0;



        for(ch=0; ch<dmdispcombconf[DMindex].NBchannel; ch++)
            cntsum += data.image[dmdispcombconf[DMindex].dmdispID[ch]].md[0].cnt0;


        
            
        if(cntsum != cntsumold)
        {
            dmdispcombconf[0].status = 3;
            cnt++;

            memcpy (data.image[IDdispt].array.F, dmdispptr_array[0], sizeof(float)*sizexy);
            for(ch=1; ch<dmdispcombconf[DMindex].NBchannel; ch++)
            {
                for(ii=0; ii<sizexy; ii++)
                    dmdispptr[ii] += dmdispcombconf[DMindex].dmdispgain[ch]*dmdispptr_array[ch][ii];
            }

            dmdispcombconf[DMindex].status = 4;

            ave = 0.0;
            if(AveMode == 1) // REMOVE DC LEVEL AND MOVE TO MEAN MOTION RANGE
                {
                    for(ii=0; ii<dmdispcombconf[DMindex].xysize; ii++)
                        ave += data.image[IDdispt].array.F[ii];
                    ave /= dmdispcombconf[DMindex].xysize;
                }
            dmdispcombconf[DMindex].status = 5;

            if(AveMode < 2)
            {
                    for(ii=0; ii<dmdispcombconf[DMindex].xysize; ii++)
                {
                    data.image[IDdispt].array.F[ii] += dmdispcombconf[DMindex].DClevel-ave;
                    if(data.image[IDdispt].array.F[ii]<0.0)
                        data.image[IDdispt].array.F[ii] = 0.0;
                }
            }
            dmdispcombconf[DMindex].status = 6;

            data.image[dmdispcombconf[DMindex].IDdisp].md[0].write = 1;
            memcpy (data.image[dmdispcombconf[DMindex].IDdisp].array.F,data.image[IDdispt].array.F, sizeof(float)*data.image[dmdispcombconf[DMindex].IDdisp].md[0].nelement);
            data.image[dmdispcombconf[DMindex].IDdisp].md[0].cnt0++;
            data.image[dmdispcombconf[DMindex].IDdisp].md[0].write = 0;            
            sem_post(data.image[dmdispcombconf[DMindex].IDdisp].semptr[0]);
 
            if(dm2dm_mode==1)
            {
                memset(data.image[IDtmpoutdm].array.F, '\0', sizeof(float)*sizexyDMout);
                for(kk=0;kk<data.image[dmdispcombconf[DMindex].IDdisp].md[0].nelement;kk++)
                    {
                        for(ii=0;ii<sizexyDMout;ii++)
                            data.image[IDtmpoutdm].array.F[ii] += data.image[dmdispcombconf[DMindex].IDdisp].array.F[kk] * data.image[dmdispcombconf[DMindex].ID_dm2dm_DMmodes].array.F[kk*sizexyDMout+ii];
                    }
                
                data.image[dmdispcombconf[DMindex].ID_dm2dm_outdisp].md[0].write = 1;
                memcpy (data.image[dmdispcombconf[DMindex].ID_dm2dm_outdisp].array.F,data.image[IDtmpoutdm].array.F, sizeof(float)*sizexyDMout);
                data.image[dmdispcombconf[DMindex].ID_dm2dm_outdisp].md[0].cnt0++;
                data.image[dmdispcombconf[DMindex].ID_dm2dm_outdisp].md[0].write = 0;            
                sem_post(data.image[dmdispcombconf[DMindex].ID_dm2dm_outdisp].semptr[0]);                
            }
            
            if(wfsrefmode==1)
            {
                memset(data.image[IDtmpoutref].array.F, '\0', sizeof(float)*sizexywfsref);
                list_image_ID();
                printf("kkmax = %ld\n", data.image[dmdispcombconf[DMindex].IDdisp].md[0].nelement);
                printf("iimax = %ld\n", sizexywfsref);
                printf("ID RespMat = %ld  (%ld)\n", dmdispcombconf[DMindex].ID_wfsref_RespMat, (data.image[dmdispcombconf[DMindex].IDdisp].md[0].nelement-1)*sizexywfsref + sizexywfsref-1);
                fflush(stdout);
                save_fits(wfsref_WFSRespMat, "!_test_wfsref_WFSRespMat.fits");
                for(kk=0;kk<data.image[dmdispcombconf[DMindex].IDdisp].md[0].nelement;kk++)
                    {
                        printf("(%ld %g) ", kk, data.image[dmdispcombconf[DMindex].IDdisp].array.F[kk]);
                        for(ii=0;ii<sizexywfsref;ii++)
                            data.image[IDtmpoutref].array.F[ii] += data.image[dmdispcombconf[DMindex].IDdisp].array.F[kk] * data.image[dmdispcombconf[DMindex].ID_wfsref_RespMat].array.F[kk*sizexywfsref+ii];
                    }
                printf("\n");
                printf("Updating Zero Point  %ld <- %ld\n", dmdispcombconf[DMindex].ID_wfsref_out, IDtmpoutref);
                fflush(stdout);   
                data.image[dmdispcombconf[DMindex].ID_wfsref_out].md[0].write = 1;
                memcpy (data.image[dmdispcombconf[DMindex].ID_wfsref_out].array.F,data.image[IDtmpoutref].array.F, sizeof(float)*sizexywfsref);
                data.image[dmdispcombconf[DMindex].ID_wfsref_out].md[0].cnt0++;
                data.image[dmdispcombconf[DMindex].ID_wfsref_out].md[0].write = 0;            
                sem_post(data.image[dmdispcombconf[DMindex].ID_wfsref_out].semptr[0]);   
                printf("Done\n");
                fflush(stdout);
            }
            
            
            dmdispcombconf[DMindex].status = 7;

            if(voltmode==1)
                AOloopControl_DM_disp2V(DMindex);

            dmdispcombconf[DMindex].status = 8;

            cntsumold = cntsum;
        }
    }

    if(voltmode==1)
        arith_image_zero(dmdispcombconf[DMindex].voltname);



    printf("LOOP STOPPED\n");
    fflush(stdout);

    free(size);
 
    delete_image_ID("_tmpoutdm");
    delete_image_ID("_tmpoutref");

    return 0;
}






int AOloopControl_DM_chan_setgain(long DMindex, int ch, float gain)
{
    
    AOloopControl_DM_loadconf();
 
    if(ch<dmdispcombconf[DMindex].NBchannel)
        dmdispcombconf[DMindex].dmdispgain[ch] = gain;

    return 0;
}





int AOloopControl_DM_dmdispcombstatus(long DMindex)
{
    long long mcnt = 0;
    int ch;

    AOloopControl_DM_loadconf();

    initscr();
    getmaxyx(stdscr, wrow, wcol);

    start_color();
    init_pair(1, COLOR_BLACK, COLOR_WHITE);
    init_pair(2, COLOR_BLACK, COLOR_RED);
    init_pair(3, COLOR_GREEN, COLOR_BLACK);
    init_pair(4, COLOR_RED, COLOR_BLACK);

    while( !kbdhit() )
    {
        usleep(dmdispcombconf[DMindex].moninterval);
        clear();
        attron(A_BOLD);
        print_header(" PRESS ANY KEY TO STOP MONITOR ", '-');
        attroff(A_BOLD);
        printw("    %ld\n", mcnt);
        printw("ON         %d\n", dmdispcombconf[DMindex].ON);
        printw("cnt       %ld\n", dmdispcombconf[DMindex].loopcnt);
        printw("updatecnt %ld\n", dmdispcombconf[DMindex].updatecnt);
        printw("busy      %d\n", dmdispcombconf[DMindex].busy);
        printw("MAXVOLT   %f\n", dmdispcombconf[DMindex].MAXVOLT);
        printw("status    %d\n",  dmdispcombconf[DMindex].status);
        printw("moninterval %d\n", dmdispcombconf[DMindex].moninterval);
        printw("\n");
        for(ch=0; ch<dmdispcombconf[DMindex].NBchannel; ch++)
        {
            printw("  %2d   %5.3f\n", ch, dmdispcombconf[DMindex].dmdispgain[ch]);
        }

        mcnt++;
        refresh();
    }
    endwin();

    return 0;
}






int AOloopControl_DM_dmdispcomboff(long DMindex)
{
    AOloopControl_DM_loadconf();
    dmdispcombconf[DMindex].ON = 0;

    return 0;
}

int AOloopControl_DM_dmtrigoff(long DMindex)
{
   AOloopControl_DM_loadconf();
    data.image[dmdispcombconf[DMindex].IDvolt].md[0].status = 101;

    return 0;
}











int AOloopControl_printDMturbconf()
{
    long DMindex;
    
    printf("ind on  ampl [um]  tint [us]  simtime [s]  wspeed [m/s]  LOcoeff\n");
    for(DMindex=0; DMindex<NB_DMindex; DMindex++)
        {
            printf("%ld  %d  %10f  %10ld  %10f %5f  %5f\n", DMindex, dmturbconf[DMindex].on, dmturbconf[DMindex].ampl, dmturbconf[DMindex].tint, dmturbconf[DMindex].simtime, dmturbconf[DMindex].wspeed, dmturbconf[DMindex].LOcoeff);
        }
    
    return 0;
}




//
// create configuration shared memory structure for moving turbulence screen for DM
// one configuration per DM
//
int AOloopControl_DMturb_createconf()
{
    int result;
    long IDc1;
    char name[200];
    long DMindex;
    char errstr[200];
    
    AOloopControl_DM_loadconf();


    if( dmturb_loaded == 0 )
    {
        printf("Create/read DMturb configuration\n");
        fflush(stdout);

        SMturbfd = open(DMTURBCONF_FILENAME, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
        if (SMturbfd == -1) {
            sprintf(errstr, "Error opening (O_RDWR | O_CREAT | O_TRUNC) file \"%s\"", name);
            perror(errstr);
            exit(EXIT_FAILURE);
        }

        result = lseek(SMturbfd, sizeof(AOLOOPCONTROL_DMTURBCONF)*NB_DMindex-1, SEEK_SET);
        if (result == -1) {
            close(SMturbfd);
            perror("Error calling lseek() to 'stretch' the file");
            exit(EXIT_FAILURE);
        }

        result = write(SMturbfd, "", 1);
        if (result != 1) {
            close(SMturbfd);
            perror("Error writing last byte of the file");
            exit(EXIT_FAILURE);
        }

        dmturbconf = (AOLOOPCONTROL_DMTURBCONF*)mmap(0, sizeof(AOLOOPCONTROL_DMTURBCONF)*NB_DMindex, PROT_READ | PROT_WRITE, MAP_SHARED, SMturbfd, 0);
        if (dmturbconf == MAP_FAILED) {
            close(SMturbfd);
            perror("Error mmapping the file");
            exit(EXIT_FAILURE);
        }


        for(DMindex=0; DMindex<NB_DMindex; DMindex++)
        {
            dmturbconf[DMindex].on = 0;
            dmturbconf[DMindex].ampl = 0.01; // [um]
            dmturbconf[DMindex].tint = 100; // [us]
            dmturbconf[DMindex].simtime = 0.0; // sec
            dmturbconf[DMindex].wspeed = 10.0; // [m/s]
            dmturbconf[DMindex].LOcoeff = 0.2;
        }
        dmturb_loaded = 1;

    }
    AOloopControl_printDMturbconf();
    
    return 0;
}








int AOloopControl_DMturb_loadconf(long DMindex)
{
    int result;
    char errstr[200];

    if( dmturb_loaded == 0 )
    {
        printf("Create/read configuration\n");

        SMturbfd = open(DMTURBCONF_FILENAME, O_RDWR, (mode_t)0600);
        if (SMturbfd == -1) {
            sprintf(errstr, "Error opening (O_RDWR) file \"%s\" in function AOloopControl_DMturb_loadconf", DMTURBCONF_FILENAME);
            perror(errstr);
            exit(EXIT_FAILURE);
        }

        dmturbconf = (AOLOOPCONTROL_DMTURBCONF*)mmap(0, sizeof(AOLOOPCONTROL_DMTURBCONF)*NB_DMindex, PROT_READ | PROT_WRITE, MAP_SHARED, SMturbfd, 0);
        if (dmturbconf == MAP_FAILED) {
            close(SMturbfd);
            printf("Error mmapping the file -> creating it\n");
            AOloopControl_DMturb_createconf();
        }
        dmturb_loaded = 1;
    }

    return 0;
}





int AOloopControl_DM_dmturboff(long DMindex)
{
    AOloopControl_DMturb_loadconf(DMindex);
    dmturbconf[DMindex].on = 0;
    AOloopControl_DM_dmturb_printstatus(DMindex);

    return 0;
}

int AOloopControl_DM_dmturb_wspeed(long DMindex, double wspeed)
{
    AOloopControl_DMturb_loadconf(DMindex);
    dmturbconf[DMindex].wspeed = wspeed;
    AOloopControl_DM_dmturb_printstatus(DMindex);

    return 0;
}

int AOloopControl_DM_dmturb_ampl(long DMindex, double ampl)
{
    AOloopControl_DMturb_loadconf(DMindex);
    dmturbconf[DMindex].ampl = ampl;
    AOloopControl_DM_dmturb_printstatus(DMindex);

    return 0;
}

int AOloopControl_DM_dmturb_LOcoeff(long DMindex, double LOcoeff)
{
    AOloopControl_DMturb_loadconf(DMindex);
    dmturbconf[DMindex].LOcoeff = LOcoeff;
    AOloopControl_DM_dmturb_printstatus(DMindex);

    return 0;
}

int AOloopControl_DM_dmturb_tint(long DMindex, long tint)
{
    AOloopControl_DMturb_loadconf(DMindex);
    dmturbconf[DMindex].tint = tint;
    AOloopControl_DM_dmturb_printstatus(DMindex);

    return 0;
}



int AOloopControl_DM_dmturb_printstatus(long DMindex)
{
    AOloopControl_DMturb_loadconf(DMindex);

    printf("Run time = %.3f sec\n", dmturbconf[DMindex].simtime);
    printf("\n");
    printf("cnt              : %ld   (ave frequ = %.2f kHz)\n", dmturbconf[DMindex].cnt, 0.001*dmturbconf[DMindex].cnt/dmturbconf[DMindex].simtime);
    printf("\n");

    if(dmturbconf[DMindex].on == 1)
        printf("LOOP IS ON\n");
    else
        printf("LOOP IS OFF\n");

    printf("ampl    =  %.2f um\n", dmturbconf[DMindex].ampl);
    printf("wspeed  =  %.2f m/s\n", dmturbconf[DMindex].wspeed);
    printf("tint    =  %ld us\n", dmturbconf[DMindex].tint);
    printf("LOcoeff =  %.2f\n", dmturbconf[DMindex].LOcoeff);
    printf("Requested uptdate frequ = %.2f kHz\n", 0.001/(1.0e-6*dmturbconf[DMindex].tint));
    printf("\n");
    printf("\n");

    return(0);
}


int AOloopControl_DM_dmturb(long DMindex)
{
    long size_sx; // screen size
    long size_sy;
    long IDs1, IDs2;
    char name[200];
    long imsize = 2048;
    struct timespec tlast;
    struct timespec tdiff;
    struct timespec tdiff1;
    double tdiff1v;

    float screen0_X;
    float screen0_Y;
    long ii, jj, ii1;
    float x, y;
    float xpix, ypix;
    float xpixf, ypixf;
    long xpix1, xpix2, ypix1, ypix2;
    float ave;

    double angle = 1.0;
    double coeff = 0.001;
    double x1, fx1;
    
    double RMSval;
    long RMSvalcnt;
    double r;
    
    float pixscale = 0.1; // [m/pix]
    // Subaru pupil ~ 80 pix diam
    // Single actuator ~7 pix
    
    long DM_Xsize, DM_Ysize;

    long IDturbs1;
    long IDturb;
    double totim;
    long IDk;


    AOloopControl_DMturb_createconf();

    IDs1 = load_fits("turbscreen1.fits", "screen1", 1);
    IDs2 = load_fits("turbscreen2.fits", "screen2", 1);
    list_image_ID();
    
    if(IDs1==-1)
    {
        make_master_turbulence_screen("screen1", "screen2", imsize, 200.0, 1.0);
        IDs1 = image_ID("screen1");
        IDk = make_gauss("kernim", imsize, imsize, 20.0, 1.0);
        totim = 0.0;
        for(ii=0;ii<imsize*imsize;ii++)
            totim += data.image[IDk].array.F[ii];
        for(ii=0;ii<imsize*imsize;ii++)
            data.image[IDk].array.F[ii] /= totim;
        IDs2 = fconvolve("screen1", "kernim", "screen2");
        delete_image_ID("kernim");
        save_fits("screen1", "!turbscreen1.fits");
        save_fits("screen2", "!turbscreen2.fits");
    }
    

    printf("ARRAY SIZE = %ld %ld\n", data.image[IDs1].md[0].size[0], data.image[IDs1].md[0].size[1]);
    size_sx = data.image[IDs1].md[0].size[0];
    size_sy = data.image[IDs1].md[0].size[1];

    clock_gettime(CLOCK_REALTIME, &dmturbconf[DMindex].tstart);
    dmturbconf[DMindex].tend = dmturbconf[DMindex].tstart;


    DM_Xsize = dmdispcombconf[DMindex].xsize;
    DM_Ysize = dmdispcombconf[DMindex].ysize;
    printf("DM %ld array size : %ld %ld\n", DMindex, DM_Xsize, DM_Ysize);
    list_image_ID();
    sprintf(name, "dm%lddisp1", DMindex);
    read_sharedmem_image(name);
    list_image_ID();
    dmturbconf[DMindex].on = 1;

    IDturbs1 = create_2Dimage_ID("turbs1", DM_Xsize, DM_Ysize);
    IDturb = create_2Dimage_ID("turbs", DM_Xsize, DM_Ysize);

    while(dmturbconf[DMindex].on == 1) // computation loop
    {
        usleep(dmturbconf[DMindex].tint);

        tlast = dmturbconf[DMindex].tend;
        clock_gettime(CLOCK_REALTIME, &dmturbconf[DMindex].tend);
        tdiff = time_diff(dmturbconf[DMindex].tstart, dmturbconf[DMindex].tend);
        tdiff1 =  time_diff(tlast, dmturbconf[DMindex].tend);
        tdiff1v = 1.0*tdiff1.tv_sec + 1.0e-9*tdiff1.tv_nsec;

        screen0_X += dmturbconf[DMindex].wspeed*tdiff1v*cos(angle); // [m]
        screen0_Y += dmturbconf[DMindex].wspeed*tdiff1v*sin(angle); // [m]


        dmturbconf[DMindex].simtime = 1.0*tdiff.tv_sec + 1.0e-9*tdiff.tv_nsec;


        for(ii=0; ii<DM_Xsize; ii++)
            for(jj=0; jj<DM_Ysize; jj++)
            {
                ii1 = jj*DM_Xsize+ii;

                x = 10.0*ii/DM_Xsize + screen0_X; // [m]
                y = 10.0*jj/DM_Ysize + screen0_Y; // [m]

                xpix = 0.5*size_sx + x/pixscale;
                ypix = 0.5*size_sy + y/pixscale;

                xpix1 = ((long) xpix)%size_sx;
                xpix2 = (xpix1+1)%size_sx;
                xpixf = xpix- (long) xpix;
                ypix1 = ((long) ypix)%size_sy;
                ypix2 = (ypix1+1)%size_sy;
                ypixf = ypix - (long) ypix;

                while(xpix1<0)
                    xpix1 = 0;
                while(xpix1>size_sx-1)
                    xpix1 = size_sx-1;

                if(ypix1<0)
                    ypix1 = 0;
                if(ypix1>size_sy-1)
                    ypix1 = size_sy-1;

                data.image[IDturbs1].array.F[ii1] = 1.0*xpix1;

                data.image[IDturb].array.F[ii1] = (1.0-xpixf)*(1.0-ypixf)*(data.image[IDs1].array.F[ypix1*size_sx+xpix1]-(1.0-dmturbconf[DMindex].LOcoeff)*data.image[IDs2].array.F[ypix1*size_sx+xpix1]);

                data.image[IDturb].array.F[ii1]  +=  (xpixf)*(1.0-ypixf)*(data.image[IDs1].array.F[ypix1*size_sx+xpix2]-(1.0-dmturbconf[DMindex].LOcoeff)*data.image[IDs2].array.F[ypix1*size_sx+xpix2]);

                data.image[IDturb].array.F[ii1]  += (1.0-xpixf)*(ypixf)*(data.image[IDs1].array.F[ypix2*size_sx+xpix1]-(1.0-dmturbconf[DMindex].LOcoeff)*data.image[IDs2].array.F[ypix2*size_sx+xpix1]);

                data.image[IDturb].array.F[ii1]  += xpixf*ypixf*(data.image[IDs1].array.F[ypix2*size_sx+xpix2]-(1.0-dmturbconf[DMindex].LOcoeff)*data.image[IDs2].array.F[ypix2*size_sx+xpix2]);
            }

        // proccess array
        ave = 0.0;
        for(ii1=0; ii1<DM_Xsize*DM_Ysize; ii1++)
            ave += data.image[IDturb].array.F[ii1];
        ave /= dmdispcombconf[DMindex].xysize;
        for(ii1=0; ii1<DM_Xsize*DM_Ysize; ii1++)
        {
            data.image[IDturb].array.F[ii1] -= ave;
            data.image[IDturb].array.F[ii1] *= coeff;
        }

        RMSval = 0.0;
        RMSvalcnt = 0;

        for(ii=0; ii<DM_Xsize; ii++)
            for(jj=0; jj<DM_Ysize; jj++)
            {
                ii1 = DM_Xsize*jj+ii;
                x = 0.5*DM_Xsize - 0.5 - ii;
                y = 0.5*DM_Ysize - 0.5 - jj;
                r = sqrt(x*x+y*y);
                if(r<DM_Xsize*0.5-1.0)
                {
                    RMSval += data.image[IDturb].array.F[ii1]*data.image[IDturb].array.F[ii1];
                    RMSvalcnt++;
                }
            }
        RMSval = sqrt(RMSval/RMSvalcnt);

        x1 = log10(RMSval/dmturbconf[DMindex].ampl);
        fx1 = 1.0 + 50.0*exp(-5.0*x1*x1);
        coeff /= pow(10.0,x1/fx1);

        
        printf("STEP 001  %f %f\n", screen0_X, screen0_Y);
        fflush(stdout);
        list_image_ID();
        
        sprintf(name, "dm%lddisp1", DMindex);
        copy_image_ID("turbs", name, 0);
        save_fits("turbs", "!turbs.fits");
        save_fits("turbs1", "!turbs1.fits");
    }


    return(0);
}



