#ifndef _AOLOOPCONTROL_H
#define _AOLOOPCONTROL_H



#define AOLOOPCONTROL_FILENAME_CONF "/tmp/aoloopcontrol.conf.shm"

int init_AOloopControl();

int AOloopControl_run();



typedef struct
{
  int ON;

  char name[80];
  

} AOLOOPCONTROL_CONF;




#endif
