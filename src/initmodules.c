#include "CLIcore.h"

extern DATA data;

int init_modules()
{

  init_COREMOD_memory();
  init_COREMOD_arith();
  init_COREMOD_iofits();
  init_COREMOD_tools();
  init_fft();
  init_info();
  init_statistic();
  init_linopt_imtools();
  init_image_gen();
  init_image_filter();
  init_image_basic();
  init_WFpropagate();
  init_ZernikePolyn();
  init_coronagraphs();
  init_OptSystProp();
  init_PIAACMCsimul();
  init_image_format();
  init_img_reduce();
  init_AOloopControl_DM();
  init_AOsystSim();
  init_AOloopControl();
  init_psf();
  init_AtmosphereModel();
  init_AtmosphericTurbulence();
  init_cudacomp();
  init_SCExAO_control();
  init_TransitLC();
 
  return 0;
}
