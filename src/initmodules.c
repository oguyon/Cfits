#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "fft/fft.h"
#include "info/info.h"
#include "statistic/statistic.h"
#include "linopt_imtools/linopt_imtools.h"
#include "image_gen/image_gen.h"
#include "image_filter/image_filter.h"
#include "image_basic/image_basic.h"
#include "WFpropagate/WFpropagate.h"
#include "ZernikePolyn/ZernikePolyn.h"
#include "coronagraphs/coronagraphs.h"
#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"
#include "image_format/image_format.h"
#include "img_reduce/img_reduce.h"
#include "AOloopControl_DM/AOloopControl_DM.h"
#include "AOsystSim/AOsystSim.h"
#include "AOloopControl/AOloopControl.h"
#include "FPAOloopControl/FPAOloopControl.h"
#include "psf/psf.h"
#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "cudacomp/cudacomp.h"
#include "SCExAO_control/SCExAO_control.h"
#include "TransitLC/TransitLC.h"
#include "linARfilterPred/linARfilterPred.h"


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
  init_FPAOloopControl();
  init_psf();
  init_AtmosphereModel();
  init_AtmosphericTurbulence();
  init_cudacomp();
  init_SCExAO_control();
  init_TransitLC();
  init_linARfilterPred();
 
  return 0;
}
