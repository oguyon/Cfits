/*^****************************************************************************
* FILE: COREMOD_fft-util.c : Module utility file
*
*
*     File naming convention: Modulename-util.c
*
*
* modules-config.h was generated at build time by module names 
* in modules-included.list
*
*******************************************************************************/
#include <sys/stat.h>
#include <time.h>
#include <Cfits.h>  // Generic  header
#include <COREMOD_fft.h>     // Header for this module


/*
* Forward references for the module glue functions. 
* Naming convention: mod_ModuleName_CommandName 
*/
PF mod_COREMOD_fft_initfft           ( struct lxrcmd *c ); 
PF mod_COREMOD_fft_permut           ( struct lxrcmd *c ); 
PF mod_COREMOD_fft_fft2d           ( struct lxrcmd *c ); 
PF mod_COREMOD_fft_fft2di           ( struct lxrcmd *c ); 
PF mod_COREMOD_fft_fft2dr           ( struct lxrcmd *c ); 
PF mod_COREMOD_fft_fft2dri           ( struct lxrcmd *c ); 
PF mod_COREMOD_fft_correl          ( struct lxrcmd *c );
PF mod_COREMOD_fft_fzoom           ( struct lxrcmd *c ); 
PF mod_COREMOD_fft_testfft           ( struct lxrcmd *c ); 

PF mod_COREMOD_fft_DFT              ( struct lxrcmd *c );



/*
* Command-control-blocks for this module
*/
struct lxrcmd mod_COREMOD_fft_cmds[] = {

    {    "initfft",
         mod_COREMOD_fft_initfft,
         "init fftw",
         "init fftw",
         "initfft",
         "%s", 
         1,
         "init_fftw_plans0",
    },
    {    "permut",
         mod_COREMOD_fft_permut,
         "quadrants permutation",
         "quadrants permutation",
         "permut im",
         "%s %s",
         2,
         "permut",
    },
    {    "fft2d",
         mod_COREMOD_fft_fft2d,
         "2d fft",
         "str1 is input, str2 is output",
         "fft2d in out",
         "%s %s %s",
         3,
         "do2dfft",
    },
    {    "fft2di",
         mod_COREMOD_fft_fft2di,
         "2d inverse fft",
         "str1 is input, str2 is output",
         "fft2di in out",
         "%s %s %s",
         3,
         "do2dffti",
    },
    {    "fft2dr",
         mod_COREMOD_fft_fft2dr,
         "2d real fft",
         "str1 is input (real), str2 is output (complex)",
         "fft2dr in out",
         "%s %s %s",
         3,
         "do2drfft",
    },
    {    "fft2dri",
         mod_COREMOD_fft_fft2dri,
         "2d inverse real fft",
         "str1 is input (complex), str2 is output (real)",
         "fft2dri in out",
         "%s %s %s",
         3,
         "do2drffti",
    },
    {    "fcorrel",
         mod_COREMOD_fft_correl,
         "fourier correlation map",
         "str1,str2 are input, str3 is output.",
         "fcorrel im1 im2 out",
         "%s %s %s %s",
         4,
         "ftt_correlation",
    },
    {    "fzoom",
         mod_COREMOD_fft_fzoom,
         "fourier zoom",
         "str1 is input, str2 is output, l1 is the zoom factor.",
         "fzoom in out 2",
         "%s %s %s %ld",
         4,
         "fftzoom",
    },
    {    "testfft",
         mod_COREMOD_fft_testfft,
         "test FFT speed",
         "test FFT speed up to size 2^n",
         "testfft",
         "%s %d",
         2,
         "test_fftspeed",
    },
    {    "custdft",
         mod_COREMOD_fft_DFT,
         "custom DFT",
         "discrete Fourier transform. str1 is input complex image, str2 is mask specifying which pixels to use in input, str3 is output result, str4 is mask specifying pixels to compute in output, f1 is zoom factor, d1 is direction",
         "custdft in inmask out outmask 1.0 1",
         "%s %s %s %s %s %f %d",
         7,
         "COREMOD_fft_DFT",
    },
};




/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_init        :  The initialization function for this module.
|                       Naming convention: mod_ModuleName_init
|                       This function is declared in the include file:
|                       COREMOD_fft.h  included via Cfits.h which was generated at build 
|                       time from the file: modules-included.list
|
|   struct module *m : This module's module-control-block
|
|
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_init ( struct module *m )
{
  int i;		    	    
  FILE *fp;
  struct stat finfo;
  char str0[200];
  char str1[200];
  

  // Set module name
  sprintf( m->name, "COREMOD_fft");
  if( Debug > 1) fprintf(stdout, "[mod_%s_init]\n",m->name);
  
  // Set module info line
  sprintf(m->info,"Fast Fourier Transform (FFT) routines with fftw library");
   
  
  // Set number of commands for this module
  m->ncommands = sizeof(mod_COREMOD_fft_cmds) / sizeof(struct lxrcmd) ;
  // Set command control block
  m->cmds = mod_COREMOD_fft_cmds;   
   
  // set module-control-block index in every command-control-block 
  for( i=0; i<m->ncommands; ++i ) {
    mod_COREMOD_fft_cmds[i].cdata_01 = m->module_number; 
  }
   

  // Set module compile time ascii entry
  sprintf(str0, "unknown");
  sprintf(str1, "%s/%s.o", OBJDIR, m->name);
  if (!stat(str1, &finfo)) {
    sprintf(str0, "%s", asctime(localtime(&finfo.st_mtime)));}
  else { fprintf(stderr,"%c[%d;%dm WARNING: cannot find file %s [ %s  %s  %d ] %c[%d;m\n",(char) 27, 1, 31, str1, __FILE__, __func__, __LINE__, (char) 27, 0);}
  strncpy( m->buildtime, str0, strlen(str0)-1);
  m->buildtime[strlen(str0)] = '\0';    
  
  return(PASS);
}






/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_initfft : command function for user command "initfft"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_initfft ( struct lxrcmd *c ) 
{

    init_fftw_plans0 ( );   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_permut : command function for user command "permut"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_permut ( struct lxrcmd *c ) 
{

    permut ( c->args[1].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_fft2d : command function for user command "fft2d"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_fft2d ( struct lxrcmd *c ) 
{

    do2dfft ( c->args[1].v.s, c->args[2].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_fft2di : command function for user command "fft2di"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_fft2di ( struct lxrcmd *c ) 
{

    do2dffti ( c->args[1].v.s, c->args[2].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_fft2dr : command function for user command "fft2dr"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_fft2dr ( struct lxrcmd *c ) 
{

    do2drfft ( c->args[1].v.s, c->args[2].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_fft2dri : command function for user command "fft2dri"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_fft2dri ( struct lxrcmd *c ) 
{

    do2drffti ( c->args[1].v.s, c->args[2].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}



PF mod_COREMOD_fft_correl ( struct lxrcmd *c ) 
{
  fft_correlation ( c->args[1].v.s, c->args[2].v.s, c->args[3].v.s );   
  c->results = NULL;   // NULL results block. Change as appropriate
    
  return(PASS);
}





/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_fzoom : command function for user command "fzoom"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_fzoom ( struct lxrcmd *c ) 
{

    fftzoom ( c->args[1].v.s, c->args[2].v.s, c->args[3].v.ld);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_fft_testfft : command function for user command "testfft"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_fft_testfft ( struct lxrcmd *c ) 
{

    test_fftspeed ( c->args[1].v.d);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}



PF mod_COREMOD_fft_DFT              ( struct lxrcmd *c )
{
  COREMOD_fft_DFT( c->args[1].v.s, c->args[2].v.s, c->args[3].v.s, c->args[4].v.s, c->args[5].v.f, c->args[5].v.d); 
  c->results = NULL;   // NULL results block. Change as appropriate	
  return(PASS);
}

