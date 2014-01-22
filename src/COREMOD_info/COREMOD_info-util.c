/*^****************************************************************************
* FILE: COREMOD_info-util.c : Module utility file
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
#include <COREMOD_info.h>     // Header for this module


/*
* Forward references for the module glue functions. 
* Naming convention: mod_ModuleName_CommandName 
*/
PF mod_COREMOD_info_imgnbpixflux           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_img_percentile           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_imgmkhistoc           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_img_mkhisto           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_img_ssquare           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_imstat           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_profile           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_profile2im           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_printpix           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_backnoise           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_tstrfunc           ( struct lxrcmd *c ); 
PF mod_COREMOD_info_strfunc           ( struct lxrcmd *c ); 


/*
* Command-control-blocks for this module
*/
struct lxrcmd mod_COREMOD_info_cmds[] = {

    {    "imgnbpixflux",
         mod_COREMOD_info_imgnbpixflux,
         "cumulative number of pixels/flux",
         "cumulative number of pixels/flux",
         "imgnbpixflux im",
         "%s %s",
         2,
         "img_nbpix_flux",
    },
    {    "img_percentile",
         mod_COREMOD_info_img_percentile,
         "image percentile",
         "image percentile",
         "img_percentile im1 0.2",
         "%s %s %f",
         3,
         "printf",
    },
    {    "imgmkhistoc",
         mod_COREMOD_info_imgmkhistoc,
         "make cumulative histogram of pixel values",
         "str1 is image name, str2 is output file",
         "imgmkhistoc im histo.dat",
         "%s %s %s",
         3,
         "img_histoc",
    },
    {    "img_mkhisto",
         mod_COREMOD_info_img_mkhisto,
         "creates an histogram",
         "min,max = f1,f2   l1=number of steps",
         "img_mkhisto im1 histo 0 10000.0 1000",
         "%s %s %s %f %f %ld",
         6,
         "make_histogram",
    },
    {    "img_ssquare",
         mod_COREMOD_info_img_ssquare,
         "sum of squared pixels",
         "sum of squared pixels",
         "img_ssquare im1",
         "%s %s",
         2,
         "printf",
    },
    {    "imstat",
         mod_COREMOD_info_imstat,
         "various info about the image",
         "various info about the image",
         "stats im1",
         "%s %s",
         2,
         "stats",
    },
    {    "profile",
         mod_COREMOD_info_profile,
         "radial profile",
         "str1 is the image file, str2 the output. center is (f1,f2). f3 is the step and n1 the number of steps.",
         "profile im1 im.prof 512 512 2.5 100",
         "%s %s %s %f %f %f %d",
         7,
         "profile",
    },
    {    "profile2im",
         mod_COREMOD_info_profile2im,
         "make 2D image from a radial profile",
         "str1 is the radial profile file name of l1 points. l2 is the size of the output image str2, of center f1,f2 and f3 is the radius of the disk.",
         "profile2im prof.dat 1000 512 200 out",
         "%s %s %ld %ld %f %f %f %s",
         8,
         "profile2im",
    },
    {    "printpix",
         mod_COREMOD_info_printpix,
         "print pixel values",
         "print pixel values of str1 in file str2",
         "printpix im im.pix",
         "%s %s %s",
         3,
         "printpix",
    },
    {    "backnoise",
         mod_COREMOD_info_backnoise,
         "prints the background noise",
         "prints the background noise",
         "backnoise im1",
         "%s %s",
         2,
         "printf",
    },
    {    "tstrfunc",
         mod_COREMOD_info_tstrfunc,
         "test the structure function",
         "str1 is the input file, l1 is the number of test points",
         "tstrfunc in 100000 out",
         "%s %s %ld %s",
         4,
         "test_structure_function",
    },
    {    "strfunc",
         mod_COREMOD_info_strfunc,
         "structure function",
         "str1 is the input file",
         "strfunc in out",
         "%s %s %s",
         3,
         "fft_structure_function",
    },
};




/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_init        :  The initialization function for this module.
|                       Naming convention: mod_ModuleName_init
|                       This function is declared in the include file:
|                       COREMOD_info.h  included via Cfits.h which was generated at build 
|                       time from the file: modules-included.list
|
|   struct module *m : This module's module-control-block
|
|
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_init ( struct module *m )					    
{										    
  int i;		    	    
  FILE *fp;  struct stat finfo;
  char str0[200];
  char str1[200];


   if( Debug > 1) fprintf(stdout, "[mod_COREMOD_info_init]\n");			    
										    
   // Set module name
   sprintf(m->name,"COREMOD_info");
   // Set module info name
   sprintf(m->info,"info about images");

   // Set number of commands for this module					    
   m->ncommands = sizeof(mod_COREMOD_info_cmds) / sizeof(struct lxrcmd) ;
   // Set command control block
   m->cmds = mod_COREMOD_info_cmds;
		 								    
   // set module-control-block index in every command-control-block 
    for( i=0; i<m->ncommands; ++i ) {
        mod_COREMOD_info_cmds[i].cdata_01 = m->module_number; 
    }

   // Set module compile time ascii entry
   sprintf(str0, "unknown");
   sprintf(str1, "%s/COREMOD_info.o", OBJDIR);
if (!stat(str1, &finfo)) {
   sprintf(str0, "%s", asctime(localtime(&finfo.st_mtime)));}
else { printf("ERROR: cannot find file %s\n",str1);}
   strncpy( m->buildtime, str0, strlen(str0)-1);
   m->buildtime[strlen(str0)] = '\0';

    return(PASS);
}






/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_imgnbpixflux : command function for user command "imgnbpixflux"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_imgnbpixflux ( struct lxrcmd *c ) 
{

  img_nbpix_flux ( c->args[1].v.s );   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_img_percentile : command function for user command "img_percentile"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_img_percentile ( struct lxrcmd *c ) 
{

    printf ( c->args[1].v.s, c->args[2].v.f);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_imgmkhistoc : command function for user command "imgmkhistoc"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_imgmkhistoc ( struct lxrcmd *c ) 
{

    img_histoc ( c->args[1].v.s, c->args[2].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_img_mkhisto : command function for user command "img_mkhisto"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_img_mkhisto ( struct lxrcmd *c ) 
{

    make_histogram ( c->args[1].v.s, c->args[2].v.s, c->args[3].v.f, c->args[4].v.f, c->args[5].v.ld);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_img_ssquare : command function for user command "img_ssquare"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_img_ssquare ( struct lxrcmd *c ) 
{

  printf ( "%f\n", ssquare(c->args[1].v.s) );   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_imstat : command function for user command "imstat"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_imstat ( struct lxrcmd *c ) 
{
  
  stats ( c->args[1].v.s, c->options);   
  
  c->results = NULL;   // NULL results block. Change as appropriate


  return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_profile : command function for user command "profile"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_profile ( struct lxrcmd *c ) 
{

    profile ( c->args[1].v.s, c->args[2].v.s, c->args[3].v.f, c->args[4].v.f, c->args[5].v.f, c->args[6].v.d);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_profile2im : command function for user command "profile2im"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_profile2im ( struct lxrcmd *c ) 
{

    profile2im ( c->args[1].v.s, c->args[2].v.ld, c->args[3].v.ld, c->args[4].v.f, c->args[5].v.f, c->args[6].v.f, c->args[7].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_printpix : command function for user command "printpix"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_printpix ( struct lxrcmd *c ) 
{

    printpix ( c->args[1].v.s, c->args[2].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_backnoise : command function for user command "backnoise"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_backnoise ( struct lxrcmd *c ) 
{

  printf ( "%f\n", background_photon_noise(c->args[1].v.s));   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_tstrfunc : command function for user command "tstrfunc"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_tstrfunc ( struct lxrcmd *c ) 
{

    test_structure_function ( c->args[1].v.s, c->args[2].v.ld, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_COREMOD_info_strfunc : command function for user command "strfunc"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_COREMOD_info_strfunc ( struct lxrcmd *c ) 
{

    fft_structure_function ( c->args[1].v.s, c->args[2].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}
