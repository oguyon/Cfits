/*^****************************************************************************
* FILE: coronagraphs-util.c : Module utility file
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
#include <coronagraphs.h>     // Header for this module


/*
* Forward references for the module glue functions. 
* Naming convention: mod_ModuleName_CommandName 
*/
PF mod_coronagraphs_cormk2Dprolate           ( struct lxrcmd *c ); 
PF mod_coronagraphs_cormk2Dprolateld          ( struct lxrcmd *c ); 
PF mod_coronagraphs_corup2Dprolate           ( struct lxrcmd *c ); 
PF mod_coronagraphs_cormk2Dprolate_CS           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimaic           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsim4qpm           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimbl8           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimrrpm           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimcpa           ( struct lxrcmd *c ); 
PF mod_coronagraphs_mksubpup           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimpiaa           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimpiaac           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimpiaacaic           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimaplc           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimpsf           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corsimtransm           ( struct lxrcmd *c ); 
PF mod_coronagraphs_corouserf           ( struct lxrcmd *c ); 
PF mod_coronagraphs_scanPIAACMC_centObs_perf  ( struct lxrcmd *c ); 

/*
* Command-control-blocks for this module
*/
struct lxrcmd mod_coronagraphs_cmds[] = {

    {    "cormk2Dprolate",
         mod_coronagraphs_cormk2Dprolate,
         "make 2D prolate",
         "f1 is mask diameter in pixel, f2 is central obs, str1 is the output",
         "cormk2Dprolate 2.0 0.0 prol2",
         "%s %f %f %s",
         4,
         "coronagraph_make_2Dprolate",
    },
    {    "cormk2Dprolateld",
         mod_coronagraphs_cormk2Dprolateld,
         "make 2D prolate, mask size in l/d",
         "f1 is mask diameter in l/d, f2 is central obs, str1 is the output",
         "cormk2Dprolateld 1.0 0.0 prol2",
         "%s %f %f %s",
         4,
         "coronagraph_make_2Dprolateld",
    },
    {    "corup2Dprolate",
         mod_coronagraphs_corup2Dprolate,
         "update 2D prolate",
         "f1 is mask diameter in l/d, f2 is central obs, f3 is DFT zoom factor",
         "corup2Dprolate 3.8 0.0 2.0",
         "%s %f %f %f",
         4,
         "coronagraph_update_2Dprolate",
    },
    {    "cormk2DprolateCS",
         mod_coronagraphs_cormk2Dprolate_CS,
         "make 2D prolate for complex shaped apertures",
         "f1 is mask diameter in pixel, str1 is the output",
         "cormk2DprolateCS 2.0 prol2",
         "%s %f %s",
         3,
         "coronagraph_make_2DprolateCS",
    },
    {    "corsimaic",
         mod_coronagraphs_corsimaic,
         "AIC coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsimiac 1.0 1.0 psf",
         "%s %f %f %s",
         4,
         "coronagraph_simul_AIC",
    },
    {    "corsim4qpm",
         mod_coronagraphs_corsim4qpm,
         "4QPM coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsim4qpm 1.0 1.0 psf",
         "%s %f %f %s",
         4,
         "coronagraph_simul_4QPM",
    },
    {    "corsimbl8",
         mod_coronagraphs_corsimbl8,
         "BL8 coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsimbl8 1.0 1.0 psf",
         "%s %f %f %s",
         4,
         "coronagraph_simul_BL8",
    },
    {    "corsimrrpm",
         mod_coronagraphs_corsimrrpm,
         "Roddier phase mask coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsimrrpm 1.0 1.0 psf",
         "%s %f %f %s",
         4,
         "coronagraph_simul_RRPM",
    },
    {    "corsimcpa",
         mod_coronagraphs_corsimcpa,
         "Classical Pupil Apodization coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsimcpa 1.0 1.0 psf",
         "%s %f %f %s",
         4,
         "coronagraph_simul_CPA",
    },
    {    "mksubpup",
         mod_coronagraphs_mksubpup,
         "make Subaru pupil",
         "make Subaru pupil",
         "mksubpup",
         "%s", 
         1,
         "coronagraphs_make_SUBARU_pupil",
    },
    {    "corsimpiaa",
         mod_coronagraphs_corsimpiaa,
         "PIAA coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsimpiaa 1.0 1.0 psf",
         "%s %f %f %s",
         4,
         "coronagraph_simul_PIAA",
    },
    {    "corsimpiaac",
         mod_coronagraphs_corsimpiaac,
         "PIAAC coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsimpiaac 1.0 1.0 psf",
         "%s %f %f %s",
         4,
         "coronagraph_simul_PIAAC",
    },
    {    "corsimpiaacaic",
         mod_coronagraphs_corsimpiaacaic,
         "PIAAC/AIC coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsimpiaacaic 1.0 1.0 psf1",
         "%s %f %f %s",
         4,
         "coronagraph_simul_AIC_PIAAC",
    },
    {    "corsimaplc",
         mod_coronagraphs_corsimaplc,
         "APLC simulations (multistep)",
         "f1,f2 is the point source offset in l/d, str1 is the output psf",
         "corsimaplc 1.0 1.0 psf1",
         "%s %f %f %s",
         4,
         "coronagraph_simul_MULTISTEP_APLC",
    },
    {    "corsimpsf",
         mod_coronagraphs_corsimpsf,
         "coronagraph simulation",
         "f1,f2 is the point source offset in l/d, str1 is the output psf, l1 is the coronagraph type (1: AIC) (2: phase mask) (3: 4 quadrants) (4: PIAA/CPA hybrid) (5: PIAAC/CPA hybrid) (6: 8 order BL, circular) (7: 8 order BL, linear) (8: AIC/PIAAC) (9: CPA).  Options include -FPMSERR xxx for focal plane mask size error.",
         "corsimpsf 1.0 1.0 psf1 1",
         "%s %f %f %s %ld",
         5,
         "coronagraph_simulPSF",
    },
    {    "corsimtransm",
         mod_coronagraphs_corsimtransm,
         "transmission",
         "str1 is output file name, l1 is the coronagraph type (1: AIC) (2: phase mask) (3: 4 quadrants) (4: PIAA/CPA hybrid) (5: PIAAC/CPA hybrid) (6: 8 order BL, circular) (7: 8 order BL, linear) (8: AIC/PIAAC) (9: CPA), f1 is logcontrast. Options include -TT xxx, where xxx is the tip-tilt error in l/d.",
         "corsimtransm transm.dat 1 10.0",
         "%s %s %ld %f",
         4,
         "coronagraph_transm",
    },
    {    "corouserf",
         mod_coronagraphs_corouserf,
         "user function",
         "user function",
         "*",
         "%s", 
         1,
         "coronagraph_userfunc",
    },
    {    "coroscanpiaacmcperf",
	 mod_coronagraphs_scanPIAACMC_centObs_perf,
	 "scan PIAACMC perf",
	 "f1 is input pupil central obstruction",
	 "coroscanpiaacmcperf 0.3",
	 "%s %f",
	 2,
	 "coronagraphs_scanPIAACMC_centObs_perf",
    },
};




/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_init        :  The initialization function for this module.
|                       Naming convention: mod_ModuleName_init
|                       This function is declared in the include file:
|                       coronagraphs.h  included via Cfits.h which was generated at build 
|                       time from the file: modules-included.list
|
|   struct module *m : This module's module-control-block
|
|
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_init ( struct module *m )					    
{										    
  int i;		    	    
  FILE *fp;  struct stat finfo;
  char str0[200];
  char str1[200];


   if( Debug > 1) fprintf(stdout, "[mod_coronagraphs_init]\n");			    
										    
    // Set number of commands for this module					    
    m->ncommands = sizeof(mod_coronagraphs_cmds) / sizeof(struct lxrcmd) ;		    
										    
    // Set command control block
    m->cmds = mod_coronagraphs_cmds;
		 								    

    // set module-control-block index in every command-control-block 
    for( i=0; i<m->ncommands; ++i ) {
        mod_coronagraphs_cmds[i].cdata_01 = m->module_number; 
    }


    strncpy( m->name, "coronagraphs", MD_MAX_MODULE_NAME_CHARS );	    
   sprintf(str0, "unknown");
   sprintf(str1, "%s/coronagraphs.o", OBJDIR);
if (!stat(str1, &finfo)) {
   sprintf(str0, "%s", asctime(localtime(&finfo.st_mtime)));}
else { printf("ERROR: cannot find file %s\n",str1);}

   strncpy( m->buildtime, str0, strlen(str0)-1);
   m->buildtime[strlen(str0)] = '\0';

   sprintf(str0,"%s/modules/coronagraphs/coronagraphs.help", SOURCEDIR);
   if((fp=fopen(str0,"r"))==NULL)
   fprintf(stderr,"ERROR: cannot open file %s\n",str0);
   else {
   fgets(str1,999,fp);
   strncpy( m->info, str1, strlen(str1)-1);
   m->info[strlen(str1)] = '\0';
   fclose(fp);
}



    return(PASS);
}






/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_cormk2Dprolate : command function for user command "cormk2Dprolate"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_cormk2Dprolate ( struct lxrcmd *c ) 
{
    coronagraph_make_2Dprolate ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   
    c->results = NULL;   // NULL results block. Change as appropriate
    return(PASS);
}

PF mod_coronagraphs_cormk2Dprolateld ( struct lxrcmd *c ) 
{
  coronagraph_make_2Dprolateld ( c->args[1].v.f, c->args[2].v.f , c->args[3].v.s);   
    c->results = NULL;   // NULL results block. Change as appropriate
    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corup2Dprolate : command function for user command "corup2Dprolate"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corup2Dprolate ( struct lxrcmd *c ) 
{

  coronagraph_update_2Dprolate ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.f);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


PF mod_coronagraphs_cormk2Dprolate_CS ( struct lxrcmd *c ) 
{
  coronagraph_make_2Dprolate_CS ( c->args[1].v.f, c->args[2].v.s);   
  c->results = NULL;   // NULL results block. Change as appropriate
  
  return(PASS);
}

/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimaic : command function for user command "corsimaic"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimaic ( struct lxrcmd *c ) 
{

    coronagraph_simul_AIC ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsim4qpm : command function for user command "corsim4qpm"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsim4qpm ( struct lxrcmd *c ) 
{

    coronagraph_simul_4QPM ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimbl8 : command function for user command "corsimbl8"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimbl8 ( struct lxrcmd *c ) 
{

    coronagraph_simul_BL8 ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimrrpm : command function for user command "corsimrrpm"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimrrpm ( struct lxrcmd *c ) 
{

    coronagraph_simul_RRPM ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimcpa : command function for user command "corsimcpa"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimcpa ( struct lxrcmd *c ) 
{

    coronagraph_simul_CPA ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_mksubpup : command function for user command "mksubpup"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_mksubpup ( struct lxrcmd *c ) 
{

    coronagraphs_make_SUBARU_pupil ( );   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimpiaa : command function for user command "corsimpiaa"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimpiaa ( struct lxrcmd *c ) 
{

    coronagraph_simul_PIAA ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimpiaac : command function for user command "corsimpiaac"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimpiaac ( struct lxrcmd *c ) 
{

    coronagraph_simul_PIAAC ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimpiaacaic : command function for user command "corsimpiaacaic"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimpiaacaic ( struct lxrcmd *c ) 
{

    coronagraph_simul_AIC_PIAAC ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimaplc : command function for user command "corsimaplc"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimaplc ( struct lxrcmd *c ) 
{

    coronagraph_simul_MULTISTEP_APLC ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimpsf : command function for user command "corsimpsf"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimpsf ( struct lxrcmd *c ) 
{

  coronagraph_simulPSF ( c->args[1].v.f, c->args[2].v.f, c->args[3].v.s, c->args[4].v.ld, c->options);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corsimtransm : command function for user command "corsimtransm"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corsimtransm ( struct lxrcmd *c ) 
{

    coronagraph_transm ( c->args[1].v.s, c->args[2].v.ld, c->args[3].v.f, c->options);   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


/*^-----------------------------------------------------------------------------
| PF
| mod_coronagraphs_corouserf : command function for user command "corouserf"
|   struct lxrcmd *c: Command-control-block calling this command function.
|
|
|
| Returns : PASS/FAIL
+-----------------------------------------------------------------------------*/
PF mod_coronagraphs_corouserf ( struct lxrcmd *c ) 
{

    coronagraph_userfunc ( );   


    c->results = NULL;   // NULL results block. Change as appropriate


    return(PASS);
}


PF mod_coronagraphs_scanPIAACMC_centObs_perf ( struct lxrcmd *c ) 
{

  CORONAGRAPHS_scanPIAACMC_centObs_perf ( (float) c->args[1].v.f );   

  c->results = NULL;   // NULL results block. Change as appropriate

  return(PASS);
}
