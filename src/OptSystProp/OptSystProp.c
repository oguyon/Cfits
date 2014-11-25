#include <fitsio.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate/WFpropagate.h"
#include "statistic/statistic.h"
#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"

#include "OptSystProp/OptSystProp.h"

extern DATA data;

#define SBUFFERSIZE 2000


// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//



int init_OptSystProp()
{
    strcpy(data.module[data.NBmodule].name, __FILE__);
    strcpy(data.module[data.NBmodule].info, "Optical propagation through system");
    data.NBmodule++;


    // add atexit functions here
   

    return 0;

}










int OptSystProp_propagateCube(OPTSYST *optsyst, long index, char *IDin_amp_name, char *IDin_pha_name, char *IDout_amp_name, char *IDout_pha_name, double zprop)
{
    int kl;
    long ii;
    long size;
    long size2;
    long IDin_amp, IDin_pha;
    long IDc_in, IDc_out;
    long IDout_amp, IDout_pha;

    double amp, pha, re, im;



    printf("propagating by %lf m\n", zprop);

    IDin_amp = image_ID(IDin_amp_name);
    IDin_pha = image_ID(IDin_pha_name);
    size = data.image[IDin_amp].md[0].size[0];
    size2 = size*size;
    IDc_in = create_2DCimage_ID("tmppropCin", size, size);

    IDout_amp = create_3Dimage_ID(IDout_amp_name, size, size, optsyst[index].nblambda);
    IDout_pha = create_3Dimage_ID(IDout_pha_name, size, size, optsyst[index].nblambda);


    for(kl=0; kl<optsyst[index].nblambda; kl++)
    {
        printf("kl = %d / %d  %f\n", kl, optsyst[index].nblambda, optsyst[index].lambdaarray[kl]);
        for(ii=0; ii<size2; ii++)
        {
            amp = data.image[IDin_amp].array.F[kl*size2+ii];
            pha = data.image[IDin_pha].array.F[kl*size2+ii];
            data.image[IDc_in].array.CF[ii].re = amp*cos(pha);
            data.image[IDc_in].array.CF[ii].im = amp*sin(pha);
        }
        Fresnel_propagate_wavefront("tmppropCin", "tmppropCout", optsyst[index].pixscale, zprop, optsyst[index].lambdaarray[kl]);

        IDc_out = image_ID("tmppropCout");
        for(ii=0; ii<size2; ii++)
        {
            re = data.image[IDc_out].array.CF[ii].re;
            im = data.image[IDc_out].array.CF[ii].im;
            amp = sqrt(re*re+im*im);
            pha = atan2(im,re);
            data.image[IDout_amp].array.F[kl*size2+ii] = amp;
            data.image[IDout_pha].array.F[kl*size2+ii] = pha;
        }
        delete_image_ID("tmppropCout");
    }


    return 0;
}


/// @param[in]  index	system index (usually 0)
/// @param[in]	elemstart	starting element index
/// @param[in]	elemend		ending element index
/// @param[in]  savedir  directory to which image results are saved
///
/// *elemkeepmem	1 if element complex amplitude should be kept in memory after use
///

int OptSystProp_run(OPTSYST *optsyst, long index, long elemstart, long elemend, char *savedir)
{
    char command[500];
    char imname[200];
    char imnameamp[200];
    char imnamepha[200];

    long IDx, IDy, IDr, IDPA;
    double x, y;
    long IDa, IDp;
    long size;
    long nblambda;
    long size2;
    long ii, jj, k;
    long elem;
    long kl;

    char fname_pupa0[200];
    char fname_pupp0[200];
    char fname[200];

    long ID;
    double proplim = 1.0e-4;
    double total;


    float beamradpix;
    long ID0, ID1, ID2;
    long size0, size1;
    long i, j;

    char imnameamp_in[200];
    char imnamepha_in[200];
    char imnameamp_out[200];
    char imnamepha_out[200];
    long emax;
    double propdist;

    long gsize, offset, ii1, jj1;
    float re, im;
    long IDre, IDim, IDre1, IDim1;
    float *convkern;
    float val, tot;
    float u, t;

    long elemstart1 = 0;
    int elemOK;

    int r;

    size = optsyst[index].size;
    size2 = size*size;
    nblambda = optsyst[index].nblambda;



    // create base complex amplitude
    sprintf(imname,"WFamp%ld", index);
    IDa = image_ID(imname);
    if(IDa==-1)
        IDa = create_3Dimage_ID(imname, size, size, nblambda);

    sprintf(imname,"WFpha%ld", index);
    IDp = image_ID(imname);
    if(IDp==-1)
        IDp = create_3Dimage_ID(imname, size, size, nblambda);

    for(ii=0; ii<size2; ii++)
        for(kl=0; kl<nblambda; kl++)
            data.image[IDa].array.F[size2*kl+ii] = 1.0;



    elemstart1 = 0;
    elemOK = 1;
    while(elemOK==1)
    {
        if(elemstart1==0)
        {
            sprintf(imnameamp_in, "WFamp%ld", index);
            sprintf(imnamepha_in, "WFpha%ld", index);
        }
        else
        {
            sprintf(imnameamp_in, "WFamp%ld_%03ld", index, elemstart1-1);
            sprintf(imnamepha_in, "WFpha%ld_%03ld", index, elemstart1-1);
        }
        if(((ID1=image_ID(imnameamp_in))!=-1)&&((ID2=image_ID(imnamepha_in))!=-1)&&(elemstart1<elemstart+1))
        {
            elemstart1++;
            elemOK = 1;
        }
        else
            elemOK = 0;

        printf("%ld/%ld %d    %s %ld   %s %ld\n", elemstart1, elemstart, elemOK, imnameamp_in, ID1, imnamepha_in, ID2);
    }
    elemstart1--;

    printf("STARTING AT ELEMENT %ld\n", elemstart1);

    emax = elemend;
    if(emax>optsyst[index].NBelem)
        emax = optsyst[index].NBelem;

    for(elem=elemstart1; elem<emax; elem++)
    {
        if(elem==0)
        {
            sprintf(imnameamp_in, "WFamp%ld", index);
            sprintf(imnamepha_in, "WFpha%ld", index);
            sprintf(imnameamp_out, "WFamp%ld_000", index);
            sprintf(imnamepha_out, "WFpha%ld_000", index);
            propdist = optsyst[index].elemZpos[0];
        }
        else
        {
            sprintf(imnameamp_in, "WFamp%ld_%03ld", index, elem-1);
            sprintf(imnamepha_in, "WFpha%ld_%03ld", index, elem-1);
            sprintf(imnameamp_out, "WFamp%ld_%03ld", index, elem);
            sprintf(imnamepha_out, "WFpha%ld_%03ld", index, elem);
            propdist = optsyst[index].elemZpos[elem]-optsyst[index].elemZpos[elem-1];
        }


        if(image_ID(imnameamp_out)!=-1)
            delete_image_ID(imnameamp_out);

        if(image_ID(imnamepha_out)!=-1)
            delete_image_ID(imnamepha_out);


        if( fabs(propdist)>proplim )
        {
            printf("Propagating to element %ld  (%lf m)\n", elem,  propdist);
            OptSystProp_propagateCube(optsyst, 0, imnameamp_in, imnamepha_in, imnameamp_out, imnamepha_out, propdist);
        }
        else
        {
            copy_image_ID(imnameamp_in, imnameamp_out);
            copy_image_ID(imnamepha_in, imnamepha_out);
        }
        IDa = image_ID(imnameamp_out);
        IDp = image_ID(imnamepha_out);

        /// discard element memory after used
        printf("*********** %ld  -> %d\n", elem-1, optsyst[index].keepMem[elem-1]);
        if(optsyst[index].keepMem[elem-1]==0)
        {
            printf("********** Deleting element %ld      %s %s\n", elem-1, imnameamp_in, imnamepha_in);
            delete_image_ID(imnameamp_in);
            delete_image_ID(imnamepha_in);
        }

        printf("Applying element %ld\n", elem);
        fflush(stdout);

        if(optsyst[index].elemtype[elem]==1)   // OPAQUE MASK
        {
            ID = optsyst[index].elemarrayindex[elem];
            printf("============= elem %ld:  Opaque mask (%s) =================\n", elem, data.image[ID].md[0].name);
            fflush(stdout);

            if(ID == -1)
            {
                printf("ERROR: ID = -1, missing mask image\n");
                exit(0);
            }

            //	save_fits(data.image[ID].md[0].name, "!opmask.fits"); //TEST
            //save_fits(data.image[IDa].md[0].name, "!opmask1.fits"); //TEST

            if(data.image[ID].md[0].size[2] != nblambda)
            {
# ifdef HAVE_LIBGOMP
                #pragma omp parallel default(shared) private(ii)
                {
                    #pragma omp for
# endif
                    for(kl=0; kl<nblambda; kl++)
                        for(ii=0; ii<size2; ii++)
                            data.image[IDa].array.F[size2*kl+ii] *= data.image[ID].array.F[ii];
# ifdef HAVE_LIBGOMP
                }
# endif
            }
            else
            {
# ifdef HAVE_LIBGOMP
                #pragma omp parallel
                {
                    #pragma omp for
# endif
                    for(ii=0; ii<size2*nblambda; ii++)
                        data.image[IDa].array.F[ii] *= data.image[ID].array.F[ii];
# ifdef HAVE_LIBGOMP
                }
# endif
            }

            //	save_fits(data.image[IDa].md[0].name, "!opmask2.fits"); //TEST


        }

        if(optsyst[index].elemtype[elem]==3)  // MIRROR OR REFRACTIVE SURFACE
        {
            printf("============= Mirror or Refractive surface =======================\n");
            fflush(stdout);
            ID = optsyst[index].ASPHSURFMarray[optsyst[index].elemarrayindex[elem]].surfID;
            printf("%d surface ID = %ld\n", optsyst[index].elemarrayindex[elem], ID);
            fflush(stdout);
# ifdef HAVE_LIBGOMP
            #pragma omp parallel default(shared) private(ii)
            {
                #pragma omp for
# endif
                for(kl=0; kl<nblambda; kl++)
                    for(ii=0; ii<size2; ii++)
                        data.image[IDp].array.F[size2*kl+ii] -= 4.0*M_PI*data.image[ID].array.F[ii]/optsyst[index].lambdaarray[kl];
# ifdef HAVE_LIBGOMP
            }
# endif
        }


        if(optsyst[index].elemtype[elem]==5)  // FOCAL PLANE MASK - MASK INPUT IS 1-MASK FOR EFFICIENT DFT
        {
            printf("============= Focal Plane Mask ==============\n");
            fflush(stdout);
            // uses 1-fpm

            // test
            //	  save_fits(imnameamp_out, "!TESTamp.fits");
            //exit(0);

            ID = mk_complex_from_amph(imnameamp_out, imnamepha_out, "_WFctmp");
            delete_image_ID(imnameamp_out);
            delete_image_ID(imnamepha_out);


            if(optsyst[index].DFTgridpad>0)
            {
                // RESAMPLE ON A SPARSE GRID TO SPEED UP DFT
                IDre = create_3Dimage_ID("dftgridre", size, size, nblambda);
                IDim = create_3Dimage_ID("dftgridim", size, size, nblambda);
                gsize = 2*optsyst[index].DFTgridpad+1; // grid size, odd number - this is the space between pixels
                offset = optsyst[index].DFTgridpad; // offset from box edge to active pixel

                ID = image_ID("_WFctmp");
                for(kl=0; kl<nblambda; kl++)
                    for(ii=0; ii<size; ii++)
                        for(jj=0; jj<size; jj++)
                        {
                            re = data.image[ID].array.CF[size2*kl+jj*size+ii].re;
                            im = data.image[ID].array.CF[size2*kl+jj*size+ii].im;
                            ii1 = offset + ((long) (ii/gsize))*gsize;
                            jj1 = offset + ((long) (jj/gsize))*gsize;
                            if((ii1<size)&&(jj1<size))
                            {
                                data.image[IDre].array.F[size2*kl+jj1*size+ii1] += re;
                                data.image[IDim].array.F[size2*kl+jj1*size+ii1] += im;
                            }
                        }
                // save_fits("dftgridre", "!dftgridre.fits");
                // save_fits("dftgridim", "!dftgridim.fits");
                mk_complex_from_reim("dftgridre", "dftgridim", "_WFctmpc");
                delete_image_ID("dftgridre");
                delete_image_ID("dftgridim");


                i = optsyst[index].elemarrayindex[elem];
                ID = optsyst[index].FOCMASKarray[i].fpmID;
                printf("focm : %s\n", data.image[ID].md[0].name);
                //	      printf("Saving to testfpm.fits\n");
                //      fflush(stdout);
                //	      list_image_ID();
                mk_amph_from_complex("piaacmcfpm", "fpma", "fpmp");

                if(optsyst[index].SAVE == 1)
                {
                    sprintf(fname, "!%s/fpm_ampl.fits", savedir);
                    save_fits("fpma", fname);
                    sprintf(fname, "!%s/fpm_pha.fits", savedir);
                    save_fits("fpmp", fname);
                }
                //	      exit(0);
                fft_DFTinsertFPM("_WFctmpc", data.image[ID].md[0].name, optsyst[index].FOCMASKarray[i].zfactor, "_WFcout");
                delete_image_ID("_WFctmpc");
                sprintf(command, "mv _DFT_foca %s/_DFT_foca_%02ld.fits", savedir, elem);
                r = system(command);
                sprintf(command, "mv _DFT_focp %s/_DFT_focp_%02ld.fits", savedir, elem);
                r = system(command);

                //mk_reim_from_complex("_WFcout", "_twfre", "_twfim");
                //save_fits("_twfre", "!_twfre.fits");
                //save_fits("_twfim", "!_twfim.fits");
                //
                // INTERPOLATE SPARSE RESULT ON CONTINUOUS GRID
                //
                convkern = (float*) malloc(sizeof(float)*(2*gsize+1)*(2*gsize+1));
                tot = 0.0;
                for(i=0; i<2*gsize+1; i++)
                    for(j=0; j<2*gsize+1; j++)
                    {
                        u = fabs(1.0*(i-gsize)/gsize);
                        t = fabs(1.0*(j-gsize)/gsize);
                        val = (1.0-u)*(1.0-t);
                        convkern[j*(2*gsize+1)+i] = val;
                        //printf("   %d %d %f\n", i, j, val);
                        tot += val;
                    }
                for(i=0; i<(2*gsize+1)*(2*gsize+1); i++)
                    convkern[i] *= gsize*gsize/tot;


                ID = image_ID("_WFcout");
                IDre1 = create_3Dimage_ID("dftgridre1", size, size, nblambda);
                IDim1 = create_3Dimage_ID("dftgridim1", size, size, nblambda);
                for(kl=0; kl<nblambda; kl++)
                    for(ii1=offset+gsize; ii1<size-gsize; ii1+=gsize)
                        for(jj1=offset+gsize; jj1<size-gsize; jj1+=gsize)
                        {
                            re = data.image[ID].array.CF[size2*kl+jj1*size+ii1].re;
                            im = data.image[ID].array.CF[size2*kl+jj1*size+ii1].im;

                            for(i=0; i<2*gsize+1; i++)
                                for(j=0; j<2*gsize+1; j++)
                                {
                                    ii = ii1+(i-gsize);
                                    jj = jj1+(j-gsize);
                                    data.image[IDre1].array.F[size2*kl+jj*size+ii] += re*convkern[j*(2*gsize+1)+i];
                                    data.image[IDim1].array.F[size2*kl+jj*size+ii] += im*convkern[j*(2*gsize+1)+i];
                                }
                        }
                //save_fits("dftgridre1", "!dftgridre1.fits");
                //save_fits("dftgridim1", "!dftgridim1.fits");
                free(convkern);
                delete_image_ID("_WFcout");
                mk_complex_from_reim("dftgridre1", "dftgridim1", "_WFcout");
                delete_image_ID("dftgridre1");
                delete_image_ID("dftgridim1");
            }
            else
            {
                i = optsyst[index].elemarrayindex[elem];
                ID = optsyst[index].FOCMASKarray[i].fpmID;
                printf("focm : %s\n", data.image[ID].md[0].name);
                fflush(stdout);
                fft_DFTinsertFPM("_WFctmp", data.image[ID].md[0].name, optsyst[index].FOCMASKarray[i].zfactor, "_WFcout");
                sprintf(command, "mv _DFT_foca %s/_DFT_foca_%02ld.fits", savedir, elem);
                r = system(command);
                sprintf(command, "mv _DFT_focp %s/_DFT_focp_%02ld.fits", savedir, elem);
                r = system(command);
            }

            i = optsyst[index].elemarrayindex[elem];

            if(optsyst[index].FOCMASKarray[i].mode == 1)
            {
                arith_image_sub_inplace("_WFctmp", "_WFcout");
                mk_amph_from_complex("_WFctmp", imnameamp_out, imnamepha_out);
            }
            else
                mk_amph_from_complex("_WFcout", imnameamp_out, imnamepha_out);


            delete_image_ID("_WFctmp");
            delete_image_ID("_WFcout");
            //  delete_image_ID("dftgrid");
        }

        IDa = image_ID(imnameamp_out);
        optsyst[index].flux[elem] = 0.0;
        for(kl=0; kl<nblambda; kl++)
            for(ii=0; ii<size2; ii++)
                optsyst[index].flux[elem] += data.image[IDa].array.F[kl*size2+ii]*data.image[IDa].array.F[kl*size2+ii];

        printf("Element %ld  [%ld %ld]  Flux = %lf\n", elem, nblambda, size2, optsyst[index].flux[elem]/nblambda);
        if(isnan(optsyst[index].flux[elem])!=0)
            exit(0);

        if(optsyst[index].SAVE == 1)
        {
            printf("Saving intermediate plane [%ld] ... ", elem);
            fflush(stdout);

            sprintf(fname, "!./%s/WFamp%ld_%03ld.fits", savedir, index, elem);
            save_fits(imnameamp_out, fname);
            sprintf(fname, "!./%s/WFpha%ld_%03ld.fits", savedir, index, elem);
            save_fits(imnamepha_out, fname);


            printf("done\n");
            fflush(stdout);
        }
    }

    if((elem==optsyst[index].NBelem)&&(optsyst[index].endmode==0)) // Compute final focal plane image
    {
        printf("COMPUTING FINAL IMAGE AS FFT OF %ld\n", elem-1);
        mk_complex_from_amph(imnameamp_out, imnamepha_out, "_WFctmp");
        permut("_WFctmp");
        sprintf(imname, "psfc%ld", index);
        do2dfft("_WFctmp", imname);
        delete_image_ID("_WFctmp");
        permut(imname);
        sprintf(imnameamp, "psfa%ld", index);
        sprintf(imnamepha, "psfp%ld", index);
        mk_amph_from_complex(imname, imnameamp, imnamepha);

        if(optsyst[index].SAVE == 1)
        {   sprintf(fname, "!%s/psfa%ld.fits", savedir, index);
            save_fits(imnameamp, fname);

            sprintf(fname, "!%s/psfp%ld.fits", savedir, index);
            save_fits(imnamepha, fname);
        }


        ID = image_ID(imnameamp);
        for(ii=0; ii<size2*nblambda; ii++)
            data.image[ID].array.F[ii] /= sqrt(size2*optsyst[index].flux[0]/nblambda);

        sprintf(imname, "psfi%ld", index);
        arith_image_mult(imnameamp, imnameamp, imname);
        total = arith_image_total(imname)/nblambda;
        printf("TOTAL = %lf\n", total);

        if(optsyst[index].SAVE == 1)
        {   sprintf(fname, "!%s/psfi%ld.fits", savedir, index);
            save_fits(imname, fname);
        }
    }


    return(0);
}


