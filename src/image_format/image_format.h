#ifndef _IMAGEFORMATMODULE_H
#define _IMAGEFORMATMODULE_H


int CR2toFITS(char *fnameCR2, char *fnameFITS);

long IMAGE_FORMAT_FITS_to_ushortintbin_lock( char *IDname, char *fname);

long IMAGE_FORMAT_FITS_to_floatbin_lock(  char *IDname, char *fname);

long IMAGE_FORMAT_read_binary16f(char *fname, long xsize, long ysize, char *IDname);

#endif

