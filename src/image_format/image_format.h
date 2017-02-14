#ifndef _IMAGEFORMATMODULE_H
#define _IMAGEFORMATMODULE_H


int_fast8_t init_image_format();


int IMAGE_FORMAT_im_to_ASCII(const char *IDname, const char *foutname);

int IMAGE_FORMAT_FITS_to_ASCII(const char *finname, const char *foutname);

int image_writeBMP_auto(const char *IDnameR, const char *IDnameG, const char *IDnameB, const char *outname);

int CR2toFITS(const char *fnameCR2, const char *fnameFITS);

long IMAGE_FORMAT_FITS_to_ushortintbin_lock( const char *IDname, const char *fname);

long IMAGE_FORMAT_FITS_to_floatbin_lock(  const char *IDname, const char *fname);

long IMAGE_FORMAT_read_binary32f(const char *fname, long xsize, long ysize, const char *IDname);

int loadCR2toFITSRGB(const char *fnameCR2, const char *fnameFITSr, const char *fnameFITSg, const char *fnameFITSb);

int image_format_extract_RGGBchan(const char *ID_name, const char *IDoutR_name, const char *IDoutG1_name, const char *IDoutG2_name, const char *IDoutB_name);

#endif

