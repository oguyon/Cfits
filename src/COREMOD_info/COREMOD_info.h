#if !defined(INFO_H)
#define INFO_H

long brighter(char *ID_name, PRECISION value);
/* number of pixels brighter than value */

int img_nbpix_flux(char *ID_name);

PRECISION img_percentile(char *ID_name, PRECISION p);

int img_histoc(char *ID_name, char *fname);

int make_histogram(char *ID_name, char *ID_out_name, PRECISION min, PRECISION max, long nbsteps);

PRECISION ssquare(char *ID_name);

PRECISION rms_dev(char *ID_name);

int stats(char *ID_name, char *options);

PRECISION img_min(char *ID_name);

PRECISION img_max(char *ID_name);

int profile(char *ID_name, char *outfile, PRECISION xcenter, PRECISION ycenter, PRECISION step, long nb_step);

int profile2im(char *profile_name, long nbpoints, long size, PRECISION xcenter, PRECISION ycenter, PRECISION radius, char *out);

int printpix(char *ID_name, char *filename);

PRECISION background_photon_noise(char *ID_name);

int test_structure_function(char *ID_name, long NBpoints, char *fname);

int full_structure_function(char *ID_name, long NBpoints, char *ID_out);

int fft_structure_function(char *ID_in, char *ID_out);

#endif
