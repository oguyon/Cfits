#ifndef _COREMODMEMORY_H
#define _COREMODMEMORY_H

/* the number of images in the data structure is kept NB_IMAGES_BUFFER above the number of used images prior to the execution of any function. It means that no function should create more than 100 images. */
#define NB_IMAGES_BUFFER 500
/* when the number of free images in the data structure is below NB_IMAGES_BUFFER, it is increased by  NB_IMAGES_BUFFER */
#define NB_IMAGES_BUFFER_REALLOC 600

/* the number of variables in the data structure is kept NB_VARIABLES_BUFFER above the number of used variables prior to the execution of any function. It means that no function should create more than 100 variables. */
#define NB_VARIABLES_BUFFER 100
/* when the number of free variables in the data structure is below NB_VARIABLES_BUFFER, it is increased by  NB_VARIABLES_BUFFER */
#define NB_VARIABLES_BUFFER_REALLOC 150


/*void print_sys_mem_info();*/


int memory_monitor(char *termttyname);





long compute_nb_image();

long compute_nb_variable();

long long compute_image_memory();

long compute_variable_memory();

long image_ID(char *name);

long image_ID_noaccessupdate(char *name);

long variable_ID(char *name);

long next_avail_image_ID();

long next_avail_variable_ID();

int delete_image_ID(char* imname);

int delete_image_ID_prefix(char *prefix);

int delete_variable_ID(char* varname);

long create_image_ID(char *name, long naxis, long *size, int atype, int shared, int nbkw);

long create_1Dimage_ID(char *ID_name, long xsize);

long create_1DCimage_ID(char *ID_name, long xsize);

long create_2Dimage_ID(char *ID_name, long xsize, long ysize);

long create_2Dimagedouble_ID(char *ID_name, long xsize, long ysize);

long create_2DCimage_ID(char *ID_name, long xsize, long ysize);

long create_3Dimage_ID(char *ID_name, long xsize, long ysize, long zsize);

long create_3Dimage_ID_double(char *ID_name, long xsize, long ysize, long zsize);

long create_3DCimage_ID(char *ID_name, long xsize, long ysize, long zsize);

long copy_image_ID(char *name, char *newname);

long create_variable_ID(char *name, double value);

int list_image_ID_ofp(FILE *fo);

int list_image_ID();

int list_image_ID_file(char *fname);

int list_variable_ID();

long chname_image_ID(char *ID_name, char *new_name);

int mk_complex_from_reim(char *re_name, char *im_name, char *out_name);

int mk_complex_from_amph(char *am_name, char *ph_name, char *out_name);

int mk_reim_from_complex(char *in_name, char *re_name, char *im_name);

int mk_amph_from_complex(char *in_name, char *am_name, char *ph_name);

int mk_reim_from_amph(char *am_name, char *ph_name, char *re_out_name, char *im_out_name);

int mk_amph_from_reim(char *re_name, char *im_name, char *am_out_name, char *ph_out_name);

int clearall();

int check_2Dsize(char *ID_name, long xsize, long ysize);

int check_3Dsize(char *ID_name, long xsize, long ysize, long zsize);

int rotate_cube(char *ID_name, char *ID_out_name, int orientation);

#endif
