#ifndef _IMGREDUCE_H
#define _IMGREDUCE_H

int init_img_reduce();

long IMG_REDUCE_cleanbadpix_fast(char *IDname, char *IDbadpix_name, char *IDoutname);

long IMG_REDUCE_cubesimplestat(char *IDin_name);
int IMG_REDUCE_cubeprocess(char *IDin_name);

#endif
