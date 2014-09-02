#ifndef _OPTICSMATERIALS_H
#define _OPTICSMATERIALS_H

double OPTICSMATERIALS_n( int material, double lambda);

// phase offset as function of mask thickness and lambda
double OPTICSMATERIALS_pha_lambda( int material, double z, double lambda );

#endif
