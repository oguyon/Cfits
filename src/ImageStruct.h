/**
 * @file    ImageStruct.h
 * @brief   Image structure definition
 * 
 * The IMAGE structure is defined here
 * Supports shared memory, low latency IPC through semaphores
 * 
 * @author  O. Guyon
 * @date    24 Jun 2017
 *
 * @bug No known bugs. 
 * 
 */

#ifndef _IMAGESTRUCT_H
#define _IMAGESTRUCT_H

#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <semaphore.h>


#define SHAREDMEMDIR "/tmp"        /**< loaction of file mapped semaphores */


#define SEMAPHORE_MAX       1000   /**< maximum number of semaphores */
#define SEMAPHORE_MAXVAL    10 	   /**< maximum value for each of the semaphore, mitigates warm-up time when processes catch up with data that has accumulated */







typedef struct
{
    float re;
    float im;
} complex_float;


typedef struct
{
    double re;
    double im;
} complex_double;




/// Data types

#define CHAR             1
#define INT              2
#define FLOAT            3
#define DOUBLE           4
#define COMPLEX_FLOAT    5
#define COMPLEX_DOUBLE   6
#define USHORT           7
#define LONG             8

/// default data type for float and complex float
#define Dtype            3  
#define CDtype           5





/** The IMAGE structure includes :
 *   - an array of IMAGE_KEWORD structures
 *   - an array of IMAGE_METADATA structures (usually only 1 element)
 */


typedef struct
{
    char name[16];         /**< keyword name */
    char type;             /**< N: unused, L: long, D: double, S: 16-char string */

    union {
        long numl;
        double numf;
        char valstr[16];
    } value;

    char comment[80];

} IMAGE_KEYWORD;




typedef struct
{
    char name[80];              /**< image name */

    long naxis;                 /**< number of axis */
    long size[3];               /**< image size */
    long nelement;              /**< number of elements in image */
    int atype;                  /**< data type code */

    double creation_time;       /**< creation time (since program start) */
    double last_access;         /**< last time the image was accessed  (since program start) */
    struct timespec wtime;

    int shared;                 /**< 1 if in shared memory */

    int write;               	/**< 1 if image is being written */
    int status;              	/**< 1 to log image (default); 0 : do not log: 2 : stop log (then goes back to 2) */
    long cnt0;               	/**< counter (incremented if image is updated) */
    long cnt1;               	/**< in 3D rolling buffer image, this is the last slice written */

    long NBkw;                  /**< number of keywords */

} IMAGE_METADATA;






typedef struct          		/**< structure used to store data arrays */
{
    int used;
    int shmfd;		 			/**< if shared memory, file descriptor */
    size_t memsize; 			/**< total size in memory if shared */

    IMAGE_METADATA *md;			

    union
    {
        char *C;
        int *I;
        long *L;
        float *F;
        double *D;
        complex_float *CF;
        complex_double *CD;
        unsigned short int *U;
    } array;                 	/**< pointer to data array */

    IMAGE_KEYWORD *kw;

    int sem; 					/**< number of semaphores in use     */
    sem_t **semptr; 			/**< semaphore array */

    sem_t *semlog; 				/**< semaphore for logging */
    
    
    char name[80]; 				/**< local name (can be different from name in shared memory) */

} IMAGE;




#endif
