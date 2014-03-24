#ifndef _00CORE_H
#define _00CORE_H


int printRED(char *string);

int printWARNING(const char *file, const char *func, int line, char *warnmessage);

int printERROR(const char *file, const char *func, int line, char *errmessage);

int setCfits_precision(int vp);

int CfitsWritePid();

struct timespec timespec_diff(struct timespec start, struct timespec end);
double timespec_diff_double(struct timespec start, struct timespec end);

#endif
