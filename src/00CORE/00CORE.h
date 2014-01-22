#ifndef _00CORE_H
#define _00CORE_H


int printRED(char *string);

int printWARNING(const char *file, const char *func, int line, char *warnmessage);

int printERROR(const char *file, const char *func, int line, char *errmessage);

int setCfits_precision(int vp);

int CfitsWritePid();

#endif
