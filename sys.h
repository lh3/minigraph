#ifndef MG_SYS_H
#define MG_SYS_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void mg_err_fputs(const char *str, FILE *fp);

double realtime(void);
double cputime(void);
long peakrss(void);

#ifdef __cplusplus
}
#endif

#endif
