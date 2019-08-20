#include <stdlib.h>
#include "mgpriv.h"
#include "ksort.h"

int mg_verbose = 1;
int mg_dbg_flag = 0;
double mg_realtime0;

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mg128_t, sort_key_128x, 8) 

KSORT_INIT_GENERIC(uint32_t)
