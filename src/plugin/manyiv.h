#ifndef MANYIV
#define MANYIV

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#define MANYIV_VERSION "0.1.0"

#define MANYIV_CHAR(cvar, len)    \
    char *cvar = new char[len]; \
    memset (cvar, '\0', sizeof(char) * len)

#define MANYIV_PWMAX(a, b) ( (a) > (b) ? (a) : (b) )
#define MANYIV_PWMIN(a, b) ( (a) > (b) ? (b) : (a) )

#endif
