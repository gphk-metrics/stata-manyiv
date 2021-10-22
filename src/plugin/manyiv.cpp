/*********************************************************************
 * Program: manyiv_d_projection.cpp
 * Created: Tue Oct 19 18:34:49 EDT 2021
 * Updated: Thu Oct 21 23:25:47 EDT 2021
 * Purpose: Stata plugin to compute the diagonal of the projection
 *          matrix of the design matrix implied by a set of FE. Given
 *          the pattern of such a design matrix leverading the sparse
 *          structure hopefully minimizes the memory requirements.
 * Version: 0.1.0
 *********************************************************************/

#include "manyiv.h"
#include "stplugin.h"
#include "sf_printf.c"
#include "sf_helpers.c"
#include "sf_d_projection.cpp"

int main()
{
    return(0);
}

int WinMain()
{
    return(0);
}

// Syntax
//     plugin call manyiv, file [benchmark]

STDLL stata_call(int argc, char * argv[])
{
    ST_retcode rc;
    uint32_t benchmark;
    if ( strcmp(argv[0], "_plugin_check") == 0 ) {
        sf_printf("(note: manyiv_plugin v" MANYIV_VERSION " successfully loaded)\n");
        return(0);
    }
    else {
        MANYIV_CHAR(fname, strlen(argv[0]) + 1);
        strcpy(fname, argv[0]);
        benchmark = strcmp(argv[1], "_plugin_bench") == 0;
        return(sf_d_projection(fname, benchmark));
    }
}
