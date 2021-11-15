/*********************************************************************
 * Program: manyiv_d_projection.cpp
 * Support: Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com>
 * Created: Tue Oct 19 18:34:49 EDT 2021
 * Updated: Sun Nov 14 16:20:43 EST 2021
 * Purpose: Stata plugin to compute the diagonal of the projection
 *          matrix of the design matrix implied by a set of FE. Given
 *          the pattern of such a design matrix leverading the sparse
 *          structure hopefully minimizes the memory requirements.
 * Version: 0.2.0
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
// Errors
//     1701 - Matrix inversion failed
//     1702 - OOM
//     1703 - Unable to export results back to Stata

STDLL stata_call(int argc, char * argv[])
{
    uint32_t benchmark;
    if ( argc ) {
        if ( strcmp(argv[0], "_plugin_check") == 0 ) {
            sf_printf("(note: manyiv_plugin v" MANYIV_VERSION " successfully loaded)\n");
            return(0);
        }
        else if ( strcmp(argv[0], "_plugin_run") == 0 ) {
            MANYIV_CHAR(fname, strlen(argv[1]) + 1);
            strcpy(fname, argv[1]);
            benchmark = strcmp(argv[2], "_plugin_bench") == 0;
            return(sf_d_projection(fname, benchmark));
        }
    }
    return(0);
}
