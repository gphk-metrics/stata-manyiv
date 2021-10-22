#include "sf_printf.h"
#include "sf_helpers.h"

void gf_count_range(
    uint32_t *groupid,
    uint32_t *index,
    uint32_t *skip,
    uint32_t nobs,
    uint32_t K,
    uint32_t *ncommon)
{
    uint32_t i;
    memset(ncommon, '\0', K * sizeof(uint32_t));
    for (i = 0; i < nobs; i++) {
        ncommon[groupid[index[i]]]++;
    }
    for (i = 0; i < K; i++) {
        if ( skip[i] ) {
            ncommon[i] = 0;
        }
    }
}

uint32_t gf_count_overlap(
    uint32_t *groupid,
    uint32_t *index,
    uint32_t *skip,
    uint32_t nobs,
    uint32_t K,
    uint32_t *ncommon)
{
    uint32_t i, nnz = 0;
    memset(ncommon, '\0', K * sizeof(uint32_t));
    for (i = 0; i < nobs; i++) {
        ncommon[groupid[index[i]]] = 1;
    }
    for (i = 0; i < K; i++) {
        if ( skip[i] ) ncommon[i] = 0;
        nnz += ncommon[i];
    }
    return(nnz);
}

ST_retcode sf_oom_error(char const *step_desc, char const *obj_desc)
{
    sf_errprintf("%s: Unable to allocate memory for object '%s'.\n", step_desc, obj_desc);
    sf_printf("See {help gcollapse##memory:help gcollapse (Out of memory)}.\n");
    return (1702);
}

void sf_running_timer (clock_t *timer, const char *msg)
{
    double diff = (double) (clock() - *timer) / CLOCKS_PER_SEC;
    sf_printf (msg);
    sf_printf (" (%.3f seconds).\n", diff);
    *timer = clock();
}
