#ifndef SF_HELPERS
#define SF_HELPERS

#include <time.h>
#include <string.h>
#include <inttypes.h>
#include "stplugin.h"
#include "sf_printf.h"

void gf_count_range(
    uint32_t *groupid,
    uint32_t *index,
    uint32_t *skip,
    uint32_t nobs,
    uint32_t K,
    uint32_t *ncommon
);

uint32_t gf_count_overlap(
    uint32_t *groupid,
    uint32_t *index,
    uint32_t *skip,
    uint32_t nobs,
    uint32_t K,
    uint32_t *ncommon 
);

ST_retcode sf_oom_error(char const *step_desc, char const *obj_desc);

void sf_running_timer (clock_t *timer, const char *msg);

#endif
