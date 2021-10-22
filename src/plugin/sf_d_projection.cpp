#include "sf_printf.h"
#include "sf_helpers.h"
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

ST_retcode sf_d_projection(char *fname, uint8_t benchmark);
ST_retcode sf_d_projection(char *fname, uint8_t benchmark)
{
    clock_t timer  = clock();
    ST_double maxr, tol, z;
    SimplicialLLT <SparseMatrix<ST_double> > solverLL;
    SparseQR<SparseMatrix<ST_double, RowMajor>, COLAMDOrdering<int32_t> > solveQR;
    uint32_t i, j, k, n, in1, in2, nobs, nabsorb, ktot, kskip, kmax, nnz, kindep;
    uint32_t offseti, offsetj, skipoffi, skipoffj;

    ST_retcode rc = 0;

    /**********************************************************************
     * 1. Read index, info, group
     **********************************************************************/

    fstream fh(fname, ios::in | ios::binary);
    if ( fh.is_open() ) {
        fh.read((char *) (&nabsorb), sizeof(nabsorb));
        fh.read((char *) (&nobs),    sizeof(nobs));

        uint32_t nlevels[nabsorb];
        uint32_t nlevelsindep[nabsorb];
        uint32_t *index[nabsorb];
        uint32_t *info[nabsorb];
        uint32_t *groupid[nabsorb];
        uint32_t *skip[nabsorb];
        uint32_t *groupmask[nabsorb];

        for (k = 0; k < nabsorb; k++) {
            fh.read((char *) (nlevels + k), sizeof(nlevels[k]));
        }

        kmax = ktot = 0;
        for (k = 0; k < nabsorb; k++) {
            ktot += nlevels[k];
            kmax  = (nlevels[k] > kmax)? nlevels[k]: kmax;

            index[k]     = (uint32_t *) calloc(nobs, sizeof(*index[k]));
            groupid[k]   = (uint32_t *) calloc(nobs, sizeof(*groupid[k]));
            info[k]      = (uint32_t *) calloc(nlevels[k] + 1, sizeof(*info[k]));
            skip[k]      = (uint32_t *) calloc(nlevels[k], sizeof(*skip[k]));
            groupmask[k] = (uint32_t *) calloc(nlevels[k], sizeof(*groupmask[k]));

            if ( index[k]     == NULL ) sf_oom_error("sf_d_projection", "index");
            if ( groupid[k]   == NULL ) sf_oom_error("sf_d_projection", "groupid");
            if ( info[k]      == NULL ) sf_oom_error("sf_d_projection", "info");
            if ( skip[k]      == NULL ) sf_oom_error("sf_d_projection", "skip");
            if ( groupmask[k] == NULL ) sf_oom_error("sf_d_projection", "groupmask");
        }

        for (k = 0; k < nabsorb; k++) {
            fh.read((char *) index[k],   sizeof(*index[k])   * nobs);
            fh.read((char *) groupid[k], sizeof(*groupid[k]) * nobs);
            fh.read((char *) info[k],    sizeof(*info[k])    * (nlevels[k] + 1));
            fh.read((char *) skip[k],    sizeof(*skip[k])    * nlevels[k]);
        }
        fh.close();

        kskip = 0;
        for (k = 0; k < nabsorb; k++) {
            nlevelsindep[k] = nlevels[k];
            for (i = 0; i < nlevels[k]; i++) {
                if ( skip[k][i] ) {
                    ktot--;
                    kskip++;
                    nlevelsindep[k]--;
                }
            }
        }

        uint32_t  *ncommon      = (uint32_t *)  calloc(kmax, sizeof(*ncommon));
        uint32_t  *indepix      = (uint32_t *)  calloc(ktot, sizeof(*indepix));
        ST_double *d_projection = (ST_double *) calloc(nobs, sizeof(*d_projection));
        SparseMatrix<ST_double> ZZ(ktot, ktot);

        if ( ncommon      == NULL ) sf_oom_error("sf_d_projection", "ncommon");
        if ( d_projection == NULL ) sf_oom_error("sf_d_projection", "d_projection");
        if ( indepix      == NULL ) sf_oom_error("sf_d_projection", "indepix");

        if ( benchmark )
            sf_running_timer(&timer, "\tPlugin Step 1: Read info and allocate memory");

    /**********************************************************************
     * 2. Compute cross product
     **********************************************************************/

        ZZ.setZero();
        nnz = 0;
        for (k = 0; k < nabsorb; k++) {
            nnz += nlevelsindep[k];
            for (n = k + 1; n < nabsorb; n++) {
                for (i = 0; i < nlevels[k]; i++) {
                    if ( skip[k][i] ) continue;
                    in1 = info[k][i];
                    in2 = info[k][i + 1];
                    nnz += 2 * gf_count_overlap(groupid[n], index[k] + in1, skip[n], in2 - in1, nlevels[n], ncommon);
                }
            }
        }
        ZZ.reserve(nnz);

        skipoffi = skipoffj = offseti = offsetj = 0;
        for (k = 0; k < nabsorb; k++) {
            skipoffi = 0;
            skipoffj = 0;
            offsetj = offseti;
            for (i = 0; i < nlevels[k]; i++) {
                if ( skip[k][i] ) {
                    skipoffi++;
                    skipoffj++;
                }
                else {
                    in1 = info[k][i];
                    in2 = info[k][i + 1];
                    ZZ.insert(i + offseti - skipoffi, i + offsetj - skipoffj) = (ST_double) (in2 - in1);
                }
            }
            for (n = k + 1; n < nabsorb; n++) {
                skipoffi = 0;
                offsetj += nlevelsindep[n - 1];
                for (i = 0; i < nlevels[k]; i++) {
                    if ( skip[k][i] ) {
                        skipoffi++;
                        continue;
                    }
                    in1 = info[k][i];
                    in2 = info[k][i + 1];
                    gf_count_range(groupid[n], index[k] + in1, skip[n], in2 - in1, nlevels[n], ncommon);
                    skipoffj = 0;
                    for (j = 0; j < nlevels[n]; j++) {
                        if ( skip[n][j] ) {
                            skipoffj++;
                        }
                        else if ( ncommon[j] ) {
                            ZZ.insert(i + offseti - skipoffi, j + offsetj - skipoffj) = (ST_double) ncommon[j];
                            ZZ.insert(j + offsetj - skipoffj, i + offseti - skipoffi) = (ST_double) ncommon[j];
                        }
                    }
                }
            }
            offseti += nlevelsindep[k];
        }

        if ( benchmark )
            sf_running_timer(&timer, "\tPlugin Step 2: Cross product");

    /**********************************************************************
     * 3. Detect collinear columns
     **********************************************************************/

        ZZ.makeCompressed();
        solveQR.compute(ZZ);

        Matrix<ST_double, Dynamic, 1> R = solveQR.matrixR().diagonal();
        Matrix<int32_t, Dynamic, 1> P = solveQR.colsPermutation().indices();
        Matrix<ST_double, Dynamic, 1> Rsort = R;

        kindep = 0;
        maxr   = 0;
        tol    = MANYIV_PWMAX(1e-12, numeric_limits<double>::epsilon() * ktot);
        for (i = 0; i < ktot; i++) {
            z = abs(R(i));
            if ( z > maxr ) maxr = z;
            Rsort(P(i)) = z;
        }

        k = 0;
        j = 0;
        offseti = 0;
        skipoffi = 0;
        for (i = 0; i < ktot + kskip; i++) {
            if ( i >= (nlevels[k] + offseti) ) {
                j = 0;
                offseti += nlevels[k++];
            }
            if ( skip[k][i - offseti] ) {
                skipoffi++;
                groupmask[k][i - offseti] = 0;
            }
            else if ( (Rsort(i - skipoffi) / maxr) < tol ) {
                nlevelsindep[k]--;
                groupmask[k][i - offseti] = 0;
            }
            else {
                indepix[kindep++] = i - skipoffi;
                groupmask[k][i - offseti] = ++j;
            }
        }

        if ( benchmark )
            sf_running_timer(&timer, "\tPlugin Step 3: Collinearity");

    /**********************************************************************
     * 4. Compute inverse after removing collinear cols
     **********************************************************************/

        SparseMatrix<ST_double> ZZindep(kindep, kindep);
        ZZindep.setZero();
        ZZindep.reserve(nnz);

        // for (i = 0; i < kindep; i++) {
        //     for (j = 0; j < kindep; j++) {
        //         if ( (z = ZZ.coeff(indepix[i], indepix[j])) ) {
        //             ZZindep.insert(i, j) = z;
        //         }
        //     }
        // }

        for (j = 0; j < kindep; j++) {
            i = 0;
            for (SparseMatrix<double>::InnerIterator it(ZZ, indepix[j]); it; ++it) {
                while ( (indepix[i] < it.row()) && (i < kindep) ) {
                    i++;
                }
                if ( indepix[i] == it.row() && (i < kindep) ) {
                    ZZindep.insert(i, j) = it.value();
                }
            }
        }

        ZZindep.makeCompressed();
        solverLL.compute(ZZindep);

        ZZ.resize(0, 0);
        ZZindep.resize(0, 0);

        if ( solverLL.info() != Success) rc = 1701;
        if ( benchmark )
            sf_running_timer(&timer, "\tPlugin Step 4: Inversion");

    /**********************************************************************
     * 5. Projection diagonal is z_i' inv(Z' Z) z_i
     **********************************************************************/

        if ( rc == 0 ) {
            SparseMatrix<ST_double> eye(kindep, kindep);
            eye.setIdentity();

            SparseMatrix<double> invZZ = solverLL.solve(eye);
            memset(d_projection, '\0', nobs * sizeof(*d_projection));
            for (n = 0; n < nobs; n++) {
                offseti = 0;
                for (i = 0; i < nabsorb; i++) {
                    offsetj = 0;
                    for (j = 0; j < nabsorb; j++) {
                        in1 = groupmask[i][groupid[i][n]];
                        in2 = groupmask[j][groupid[j][n]];
                        if ( in1 && in2 ) {
                            d_projection[n] += invZZ.coeff(--in1 + offseti, --in2 + offsetj);
                        }
                        offsetj += nlevelsindep[j];
                    }
                    offseti += nlevelsindep[i];
                }
            }
            invZZ.resize(0, 0);

             fstream fh(fname, ios::out | ios::binary);
             if ( fh.is_open() ) {
                 fh.write((char*) d_projection, nobs * sizeof(*d_projection));
                 for (k = 0; k < nabsorb; k++) {
                     fh.write((char*) groupmask[k], nlevels[k] * sizeof(*groupmask[k]));
                 }
                 fh.close();
             }
             else {
                 rc = 1703;
             }
        }

        if ( benchmark )
            sf_running_timer(&timer, "\tPlugin Step 5: Diagonal projection");

    /**********************************************************************
     * 6. Free stuff
     **********************************************************************/

        for (k = 0; k < nabsorb; k++) {
            free(index[k]);
            free(groupid[k]);
            free(info[k]);
            free(skip[k]);
            free(groupmask[k]);
        }

        free(ncommon);
        free(indepix);
        free(d_projection);
    }

    return (rc);
}
