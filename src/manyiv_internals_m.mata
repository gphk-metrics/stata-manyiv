cap mata mata drop ManyIVStats()
cap mata mata drop ManyIVreg_IM()

cap mata mata drop sf_helper_epsilon()
cap mata mata drop sf_helper_sig()
cap mata mata drop sf_helper_annihilator()
cap mata mata drop sf_helper_solve()
cap mata mata drop sf_helper_tsolve()

mata
struct ManyIVStats {
    real scalar F
    real matrix Omega
    real matrix Xi
    real vector Sargan
    real vector CD
}

class ManyIVreg_IM
{
    real rowvector beta
    real matrix se
    struct ManyIVStats scalar stats
    real scalar n, K, L, F

    real scalar clustered
    real scalar estimatese
    real scalar estimatestats

    void results()
    void print()
    void fit()

    string colvector betaLabels
    string colvector seLabels
}

void function ManyIVreg_IM::fit(
    real colvector y,
    real colvector T,
    real matrix Z,
    real matrix W,
    real scalar estimatese,
    real scalar estimatestats, |real colvector cluster)
{
    real scalar i, j
    real vector overid, pvalue
    real matrix Xi

    real matrix Yp, Zp, RFS, YY, YPY, YMY, ZW, Cinfo
    real matrix epsilon, epsilon_i, hatP, hatP_i
    real vector k, ei, sec, DZW, DW, iIDZW, iIDW
    real vector hatTjive, hatTujive, hatPjive, hatPujive

    n = rows(Z)
    K = cols(Z)
    L = cols(W)
    betaLabels = ("OLS", "TSLS", "LIML", "MBTSLS", "JIVE", "UJIVE", "RTSLS")'
    seLabels   = ("Homoskedastic", "Heteroscedastic", "Cluster", "ManyIV", "ManyIV")'

    // 2. Point estimates
    // Note: These apply Frisch–Waugh–Lovell; project the covariate
    // of interest and the instrument onto W's null space then use
    // univariate formulas.

    Yp  = sf_helper_annihilator(W, (y, T))   // [y_⊥ T_⊥] = M_W [y T]
    Zp  = sf_helper_annihilator(W, Z)        // Z_⊥       = M_W Z
    YY  = (Yp' * Yp)                         // [y_⊥ T_⊥]' [y_⊥ T_⊥] = [y T]' M_W [y T]
    RFS = sf_helper_solve(Zp, Yp)            // [solve(Z_⊥, y_⊥) solve(Z_⊥, T_⊥)] = Reduced form and First stage
    YPY = (Yp' * Zp) * RFS                   // [y_⊥ T_⊥]' H_{Z_⊥} [y_⊥ T_⊥]
    YMY = YY - YPY                           // [y_⊥ T_⊥]' M_{Z_⊥} [y_⊥ T_⊥]

    // 2.1 k-class: OLS, TSLS, LIML, MBTLS
    // Note: These are all coded as
    //
    //     (T_⊥' (I - k W_{Z_⊥}) y_⊥) / (T_⊥' (I - k W_{Z_⊥}) T_⊥)
    //
    // So different values of k give different estimands.
    //
    // - k = 0 -> YY[1, 2]  / YY[2, 2]  = (y_⊥' T_⊥) / (T_⊥' T_⊥)
    // - k = 1 -> YPY[1, 2] / YPY[2, 2] = (y_⊥' H_{Z_⊥} T_⊥) / (T_⊥' H_{Z_⊥} T_⊥)
    // - The other two give liml and mbtsls

    eigensystem(invsym(YMY) * YY, ., ei=.)
    k     = (0, 1, min(Re(ei)), (1 - L/n) / (1 - (K - 1) / n - L/n))
    beta = (YY[1, 2] :- k :* YMY[1, 2]) :/ (YY[2, 2] :- k :* YMY[2, 2])

    // 2.2 JIVE, UJIVE
    ZW     = (W, Z)
    DZW    = rowsum((ZW * invsym(ZW' * ZW)) :* ZW) // D_{Z W} = diag(H_{Z W}) as a vector
    iIDZW  = 1 :/ (1 :- DZW)                       // (I - D_{Z W})^{-1} as a vector
    DW     = rowsum((W * invsym(W'*W)) :* W)       // D_W = diag(H_W) as a vector
    iIDW   = 1 :/ (1 :- DW)                        // (I - D_W)^{-1} as a vector

    hatTujive = T :- iIDZW :* sf_helper_annihilator(ZW, T) //     (I - (I - D_{Z W})^{-1} M_{Z W}) T
    hatPjive  = sf_helper_annihilator(W, hatTujive)        // M_W (I - (I - D_{Z W})^{-1} M_{Z W}) T
    hatTjive  = T :- iIDW :* sf_helper_annihilator(W, T)   // (I - (I - D_W)^{-1} M_W) T
    hatPujive = hatTujive - hatTjive                       // (I - D_W)^{-1} M_W T - (I - D_{Z W})^{-1} M_{Z W} T
                                                           // = (I - D_W)^{-1} (I - H_W) T - (I - D_{Z W})^{-1} (I - H_{Z W}) T
                                                           // = ((I - D_{Z W})^{-1} (H_{Z W} - I) - (I - D_W)^{-1} (H_W - I) T) T
                                                           // = (
                                                           //     (I - D_{Z W})^{-1} (H_{Z W} - D_{Z W} - (I - D_{Z W})) -
                                                           //     (I - D_W)^{-1} (H_W - D_W - (I - D_W))
                                                           // ) T
                                                           // = ((I - D_{Z W})^{-1} (H_{Z W} - D_{Z W}) - (I - D_W)^{-1} (H_W - D_W)) T

    beta =  (
        beta,
        (hatPjive'  * y) / (hatPjive'  * T),
        (hatPujive' * y) / (hatPujive' * T),
        YPY[1, 1] / YPY[1, 2]
    )

    // ------------------
    // 5. Standard Errors
    // ------------------

    Sp = YMY/(n-K-L)  // S_{perp}
    se = J(5, 7, .)
    if ( estimatese ) {

        // -----------------
        // 5.1 Homoscedastic
        // -----------------

        se[1, 1::6] = sqrt((
            sf_helper_sig(Yp, beta[1]) / (Yp[.,2]'*Yp[.,2]),
            (
                sf_helper_sig(Yp, beta[2]), sf_helper_sig(Yp, beta[3]), sf_helper_sig(Yp, beta[4])
            ) / YPY[2,2],
            sf_helper_sig(Yp, beta[5]) * (hatPjive'  * hatPjive)  / (hatPjive'  * T)^2,
            sf_helper_sig(Yp, beta[6]) * (hatPujive' * hatPujive) / (hatPujive' * T)^2
       ))

        // -------------------
        // 5.2 Heteroscedastic
        // -------------------

        // ols, tsls, liml, mbtsls, jive, ujive
        // Note for asymptotics k -> 1 for the various estimators.
        hatP    = Yp[.,2], J(1, 3, Zp * sf_helper_solve(Zp, Yp[.,2])), hatPjive, hatPujive
        epsilon = sf_helper_epsilon(Yp, beta[1::6])
        se[2, 1::6] = sqrt(colsum((epsilon :* hatP):^2)) :/ (T' * hatP)

        // -------------------------
        // 5.4 Cluster, if requested
        // -------------------------

        clustered = missing(cluster) == 0
        if ( clustered ) {
            sec   = J(1, 6, 0)
            Cinfo = panelsetup(cluster, 1)
            for(i = 1; i <= rows(Cinfo); i++) {
                hatP_i    = panelsubmatrix(hatP,    i, Cinfo)
                epsilon_i = panelsubmatrix(epsilon, i, Cinfo)
                for (j = 1; j <= 6; j++) {
                    sec[j] = sec[j] + hatP_i[.,j]' * ((epsilon_i[.,j] * epsilon_i[.,j]') * hatP_i[.,j])
                }
            }
            se[3, 1::6] = sqrt(sec) :/ (T' * hatP)
        }

        // --------------------
        // 5.3 Many instruments
        // --------------------

        // Notation
        S    = YPY/n
        eigensystem(sf_helper_solve(Sp, S), ., ei=.)
        mmin = min(Re(ei))

        // Hessian of random-effects
        lamre = max(Re(ei)) - K/n;
        a     = beta[3] \ 1
        b     = 1 \ -beta[3]
        Omre  = (n-K-L) * Sp/(n-L) + n * (S :- lamre * sf_helper_tsolve((a*a'), sf_helper_tsolve(a', Sp) * a)) / (n-L);
        Qs    = (b' * S * b) / (b' * Omre * b)
        c     = lamre * Qs / ((1-L/n) * (K/n+lamre))

        se[4, 3] = sqrt(
            -b'*Omre*b / (n*lamre) * (lamre+K/n) /
            (Qs*Omre[2,2] - S[2,2] + (c/(1-c)) * Qs / (sf_helper_tsolve(a', Omre) * a))
        )

        // mbtsls, using maximum URE likelihood plug-in estimator
        b = 1 \ -beta[4] // b_mbtsls
        Lam11 = max((0, b' * (S - K/n * Sp) * b))

        if ( mmin > K/n ) {
            Lam22 = S[2, 2] - K/n * Sp[2, 2]
            Omure = Sp
        }
        else {
            Lam22 = lamre / (sf_helper_tsolve(a', Omre) * a)
            Omure = Omre
        }

        Gamma = (1, 0) \ (-beta[4], 1)
        Sig   = Gamma' * Omure * Gamma
        h     = ((1-L/n) * (K-1)/n) / (1 - L/n - (K-1)/n)

        Vvalid   = Sig[1,1] / Lam22 + h * (Sig[1,1] * Sig[2,2] + Sig[1,2]^2) / Lam22^2
        Vinvalid = Vvalid + (Lam11 * Omure[2,2] + Lam11 * Lam22 * n/K) / Lam22^2

        se[4::5, 4] = sqrt((Vvalid \ Vinvalid) / n)
    }

    F = YPY[2, 2] / (K * Sp[2,2]) // First-stage F
    if ( estimatestats ) {
        Xi = YPY/n - (K/n) * Sp // % Xi

        overid = J(2, 1, .)
        pvalue = J(2, 1, .)
        if ( K > 1 ) {
            overid[1] = n * mmin / (1 - K/n - L/n + mmin) // n* J_sargan
            pvalue[1] = chi2tail(K - 1, overid[1])        // p-value for Sargan

            overid[2] = n * mmin  // Cragg-Donald
            pvalue[2] = 1 - normal(sqrt((n-K-L)/(n-L)) * invnormal(chi2(K - 1, overid[2])))
        }

        stats.F      = F
        stats.Omega  = Sp
        stats.Xi     = Xi
        stats.Sargan = overid[1], pvalue[1]
        stats.CD     = overid[2], pvalue[2]
    }
}

void function ManyIVreg_IM::results(string scalar bname, string scalar sename)
{

    if ( (bname != "") ) {
        st_matrix(bname, beta)
        st_matrixcolstripe(bname, (J(rows(betaLabels), 1, ""), betaLabels))
    }

    if ( (sename != "") ) {
        st_matrix(sename, se)
        st_matrixcolstripe(sename, (J(rows(betaLabels), 1, ""), betaLabels))
        st_matrixrowstripe(sename, (J(rows(seLabels),   1, ""), seLabels))
    }
}

void function ManyIVreg_IM::print()
{
    real scalar j, maxl
    maxl = max(strlen(betaLabels)) + 1
    if ( estimatese ) {
        printf(sprintf("%%%gs %%9s %%11s\n", maxl), "", "Coef.", clustered? "Cluster": "Hetero")
        printf(sprintf("%%%gs %%9s %%11s\n", maxl), maxl * " ", 5 * "-", (clustered? 7: 6) * "-")
        for(j = 1; j <= length(beta); j++) {
            printf(sprintf("%%%gs %%9.4f (%%9.4f)\n", maxl), betaLabels[j], beta[j], se[2 :+ clustered, j])
        }
    }
    else {
        printf(sprintf("%%%gs %%9s\n", maxl), "", "Coef.")
        printf(sprintf("%%%gs %%9s\n", maxl), maxl * " ", 5 * "-")
        for(j = 1; j <= length(beta); j++) {
            printf(sprintf("%%%gs %%9.4f\n", maxl), betaLabels[j], beta[j])
        }
    }
    printf("\n%g observations, ", n)
    printf("%g instrument%s, ", K, (K > 1)? "s": "")
    printf("%g covariate%s, ", L, (L > 1)? "s": "")
    printf("first-stage F = %5.1fc\n", F)
}

real matrix function sf_helper_epsilon(real matrix Yp, real rowvector beta)
{
    return(Yp[., 1] :- Yp[., 2] * beta)
}

real scalar function sf_helper_sig(real matrix Yp, real scalar beta)
{
    real colvector e
    e = sf_helper_epsilon(Yp, beta)
    return(e' * e / length(e))
}

real matrix function sf_helper_annihilator(real matrix X, real matrix Y)
{
    return(Y - X * sf_helper_solve(X, Y))
}

real matrix function sf_helper_solve(real matrix X, real matrix Y)
{
    return(invsym(cross(X, X)) * cross(X, Y))
}

real matrix function sf_helper_tsolve(real matrix X, real matrix Y)
{
    return((invsym(cross(Y', Y')) * cross(Y', X'))')
}
end
