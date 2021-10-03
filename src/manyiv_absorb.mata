cap mata mata drop ManyIVreg_Absorb()
cap mata mata drop New_ManyIVreg_Absorb()
mata
class ManyIVreg_Absorb
{
    transmorphic absorbinfo
    string vector absorbvars
    string scalar touse
    real scalar nobs
    real scalar nabsorb
    real scalar hdfetol
    real scalar maxiter
    real scalar encoded
    real vector nlevels
    real vector total_nlevels

    void new()
    void init()
    real scalar makepanel()
    real colvector index()
    real colvector select()
    real matrix info()
    real colvector nj()
    real colvector groupid()
    void encode()
    real scalar _demean()
    real matrix hdfe()
    void _hdfe()
}

class ManyIVreg_Absorb scalar function New_ManyIVreg_Absorb(string vector _absorbvars, string scalar _touse)
{
    class ManyIVreg_Absorb scalar Absorb 
    Absorb = ManyIVreg_Absorb()
    Absorb.init(_absorbvars, _touse)
    return(Absorb)
}

void function ManyIVreg_Absorb::new()
{
    this.hdfetol = 1e-8
    this.maxiter = 1000
    this.encoded = 0
}

void function ManyIVreg_Absorb::init(string vector _absorbvars, string scalar _touse)
{
    absorbvars = _absorbvars
    touse = _touse
    nabsorb = length(absorbvars)
    absorbinfo = asarray_create()
    total_nlevels = 0
    nlevels = J(1, nabsorb, .)
    for(j = 1; j <= nabsorb; j++) {
        nlevels[j] = makepanel(absorbvars[j])
        total_nlevels = total_nlevels + nlevels[j]
    }
}

real scalar function ManyIVreg_Absorb::makepanel(string scalar var)
{

    string colvector svar
    real colvector varindex, rvar
    real matrix varinfo

    if ( strpos(st_vartype(var), "str") ) {
        svar     = st_sdata(., var, touse)
        varindex = order(svar, 1)
        varinfo  = panelsetup(svar[varindex], 1)
    }
    else {
        rvar     = st_data(., var, touse)
        varindex = order(rvar, 1)
        varinfo  = panelsetup(rvar[varindex], 1)
    }
    asarray(absorbinfo, var + ".index", varindex)
    asarray(absorbinfo, var + ".info",  varinfo)
    asarray(absorbinfo, var + ".nj",    varinfo[., 2] :- varinfo[., 1] :+ 1)

    nobs = length(varindex)
    return(rows(varinfo))
}

real matrix function ManyIVreg_Absorb::info(real scalar j)
{
    return(asarray(absorbinfo, absorbvars[j] + ".info"))
}

real colvector function ManyIVreg_Absorb::index(real scalar j)
{
    return(asarray(absorbinfo, absorbvars[j] + ".index"))
}

real colvector function ManyIVreg_Absorb::nj(real scalar j)
{
    return(asarray(absorbinfo, absorbvars[j] + ".nj"))
}

real colvector function ManyIVreg_Absorb::groupid(real scalar j)
{
    return(asarray(absorbinfo, absorbvars[j] + ".groupid"))
}

real colvector function ManyIVreg_Absorb::select(real scalar j, real scalar i)
{
    return(panelsubmatrix(index(j), i, info(j)))
}

real matrix function ManyIVreg_Absorb::hdfe(real matrix X, | real scalar base)
{
    _hdfe(X, base)
    return(X)
}

void function ManyIVreg_Absorb::_hdfe(real matrix X, | real scalar base)
{
    real scalar i, dev
    if ( nabsorb == 1 ) {
        (void) _demean(X, 1, base)
    }
    else if ( nabsorb > 1 ) {
        i       = 0
        dev     = 1
        while ( i++ < maxiter ) {
            for (j = 1; j <= nabsorb; j++) {
                dev = _demean(X, j, base)
                if ( dev < hdfetol ) break
            }
            if ( dev < hdfetol ) break
        }
        if ( i > maxiter ) {
            errprintf("maximum number of hdfe iterations exceeded (%g)\n", maxiter)
            error(1234)
        }
        else {
            printf("hdfe convergence after %g projections (error = %5.3g)", i * nabsorb + j, dev)
        }
    }
}

real scalar function ManyIVreg_Absorb::_demean(real matrix X, real scalar j, | real scalar base)
{
    real colvector avg
    real colvector sel
    real matrix X_i
    real scalar i, dev

    dev = 0
    if ( args() < 3 ) base = 0
    for(i = 1; i <= nlevels[j]; i++) {
        if ( i == base ) continue
        sel = panelsubmatrix(index(j), i, info(j))
        X_i = X[sel, .]
        avg = (colsum(X_i) :/ length(sel))
        X[sel, .] = X_i :- avg
        dev = max((abs(avg), dev))
    }

    return(dev)
}

void function ManyIVreg_Absorb::encode()
{
    real colvector sel, groupid
    real scalar i, j

    if ( encoded == 0 ) {
        for(j = 1; j <= nabsorb; j++) {
            groupid = J(nobs, 1, .)
            for(i = 1; i <= nlevels[j]; i++) {
                sel = panelsubmatrix(index(j), i, info(j))
                groupid[sel] = J(length(sel), 1, i)
            }
            asarray(absorbinfo, absorbvars[j] + ".groupid", groupid)
        }
    }
    encoded = 1
}
end
