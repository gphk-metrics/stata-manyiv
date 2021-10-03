*! version 0.2.0 30Sep2021
*! Instrumental variables regression (OLS, TSLS, LIML, MBTSLS, JIVE, UJIVE, RTSLS)
*! Adapted for Stata from ivreg.m by Michal Kolesár <kolesarmi@googlemail dotcom>

* TODO: allow absorbiv with no instruments
* TODO: absorb and absorbiv don't account for collinearity atm
* TODO: number of instruments check does not account for absorbed instruments
* TODO: detect repeated variables between absorb and varist (z or w varlist)
capture program drop manyiv
program manyiv, eclass sortpreserve
    syntax anything(equalok) /// dependent_var covariates
           [if] [in] ,       /// subset
    [                        ///
        absorb(str)          /// Absorb controls
        absorbiv(str)        /// Absorb instruments
                             ///
        internals(str)       ///
        SAVEresults(str)     /// name of mata object to store results
                             ///
        /// noCONStant       /// omit constant; a bit of a pain to account for
        noPRINTtable         /// do not print table
        nose                 /// omit standard errors
        nostats              /// omit additional statistics
        cluster(varname)     /// clustered SE
                             ///
    ]

    local estimatese    = ("`estimatese'"    == "")
    local estimatestats = ("`estimatestats'" == "")

    ***********************************************************************
    *                           Parse IV Syntax                           *
    ***********************************************************************

    if !regexm(`"`anything'"', ".+\((.+=.+)\)") {
        disp as err "Unable to parse IV syntax: depvar (endogenous = instrument) [exogenous]"
        exit 198
    }

    local iveq   = regexr(regexs(1), "\(|\)", "")
    local ivexog = trim(regexr("`anything'", "\(.+=.+\)", ""))

    gettoken ivdepvar ivexog: ivexog
    gettoken ivendog ivinst: iveq, p(=)
    gettoken _ ivinst: ivinst

    helper_syntax_error `ivdepvar', msg(dependent variable) plain
    helper_syntax_error `ivinst',   msg(instrument variables)
    helper_syntax_error `ivendog',  msg(endogenous variables)

    ***********************************************************************
    *                           Parse varlists                            *
    ***********************************************************************

    local varlist `ivdepvar' `ivendog' `ivinst' `ivexog'
    marksample touse

    tempname Y X Z W C A IV
    fvexpand `ivexog'
    local ivexog `r(varlist)'

    fvexpand `ivinst'
    local ivinst `r(varlist)'
    local ivinst: list ivinst - ivexog

    mata `Y' = st_data(., "`ivdepvar'", "`touse'")
    mata `X' = st_data(., "`ivendog'",  "`touse'")

    helper_strip_omitted `ivinst' if `touse',
    mata `Z' = select(st_data(., "`ivinst'", "`touse'"), !st_matrix("r(omit)"))
    local ivinst `r(varlist)'

    if ( `:list sizeof ivexog' ) {
        helper_strip_omitted `ivexog' if `touse', `constant'
        mata `W' = select(st_data(., "`ivexog'", "`touse'"), !st_matrix("r(omit)")), J(rows(`Y'), 1, 1)
        local ivexog `r(varlist)'
    }
    else {
        mata `W' = J(rows(`Y'), 1, 1)
    }

    tempname kendog kinst
    mata: st_numscalar("`kendog'", cols(st_data(1, "`ivendog'")))
    mata: st_numscalar("`kinst'",  cols(st_data(1, "`ivinst'")))

    if ( scalar(`kinst') < scalar(`kendog') ) {
        disp as error "Need at least as many instruments as endogenous variables"
        exit 198
    }

    if ( scalar(`kendog') != 1 ) {
        disp as error "Only one endogenous variabe allowed (currently hard-coded)"
        exit 198
    }

    * TODO: Decide whether to make hard errors or leave as warnings
    * -------------------------------------------------------------
    * Finally, flag variables that are both
    *
    * - dependent variable _and_ instrumented
    * - dependent variable _and_ instrument
    * - dependent variable _and_ exogenous
    * - instrumented _and_ instrument
    * - instrumented _and_ exogenous
    * - instrument _and_ exogenous

    local problems: list ivdepvar & ivendog
    if ( `"`problems'"' != `""' ) {
        disp as error "{bf:warning:} `problems' included as both regressand and endogenous variable"
    }

    local problems: list ivdepvar & ivexog
    if ( `"`problems'"' != `""' ) {
        disp as error "{bf:warning:} `problems' included as both regressand and exogenous variable"
    }

    local problems: list ivdepvar & ivinst
    if ( `"`problems'"' != `""' ) {
        disp as error "{bf:warning:} `problems' included as both regressand and instrument"
    }

    local problems: list ivendog  & ivexog
    if ( `"`problems'"' != `""' ) {
        disp as error "{bf:warning:} included as both an endogenous and exogenous variable: `problems'"
    }

    local problems: list ivendog  & ivinst
    if ( `"`problems'"' != `""' ) {
        disp as error "{bf:warning:} included as both an endogenous variable and an instrument: `problems'"
    }

    local problems: list ivexog   & ivinst
    if ( `"`problems'"' != `""' ) {
        disp as error "{bf:warning:} included as both an exogenous variable and an instrument: `problems'"
    }

    ***********************************************************************
    *                             Estimation                              *
    ***********************************************************************

    local ManyIVreg:  copy local saveresults
    local ManyIVData: copy local savedata

    tempname results beta se
    if "`ManyIVData'" == "" tempname ManyIVData
    if "`ManyIVreg'"  == "" tempname ManyIVreg

    mata `A'  = New_ManyIVreg_Absorb(tokens(st_local("absorb")),   "`touse'")
    mata `IV' = New_ManyIVreg_Absorb(tokens(st_local("absorbiv")), "`touse'")

    if ("`cluster'" != "") {
        tempvar clusterid
        sort `cluster' `touse'
        by `cluster' `touse': gen long `clusterid' = (_n == 1) & `touse'
        replace `clusterid' = sum(`clusterid')
        replace `clusterid' = . if !`touse'
        mata `C' = st_data(., "`clusterid'",  "`touse'")
    }
    else mata `C' = .

    mata `ManyIVreg' = ManyIVreg_IM()
    mata `ManyIVreg'.fit(`Y', `X', `Z', `W', `estimatese', `estimatestats', `C', `A', `IV')
    mata `ManyIVreg'.results("`beta'", "`se'")

    if ( "`printtable'" != "noprinttable" ) {
        mata `ManyIVreg'.print()
    }

    ereturn post `beta', esample(`touse')
    if ( `estimatese' ) {
        ereturn matrix se = `se'
    }

    tempname F Omega Xi Sargan CD
    mata st_numscalar("`F'", `ManyIVreg'.stats.F)
    ereturn scalar F = `F'
    if ( `estimatestats' ) {
        mata st_matrix("`Omega'",  `ManyIVreg'.stats.Omega)
        mata st_matrix("`Xi'",     `ManyIVreg'.stats.Xi)
        mata st_matrix("`Sargan'", `ManyIVreg'.stats.Sargan)
        mata st_matrix("`CD'",     `ManyIVreg'.stats.CD)

        ereturn matrix Omega  = `Omega'
        ereturn matrix Xi     = `Xi'
        ereturn matrix Sargan = `Sargan'
        ereturn matrix CD     = `CD'
    }
end

capture program drop helper_syntax_error
program helper_syntax_error
    syntax anything(equalok), [msg(str) plain]
    if ( `:list sizeof anything' == 0 ) {
        disp as err "Unable to parse IV syntax: no `msg' variable detected"
        exit 198
    }

    loca 0: copy local anything
    if ( "`plain'" == "pain" ) {
        cap noi syntax varlist(numeric)
    }
    else {
        cap noi syntax varlist(numeric fv ts)
    }

    if ( _rc ) {
        disp as err "Unable to parse IV syntax (`msg')"
        exit _rc
    }
end

capture program drop helper_strip_omitted
program helper_strip_omitted, rclass
    syntax anything(equalok) [if], [*]
    _rmcoll `anything' `if', expand `options'
    local expanded `r(varlist)'

    tempname b omit
    matrix `b' = J(1, `:list sizeof expanded', .)
    matrix colnames `b' = `expanded'
    _ms_omit_info `b'
    matrix `omit' = r(omit)
    mata st_local("varlist", invtokens(select(tokens("`expanded'"), !st_matrix("r(omit)"))))

    return local expanded: copy local expanded
    return local varlist:  copy local varlist
    return matrix omit =  `omit'
end

// cap findfile manyiv_absorb.mata
// cap do `"`r(fn)'"'
// 
// cap findfile manyiv_internals_m.mata
// cap do `"`r(fn)'"'
