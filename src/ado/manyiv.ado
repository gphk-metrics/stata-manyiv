*! version 0.6.4 29Jul2022
*! Instrumental variables regression (OLS, TSLS, LIML, MBTSLS, JIVE, UJIVE, RTSLS)
*! Based on ivreg.m by Michal Kolesár <kolesarmi@googlemail dotcom>
*! Adapted for Stata by Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com>

* Idea: xx Add a token to start and end of C++ file as a check;
*       if token at start and end doesn't match then spit out error
*       "unable to compute d_projection using plugin; temporary
*       file corrupted". Internal element can be called d_token or
*       similar. String scalar.

* TODO: xx after the above, add unit tests, including error checking
*       and expected error exits and CLT sims for SEs.

* TODO: xx It's not obvious how to detect collinearity (numerically);
*       ask if this algorithm is OK:
*       - QR decomposition
*       - Diagonal entries of R matrix
*       - Scale by largest element (if > 1)
*       - Drop if below tolerance (max of 1e-12 and machine epsilon * nobs)

capture program drop manyiv
program manyiv, eclass
    version 14.1

    ***********************************************************************
    *                            Plugin Check                             *
    ***********************************************************************

    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    local _plugin_file manyiv_`c_os_'.plugin
    cap plugin call manyiv_plugin
    if ( _rc ) {
        local _plugin_notloaded _plugin_notloaded
        cap findfile "`_plugin_file'"
        if ( _rc ) {
            local _plugin_filenotfound _plugin_filenotfound
        }
    }

    if ( "`0'" == "_plugin_check" ) {
        cap noi plugin call manyiv_plugin, `"_plugin_check"'
        exit _rc
    }
    else {
        cap plugin call manyiv_plugin, `"_plugin_check"'
        if ( _rc ) {
            local _plugin_skip_override _plugin_skip_override
        }
    }

    ***********************************************************************
    *                           Actual Program                            *
    ***********************************************************************

    if _N == 0 error 2000

    FreeTimer
    local t99: copy local FreeTimer
    manyiv_timer on `t99'

    syntax anything(equalok) /// dependent_var covariates
           [if] [in] ,       /// subset
    [                        ///
        absorb(str)          /// Absorb controls
        absorbiv(str)        /// Absorb instruments
        skipsingletons       /// skip singleton groups
        keepsingletons       /// do not drop singleton groups
        forcejive            /// make sure jive/ujive will run
        cluster(varname)     /// clustered SE
                             ///
        internals(str)       ///
        SAVEresults(str)     /// name of mata object to store results
        method(str)          /// SQUAREM (default), CG (conjugate gradient), FPI (fixed point iteration)
                             ///
        noConstant           /// omit constant
        nosmall              /// omit small-sample adjustments
        noPRINTtable         /// do not print table
        nose                 /// omit standard errors
        nostats              /// omit additional statistics
                             ///
        _plugin_skip         ///
        _plugin_bench        ///
    ]

    if ( "`forcejive'" != "" ) {
        foreach opt in skipsingletons keepsingletons {
            if ( "``opt''" != "" ) {
                disp as txt "Option -`opt'- ignored with -forcejive-"
                local `opt'
            }
        }
    }

    local benchmark     = ("`_plugin_bench'" != "")
    local estimatese    = ("`se'"            == "")
    local estimatestats = ("`stats'"         == "")
    local small         = ("`small'"         == "")
    local cons          = ("`constant'"      == "")

    if ( "`method'"   == "" ) local method squarem
    if ( "`absorb'"   != "" ) unab absorb:   `absorb'
    if ( "`absorbiv'" != "" ) unab absorbiv: `absorbiv'

    if ( lower(`"`method'"') == "fpi" ) {
        local method_code 1
    }
    else if ( lower(`"`method'"') == "squarem" ) {
        local method_code 2
    }
    else if ( inlist(lower(`"`method'"'), "conjugate gradient", "conjugate_gradient", "cg") ) {
        local method_code 3
    }
    else {
        disp as err "method() must be one of: cg, squarem, fpi"
        exit 198
    }

    local problems: list absorb & absorbiv
    if ( `"`problems'"' != `""' ) {
        disp as error "{bf:warning:} included as both an exogenous variable and an instrument: `problems'"
    }
    local absorbiv: list absorbiv - absorb

    if ( "`_plugin_notloaded'`_plugin_skip_override'" != "" ) {
        local _plugin_skip _plugin_skip
        if ( `:list sizeof absorb' + `:list sizeof absorbiv' > 2 ) {
            disp as err "{bf:warning:} unable to load helper plugin; jive/ujive skipped with 3+ absorb levels"
            if ( "`_plugin_filenotfound'" != "" ) {
                disp as err "(plugin file '`_plugin_file'' not found; please check your installation)"
            }
        }
    }

    ***********************************************************************
    *                           Parse varlists                            *
    ***********************************************************************

    * 0. Parse IV Syntax
    * ------------------

    ParseIVSyntax `anything', absorbiv(`absorbiv')

    * 1. Touse from varlist
    * ---------------------

    tempname Y X Z W C A IV singlecons
    local varlist `ivdepvar' `ivendog' `ivinst' `ivexog'
    marksample touse

    * 2. Drop singletons
    * ------------------

    DropSingletons `touse' `A' `IV' `method_code', ///
        `keepsingletons' `skipsingletons' absorb(`absorb') absorbiv(`absorbiv')
    scalar `singlecons' = r(singlecons)

    manyiv_timer info `t99' `"Steps 1 & 2. Parse variables; drop singletons"', prints(`benchmark')

    * 3. Collinear fixed effects and diagonal of projection matrix
    * ------------------------------------------------------------

    DiagonalProjection `A' `IV', `_plugin_skip' `_plugin_bench'

    manyiv_timer info `t99' `"Step 3. Compute diagonal projection for FEs"', prints(`benchmark')

    * 4. Variables to use in regressions
    * ----------------------------------

    ReadVarlist `touse' `C' `Y' `X' `W' `Z' `method_code' `singlecons' `constant', ///
        ivcall(`anything') cluster(`cluster') absorb(`absorb') absorbiv(`absorbiv')
    local cons = `r(cons)'

    manyiv_timer info `t99' `"Step 4. Variables to use in regression"', prints(`benchmark')

    ***********************************************************************
    *                             Estimation                              *
    ***********************************************************************

    local ManyIVreg:  copy local saveresults
    local ManyIVData: copy local savedata

    tempname results beta se rc
    if "`ManyIVData'" == "" tempname ManyIVData
    if "`ManyIVreg'"  == "" tempname ManyIVreg

    mata `ManyIVreg' = ManyIVreg_IM()
    mata `ManyIVreg'.loadvars(`Y', `X', `Z', `W', `cons', `A', `IV')

    * Very inefficient recursive loop; only run if requested
    if ( "`forcejive'" != "" ) {
        local checkjive 0
        mata `ManyIVreg'.checkjive("`touse'", `A', `IV')
        qui count if `touse'
        if ( `r(N)' == 0 ) error 2000
        if ( `checkjive' ) {
            disp
            disp as err "{bf:warning:} option -forcejive- will drop observations and/or covariates"
            disp as err "until jive/ujive can be computed. However, it is preferable for the user"
            disp as err "to investigate why jive/ujive fails and manually drop the responsible"
            disp as err "variables or observations. The root issue is a covariate or fixed effect"
            disp as err "group that identifies a single observation."
        }
        qui while ( `checkjive' ) {
            DropSingletons `touse' `A' `IV' `method_code', ///
                `keepsingletons' `skipsingletons' absorb(`absorb') absorbiv(`absorbiv')
            qui count if `touse'
            if ( `r(N)' == 0 ) error 2000
            scalar singlecons = r(singlecons)
            DiagonalProjection `A' `IV', `_plugin_skip' `_plugin_bench'
            ReadVarlist `touse' `C' `Y' `X' `W' `Z' `method_code' `singlecons' `constant', ///
                ivcall(`anything') cluster(`cluster') absorb(`absorb') absorbiv(`absorbiv')
            local cons = `r(cons)'
            mata `ManyIVreg'.loadvars(`Y', `X', `Z', `W', `cons', `A', `IV')
            mata `ManyIVreg'.checkjive("`touse'", `A', `IV')
        }
        manyiv_timer info `t99' `"Step X. Force JIVE/UJIVE"', prints(`benchmark')
    }

    * Finally, fit the model
    mata `ManyIVreg'.fit(`Y', `X', `Z', `W', `estimatese', `estimatestats', `small', `cons', `C', `A', `IV')
    mata `ManyIVreg'.dropvars()

    mata st_numscalar("`rc'", `ManyIVreg'.rc)
    if ( `=scalar(`rc')' ) exit `=scalar(`rc')'
    mata `ManyIVreg'.results("`beta'", "`se'")

    manyiv_timer info `t99' `"Step 5. ManyIV fit"', prints(`benchmark') off

    if ( "`printtable'" != "noprinttable" ) {
        mata `ManyIVreg'.print()
    }

    tempname F Omega Xi Sargan CD rf fs jive kinst
    ereturn post `beta', esample(`touse')
    if ( `estimatese' ) {
        ereturn matrix se    = `se'
        ereturn scalar small = `small'
    }
    mata st_numscalar("`jive'", `ManyIVreg'.jive)
    ereturn scalar jive = `jive'

    ereturn local depvar:       copy local ivdepvar
    ereturn local instrumented: copy local ivendog
    ereturn local instruments:  copy local ivinst
    ereturn local exogenous:    copy local ivexog

    if ( `:list sizeof ivinst' ) {
        mata st_numscalar("`kinst'", rows(`ManyIVreg'.RFS))
        if ( `=scalar(`kinst')' ) {
            mata st_matrix("`rf'", `ManyIVreg'.RFS[.,1])
            mata st_matrix("`fs'", `ManyIVreg'.RFS[.,2])
            ereturn matrix rf = `rf'
            ereturn matrix fs = `fs'
        }
        else {
            disp as txt "{bf:warning:} unable to recover reduced form or first stage"
            if ( `:list sizeof absorb' + `:list sizeof absorbiv' ) {
                disp as txt "(instruments likely collinear with absorb levels)"
            }
        }
    }

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

***********************************************************************
*                                                                     *
*                          Modularized Steps                          *
*                                                                     *
***********************************************************************

capture program drop ParseIVSyntax
program ParseIVSyntax
    syntax anything(equalok), [absorbiv(str)]

    if !regexm(`"`anything' "', ".+\((.+=.+)\)[^\.]") {
        disp as err "Unable to parse IV syntax: depvar (endogenous = instrument) [exogenous]"
        exit 198
    }

    local iveq   = regexs(1)
    local ivexog = trim(regexr("`anything' ", "\(.+=.+\)[^\.]", " "))

    gettoken ivdepvar ivexog: ivexog
    gettoken ivendog ivinst: iveq, p(=)
    gettoken _ ivinst: ivinst, p(=)
    local ivinst `ivinst'

    helper_syntax_error `ivdepvar', msg(dependent variable) plain
    helper_syntax_error `ivendog',  msg(endogenous variables)
    if ( inlist("`ivinst'", ".", "") ) {
        if ( "`absorbiv'" == "" ) {
            disp as err "no instruments specified"
            exit 198
        }
        local ivinst ""
    }
    else helper_syntax_error `ivinst', msg(instrument variables)

    c_local ivdepvar: copy local ivdepvar
    c_local ivendog:  copy local ivendog
    c_local ivinst:   copy local ivinst
    c_local ivexog:   copy local ivexog
end

capture program drop DropSingletons
program DropSingletons, rclass
    syntax anything, [skipsingletons keepsingletons absorb(str) absorbiv(str)]
    tokenize `anything'
    local touse:       copy local 1
    local A:           copy local 2
    local IV:          copy local 3
    local method_code: copy local 4

    tempname singlecons nsingletons addnsingletons singleix onesingleton tag

    * Note: Need this here to drop singletons
    mata `A'  = ManyIVreg_Absorb_New(tokens(st_local("absorb")),   "`touse'", `method_code')
    mata `IV' = ManyIVreg_Absorb_New(tokens(st_local("absorbiv")), "`touse'", `method_code')
    mata `IV'.append(`A')

    scalar `nsingletons' = 0
    scalar `singlecons'  = 0
    if ( "`keepsingletons'`skipsingletons'" == "" ) {
        scalar `addnsingletons' = 1
        while ( `=scalar(`addnsingletons')' ) {
            * method 1 drops singleton groups across absorb variables
            * and re-indexes the internal index and encoded group IDs.
            if ( `:list sizeof absorbiv' ) {
                mata `IV'.dropsingletons(`singleix'=J(0, 1, .), 1)
                mata `A'.dropfromindex(`singleix')
            }
            else {
                mata `A'.dropsingletons(`singleix'=J(0, 1, .), 1)
            }

            mata st_numscalar("`addnsingletons'", rows(`singleix'))
            scalar `nsingletons' = `nsingletons' + `addnsingletons'
            if ( `=scalar(`addnsingletons')' ) {
                mata `tag' = st_data(., "`touse'", "`touse'")
                mata `tag'[`singleix'] = J(rows(`singleix'), 1, 0)
                mata st_store(., "`touse'", "`touse'", `tag')
            }
        }
        scalar `addnsingletons' = 0
    }
    else if ( "`skipsingletons'" != "" ) {
        * method 2 flags singleton groups to be skipped
        mata `IV'.dropsingletons(., 2)
        mata `A'.dropsingletons(., 2)
        mata st_numscalar("`onesingleton'", `IV'.onesingleton() | `A'.onesingleton())
        mata st_numscalar("`singlecons'", `IV'.allhaveskip() & `A'.allhaveskip())
        if ( `=scalar(`onesingleton')' ) {
            disp as err "{bf:warning} jive/ujive unstable with -skipsingletons- and an absorb group has only one singleton"
        }
    }
    else if ( "`keepsingletons'" != "" ) {
        * method 3 only counts number of singletons
        mata `IV'.dropsingletons(., 3)
        mata `A'.dropsingletons(., 3)
    }

    if ( `=scalar(`nsingletons')' ) {
        if ( `=scalar(`nsingletons')' > 1 ) local s s
        local strnsingletons = trim("`:disp %21.0fc scalar(`nsingletons')'")
        disp as txt "{bf:note: will drop `strnsingletons' singleton group`s'}"
    }

    return scalar singlecons = `singlecons'
end

capture program drop DiagonalProjection
program DiagonalProjection
    syntax anything, [_plugin_skip _plugin_bench]

    tokenize `anything'
    local A:  copy local 1
    local IV: copy local 2

    tempfile info
    tempname nabsorb
    mata st_numscalar("`nabsorb'", `IV'.nabsorb)
    if ( (`=scalar(`nabsorb')' > 1) & ("`_plugin_skip'" == "") ) {
        mata `IV'.exportc("`info'")
        cap noi plugin call manyiv_plugin, `"_plugin_run"' `"`info'"' `"`_plugin_bench'"'
        plugin_error_dispatcher `=_rc'
        mata `IV'.importc("`info'")
    }

    mata st_numscalar("`nabsorb'", `A'.nabsorb)
    if ( (`=scalar(`nabsorb')' > 1) & ("`_plugin_skip'" == "") ) {
        mata `A'.exportc("`info'")
        cap noi plugin call manyiv_plugin, `"_plugin_run"' `"`info'"' `"`_plugin_bench'"'
        plugin_error_dispatcher `=_rc'
        mata `A'.importc("`info'")
    }
    cap erase `"`info'"'
end

capture program drop ReadVarlist
program ReadVarlist, rclass
    syntax anything, [ivcall(str) cluster(str) absorb(str) absorbiv(str)]

    tokenize `anything'
    local touse:       copy local 1
    local C:           copy local 2
    local Y:           copy local 3
    local X:           copy local 4
    local W:           copy local 5
    local Z:           copy local 6
    local method_code: copy local 7
    local singlecons:  copy local 8
    local constant:    copy local 9

    ParseIVSyntax `ivcall', absorbiv(`absorbiv')

    local cons = ("`constant'" == "")

    mata `C' = ManyIVreg_Absorb_New(tokens(st_local("cluster")), "`touse'", `method_code')
    mata `Y' = st_data(., "`ivdepvar'", "`touse'")
    mata `X' = st_data(., "`ivendog'",  "`touse'")

    fvexpand `ivexog'
    local ivexog `r(varlist)'

    if ( `:list sizeof ivexog' ) {
        helper_strip_omitted `ivexog' if `touse', `constant'
        mata `W' = select(st_data(., "`ivexog'", "`touse'"), !st_matrix("r(omit)"))
        local ivexog `r(varlist)'
    }
    else {
        mata `W' = J(rows(`Y'), 0, .)
    }

    if ( `:list sizeof ivinst' ) {
        fvexpand `ivinst'
        local ivinst `r(varlist)'
        local ivinst: list ivinst - ivexog

        helper_strip_omitted `ivinst' if `touse', extra(`ivexog') `constant'
        mata `Z' = select(st_data(., "`ivinst'", "`touse'"), !st_matrix("r(omit)"))
        local ivinst `r(varlist)'
    }
    else {
        mata `Z' = J(rows(`Y'), 0, .)
    }

    if ( `cons' & (("`absorb'" == "") | `=scalar(`singlecons')') ) {
        mata `W' = `W', J(rows(`Y'), 1, 1)
    }
    else local cons = 0

    tempname kendog
    mata: st_numscalar("`kendog'", cols(st_data(1, "`ivendog'")))

    if ( `=scalar(`kendog')' != 1 ) {
        disp as error "Only one endogenous variabe allowed"
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

    local problems
    local ivallvars = `" `ivdepvar' `ivendog' `ivinst' `ivexog' "'
    local ivallvars: subinstr local ivallvars " " "  "
    foreach absvar in `absorb' `absorbiv' {
        if ( strpos("`ivallvars'", " `absvar' ") ) {
            local problems `problems' `absvar'
        }
    }

    if ( `"`problems'"' != `""' ) {
        disp as error "{bf:warning:} included as both absorb and non-absorb: `problems'"
    }

    return local cons: copy local cons

    c_local ivdepvar: copy local ivdepvar
    c_local ivendog:  copy local ivendog
    c_local ivinst:   copy local ivinst
    c_local ivexog:   copy local ivexog
end

***********************************************************************
*                                                                     *
*                               Helpers                               *
*                                                                     *
***********************************************************************

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
    syntax anything(equalok) [if], [extra(str) *]
    _rmcoll `anything' `extra' `if', expand `options'
    local expanded `r(varlist)'

    tempname b omit final
    matrix `b' = J(1, `:list sizeof expanded', .)
    matrix colnames `b' = `expanded'
    matrix `b' = `b'[1,1..(`:list sizeof expanded' - `:list sizeof extra')]

    _ms_omit_info `b'
    matrix `omit' = r(omit)
    mata `final' = select(st_matrixcolstripe("`b'")[., 2]', !st_matrix("r(omit)"))
    mata st_local("varlist", invtokens(cols(`final')? `final': ""))

    return local expanded: copy local expanded
    return local varlist:  copy local varlist
    return matrix omit =  `omit'
end

capture program drop FreeTimer
program FreeTimer
    qui {
        timer list
        local i = 99
        while ( (`i' > 0) & ("`r(t`i')'" != "") ) {
            local --i
        }
    }
    c_local FreeTimer `i'
end

capture program drop manyiv_timer
program manyiv_timer, rclass
    syntax anything, [prints(int 0) end off]
    tokenize `"`anything'"'
    local what  `1'
    local timer `2'
    local msg   `"`3'; "'

    * If timer is 0, then there were no free timers; skip this benchmark
    if ( `timer' == 0 ) exit 0

    if ( inlist("`what'", "start", "on") ) {
        cap timer off `timer'
        cap timer clear `timer'
        timer on `timer'
    }
    else if ( inlist("`what'", "info") ) {
        timer off `timer'
        qui timer list
        return scalar t`timer' = `r(t`timer')'
        return local pretty`timer' = trim("`:di %21.4gc r(t`timer')'")
        if ( `prints' ) di `"`msg'`:di trim("`:di %21.4gc r(t`timer')'")' seconds"'
        timer off `timer'
        timer clear `timer'
        timer on `timer'
    }

    if ( "`end'`off'" != "" ) {
        timer off `timer'
        timer clear `timer'
    }
end

capture program drop plugin_error_dispatcher
program plugin_error_dispatcher
    args rc
    if `rc' == 1701 disp as err "plugin error: fixed effects matrix inversion failed"
    if `rc' == 1702 disp as err "plugin error: unable to allocate plugin obbject (out of memory)"
    if `rc' == 1703 disp as err "plugin error: unable to export results back to Stata from plugin"
    if `rc' disp as err "(performance warning: plugin failed to run; will fall back on pure Stata code)"
    exit `rc'
end

* cap findfile manyiv_absorb.mata
* cap do `"`r(fn)'"'
*
* cap findfile manyiv_internals.mata
* cap do `"`r(fn)'"'
*
* cap findfile "manyiv_`c_os_'.plugin"
* cap do `"`r(fn)'"'

if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
else local c_os_: di lower("`c(os)'")

cap program drop manyiv_plugin
cap program manyiv_plugin, plugin using("manyiv_`c_os_'.plugin")
