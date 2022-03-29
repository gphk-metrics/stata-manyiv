* mata mata clear
* cap noi net uninstall manyiv
* net install manyiv, from(`c(pwd)'/../src/)
*
* mata mata clear
* qui include ../src/mata/manyiv_absorb.mata
* qui include ../src/mata/manyiv_internals.mata
* qui include ../src/ado/manyiv.ado

capture program drop main
program main
    load_data
    replicate

    manyoutcomes
    manyoutcomes_export
end

capture program drop load_data
program load_data
    * global dropbox ~/Dropbox/Papers/GPH_ExaminerDesign/Applications/DGPY
    global dropbox ~/Dropbox/GPH_ExaminerDesign/Applications/DGPY
    use "${dropbox}/Data/paulData",clear
    gegen office_id = group(office)
    gegen judge_id  = group(judge)
    qui prep_data
    global controls r_a_roundag_* homeflag2_neg1  income findex_std_neg1  revtotbl_tc_neg1  agcolthb_tc_neg1 mortin12_ind_neg1 trdbalnm_tc_neg1 autopv6_ind_neg1 autbalt_tc_neg1 revratio_tc_neg1  nonmortinq6_tc_neg1  missing_age_0 missing_neg1  missing_zip_income
    gegen absorbid2_judge = group(absorbid2 judge)
end

capture program drop replicate
program replicate
    local var findex_0_4_std
    set rmsg on

    ** Table 3, Column 2 and 3 of Row 1 (Panel A)
    ** Online Appendix Table 7, Column 5
    ** (TSLS matches the commented out results below)
    manyiv `var'           (discharge = z_ij_pooled) if sample == 1,  absorb(absorbid2) cluster(office_id)
    manyiv `var' $controls (discharge = z_ij_pooled) if sample == 1,  absorb(absorbid2) cluster(office_id)
    manyiv `var' $controls (discharge = .)           if sample == 1,  absorb(absorbid2) cluster(office_id) absorbiv(judge)
    manyiv `var' $controls (discharge = i.judge_id)  if sample == 1,  absorb(absorbid2) cluster(office_id)

    * ** Table 3, Column 2 and 3 of Row 1 (Panel A)
    * qui reghdfe `var' (discharge = z_ij_pooled) if sample == 1,  absorb(absorbid2) vce(cluster office_id) old
    * disp _b[discharge], _se[discharge]
    * qui reghdfe `var' $controls (discharge = z_ij_pooled) if sample == 1  ,   absorb(absorbid2) vce(cluster office_id)  old
    * disp _b[discharge], _se[discharge]
    *
    * ** Table 2, Column 1: First stage
    * qui reghdfe discharge  z_ij_pooled if sample == 1,  absorb(absorbid2) vce(cluster office_id) old
    * disp _b[z_ij_pooled], _se[z_ij_pooled]
    *
    * * As the replication code suggests the SE are off because of
    * * reghdfe; but here you can just use areg
    * qui areg discharge  z_ij_pooled if sample == 1,  absorb(absorbid2) vce(cluster office_id)
    * disp _b[z_ij_pooled], _se[z_ij_pooled]
    *
    * ** Online Appendix Table 7, Column 5
    * xi, pre(j_) i.judge
    * local var  findex_0_4_std
    * qui reghdfe `var' $controls (discharge = j_*) if sample == 1  ,   absorb(absorbid2) vce(cluster office_id)  old
    * disp _b[discharge], _se[discharge]
end

* ---------------------------------------------------------------------
* Coefficient comparison

capture program drop manyoutcomes
program manyoutcomes
    keep if sample == 1
    gegen leniency_jive_sum    = sum(discharge),    by(judge)
    gegen leniency_jive_n      = count(discharge),  by(judge)
    gen leniency_iv            = leniency_jive_sum / leniency_jive_n
    gen leniency_jive          = (leniency_jive_sum - discharge) / (leniency_jive_n - 1)
    gegen leniency_oym_sum     = sum(discharge),   by(office year file_month)
    gegen leniency_oym_n       = count(discharge), by(office year file_month)
    gen leniency_ujive         = leniency_jive - (leniency_oym_sum - discharge) / (leniency_oym_n - 1)
    gegen leniency_jive_w_sum  = sum(discharge),   by(judge office file_year file_month)
    gegen leniency_jive_w_n    = count(discharge), by(judge office file_year file_month)
    gen leniency_iv_w          = leniency_jive_w_sum / leniency_jive_w_n
    gen leniency_jive_w        = (leniency_jive_w_sum - discharge) / (leniency_jive_w_n - 1)
    gen leniency_ujive_w       = leniency_jive_w - (leniency_oym_sum - discharge) / (leniency_oym_n - 1)

    local vars             ///
        findex_0_4_std     ///
        cur2gn12_ind_neg1  ///
        collec12_ind_neg1  ///
        chgoff12_ind_neg1  ///
        bankin12_ind_neg1  ///
        forecl12_ind_neg1  ///
        jdgst12h_ind_neg1  ///
        lien12_ind_neg1    ///
        repo12_ind_neg1    ///
        revtotbl_0_4_tc    ///
        agcolthb_0_4_tc    ///
        mortin12_ind_0_4   ///
        trdbalnm_0_4_tc    ///
        autopv6_ind_0_4    ///
        autbalt_0_4_tc     ///
        revratio_0_4_tc    ///
        nonmortinq6_0_4_tc ///
        score_0_4

    local vars findex_0_4_std revtotbl_0_4_tc agcolthb_0_4_tc mortin12_ind_0_4 trdbalnm_0_4_tc autopv6_ind_0_4 autbalt_0_4_tc  revratio_0_4_tc  nonmortinq6_0_4_tc
    local appendix_vars deadflag2_cum_4 rellien_0_4 paidjudg_0_4 unpdjudg_0_4 paidcol_0_4 unpdcol_0_4 colm12h_0_4 cmop4up_0_4 cmop5up_0_4 studact2_0_4 studdef2_0_4

    local opts       absorb(absorbid2) cluster(office)
    local judge_fe   `opts' absorbiv(judge)
    local judge_fe_w `opts' absorbiv(absorbid2_judge judge)
    mata out = J(0, 28, .)
    foreach var of local vars {
        disp "`var'"
        qui manyiv `var' (discharge = leniency_iv)      if sample == 1, save(res_iv)      `opts'
        qui manyiv `var' (discharge = z_ij_pooled)      if sample == 1, save(res_ujive2)      `opts'
        qui manyiv `var' (discharge = leniency_jive)    if sample == 1, save(res_jive)    `opts'
        qui manyiv `var' (discharge = leniency_ujive)   if sample == 1, save(res_ujive)   `opts'
        qui manyiv `var' (discharge = .)                if sample == 1, save(res_fe)      `judge_fe'
        qui manyiv `var' (discharge = leniency_iv_w)    if sample == 1, save(res_iv_w)    `opts'
        qui manyiv `var' (discharge = leniency_jive_w)  if sample == 1, save(res_jive_w)  `opts'
        qui manyiv `var' (discharge = leniency_ujive_w) if sample == 1, save(res_ujive_w) `opts'
        qui manyiv `var' (discharge = .)                if sample == 1, save(res_fe_w)    `judge_fe_w'

        mata row =
             res_iv.beta[1]      , res_iv.se[3, 1]     ,
             res_ujive2.beta[2]  , res_ujive2.se[3, 2] ,
             res_iv.beta[2]      , res_iv.se[3, 2]     ,
             res_jive.beta[2]    , res_jive.se[3, 2]   ,
             res_ujive.beta[2]   , res_ujive.se[3, 2]  ,
             res_fe.beta[2]      , res_fe.se[3, 2]     ,
             res_fe.beta[5]      , res_fe.se[3, 5]     ,
             res_fe.beta[6]      , res_fe.se[3, 6]     ,
             res_iv_w.beta[2]    , res_iv_w.se[3, 2]   ,
             res_jive_w.beta[2]  , res_jive_w.se[3, 2] ,
             res_ujive_w.beta[2] , res_ujive_w.se[3, 2],
             res_fe_w.beta[2]    , res_fe_w.se[3, 2]   ,
             res_fe_w.beta[5]    , res_fe_w.se[3, 5]   ,
             res_fe_w.beta[6]    , res_fe_w.se[3, 6]
        mata out = out \ row
    }

    gstats tab `vars', s(sd) mata(sd)
    mata reix = colshape(rowshape(1::(2 * `:list sizeof vars')', 2)', 1)
    mata save_table("../misc/tables.txt", "varegs", colshape(colshape(out, 2)', 14)[reix, .])
end

capture program drop manyoutcomes_export
program manyoutcomes_export
    preserve
        clear
        mata stata(sprintf("set obs %g", rows(out)))

        local varname `"."'                      ///
                      `"Financial Strain Index"' ///
                      `"Revolving Balance"'      ///
                      `"Collection Balance"'     ///
                      `"Have a Mortgage"'        ///
                      `"Mortgage Balance"'       ///
                      `"Have an Auto Loan"'      ///
                      `"Auto Balance"'           ///
                      `"Revolving Utilization"'  ///
                      `"Non-Mortgage Inquiries"'

        mata varname = tokens(st_local("varname"))[2::`:list sizeof varname']
        mata st_addvar(sprintf("str%g", max(strlen(varname))), "varname")
        mata (void) st_sstore(., "varname", varname')

        mata st_addvar("`:set type'", "sd")
        mata (void) st_store(., "sd", sd.output)

        local types iv jive ujive iv_fe jive_fe ujive_fe
        local varnames
        foreach est_type in ols paper_no_w {
            local varnames `varnames' beta_`est_type' se_`est_type'
        }

        foreach w in no_w w {
            foreach est_type of local types {
                local varnames `varnames' beta_`est_type'_`w' se_`est_type'_`w'
            }
        }

        mata (void) st_addvar(J(1, cols(out), "`:set type'"), tokens("`varnames'"))
        mata st_store(., tokens("`varnames'"), out)

        foreach est_type in ols paper_no_w {
            gen t_`est_type' = beta_`est_type' / se_`est_type'
        }

        foreach w in no_w w {
            foreach est_type of local types {
                gen t_`est_type'_`w' = beta_`est_type'_`w' / se_`est_type'_`w'
            }
        }
        outsheet using "../misc/graph_input.csv", comma replace names
    restore
end

* ---------------------------------------------------------------------
* Replication

capture program drop prep_data
program define prep_data
    gen roundage_0 = round(age_0, 5)
    xi, pre(r_a_) i.roundage_0

    ** creating baseline financial strain index with 8 components
    foreach var in forecl12_ind_neg1 repo12_ind_neg1 chgoff12_ind_neg1 lien12_ind_neg1 collec12_ind_neg1 jdgst12h_ind_neg1 bankin12_ind_neg1 cur2gn12_ind_neg1{
        sum `var' if discharge==0 & missing_neg1==0
        local mean_`var' = string(r(mean), "%12.3f")
        local sd_`var' = string(r(sd), "%12.3f")
    }

    foreach var in forecl12_ind_neg1 repo12_ind_neg1 chgoff12_ind_neg1 lien12_ind_neg1 collec12_ind_neg1 jdgst12h_ind_neg1 bankin12_ind_neg1 cur2gn12_ind_neg1 {
        gen `var'_std = (`var' - `mean_`var'')/(`sd_`var'') if missing_neg1==0
    }

    gen findex_neg1 = (forecl12_ind_neg1_std +  repo12_ind_neg1_std +  chgoff12_ind_neg1_std +  lien12_ind_neg1_std +  collec12_ind_neg1_std +  jdgst12h_ind_neg1_std +  bankin12_ind_neg1_std +  cur2gn12_ind_neg1_std)/8
    sum findex_neg1 if discharge==0 & missing_neg1==0
        local mean_findex_neg1 = string(r(mean), "%12.3f")
        local sd_findex_neg1 = string(r(sd), "%12.3f")
    gen findex_std_neg1 = (findex_neg1 - `mean_findex_neg1')/(`sd_findex_neg1')

    ** creating baseline financial strain index with 6 components
    gen findex_neg1_6 = (forecl12_ind_neg1_std + repo12_ind_neg1_std + lien12_ind_neg1_std + collec12_ind_neg1_std +  jdgst12h_ind_neg1_std + cur2gn12_ind_neg1_std)/6
    sum findex_neg1_6 if discharge==0 & missing_neg1==0
        local mean_findex_neg1_6 = string(r(mean), "%12.3f")
        local sd_findex_neg1_6 = string(r(sd), "%12.3f")
    gen findex_std_neg1_6 = (findex_neg1_6 - `mean_findex_neg1_6')/(`sd_findex_neg1_6')

    foreach var in findex_std  {
        replace `var'_neg1 = 0 if missing_neg1==1
    }
    foreach var in forecl12_ind_neg1 repo12_ind_neg1 chgoff12_ind_neg1 lien12_ind_neg1 collec12_ind_neg1 jdgst12h_ind_neg1 bankin12_ind_neg1 cur2gn12_ind_neg1{
        sum `var' if discharge==0 & missing_neg1==0
        local mean_`var' = string(r(mean), "%12.3f")
        local sd_`var' = string(r(sd), "%12.3f")
    }

    foreach var in forecl12_ind repo12_ind  lien12_ind collec12_ind jdgst12h_ind cur2gn12_ind chgoff12_ind bankin12_ind {
        foreach timeperiod in neg2 0 1 2 3 4 5 6 7 {
            gen `var'_`timeperiod'_std = (`var'_`timeperiod' - `mean_`var'_neg1')/(`sd_`var'_neg1')
        }
    }

    gen findex_0 = (forecl12_ind_0_std + repo12_ind_0_std + lien12_ind_0_std + collec12_ind_0_std + jdgst12h_ind_0_std +  cur2gn12_ind_0_std)/6
    sum findex_neg1_6 if discharge==0 & missing_neg1==0
        local mean_findex_neg1_6 = string(r(mean), "%12.3f")
        local sd_findex_neg1_6 = string(r(sd), "%12.3f")
    gen findex_0_std = (findex_0 - `mean_findex_neg1_6')/(`sd_findex_neg1_6')

    foreach timeperiod in neg2 1 2 3 4 5 6 7 {
        gen findex_`timeperiod' = (forecl12_ind_`timeperiod'_std +  repo12_ind_`timeperiod'_std +  chgoff12_ind_`timeperiod'_std +  lien12_ind_`timeperiod'_std +  collec12_ind_`timeperiod'_std +  jdgst12h_ind_`timeperiod'_std +  bankin12_ind_`timeperiod'_std +  cur2gn12_ind_`timeperiod'_std)/8
        sum findex_neg1 if discharge==0 & missing_neg1==0
            local mean_findex_neg1 = string(r(mean), "%12.3f")
            local sd_findex_neg1 = string(r(sd), "%12.3f")
        gen findex_`timeperiod'_std = (findex_`timeperiod' - `mean_findex_neg1')/(`sd_findex_neg1')
    }
    egen findex_0_4_std = rowmean(findex_0_std findex_1_std findex_2_std findex_3_std findex_4_std)
    egen findex_5_7_std = rowmean(findex_5_std findex_6_std findex_7_std)
end

capture mata: mata drop save_table()
mata:
void function save_table(
    string scalar output,
    string scalar tag,
    real matrix M)
{
    real scalar i, j, fh
    string scalar fmt

    fmt = "\t%21.9f"
    fh  = fopen(output, "a")
    fwrite(fh, sprintf("<tab:%s>\n", tag))
    for(i = 1; i <= rows(M); i++) {
        for(j = 1; j <= cols(M); j++) {
            fwrite(fh, sprintf(fmt, M[i, j]))
        }
        fwrite(fh, sprintf("\n"))
    }
    fclose(fh)
}
end

main
