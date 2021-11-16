* mata mata clear
* cap noi net uninstall manyiv
* net install manyiv, from(`c(pwd)'/../src/)
*
* mata mata clear
* qui include ../src/mata/manyiv_absorb.mata
* qui include ../src/mata/manyiv_internals_m.mata
* qui include ../src/ado/manyiv.ado

capture program drop main
program main
    load_data
    replicate
    manyoutcomes
end

capture program drop load_data
program load_data
    global dropbox ~/Dropbox/GPH_ExaminerDesign/Applications/DGPY
    use "${dropbox}/Data/paulData",clear
    egen office_id = group(office)
    egen judge_id  = group(judge)
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

capture program drop manyoutcomes
program manyoutcomes
    tempvar sumj obsj sumc obsc
    gegen `sumj' = sum(discharge),   by(judge office)
    gegen `obsj' = count(discharge), by(judge office)
    gegen `sumc' = sum(discharge),   by(office) replace
    gegen `obsc' = count(discharge), by(office) replace

    * cap drop z_ij_pooled
    gen z_ij_pooled_JIVE = (`sumj' - discharge) / (`obsj' - 1)
    gen z_ij_pooled_IV   = `sumj' / `obsj'
    * gen z_ij_pooled      = z_ij_pooled_JIVE - (`sumc' - discharge) / (`obsc' - 1)

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

    local optsfe absorb(absorbid2) cluster(office_id) absorbiv(absorbid2_judge judge)
    local opts   absorb(absorbid2) cluster(office_id)
    mata out_pooled  = J(0, 6, .)
    mata out_judgefe = J(0, 6, .)
    mata out_iv      = J(0, 6, .)
    mata out_jiv     = J(0, 6, .)
    foreach var of local vars {
        disp "`var'"
        qui manyiv `var' (discharge = z_ij_pooled)      if sample == 1, save(res_pooled)  `opts'
        qui manyiv `var' (discharge = .)                if sample == 1, save(res_judgefe) `optsfe'
        qui manyiv `var' (discharge = z_ij_pooled_JIVE) if sample == 1, save(res_iv)      `opts'
        qui manyiv `var' (discharge = z_ij_pooled_IV)   if sample == 1, save(res_jiv)     `opts'
        // qui manyiv `var' $controls (discharge = z_ij_pooled)      if sample == 1, save(res_pooled)  `opts'
        // qui manyiv `var' $controls (discharge = .)                if sample == 1, save(res_judgefe) `optsfe'
        // qui manyiv `var' $controls (discharge = z_ij_pooled_JIVE) if sample == 1, save(res_iv)      `opts'
        // qui manyiv `var' $controls (discharge = z_ij_pooled_IV)   if sample == 1, save(res_jiv)     `opts'
        mata out_pooled  = out_pooled  \ rowshape((res_pooled.beta[(2 \ 5 \ 6)]', res_pooled.se[3, (2 \ 5 \ 6)]'), 1)
        mata out_judgefe = out_judgefe \ rowshape((res_judgefe.beta[(2 \ 5 \ 6)]', res_judgefe.se[3, (2 \ 5 \ 6)]'), 1)
        mata out_iv      = out_iv      \ rowshape((res_iv.beta[(2 \ 5 \ 6)]', res_iv.se[3, (2 \ 5 \ 6)]'), 1)
        mata out_jiv     = out_jiv     \ rowshape((res_jiv.beta[(2 \ 5 \ 6)]', res_jiv.se[3, (2 \ 5 \ 6)]'), 1)
    }

    mata save_table("../misc/tables.txt", "varegs", (out_pooled, out_judgefe) \ (out_iv, out_jiv))
end

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
