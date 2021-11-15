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
    global dropbox ~/Dropbox/GPH_ExaminerDesign/Applications/DGPY
    use "${dropbox}/Data/paulData",clear
    egen office_id = group(office)
    egen judge_id  = group(judge)
    qui prep_data

    global controls r_a_roundag_* homeflag2_neg1  income findex_std_neg1  revtotbl_tc_neg1  agcolthb_tc_neg1 mortin12_ind_neg1 trdbalnm_tc_neg1 autopv6_ind_neg1 autbalt_tc_neg1 revratio_tc_neg1  nonmortinq6_tc_neg1  missing_age_0 missing_neg1  missing_zip_income
    local var findex_0_4_std
    set rmsg on

    ** Table 3, Column 2 and 3 of Row 1 (Panel A)
    ** Online Appendix Table 7, Column 5
    ** (TSLS matches the commented out results below)
    manyiv `var'           (discharge = z_ij_pooled) if sample == 1,  absorb(absorbid2) cluster(office_id)
    manyiv `var' $controls (discharge = z_ij_pooled) if sample == 1,  absorb(absorbid2) cluster(office_id)
    manyiv `var' $controls (discharge = .)           if sample == 1,  absorb(absorbid2) cluster(office_id) absorbiv(judge)
    manyiv `var' $controls (discharge = i.judge_id)  if sample == 1,  absorb(absorbid2) cluster(office_id)

    * Other results
    * manyiv score_0_4           (discharge = z_ij_pooled) if sample == 1,  absorb(absorbid2) cluster(office_id) absorbiv(judge)
    * manyiv score_0_4 $controls (discharge = z_ij_pooled) if sample == 1,  absorb(absorbid2) cluster(office_id) absorbiv(judge)
    * manyiv score_0_4 $controls (discharge = .)           if sample == 1,  absorb(absorbid2) cluster(office_id) absorbiv(judge)
    * manyiv score_0_4           (discharge = .)           if sample == 1,  absorb(absorbid2) cluster(office_id) absorbiv(judge)

    * TODO xx you are here; make sepparate test file for new task I think? meh, but def have output be more automated
    *
    * Z = judge
    * W = absorbid2
    * excluded, Z, Z:W
    * included, W
    * gegen absorbid2_judge = group(absorbid2 judge)
    * manyiv `var' $controls (discharge = .) if sample == 1, absorb(absorbid2) cluster(office_id) absorbiv(judge absorbid2_judge)

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

main
