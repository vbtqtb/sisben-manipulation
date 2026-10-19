*******************************************************************************

* PART 0: PRELIMINARY 

*******************************************************************************
version 17
clear all                
set more off             
set linesize 120         
capture log close        

* PATHS
* Set ROOT to the folder where you cloned this repository.
* Put the original data (sisben_aejep.dta) into data/ (see README).

global ROOT "."
global DATA "$ROOT/data"
global OUT  "$ROOT/output"
capture mkdir "$OUT"


* DECIDE WHICH TO RUN BC OF COMPUTATIONAL LIMIT 
* 0: dont run, 1: run

global RUN_A       1   // slow
global RUN_YEARLY  1   // Part C: yearly tests (a few minutes per year).
global RUN_CELLS   0   // one CJM test per municipality-year, slow

					   
* estimation settings
					   
global CUT        48   // eligibility cutoff on continuous score
global MIN_CELL  250   // only estimate cells with at least MIN_CELL
                       //   HHs (estimator needs enough data)
global CELL_CAMPAIGN 1 // 1 = use only interviews from May-Oct,
                       //   0 = use all 12 months
				
* rddensity, lpdensity, estout are needed for the estimations
* automatically installs missing packages

foreach pkg in rddensity lpdensity estout {
    capture which `pkg'
    if _rc != 0 {
        display as text "Install `pkg'"
        ssc install `pkg', replace
    }
}



log using "$OUT/cc_partial_replication.log", replace text

*******************************************************************************

* PART 1: SAMPLE CREATION & COLLAPSE TO FAMILY LEVEL
* (slow)
*******************************************************************************

if $RUN_A == 1 {

    use "$DATA/sisben_aejep.dta", clear

    * construct score from components and track difference to recorded score
    * authors: table6 do file
    gen double score_cc = facviv + facedu + facpat + facdem
    gen double diff     = score_cc - puntaje

	* should be in 0-1 if puntaje is rounded down
    summarize diff, detail

	* drop anomaly records, 4900 records ca. 0.03%
    count if abs(diff) > 1
    drop if abs(diff) > 1

    * uhhid2 is assumed to be family ID, maybe other definitions in diff dataset
	* collapse to family level for score 
    bysort uhhid2: gen byte fam_tag = (_n == 1) // first person in family
    keep if fam_tag == 1
    drop fam_tag
    count                               

    gen int  year  = year(fencuesta)
    gen byte month = month(fencuesta)
    keep if inrange(year, 1994, 2003)

    * for computational runtime
    keep if inrange(score_cc, 33, 63)

    keep codunico year month score_cc puntaje
    label var score_cc "Constructed poverty index score"
    label var puntaje  "Recorded poverty index score"
    compress                            
    save "$OUT/sample.dta", replace
}



*******************************************************************************

* PART 2: PLOT SCORE DISTRIBUTION FOR ALL YEARS 

*******************************************************************************

use "$OUT/sample.dta", clear

* check data by family per year
tab year

* check if puntaje is rounded down score_cc
count if puntaje != floor(score_cc)

*only 33 - 63 score window is plotted
histogram puntaje if inrange(puntaje, 33, 62), discrete percent        ///
    by(year, cols(5) note("")                                          ///
       title("Poverty index score by interview year"))                 ///
    xline(47.5, lcolor(green))                                           ///
    xtitle("Poverty index score (integer)")                            ///
    ytitle("Percent of families")
graph export "$OUT/fig1_score_by_year.png", replace width(2000)

*******************************************************************************

* PART 3: REPLICATION OF DISCONTINUITY ESTIMATE AND EXTENSION WITH CJM

*******************************************************************************

if $RUN_YEARLY == 1 {

*******************************************************************************

* PART 3.1: REPLICATION OF DISCONTINUITY ESTIMATE 

*******************************************************************************

* counts per year and score, each for person-level and household-level
use uhhid2 puntaje fencuesta zona estrato using "$DATA/sisben_aejep.dta", clear
keep if zona == 1 & estrato <= 3           
gen int year = year(fencuesta)
keep if inrange(year, 1994, 2003)
drop if missing(puntaje)

preserve
	* how many people with one score per year
    contract year puntaje, freq(n_ind)
    tempfile ind_counts
    save `ind_counts'
restore

* collapse to family level
bysort uhhid2: keep if _n == 1
contract year puntaje, freq(n_fam)

* Put both counts side by side (1:1 = each year x score appears once).
merge 1:1 year puntaje using `ind_counts', nogenerate
tempfile counts
save `counts'

* create author's estimator per unit and per year
tempname t3
tempfile t3res
postfile `t3' str3 unit int year double est double se double bw ///
    using `t3res', replace

foreach u in ind fam {

    use `counts', clear
    keep year puntaje n_`u'
    drop if missing(n_`u')           // year x score cells with no one
    rename (puntaje n_`u') (score survey)

    * wide format like authors: one row per score, columns survey1994 etc
    reshape wide survey, i(score) j(year)

    * bar heights: percent of the year's total at each score
    forvalues i = 1994/2003 {
        egen sisses`i' = pc(survey`i')
    }

    * running variable and polynoms
    gen align      = score - 47
    gen jump       = (score <= 47)        // 1 = eligible side
    gen jump_align = jump * align
    gen align2     = align^2
    gen align3     = align^3
    gen align4     = align^4

    * distance from cutoff to the lowest / highest score in the data
    quietly summarize align if align <= 0
    local predist = -r(min)
    quietly summarize align if align > 0
    local postdist = r(max)

    forvalues i = 1994/2003 {

        capture drop presd presd2 postsd postsd2 k

        * rule-of-thumb bandwidth, left side 
        quietly regress sisses`i' align align2 align3 align4 if align <= 0
        local mse_pre = e(rmse)^2
        * 2nd derivative of quartic term at every left-side score, squared,
        * then summed (r(sum) is the total of the variable)
        gen presd  = 2*_b[align2] + 6*_b[align3]*align ///
                   + 12*_b[align4]*align2 if align <= 0
        gen presd2 = presd^2
        quietly summarize presd2
        local pressd = r(sum)
        local prebw  = 3.348 * ((`mse_pre' * `predist') / `pressd')^0.2

        * right side
        quietly regress sisses`i' align align2 align3 align4 if align > 0
        local mse_post = e(rmse)^2
        gen postsd  = 2*_b[align2] + 6*_b[align3]*align ///
                    + 12*_b[align4]*align2 if align > 0
        gen postsd2 = postsd^2
        quietly summarize postsd2
        local postssd = r(sum)
        local postbw  = 3.348 * ((`mse_post' * `postdist') / `postssd')^0.2

        * final bandwidth = average of both sides
        local h = float(0.5 * (`prebw' + `postbw'))

        * table 3
        * triangular kernel, ie weight 1 at score 47, 0 beyond bandwidth
        gen k = max(0, 1 - abs(align/`h'))
        quietly regress sisses`i' jump align jump_align [aweight = k], vce(robust)

        display as text "Table 3 method, unit `u', `i': b = " ///
            as result %6.3f _b[jump] as text "  (se " ///
            as result %5.3f _se[jump] as text ")  bw = " as result %5.2f `h'

        post `t3' ("`u'") (`i') (_b[jump]) (_se[jump]) (`h')
    }
}
postclose `t3'

use `t3res', clear
label var est "Table 3 method: jump at 47, percentage points"
label var se  "Robust standard error"
label var bw  "Rule-of-thumb bandwidth (score points)"
save "$OUT/yearly_table3.dta", replace

*******************************************************************************

* PART 3.2: CJM ESTIMATION

*******************************************************************************

use "$OUT/sample.dta", clear

tempname cj
tempfile cjres
postfile `cj' int year int cutoff byte bw                            ///
    double f_left double f_right double se_left double se_right      ///
    double h_left double h_right double n_left double n_right        ///
    double logratio double se_logratio double pv_official            ///
    using `cjres', replace

* outer loop: cutoffs
* middle loop: bandwidth (0 = automatic, 3 = fixed).
* inner loop: years
* with artificial thresholds

foreach c of numlist $CUT 42 54 {
    foreach bw of numlist 0 3 {

        * fixed bandwidth only for 47
        if `c' != $CUT & `bw' == 3 continue

        * hopt saves bandwidth option, either empty (automatic) or h(3 3)
        if `bw' == 0 local hopt ""
        else         local hopt "h(`bw' `bw')"

        forvalues y = 1994/2003 {

            foreach v in fl fr sl sr hl hr nl nr lr slr pv {
                local `v' = .
            }

            ereturn clear      

            * _rc is 0 if the command worked
            capture quietly rddensity score_cc if year == `y', c(`c') `hopt'

            if _rc == 0 {
                local fl = e(f_ql)
                local fr = e(f_qr)
                local sl = e(se_ql)
                local sr = e(se_qr)
                local hl = e(h_l)
                local hr = e(h_r)
                local nl = e(N_h_l)
                local nr = e(N_h_r)
                local pv = e(pv_q)

                * missing value counts as larger than any number, so check "< ." (= not missing)
                if `fl' > 0 & `fr' > 0 & `fl' < . & `fr' < . {
                    local lr  = ln(`fl') - ln(`fr')
                    local slr = sqrt( (`sl'/`fl')^2 + (`sr'/`fr')^2 )
                }
            }
            else {
                display as error "rddensity fail: year `y' cutoff `c' bw `bw'"
            }

            display as text "CJM  year `y'  cutoff `c'  bw `bw':  " ///
                as result %7.3f `lr' as text "  (se " as result %6.3f `slr' as text ")"

            post `cj' (`y') (`c') (`bw') (`fl') (`fr') (`sl') (`sr') ///
                (`hl') (`hr') (`nl') (`nr') (`lr') (`slr') (`pv')
        }
    }
}
postclose `cj'

use `cjres', clear

* 95% confidence interval and p-value for the log ratio
gen double ci_lo       = logratio - 1.96 * se_logratio
gen double ci_hi       = logratio + 1.96 * se_logratio
gen double pv_logratio = 2 * normal( -abs(logratio / se_logratio) )

label var bw          "Bandwidth: 0 = automatic, 3 = fixed 3 points"
label var logratio    "CJM: log f(below) - log f(above)"
label var se_logratio "Delta method standard error of the log ratio"
label var pv_official "p-value, rddensity's own test"
label var pv_logratio "p-value, log ratio = 0"
save "$OUT/yearly_cjm.dta", replace


*******************************************************************************

* PART 3.3: MAIN TABLE FOR SEMINAR PAPER FOR COMPARISON

*******************************************************************************

* CJM with fixed bandwidth
use "$OUT/yearly_cjm.dta", clear
keep if cutoff == $CUT & bw == 3
keep year logratio se_logratio
rename (logratio se_logratio) (cjm_h3 cjm_h3_se)
tempfile fix3
save `fix3'

* artificial cutoffs, one column each
use "$OUT/yearly_cjm.dta", clear
keep if cutoff != $CUT
keep year cutoff logratio
reshape wide logratio, i(year) j(cutoff)       // -> logratio42, logratio54
rename (logratio42 logratio54) (placebo42 placebo54)
tempfile plac
save `plac'

* authors' method, persons and HHs
use "$OUT/yearly_table3.dta", clear
keep unit year est se
reshape wide est se, i(year) j(unit) string   
rename (estind seind estfam sefam) (t3_ind t3_ind_se t3_fam t3_fam_se)
tempfile t3wide
save `t3wide'

* CJM results with above
use "$OUT/yearly_cjm.dta", clear
keep if cutoff == $CUT & bw == 0
keep year logratio se_logratio pv_logratio h_left h_right
rename (logratio se_logratio pv_logratio) (cjm_auto cjm_auto_se cjm_auto_p)
merge 1:1 year using `fix3', nogenerate
merge 1:1 year using `plac', nogenerate
merge 1:1 year using `t3wide', nogenerate

* authors' Table 3 
local paper_est 0.033 0.080 0.008 0.024 0.868 1.209 1.422 1.683 1.565 1.547
local paper_se  0.086 0.083 0.121 0.097 0.119 0.145 0.154 0.150 0.132 0.132
gen double paper_t3    = .
gen double paper_t3_se = .
forvalues i = 1/10 {
    local yr = 1993 + `i'                      // i = 1 is 1994, etc.
    local e : word `i' of `paper_est'          // number in position i in list
    local s : word `i' of `paper_se'
    replace paper_t3    = `e' if year == `yr'
    replace paper_t3_se = `s' if year == `yr'
}


order year paper_t3 paper_t3_se t3_ind t3_ind_se t3_fam t3_fam_se ///
      cjm_auto cjm_auto_se cjm_auto_p h_left h_right cjm_h3 cjm_h3_se ///
      placebo42 placebo54
format paper_t3-placebo54 %7.3f
list, noobs separator(0) abbreviate(12)

save "$OUT/table_yearly_discontinuity.dta", replace
export delimited using "$OUT/table_yearly_discontinuity.csv", replace


*******************************************************************************

* PART 3.4: CJM ESTIMATION WITH YEARS AND ARTIFICIAL THRESHOLDS

*******************************************************************************

use "$OUT/yearly_cjm.dta", clear
keep if bw == 0                                 // automatic bandwidth only

twoway ///
    (rcap ci_hi ci_lo year if cutoff == $CUT, lcolor(navy))                 ///
    (connected logratio year if cutoff == $CUT, lcolor(navy) mcolor(navy))  ///
    (connected logratio year if cutoff == 42, lpattern(dash)                ///
        lcolor(gs9) mcolor(gs9) msymbol(Th))                                ///
    (connected logratio year if cutoff == 54, lpattern(shortdash)           ///
        lcolor(gs6) mcolor(gs6) msymbol(Sh)),                               ///
    yline(0, lcolor(black))                                                 ///
    xline(1997.5, lpattern(dash) lcolor(green))                               ///
    xlabel(1994(1)2003)                                                     ///
    xtitle("Interview year")                                                ///
    ytitle("log density below - log density above")                         ///
    legend(order(2 "Eligibility cutoff" 1 "95% CI"                          ///
                 3 "Placebo cutoff 42" 4 "Placebo cutoff 54")               ///
           position(6) rows(2))                                             ///
    xsize(8) ysize(5)
graph export "$OUT/fig2_cjm_by_year.png", replace width(2000)

* end of: if $RUN_YEARLY == 1
}



*******************************************************************************

* PART 4:  MUNICIPALITY YEAR MEASURES

*******************************************************************************

*******************************************************************************

* PART 4.1:  BIN DISCONTINUITY MEASURE BY AUTHORS

*******************************************************************************

use codunico fencuesta estrato zona puntaje using "$DATA/sisben_aejep.dta", clear

rename puntaje score
gen int  doi_yr = year(fencuesta)      // interview year
gen int  doi_m  = month(fencuesta)     // interview month
gen byte count  = 1                    // each row counts as one interview

keep if inrange(doi_yr, 1994, 2003)

* count interviews per municipality x stratum x zone x year x month x score
collapse (sum) survey = count, by(codunico estrato zona doi_yr doi_m score)

* binary variables, looks if score in bin
gen byte d45_47 = inrange(score, 45, 47)
gen byte d48_50 = inrange(score, 48, 50)
gen byte d43_47 = inrange(score, 43, 47)
gen byte d48_52 = inrange(score, 48, 52)

*may to october range
keep if inrange(doi_m, 5, 10)

* counts of interviews in each bin (0 if score is outside bin)
gen double s45_47 = d45_47 * survey
gen double s48_50 = d48_50 * survey
gen double s43_47 = d43_47 * survey
gen double s48_52 = d48_52 * survey

tempfile empirics
save `empirics'

* before 1998 manipulation, per munic
* same as Table6.do (line 50-91)

keep if inrange(doi_yr, 1994, 1997)
gen byte enc = 1
collapse (sum) d45_47 d48_50 d43_47 d48_52 enc, by(codunico)
rename enc total
gen double jump3_base = d45_47/total - d48_50/total
gen double jump5_base = d43_47/total - d48_52/total
keep codunico jump3_base jump5_base
tempfile basejump
save `basejump'

* discontinuity per year and munic
use `empirics', clear
collapse (sum) survey s45_47 s48_50 s43_47 s48_52, by(codunico doi_yr)

foreach v in s45_47 s48_50 s43_47 s48_52 {
    replace `v' = . if `v' == 0
}

gen double jump3 = s45_47/survey - s48_50/survey
gen double jump5 = s43_47/survey - s48_52/survey

rename doi_yr year
merge m:1 codunico using `basejump', keep(master match) nogenerate

gen double jumpd3 = (jump3 - jump3_base) / 3
gen double jumpd5 = (jump5 - jump5_base) / 5

label var jump3  "Bin difference +/-3 (share below - share above)"
label var jump5  "Bin difference +/-5"
label var jumpd3 "Bin difference +/-3, net of 1994-97 baseline, /3"
label var jumpd5 "Bin difference +/-5, net of 1994-97 baseline, /5"

keep codunico year survey jump3 jump5 jumpd3 jumpd5
save "$OUT/measure_bindiff.dta", replace


*******************************************************************************

* PART 4.2: CJM ESTIMATOR 
* (very slow)
*******************************************************************************

* run in cells bc otherwise runtime is too long

if $RUN_CELLS == 1 {

    use "$OUT/sample.dta", clear

    if $CELL_CAMPAIGN == 1 keep if inrange(month, 5, 10)

    egen long cell = group(codunico year)

    sort cell                          
    gen long rowid = _n                 

    by cell: gen long start_i = rowid[1]    // first row of the cell's block
    by cell: gen long end_i   = rowid[_N]   // last row of the cell's block
    by cell: gen long n_cell  = _N          // number of families in the cell

    * each block must be as long as its cell
    assert end_i - start_i + 1 == n_cell

    * copy one row per cell into a matrix,
    * so loop below can read cell ID, municipality, year and range
    preserve
        by cell: keep if _n == 1
        keep if n_cell >= $MIN_CELL
        keep cell codunico year start_i end_i n_cell
        mkmat cell codunico year start_i end_i n_cell, matrix(CELLS)
        local ncells = rowsof(CELLS)
        display as result "Cells to estimate: `ncells'"
    restore

    tempname pf
    tempfile results
    postfile `pf' long cell long codunico int year long n_cell           ///
        double f_left double f_right double se_left double se_right      ///
        double logratio double se_logratio double h_left double h_right  ///
        double n_eff_l double n_eff_r                                    ///
        double logratio_fix double se_logratio_fix                       ///
        using `results', replace

    forvalues i = 1/`ncells' {

        local c   = CELLS[`i', 1]
        local mun = CELLS[`i', 2]
        local yr  = CELLS[`i', 3]
        local a   = CELLS[`i', 4]
        local b   = CELLS[`i', 5]
        local n   = CELLS[`i', 6]

        * the first and last row of the range must belong to cell c
        if cell[`a'] != `c' | cell[`b'] != `c' {
            display as error "Range `a'/`b' does not belong to cell `c'"
            exit 459
        }

        * progression text to display
        if mod(`i', 50) == 0 display as text "`i' of `ncells'"


        foreach v in fl fr sl sr lr slr hl hr nl nr lrf slf {
            local `v' = .
        }

        * case 1 : automatic bandwidth
        ereturn clear
        capture quietly rddensity score_cc in `a'/`b', c($CUT)
        if _rc == 0 {
            local fl = e(f_ql)
            local fr = e(f_qr)
            local sl = e(se_ql)
            local sr = e(se_qr)
            local hl = e(h_l)
            local hr = e(h_r)
            local nl = e(N_h_l)
            local nr = e(N_h_r)
            if `fl' > 0 & `fr' > 0 & `fl' < . & `fr' < . {
                local lr  = ln(`fl') - ln(`fr')
                local slr = sqrt( (`sl'/`fl')^2 + (`sr'/`fr')^2 )
            }
        }

        * case 2 : fixed bandwidth, 3 points on each side
        ereturn clear
        capture quietly rddensity score_cc in `a'/`b', c($CUT) h(3 3)
        if _rc == 0 {
            local f1 = e(f_ql)
            local f2 = e(f_qr)
            local s1 = e(se_ql)
            local s2 = e(se_qr)
            if `f1' > 0 & `f2' > 0 & `f1' < . & `f2' < . {
                local lrf = ln(`f1') - ln(`f2')
                local slf = sqrt( (`s1'/`f1')^2 + (`s2'/`f2')^2 )
            }
        }

        post `pf' (`c') (`mun') (`yr') (`n') (`fl') (`fr') (`sl') (`sr') ///
            (`lr') (`slr') (`hl') (`hr') (`nl') (`nr') (`lrf') (`slf')
    }
    postclose `pf'

    use `results', clear
    label var logratio     "CJM log ratio, automatic bandwidth"
    label var logratio_fix "CJM log ratio, bandwidth 3"
    label var n_cell       "Families in municipality-year cell"
    save "$OUT/measure_density_cells.dta", replace
}



*******************************************************************************

* PART 4.3:  COMPARISON OF TWO MEASURES

*******************************************************************************

* needs the output of Part 4.2 (RUN_CELLS == 1 at least once)
capture confirm file "$OUT/measure_density_cells.dta"
if _rc == 0 {

use "$OUT/measure_density_cells.dta", clear
merge 1:1 codunico year using "$OUT/measure_bindiff.dta", keep(match) nogenerate

gen byte post1998 = (year >= 1998)
label define post 0 "1994-1997" 1 "1998-2003"
label values post1998 post

tabstat logratio logratio_fix jump3 jump5, by(year) statistics(mean median n) ///
    format(%7.3f)

pwcorr logratio logratio_fix jump3 jump5, obs sig
spearman logratio_fix jump3
bysort post1998: pwcorr logratio_fix jump3, obs sig

twoway (scatter logratio_fix jump3, msymbol(oh) mcolor(navy%40)),        ///
    by(post1998, note(""))                                               ///
    xtitle("C&C measure: share 45-47 minus share 48-50")               ///
    ytitle("CJM log ratio (bandwidth 3)")
graph export "$OUT/fig3_measure_comparison.png", replace width(2000)

save "$OUT/measures_compared.dta", replace

}
else display as error "Part 4.3 skipped: run Part 4.2 first (set RUN_CELLS to 1)"


log close

*******************************************************************************

* PART 5:  EXPLANATION ARTIFICIAL CUTOFFS

*******************************************************************************

use "$OUT/sample.dta", clear
histogram score_cc if inrange(score_cc, 39, 51), width(0.05) frequency ///
    xline(42 48 , lcolor(green)) xtitle("Continuous score") ///
    title("Score distribution, all years")
graph export "$OUT/fig_fine_hist.png", replace width(2000)


*separate for each period

use "$OUT/sample.dta", clear
gen byte post = (year >= 1998)
label define per 0 "1994-1997" 1 "1998-2003"
label values post per

* percent (not frequency): the periods differ in size
histogram score_cc if inrange(score_cc, 39, 57), width(0.05) percent ///
    by(post, cols(1) note("")) xline(42 48 54, lcolor(green))           ///
    xtitle("Continuous score")
graph export "$OUT/fig_fine_hist_by_period.png", replace width(2000)





