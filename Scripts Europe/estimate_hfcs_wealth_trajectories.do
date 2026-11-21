/*
HFCS wealth trajectories

Cross-country adaptation of Wealth_Parent_Traj_2026-11-09.do and
Wealth_Child_Traj_2026-11-09.do. I use repeated HFCS cross-sections rather
than a panel. Quartiles refer to age-specific current wealth, not a fixed
household classification. Estimates are produced separately in the five
imputations and then averaged. Fan bands are household dispersion, not CIs.

Outputs are written to HFCS/wealth_trajectories/{country}/{role}_{sample}/.
Run from the repository root, or change repo below.

Requires Stata 16+ with `pq` installed (ssc install pq), or a current
StataNow supporting `import parquet`. The dataset must fit in Stata memory.
*/

version 16.0
clear all
set more off
set varabbrev off

* Keep figure windows closed while writing PDF and PNG exports.
set graphics on

* Settings
local repo "`c(pwd)'"
local data "`repo'/HFCS/analysis_ready/hfcs_analysis.parquet"
capture confirm file `"`data'"'
if _rc {
    local repo "`c(pwd)'/.."
    local data "`repo'/HFCS/analysis_ready/hfcs_analysis.parquet"
}
capture confirm file `"`data'"'
if _rc {
    di as error "Cannot find hfcs_analysis.parquet. Start Stata in the repository or set local repo."
    exit 601
}
* Use short paths relative to the repository. Windows can reject long Excel paths.
cd `"`repo'"'
local data "HFCS/analysis_ready/hfcs_analysis.parquet"
local out "HFCS/wealth_trajectories"
capture mkdir `"`out'"'

local run_comparisons 1   // 0: Parent/all and Child/no_inheritance only
local min_age_n 20
local min_cell_n 8
local min_resid_n 8
local min_lowess_ages 6
local min_imp 3
local tail_years 10
local min_tail_ages 4
local make_png 1

local need "country country_raw wave imputation age_years age_lower net_wealth_real hh_weight eligible_all eligible_no_inheritance wealth_quartile_all wealth_quartile_no_inheritance"

* StataNow 19 provides import parquet. The pq package works with Stata 16+.
capture noisily import parquet `need' using `"`data'"', clear
if _rc {
    capture which pq
    if _rc {
        di as error "Parquet importer unavailable. Run: ssc install pq"
        exit 499
    }
    pq use `need' using `"`data'"', clear
}
keep if country_raw != "E1"
assert inrange(imputation, 1, 5)
assert inlist(wave, 2010, 2014, 2017, 2021, 2023)
compress

tempfile master allcurves audit_data binaudit_data
save `master', replace
levelsof country, local(countries)
* To limit the run: local countries "AT DE FR"

* One country, one specification, one imputation at a time.
capture program drop _hfcs_single_imputation
program define _hfcs_single_imputation
    version 16.0
    syntax, MINRAW(integer) MINCELL(integer) MINRESID(integer) MINAGES(integer)
    keep age_years net_wealth_real hh_weight wealth_quar
    rename age_years age
    rename net_wealth_real wealth
    drop if missing(age, wealth, hh_weight, wealth_quar) | hh_weight <= 0
    keep if inrange(wealth_quar, 1, 4)
    if _N == 0 exit
    tempfile individuals cells

    * Match the original unweighted p1/p99 trimming within exact age.
    bysort age: egen long n_age = count(wealth)
    bysort age: egen double p1_age = pctile(wealth), p(1)
    bysort age: egen double p99_age = pctile(wealth), p(99)
    gen double wealth_graph = wealth if n_age >= `minraw' & ///
        inrange(wealth, p1_age, p99_age)
    save `individuals', replace

    keep if !missing(wealth_graph)
    if _N == 0 {
        clear
        exit
    }
    gen double weighted_wealth = wealth_graph * hh_weight
    collapse (sum) weighted_wealth hh_weight (count) n = wealth_graph, ///
        by(age wealth_quar)
    gen double wealth_graph = weighted_wealth / hh_weight if ///
        n >= `mincell' & hh_weight > 0
    keep age wealth_quar wealth_graph n
    gen double lp_y = .
    forvalues q = 1/4 {
        quietly count if wealth_quar == `q' & !missing(wealth_graph)
        if r(N) >= `minages' {
            quietly lowess wealth_graph age if wealth_quar == `q' & ///
                !missing(wealth_graph), gen(smooth_q) nograph bwidth(.8)
            replace lp_y = smooth_q if wealth_quar == `q' & !missing(smooth_q)
            drop smooth_q
        }
    }
    save `cells', replace

    use `individuals', clear
    rename wealth_graph wealth_trimmed
    merge m:1 age wealth_quar using `cells', ///
        keep(master match) keepusing(wealth_graph lp_y n) nogen
    keep if !missing(lp_y)
    if _N == 0 {
        clear
        exit
    }
    gen double resid = wealth_trimmed - lp_y
    bysort wealth_quar: egen long n_resid_q = count(resid)
    bysort wealth_quar: egen double resid_p5_q = pctile(resid), p(5)
    bysort wealth_quar: egen double resid_p95_q = pctile(resid), p(95)
    replace resid = . if n_resid_q < 20 | ///
        (!missing(resid) & !inrange(resid, resid_p5_q, resid_p95_q))
    bysort age wealth_quar: egen long n_resid_age = count(resid)

    foreach p in 5 10 20 30 40 {
        local other = 100 - `p'
        bysort age wealth_quar: egen double rlo`p' = pctile(resid), p(`p')
        bysort age wealth_quar: egen double rhi`p' = pctile(resid), p(`other')
        gen double lo`p' = lp_y + rlo`p' if n_resid_age >= `minresid'
        gen double hi`p' = lp_y + rhi`p' if n_resid_age >= `minresid'
    }
    collapse (mean) wealth_graph lp_y lo5 hi5 lo10 hi10 ///
        lo20 hi20 lo30 hi30 lo40 hi40 n, by(age wealth_quar)
    keep if !missing(lp_y)
end

* Each band and the LOWESS estimate are extrapolated separately.
capture program drop _hfcs_extrapolate
program define _hfcs_extrapolate
    version 16.0
    syntax, TAIL(integer) MINPOINTS(integer)
    tempfile projected
    tempname ph
    local series "lp_y lo5 hi5 lo10 hi10 lo20 hi20 lo30 hi30 lo40 hi40"
    postfile `ph' byte wealth_quar int age byte extrapolated ///
        double wealth_graph n_avg n_imp `series' using `projected', replace

    forvalues q = 1/4 {
        quietly summarize age if wealth_quar == `q' & !missing(lp_y), meanonly
        if r(N) == 0 continue
        local last = r(max)
        if `last' >= 100 continue

        foreach v of local series {
            local slope_`v' = .
            local anchor_`v' = .
            local anchorage_`v' = .
            quietly count if wealth_quar == `q' & ///
                inrange(age, `last' - `tail', `last') & !missing(`v')
            if r(N) >= `minpoints' {
                capture quietly regress `v' age if wealth_quar == `q' & ///
                    inrange(age, `last' - `tail', `last') & !missing(`v')
                if !_rc {
                    local slope_`v' = _b[age]
                    quietly summarize age if wealth_quar == `q' & !missing(`v'), meanonly
                    local anchorage_`v' = r(max)
                    quietly summarize `v' if wealth_quar == `q' & ///
                        age == `anchorage_`v'', meanonly
                    local anchor_`v' = r(mean)
                }
            }
        }
        local first = `last' + 1
        forvalues a = `first'/100 {
            local values ""
            foreach v of local series {
                local prediction = .
                if !missing(`slope_`v'', `anchor_`v'', `anchorage_`v'') {
                    local prediction = `anchor_`v'' + (`a' - `anchorage_`v'') * `slope_`v''
                }
                local values `"`values' (`prediction')"'
            }
            post `ph' (`q') (`a') (1) (.) (.) (.) `values'
        }
    }
    postclose `ph'
    gen byte extrapolated = 0
    append using `projected'
    sort wealth_quar age extrapolated
end

* The binned export also covers countries where exact age is suppressed.
capture program drop _hfcs_binned
program define _hfcs_binned
    version 16.0
    syntax, MINCELL(integer) MINIMP(integer)
    rename age_lower age_bin
    rename net_wealth_real wealth
    keep if !missing(age_bin, wealth, hh_weight, wealth_quar) & hh_weight > 0
    if _N == 0 exit
    gen double weighted_wealth = wealth * hh_weight
    collapse (sum) weighted_wealth hh_weight (count) n = wealth, ///
        by(imputation wealth_quar age_bin)
    gen double yhat_xbd = weighted_wealth / hh_weight if n >= `mincell'
    keep if !missing(yhat_xbd)
    if _N == 0 exit
    collapse (mean) yhat_xbd n (count) n_imp = yhat_xbd, by(wealth_quar age_bin)
    keep if n_imp >= `minimp'
    sort wealth_quar age_bin
end

* Build an audit in addition to the original-format Excel files.
tempname audit bin_audit
postfile `audit' str2 country str6 role str18 sample str12 status ///
    long n_sample int n_exact int n_quartiles int observed_points int projected_points ///
    using `audit_data', replace
postfile `bin_audit' str2 country str6 role str18 sample ///
    str12 status int observed_bins projected_bins using `binaudit_data', replace
local first_curves = 1

foreach c of local countries {
    capture mkdir `"`out'/`c'"'
    foreach role in Parent Child {
        foreach sample in all no_inheritance {
            if `run_comparisons' == 0 & ///
                (("`role'" == "Parent" & "`sample'" == "no_inheritance") | ///
                 ("`role'" == "Child" & "`sample'" == "all")) continue

            local folder `"`out'/`c'/`role'_`sample'"'
            capture mkdir `"`folder'"'
            capture mkdir `"`folder'/Coefficients"'
            capture mkdir `"`folder'/Figures"'
            if "`sample'" == "all" {
                local elig "eligible_all"
                local quart "wealth_quartile_all"
            }
            else {
                local elig "eligible_no_inheritance"
                local quart "wealth_quartile_no_inheritance"
            }
            if "`role'" == "Parent" {
                local agemin 50
                local agemax 100
            }
            else {
                local agemin 20
                local agemax 65
            }

            * Exact-age sample and one-imputation LOWESS fits.
            use `master', clear
            keep if country == "`c'" & `elig' == 1 & ///
                inrange(age_years, `agemin', `agemax') & ///
                !missing(net_wealth_real, hh_weight, `quart') & ///
                hh_weight > 0 & inrange(`quart', 1, 4)
            gen byte wealth_quar = `quart'
            quietly count if imputation == 1
            local n_sample = r(N)
            quietly levelsof age_years, local(exactages)
            local n_exact : word count `exactages'
            tempfile exactsample estimates
            save `exactsample', replace
            local impnums ""
            forvalues m = 1/5 {
                use `exactsample', clear
                keep if imputation == `m'
                if _N == 0 continue
                quietly _hfcs_single_imputation, minraw(`min_age_n') ///
                    mincell(`min_cell_n') minresid(`min_resid_n') ///
                    minages(`min_lowess_ages')
                if _N == 0 continue
                gen byte imputation = `m'
                tempfile i`m'
                save `i`m'', replace
                local impnums "`impnums' `m'"
            }

            local usable 0
            if "`impnums'" != "" {
                local first = 1
                foreach m of local impnums {
                    if `first' {
                        use `i`m'', clear
                        local first = 0
                    }
                    else append using `i`m''
                }
                foreach band in lo5 hi5 lo10 hi10 lo20 hi20 lo30 hi30 lo40 hi40 {
                    bysort age wealth_quar: egen byte count_`band' = count(`band')
                }
                collapse (mean) wealth_graph lp_y lo5 hi5 lo10 hi10 ///
                    lo20 hi20 lo30 hi30 lo40 hi40 (mean) n_avg=n ///
                    (count) n_imp=lp_y ///
                    (max) count_lo5 count_hi5 count_lo10 count_hi10 ///
                    count_lo20 count_hi20 count_lo30 count_hi30 count_lo40 count_hi40, ///
                    by(age wealth_quar)
                keep if n_imp >= `min_imp' & !missing(lp_y)
                foreach band in lo5 hi5 lo10 hi10 lo20 hi20 lo30 hi30 lo40 hi40 {
                    replace `band' = . if count_`band' < `min_imp'
                    drop count_`band'
                }
                quietly levelsof wealth_quar, local(qavailable)
                local n_q : word count `qavailable'
                local observed_n = _N
                if `n_q' == 4 & `observed_n' > 0 local usable 1
            }
            if !`usable' {
                post `audit' ("`c'") ("`role'") ("`sample'") ("skipped") ///
                    (`n_sample') (`n_exact') (0) (0) (0)
                di as text "SKIP `c' / `role' / `sample': insufficient exact-age data"
            }
            else {
                tempfile observed extended
                save `observed', replace
                quietly _hfcs_extrapolate, tail(`tail_years') ///
                    minpoints(`min_tail_ages')
                save `extended', replace
                quietly count if extrapolated == 1
                local projected_n = r(N)
                post `audit' ("`c'") ("`role'") ("`sample'") ("ok") ///
                    (`n_sample') (`n_exact') (`n_q') (`observed_n') (`projected_n')

                * ABM exports: original Stata column names and Sheet1 layout.
                gen str13 quartile_label = ""
                replace quartile_label = "1st Quartile" if wealth_quar == 1
                replace quartile_label = "2nd Quartile" if wealth_quar == 2
                replace quartile_label = "3rd Quartile" if wealth_quar == 3
                replace quartile_label = "4th Quartile" if wealth_quar == 4
                sort age wealth_quar extrapolated
                if "`role'" == "Parent" {
                    export excel quartile_label wealth_quar age lp_y ///
                        lo5 hi5 lo10 hi10 lo20 hi20 lo30 hi30 lo40 hi40 ///
                        extrapolated using ///
                        `"`folder'/Parent_wealth_byage_lowess_fan_extrapolated100_2026_26_05.xlsx"', ///
                        sheet("Sheet1") firstrow(variables) replace
                    preserve
                    rename lp_y lowess_y
                    export excel wealth_quar age lowess_y extrapolated using ///
                        `"`folder'/Coefficients/Parent_wealth_lowess_trimmed_extrapolated100_last5_2026_20_05.xlsx"', ///
                        sheet("Sheet1") firstrow(variables) replace
                    restore
                }
                else {
                    preserve
                    keep if extrapolated == 0
                    export excel quartile_label wealth_quar age lp_y ///
                        lo5 hi5 lo10 hi10 lo20 hi20 lo30 hi30 lo40 hi40 using ///
                        `"`folder'/Child_wealth_byage_lowess_fan_2026_26_05.xlsx"', ///
                        sheet("Sheet1") firstrow(variables) replace
                    restore
                    export excel quartile_label wealth_quar age lp_y ///
                        lo5 hi5 lo10 hi10 lo20 hi20 lo30 hi30 lo40 hi40 ///
                        extrapolated using ///
                        `"`folder'/Child_wealth_byage_lowess_fan_extrapolated100_HFCS.xlsx"', ///
                        sheet("Sheet1") firstrow(variables) replace
                    preserve
                    rename lp_y lowess_y
                    export excel age wealth_quar wealth_graph lowess_y extrapolated using ///
                        `"`folder'/Coefficients/Child_wealth_lowess_trimmed_extrapolated100_last10_2026_08_05.xlsx"', ///
                        sheet("Sheet1") firstrow(variables) replace
                    restore
                }
                gen str2 country = "`c'"
                gen str6 role = "`role'"
                gen str18 sample = "`sample'"
                export delimited using `"`folder'/trajectory_numeric.csv"', replace
                preserve
                keep if extrapolated == 0
                export delimited using `"`folder'/trajectory_observed.csv"', replace
                restore

                * Raw, LOWESS and fan figures, in the original four-panel layout.
                forvalues q = 1/4 {
                    local col navy
                    if `q' == 2 local col cranberry
                    if `q' == 3 local col forest_green
                    if `q' == 4 local col orange_red
                    local title "Low: 1st Quartile"
                    if `q' == 2 local title "Low-Middle: 2nd Quartile"
                    if `q' == 3 local title "High-Middle: 3rd Quartile"
                    if `q' == 4 local title "High: 4th Quartile"
                    twoway (line wealth_graph age if wealth_quar == `q' & extrapolated == 0, ///
                        sort lcolor(`col') lwidth(medthick)), ///
                        title("`title'") legend(off) xtitle("Age") ///
                        ytitle("Net wealth (2023 EUR)") name(raw`q', replace)
                    twoway ///
                        (line lp_y age if wealth_quar == `q' & extrapolated == 0, ///
                            sort lcolor(`col') lwidth(medthick)) ///
                        (line lp_y age if wealth_quar == `q' & extrapolated == 1, ///
                            sort lcolor(`col') lpattern(dash) lwidth(medthick)), ///
                        title("`title'") legend(off) xtitle("Age") ///
                        ytitle("Net wealth (2023 EUR)") name(line`q', replace)
                    twoway ///
                        (rarea lo5 hi5 age if wealth_quar == `q' & extrapolated == 0, sort color(`col'%12) lwidth(none)) ///
                        (rarea lo10 hi10 age if wealth_quar == `q' & extrapolated == 0, sort color(`col'%15) lwidth(none)) ///
                        (rarea lo20 hi20 age if wealth_quar == `q' & extrapolated == 0, sort color(`col'%20) lwidth(none)) ///
                        (rarea lo30 hi30 age if wealth_quar == `q' & extrapolated == 0, sort color(`col'%25) lwidth(none)) ///
                        (rarea lo40 hi40 age if wealth_quar == `q' & extrapolated == 0, sort color(`col'%30) lwidth(none)) ///
                        (rarea lo5 hi5 age if wealth_quar == `q' & extrapolated == 1, sort color(`col'%6) lwidth(none)) ///
                        (rarea lo10 hi10 age if wealth_quar == `q' & extrapolated == 1, sort color(`col'%8) lwidth(none)) ///
                        (rarea lo20 hi20 age if wealth_quar == `q' & extrapolated == 1, sort color(`col'%10) lwidth(none)) ///
                        (rarea lo30 hi30 age if wealth_quar == `q' & extrapolated == 1, sort color(`col'%12) lwidth(none)) ///
                        (rarea lo40 hi40 age if wealth_quar == `q' & extrapolated == 1, sort color(`col'%15) lwidth(none)) ///
                        (line lp_y age if wealth_quar == `q' & extrapolated == 0, sort lcolor(`col') lwidth(medthick)) ///
                        (line lp_y age if wealth_quar == `q' & extrapolated == 1, sort lcolor(`col') lwidth(medthick) lpattern(dash)), ///
                        title("`title'") legend(off) xtitle("Age") ///
                        ytitle("Net wealth (2023 EUR)") name(fan`q', replace)
                }
                graph combine raw1 raw2 raw3 raw4, rows(2) cols(2) ///
                    name(raw_all, replace)
                graph export `"`folder'/Figures/`role'_wealth_byage_trimmed_check_HFCS.pdf"', replace
                if `make_png' graph export `"`folder'/Figures/`role'_wealth_byage_trimmed_check_HFCS.png"', replace
                graph combine line1 line2 line3 line4, rows(2) cols(2) ///
                    name(line_all, replace)
                graph export `"`folder'/Figures/`role'_wealth_byage_lowess_trimmed_extrapolated100_last10_HFCS.pdf"', replace
                if `make_png' graph export `"`folder'/Figures/`role'_wealth_byage_lowess_trimmed_extrapolated100_last10_HFCS.png"', replace
                graph combine fan1 fan2 fan3 fan4, rows(2) cols(2) ///
                    name(fan_all, replace)
                graph export `"`folder'/Figures/`role'_wealth_byage_lowess_fan_extrapolated100_HFCS.pdf"', replace
                if `make_png' graph export `"`folder'/Figures/`role'_wealth_byage_lowess_fan_extrapolated100_HFCS.png"', replace
                if "`role'" == "Child" {
                    preserve
                    keep if extrapolated == 0
                    forvalues q = 1/4 {
                        local col navy
                        if `q' == 2 local col cranberry
                        if `q' == 3 local col forest_green
                        if `q' == 4 local col orange_red
                        twoway ///
                            (rarea lo5 hi5 age if wealth_quar == `q', sort color(`col'%12) lwidth(none)) ///
                            (rarea lo10 hi10 age if wealth_quar == `q', sort color(`col'%15) lwidth(none)) ///
                            (rarea lo20 hi20 age if wealth_quar == `q', sort color(`col'%20) lwidth(none)) ///
                            (rarea lo30 hi30 age if wealth_quar == `q', sort color(`col'%25) lwidth(none)) ///
                            (rarea lo40 hi40 age if wealth_quar == `q', sort color(`col'%30) lwidth(none)) ///
                            (line lp_y age if wealth_quar == `q', sort lcolor(`col') lwidth(medthick)), ///
                            title("Quartile `q'") legend(off) xtitle("Age") ///
                            ytitle("Net wealth (2023 EUR)") name(obsfan`q', replace)
                    }
                    graph combine obsfan1 obsfan2 obsfan3 obsfan4, rows(2) cols(2) ///
                        name(obsfan_all, replace)
                    graph export `"`folder'/Figures/Child_wealth_byage_lowess_fan_observed_HFCS.pdf"', replace
                    if `make_png' graph export `"`folder'/Figures/Child_wealth_byage_lowess_fan_observed_HFCS.png"', replace
                    restore
                }
                * Collect the curves across countries for the ABM.
                if `first_curves' {
                    save `allcurves', replace
                    local first_curves = 0
                }
                else {
                    append using `allcurves'
                    save `allcurves', replace
                }
                di as text "OK `c' / `role' / `sample': `observed_n' observed cells"
            }

            * Five-year estimates also use genuine bracket-only ages.
            use `master', clear
            keep if country == "`c'" & `elig' == 1 & ///
                inrange(age_lower, `agemin', 84) & ///
                !missing(net_wealth_real, hh_weight, `quart') & ///
                hh_weight > 0 & inrange(`quart', 1, 4)
            gen byte wealth_quar = `quart'
            quietly _hfcs_binned, mincell(`min_cell_n') minimp(`min_imp')
            if _N == 0 {
                post `bin_audit' ("`c'") ("`role'") ("`sample'") ("skipped") (0) (0)
                continue
            }
            local n_bins = _N
            export excel wealth_quar age_bin yhat_xbd n using ///
                `"`folder'/Coefficients/`role'_wealth_byage_binned_trimmed_HFCS.xlsx"', ///
                sheet("Sheet1") firstrow(variables) replace
            export delimited using `"`folder'/age_bin_5y.csv"', replace

            tempfile binproject
            tempname bh
            postfile `bh' byte wealth_quar int age_bin byte extrapolated ///
                double yhat_xbd n n_imp using `binproject', replace
            forvalues q = 1/4 {
                quietly summarize age_bin if wealth_quar == `q' & !missing(yhat_xbd), meanonly
                if r(N) < 3 continue
                local last = r(max)
                quietly count if wealth_quar == `q' & ///
                    inrange(age_bin, `last' - `tail_years', `last') & !missing(yhat_xbd)
                if r(N) < 3 continue
                quietly regress yhat_xbd age_bin if wealth_quar == `q' & ///
                    inrange(age_bin, `last' - `tail_years', `last')
                local slope = _b[age_bin]
                quietly summarize yhat_xbd if wealth_quar == `q' & age_bin == `last', meanonly
                local anchor = r(mean)
                local firstbin = `last' + 5
                forvalues a = `firstbin'(5)100 {
                    local prediction = `anchor' + `slope' * (`a' - `last')
                    post `bh' (`q') (`a') (1) (`prediction') (.) (.)
                }
            }
            postclose `bh'
            gen byte extrapolated = 0
            append using `binproject'
            sort wealth_quar age_bin extrapolated
            quietly count if extrapolated == 1
            local n_projected_bins = r(N)
            if "`role'" == "Parent" {
                export excel wealth_quar age_bin yhat_xbd extrapolated using ///
                    `"`folder'/Parent_wealth_binned_trimmed_extrapolated100_HFCS.xlsx"', ///
                    sheet("Sheet1") firstrow(variables) replace
            }
            else {
                export excel wealth_quar age_bin yhat_xbd n extrapolated using ///
                    `"`folder'/Child_wealth_binned_trimmed_extrapolated100_HFCS.xlsx"', ///
                    sheet("Sheet1") firstrow(variables) replace
            }
            export delimited using `"`folder'/age_bin_5y_extrapolated.csv"', replace
            quietly levelsof wealth_quar, local(binquartiles)
            local n_bin_q : word count `binquartiles'
            if `n_bin_q' == 4 {
                forvalues q = 1/4 {
                    local col navy
                    if `q' == 2 local col cranberry
                    if `q' == 3 local col forest_green
                    if `q' == 4 local col orange_red
                    twoway ///
                        (line yhat_xbd age_bin if wealth_quar == `q' & extrapolated == 0, sort lcolor(`col') lwidth(medthick)) ///
                        (line yhat_xbd age_bin if wealth_quar == `q' & extrapolated == 1, sort lcolor(`col') lpattern(dash) lwidth(medthick)), ///
                        title("Quartile `q'") legend(off) xtitle("Age bin") ///
                        ytitle("Net wealth (2023 EUR)") name(bin`q', replace)
                }
                graph combine bin1 bin2 bin3 bin4, rows(2) cols(2) name(bin_all, replace)
                graph export `"`folder'/Figures/`role'_wealth_byage_binned_extrapolated100_HFCS.pdf"', replace
                if `make_png' graph export `"`folder'/Figures/`role'_wealth_byage_binned_extrapolated100_HFCS.png"', replace
            }
            post `bin_audit' ("`c'") ("`role'") ("`sample'") ///
                ("ok") (`n_bins') (`n_projected_bins')
        }
    }
}
postclose `audit'
postclose `bin_audit'

use `audit_data', clear
export delimited using `"`out'/trajectory_sample_audit.csv"', replace
use `binaudit_data', clear
export delimited using `"`out'/binned_sample_audit.csv"', replace
if !`first_curves' {
    use `allcurves', clear
    sort country role sample age wealth_quar extrapolated
    export delimited using `"`out'/ABM_wealth_trajectories_all_countries.csv"', replace
    export delimited country role sample quartile_label wealth_quar age lp_y ///
        lo5 hi5 lo10 hi10 lo20 hi20 lo30 hi30 lo40 hi40 extrapolated ///
        using `"`out'/ABM_wealth_fans_all_countries.csv"', replace
    gen byte crossing = 0
    foreach p in 5 10 20 30 40 {
        replace crossing = 1 if lo`p' > hi`p' & !missing(lo`p', hi`p')
    }
    preserve
    keep if crossing == 1
    export delimited using `"`out'/fan_band_crossings_review.csv"', replace
    restore
}
set graphics on
di as result "Done. Outputs: `out'"
