*! version 0.1.0 19mar2026
program xt2denoise, eclass
version 18.0

syntax varname(numeric) [if], z(varname numeric) treatment(varname numeric) control(varname numeric) [, pre(integer 1) post(integer 3) baseline(string) cluster(varname) graph]

* read panel structure
xtset
local group = r(panelvar)
local time = r(timevar)
local y `varlist'

* set default values
if ("`baseline'" == "") {
    local baseline "-1"
}
if ("`cluster'" == "") {
    local cluster "`group'"
}
local K = `pre' + `post' + 1

marksample touse

tempvar evert everc time_g eventtime

* identify ever-treated and ever-control groups
quietly egen `evert' = max(cond(`touse', `treatment', 0)), by(`group')
quietly egen `everc' = max(cond(`touse', `control', 0)), by(`group')

* no group can receive both treatments
capture assert !(`evert' & `everc') if `touse'
if (_rc) {
    display as error "Some groups receive both treatments"
    inspect `group' if `evert' & `everc' & `touse'
    error 459
}

* compute treatment timing for each group
quietly egen `time_g' = min(cond(`treatment' | `control', `time', .)) if `touse', by(`group')

* everyone must receive a treatment
capture assert !missing(`time_g') if `touse'
if (_rc) {
    display as error "Some groups do not receive any treatment"
    inspect `group' if !`evert' & !`everc' & `touse'
    error 459
}

* compute event time relative to treatment
quietly generate `eventtime' = `time' - `time_g' if `touse'

* validate baseline option
if ("`baseline'" != "average" & "`baseline'" != "atet") {
    if (!inrange(`baseline', -`pre', -1)) {
        display as error "Baseline must be between -`pre' and -1, or 'average' or 'atet'"
        error 198
    }
}

***** STEP 1: Compute quality change dz = mean(z after) - mean(z before) for each group
tempvar z_before z_after dz

* mean of z before event (eventtime < 0)
quietly egen `z_before' = mean(cond(`eventtime' < 0, `z', .)) if `touse', by(`group')
* mean of z after event (eventtime >= 0)
quietly egen `z_after' = mean(cond(`eventtime' >= 0, `z', .)) if `touse', by(`group')
* quality change
quietly generate `dz' = `z_after' - `z_before' if `touse'

***** STEP 2: Compute yt - yg to remove group fixed effect
tempvar yg dy

* use baseline period (eventtime == -1) mean as group effect
quietly egen `yg' = mean(cond(`eventtime' == -1, `y', .)) if `touse', by(`group')
quietly generate `dy' = `y' - `yg' if `touse'

***** STEP 3: Remove event time X treatment group specific mean from dy
* Save number of treated and control units in each event time
tempvar n1e n0e dy_mean dy_demean

* count treated units by event time
quietly egen `n1e' = total(cond(`evert', 1, 0)) if `touse', by(`eventtime')
* count control units by event time
quietly egen `n0e' = total(cond(`everc', 1, 0)) if `touse', by(`eventtime')

* compute event time X treatment group specific mean
quietly egen `dy_mean' = mean(`dy') if `touse', by(`eventtime' `evert')
* demean dy
quietly generate `dy_demean' = `dy' - `dy_mean' if `touse'

***** STEP 4: For each event time, compute variance and covariance using regression
* var(dy_demean), cov(dy_demean, dz), var(dz)
* OLS of y^2 on constant gives variance estimate

quietly levelsof `eventtime' if `touse' & inrange(`eventtime', -`pre', `post'), local(eventtimes)

tempvar dy2 dz2 dydz
quietly generate `dy2' = `dy_demean'^2 if `touse'
quietly generate `dz2' = `dz'^2 if `touse'
quietly generate `dydz' = `dy_demean' * `dz' if `touse'

* store results
local K = `pre' + `post' + 1
tempname var_y cov_yz var_z se_var_y se_cov_yz se_var_z n1_vec n0_vec

matrix `var_y' = J(1, `K', .)
matrix `cov_yz' = J(1, `K', .)
matrix `var_z' = J(1, `K', .)
matrix `se_var_y' = J(1, `K', .)
matrix `se_cov_yz' = J(1, `K', .)
matrix `se_var_z' = J(1, `K', .)
matrix `n1_vec' = J(1, `K', .)
matrix `n0_vec' = J(1, `K', .)

foreach e of local eventtimes {
    local col = `e' + `pre' + 1

    * variance of dy_demean (regress dy^2 on constant)
    quietly regress `dy2' if `touse' & `eventtime' == `e', cluster(`cluster')
    matrix `var_y'[1, `col'] = _b[_cons]
    matrix `se_var_y'[1, `col'] = _se[_cons]

    * covariance of dy_demean with dz (regress dy*dz on constant)
    quietly regress `dydz' if `touse' & `eventtime' == `e', cluster(`cluster')
    matrix `cov_yz'[1, `col'] = _b[_cons]
    matrix `se_cov_yz'[1, `col'] = _se[_cons]

    * variance of dz (regress dz^2 on constant)
    quietly regress `dz2' if `touse' & `eventtime' == `e', cluster(`cluster')
    matrix `var_z'[1, `col'] = _b[_cons]
    matrix `se_var_z'[1, `col'] = _se[_cons]

    * store sample sizes
    quietly summarize `n1e' if `touse' & `eventtime' == `e', meanonly
    matrix `n1_vec'[1, `col'] = r(mean)
    quietly summarize `n0e' if `touse' & `eventtime' == `e', meanonly
    matrix `n0_vec'[1, `col'] = r(mean)
}

* label columns
local colnames ""
forvalues t = -`pre'/`post' {
    local colnames `colnames' `t'
}
matrix colname `var_y' = `colnames'
matrix colname `cov_yz' = `colnames'
matrix colname `var_z' = `colnames'
matrix colname `n1_vec' = `colnames'
matrix colname `n0_vec' = `colnames'

* display results
display
display as text "Second moment estimates by event time"
display as text "{hline 70}"
display as text "Event time" _col(15) "Var(dy)" _col(30) "Cov(dy,dz)" _col(45) "Var(dz)" _col(60) "N1" _col(67) "N0"
display as text "{hline 70}"

forvalues t = -`pre'/`post' {
    local col = `t' + `pre' + 1
    display as result %9.0f `t' _col(15) %9.4f `var_y'[1, `col'] _col(30) %9.4f `cov_yz'[1, `col'] _col(45) %9.4f `var_z'[1, `col'] _col(60) %5.0f `n1_vec'[1, `col'] _col(67) %5.0f `n0_vec'[1, `col']
}
display as text "{hline 70}"

* return results
ereturn clear
ereturn matrix var_y = `var_y'
ereturn matrix cov_yz = `cov_yz'
ereturn matrix var_z = `var_z'
ereturn matrix se_var_y = `se_var_y'
ereturn matrix se_cov_yz = `se_cov_yz'
ereturn matrix se_var_z = `se_var_z'
ereturn matrix n1 = `n1_vec'
ereturn matrix n0 = `n0_vec'
ereturn local cmd "xt2denoise"
ereturn local cmdline "xt2denoise `0'"
ereturn local depvar "`y'"

end
