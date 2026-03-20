*! version 0.3.0 20mar2026
program xt2denoise, eclass
version 18.0

syntax varname(numeric) [if], z(varname numeric) treatment(varname numeric) control(varname numeric) [, pre(integer 1) post(integer 3) baseline(string) cluster(varname) graph detail COVariance]

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

quietly egen `z_before' = mean(cond(`eventtime' < 0, `z', .)) if `touse', by(`group')
quietly egen `z_after' = mean(cond(`eventtime' >= 0, `z', .)) if `touse', by(`group')
quietly generate `dz' = `z_after' - `z_before' if `touse'

***** STEP 2: Compute yt - yg to remove group fixed effect
tempvar yg dy

quietly egen `yg' = mean(cond(`eventtime' == -1, `y', .)) if `touse', by(`group')
quietly generate `dy' = `y' - `yg' if `touse'

***** STEP 3: Remove event time X treatment group specific mean from dy
tempvar dy_mean dy_demean

quietly egen `dy_mean' = mean(`dy') if `touse', by(`eventtime' `evert')
quietly generate `dy_demean' = `dy' - `dy_mean' if `touse'

***** STEP 4: Compute variance and covariance separately for treated and control
tempvar dy2 dz2 dydz eventtime100 postvar

quietly generate `dy2' = `dy_demean'^2 if `touse'
quietly generate `dz2' = `dz'^2 if `touse'
quietly generate `dydz' = `dy_demean' * `dz' if `touse'
quietly generate `eventtime100' = `eventtime' + 100 if `touse'
quietly generate `postvar' = (`eventtime' >= 0) if `touse'

* determine grouping variable and column names based on baseline
if ("`baseline'" == "atet") {
    local groupvar "`postvar'"
    local K = 1
    local colnames "ATET"
    local Kreg = 2
}
else {
    local groupvar "`eventtime100'"
    local K = `pre' + `post' + 1
    local Kreg = `K'
    local colnames ""
    forvalues t = -`pre'/`post' {
        local colnames `colnames' `t'
    }
}

* run regressions and extract coefficients using helper program
tempname cov1 cov0 var_z1 var_z0 se_cov1 se_cov0 se_var_z1 se_var_z0

_xt2denoise_regcoef `dydz' `groupvar' if `touse' & `evert' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') k(`Kreg') pre(`pre')
matrix `cov1' = r(coef)
matrix `se_cov1' = r(se)

_xt2denoise_regcoef `dydz' `groupvar' if `touse' & `everc' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') k(`Kreg') pre(`pre')
matrix `cov0' = r(coef)
matrix `se_cov0' = r(se)

_xt2denoise_regcoef `dz2' `groupvar' if `touse' & `evert' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') k(`Kreg') pre(`pre')
matrix `var_z1' = r(coef)
matrix `se_var_z1' = r(se)

_xt2denoise_regcoef `dz2' `groupvar' if `touse' & `everc' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') k(`Kreg') pre(`pre')
matrix `var_z0' = r(coef)
matrix `se_var_z0' = r(se)

* compute true variance, beta, and standard errors
tempname true_var_z beta se_beta beta_naive se_beta_naive n1_vec n0_vec

matrix `true_var_z' = `var_z1' - `var_z0'

_xt2denoise_compute_beta `cov1' `cov0' `true_var_z' `se_cov1' `se_cov0'
matrix `beta' = r(beta)
matrix `se_beta' = r(se)

* naive beta = cov1 / var_z1 (no differencing)
matrix `beta_naive' = J(1, `Kreg', .)
matrix `se_beta_naive' = J(1, `Kreg', .)
forvalues i = 1/`Kreg' {
    if `var_z1'[1, `i'] != 0 & `var_z1'[1, `i'] != . {
        matrix `beta_naive'[1, `i'] = `cov1'[1, `i'] / `var_z1'[1, `i']
        matrix `se_beta_naive'[1, `i'] = `se_cov1'[1, `i'] / abs(`var_z1'[1, `i'])
    }
}

* for ATET: compute difference (post - pre)
if ("`baseline'" == "atet") {
    tempname beta_atet se_atet beta_naive_atet se_naive_atet
    scalar `beta_atet' = `beta'[1, 2] - `beta'[1, 1]
    scalar `se_atet' = sqrt(`se_beta'[1, 1]^2 + `se_beta'[1, 2]^2)
    scalar `beta_naive_atet' = `beta_naive'[1, 2] - `beta_naive'[1, 1]
    scalar `se_naive_atet' = sqrt(`se_beta_naive'[1, 1]^2 + `se_beta_naive'[1, 2]^2)

    matrix `beta' = J(1, 1, `beta_atet')
    matrix `se_beta' = J(1, 1, `se_atet')
    matrix `beta_naive' = J(1, 1, `beta_naive_atet')
    matrix `se_beta_naive' = J(1, 1, `se_naive_atet')
}

* count sample sizes
matrix `n1_vec' = J(1, `K', .)
matrix `n0_vec' = J(1, `K', .)
if ("`baseline'" == "atet") {
    quietly count if `touse' & `evert' & inrange(`eventtime', -`pre', `post')
    matrix `n1_vec'[1, 1] = r(N)
    quietly count if `touse' & `everc' & inrange(`eventtime', -`pre', `post')
    matrix `n0_vec'[1, 1] = r(N)
}
else {
    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        quietly count if `touse' & `eventtime' == `t' & `evert'
        matrix `n1_vec'[1, `col'] = r(N)
        quietly count if `touse' & `eventtime' == `t' & `everc'
        matrix `n0_vec'[1, `col'] = r(N)
    }
}

* build variance matrices (diagonal)
tempname V V_naive
_xt2denoise_build_vmat `se_beta'
matrix `V' = r(V)

_xt2denoise_build_vmat `se_beta_naive'
matrix `V_naive' = r(V)

* apply column names
matrix colname `beta' = `colnames'
matrix colname `beta_naive' = `colnames'
matrix colname `V' = `colnames'
matrix rowname `V' = `colnames'
matrix colname `V_naive' = `colnames'
matrix rowname `V_naive' = `colnames'
matrix colname `cov1' = `colnames'
matrix colname `cov0' = `colnames'
matrix colname `se_cov1' = `colnames'
matrix colname `se_cov0' = `colnames'
matrix colname `var_z1' = `colnames'
matrix colname `var_z0' = `colnames'
matrix colname `true_var_z' = `colnames'
matrix colname `n1_vec' = `colnames'
matrix colname `n0_vec' = `colnames'

* count observations and post results
quietly count if `touse' & inrange(`eventtime', -`pre', `post')
local Nobs = r(N)

tempvar esample
quietly generate `esample' = `touse' & inrange(`eventtime', -`pre', `post')

ereturn post `beta' `V', obs(`Nobs') esample(`esample')
ereturn local depvar "`y'"
ereturn local cmd "xt2denoise"
ereturn local cmdline "xt2denoise `0'"

* store additional matrices
ereturn matrix cov1 = `cov1'
ereturn matrix cov0 = `cov0'
ereturn matrix se_cov1 = `se_cov1'
ereturn matrix se_cov0 = `se_cov0'
ereturn matrix var_z1 = `var_z1'
ereturn matrix var_z0 = `var_z0'
ereturn matrix true_var_z = `true_var_z'
ereturn matrix n1 = `n1_vec'
ereturn matrix n0 = `n0_vec'
ereturn matrix b_naive = `beta_naive'
ereturn matrix V_naive = `V_naive'

* compute covariance matrices for display
tempname cov_debiased V_cov_debiased V_cov_naive
matrix `cov_debiased' = e(cov1) - e(cov0)

local Kcov = colsof(e(cov1))
matrix `V_cov_debiased' = J(`Kcov', `Kcov', 0)
matrix `V_cov_naive' = J(`Kcov', `Kcov', 0)
forvalues i = 1/`Kcov' {
    matrix `V_cov_debiased'[`i', `i'] = e(se_cov1)[1, `i']^2 + e(se_cov0)[1, `i']^2
    matrix `V_cov_naive'[`i', `i'] = e(se_cov1)[1, `i']^2
}
matrix colname `cov_debiased' = `colnames'
matrix colname `V_cov_debiased' = `colnames'
matrix rowname `V_cov_debiased' = `colnames'
matrix colname `V_cov_naive' = `colnames'
matrix rowname `V_cov_naive' = `colnames'

ereturn matrix cov_debiased = `cov_debiased'
ereturn matrix V_cov_debiased = `V_cov_debiased'
ereturn matrix V_cov_naive = `V_cov_naive'

* display results
_xt2denoise_display, baseline(`baseline') depvar(`y') `detail' `covariance'

end

***** Helper program: run regression and extract coefficients
program _xt2denoise_regcoef, rclass
    syntax varlist(min=2 max=2) [if], cluster(varname) k(integer) pre(integer)

    gettoken depvar groupvar : varlist

    quietly regress `depvar' ibn.`groupvar' `if', cluster(`cluster') nocons

    tempname coef se
    matrix `coef' = J(1, `k', .)
    matrix `se' = J(1, `k', .)

    * check if this is ATET (k=2 with postvar) or event study
    if (`k' == 2) {
        * ATET: groupvar is postvar (0/1)
        capture matrix `coef'[1, 1] = _b[0.`groupvar']
        capture matrix `coef'[1, 2] = _b[1.`groupvar']
        capture matrix `se'[1, 1] = _se[0.`groupvar']
        capture matrix `se'[1, 2] = _se[1.`groupvar']
    }
    else {
        * event study: groupvar is eventtime100
        forvalues i = 1/`k' {
            local t = `i' - `pre' - 1
            local e100 = `t' + 100
            capture matrix `coef'[1, `i'] = _b[`e100'.`groupvar']
            capture matrix `se'[1, `i'] = _se[`e100'.`groupvar']
        }
    }

    return matrix coef = `coef'
    return matrix se = `se'
end

***** Helper program: compute beta = (cov1 - cov0) / true_var
program _xt2denoise_compute_beta, rclass
    args cov1 cov0 true_var se_cov1 se_cov0

    local K = colsof(`cov1')
    tempname beta se
    matrix `beta' = J(1, `K', .)
    matrix `se' = J(1, `K', .)

    forvalues i = 1/`K' {
        if `true_var'[1, `i'] != 0 & `true_var'[1, `i'] != . {
            matrix `beta'[1, `i'] = (`cov1'[1, `i'] - `cov0'[1, `i']) / `true_var'[1, `i']
            matrix `se'[1, `i'] = sqrt(`se_cov1'[1, `i']^2 + `se_cov0'[1, `i']^2) / abs(`true_var'[1, `i'])
        }
    }

    return matrix beta = `beta'
    return matrix se = `se'
end

***** Helper program: build diagonal variance matrix from SE vector
program _xt2denoise_build_vmat, rclass
    args se_vec

    local K = colsof(`se_vec')
    tempname V
    matrix `V' = J(`K', `K', 0)

    forvalues i = 1/`K' {
        if `se_vec'[1, `i'] != . {
            matrix `V'[`i', `i'] = `se_vec'[1, `i']^2
        }
    }

    return matrix V = `V'
end

***** Helper program: display regression tables
program _xt2denoise_display
    syntax, baseline(string) depvar(string) [detail COVariance]

    if ("`detail'" == "detail") {
        if ("`covariance'" == "covariance") {
            _coef_table_header, title(Naive covariance (Cov1)) width(62)
            display
            _coef_table, bmat(e(cov1)) vmat(e(V_cov_naive)) level(95) depname(`depvar') coeftitle(cov)
        }
        else {
            _coef_table_header, title(Naive (biased) estimate) width(62)
            display
            _coef_table, bmat(e(b_naive)) vmat(e(V_naive)) level(95) depname(`depvar') coeftitle(beta)
        }
        display
    }

    if ("`covariance'" == "covariance") {
        _coef_table_header, title(Debiased covariance (Cov1 - Cov0) relative to `baseline') width(62)
        display
        _coef_table, bmat(e(cov_debiased)) vmat(e(V_cov_debiased)) level(95) depname(`depvar') coeftitle(cov)
    }
    else {
        _coef_table_header, title(Denoised event study relative to `baseline') width(62)
        display
        _coef_table, bmat(e(b)) vmat(e(V)) level(95) depname(`depvar') coeftitle(beta)
    }
end
