*! version 0.5.0 20mar2026
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
local baseline = lower("`baseline'")
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
tempvar dz_mean dz_demean dz_naive dz_naive_mean dz_naive_demean dy2 dz2 dydz dz2_naive dydz_naive eventtime100 postvar

* demean dz by event time and treatment group for debiased estimator
quietly egen `dz_mean' = mean(`dz') if `touse', by(`eventtime' `evert')
quietly generate `dz_demean' = `dz' - `dz_mean' if `touse'

* for naive estimator: dz = 0 for controls, demean by eventtime only (full sample)
quietly generate `dz_naive' = `dz' * `evert' if `touse'
quietly egen `dz_naive_mean' = mean(`dz_naive') if `touse', by(`eventtime')
quietly generate `dz_naive_demean' = `dz_naive' - `dz_naive_mean' if `touse'

quietly generate `dy2' = `dy_demean'^2 if `touse'
quietly generate `dz2' = `dz_demean'^2 if `touse'
quietly generate `dydz' = `dy_demean' * `dz_demean' if `touse'
quietly generate `dz2_naive' = `dz_naive_demean'^2 if `touse'
quietly generate `dydz_naive' = `dy_demean' * `dz_naive_demean' if `touse'
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
tempname cov1 var_z1 se_cov1 se_var_z1 V_cov1 V_var_z1
tempname cov_diff var_z_diff se_cov_diff se_var_z_diff V_cov_diff V_var_z_diff

* estimate cov1 and var_z1 for naive estimator (full sample, dz=0 for controls)
_xt2denoise_regcoef `dydz_naive' `groupvar' if `touse' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') k(`Kreg') pre(`pre')
matrix `cov1' = r(coef)
matrix `se_cov1' = r(se)
matrix `V_cov1' = r(V)

_xt2denoise_regcoef `dz2_naive' `groupvar' if `touse' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') k(`Kreg') pre(`pre')
matrix `var_z1' = r(coef)
matrix `se_var_z1' = r(se)
matrix `V_var_z1' = r(V)

* estimate differences directly using areg with interaction (all groups)
_xt2denoise_regdiff `dydz' `groupvar' `evert' if `touse' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') k(`Kreg') pre(`pre')
matrix `cov_diff' = r(coef)
matrix `se_cov_diff' = r(se)
matrix `V_cov_diff' = r(V)

_xt2denoise_regdiff `dz2' `groupvar' `evert' if `touse' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') k(`Kreg') pre(`pre')
matrix `var_z_diff' = r(coef)
matrix `se_var_z_diff' = r(se)
matrix `V_var_z_diff' = r(V)

* compute true variance, beta, and standard errors
tempname true_var_z beta se_beta V_beta beta_naive se_beta_naive n1_vec n0_vec

* true_var_z is now estimated directly as the difference
matrix `true_var_z' = `var_z_diff'

_xt2denoise_compute_beta `cov_diff' `true_var_z' `V_cov_diff'
matrix `beta' = r(beta)
matrix `se_beta' = r(se)
matrix `V_beta' = r(V)

* naive beta = cov1 / var_z1 (no differencing)
matrix `beta_naive' = J(1, `Kreg', 0)
matrix `se_beta_naive' = J(1, `Kreg', 0)
forvalues i = 1/`Kreg' {
    if `var_z1'[1, `i'] != 0 & `var_z1'[1, `i'] != . {
        matrix `beta_naive'[1, `i'] = `cov1'[1, `i'] / `var_z1'[1, `i']
        matrix `se_beta_naive'[1, `i'] = `se_cov1'[1, `i'] / abs(`var_z1'[1, `i'])
    }
    * if var_z1 == 0, beta_naive and se stay at 0 (baseline case)
}

* for ATET: compute difference (post - pre) with proper covariance
if ("`baseline'" == "atet") {
    tempname beta_atet se_atet beta_naive_atet se_naive_atet
    scalar `beta_atet' = `beta'[1, 2] - `beta'[1, 1]
    * SE(post - pre) = sqrt(Var(post) + Var(pre) - 2*Cov(pre, post))
    scalar `se_atet' = sqrt(`V_beta'[1, 1] + `V_beta'[2, 2] - 2 * `V_beta'[1, 2])
    scalar `beta_naive_atet' = `beta_naive'[1, 2] - `beta_naive'[1, 1]
    * naive uses V_cov1 covariance (same regression)
    scalar `se_naive_atet' = sqrt(`V_cov1'[1, 1] / `var_z1'[1, 1]^2 + `V_cov1'[2, 2] / `var_z1'[1, 2]^2 - 2 * `V_cov1'[1, 2] / (`var_z1'[1, 1] * `var_z1'[1, 2]))

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
matrix colname `cov_diff' = `colnames'
matrix colname `se_cov1' = `colnames'
matrix colname `se_cov_diff' = `colnames'
matrix colname `var_z1' = `colnames'
matrix colname `var_z_diff' = `colnames'
matrix colname `true_var_z' = `colnames'
matrix colname `n1_vec' = `colnames'
matrix colname `n0_vec' = `colnames'
matrix colname `V_cov_diff' = `colnames'
matrix rowname `V_cov_diff' = `colnames'

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
ereturn matrix cov_diff = `cov_diff'
ereturn matrix se_cov1 = `se_cov1'
ereturn matrix se_cov_diff = `se_cov_diff'
ereturn matrix V_cov_diff = `V_cov_diff'
ereturn matrix var_z1 = `var_z1'
ereturn matrix var_z_diff = `var_z_diff'
ereturn matrix true_var_z = `true_var_z'
ereturn matrix n1 = `n1_vec'
ereturn matrix n0 = `n0_vec'
ereturn matrix b_naive = `beta_naive'
ereturn matrix V_naive = `V_naive'

* compute covariance matrices for display (cov_diff is already the debiased covariance)
tempname V_cov_naive cov_debiased_copy V_cov_debiased_copy
local Kcov = colsof(e(cov1))
matrix `V_cov_naive' = J(`Kcov', `Kcov', 0)
forvalues i = 1/`Kcov' {
    matrix `V_cov_naive'[`i', `i'] = e(se_cov1)[1, `i']^2
}
matrix colname `V_cov_naive' = `colnames'
matrix rowname `V_cov_naive' = `colnames'

* copy matrices for cov_debiased (cov_diff already stored, just make aliases)
matrix `cov_debiased_copy' = e(cov_diff)
matrix `V_cov_debiased_copy' = e(V_cov_diff)

ereturn matrix cov_debiased = `cov_debiased_copy'
ereturn matrix V_cov_debiased = `V_cov_debiased_copy'
ereturn matrix V_cov_naive = `V_cov_naive'

* display results
_xt2denoise_display, baseline(`baseline') depvar(`y') `detail' `covariance'

end

***** Helper program: run regression and extract coefficients
program _xt2denoise_regcoef, rclass
    syntax varlist(min=2 max=2) [if], cluster(varname) k(integer) pre(integer)

    gettoken depvar groupvar : varlist

    quietly regress `depvar' ibn.`groupvar' `if', cluster(`cluster') nocons

    tempname coef se V
    matrix `coef' = J(1, `k', .)
    matrix `se' = J(1, `k', .)
    matrix `V' = J(`k', `k', .)

    * check if this is ATET (k=2 with postvar) or event study
    if (`k' == 2) {
        * ATET: groupvar is postvar (0/1)
        capture matrix `coef'[1, 1] = _b[0.`groupvar']
        capture matrix `coef'[1, 2] = _b[1.`groupvar']
        capture matrix `se'[1, 1] = _se[0.`groupvar']
        capture matrix `se'[1, 2] = _se[1.`groupvar']
        * extract full variance-covariance matrix
        capture matrix `V'[1, 1] = e(V)[1, 1]
        capture matrix `V'[1, 2] = e(V)[1, 2]
        capture matrix `V'[2, 1] = e(V)[2, 1]
        capture matrix `V'[2, 2] = e(V)[2, 2]
    }
    else {
        * event study: groupvar is eventtime100
        forvalues i = 1/`k' {
            local t = `i' - `pre' - 1
            local e100 = `t' + 100
            capture matrix `coef'[1, `i'] = _b[`e100'.`groupvar']
            capture matrix `se'[1, `i'] = _se[`e100'.`groupvar']
        }
        * extract full variance-covariance matrix
        forvalues i = 1/`k' {
            forvalues j = 1/`k' {
                local ti = `i' - `pre' - 1
                local tj = `j' - `pre' - 1
                local e100i = `ti' + 100
                local e100j = `tj' + 100
                capture matrix `V'[`i', `j'] = e(V)["`e100i'.`groupvar'", "`e100j'.`groupvar'"]
            }
        }
    }

    return matrix coef = `coef'
    return matrix se = `se'
    return matrix V = `V'
end

***** Helper program: compute beta = cov_diff / true_var using direct difference estimates
program _xt2denoise_compute_beta, rclass
    args cov_diff true_var V_cov_diff

    local K = colsof(`cov_diff')
    tempname beta se V
    matrix `beta' = J(1, `K', 0)
    matrix `se' = J(1, `K', 0)
    matrix `V' = J(`K', `K', 0)

    * compute beta and full variance-covariance matrix
    forvalues i = 1/`K' {
        forvalues j = 1/`K' {
            if `true_var'[1, `i'] != 0 & `true_var'[1, `i'] != . & `true_var'[1, `j'] != 0 & `true_var'[1, `j'] != . {
                * Cov(beta[i], beta[j]) = V_cov_diff[i,j] / (true_var[i] * true_var[j])
                matrix `V'[`i', `j'] = `V_cov_diff'[`i', `j'] / (`true_var'[1, `i'] * `true_var'[1, `j'])
            }
        }
    }

    forvalues i = 1/`K' {
        if `true_var'[1, `i'] != 0 & `true_var'[1, `i'] != . {
            matrix `beta'[1, `i'] = `cov_diff'[1, `i'] / `true_var'[1, `i']
            matrix `se'[1, `i'] = sqrt(`V'[`i', `i'])
        }
        * if true_var == 0, beta and se stay at 0 (baseline case)
    }

    return matrix beta = `beta'
    return matrix se = `se'
    return matrix V = `V'
end

***** Helper program: estimate difference (treated - control) using areg with interaction
program _xt2denoise_regdiff, rclass
    syntax varlist(min=3 max=3) [if], cluster(varname) k(integer) pre(integer)

    gettoken depvar rest : varlist
    gettoken groupvar evert : rest
    local evert = strtrim("`evert'")

    * estimate difference directly: absorb event-time FE, interact with evert
    quietly areg `depvar' c.`evert'#ibn.`groupvar' `if', absorb(`groupvar') cluster(`cluster')

    tempname coef se V b_full V_full
    matrix `coef' = J(1, `k', .)
    matrix `se' = J(1, `k', .)
    matrix `V' = J(`k', `k', .)

    * extract coefficients by position (areg uses tempvar names, so we can't use named references)
    matrix `b_full' = e(b)
    matrix `V_full' = e(V)

    * coefficients are in order: level1#c.evert, level2#c.evert, ..., _cons
    * for ATET (k=2): positions 1 (pre/0) and 2 (post/1)
    * for event study: positions 1 to k (eventtime -pre to +post)
    forvalues i = 1/`k' {
        matrix `coef'[1, `i'] = `b_full'[1, `i']
        matrix `se'[1, `i'] = sqrt(`V_full'[`i', `i'])
    }

    * extract full variance-covariance matrix
    forvalues i = 1/`k' {
        forvalues j = 1/`k' {
            matrix `V'[`i', `j'] = `V_full'[`i', `j']
        }
    }

    return matrix coef = `coef'
    return matrix se = `se'
    return matrix V = `V'
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
