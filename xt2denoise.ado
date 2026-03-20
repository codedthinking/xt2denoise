*! version 0.2.0 19mar2026
program xt2denoise, eclass
version 18.0

syntax varname(numeric) [if], z(varname numeric) treatment(varname numeric) control(varname numeric) [, pre(integer 1) post(integer 3) baseline(string) cluster(varname) graph detail]

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

***** STEP 4: Compute variance and covariance separately for treated and control
* Use factor variables: eventtime+100 for event study, or post (0/1) for ATET
* OLS of y^2 on constant gives variance estimate

tempvar dy2 dz2 dydz eventtime100 postvar
quietly generate `dy2' = `dy_demean'^2 if `touse'
quietly generate `dz2' = `dz'^2 if `touse'
quietly generate `dydz' = `dy_demean' * `dz' if `touse'
quietly generate `eventtime100' = `eventtime' + 100 if `touse'
quietly generate `postvar' = (`eventtime' >= 0) if `touse'

* determine number of coefficients based on baseline
if ("`baseline'" == "atet") {
    local K = 1
    local colnames "ATET"
}
else {
    local K = `pre' + `post' + 1
    local colnames ""
    forvalues t = -`pre'/`post' {
        local colnames `colnames' `t'
    }
}

* store results separately for treated (1) and control (0)
tempname cov1 cov0 var_z1 var_z0 se_cov1 se_cov0 se_var_z1 se_var_z0 n1_vec n0_vec
tempname true_var_z beta b V

if ("`baseline'" == "atet") {
    * ATET: use post (0/1) as the grouping variable
    matrix `cov1' = J(1, 2, .)
    matrix `cov0' = J(1, 2, .)
    matrix `var_z1' = J(1, 2, .)
    matrix `var_z0' = J(1, 2, .)
    matrix `se_cov1' = J(1, 2, .)
    matrix `se_cov0' = J(1, 2, .)
    matrix `se_var_z1' = J(1, 2, .)
    matrix `se_var_z0' = J(1, 2, .)

    * covariance for treated
    quietly regress `dydz' ibn.`postvar' if `touse' & `evert' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') nocons
    matrix `cov1'[1, 1] = _b[0.`postvar']
    matrix `cov1'[1, 2] = _b[1.`postvar']
    matrix `se_cov1'[1, 1] = _se[0.`postvar']
    matrix `se_cov1'[1, 2] = _se[1.`postvar']

    * covariance for control
    quietly regress `dydz' ibn.`postvar' if `touse' & `everc' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') nocons
    matrix `cov0'[1, 1] = _b[0.`postvar']
    matrix `cov0'[1, 2] = _b[1.`postvar']
    matrix `se_cov0'[1, 1] = _se[0.`postvar']
    matrix `se_cov0'[1, 2] = _se[1.`postvar']

    * variance of dz for treated
    quietly regress `dz2' ibn.`postvar' if `touse' & `evert' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') nocons
    matrix `var_z1'[1, 1] = _b[0.`postvar']
    matrix `var_z1'[1, 2] = _b[1.`postvar']
    matrix `se_var_z1'[1, 1] = _se[0.`postvar']
    matrix `se_var_z1'[1, 2] = _se[1.`postvar']

    * variance of dz for control
    quietly regress `dz2' ibn.`postvar' if `touse' & `everc' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') nocons
    matrix `var_z0'[1, 1] = _b[0.`postvar']
    matrix `var_z0'[1, 2] = _b[1.`postvar']
    matrix `se_var_z0'[1, 1] = _se[0.`postvar']
    matrix `se_var_z0'[1, 2] = _se[1.`postvar']

    * compute true variance and beta for pre and post
    tempname true_var_pre true_var_post beta_pre beta_post
    scalar `true_var_pre' = `var_z1'[1, 1] - `var_z0'[1, 1]
    scalar `true_var_post' = `var_z1'[1, 2] - `var_z0'[1, 2]

    if `true_var_pre' != 0 & `true_var_pre' != . {
        scalar `beta_pre' = (`cov1'[1, 1] - `cov0'[1, 1]) / `true_var_pre'
    }
    else {
        scalar `beta_pre' = 0
    }

    if `true_var_post' != 0 & `true_var_post' != . {
        scalar `beta_post' = (`cov1'[1, 2] - `cov0'[1, 2]) / `true_var_post'
    }
    else {
        scalar `beta_post' = .
    }

    * ATET = beta_post - beta_pre
    matrix `beta' = J(1, 1, .)
    matrix `beta'[1, 1] = `beta_post' - `beta_pre'

    * SE for ATET (simplified delta method)
    tempname se_pre se_post
    scalar `se_pre' = sqrt(`se_cov1'[1, 1]^2 + `se_cov0'[1, 1]^2) / abs(`true_var_pre')
    scalar `se_post' = sqrt(`se_cov1'[1, 2]^2 + `se_cov0'[1, 2]^2) / abs(`true_var_post')

    matrix `V' = J(1, 1, .)
    matrix `V'[1, 1] = `se_pre'^2 + `se_post'^2

    matrix colname `beta' = `colnames'
    matrix colname `V' = `colnames'
    matrix rowname `V' = `colnames'

    * compute naive ATET = (Cov1_post / Var_z1_post) - (Cov1_pre / Var_z1_pre)
    tempname beta_naive_pre beta_naive_post beta_naive V_naive se_naive_pre se_naive_post
    if `var_z1'[1, 1] != 0 & `var_z1'[1, 1] != . {
        scalar `beta_naive_pre' = `cov1'[1, 1] / `var_z1'[1, 1]
        scalar `se_naive_pre' = `se_cov1'[1, 1] / abs(`var_z1'[1, 1])
    }
    else {
        scalar `beta_naive_pre' = 0
        scalar `se_naive_pre' = 0
    }
    if `var_z1'[1, 2] != 0 & `var_z1'[1, 2] != . {
        scalar `beta_naive_post' = `cov1'[1, 2] / `var_z1'[1, 2]
        scalar `se_naive_post' = `se_cov1'[1, 2] / abs(`var_z1'[1, 2])
    }
    else {
        scalar `beta_naive_post' = .
        scalar `se_naive_post' = .
    }
    matrix `beta_naive' = J(1, 1, .)
    matrix `beta_naive'[1, 1] = `beta_naive_post' - `beta_naive_pre'
    matrix `V_naive' = J(1, 1, .)
    matrix `V_naive'[1, 1] = `se_naive_pre'^2 + `se_naive_post'^2
    matrix colname `beta_naive' = `colnames'
    matrix colname `V_naive' = `colnames'
    matrix rowname `V_naive' = `colnames'

    * sample sizes
    matrix `n1_vec' = J(1, 1, .)
    matrix `n0_vec' = J(1, 1, .)
    quietly count if `touse' & `evert' & inrange(`eventtime', -`pre', `post')
    matrix `n1_vec'[1, 1] = r(N)
    quietly count if `touse' & `everc' & inrange(`eventtime', -`pre', `post')
    matrix `n0_vec'[1, 1] = r(N)

    matrix `true_var_z' = J(1, 1, .)
    matrix `true_var_z'[1, 1] = `true_var_post' - `true_var_pre'
}
else {
    * Event study: use eventtime as the grouping variable
    matrix `cov1' = J(1, `K', .)
    matrix `cov0' = J(1, `K', .)
    matrix `var_z1' = J(1, `K', .)
    matrix `var_z0' = J(1, `K', .)
    matrix `se_cov1' = J(1, `K', .)
    matrix `se_cov0' = J(1, `K', .)
    matrix `se_var_z1' = J(1, `K', .)
    matrix `se_var_z0' = J(1, `K', .)
    matrix `n1_vec' = J(1, `K', .)
    matrix `n0_vec' = J(1, `K', .)
    matrix `true_var_z' = J(1, `K', .)
    matrix `beta' = J(1, `K', .)

    * covariance for treated (evert == 1)
    quietly regress `dydz' ibn.`eventtime100' if `touse' & `evert' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') nocons
    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        local e100 = `t' + 100
        capture matrix `cov1'[1, `col'] = _b[`e100'.`eventtime100']
        capture matrix `se_cov1'[1, `col'] = _se[`e100'.`eventtime100']
    }

    * covariance for control (everc == 1)
    quietly regress `dydz' ibn.`eventtime100' if `touse' & `everc' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') nocons
    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        local e100 = `t' + 100
        capture matrix `cov0'[1, `col'] = _b[`e100'.`eventtime100']
        capture matrix `se_cov0'[1, `col'] = _se[`e100'.`eventtime100']
    }

    * variance of dz for treated
    quietly regress `dz2' ibn.`eventtime100' if `touse' & `evert' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') nocons
    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        local e100 = `t' + 100
        capture matrix `var_z1'[1, `col'] = _b[`e100'.`eventtime100']
        capture matrix `se_var_z1'[1, `col'] = _se[`e100'.`eventtime100']
    }

    * variance of dz for control
    quietly regress `dz2' ibn.`eventtime100' if `touse' & `everc' & inrange(`eventtime', -`pre', `post'), cluster(`cluster') nocons
    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        local e100 = `t' + 100
        capture matrix `var_z0'[1, `col'] = _b[`e100'.`eventtime100']
        capture matrix `se_var_z0'[1, `col'] = _se[`e100'.`eventtime100']
    }

    * count sample sizes
    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        quietly count if `touse' & `eventtime' == `t' & `evert'
        matrix `n1_vec'[1, `col'] = r(N)
        quietly count if `touse' & `eventtime' == `t' & `everc'
        matrix `n0_vec'[1, `col'] = r(N)
    }

    * compute true variance and beta (debiased)
    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        * true variance of dz = Var(dz1) - Var(dz0)
        matrix `true_var_z'[1, `col'] = `var_z1'[1, `col'] - `var_z0'[1, `col']
        * beta = (Cov1 - Cov0) / true_Var(dz)
        if `true_var_z'[1, `col'] != 0 & `true_var_z'[1, `col'] != . {
            matrix `beta'[1, `col'] = (`cov1'[1, `col'] - `cov0'[1, `col']) / `true_var_z'[1, `col']
        }
    }

    * compute naive beta = Cov1 / Var_z1 (biased, without denoising)
    tempname beta_naive se_beta_naive V_naive
    matrix `beta_naive' = J(1, `K', .)
    matrix `se_beta_naive' = J(1, `K', .)

    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        if `var_z1'[1, `col'] != 0 & `var_z1'[1, `col'] != . {
            matrix `beta_naive'[1, `col'] = `cov1'[1, `col'] / `var_z1'[1, `col']
            matrix `se_beta_naive'[1, `col'] = `se_cov1'[1, `col'] / abs(`var_z1'[1, `col'])
        }
    }

    matrix colname `beta' = `colnames'
    matrix colname `beta_naive' = `colnames'
    matrix colname `cov1' = `colnames'
    matrix colname `cov0' = `colnames'
    matrix colname `var_z1' = `colnames'
    matrix colname `var_z0' = `colnames'
    matrix colname `true_var_z' = `colnames'
    matrix colname `n1_vec' = `colnames'
    matrix colname `n0_vec' = `colnames'

    ***** Compute standard errors for debiased beta using delta method
    * beta = (cov1 - cov0) / (var_z1 - var_z0)
    * SE(beta) approx = sqrt( (se_cov1^2 + se_cov0^2) / true_var_z^2 )
    tempname se_beta
    matrix `se_beta' = J(1, `K', .)

    forvalues t = -`pre'/`post' {
        local col = `t' + `pre' + 1
        if `true_var_z'[1, `col'] != 0 & `true_var_z'[1, `col'] != . {
            matrix `se_beta'[1, `col'] = sqrt(`se_cov1'[1, `col']^2 + `se_cov0'[1, `col']^2) / abs(`true_var_z'[1, `col'])
        }
    }
    matrix colname `se_beta' = `colnames'
    matrix colname `se_beta_naive' = `colnames'

    * build variance matrix for debiased beta (diagonal)
    matrix `V' = J(`K', `K', 0)
    forvalues i = 1/`K' {
        if `se_beta'[1, `i'] != . {
            matrix `V'[`i', `i'] = `se_beta'[1, `i']^2
        }
    }
    matrix colname `V' = `colnames'
    matrix rowname `V' = `colnames'

    * build variance matrix for naive beta (diagonal)
    matrix `V_naive' = J(`K', `K', 0)
    forvalues i = 1/`K' {
        if `se_beta_naive'[1, `i'] != . {
            matrix `V_naive'[`i', `i'] = `se_beta_naive'[1, `i']^2
        }
    }
    matrix colname `V_naive' = `colnames'
    matrix rowname `V_naive' = `colnames'
}

* count observations
quietly count if `touse' & inrange(`eventtime', -`pre', `post')
local Nobs = r(N)

tempvar esample
quietly generate `esample' = `touse' & inrange(`eventtime', -`pre', `post')

* post results
ereturn post `beta' `V', obs(`Nobs') esample(`esample')
ereturn local depvar "`y'"
ereturn local cmd "xt2denoise"
ereturn local cmdline "xt2denoise `0'"

* store additional matrices
ereturn matrix cov1 = `cov1'
ereturn matrix cov0 = `cov0'
ereturn matrix var_z1 = `var_z1'
ereturn matrix var_z0 = `var_z0'
ereturn matrix true_var_z = `true_var_z'
ereturn matrix n1 = `n1_vec'
ereturn matrix n0 = `n0_vec'

* store naive beta matrices
ereturn matrix b_naive = `beta_naive'
ereturn matrix V_naive = `V_naive'

* display regression table(s)
if ("`detail'" == "detail") {
    * show naive (biased) beta first
    _coef_table_header, title(Naive (biased) estimate) width(62)
    display
    _coef_table, bmat(e(b_naive)) vmat(e(V_naive)) level(95) depname(`y') coeftitle(beta)
    display
}

_coef_table_header, title(Denoised event study relative to `baseline') width(62)
display
_coef_table, bmat(e(b)) vmat(e(V)) level(95) depname(`y') coeftitle(beta)

end
