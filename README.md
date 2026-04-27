---
author: Koren, Miklós (https://koren.dev)
date: 2026-04-27
version: 0.8.0
title: XT2DENOISE - denoise second moments in panel event studies
description: |
    Denoises second moments in panel event study estimates.
url: https://github.com/codedthinking/xt2denoise
requires: Stata version 18
---
# `xt2denoise` denoises second moments in panel event studies

This package implements the debiasing estimator of Koren, Orbán, and Telegdy (2025) to correct for small-sample bias in fixed-effect estimates of second moments in panel event studies.

# Syntax

- `xt2denoise` varname(numeric) [*if*], **z**(varname numeric) **treatment**(varname numeric) **control**(varname numeric), [**pre**(#) **post**(#) **baseline**(*string*) **cluster**(varname) **graph** **detail** **cov** **excessvariance** **includenonchangers**]

The package can be installed with
```
net install xt2denoise, from(https://raw.githubusercontent.com/codedthinking/xt2denoise/main/) replace
```

# Options
Option | Description
-------|------------
**z** | Variable defining spell quality, used as a right-hand-side variable for denoising second moments.
**treatment** | Dummy variable indicating the treatment of interest.
**control** | Dummy variable indicating the control treatment.
**pre** | Number of periods before treatment to include in the estimation (default 1)
**post** | Number of periods after treatment to include in the estimation (default 3)
**baseline** | Either a negative number between `-pre` and `-1` or `average`, or `atet`. If `-k`, the baseline is the kth period before the treatment. If `average`, the baseline is the average of the pre-treatment periods. If `atet`, the regression table reports the average of the post-treatment periods minus the average of the pre-treatment periods. Default is `-1`.
**cluster** | Variable to cluster standard errors. If not specified, the panel identifier is used.
**graph** (optional) | Plot the event study graph with the default settings of `hetdid_coefplot`.
**detail** (optional) | Display both the naive (biased) and denoised estimates. The naive estimate uses only treated group moments (Cov1/Var_z1), while the denoised estimate differences out the control group.
**cov** (optional) | Report covariance instead of beta coefficients. When specified, displays Cov(dy, dz) for naive (Cov1) and debiased (Cov1 - Cov0) estimates. Can be abbreviated as `cov` (minimum abbreviation of `covariance`).
**excessvariance** (optional) | Apply excess variance correction. Use when the control group has different variance than the treatment group due to compositional differences. See below for details.
**includenonchangers** (optional) | Include control group observations (with dz=0) in the naive estimator. By default, the naive estimator uses only treated observations. This option restores the old behavior of including controls with zero change in the naive variance and covariance calculations.

# Examples

Using simulated data with a known treatment effect (beta = 0 pre-treatment, beta = 1 post-treatment):

```
use "test/xt2denoise_testdata.dta", clear
xtset fake_id year
xt2denoise lnR, z(z) treatment(treatment) control(control) pre(3) post(3)
```

Output:
```
Denoised event study relative to -1        Number of obs = 370

------------------------------------------------------------------------------
         lnR |       beta   Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
          -3 |      .0788   .0728743     1.08   0.280     -.064031    .2216309
          -2 |   .0864862   .0636792     1.36   0.174    -.0383228    .2112951
          -1 |          0  (omitted)
           0 |   1.026701   .2403755     4.27   0.000     .5555735    1.497828
           1 |   .8458092   .2522092     3.35   0.001     .3514882     1.34013
           2 |   .8696062   .2358194     3.69   0.000     .4074088    1.331804
           3 |   .8820123   .2475426     3.56   0.000     .3968377    1.367187
------------------------------------------------------------------------------
```

The estimator correctly recovers the true parameters: approximately 0 in pre-treatment periods and approximately 1 in post-treatment periods.

# Background

In panel event studies, second moment estimates (variances, covariances) are biased due to small-sample noise in fixed-effect estimation. This package implements a debiasing estimator that uses a placebo control group to difference out the noise component.

The key insight is that the "true" variance of the quality change can be identified as:

    true_Var(dz) = Var(dz | treated) - Var(dz | control)

The debiased coefficient is then:

    beta = (Cov(dy, dz | treated) - Cov(dy, dz | control)) / true_Var(dz)

## Excess variance correction

The baseline debiasing method assumes that the noise variance is the same in the treatment and control groups. This holds when both groups are random samples from the same population. In practice, however, the control group may differ systematically from the treatment group. For example, if the control group consists of smaller firms, they may have higher measurement noise.

The `excessvariance` option corrects for this by estimating the ratio of variances between the two groups from the pre-treatment period, when there should be no treatment effect. Let

    c_z = Var(dz | treated, pre) / Var(dz | control, pre)
    c_y = Var(dy | treated, pre) / Var(dy | control, pre)

Under the assumption that variance differences are proportional (i.e., if the control group has twice the noise in z, it also has twice the noise in y), the estimator scales the control group moments by sqrt(c_z) and sqrt(c_y) before differencing. This ensures that the noise component is properly removed even when the two groups have different baseline variances.

# Remarks

The command returns coefficients and standard errors in `e(b)` and `e(V)`. Additional matrices are stored:

- `e(b_naive)`, `e(V_naive)`: Naive (biased) estimates and variance-covariance matrix
- `e(cov1)`, `e(V_cov_naive)`: Naive covariance and its variance matrix
- `e(cov_diff)`, `e(V_cov_diff)`: Debiased covariance (Cov1 - Cov0) and its variance matrix
- `e(var_z1)`, `e(var_z_diff)`: Variance of dz (naive and differenced)
- `e(n1)`, `e(n0)`: Sample sizes by event time

Standard post-estimation commands can be used, such as `coefplot`, `esttab`, or `outreg2`.

# Authors
- Miklós Koren (Central European University, https://koren.mk), *maintainer*

# License and Citation
You are free to use this package under the terms of its [license](https://github.com/codedthinking/xt2denoise/blob/main/LICENSE). If you use it, please cite the software package in your work:

- Koren, Miklós. (2026). XT2DENOISE - denoise second moments in panel event studies [Computer software]. Available at https://github.com/codedthinking/xt2denoise

Please also cite the paper describing the methodology:

- Koren, M., Orbán, K., & Telegdy, Á. (2025). CEO replacements and firm growth: Correcting for small-sample bias in fixed-effect estimates. Unpublished manuscript. https://koren.mk/publications/debiasing
