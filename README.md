---
author: Koren, Miklós (https://koren.dev)
date: 2026-03-19
version: 0.1.0
title: XT2DENOISE - denoise second moments in panel event studies
description: |
    Denoises second moments in panel event study estimates.
url: https://github.com/codedthinking/xt2denoise
requires: Stata version 18
---
# `xt2denoise` denoises second moments in panel event studies

This package implements the debiasing estimator of Koren, Orbán, and Telegdy (2025) to correct for small-sample bias in fixed-effect estimates of second moments in panel event studies.

# Syntax

- `xt2denoise` varname(numeric) [*if*], **z**(varname numeric) **treatment**(varname numeric) **control**(varname numeric), [**pre**(#) **post**(#) **baseline**(*string*) **cluster**(varname) **graph**]

The package can be installed with
```
net install xt2denoise, from(https://raw.githubusercontent.com/codedthinking/xt2denoise/main/) replace
```

# Options
## Options
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

# Remarks

The command returns coefficients and standard errors in `e(b)` and `e(V)`. Additional matrices are stored:

- `e(cov1)`, `e(cov0)`: Covariances for treated and control groups
- `e(var_z1)`, `e(var_z0)`: Variance of dz for treated and control groups
- `e(true_var_z)`: True variance of dz (difference)
- `e(n1)`, `e(n0)`: Sample sizes by event time

Standard post-estimation commands can be used, such as `coefplot`, `esttab`, or `outreg2`.

# Authors
- Miklós Koren (Central European University, https://koren.mk), *maintainer*

# License and Citation
You are free to use this package under the terms of its [license](https://github.com/codedthinking/xt2denoise/blob/main/LICENSE). If you use it, please cite the software package in your work:

- Koren, Miklós. (2026). XT2DENOISE - denoise second moments in panel event studies [Computer software]. Available at https://github.com/codedthinking/xt2denoise

Please also cite the paper describing the methodology:

- Koren, M., Orbán, K., & Telegdy, Á. (2025). CEO replacements and firm growth: Correcting for small-sample bias in fixed-effect estimates. Unpublished manuscript. https://koren.mk/publications/debiasing
