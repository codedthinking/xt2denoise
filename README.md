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

- `xt2denoise` varname(numeric) [*if*], **treatment**(varname numeric) **control**(varname numeric), [**pre**(#) **post**(#) **baseline**(*string*) **weighting**(string) **cluster**(varname) **graph**]

The package can be installed with
```
net install xt2denoise, from(https://raw.githubusercontent.com/codedthinking/xt2denoise/main/) replace
```

# Options
## Options
Option | Description
-------|------------
**treatment** | Dummy variable indicating the treatment of interest.
**control** | Dummy variable indicating the control treatment.
**pre** | Number of periods before treatment to include in the estimation (default 1)
**post** | Number of periods after treatment to include in the estimation (default 3)
**baseline** | Either a negative number between `-pre` and `-1` or `average`, or `atet`. If `-k`, the baseline is the kth period before the treatment. If `average`, the baseline is the average of the pre-treatment periods. If `atet`, the regression table reports the average of the post-treatment periods minus the average of the pre-treatment periods. Default is `-1`.
**weighting** | Method to weight different cohorts in the estimation.
**cluster** | Variable to cluster standard errors. If not specified, the panel identifier is used.
**graph** (optional) | Plot the event study graph with the default settings of `hetdid_coefplot`.

# Examples
```
TBD
```

# Background
TBD

# Remarks
TBD

# Authors
- Miklós Koren (Central European University, https://koren.mk), *maintainer*

# License and Citation
You are free to use this package under the terms of its [license](https://github.com/codedthinking/xt2denoise/blob/main/LICENSE). If you use it, please cite the software package in your work:

- Koren, Miklós. (2026). XT2DENOISE - denoise second moments in panel event studies [Computer software]. Available at https://github.com/codedthinking/xt2denoise

Please also cite the paper describing the methodology:

- Koren, M., Orbán, K., & Telegdy, Á. (2025). CEO replacements and firm growth: Correcting for small-sample bias in fixed-effect estimates. Unpublished manuscript. https://koren.mk/publications/debiasing
