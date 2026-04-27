{smcl}


{marker xt2denoise-denoises-second-moments-in-panel-event-studies}{...}
{title:{cmd:xt2denoise} denoises second moments in panel event studies}

{pstd}This package implements the debiasing estimator of Koren, Orbán, and Telegdy (2025) to correct for small-sample bias in fixed-effect estimates of second moments in panel event studies.{p_end}


{marker syntax}{...}
{title:Syntax}

{text}{phang2}{cmd:xt2denoise} varname(numeric) [{it:if}], {bf:z}(varname numeric) {bf:treatment}(varname numeric) {bf:control}(varname numeric), [{bf:pre}(#) {bf:post}(#) {bf:baseline}({it:string}) {bf:cluster}(varname) {bf:graph} {bf:detail} {bf:cov} {bf:excessvariance} {bf:includenonchangers}]{p_end}


{pstd}The package can be installed with{p_end}

{p 8 16 2}net install xt2denoise, from(https://raw.githubusercontent.com/codedthinking/xt2denoise/main/) replace


{marker options}{...}
{title:Options}

{synoptset tabbed}{...}
{synopthdr:Option}
{synoptline}
{synopt:{bf:z}}Variable defining spell quality, used as a right-hand-side variable for denoising second moments.{p_end}
{synopt:{bf:treatment}}Dummy variable indicating the treatment of interest.{p_end}
{synopt:{bf:control}}Dummy variable indicating the control treatment.{p_end}
{synopt:{bf:pre}}Number of periods before treatment to include in the estimation (default 1){p_end}
{synopt:{bf:post}}Number of periods after treatment to include in the estimation (default 3){p_end}
{synopt:{bf:baseline}}Either a negative number between {cmd:-pre} and {cmd:-1} or {cmd:average}, or {cmd:atet}. If {cmd:-k}, the baseline is the kth period before the treatment. If {cmd:average}, the baseline is the average of the pre-treatment periods. If {cmd:atet}, the regression table reports the average of the post-treatment periods minus the average of the pre-treatment periods. Default is {cmd:-1}.{p_end}
{synopt:{bf:cluster}}Variable to cluster standard errors. If not specified, the panel identifier is used.{p_end}
{synopt:{bf:graph} (optional)}Plot the event study graph with the default settings of {cmd:hetdid_coefplot}.{p_end}
{synopt:{bf:detail} (optional)}Display both the naive (biased) and denoised estimates. The naive estimate uses only treated group moments (Cov1/Var_z1), while the denoised estimate differences out the control group.{p_end}
{synopt:{bf:cov} (optional)}Report covariance instead of beta coefficients. When specified, displays Cov(dy, dz) for naive (Cov1) and debiased (Cov1 - Cov0) estimates. Can be abbreviated as {cmd:cov} (minimum abbreviation of {cmd:covariance}).{p_end}
{synopt:{bf:excessvariance} (optional)}Apply excess variance correction. Use when the control group has different variance than the treatment group due to compositional differences. See below for details.{p_end}
{synopt:{bf:includenonchangers} (optional)}Include control group observations (with dz=0) in the naive estimator. By default, the naive estimator uses only treated observations. This option restores the old behavior of including controls with zero change in the naive variance and covariance calculations.{p_end}
{synoptline}


{marker examples}{...}
{title:Examples}

{pstd}Using simulated data with a known treatment effect (beta = 0 pre-treatment, beta = 1 post-treatment):{p_end}

{p 8 16 2}use "test/xt2denoise_testdata.dta", clear
xtset fake_id year
xt2denoise lnR, z(z) treatment(treatment) control(control) pre(3) post(3)

{pstd}Output:{p_end}

{p 8 16 2}Denoised event study relative to -1        Number of obs = 370

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

{pstd}The estimator correctly recovers the true parameters: approximately 0 in pre-treatment periods and approximately 1 in post-treatment periods.{p_end}


{marker background}{...}
{title:Background}

{pstd}In panel event studies, second moment estimates (variances, covariances) are biased due to small-sample noise in fixed-effect estimation. This package implements a debiasing estimator that uses a placebo control group to difference out the noise component.{p_end}

{pstd}The key insight is that the "true" variance of the quality change can be identified as:{p_end}

{p 8 16 2}true_Var(dz) = Var(dz | treated) - Var(dz | control)

{pstd}The debiased coefficient is then:{p_end}

{p 8 16 2}beta = (Cov(dy, dz | treated) - Cov(dy, dz | control)) / true_Var(dz)


{marker excess-variance-correction}{...}
{dlgtab:Excess variance correction}

{pstd}The baseline debiasing method assumes that the noise variance is the same in the treatment and control groups. This holds when both groups are random samples from the same population. In practice, however, the control group may differ systematically from the treatment group. For example, if the control group consists of smaller firms, they may have higher measurement noise.{p_end}

{pstd}The {cmd:excessvariance} option corrects for this by estimating the ratio of variances between the two groups from the pre-treatment period, when there should be no treatment effect. Let{p_end}

{p 8 16 2}c_z = Var(dz | treated, pre) / Var(dz | control, pre)
c_y = Var(dy | treated, pre) / Var(dy | control, pre)

{pstd}Under the assumption that variance differences are proportional (i.e., if the control group has twice the noise in z, it also has twice the noise in y), the estimator scales the control group moments by sqrt(c_z) and sqrt(c_y) before differencing. This ensures that the noise component is properly removed even when the two groups have different baseline variances.{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}The command returns coefficients and standard errors in {cmd:e(b)} and {cmd:e(V)}. Additional matrices are stored:{p_end}

{text}{phang2}{cmd:e(b_naive)}, {cmd:e(V_naive)}: Naive (biased) estimates and variance-covariance matrix{p_end}
{phang2}{cmd:e(cov1)}, {cmd:e(V_cov_naive)}: Naive covariance and its variance matrix{p_end}
{phang2}{cmd:e(cov_diff)}, {cmd:e(V_cov_diff)}: Debiased covariance (Cov1 - Cov0) and its variance matrix{p_end}
{phang2}{cmd:e(var_z1)}, {cmd:e(var_z_diff)}: Variance of dz (naive and differenced){p_end}
{phang2}{cmd:e(n1)}, {cmd:e(n0)}: Sample sizes by event time{p_end}


{pstd}Standard post-estimation commands can be used, such as {cmd:coefplot}, {cmd:esttab}, or {cmd:outreg2}.{p_end}


{marker authors}{...}
{title:Authors}

{text}{phang2}Miklós Koren (Central European University, {browse "https://koren.mk"}), {it:maintainer}{p_end}



{marker license-and-citation}{...}
{title:License and Citation}

{pstd}You are free to use this package under the terms of its {browse "https://github.com/codedthinking/xt2denoise/blob/main/LICENSE"}. If you use it, please cite the software package in your work:{p_end}

{text}{phang2}Koren, Miklós. (2026). XT2DENOISE - denoise second moments in panel event studies [Computer software]. Available at {browse "https://github.com/codedthinking/xt2denoise"}{p_end}


{pstd}Please also cite the paper describing the methodology:{p_end}

{text}{phang2}Koren, M., Orbán, K., & Telegdy, Á. (2025). CEO replacements and firm growth: Correcting for small-sample bias in fixed-effect estimates. Unpublished manuscript. {browse "https://koren.mk/publications/debiasing"}{p_end}
