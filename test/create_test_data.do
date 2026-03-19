*! Create test data for xt2denoise
*! Based on Monte Carlo simulation for placebo-controlled event study

clear all
set seed 2191

* Small test parameters
local rho0 = 0.5
local rho1 = 0.5
local sigma_epsilon0 = 0.1
local sigma_epsilon1 = 0.1
local hazard = 0.2
local T_max = 10
local N_changes = 50
local sigma_z = 0.3
local control_treated_ratio = 1
local delta_trend = 0.01

* Generate CEO changes
set obs `N_changes'
generate frame_id_numeric = _n

generate T1 = invexponential(1/`hazard', uniform())
generate T2 = invexponential(1/`hazard', uniform())
replace T1 = ceil(T1)
replace T2 = ceil(T2)

keep if T1 <= `T_max' & T2 <= `T_max'

* Construct placebo pairs
expand 1 + `control_treated_ratio', generate(placebo)
bysort frame_id_numeric (placebo): generate index = _n
egen fake_id = group(frame_id_numeric index)

* Add time dimension
expand T1 + T2
bysort fake_id: generate year = _n
generate byte ceo_spell = cond(year <= T1, 1, 2)

xtset fake_id year
generate change_year = T1 + 1

* Generate outcome with AR(1) process
generate trend_t = (year - change_year) * `delta_trend'
generate dlnR = rnormal(0, cond(placebo == 1, `sigma_epsilon0', `sigma_epsilon1'))
bysort fake_id (year): generate lnR = trend_t if _n == 1
bysort fake_id (year): replace lnR = trend_t + `rho1' * (lnR[_n-1] - trend_t[_n-1]) + dlnR if _n > 1 & placebo == 1
bysort fake_id (year): replace lnR = `rho0' * (lnR[_n-1] - trend_t[_n-1]) + dlnR if _n > 1 & placebo == 0

* z is lnR (spell quality measured by outcome)
generate z = lnR

* Measured manager skill
egen manager_skill = mean(lnR), by(fake_id ceo_spell)
summarize manager_skill if placebo == 0, meanonly
replace manager_skill = manager_skill - r(mean)

* Create treatment indicators
generate byte treatment = (placebo == 0) & (year == change_year)
generate byte control = (placebo == 1) & (year == change_year)

* Keep required variables
keep frame_id_numeric year lnR ceo_spell manager_skill change_year placebo fake_id treatment control z

* Save test data
save "test/xt2denoise_testdata.dta", replace

* Display summary
describe
summarize
tabulate placebo treatment
