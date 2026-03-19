*! Test xt2denoise with simulated data

use "test/xt2denoise_testdata.dta", clear

xtset fake_id year

xt2denoise lnR, z(z) treatment(treatment) control(control) pre(3) post(3)

ereturn list
