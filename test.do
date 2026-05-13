net install xt2denoise, from(`=c(pwd)') replace
use "test/xt2denoise_testdata.dta", clear
xt2denoise lnR, z(lnR) treatment(treatment) control(control) baseline(atet) pre(3) post(3) cov detail 