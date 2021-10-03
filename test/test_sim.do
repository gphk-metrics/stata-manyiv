set seed 42
set linesize 112
mata mata clear
* cap noi net uninstall manyiv
* net install manyiv, from(`c(pwd)'/../src/)
qui include ../src/manyiv_absorb.mata
qui include ../src/manyiv_internals_m.mata
qui include ../src/manyiv.ado

clear
set obs 1000
gen u  = rnormal()
gen z1 = rnormal()
gen z2 = rnormal()
gen e  = rnormal() + u
gen c  = int(runiform() * 10)
gen fe = int(runiform() * 15)
gen iv = int(runiform() * 8)
gen w  = rnormal()
gen x  = 1 + 0.1 * z1 - 0.2 * z2 - 1/(1 + iv) + u
gen y  = 1 + x + w + 1/(1 + fe) + e

manyiv y (x = z1 z2 i.iv) w
matrix list e(b)
matrix list e(se)
manyiv y (x = z1 z2) w, absorbiv(iv)
matrix list e(b)
matrix list e(se)

manyiv y (x = z1 z2) w
manyiv y (x = z1 z2) w i.fe
matrix list e(b)
matrix list e(se)
manyiv y (x = z1 z2) w, absorb(fe)
matrix list e(b)
matrix list e(se)

manyiv y (x = z1 z2) w i.fe, cluster(c)
matrix list e(b)
matrix list e(se)
manyiv y (x = z1 z2) w, absorb(fe) cluster(c)
matrix list e(b)
matrix list e(se)

tab c,  gen(_c)
tab fe, gen(_fe)
jive y (x = z1 z2) w
jive y (x = z1 z2) w _fe*
jive y (x = z1 z2) w, r
jive y (x = z1 z2) w _fe*, r
