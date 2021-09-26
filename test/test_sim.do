set seed 42
set linesize 112
mata mata clear
qui include ../src/manyiv_internals_m.mata
qui include ../src/manyiv.ado

clear
set obs 1000
gen u  = rnormal()
gen z1 = rnormal()
gen z2 = rnormal()
gen e  = rnormal() + u
gen c  = mod(_n, 10)
gen w  = rnormal()
gen x  = 1 + z1 - z2 + u
gen y  = 1 + x + w + e
reg y x w
ivreg y (x = z1 z2) w
manyiv y (x = z1 z2) w, save(rsim)
manyiv y (x = z1 z2) w, save(rsim) cluster(c)
