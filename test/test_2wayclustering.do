set seed 42
set linesize 112
tempfile info

clear
set obs 1000
gen u   = rnormal()
gen z1  = rnormal()
gen z2  = rnormal()
gen e   = rnormal() + u
gen c   = int(runiform() * 10)
gen c2  = int(runiform() * 8)
gen fe  = int(runiform() * 15)
gen iv  = int(runiform() * 8)
gen fe2 = int(runiform() * 15)
gen iv2 = int(runiform() * 8)
gen w   = rnormal()
gen w2  = rnormal()
gen x   = 1 + 0.1 * z1 - 0.2 * z2 - 1/(1 + iv * iv2) + u
gen y   = 1 + x + w + 1/(1 + fe) - 1/(1 + fe2) + e
gen byte touse = _n > 2

mata mata clear
qui include ../src/mata/manyiv_absorb.mata
qui include ../src/mata/manyiv_internals.mata
qui include ../src/ado/manyiv.ado
manyiv y (x = z1 z2), cluster(c c2) nosmall
ivreghdfe y (x = z1 z2), cluster(c c2)
manyiv y (x = z1 z2), cluster(c c2)
ivreghdfe y (x = z1 z2), cluster(c c2) small
