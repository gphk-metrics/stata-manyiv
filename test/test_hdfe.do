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
gen fe  = int(runiform() * 15)
gen iv  = int(runiform() * 8)
gen fe2 = int(runiform() * 15)
gen iv2 = int(runiform() * 8)
gen w   = rnormal()
gen w2  = rnormal()
gen x   = 1 + 0.1 * z1 - 0.2 * z2 - 1/(1 + iv * iv2) + u
gen y   = 1 + x + w + 1/(1 + fe) - 1/(1 + fe2) + e
gen byte touse = _n > 2

* mata mata clear
* qui do ../src/mata/manyiv_absorb.mata
mata A  = ManyIVreg_Absorb_New(tokens("fe iv"), "touse", 1)
mata A.dropsingletons(.)
mata A.encode()
mata A.flagredundant()
mata maxiter = A.maxiter
mata hdfetol = A.hdfetol
mata X = st_data(., "z1 z2", "touse")
mata A._hdfe_squarem(X1 = X)
mata A._hdfe_fpi(X2 = X)
mata A._hdfe_cg(X3 = X)
mata colmax(abs(X1 :- X3))
mata colmax(abs(X1 :- X2))
