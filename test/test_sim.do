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

manyiv y (x = z1 z2 i.iv) w, cluster(c)
mata b1  = st_matrix("e(b)")
mata se1 = st_matrix("e(se)")
mata F1  = st_numscalar("e(F)")
manyiv y (x = z1 z2) w, absorbiv(iv) cluster(c)
mata b2  = st_matrix("e(b)")
mata se2 = st_matrix("e(se)")
mata F2  = st_numscalar("e(F)")
mata max(abs(b2  :- b1))
mata max(abs(se2 :- se1))
mata max(abs(F2  :- F1))

manyiv y (x = z1 z2) w i.fe,  cluster(c)
mata b1  = st_matrix("e(b)")
mata se1 = st_matrix("e(se)")
mata F1  = st_numscalar("e(F)")
manyiv y (x = z1 z2) w, absorb(fe) cluster(c)
mata b2  = st_matrix("e(b)")
mata se2 = st_matrix("e(se)")
mata F2  = st_numscalar("e(F)")
mata max(abs(b2  :- b1))
mata max(abs(se2 :- se1))
mata max(abs(F2  :- F1))

manyiv y (x = z1 z2 i.iv) w i.fe, cluster(c)
mata b1  = st_matrix("e(b)")
mata se1 = st_matrix("e(se)")
mata F1  = st_numscalar("e(F)")
manyiv y (x = z1 z2 i.iv) w, absorb(fe) cluster(c)
mata b2  = st_matrix("e(b)")
mata se2 = st_matrix("e(se)")
mata F2  = st_numscalar("e(F)")
manyiv y (x = z1 z2) i.fe w, absorbiv(iv) cluster(c)
mata b3  = st_matrix("e(b)")
mata se3 = st_matrix("e(se)")
mata F3  = st_numscalar("e(F)")
manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c)
mata b4  = st_matrix("e(b)")
mata se4 = st_matrix("e(se)")
mata F4  = st_numscalar("e(F)")

mata max(abs(b2  :- b1))
mata max(abs(se2 :- se1))
mata max(abs(F2  :- F1))
mata max(abs(b3  :- b1))
mata max(abs(se3 :- se1))
mata max(abs(F3  :- F1))
mata max(abs(b4  :- b1))
mata max(abs(se4 :- se1))
mata max(abs(F4  :- F1))

tab c,  gen(_c)
tab fe, gen(_fe)
tab iv, gen(_iv)
jive y (x = z1 z2) w
jive y (x = z1 z2) w _fe*
jive y (x = z1 z2) w, r
jive y (x = z1 z2) w _fe*, r
jive y (x = z1 z2 _iv*) w, r
jive y (x = z1 z2 _iv*) w _fe*, r
