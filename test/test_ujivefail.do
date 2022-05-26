clear
input z1 z2 fe iv
-1 -1 0 0
1  1 0 0
3  3 1 0
-3 -3 1 1
0  0 2 1
0  0 2 1
end
gen z3 = _n^2
gen x  = z1 + z2 + z3 + rnormal()
gen y  = x + rnormal()

manyiv y (x = z1 z2), absorb(fe)
manyiv y (x = z1 z2), absorbiv(iv)
manyiv y (x = z3),    absorb(fe) absorbiv(iv)
manyiv y (x = z1 z2), absorb(fe) absorbiv(iv)
cap noi manyiv y (x = z1 z2), absorb(fe) absorbiv(iv) forcejive
assert _rc == 481

mata A = ManyIVreg_Absorb_New(tokens("iv fe"), "")
mata Z = A.hdfe(st_data(., "z1 z2", ""))
mata D = 1 :- (rowsum((Z * invsym(Z' * Z)) :* Z) + A.d_projection())
mata D

set obs 8
replace z1 = 1 in 7
replace z1 = 2 in 8
replace z2 = 2 in 7
replace z2 = 3 in 8
replace fe = 3 in 7/8
replace iv = 2 in 7/8
replace z3 = _n^2
replace x  = z1 + z2 + z3 + rnormal() in 7/8
replace y  = x + rnormal() in 7/8

manyiv y (x = z1 z2), absorb(fe) absorbiv(iv)
manyiv y (x = z1 z2), absorb(fe) absorbiv(iv) forcejive

clear
set obs 4
input z
-1
1
0
0
end
gen fe = _n > 2
gen x  = z + rnormal()
gen y  = x + rnormal()
manyiv y (x = z), absorb(fe)
cap noi manyiv y (x = z), absorb(fe) forcejive
assert _rc == 481

mata A = ManyIVreg_Absorb_New(tokens("fe"), "")
mata Z = A.hdfe(st_data(., "z", ""))
mata D = 1 :- (rowsum((Z * invsym(Z' * Z)) :* Z) + A.d_projection())
mata D

set obs 6
replace z = 1 in 5
replace z = 2 in 6
replace fe = _n > 2
replace x = z + rnormal() in 5/6
replace y = x + rnormal() in 5/6
