clear
foreach var in y x z fe {
    gen `var' = .
}
cap noi manyiv y (x = .), absorbiv(z)
assert _rc == 2000
cap noi manyiv y (x = z), absorb(fe)
assert _rc == 2000
cap noi manyiv y (x = .), absorb(fe) absorbiv(z)
assert _rc == 2000

clear
set obs 1
foreach var in y x z fe {
    gen `var' = 1
}
cap noi manyiv y (x = .), absorbiv(z)
assert _rc == 2000
cap noi manyiv y (x = z), absorb(fe)
assert _rc == 2000
cap noi manyiv y (x = .), absorb(fe) absorbiv(z)
assert _rc == 2000

clear
sysuse auto
manyiv price (mpg = foreign), absorb(rep78)
cap noi manyiv price (mpg = foreign) if _n == 0, absorb(rep78)
assert _rc == 2000
replace price = .
cap noi manyiv price (mpg = foreign), absorb(rep78)
assert _rc == 2000
