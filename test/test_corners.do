clear
foreach var in y x z fe {
    gen `var' = .
}
cap manyiv y (x = .), absorbiv(z)
assert _rc == 1234
cap manyiv y (x = z), absorb(fe)
assert _rc == 1234
cap manyiv y (x = .), absorb(fe) absorbiv(z)
assert _rc == 1234
