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

replace fe = 99 in 10/11
replace iv = 98 in 10
replace fe = 5.5 in 100/101
replace iv = 5.5 in 100
replace fe = -9999 in 200
replace iv = -9998 in 300
replace fe = 1234 if iv == 3 & _n == 95

* cd ..
* shell make all
* cd test

if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
else local c_os_: di lower("`c(os)'")

cap program drop manyiv_plugin
program manyiv_plugin, plugin using("../src/build/manyiv_`c_os_'.plugin")

mata mata clear
qui do ../src/mata/manyiv_absorb.mata
mata A  = ManyIVreg_Absorb_New(tokens("fe"), "touse", 1)
mata IV = ManyIVreg_Absorb_New(tokens("iv"), "touse", 1)
mata IV.append(A)
* mata IV.dropsingletons(singleix=J(0, 1, .), 1)
* mata A.dropfromindex(singleix)
mata IV.dropsingletons(., 2)
mata A.dropsingletons(., 2)
mata IV.encode()
mata A.encode()

plugin call manyiv_plugin, `"_plugin_check"'

mata IV.exportc("`info'")
mata dtest = IV.d_projection()
plugin call manyiv_plugin, `"`info'"'
mata IV.importc("`info'")
mata max(abs(dtest - IV.d_projection))

mata A.exportc("`info'")
mata dtest = A.d_projection()
plugin call manyiv_plugin, `"`info'"'
mata A.importc("`info'")
mata max(abs(dtest - A.d_projection))

mata B = ManyIVreg_Absorb_New(tokens("iv iv2 fe fe2"), "touse", 1)
mata B.encode()
mata B.dropsingletons(., 2)
mata D1 = select(designmatrix(B.groupid(1)), !B.skip(1)')
mata D2 = select(designmatrix(B.groupid(2)), !B.skip(2)')
mata D3 = select(designmatrix(B.groupid(3)), !B.skip(3)')
mata D4 = select(designmatrix(B.groupid(4)), !B.skip(4)')
mata DD = (D1, D2, D3, D4)
mata ZZ = DD' * DD
mata dd = rowsum((DD * invsym(ZZ)) :* DD)

mata B.exportc("`info'")
plugin call manyiv_plugin, `"`info'"'
mata B.importc("`info'")
mata max(abs(dd - B.d_projection))
