set seed 42
set linesize 112
* mata mata clear
* cap noi net uninstall manyiv
* net install manyiv, from(`c(pwd)'/../src/)
*
* mata mata clear
* qui include ../src/manyiv_absorb.mata
* qui include ../src/manyiv_internals_m.mata
* qui include ../src/manyiv.ado

capture program drop main
program main
    consistency
end

***********************************************************************
*                                                                     *
*                        Internal Consistency                         *
*                                                                     *
***********************************************************************

capture program drop consistency
program consistency
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

    manyiv y (x = z1 z2), cluster(c)
    manyiv y (x = z1 z2), cluster(c) noc
    mata errors = J(0, 4, .)
    mata errow  = 0

    manyiv y (x = z1 z2) i.fe, cluster(c)
    mata b1  = st_matrix("e(b)")
    mata se1 = st_matrix("e(se)")
    mata F1  = st_numscalar("e(F)")
    manyiv y (x = z1 z2), absorb(fe) cluster(c)
    mata b2  = st_matrix("e(b)")
    mata se2 = st_matrix("e(se)")
    mata F2  = st_numscalar("e(F)")
    mata max(abs(b2  :- b1))
    mata max(abs(se2 :- se1))
    mata max(abs(F2  :- F1))
    mata errors = errors \ (++errow, max(abs(b2  :- b1)), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

    manyiv y (x = z1 z2) i.fe i.fe2, cluster(c)
    mata b1  = st_matrix("e(b)")
    mata se1 = st_matrix("e(se)")
    mata F1  = st_numscalar("e(F)")
    manyiv y (x = z1 z2), absorb(fe fe2) cluster(c)
    mata b2  = st_matrix("e(b)")
    mata se2 = st_matrix("e(se)")
    mata F2  = st_numscalar("e(F)")
    mata max(abs(b2  :- b1))
    mata max(abs(se2 :- se1))
    mata max(abs(F2  :- F1))
    mata errors = errors \ (++errow, max(abs(b2  :- b1)), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

    manyiv y (x = z1 z2 i.iv i.iv2), cluster(c)
    mata b1  = st_matrix("e(b)")
    mata se1 = st_matrix("e(se)")
    mata F1  = st_numscalar("e(F)")
    manyiv y (x = z1 z2), absorbiv(iv iv2) cluster(c)
    mata b2  = st_matrix("e(b)")
    mata se2 = st_matrix("e(se)")
    mata F2  = st_numscalar("e(F)")
    mata max(abs(b2  :- b1))
    mata max(abs(se2 :- se1))
    mata max(abs(F2  :- F1))
    mata errors = errors \ (++errow, max(abs(b2  :- b1)), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

    manyiv y (x = z1 z2 i.iv i.iv2) i.fe i.fe2, cluster(c)
    mata b1  = st_matrix("e(b)")
    mata se1 = st_matrix("e(se)")
    mata F1  = st_numscalar("e(F)")
    manyiv y (x = z1 z2), absorbiv(iv iv2) absorb(fe fe2) cluster(c)
    mata b2  = st_matrix("e(b)")
    mata se2 = st_matrix("e(se)")
    mata F2  = st_numscalar("e(F)")
    mata max(abs(editvalue(b2, 0, .) :- b1))
    mata max(abs(se2 :- se1))
    mata max(abs(F2  :- F1))
    mata errors = errors \ (++errow, max(abs(editvalue(b2, 0, .) :- b1)), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

    manyiv y (x = i.iv), cluster(c)
    mata b1  = st_matrix("e(b)")
    mata se1 = st_matrix("e(se)")
    mata F1  = st_numscalar("e(F)")
    manyiv y (x = .), absorbiv(iv) cluster(c)
    mata b2  = st_matrix("e(b)")
    mata se2 = st_matrix("e(se)")
    mata F2  = st_numscalar("e(F)")
    mata max(abs(b2  :- b1))
    mata max(abs(se2 :- se1))
    mata max(abs(F2  :- F1))
    mata errors = errors \ (++errow, max(abs(b2  :- b1)), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

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
    mata errors = errors \ (++errow, max(abs(b2  :- b1)), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

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
    mata errors = errors \ (++errow, max(abs(b2  :- b1)), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

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
    mata errors = errors \ (++errow, max(abs(b2  :- b1)), max(abs(se2 :- se1)), max(abs(F2  :- F1)))
    mata errors = errors \ (++errow, max(abs(b3  :- b1)), max(abs(se3 :- se1)), max(abs(F3  :- F1)))
    mata errors = errors \ (++errow, max(abs(b4  :- b1)), max(abs(se4 :- se1)), max(abs(F4  :- F1)))

    * Singletons
    replace fe = 99 in 10/11
    replace iv = 98 in 10
    replace fe = 5.5 in 100/101
    replace iv = 5.5 in 100
    replace fe = -9999 in 200
    replace iv = -9998 in 300
    replace fe = 1234 if iv == 3 & _n == 95

    cap drop feid
    cap drop ivid
    cap drop nfeid
    cap drop nivid
    egen feid = group(fe)
    egen ivid = group(iv)
    egen nfeid = count(1), by(feid)
    egen nivid = count(1), by(ivid)
    levelsof feid if nfeid == 1, loc(ofeid)
    levelsof ivid if nivid == 1, loc(oivid)

    ivregress 2sls y (x = z1 z2 io(`oivid').ivid) w io(`ofeid').feid, cluster(c) small
    manyiv y (x = z1 z2 io(`oivid').ivid) w io(`ofeid').feid, cluster(c)
    mata b1  = st_matrix("e(b)")
    mata se1 = st_matrix("e(se)")
    mata F1  = st_numscalar("e(F)")
    manyiv y (x = z1 z2) w , absorbiv(iv) absorb(fe) cluster(c) skipsingletons _plugin_skip
    * manyiv y (x = z1 z2) w , absorbiv(iv) absorb(fe) cluster(c) skipsingletons
    mata b2  = st_matrix("e(b)")
    mata se2 = st_matrix("e(se)")
    mata F2  = st_numscalar("e(F)")
    mata max(abs(editvalue(b2, 0, .) :- editvalue(b1, 0, .)))
    mata (any((b2, b1) :== 0))? b2 \ b1: .
    mata max(abs(se2 :- se1))
    mata max(abs(F2  :- F1))
    mata errors = errors \ (++errow, max(abs(editvalue(b2, 0, .) :- editvalue(b1, 0, .))), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

    manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c)
    manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c) _plugin_skip
    manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c) skipsingletons
    manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c) skipsingletons _plugin_skip
    manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c) keepsingletons
    manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c) keepsingletons _plugin_skip
    manyiv y (x = z1 z2) w if (fe != 8) | (_n == 1), absorbiv(iv) absorb(fe) cluster(c)
    manyiv y (x = z1 z2) w if (fe != 8) | (_n == 1), absorbiv(iv) absorb(fe) cluster(c) _plugin_skip
    manyiv y (x = z1 z2) w if (fe != 8) | (_n == 1), absorbiv(iv) absorb(fe) cluster(c) skipsingletons
    manyiv y (x = z1 z2) w if (fe != 8) | (_n == 1), absorbiv(iv) absorb(fe) cluster(c) skipsingletons _plugin_skip
    manyiv y (x = z1 z2) w if (fe != 8) | (_n == 1), absorbiv(iv) absorb(fe) cluster(c) keepsingletons
    manyiv y (x = z1 z2) w if (fe != 8) | (_n == 1), absorbiv(iv) absorb(fe) cluster(c) keepsingletons _plugin_skip

    * Redundant
    replace iv = 10.5 if fe == 7
    cap drop feid
    cap drop ivid
    cap drop nfeid
    cap drop nivid
    egen feid = group(fe)
    egen ivid = group(iv)
    egen nfeid = count(1), by(feid)
    egen nivid = count(1), by(ivid)
    levelsof feid if nfeid == 1, loc(ofeid)
    levelsof ivid if nivid == 1, loc(oivid)

    ivregress 2sls y (x = z1 z2 io(`oivid').ivid) w io(`ofeid').feid, cluster(c) small
    manyiv y (x = z1 z2 io(`oivid').ivid) w io(`ofeid').feid, cluster(c)
    mata b1  = st_matrix("e(b)")
    mata se1 = st_matrix("e(se)")
    mata F1  = st_numscalar("e(F)")
    manyiv y (x = z1 z2) w , absorbiv(iv) absorb(fe) cluster(c) skipsingletons _plugin_skip
    * manyiv y (x = z1 z2) w , absorbiv(iv) absorb(fe) cluster(c) skipsingletons
    mata b2  = st_matrix("e(b)")
    mata se2 = st_matrix("e(se)")
    mata F2  = st_numscalar("e(F)")
    mata max(abs(editvalue(b2, 0, .) :- editvalue(b1, 0, .)))
    mata (any((b2, b1) :== 0))? b2 \ b1: .
    mata max(abs(se2 :- se1))
    mata max(abs(F2  :- F1))
    mata errors = errors \ (++errow, max(abs(editvalue(b2, 0, .) :- editvalue(b1, 0, .))), max(abs(se2 :- se1)), max(abs(F2  :- F1)))

    misc_checks
    mata errors
end

***********************************************************************
*                                                                     *
*                             Misc Checks                             *
*                                                                     *
***********************************************************************

capture program drop misc_checks
program misc_checks
    qui tab c,  gen(_c)
    qui tab fe, gen(_fe)
    qui tab iv, gen(_iv)
    qui tab fe2, gen(_2fe)
    qui tab iv2, gen(_2iv)

    ivregress 2sls y (x = z1 z2 i.ivid) w i.feid if fe != 1234, cluster(c) small
    manyiv y (x = z1 z2 i.ivid) w i.feid if fe != 1234, cluster(c)
    manyiv y (x = z1 z2) w if fe != 1234, absorbiv(iv) absorb(fe) cluster(c) keepsingletons
    manyiv y (x = z1 z2) w if fe != 1234, absorbiv(iv) absorb(fe) cluster(c) skipsingletons
    cap noi jive y (x = z1 z2 _iv*) w _fe* if fe != 1234, r

    ivregress 2sls y (x = z1 z2 i.ivid) w i.feid, cluster(c) small
    manyiv y (x = z1 z2 i.ivid) w i.feid, cluster(c)
    manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c) keepsingletons
    manyiv y (x = z1 z2) w, absorbiv(iv) absorb(fe) cluster(c) skipsingletons
    cap noi jive y (x = z1 z2 _iv*) w _fe*, r

    cap noi jive y (x = z1 z2) w
    cap noi jive y (x = z1 z2) w _fe*
    cap noi jive y (x = z1 z2) w, r
    cap noi jive y (x = z1 z2) w _fe*, r
    cap noi jive y (x = z1 z2 _iv*) w, r
    cap noi jive y (x = z1 z2 _iv*) w _fe*, r
    cap noi jive y (x = _iv*), r
end

main
