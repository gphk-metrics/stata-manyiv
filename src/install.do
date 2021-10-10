cap noi net uninstall manyiv
mata: mata clear
mata: mata set matastrict on
mata: mata set mataoptimize on
cap noi erase src/lmanyiv.mlib
do src/manyiv_absorb.mata
do src/manyiv_internals_m.mata
mata: mata mlib create lmanyiv, dir("src") replace
mata: mata mlib add lmanyiv ManyIV*() sf_helper*(), dir("src") complete
net install manyiv, from(`c(pwd)'/src)
