shell make all
cap noi ado uninstall manyiv
mata: mata clear
mata: mata set matastrict on
mata: mata set mataoptimize on
cap noi erase src/build/lmanyiv.mlib
do src/mata/manyiv_absorb.mata
do src/mata/manyiv_internals.mata
mata: mata mlib create lmanyiv, dir("src/build") replace
mata: mata mlib add lmanyiv ManyIV*() sf_helper*(), dir("src/build") complete
net install manyiv, from(`c(pwd)') replace
