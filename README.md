ManyIV
======

Various instrumental variables regressions (OLS, TSLS, LIML, MBTSLS,
JIVE, UJIVE, RTSLS) with absorbed instruments and controls.

`version 0.6.1 26May2022` | [Installation](#installation) | [Usage](#usage) | [Examples](#examples) | [Compiling](#compiling)

### Installation

From the command line:

```
git clone git@github.com:gphk-metrics/stata-manyiv
```

(or download the code manually and unzip). From Stata:

```
cap noi net uninstall manyiv
net install manyiv, from(`c(pwd)'/stata-manyiv)
```

(Change `stata-manyiv` if you download the package to a different
folder; e.g. `stata-manyiv-main`.) To check that `manyiv` was
installed correctly, please run

```
manyiv _plugin_check
```

If you get an error, `manyiv` may experience issues in some cases.
In particular, without the plugin the number of absorb variables is
capped at 2, and when the absorb variables have many levels the
function's performance will be very slow. If the plugin fails to
load, please see the [compiling](#compiling) section below.

NB: if the repo were public, this would be installed from Stata via
```
local github "https://raw.githubusercontent.com"
net install manyiv, from(`github'/gphk-metrics/stata-manyiv/master/)
```

### Usage

```
manyiv depvar (endogenous = instrument) [exogenous], options
help manyiv
```

### Examples

```stata
clear
set seed 1729
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

manyiv y (x = z1 z2) w
manyiv y (x = z1 z2) w, cluster(c)
manyiv y (x = z1 z2) w, absorb(fe) cluster(c)
manyiv y (x = z1 z2) w, absorbiv(iv) cluster(c)
manyiv y (x = z1 z2) w, absorb(fe) absorbiv(iv) cluster(c)
manyiv y (x = .)     w, absorb(fe) absorbiv(iv) cluster(c)
```

There are instances where jive/ujive will fail because they are
leave-one-out estimators and a covariate or fixed effect group will be
collinear with a given observation. If dealing with the issue manually
is otherwise prohibitive, the `forcejive` option will prompt `manyiv` to
drop observations and/or covariates until jive/ujive can run.

```stata
clear
input z1 z2 fe iv
-1 -1 0 0
1  1 0 0
3  3 1 0
-3 -3 1 1
0  0 2 1
0  0 2 1
1 2 3 2
2 3 3 2
end
gen z3 = _n^2
gen x  = z1 + z2 + z3 + rnormal()
gen y  = x + rnormal()
manyiv y (x = z1 z2), absorb(fe) absorbiv(iv)
manyiv y (x = z1 z2), absorb(fe) absorbiv(iv) forcejive
```

### Compiling

To compute the jive and ujive estimators with multiple absorb levels, `manyiv` uses a plugin compiled with the [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) library. Pre-compiled binaries are provided in [this reposirtory](src/build). To re-compile the plugin from source, download Eigen3 using your system's package manager (e.g. `brew`, `apt`, etc.) or from [their website](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) and run:

```
git clone git@github.com:gphk-metrics/stata-manyiv
cd stata-manyiv
make all EIGEN=/path/to/eigen3
```

(NB: On Windows, this needs to be compiled using MinGW.) If your system's
libraries are not compatible with the latest version of Eigen3, you can
still compile `manyiv` using an older version. For example,

```
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
git checkout 3.3.9
```

Replace 3.3.9 with the version you want to try. (Run `git tag`
to see all available versions.) Then you can run `make` in the
`stata-manyiv` folder, setting `EIGEN` to the folder you just downloaded
from `gitlab`. Once the plugin is compiled, open Stata from the
`stata-manyiv` directory and run:

```
cap ado uninstall manyiv
cap net uninstall manyiv
net install manyiv, from(`c(pwd)')
```
