ManyIV
======

Various instrumental variables regressions (OLS, TSLS, LIML, MBTSLS, JIVE, UJIVE, RTSLS)

`version 0.4.0 21Oct2021` | [Installation](#installation) | [Usage](#usage) | [Examples](#examples) | [Compiling](#compiling)

### Installation

From the command line:

```
git clone git@github.com:mcaceresb/stata-manyiv
```

(or download the code manually and unzip). From Stata:

```
cap noi net uninstall manyiv
net install manyiv, from(`c(pwd)'/stata-manyiv/src)
```

If the repo were public, this can be installed directly from Stata:

```
local github "https://raw.githubusercontent.com"
net install manyiv, from(`github'/mcaceresb/stata-manyiv/master/src/)
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

### Compiling

To compute the jive and ujive estimators with multiple absorb levels, `manyiv` uses a plugin compiled with the [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) library. To compile the plugin from source, install Eigen3 (using your system's package manager or from [their website](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)) and run:

```
git clone git@github.com:mcaceresb/stata-manyiv
cd stata-manyiv
make all EIGEN=/path/to/eigen3
```

Then open Stata from the current directory and run:

```
cap ado uninstall manyiv
net install manyiv, from(`c(pwd)'/src)
```
