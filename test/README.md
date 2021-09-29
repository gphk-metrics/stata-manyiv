ManyIV
======

Various instrumental variables regressions (OLS, TSLS, LIML, MBTSLS, JIVE, UJIVE, RTSLS)

`version 0.1.0 27Sep2021` | [Installation](#installation) | [Usage](#usage) | [Examples](#examples)

### Installation

From the command line

```
git clone https://github.com/mcaceresb/stata-manyiv
```

From Stata

```
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
gen c  = mod(_n, 10)
gen w  = rnormal()
gen x  = 1 + z1 - z2 + u
gen y  = 1 + x + w + e

manyiv y (x = z1 z2) w
manyiv y (x = z1 z2) w, cluster(c)
```
