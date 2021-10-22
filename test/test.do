***********************************************************************
*                              R Package                              *
***********************************************************************

library(ManyIV)
formula <- lwage ~ education + as.factor(yob) | as.factor(qob) * as.factor(yob)
r1 <- IVreg(formula, data = ak80, inference = c("standard", "re", "il", "lil"))
r2 <- IVreg(formula, data = ak80, inference = c("md"))
# write_dta(ak80, "/tmp/tmp.dta")

***********************************************************************
*                             Stata Test                              *
***********************************************************************

set linesize 112
cd ~/projects/ra/hull/manyiv/stata-manyiv/src
mata mata clear
qui do ../src/mata/manyiv_absorb.mata
qui do ../src/mata/manyiv_internals.mata
qui do ../src/ado/manyiv.ado
use "/tmp/tmp.dta", clear
manyiv lwage (education = i.yob##i.qob) i.yob, save(r0)

***********************************************************************
*                            Malab Export                             *
***********************************************************************

mata
fh = fopen("/tmp/tmpY", "rw")
C = bufio()
fbufput(C, fh, "%8z", `Y'')
fclose(fh)

fh = fopen("/tmp/tmpX", "rw")
C = bufio()
fbufput(C, fh, "%8z", `X'')
fclose(fh)

fh = fopen("/tmp/tmpZ", "rw")
C = bufio()
fbufput(C, fh, "%8z", `Z'')
fclose(fh)

fh = fopen("/tmp/tmpW", "rw")
C = bufio()
fbufput(C, fh, "%8z", `W'')
fclose(fh)

fh = fopen("/tmp/tmpD", "rw")
C = bufio()
fbufput(C, fh, "%8z", rows(`Y'))
fbufput(C, fh, "%8z", cols(`Y'))
fbufput(C, fh, "%8z", rows(`X'))
fbufput(C, fh, "%8z", cols(`X'))
fbufput(C, fh, "%8z", rows(`Z'))
fbufput(C, fh, "%8z", cols(`Z'))
fbufput(C, fh, "%8z", rows(`W'))
fbufput(C, fh, "%8z", cols(`W'))
fclose(fh)
end

***********************************************************************
*                             Matlab Test                             *
***********************************************************************

fileID = fopen('/tmp/tmpD', 'r');
dims = fread(fileID, 'double');
fclose(fileID)

fileID = fopen('/tmp/tmpY', 'r');
Y = fread(fileID, dims(1:2)', 'double');
fclose(fileID)

fileID = fopen('/tmp/tmpX', 'r');
X = fread(fileID, dims(3:4)', 'double');
fclose(fileID)

fileID = fopen('/tmp/tmpZ', 'r');
Z = fread(fileID, dims(5:6)', 'double');
fclose(fileID)

fileID = fopen('/tmp/tmpW', 'r');
W = fread(fileID, dims(5:6)', 'double');
fclose(fileID)

cd ~/projects/ra/hull/manyiv/kolesarm-ivreg
[b, se, stats] = ivreg(Y, X, Z, W)
