#!/bin/bash

. ../MasterTest.sh
CleanFiles corr.in corr.dat cross.dat
TESTNAME='Correlation test'
Requires netcdf
INPUT="-i corr.in"
cat > corr.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

distance d1 :2 :12
radgyr rg :1-3 mass nomax 

multidihedral TZ180 phi psi resrange 7
multidihedral TZ360 phi psi resrange 7 range360

corr d1 d1 out corr.dat 
corr d1 rg out cross.dat 
run
writedata tz2.phipsi.agr TZ180[*] TZ360[*]
set TYPE = psi
runanalysis corr TZ180[\$TYPE]:7 out dih.corr.agr name \$TYPE.7.180
runanalysis corr TZ360[\$TYPE]:7 out dih.corr.agr name \$TYPE.7.360
runanalysis corr TZ180[\$TYPE]:7 out dih.corr.agr name c\$TYPE.7.180 cosphi
runanalysis corr TZ360[\$TYPE]:7 out dih.corr.agr name c\$TYPE.7.360 cosphi

runanalysis corr TZ180[phi]:7 TZ180[psi]:7 out dih.cross.agr name phipsi.7.180
runanalysis corr TZ360[phi]:7 TZ360[psi]:7 out dih.cross.agr name phipsi.7.360

dataset wrap TZ180[*] TZ360[*]
writedata tz2.phipsi.wrap.agr TZ180[*] TZ360[*]

runanalysis corr TZ180[\$TYPE]:7 out dih.corr.agr name w\$TYPE.7.180
runanalysis corr TZ360[\$TYPE]:7 out dih.corr.agr name w\$TYPE.7.360

dataset mode distance  TZ*
runanalysis corr TZ180[phi]:7 TZ180[psi]:7 out dih.cross.agr name wphipsi.7.180
runanalysis corr TZ360[phi]:7 TZ360[psi]:7 out dih.cross.agr name wphipsi.7.360

writedata dih.cross.dat phipsi.7.180 phipsi.7.360 wphipsi.7.180 wphipsi.7.360
list
EOF
RunCpptraj "$TESTNAME"
DoTest corr.dat.save corr.dat
DoTest cross.dat.save cross.dat

cat > corr.in <<EOF
for FILE in tz2.phipsi.agr,tz2.phipsi.wrap.agr NAME in NoWrap,Wrap
  readdata \$FILE
  list
  runanalysis corr \$FILE:0 out \$NAME.corr.agr name \$NAME.phi.7.180
  runanalysis corr \$FILE:1 out \$NAME.corr.agr name \$NAME.psi.7.180
  runanalysis corr \$FILE:2 out \$NAME.corr.agr name \$NAME.phi.7.360
  runanalysis corr \$FILE:3 out \$NAME.corr.agr name \$NAME.psi.7.360
done
EOF
RunCpptraj "Corr test, read dihedral data"

EndTest

exit 0
