#!/bin/bash

. ../MasterTest.sh

CleanFiles rama.in rama.dat total.dat dih.dat sum.agr

INPUT='-i rama.in'

cat > rama.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc
#rama DPDP out rama.dat totalout total.dat :6-7
rama usechars sumout sum.agr DPDP out rama.dat totalout total.dat
multidihedral DIH phi psi resrange 6-7 out rama.dat
EOF
RunCpptraj "Ramachandran test."
DoTest sum.agr.save sum.agr

EndTest
exit 0
