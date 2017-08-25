#!/bin/bash

. ../MasterTest.sh

CleanFiles rama.in rama.dat total.dat dih.dat

INPUT='-i rama.in'

cat > rama.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc
rama DPDP out rama.dat totalout total.dat :6-7
multidihedral DIH phi psi resrange 6-7 out dih.dat
EOF
RunCpptraj "Ramachandran test."

EndTest
exit 0
