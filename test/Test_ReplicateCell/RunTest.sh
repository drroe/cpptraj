#!/bin/bash

. ../MasterTest.sh

CleanFiles replicate.in cell.mol2 MyCoords.mol2

TESTNAME='Replicate cell test'
Requires netcdf maxthreads 1

INPUT='-i replicate.in'

cat > replicate.in <<EOF
parm ../Test_SymmRmsd/TYR.parm7
trajin ../Test_SymmRmsd/TYR.nc 1 1
box x 20 y 20 z 20 beta 90
replicatecell out cell.mol2 dir 100 dir 0-10 dir 001 name MyCoords
run
crdout MyCoords MyCoords.mol2
EOF
RunCpptraj "$TESTNAME"
DoTest cell.mol2.save cell.mol2
DoTest cell.mol2.save MyCoords.mol2

EndTest
exit 0
