#!/bin/bash

. ../MasterTest.sh

CleanFiles proj.in project.dat

TESTNAME='PC projection test'
Requires netcdf

CPPTRAJ_TOP=../tz2.parm7
TRJ=../tz2.nc

CPPTRAJ_INPUT="proj.in"

cat > proj.in <<EOF
trajin $TRJ
readdata tz2.modes.mwcovar.dat
projection modes tz2.modes.mwcovar.dat @CA out project.dat
EOF
RunCpptraj "Projection"
DoTest project.dat.save project.dat

EndTest

exit 0
