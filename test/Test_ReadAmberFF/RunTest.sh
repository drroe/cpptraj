#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in tip3pf.parm.dat toyrna.parm.dat

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
readdata frcmod.tip3pf as frcmod name PARM
writedata tip3pf.parm.dat PARM

readdata toyrna.dat as amberff name TOYRNA
writedata toyrna.parm.dat TOYRNA
EOF
RunCpptraj "Test read Amber FF and force modification files."
DoTest tip3pf.parm.dat.save tip3pf.parm.dat
DoTest toyrna.parm.dat.save toyrna.parm.dat

EndTest
