#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in tip3pf.parm.dat

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
readdata frcmod.tip3pf as frcmod name PARM
writedata tip3pf.parm.dat PARM
EOF
RunCpptraj "Test read Amber force modification file."
DoTest tip3pf.parm.dat.save tip3pf.parm.dat

EndTest
