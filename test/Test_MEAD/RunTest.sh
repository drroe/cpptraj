#!/bin/bash

. ../MasterTest.sh

TESTNAME='MEAD tests'

CleanFiles cpptraj.in potential.dat

INPUT='-i cpptraj.in'

UNITNAME='MEAD potential test'
# The MEAD potential output for these coords is:
# 0.00249097
# 0.0203602
cat > cpptraj.in <<EOF
parm tz2.pqr pqr
loadcrd tz2.pqr name TZ2
mead \
  ogm 31,1.0 \
  crdset TZ2 \
  out potential.dat \
  verbose 2 \
  potential epsin 1 epsext 80 \
    fpt 0.0,0.0,0.0 fpt 2.0,0.0,0.0
EOF
RunCpptraj "$UNITNAME"
DoTest potential.dat.save potential.dat

UNITNAME='MEAD solvate test'
cat > cpptraj.in <<EOF
parm sphere.pqr pqr
loadcrd sphere.pqr name SPHERE
mead \
  ogm 41,1.0 \
  ogm 41,0.25 \
  verbose 2 \
  solvate epsin 1
EOF
RunCpptraj "$UNITNAME"

EndTest
