#!/bin/bash

. ../MasterTest.sh

TESTNAME='MEAD tests'

CleanFiles cpptraj.in potential.dat test.sphere.pqr solvate.dat bounds.dat MyGrid.dx

INPUT='-i cpptraj.in'

Potential() {
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
}

Solvate() {
  UNITNAME='MEAD solvate test'
  cat > cpptraj.in <<EOF
parm sphere.pqr pqr
loadcrd sphere.pqr name SPHERE
crdaction SPHERE bounds out bounds.dat name MyGrid dx 1.0 offset 5
crdout SPHERE test.sphere.pqr pdb dumpq
mead \
  ogm 41,1.0 \
  ogm 41,0.25 \
  crdset SPHERE \
  out solvate.dat \
  name EPS1 \
  verbose 2 \
  solvate epsin 1 rxngrid MyGrid
writedata MyGrid.dx MyGrid
mead \
  ogm 41,1.0 \
  ogm 41,0.25 \
  crdset SPHERE \
  out solvate.dat \
  name EPS4 \
  verbose 0 \
  solvate epsin 4
EOF
  RunCpptraj "$UNITNAME"
  DoTest solvate.dat.save solvate.dat
  DoTest MyGrid.dx.save MyGrid.dx
}

Potential
Solvate

EndTest
