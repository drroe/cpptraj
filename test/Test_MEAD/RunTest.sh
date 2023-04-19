#!/bin/bash

. ../MasterTest.sh

TESTNAME='MEAD tests'

CleanFiles cpptraj.in potential.dat test.sphere.pqr solvate.dat bounds.dat MyGrid.dx \
           tz2.dat tz2.ssi.dat tz2.pkint tz2.summ tz2.g

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

MultiFlex() {
  UNITNAME='MEAD multiflex test'
  cat > cpptraj.in <<EOF
parm tz2.pqr pqr
loadcrd tz2.pqr name TZ2
mead multiflex crdset TZ2 \
     ionicstr 0.15 solrad 1.4 epsin 4 epssol 80 \
     sites tz2.sites \
     out tz2.dat ssiout tz2.ssi.dat name TZ2 \
     pkint tz2.pkint summ tz2.summ gfile tz2.g
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.pkint.save tz2.pkint -a 0.0001
  DoTest tz2.summ.save tz2.summ -a 0.0001
  DoTest tz2.g.save tz2.g -a 0.000000001
  DoTest tz2.dat.save tz2.dat -a 0.0002
  DoTest tz2.ssi.dat.save tz2.ssi.dat
}

Potential
Solvate
MultiFlex

EndTest
