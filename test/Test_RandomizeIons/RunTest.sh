#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in random.crd

TESTNAME='randomizeions test'
Requires maxthreads 1 zlib

CPPTRAJ_INPUT="ptraj.in"
CPPTRAJ_TOP="adh206.ff10.tip3p.parm7.gz"
cat > ptraj.in <<EOF
trajin adh206.tip3p.rst7.gz
randomizeions @Na+ around :1-16 by 5.0 overlap 3.0 seed 113698
trajout random.crd title "Test"
EOF
RunCpptraj "randomizeions test"
DoTest random.crd.save random.crd

EndTest

exit 0
