#!/bin/bash

. ../MasterTest.sh

TESTNAME='MCTI tests'

CleanFiles cpptraj.in mc.log mc.pkout

INPUT='-i cpptraj.in'

MCTI() {
  UNITNAME='MCTI test'
  cat > cpptraj.in <<EOF
random setdefault drand48
mcti nmcsteps 1000 startph 0 stopph 12 phincr 0.2 \
     iseed 1 mclog mc.log mcprefix ../Test_MEAD/tz2 \
     mcmode reduced pkout mc.pkout
EOF
  RunCpptraj "$UNITNAME"
  DoTest mc.pkout.save mc.pkout
}

MCTI

EndTest
