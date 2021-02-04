#!/bin/bash

. ../MasterTest.sh
# Test that command line arguments are working properly.
# Note that these tests are a little different in that
# in some cases we want to compare the output. Therefore,
# the test output file will be directed to a new save
# each test.
TESTNAME='Command line arguments tests'

CleanFiles out.crd test.traj.out test.tl.out

# Test that -p, -y, and -x work. Also uses -o.
#INPUT='-p ../tz2.parm7 -y ../tz2.crd -x out.crd -o test.traj.out'
#RunCpptraj "Test command-line trajectory conversion"

# Test -tl
INPUT='-p ../tz2.parm7 -y ../tz2.nc -y ../tz2.crd -y ../tz2.pdb -tl -o test.tl.out'
RunCpptraj "Test command-line trajectory length"
# For old version of cpptraj TODO deprecate
if [ ! -f 'test.tl.out' -a -f 'test.out' ] ; then
  mv test.out test.tl.out
fi
DoTest test.tl.out.save test.tl.out

EndTest
