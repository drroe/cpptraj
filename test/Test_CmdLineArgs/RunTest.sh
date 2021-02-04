#!/bin/bash

. ../MasterTest.sh
# Test that command line arguments are working properly.
# Note that these tests are a little different in that
# in some cases we want to compare the output. Therefore,
# the test output file will be directed to a new save
# each test.
TESTNAME='Command line arguments tests'

CleanFiles out.crd test.traj.out test.tl.out test.ms.out

# For old version of cpptraj that did not redirect output for command-line
# arguments when -o given, this moves current test.out to target file.
# TODO deprecate
MoveTestOut() {
  TGT=$1
  if [ ! -f "$TGT" -a -f 'test.out' ] ; then
    mv test.out $TGT
  fi
}

# Test that -p, -y, and -x work. Also uses -o.
#INPUT='-p ../tz2.parm7 -y ../tz2.crd -x out.crd -o test.traj.out'
#RunCpptraj "Test command-line trajectory conversion"

# Test -tl
INPUT='-p ../tz2.parm7 -y ../tz2.nc -y ../tz2.crd -y ../tz2.pdb -tl -o test.tl.out'
RunCpptraj "Test command-line trajectory length"
# TODO deprecate
MoveTestOut test.tl.out
DoTest test.tl.out.save test.tl.out

# Test -ms
INPUT='-p ../tz2.parm7 -ms @CA -o test.ms.out'
RunCpptraj "Test command-line select atom numbers by mask"
MoveTestOut test.ms.out
DoTest test.ms.out.save test.ms.out

EndTest
