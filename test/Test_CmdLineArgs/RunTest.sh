#!/bin/bash

. ../MasterTest.sh
# Test that command line arguments are working properly.
# Note that these tests are a little different in that
# in some cases we want to compare the output. Therefore,
# the test output file will be directed to a new save
# each test.
TESTNAME='Command line arguments tests'

CleanFiles out.crd test.traj.out test.tl.out test.ms.out test.mr.out \
           test.mask.out test.resmask.out

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
# NOTE: mmpbsa.py relies on this.
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

# Test -mr
INPUT='-p ../tz2.parm7 -mr :TRP -o test.mr.out'
RunCpptraj "Test command-line select residue numbers by mask"
MoveTestOut test.mr.out
DoTest test.mr.out.save test.mr.out

# Test --mask
INPUT='-p ../tz2.parm7 --mask @CA -o test.mask.out'
RunCpptraj "Test command-line select atoms by mask"
MoveTestOut test.mask.out
DoTest test.mask.out.save test.mask.out

# Test --resmask
INPUT='-p ../tz2.parm7 --resmask :TRP -o test.resmask.out'
RunCpptraj "Test command-line select residues by mask"
MoveTestOut test.resmask.out
DoTest test.resmask.out.save test.resmask.out

EndTest
