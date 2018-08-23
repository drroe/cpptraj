#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.in rms.save gzrms.dat bzrms.dat

TESTNAME='Compressed trajectory read tests'
INPUT='-i cpptraj.in'
Requires maxthreads 10 zlib

# Test 1
UNITNAME='Uncompressed trajectory read'
cat > cpptraj.in <<EOF
parm ../Test_TrajinOffset/ala2.99sb.mbondi2.parm7
trajin ../Test_TrajinOffset/rem.crd.000
rms first !@H= mass out rms.save
EOF
RunCpptraj "$UNITNAME."

UNITNAME='Gzipped trajectory read'
cat > cpptraj.in <<EOF
parm ../Test_TrajinOffset/ala2.99sb.mbondi2.parm7
trajin ../Test_TrajinOffset/rem.crd.000.gz
rms first !@H= mass out gzrms.dat
EOF
RunCpptraj "$UNITNAME."
DoTest rms.save gzrms.dat 

UNITNAME='Bzip2ed trajectory read'
CheckFor bzlib
if [ $? -eq 0 ] ; then
  cat > cpptraj.in <<EOF
parm ../Test_TrajinOffset/ala2.99sb.mbondi2.parm7
trajin ../Test_TrajinOffset/rem.crd.000.bz2
rms first !@H= mass out bzrms.dat
EOF
  RunCpptraj "$UNITNAME."
  DoTest rms.save bzrms.dat 
fi

EndTest

exit 0
