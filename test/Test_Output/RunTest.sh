#!/bin/bash

. ../MasterTest.sh

CleanFiles out.in mytest.out

INPUT='-i out.in'

cat > out.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
distance d1 :1 :12
run
output to mytest.out
printdata d1
output reset
quit
EOF
RunCpptraj "Output redirection test."
DoTest mytest.out.save mytest.out

EndTest
