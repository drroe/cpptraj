#!/bin/bash

. ../MasterTest.sh

TESTNAME='MDANCE tests'
Requires netcdf

INPUT='mdance.in'

CleanFiles mdance.in cnumvtime.dat noh.cnumvtime.dat

UNITNAME='MDANCE Kmeans'
cat > mdance.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
align first
trajout tz2.aligned.nc
createcrd MyCrd

mdance crdset MyCrd clusters 5 out cnumvtime.dat
EOF
RunCpptraj "$UNITNAME"
DoTest cnumvtime.dat.save cnumvtime.dat

UNITNAME='MDANCE Kmeans with atom selection'
cat > mdance.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
align first
trajout tz2.aligned.nc
createcrd MyCrd

mdance crdset MyCrd clusters 5 out noh.cnumvtime.dat mask !@H=
EOF
RunCpptraj "$UNITNAME"
DoTest noh.cnumvtime.dat.save noh.cnumvtime.dat

EndTest

