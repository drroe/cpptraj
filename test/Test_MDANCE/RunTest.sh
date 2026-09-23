#!/bin/bash

. ../MasterTest.sh

TESTNAME='MDANCE tests'

INPUT='mdance.in'

CleanFiles mdance.in cnumvtime.dat noh.cnumvtime.dat cpptraj.json

UNITNAME='MDANCE CSV data test'
cat > mdance.in <<EOF
readdata sim.csv as coords name MyCrd
mdance crdset MyCrd clusters 10 nosort json cpptraj.json kseed 1
EOF
RunCpptraj "$UNITNAME"
DoTest result.json cpptraj.json

UNITNAME='MDANCE Kmeans'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > mdance.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
align first
trajout tz2.aligned.nc
createcrd MyCrd

mdance crdset MyCrd clusters 5 out cnumvtime.dat kseed 1
EOF
  RunCpptraj "$UNITNAME"
  DoTest cnumvtime.dat.save cnumvtime.dat
fi

UNITNAME='MDANCE Kmeans with atom selection'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > mdance.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
align first
trajout tz2.aligned.nc
createcrd MyCrd

mdance crdset MyCrd clusters 5 out noh.cnumvtime.dat mask !@H= kseed 1
EOF
  RunCpptraj "$UNITNAME"
  DoTest noh.cnumvtime.dat.save noh.cnumvtime.dat
fi

EndTest

