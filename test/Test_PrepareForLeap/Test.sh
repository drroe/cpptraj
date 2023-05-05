#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in 1qos.cpptraj.pdb leap.1qos.in \
           leap.4zzw.in 4zzw.cpptraj.pdb 1qos.cpptraj.pdb.2 \
           leap.1qos.in.2 \
           leap.log leap.ff.in \
           m.parm7 m.rst7 leap.6abz.in 6abz.cpptraj.pdb cpptraj.6abz.summ \
           n.parm7 n.rst7 leap.1opd.in 1opd.cpptraj.pdb cpptraj.1opd.summ

INPUT='-i cpptraj.in'

TESTNAME='PrepareForLeap tests'

# Just in case CPPTRAJHOME is not set, see if we can find the map file
RESMAPFILE=''
if [ -z "$CPPTRAJHOME" ] ; then
  if [ -e '../../dat/Carbohydrate_PDB_Glycam_Names.txt' ] ; then
    RESMAPFILE='resmapfile ../../dat/Carbohydrate_PDB_Glycam_Names.txt'
  fi
fi

Run1qos() {
cat > cpptraj.in <<EOF
parm 1qos.pdb
loadcrd 1qos.pdb name MyCrd

prepareforleap \
  crdset MyCrd \
  name Final \
  out leap.1qos.in \
  leapunitname m \
  pdbout 1qos.cpptraj.pdb \
  nowat noh $RESMAPFILE
EOF
RunCpptraj "Prepare PDB 1qos for LEaP"
DoTest leap.1qos.in.save leap.1qos.in
DoTest 1qos.cpptraj.pdb.save 1qos.cpptraj.pdb
}

Run1qosGlycam() {
cat > cpptraj.in <<EOF
parm 1qos.cpptraj.pdb.save
loadcrd 1qos.cpptraj.pdb.save name MyCrd
prepareforleap \
  hasglycam \
  crdset MyCrd \
  name Final \
  cysmask :CYX@SG \
  out leap.1qos.in.2 \
  leapunitname m \
  pdbout 1qos.cpptraj.pdb.2 \
  nowat noh $RESMAPFILE
EOF
RunCpptraj "Prepare PDB with existing glycam names (1qos)."
DoTest 1qos.cpptraj.pdb.save 1qos.cpptraj.pdb.2
}

Run4zzw() {
cat > cpptraj.in <<EOF
parm 4zzw.pdb
loadcrd 4zzw.pdb name MyCrd

prepareforleap \
  crdset MyCrd \
  name Final \
  out leap.4zzw.in \
  leapunitname m \
  pdbout 4zzw.cpptraj.pdb \
  nowat noh keepaltloc highestocc
EOF
RunCpptraj "Prepare PDB 4zzw for LEaP"
DoTest leap.4zzw.in.save leap.4zzw.in
DoTest 4zzw.cpptraj.pdb.save 4zzw.cpptraj.pdb
}

Run6abz() {
  UNITNAME='Prepare PDB 6abz for LEaP with protonation state calculation'
  CheckFor inpath tleap
  if [ $? -eq 0 ] ; then
    echo "source leaprc.protein.ff14SB" > leap.ff.in
    cat > cpptraj.in <<EOF
parm 6ABZ.pdb
loadcrd 6ABZ.pdb name MyCrd

prepareforleap \
  crdset MyCrd \
  name Final \
  out leap.6abz.in \
  leapunitname m runleap leap.ff.in \
  pdbout 6abz.cpptraj.pdb \
  nowat noh keepaltloc highestocc stripmask :PGO,NA,CL \
  doprot summ cpptraj.6abz.summ
list
EOF
    RunCpptraj "$UNITNAME"
    DoTest 6abz.cpptraj.pdb.save 6abz.cpptraj.pdb
  fi
}

Run1opd() {
  UNITNAME='Prepare PDB 1opd for LEaP with protonation state calculation'
  CheckFor inpath tleap
  if [ $? -eq 0 ] ; then
    cat > leap.ff.in <<EOF
source leaprc.protein.ff19SB
set default PBradii mbondi2
EOF
    cat > cpptraj.in <<EOF
#parm 1OPD.pdb
#loadcrd 1OPD.pdb name MyCrd
parm H++_1opd_noReduce2/tmp.pqr relaxfmt pqr
loadcrd H++_1opd_noReduce2/tmp.pqr name MyCrd

prepareforleap \
  crdset MyCrd \
  name Final \
  out leap.1opd.in \
  leapunitname n runleap leap.ff.in \
  pdbout 1opd.cpptraj.pdb \
  nowat keepaltloc highestocc stripmask :SO4 \
  doprot summ cpptraj.1opd.summ
list
EOF
    RunCpptraj "$UNITNAME"
  fi
}

#Run6abz
Run1opd
 
EndTest
