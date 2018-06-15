#!/bin/bash

. ../MasterTest.sh

CleanFiles tpd.in tpd.dat tpd.gnu tpd.*.gnu

INPUT='-i tpd.in'

cat > tpd.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
twopartdiff TPD out tpd.gnu
EOF
RunCpptraj "Two particle diffusion test"
DoTest tpd.drr.gnu.save tpd.drr.gnu
DoTest tpd.dtt.gnu.save tpd.dtt.gnu

EndTest
exit 0
