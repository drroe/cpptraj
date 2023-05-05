#!/bin/bash

cpptraj > cpptraj.out <<EOF
#set DIR = H++_1opd_noReduce2
set PREFIX = MyCrd

parm tmp.pqr relaxfmt pqr
atoms :1
loadcrd tmp.pqr relaxfmt name \$PREFIX

#     sites 0.15_80_4_pH6.5_1OPD.sites \

random setdefault drand48
debug 10
mead multiflex verbose 3 crdset \$PREFIX \
     ionicstr 0.15 solrad 1.4 epsin 4 epssol 80 \
     out cpptraj.\$PREFIX.dat ssiout cpptraj.ssi.dat name \$PREFIX \
     pkint cpptraj.pkint summ cpptraj.summ gfile cpptraj.g 

mcti setname \$PREFIX nmcsteps 1000 startph 0 stopph 12 phincr 0.2 \
     iseed 1 mclog cpptraj.mc.log mcmode reduced pkout cpptraj.pkout
EOF
