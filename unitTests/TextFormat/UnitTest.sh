#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o TextFormat.o a.out StringRoutines.o CpptrajStdio.o

UNITSOURCES='TextFormat.cpp StringRoutines.cpp CpptrajStdio.cpp'

CreateMakefile

RunMake "TextFormat class unit test."

EndTest
