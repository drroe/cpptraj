#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o Range.o a.out CpptrajStdio.o ArgList.o StringRoutines.o

UNITSOURCES='Range.cpp CpptrajStdio.cpp ArgList.cpp StringRoutines.cpp'

CreateMakefile

RunMake "Range class unit test."

EndTest
