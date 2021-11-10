#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o a.out Counter.o Counter_Regular.o Counter_Array.o CpptrajStdio.o StringRoutines.o

UNITSOURCES='ArgList.cpp Counter.cpp Counter_Regular.cpp Counter_Array.cpp CpptrajStdio.cpp StringRoutines.cpp TrajFrameCounter.cpp'

CreateMakefile

RunMake "TrajFrameIndex template class unit test."

EndTest
