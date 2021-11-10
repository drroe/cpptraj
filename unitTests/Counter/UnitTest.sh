#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o a.out Counter.o Counter_Regular.o Counter_Array.o CpptrajStdio.o StringRoutines.o

UNITSOURCES='Counter.cpp Counter_Regular.cpp Counter_Array.cpp CpptrajStdio.cpp StringRoutines.cpp'

CreateMakefile

RunMake "Counter classes unit test."

EndTest
