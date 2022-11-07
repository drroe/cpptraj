#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o 

#UNITSOURCES='NameType.cpp CpptrajStdio.cpp'

CreateMakefile

RunMake "Top/Param classes unit test."

EndTest
