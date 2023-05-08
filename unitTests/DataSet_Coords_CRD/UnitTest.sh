#!/bin/bash

. ../UnitMaster.sh

CleanFiles a.out Makefile *.o 

UNITSOURCES='DataSet_Coords_CRD.cpp DataSet_Coords.cpp Topology.cpp Frame.cpp CpptrajStdio.cpp Box.cpp Vec3.cpp NameType.cpp Matrix_3x3.cpp Atom.cpp Residue.cpp AtomMask.cpp CharMask.cpp CoordinateInfo.cpp StringRoutines.cpp MaskToken.cpp DistRoutines.cpp FileName.cpp Range.cpp CompactFrameArray.cpp ArgList.cpp DataSet.cpp TextFormat.cpp MetaData.cpp ParameterSet.cpp CpptrajFile.cpp FileIO_Std.cpp FileIO_Bzip2.cpp FileIO_Gzip.cpp'

CreateMakefile

RunMake "DataSet_Coords_CRD unit test."

EndTest
