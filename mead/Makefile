# Makefile for MEAD bundled with cpptraj.
include ../config.h
DEL_FILE      = /bin/rm -f
AR            = ar cqs
TARGET        = libmead.a

SOURCES       = AnalyMaker.cc AnalySlab.cc AnalySphere.cc Angle.cc Atom.cc AtomChargeSet.cc AtomID.cc AtomSet.cc AvsScalarField.cc AvsScalarField_input.cc AvsSF_read_header.cc Bigmem.cc ChargeCubeRep.cc ChargeDist.cc CLShell.cc CubeLatSpec.cc Debye.cc DielByAtoms.cc DielCubeRep.cc DielectricEnvironment.cc DielectricSlab.cc DielectricSphere.cc DielMembAtoms.cc ElectrolyteByAtoms.cc ElectrolyteEnvironment.cc Elstat_list_init.cc ElstatMaker.cc ElstatPot.cc ElstatPotCombination.cc ElySphere.cc FDChargeIterator.cc FDElstatMaker.cc FDGridLevel.cc FinDiffElstatPot.cc FinDiffMethod.cc globals.cc ManyPointCharge.cc MomentAnalysis.cc multPotxCharge.cc OnePointCharge.cc Pair.cc PairGenBorn.cc PhysCond.cc Polynomial.cc Potat.cc Sausage.cc SAVanal_calc.cc SAVcalc_cuberep.cc SAVtag_points.cc Shell.cc SolvAccVol.cc sor.cc SphericalHarmonic.cc TPShell.cc UniformDielectric.cc UniformElectrolyte.cc

OBJECTS=$(SOURCES:.cc=.o)

.cc.o:
	$(VB)echo [MEAD] CXX $<
	$(VB)$(CXX) $(DIRECTIVES) -DUSE_EXCEPTIONS -I. $(CXXFLAGS) -c -o $@ $<

all: $(TARGET)

$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

uninstall: clean

clean:
	$(DEL_FILE) *.o $(TARGET)

# Dependencies
depend: ../src/findDepend
	../src/findDepend $(SOURCES) > meaddepend

include meaddepend
