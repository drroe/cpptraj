# Test makefile for cpptraj
test: test.all 

test.general:
	@-cd Test_General && ./RunTest.sh $(OPT)
	@-cd Test_Parallel && ./RunTest.sh $(OPT)

test.strip:
	@-cd Test_Strip && ./RunTest.sh $(OPT)

test.brokentraj:
	@-cd Test_BrokenTraj && ./RunTest.sh $(OPT)

test.trajinoffset:
	@-cd Test_TrajinOffset && ./RunTest.sh $(OPT)

test.remdnetcdf:
	@-cd Test_RemdTraj && ./RunTest.sh $(OPT)

test.framerange:
	@-cd Test_trajout_framerange && ./RunTest.sh $(OPT)

test.perres:
	@-cd Test_PerResRMSD && ./RunTest.sh $(OPT)

test.dssp:
	@-cd Test_DSSP && ./RunTest.sh $(OPT)

test.center:
	@-cd Test_Center && ./RunTest.sh $(OPT)

test.hbond:
	@-cd Test_Hbond && ./RunTest.sh $(OPT)
	@-cd Test_Hbond_Split && ./RunTest.sh $(OPT)
	@-cd Test_Hbond_Matrix && ./RunTest.sh $(OPT)

test.image:
	@-cd Test_Image && ./RunTest.sh $(OPT)

test.netcdftraj:
	@-cd Test_NetcdfTraj && ./RunTest.sh $(OPT)

test.surf:
	@-cd Test_Surf && ./RunTest.sh $(OPT)

test.radgyr:
	@-cd Test_Radgyr && ./RunTest.sh $(OPT)

test.mask:
	@-cd Test_Mask && ./RunTest.sh $(OPT)

test.closest:
	@-cd Test_Closest && ./RunTest.sh $(OPT)

test.rmsd:
	@-cd Test_RMSD && ./RunTest.sh $(OPT)

test.outtraj:
	@-cd Test_Outtraj && ./RunTest.sh $(OPT)

test.datafile:
	@-cd Test_Datafile && ./RunTest.sh $(OPT)

test.avgref:
	@-cd Test_AvgRef && ./RunTest.sh $(OPT)

test.mol2:
	@-cd Test_Mol2 && ./RunTest.sh $(OPT)

test.nastruct:
	@-cd Test_NAstruct && ./RunTest.sh $(OPT)
	@-cd Test_NAstruct_8OG && ./RunTest.sh $(OPT)

test.pucker:
	@-cd Test_Pucker && ./RunTest.sh $(OPT)

test.atommap:
	@-cd Test_AtomMap && ./RunTest.sh $(OPT)

test.2drms:
	@-cd Test_2DRMS && ./RunTest.sh $(OPT)

test.average:
	@-cd Test_Average && ./RunTest.sh $(OPT)

test.radial:
	@-cd Test_Radial && ./RunTest.sh $(OPT)

test.hist:
	@-cd Test_Hist && ./RunTest.sh $(OPT)

test.ncrestart:
	@-cd Test_NCrestart && ./RunTest.sh $(OPT)

test.drmsd:
	@-cd Test_DRMSD && ./RunTest.sh $(OPT)

test.jcoupling:
	@-cd Test_Jcoupling && ./RunTest.sh $(OPT)

test.cluster:
	@-cd Test_Cluster && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_Sieve && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_RandomSieve && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_DataSet && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_DBSCAN && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_Kmeans && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_SymmRMSD && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_TrajWrites && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_Dpeaks && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_Nreps && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_AssignRefs && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_Cache && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_ReadInfo && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_CoordsAndData && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_QuatRMSD && ./RunTest.sh $(OPT)

test.ired:
	@-cd Test_IRED && ./RunTest.sh $(OPT)

test.molsurf:
	@-cd Test_Molsurf && ./RunTest.sh $(OPT)

test.charmm:
	@-cd Test_Charmm && ./RunTest.sh $(OPT)
	@-cd Test_CharmmParams && ./RunTest.sh $(OPT)

test.rotdif:
	@-cd Test_Rotdif && ./RunTest.sh $(OPT)

test.corr:
	@-cd Test_Corr && ./RunTest.sh $(OPT)

test.checkstructure:
	@-cd Test_CheckStructure && ./RunTest.sh $(OPT)

test.runavg:
	@-cd Test_Runavg && ./RunTest.sh $(OPT)

test.analrunavg:
	@-cd Test_AnalysisRunAvg && ./RunTest.sh $(OPT)

test.rmsavgcorr:
	@-cd Test_RmsAvgCorr && ./RunTest.sh $(OPT)

test.contacts:
	@-cd Test_Contacts && ./RunTest.sh $(OPT)

test.principal:
	@-cd Test_Principal && ./RunTest.sh $(OPT)

test.clusterdihedral:
	@-cd Test_ClusterDihedral && ./RunTest.sh $(OPT)

test.vector:
	@-cd Test_Vector && ./RunTest.sh $(OPT)

test.matrix:
	@-cd Test_Matrix && ./RunTest.sh $(OPT)

test.atomicfluct:
	@-cd Test_AtomicFluct && ./RunTest.sh $(OPT)

test.grid:
	@-cd Test_Grid && ./RunTest.sh $(OPT)
	@-cd Test_Grid_Rotate_Simple && ./RunTest.sh $(OPT)

test.unwrap:
	@-cd Test_Unwrap && ./RunTest.sh $(OPT)

test.project:
	@-cd Test_Project && ./RunTest.sh $(OPT)

test.analyzemodes:
	@-cd Test_Analyze_Modes && ./RunTest.sh $(OPT)

test.watershell:
	@-cd Test_Watershell && ./RunTest.sh $(OPT)

test.diffusion:
	@-cd Test_Diffusion && ./RunTest.sh $(OPT)

test.volmap:
	@-cd Test_Volmap && ./RunTest.sh $(OPT)

test.lie:
	@-cd Test_LIE && ./RunTest.sh $(OPT)

test.spam:
	@-cd Test_SPAM && ./RunTest.sh $(OPT)

test.randomizeions:
	@-cd Test_RandomizeIons && ./RunTest.sh $(OPT)

test.density:
	@-cd Test_Density && ./RunTest.sh $(OPT)

test.lipidorder:
	@-cd Test_LipidOrder && ./RunTest.sh $(OPT)
	@-cd Test_LipidSCD && ./RunTest.sh $(OPT)

test.pairdist:
	@-cd Test_PairDist && ./RunTest.sh $(OPT)

test.gist:
	@-cd Test_GIST && ./RunTest.sh $(OPT)
	@-cd Test_GIST_Order && ./RunTest.sh $(OPT)

test.readdata:
	@-cd Test_ReadData && ./RunTest.sh $(OPT)

test.systemVF:
	@-cd Test_systemVF && ./RunTest.sh $(OPT)

test.atomiccorr:
	@-cd Test_AtomicCorr && ./RunTest.sh $(OPT)

test.autocorr:
	@-cd Test_AutoCorr && ./RunTest.sh $(OPT)

test.box:
	@-cd Test_Box && ./RunTest.sh $(OPT)

test.atomtypeselect:
	@-cd Test_AtomTypeSelect && ./RunTest.sh $(OPT)

test.crankstat:
	@-cd Test_Crank_Stat && ./RunTest.sh $(OPT)

test.createcrd:
	@-cd Test_CreateCrd && ./RunTest.sh $(OPT)

test.autoimage:
	@-cd Test_AutoImage && ./RunTest.sh $(OPT)

test.loadcrd:
	@-cd Test_LoadCrd && ./RunTest.sh $(OPT)

test.multivector:
	@-cd Test_MultiVector && ./RunTest.sh $(OPT)

test.lifetime:
	@-cd Test_Lifetime && ./RunTest.sh $(OPT)

test.fft:
	@-cd Test_FFT && ./RunTest.sh $(OPT)

test.timecorr:
	@-cd Test_Timecorr && ./RunTest.sh $(OPT)

test.crosscorr:
	@-cd Test_CrossCorr && ./RunTest.sh $(OPT)

test.symmrmsd:
	@-cd Test_SymmRmsd && ./RunTest.sh $(OPT)

test.multidihedral:
	@-cd Test_Multidihedral && ./RunTest.sh $(OPT)

test.mremd:
	@-cd Test_Ensemble_MREMD && ./RunTest.sh $(OPT)

test.makestructure:
	@-cd Test_MakeStructure && ./RunTest.sh $(OPT)

test.rotate:
	@-cd Test_Rotate && ./RunTest.sh $(OPT)

test.integrate:
	@-cd Test_Integrate && ./RunTest.sh $(OPT)

test.nativecontacts:
	@-cd Test_NativeContacts && ./RunTest.sh $(OPT)

test.checkchirality:
	@-cd Test_CheckChiral && ./RunTest.sh $(OPT)

test.energy:
	@-cd Test_Energy && ./RunTest.sh $(OPT)

test.lessplit:
	@-cd Test_LESsplit && ./RunTest.sh $(OPT)

test.permutedih:
	@-cd Test_PermuteDihedrals && ./RunTest.sh $(OPT)

test.ccp4:
	@-cd Test_CCP4 && ./RunTest.sh $(OPT)

test.dihcovar:
	@-cd Test_DihCovar && ./RunTest.sh $(OPT)

test.curvefit:
	@-cd Test_CurveFit && ./RunTest.sh $(OPT)
	@-cd Test_Multicurve && ./RunTest.sh $(OPT)

test.datasetavg:
	@-cd Test_DataSetAverage && ./RunTest.sh $(OPT)

test.vectormath:
	@-cd Test_VectorMath && ./RunTest.sh $(OPT)

test.esander:
	@-cd Test_Esander && ./RunTest.sh $(OPT)
	@-cd Test_Esander_Short && ./RunTest.sh $(OPT)

test.temperature:
	@-cd Test_Temperature && ./RunTest.sh $(OPT)

test.spline:
	@-cd Test_Spline && ./RunTest.sh $(OPT)

test.kde:
	@-cd Test_KDE && ./RunTest.sh $(OPT)

test.wavelet:
	@-cd Test_Wavelet && ./RunTest.sh $(OPT)

test.cif:
	@-cd Test_CIF && ./RunTest.sh $(OPT)

test.pdb:
	@-cd Test_PDB && ./RunTest.sh $(OPT)

test.align:
	@-cd Test_Align && ./RunTest.sh $(OPT)

test.velocityautocorr:
	@-cd Test_VelocityAutoCorr && ./RunTest.sh $(OPT)

test.velfrc:
	@-cd Test_VelFrc && ./RunTest.sh $(OPT)

test.rotatedihedral:
	@-cd Test_RotateDihedral && ./RunTest.sh $(OPT)

test.dataset:
	@-cd Test_DataSetCmd && ./RunTest.sh $(OPT)
	@-cd Test_Parallel_Dataset && ./RunTest.sh $(OPT)

test.remap:
	@-cd Test_Remap && ./RunTest.sh $(OPT)

test.comparetop:
	@-cd Test_CompareTop && ./RunTest.sh $(OPT)

test.filter:
	@-cd Test_Filter && ./RunTest.sh $(OPT)

test.combinecrd:
	@-cd Test_CombineCrd && ./RunTest.sh $(OPT)

test.topinfo:
	@-cd Test_TopInfo && ./RunTest.sh $(OPT)

test.gromacs:
	@-cd Test_GromacsTrr && ./RunTest.sh $(OPT)
	@-cd Test_GromacsXtc && ./RunTest.sh $(OPT)
	@-cd Test_GromacsTNG && ./RunTest.sh $(OPT)

test.ewald:
	@-cd Test_Ewald && ./RunTest.sh $(OPT)
	@-cd Test_LibPME && ./RunTest.sh $(OPT)
	@-cd Test_LJ_PME && ./RunTest.sh $(OPT)

test.ti:
	@-cd Test_TI && ./RunTest.sh $(OPT)

test.fiximagedbonds:
	@-cd Test_FixImagedBonds && ./RunTest.sh $(OPT)

test.distance:
	@-cd Test_Distance && ./RunTest.sh $(OPT)

test.volume:
	@-cd Test_Volume && ./RunTest.sh $(OPT)

test.setvelocity:
	@-cd Test_SetVelocity && ./RunTest.sh $(OPT)

test.remlog:
	@-cd Test_Remlog && ./RunTest.sh $(OPT)

test.areapermol:
	@-cd Test_AreaPerMol && ./RunTest.sh $(OPT)

test.distbasedmask:
	@-cd Test_DistBasedMask && ./RunTest.sh $(OPT)

test.calc:
	@-cd Test_Calc && ./RunTest.sh $(OPT)

test.calcstate:
	@-cd Test_CalcState && ./RunTest.sh $(OPT)

test.change:
	@-cd Test_Change && ./RunTest.sh $(OPT)

test.createreservoir:
	@-cd Test_CreateReservoir && ./RunTest.sh $(OPT)

test.fixatomorder:
	@-cd Test_FixAtomOrder && ./RunTest.sh $(OPT)

test.lowestcurve:
	@-cd Test_LowestCurve && ./RunTest.sh $(OPT)

test.minimage:
	@-cd Test_MinImage && ./RunTest.sh $(OPT)

test.regress:
	@-cd Test_Regress && ./RunTest.sh $(OPT)

test.rst:
	@-cd Test_Rst && ./RunTest.sh $(OPT)

test.pairwise:
	@-cd Test_Pairwise && ./RunTest.sh $(OPT)

test.parmbox:
	@-cd Test_ParmBox && ./RunTest.sh $(OPT)

test.solvent:
	@-cd Test_Solvent && ./RunTest.sh $(OPT)

test.translate:
	@-cd Test_Translate && ./RunTest.sh $(OPT)

test.replicatecell:
	@-cd Test_ReplicateCell && ./RunTest.sh $(OPT)

test.scaledihedralk:
	@-cd Test_ScaleDihedralK && ./RunTest.sh $(OPT)

test.phipsi:
	@-cd Test_PhiPsi && ./RunTest.sh $(OPT)

test.control:
	@-cd Test_Control && ./RunTest.sh $(OPT)

test.infraredspec:
	@-cd Test_InfraredSpec && ./RunTest.sh $(OPT)

test.cphstats:
	@-cd Test_Cphstats && ./RunTest.sh $(OPT)

test.cmdline:
	@-cd Test_CmdLine && ./RunTest.sh $(OPT)

test.xyzfmt:
	@-cd Test_XYZformat && ./RunTest.sh $(OPT)

test.hausdorff:
	@-cd Test_Hausdorff && ./RunTest.sh $(OPT)

test.catcrd:
	@-cd Test_CatCrd && ./RunTest.sh $(OPT)

test.xtalsymm:
	@-cd Test_XtalSymm && ./RunTest.sh $(OPT)

test.time:
	@-cd Test_Time && ./RunTest.sh $(OPT)

test.splitcoords:
	@-cd Test_SplitCoords && ./RunTest.sh $(OPT)

test.dihedralrms:
	@-cd Test_DihedralRMS && ./RunTest.sh $(OPT)

test.slope:
	@-cd Test_Slope && ./RunTest.sh $(OPT)

test.diagmatrix:
	@-cd Test_Diagmatrix && ./RunTest.sh $(OPT)

test.emin:
	@-cd Test_Emin && ./RunTest.sh $(OPT)

test.evalplateau:
	@-cd Test_EvalPlateau && ./RunTest.sh $(OPT)

test.graft:
	@-cd Test_Graft && ./RunTest.sh $(OPT)

test.flatten:
	@-cd Test_Flatten && ./RunTest.sh $(OPT)

test.random:
	@-cd Test_Random && ./RunTest.sh $(OPT)

test.multipucker:
	@-cd Test_MultiPucker && ./RunTest.sh $(OPT)

test.ljc:
	@-cd Test_LJC && ./RunTest.sh $(OPT)

test.hmassrepartition:
	@-cd Test_HmassRepartition && ./RunTest.sh $(OPT)

test.prepareforleap:
	@-cd Test_PrepareForLeap && ./RunTest.sh $(OPT)

test.keep:
	@-cd Test_Keep && ./RunTest.sh $(OPT)

test.hremd:
	@-cd Test_Ensemble_HREMD && ./RunTest.sh $(OPT)

test.ndiff:
	@-cd Test_NDiff && ./RunTest.sh $(OPT)

test.netcdfdata:
	@-cd Test_NetcdfData && ./RunTest.sh $(OPT)

test.compareenergy:
	@-cd Test_CompareEnergy && ./RunTest.sh $(OPT)

test.avgbox:
	@-cd Test_AvgBox && ./RunTest.sh $(OPT)

test.tordiff:
	@-cd Test_TorDiff && ./RunTest.sh $(OPT)

test.calcdiffusion:
	@-cd Test_CalcDiffusion && ./RunTest.sh $(OPT)

test.parsetiming:
	@-cd Test_ParseTiming && ./RunTest.sh $(OPT)

test.readprep:
	@-cd Test_ReadPrep && ./RunTest.sh $(OPT)

test.zmatrix:
	@-cd Test_Zmatrix && ./RunTest.sh $(OPT)

test.readoff:
	@-cd Test_ReadOFF && ./RunTest.sh $(OPT)

test.sequence:
	@-cd Test_Sequence && ./RunTest.sh $(OPT)

test.tica:
	@-cd Test_TICA && ./RunTest.sh $(OPT)

test.enedecomp:
	@-cd Test_EneDecomp && ./RunTest.sh $(OPT)

# Only GPU-specific tests should go here
GPUTESTS=test.closest \
         test.solvent \
         test.watershell \
         test.gist \
         test.radial

# Every test target should go here.
COMPLETETESTS=test.general \
              test.strip \
              test.brokentraj \
              test.trajinoffset \
              test.remdnetcdf \
              test.framerange \
              test.perres \
              test.dssp \
              test.center \
              test.hbond \
              test.ncrestart \
              test.image \
              test.surf \
              test.radgyr \
              test.mask \
              test.jcoupling \
              test.diffusion \
              test.outtraj \
              test.charmm \
              test.rotdif \
              test.netcdftraj \
              test.runavg \
              test.rmsd \
              test.2drms \
              test.mol2 \
              test.nastruct \
              test.average \
              test.pucker \
              test.atommap \
              test.datafile \
              test.radial \
              test.hist \
              test.closest \
              test.drmsd \
              test.cluster \
              test.ired \
              test.randomizeions \
              test.molsurf \
              test.corr \
              test.checkstructure \
              test.rmsavgcorr \
              test.contacts \
              test.principal \
              test.clusterdihedral \
              test.vector \
              test.matrix \
              test.atomicfluct \
              test.grid \
              test.unwrap \
              test.project \
              test.analyzemodes \
              test.watershell \
              test.analrunavg \
              test.volmap \
              test.spam \
              test.lie \
              test.density \
              test.lipidorder \
              test.pairdist \
              test.readdata \
              test.gist \
              test.systemVF \
              test.atomiccorr \
              test.autocorr \
              test.box \
              test.atomtypeselect \
              test.avgref \
              test.crankstat \
              test.createcrd \
              test.autoimage \
              test.loadcrd \
              test.multivector \
              test.lifetime \
              test.fft \
              test.timecorr \
              test.crosscorr \
              test.symmrmsd \
              test.multidihedral \
              test.mremd \
              test.hremd \
              test.makestructure \
              test.rotate \
              test.integrate \
              test.nativecontacts \
              test.checkchirality \
              test.energy \
              test.lessplit \
              test.permutedih \
              test.ccp4 \
              test.dihcovar \
              test.curvefit \
              test.datasetavg \
              test.vectormath \
              test.esander \
              test.temperature \
              test.spline \
              test.kde \
              test.wavelet \
              test.cif \
              test.pdb \
              test.align \
              test.velocityautocorr \
              test.velfrc \
              test.rotatedihedral \
              test.dataset \
              test.remap \
              test.comparetop \
              test.filter \
              test.combinecrd \
              test.topinfo \
              test.gromacs \
              test.ewald \
              test.ti \
              test.fiximagedbonds \
              test.distance \
              test.volume \
              test.setvelocity \
              test.remlog \
              test.areapermol \
              test.distbasedmask \
              test.calc \
              test.calcstate \
              test.change \
              test.createreservoir \
              test.fixatomorder \
              test.lowestcurve \
              test.minimage \
              test.regress \
              test.rst \
              test.pairwise \
              test.parmbox \
              test.solvent \
              test.translate \
              test.replicatecell \
              test.scaledihedralk \
              test.phipsi \
              test.control \
              test.infraredspec \
              test.cphstats \
              test.cmdline \
              test.xyzfmt \
              test.hausdorff \
              test.catcrd \
              test.xtalsymm \
              test.time \
              test.splitcoords \
              test.dihedralrms \
              test.slope \
              test.diagmatrix \
              test.emin \
              test.evalplateau \
              test.graft \
              test.flatten \
              test.random \
              test.multipucker \
              test.ljc \
              test.hmassrepartition \
              test.prepareforleap \
              test.keep \
              test.ndiff \
              test.netcdfdata \
              test.compareenergy \
              test.avgbox \
              test.tordiff \
              test.calcdiffusion \
              test.parsetiming \
              test.readprep \
              test.zmatrix \
              test.readoff \
              test.sequence \
              test.tica \
              test.enedecomp

test.all:
	$(MAKE) test.complete summary

test.mpi:
	$(MAKE) test.complete summary OPT="$(OPT) mpi"

test.openmp:
	$(MAKE) test.complete summary OPT="$(OPT) openmp"

test.cuda:
	$(MAKE) test.gpu summary OPT="$(OPT) cuda"

test.hip:
	$(MAKE) test.gpu summary OPT="$(OPT) hip"

test.libcpptraj:
	@cd Test_Libcpptraj && ./RunTest.sh $(OPT)

test.complete: CpptrajTest.sh MasterTest.sh
	@./CpptrajTest.sh --target test.cpptraj $(OPT)

test.gpu: CpptrajTest.sh MasterTest.sh
	@./CpptrajTest.sh --target test.cpptraj.gpu $(OPT)

test.cpptraj: $(COMPLETETESTS)

test.cpptraj.gpu: $(GPUTESTS)

test.showerrors:
	$(MAKE) test.complete summary OPT="$(OPT) showerrors long"

summary: CpptrajTest.sh MasterTest.sh
	@./CpptrajTest.sh summary $(OPT)

test.vg:
	$(MAKE) test.complete summary OPT="$(OPT) vg"

clean:
	@echo "Cleaning all test directories:"
	$(MAKE) test.complete test.libcpptraj OPT=clean

cleanprofile:
	@echo "Removing GNU profile files:"
	(find . -name 'libstdcxx-profile.*' | \
        while read dif ;\
        do \
            rm -f $$dif ;\
        done ;\
        )

