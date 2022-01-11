#include "Exec_Replicate.h"
#include "CpptrajStdio.h"
#include "TopInfo.h" // DEBUG

// Exec_Replicate::Help()
void Exec_Replicate::Help() const
{
  //mprintf("\t<mask> [nrep <#>] [%s]\n", DataSetList::TopIdxArgs);
  mprintf("\tcrdset <coords set> [useframes <range>|frame <#>] [name <out coords>]\n"
          "\t[<mask>] [nrep <#>]\n");
}

// Exec_Replicate::Execute()
/** Replicate all or part of a system in a COORDS set.
  */
Exec::RetType Exec_Replicate::Execute(CpptrajState& State, ArgList& argIn)
{
  int nrep = argIn.getKeyInt("nrep", -1);

  // Get input COORDS set
  std::string setname = argIn.GetStringKey("crdset");
  if (setname.empty()) {
    mprinterr("Error: Specify COORDS dataset name with 'crdset'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: Could not find COORDS set '%s'\n", setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("    REPLICATE: Using COORDS '%s'\n", CRD->legend());
  if (CRD->Size() < 1) {
    mprinterr("Error: COORDS set is empty.\n");
    return CpptrajState::ERR;
  }

  Range useFrames;
  std::string useFramesArg = argIn.GetStringKey("useframes");
  if (!useFramesArg.empty()) {
    if (useFrames.SetRange( useFramesArg, Range::UNSORTED )) {
      mprinterr("Error: Could not set range '%s'\n", useFramesArg.c_str());
      return CpptrajState::ERR;
    }
    // User frame args start from 1
    useFrames.ShiftBy(-1);
    if (useFrames.Size() < 2) {
      mprinterr("Error: Need to specify at least 2 frames for 'useframes'\n");
      return CpptrajState::ERR;
    }
    if (nrep == -1)
      nrep = useFrames.Size() - 1;
    else if (nrep != useFrames.Size()) {
      mprinterr("Error: # reps (%i) not equal to # frames to use (%i).\n", nrep, useFrames.Size());
      return CpptrajState::ERR;
    }
  } else {
    useFrames = Range( argIn.getKeyInt("frame", 1) - 1 );
    if (nrep == -1)
      nrep = 1;
  }
  mprintf("\t%i replicates\n", nrep);
  // Get the first frame
  Range::const_iterator frmnum = useFrames.begin();
  Frame frameIn = CRD->AllocateFrame();
  CRD->GetFrame(*(frmnum++), frameIn);

  
  // Create output COORDS set if necessary
  DataSet_Coords* OUT = 0;
//  int outframe = 0;
  std::string outname = argIn.GetStringKey("name");
  if (outname.empty()) {
    // This will not work for TRAJ data sets
    if (CRD->Type() == DataSet::TRAJ) {
      mprinterr("Error: Using TRAJ as input set requires use of 'name' keyword for output.\n");
      return CpptrajState::ERR;
    }
    //OUT = CRD;
    //outframe = tgtframe;
  } else {
    // Create new output set with 1 empty frame.
    OUT = (DataSet_Coords*)State.DSL().AddSet( DataSet::COORDS, outname );
    if (OUT == 0) return CpptrajState::ERR;
    OUT->Allocate( DataSet::SizeArray(1, 1) );
    //OUT->CoordsSetup( CRD->Top(), CRD->CoordsInfo() );
    //OUT->AddFrame( CRD->AllocateFrame() );
    mprintf("\tOutput to set '%s'\n", OUT->legend());
  }
  
  AtomMask tempMask( argIn.GetMaskNext() );
  // -------------------------

  // Topology from input coords will be modified.
  Topology topIn = CRD->Top();

  if (topIn.SetupIntegerMask( tempMask )) return CpptrajState::ERR;
  mprintf("\tReplicating atoms in topology '%s' %i times.\n", topIn.c_str(), nrep);
  tempMask.MaskInfo();

  if (topIn.ReplicateAtoms( tempMask, nrep )) {
    mprinterr("Error: Replication failed for topology.\n");
    return CpptrajState::ERR;
  }

  TopInfo tinfo( &topIn );
  tinfo.PrintReplicateInfo("*");

  if (frameIn.ReplicateFrameAtoms( tempMask, nrep )) {
    mprinterr("Error: Replication failed for frame.\n");
    return CpptrajState::ERR;
  }

  // Set

  // Assign to output coords
  CoordinateInfo cinfo = CRD->CoordsInfo();
  if (OUT == 0) {
    // Want to replace set CRD with new top/coords.
    MetaData md = CRD->Meta();
    State.DSL().RemoveSet( CRD );
    CRD = 0;
    OUT = (DataSet_Coords*)State.DSL().AddSet( DataSet::COORDS, md );
  }
  OUT->CoordsSetup( topIn, cinfo );
  OUT->AddFrame( frameIn );


  return CpptrajState::OK;
}
