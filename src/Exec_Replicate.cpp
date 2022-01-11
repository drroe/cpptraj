#include "Exec_Replicate.h"
#include "CpptrajStdio.h"
#include "TopInfo.h" // DEBUG

// Exec_Replicate::Help()
void Exec_Replicate::Help() const
{
  //mprintf("\t<mask> [nrep <#>] [%s]\n", DataSetList::TopIdxArgs);
  mprintf("\tcrdset <coords set> [name <out coords>] [<mask>] [nrep <#>]\n");
}

// Exec_Replicate::Execute()
/** Replicate all or part of a system in a COORDS set. If there is more than
  * 1 frame in the incoming COORDS set, assume we want that many replicates.
  */
Exec::RetType Exec_Replicate::Execute(CpptrajState& State, ArgList& argIn)
{
  int nrep = argIn.getKeyInt("nrep", 1);
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

  // Get the first frame
  Frame frameIn = CRD->AllocateFrame();
  CRD->GetFrame(0, frameIn);

/*  int tgtframe = argIn.getKeyInt("frame", 0);
  if (tgtframe < 0 || tgtframe >= (int)CRD->Size()) {
    mprinterr("Error: Specified frame %i is out of range.\n", tgtframe+1);
    return CpptrajState::ERR;
  }
  Frame frameIn = CRD->AllocateFrame();
  CRD->GetFrame(tgtframe, frameIn);
  mprintf("    REPLICATE: Using COORDS '%s', frame %i\n", CRD->legend(), tgtframe+1);*/
  // ----------
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
  // ----------
  
  // Topology from input coords will be modified.
  Topology topIn = CRD->Top();

  /*Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  // Check if this topology has already been used to set up an input
  // trajectory, as this will break the traj read.
  std::string fname = State.TopUsedInInputTraj(parm);
  if (!fname.empty()) {
    mprinterr("Error: Topology '%s' has already been used to set up trajectory '%s'.\n"
              "Error:   To strip this topology use the 'strip' action.\n",
              parm->c_str(), fname.c_str());
    return CpptrajState::ERR;
  }*/
  AtomMask tempMask( argIn.GetMaskNext() );

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
