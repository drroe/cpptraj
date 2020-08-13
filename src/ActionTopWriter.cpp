#include "ActionTopWriter.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "Topology.h"
#include "ParmFile.h"
#include "FileName.h"
#include "DataSetList.h"
#include "DataSet_Topology.h"

static const char* keywords_ = "\t[outprefix <prefix>] [parmout <filename>]\n"
                               "\t[parmopts <comma-separated-list>]\n";
static const char* options_ = 
          "    outprefix <prefix> : Write modified topology to <prefix>.<originalname>\n"
          "    parmout <filename> : Write modified topology to <filename>\n"
          "    parmopts <list>    : Options for writing topology file\n";

const char* ActionTopWriter::Keywords() { return keywords_; }
const char* ActionTopWriter::Options() { return options_; }

/** CONSTRUCTOR */
ActionTopWriter::ActionTopWriter() :
  debug_(0)
{}

/** Parse arguments. */
int ActionTopWriter::InitTopWriter(ArgList& argIn, const char* typeStrIn, int debugIn,
                                   DataSetList* DSLin)
{
  if (DSLin == 0) {
    mprinterr("Internal Error: InitTopWriter(): DSLin is null.\n");
    return 1;
  }
  debug_ = debugIn;
  prefix_ = argIn.GetStringKey("outprefix");
  parmoutName_ = argIn.GetStringKey("parmout");
  parmOpts_ = argIn.GetStringKey("parmopts");
  if (typeStrIn==0) {
    mprinterr("Internal Error: No type string given for ActionTopWriter::InitTopWriter\n");
    return 1;
  }
  typeStr_.assign(typeStrIn);
  masterDSL_ = DSLin;
  return 0;
}

/** Print arguments to stdout. */
void ActionTopWriter::PrintOptions() const {
  if (!prefix_.empty())
    mprintf("\tWriting '%s' topology file with prefix '%s'\n", typeStr_.c_str(), prefix_.c_str());
  if (!parmoutName_.empty())
    mprintf("\tWriting '%s' topology file with name '%s'\n", typeStr_.c_str(), parmoutName_.c_str());
  if (!parmOpts_.empty())
    mprintf("\tUsing the following write options: %s\n", parmOpts_.c_str());
}

/** Create a DataSet_Topology based on the given names. */
DataSet_Topology* ActionTopWriter::CreateTopSet(Topology const& oldTop)
{
  DataSet_Topology* topSet = 0;
  // NOTE: We identify "old" (i.e. non-modified) Topology by memory address,
  //       not Pindex(), since modified Topologies are given the same
  //       Pindex as the Topology they are based on.
  // Has this topology already been modified?
  Imap::const_iterator it = modifiedTopMap_.find( &oldTop );
  if (it != modifiedTopMap_.end()) {
    mprintf("\tTopology %s has been previously modified.\n", oldTop.c_str());
    // This topology has been previously modified here. Check and return.
    topSet = it->second;
    if (topSet->TopPtr() == 0) {
      mprinterr("Internal Error: CreateTopSet(): Previously modified Topology set is null.\n");
      return 0;
    }
  } else {
    // Create a new set.
    MetaData md;
    // Determine the filename for the DataSet TODO tag
    FileName fname;
    if (!parmoutName_.empty())
      fname = parmoutName_;
    else if (!prefix_.empty())
      fname = oldTop.OriginalFilename().PrependFileName(prefix_);
    if (fname.empty()) {
      // No file name. Use generic name
      md = MetaData( masterDSL_->GenerateDefaultName("TOP") );
      mprintf("\tModified topology data set name: %s\n", md.Name().c_str());
    } else {
      mprintf("\tModified topology data set name: %s\n", fname.full());
      md = MetaData(fname, "", -1);
    }
    topSet = (DataSet_Topology*)masterDSL_->AddSet(DataSet::TOPOLOGY, md);
    if (topSet == 0) {
      mprinterr("Internal Error: CreateTopSet(): Could not allocate Topology DataSet.\n");
      return 0;
    }
    modifiedTopMap_.insert( Ipair(&oldTop, topSet) );
  }
  return topSet;
}
  
/** Write Topology to topology file(s).
  * \param topIn Topology to write.
  * \return 0 if OK, 1 if prefix write failed, 2 if parmout write failed, 3 if both failed.
  */
int ActionTopWriter::WriteTops(Topology const& topIn) const {
  ArgList writeOpts;
  int err = 0;
  if (!parmOpts_.empty())
    writeOpts.SetList(parmOpts_, ",");
  if (!prefix_.empty()) {
    ParmFile pfile;
    if ( pfile.WritePrefixTopology(topIn, prefix_, writeOpts, ParmFile::UNKNOWN_PARM, debug_) )
    {
      mprinterr("Error: Could not write out '%s' topology file with prefix '%s'.\n",
                typeStr_.c_str(), prefix_.c_str());
      err = 1;
    }
  }
  if (!parmoutName_.empty()) {
    ParmFile pfile;
    if (pfile.WriteTopology(topIn, parmoutName_, writeOpts, ParmFile::UNKNOWN_PARM, debug_))
    {
      mprinterr("Error: Could not write out '%s' topology file '%s'\n",
                typeStr_.c_str(), parmoutName_.c_str());
      err += 2;
    }
  }
  return err;
}
