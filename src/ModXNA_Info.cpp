#include "ModXNA_Info.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
ModXNA_Info::ModXNA_Info() {}

/** Given a ModXNA metadata string, return an array of just the
  * atom names.
  */
ModXNA_Info::Sarray ModXNA_Info::atomNameList(std::string const& mheader, int& errval)
{
  Sarray ListOut;
  if (mheader.empty()) return ListOut;
  ArgList anames( mheader, "@," );
  for (int iarg = 0; iarg < anames.Nargs(); iarg++)
    ListOut.push_back( anames[iarg] );
  if (ListOut.empty()) {
    mprinterr("Error: ModXNA header string '%s' has no valid atom names.\n", mheader.c_str());
    errval++;
  }
  return ListOut;
}

/** Parse ModXNA info from given string. */
ModXNA_Info::StatType ModXNA_Info::ParseModxnaStr(std::string const& line) {
  if (line.empty()) return NOT_MODXNA;
  ArgList args(line, ":");
  //args.PrintDebug(); // DEBUG

  if (args.GetStringNext() != "modXNA-fragment") return NOT_MODXNA;

  if (args.Nargs() < 2) {
    mprinterr("Error: modXNA string is too short.\n"
              "Error: '%s'\n", line.c_str());
    return ERR;
  }
  fragmentName_ = args[1];

  int errval = 0;
  head_ = atomNameList(args.GetStringKey("HEAD01"), errval);
  headStrip_ = atomNameList(args.GetStringKey("HEAD01STRIP"), errval);
  tail_ = atomNameList(args.GetStringKey("TAIL01"), errval);
  tailStrip_ = atomNameList(args.GetStringKey("TAIL01STRIP"), errval);
  anchor_ = atomNameList(args.GetStringKey("HEAD03"), errval);
  anchorStrip_ = atomNameList(args.GetStringKey("HEAD03STRIP"), errval);

  if (errval > 0) return ERR;
  return OK;
}

/** Generate mask string from atom names. */
std::string ModXNA_Info::atomMaskStr(Sarray const& anames)
{
  if (anames.empty()) return std::string("");
  std::string maskstr("@" + anames[0]);
  for (unsigned int idx = 1; idx < anames.size(); idx++)
    maskstr.append("," + anames[idx]);
  return maskstr;
}

/** Print name array to stdout. */
void ModXNA_Info::printmxna(const char* desc, Sarray const& anames) {
  if (!anames.empty()) {
    std::string str = atomMaskStr(anames);
    if (!str.empty()) mprintf("%s: %s\n", desc, str.c_str());
  }
}

/** Print info to stdout */
void ModXNA_Info::PrintModxna() const {
  if (!HasModxna()) return;
  mprintf("\tFragment: %s\n", fragmentName_.c_str());
  printmxna("\tHead        ", head_);
  printmxna("\tHeadStrip   ", headStrip_);
  printmxna("\tTail        ", tail_);
  printmxna("\tTailStrip   ", tailStrip_);
  printmxna("\tAnchor      ", anchor_);
  printmxna("\tAnchorStrip ", anchorStrip_);
}

/** \return specified array of names */
ModXNA_Info::Sarray const& ModXNA_Info::AtomNames(MaskType typeIn) const {
  switch (typeIn) {
    case HEAD        : return head_;
    case HEADSTRIP   : return headStrip_;
    case TAIL        : return tail_;
    case TAILSTRIP   : return tailStrip_;
    case ANCHOR      : return anchor_;
    case ANCHORSTRIP : return anchorStrip_;
  }
  return head_;
}

/** \return Atom mask string for the specified mask type */
std::string ModXNA_Info::AtomMaskString(MaskType typeIn) const {
  switch (typeIn) {
    case HEAD        : return atomMaskStr(head_);
    case HEADSTRIP   : return atomMaskStr(headStrip_);
    case TAIL        : return atomMaskStr(tail_);
    case TAILSTRIP   : return atomMaskStr(tailStrip_);
    case ANCHOR      : return atomMaskStr(anchor_);
    case ANCHORSTRIP : return atomMaskStr(anchorStrip_);
  }
  return std::string("");
}

/** Generate ModXNA metadata string */
std::string ModXNA_Info::GenMetadataString() const {
  if (fragmentName_.empty()) return std::string("");
  std::string out("modXNA-fragment:" + fragmentName_);
  // Base
  out.append(":HEAD01:" + atomMaskStr(head_));
  out.append(":HEAD01STRIP:" + atomMaskStr(headStrip_));
  if (!tail_.empty()) {
    // BB
    out.append(":TAIL01:" + atomMaskStr(tail_));
    out.append(":TAIL01STRIP:" + atomMaskStr(tailStrip_));
    if (!anchor_.empty()) {
      // Sugar
      out.append(":HEAD03:" + atomMaskStr(anchor_));
      out.append(":HEAD03STRIP:" + atomMaskStr(anchorStrip_));
    }
  }
  return out;
}


