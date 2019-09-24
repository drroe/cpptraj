#include "Residue.h"
#include <cctype> // tolower

const char Residue::BLANK_CHAINID_ = '\0';

const char Residue::DEFAULT_CHAINID_ = 'Z';

char Residue::ConvertResName(std::string const& r) {
  if (r.compare(0,3,"ALA")==0) return 'A';
  if (r.compare(0,3,"ARG")==0) return 'R';
  if (r.compare(0,3,"ASN")==0) return 'N';
  if (r.compare(0,3,"ASP")==0) return 'D';
  if (r.compare(0,3,"ASH")==0) return 'D'; // Protonated ASP
  if (r.compare(0,3,"AS4")==0) return 'D'; // Constant pH ASP
  if (r.compare(0,3,"CYS")==0) return 'C';
  if (r.compare(0,3,"CYM")==0) return 'C'; // Deprotonated CYS
  if (r.compare(0,3,"CYX")==0) return 'C';
  if (r.compare(0,3,"GLN")==0) return 'Q';
  if (r.compare(0,3,"GLU")==0) return 'E';
  if (r.compare(0,3,"GLH")==0) return 'E'; // Protonated GLU
  if (r.compare(0,3,"GL4")==0) return 'E'; // Constant pH GLU 
  if (r.compare(0,3,"GLY")==0) return 'G';
  if (r.compare(0,3,"HIS")==0) return 'H';
  if (r.compare(0,3,"HIE")==0) return 'H'; // NE-protonated (HIS)
  if (r.compare(0,3,"HID")==0) return 'H'; // ND-protonated
  if (r.compare(0,3,"HIP")==0) return 'H'; // NE/ND protonated
  if (r.compare(0,3,"ILE")==0) return 'I';
  if (r.compare(0,3,"LEU")==0) return 'L';
  if (r.compare(0,3,"LYS")==0) return 'K';
  if (r.compare(0,3,"LYN")==0) return 'K'; // Deprotonated (neutral) LYS 
  if (r.compare(0,3,"MET")==0) return 'M';
  if (r.compare(0,3,"PHE")==0) return 'F';
  if (r.compare(0,3,"PRO")==0) return 'P';
  if (r.compare(0,3,"SER")==0) return 'S';
  if (r.compare(0,3,"THR")==0) return 'T';
  if (r.compare(0,3,"TRP")==0) return 'W';
  if (r.compare(0,3,"TYR")==0) return 'Y';
  if (r.compare(0,3,"VAL")==0) return 'V';
  // Nucleic acids
  if (r.compare(0,2,"DA")==0 || r.compare(0,1,"A")==0) return 'A';
  if (r.compare(0,2,"DG")==0 || r.compare(0,1,"G")==0) return 'G';
  if (r.compare(0,2,"DC")==0 || r.compare(0,1,"C")==0) return 'C';
  if (r.compare(0,2,"DT")==0 || r.compare(0,1,"T")==0) return 'T';
  if (r.compare(0,1,"U")==0) return 'U';
  // Make lower case letter when unrecognized.
  if (!r.empty()) return tolower(r[0]);
  return ' ';
}

const char* Residue::ConvertResName(char letter) {
  switch (letter) {
      case 'A': return "ALA";
      case 'R': return "ARG";
      case 'N': return "ASN";
      case 'D': return "ASP";
  //if (r.compare(0,3,"ASH")==0) return 'D'; // Protonated ASP
      case 'C': return "CYS";
  //if (r.compare(0,3,"CYM")==0) return 'C'; // Deprotonated CYS
  //if (r.compare(0,3,"CYX")==0) return 'C';
      case 'Q': return "GLN";
      case 'E': return "GLU";
  //if (r.compare(0,3,"GLH")==0) return 'E'; // Protonated GLU
      case 'G': return "GLY";
      case 'H': return "HIS";
  //if (r.compare(0,3,"HIE")==0) return 'H'; // NE-protonated (HIS)
  //if (r.compare(0,3,"HID")==0) return 'H'; // ND-protonated
  //if (r.compare(0,3,"HIP")==0) return 'H'; // NE/ND protonated
      case 'I': return "ILE";
      case 'L': return "LEU";
      case 'K': return "LYS";
  //if (r.compare(0,3,"LYN")==0) return 'K'; // Deprotonated (neutral) LYS 
      case 'M': return "MET";
      case 'F': return "PHE";
      case 'P': return "PRO";
      case 'S': return "SER";
      case 'T': return "THR";
      case 'W': return "TRP";
      case 'Y': return "TYR";
      case 'V': return "VAL";
  }
  return 0;
}

const char* Residue::resTypeStr_[] = {
  "protein", "nucleic", "lipid", "solvent", "other"
};

const char* Residue::ResTypeStr(ResidueType rt) {
  return resTypeStr_[rt];
}

Residue::ResNameMapType Residue::resNameMap_ = ResNameMapType();

void Residue::InitResNameMap() {
  // PROTEIN RESIDUES
  resNameMap_.insert( ResNamePairType("ACE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("ALA", PROTEIN) );
  resNameMap_.insert( ResNamePairType("ARG", PROTEIN) );
  resNameMap_.insert( ResNamePairType("ASH", PROTEIN) );
  resNameMap_.insert( ResNamePairType("ASN", PROTEIN) );
  resNameMap_.insert( ResNamePairType("ASP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CALA", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CARG", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CASN", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CASP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CCYS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CCYX", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CGLN", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CGLU", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CGLY", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CHID", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CHIE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CHIP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CHIS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CHYP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CILE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CLEU", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CLYS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CMET", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CPHE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CPRO", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CSER", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CTHR", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CTRP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CTYR", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CVAL", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CYM", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CYS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("CYX", PROTEIN) );
  resNameMap_.insert( ResNamePairType("GLH", PROTEIN) );
  resNameMap_.insert( ResNamePairType("GLN", PROTEIN) );
  resNameMap_.insert( ResNamePairType("GLU", PROTEIN) );
  resNameMap_.insert( ResNamePairType("GLY", PROTEIN) );
  resNameMap_.insert( ResNamePairType("HID", PROTEIN) );
  resNameMap_.insert( ResNamePairType("HIE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("HIP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("HIS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("HYP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("ILE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("LEU", PROTEIN) );
  resNameMap_.insert( ResNamePairType("LYN", PROTEIN) );
  resNameMap_.insert( ResNamePairType("LYS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("MET", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NALA", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NARG", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NASN", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NASP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NCYS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NCYX", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NGLN", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NGLU", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NGLY", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NHE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NHID", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NHIE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NHIP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NHIS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NILE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NLEU", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NLYS", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NME", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NMET", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NPHE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NPRO", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NSER", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NTHR", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NTRP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NTYR", PROTEIN) );
  resNameMap_.insert( ResNamePairType("NVAL", PROTEIN) );
  resNameMap_.insert( ResNamePairType("PHE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("PRO", PROTEIN) );
  resNameMap_.insert( ResNamePairType("SER", PROTEIN) );
  resNameMap_.insert( ResNamePairType("THR", PROTEIN) );
  resNameMap_.insert( ResNamePairType("TRP", PROTEIN) );
  resNameMap_.insert( ResNamePairType("TYR", PROTEIN) );
  resNameMap_.insert( ResNamePairType("VAL", PROTEIN) );
  // CHARMM PROTEIN RESIDUES
  resNameMap_.insert( ResNamePairType("HSD", PROTEIN) );
  resNameMap_.insert( ResNamePairType("HSE", PROTEIN) );
  resNameMap_.insert( ResNamePairType("HSP", PROTEIN) );
  // DNA RESIDUES
  resNameMap_.insert( ResNamePairType("DA", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DA3", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DA5", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DAN", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DC", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DC3", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DC5", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DCN", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DG", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DG3", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DG5", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DGN", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DT", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DT3", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DT5", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("DTN", NUCLEIC) );
  // RNA RESIDUES
  resNameMap_.insert( ResNamePairType("A", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("A3", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("A5", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("AMP", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("AN", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("C", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("C3", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("C5", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("CMP", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("CN", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("G", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("G3", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("G5", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("GMP", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("GN", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("OHE", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("U", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("U3", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("U5", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("UMP", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("UN", NUCLEIC) );
  // CHARMM NA RESIDUES
  resNameMap_.insert( ResNamePairType("GUA", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("ADE", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("CYT", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("THY", NUCLEIC) );
  resNameMap_.insert( ResNamePairType("URA", NUCLEIC) );
  // LIPID RESIDUES
  resNameMap_.insert( ResNamePairType("AR", LIPID) );
  resNameMap_.insert( ResNamePairType("CHL", LIPID) );
  resNameMap_.insert( ResNamePairType("DHA", LIPID) );
  resNameMap_.insert( ResNamePairType("LAL", LIPID) );
  resNameMap_.insert( ResNamePairType("MY", LIPID) );
  resNameMap_.insert( ResNamePairType("OL", LIPID) );
  resNameMap_.insert( ResNamePairType("OL2", LIPID) );
  resNameMap_.insert( ResNamePairType("PA", LIPID) );
  resNameMap_.insert( ResNamePairType("PC", LIPID) );
  resNameMap_.insert( ResNamePairType("PE", LIPID) );
  resNameMap_.insert( ResNamePairType("PGR", LIPID) );
  resNameMap_.insert( ResNamePairType("PH-", LIPID) );
  resNameMap_.insert( ResNamePairType("PS", LIPID) );
  resNameMap_.insert( ResNamePairType("ST", LIPID) );
  // CHARMM LIPID RESIDUES
  resNameMap_.insert( ResNamePairType("LPPC", LIPID) );
  resNameMap_.insert( ResNamePairType("DLPC", LIPID) );
  resNameMap_.insert( ResNamePairType("DLPE", LIPID) );
  resNameMap_.insert( ResNamePairType("DLPS", LIPID) );
  resNameMap_.insert( ResNamePairType("DLPA", LIPID) );
  resNameMap_.insert( ResNamePairType("DLPG", LIPID) );
  resNameMap_.insert( ResNamePairType("DMPC", LIPID) );
  resNameMap_.insert( ResNamePairType("DMPE", LIPID) );
  resNameMap_.insert( ResNamePairType("DMPS", LIPID) );
  resNameMap_.insert( ResNamePairType("DMPA", LIPID) );
  resNameMap_.insert( ResNamePairType("DMPG", LIPID) );
  resNameMap_.insert( ResNamePairType("DPPC", LIPID) );
  resNameMap_.insert( ResNamePairType("DPPE", LIPID) );
  resNameMap_.insert( ResNamePairType("DPPS", LIPID) );
  resNameMap_.insert( ResNamePairType("DPPA", LIPID) );
  resNameMap_.insert( ResNamePairType("DPPG", LIPID) );
  resNameMap_.insert( ResNamePairType("DSPC", LIPID) );
  resNameMap_.insert( ResNamePairType("DSPE", LIPID) );
  resNameMap_.insert( ResNamePairType("DSPS", LIPID) );
  resNameMap_.insert( ResNamePairType("DSPA", LIPID) );
  resNameMap_.insert( ResNamePairType("DSPG", LIPID) );
  resNameMap_.insert( ResNamePairType("DOPC", LIPID) );
  resNameMap_.insert( ResNamePairType("DOPE", LIPID) );
  resNameMap_.insert( ResNamePairType("DOPS", LIPID) );
  resNameMap_.insert( ResNamePairType("DOPA", LIPID) );
  resNameMap_.insert( ResNamePairType("DOPG", LIPID) );
  resNameMap_.insert( ResNamePairType("POPC", LIPID) );
  resNameMap_.insert( ResNamePairType("POPE", LIPID) );
  resNameMap_.insert( ResNamePairType("POPS", LIPID) );
  resNameMap_.insert( ResNamePairType("POPA", LIPID) );
  resNameMap_.insert( ResNamePairType("POPG", LIPID) );
  resNameMap_.insert( ResNamePairType("SAPC", LIPID) );
  resNameMap_.insert( ResNamePairType("SDPC", LIPID) );
  resNameMap_.insert( ResNamePairType("SOPC", LIPID) );
  resNameMap_.insert( ResNamePairType("DAPC", LIPID) );
  // SOLVENT RESIDUES
  resNameMap_.insert( ResNamePairType("WAT", SOLVENT) );
  resNameMap_.insert( ResNamePairType("HOH", SOLVENT) );
  resNameMap_.insert( ResNamePairType("H2O", SOLVENT) );
  resNameMap_.insert( ResNamePairType("TIP3", SOLVENT) );
  resNameMap_.insert( ResNamePairType("TIP4", SOLVENT) );
  resNameMap_.insert( ResNamePairType("TIP5", SOLVENT) );
  resNameMap_.insert( ResNamePairType("SOL", SOLVENT) );
}

Residue::ResidueType Residue::GetTypeFromName(NameType const& nameIn) {
  ResNameMapType::iterator it = resNameMap_.find( nameIn );
  if (it == resNameMap_.end())
    return UNKNOWN;
  else
    return it->second;
}
