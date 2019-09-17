#include "Residue.h"
#include <cctype> // tolower

const char Residue::BLANK_CHAINID_ = '\0';

const char Residue::DEFAULT_CHAINID_ = 'Z';

Residue::ResidueType Residue::TypeFromName(NameType const& nameIn) {
  switch (nameIn.len()) {
    case 4 :
      // 4 letters. Potentially Amber terminal protein or Charmm lipid
      switch (nameIn[0]) {
        case 'C' :
          switch (nameIn[1]) {
            case 'A' :
              if (nameIn[2] == 'L' && nameIn[3] == 'A') return PROTEIN; // C-Alanine
              if (nameIn[2] == 'R' && nameIn[3] == 'G') return PROTEIN; // C-Arginine
              if (nameIn[2] == 'S' && nameIn[3] == 'N') return PROTEIN; // C-Asparagine
              if (nameIn[2] == 'S' && nameIn[3] == 'P') return PROTEIN; // C-Aspartic acid
            break; // END case nameIn[1]==A
            case 'C' :
              if (nameIn[2] == 'Y' && nameIn[3] == 'S') return PROTEIN; // Protonated (normal) C-CYS
              if (nameIn[2] == 'Y' && nameIn[3] == 'X') return PROTEIN; // Disulphide C-CYS
            break; // END case nameIn[1]==C
            case 'G' :
              if (nameIn[2] == 'L' && nameIn[3] == 'N') return PROTEIN; // C-Glutamine
              if (nameIn[2] == 'L' && nameIn[3] == 'U') return PROTEIN; // C-Glutamic acid
              if (nameIn[2] == 'L' && nameIn[3] == 'Y') return PROTEIN; // C-Glycine
            break; // END case nameIn[1]==G
            case 'H' :
              if (nameIn[2] == 'I' && nameIn[3] == 'D') return PROTEIN;
              if (nameIn[2] == 'I' && nameIn[3] == 'E') return PROTEIN;
              if (nameIn[2] == 'I' && nameIn[3] == 'P') return PROTEIN;
              if (nameIn[2] == 'I' && nameIn[3] == 'S') return PROTEIN;
              if (nameIn[2] == 'Y' && nameIn[3] == 'P') return PROTEIN;
            break; // END case nameIn[1]==H
            case 'I' :
              if (nameIn[2] == 'L' && nameIn[3] == 'E') return PROTEIN;
            break; // END case nameIn[1]==I
            case 'L' :
              if (nameIn[2] == 'E' && nameIn[3] == 'U') return PROTEIN;
              if (nameIn[2] == 'Y' && nameIn[3] == 'S') return PROTEIN;
            break; // END case nameIn[1]==L
            case 'M' :
              if (nameIn[2] == 'E' && nameIn[3] == 'T') return PROTEIN;
            break; // END case nameIn[1]==M
            case 'P' :
              if (nameIn[2] == 'H' && nameIn[3] == 'E') return PROTEIN;
              if (nameIn[2] == 'R' && nameIn[3] == 'O') return PROTEIN;
            break; // END case nameIn[1]==P
            case 'S' :
              if (nameIn[2] == 'E' && nameIn[3] == 'R') return PROTEIN;
            break; // END case nameIn[1]==S
            case 'T' :
              if (nameIn[2] == 'H' && nameIn[3] == 'R') return PROTEIN;
              if (nameIn[2] == 'R' && nameIn[3] == 'P') return PROTEIN;
              if (nameIn[2] == 'Y' && nameIn[3] == 'R') return PROTEIN;
            break; // END case nameIn[1]==T
            case 'V' :
              if (nameIn[2] == 'A' && nameIn[3] == 'L') return PROTEIN;
            break; // END case nameIn[1]==V
          } // END switch nameIn[1]
        break; // END case nameIn[0]==C

        case 'N' :
          switch (nameIn[1]) {
            case 'A' :
              if (nameIn[2] == 'L' && nameIn[3] == 'A') return PROTEIN; // N-Alanine
              if (nameIn[2] == 'R' && nameIn[3] == 'G') return PROTEIN; // N-Arginine
              if (nameIn[2] == 'S' && nameIn[3] == 'N') return PROTEIN; // N-Asparagine
              if (nameIn[2] == 'S' && nameIn[3] == 'P') return PROTEIN; // N-Aspartic acid
            break; // END case nameIn[1]==A
            case 'C' :
              if (nameIn[2] == 'Y' && nameIn[3] == 'S') return PROTEIN; // Protonated (normal) N-CYS
              if (nameIn[2] == 'Y' && nameIn[3] == 'X') return PROTEIN; // Disulphide N-CYS
            break; // END case nameIn[1]==C
            case 'G' :
              if (nameIn[2] == 'L' && nameIn[3] == 'N') return PROTEIN; // N-Glutamine
              if (nameIn[2] == 'L' && nameIn[3] == 'U') return PROTEIN; // N-Glutamic acid
              if (nameIn[2] == 'L' && nameIn[3] == 'Y') return PROTEIN; // N-Glycine
            break; // END case nameIn[1]==G
            case 'H' :
              if (nameIn[2] == 'I' && nameIn[3] == 'D') return PROTEIN;
              if (nameIn[2] == 'I' && nameIn[3] == 'E') return PROTEIN;
              if (nameIn[2] == 'I' && nameIn[3] == 'P') return PROTEIN;
              if (nameIn[2] == 'I' && nameIn[3] == 'S') return PROTEIN;
            break; // END case nameIn[1]==H
            case 'I' :
              if (nameIn[2] == 'L' && nameIn[3] == 'E') return PROTEIN;
            break; // END case nameIn[1]==I
            case 'L' :
              if (nameIn[2] == 'E' && nameIn[3] == 'U') return PROTEIN;
              if (nameIn[2] == 'Y' && nameIn[3] == 'S') return PROTEIN;
            break; // END case nameIn[1]==L
            case 'M' :
              if (nameIn[2] == 'E' && nameIn[3] == 'T') return PROTEIN;
            break; // END case nameIn[1]==M
            case 'P' :
              if (nameIn[2] == 'H' && nameIn[3] == 'E') return PROTEIN;
              if (nameIn[2] == 'R' && nameIn[3] == 'O') return PROTEIN;
            break; // END case nameIn[1]==P
            case 'S' :
              if (nameIn[2] == 'E' && nameIn[3] == 'R') return PROTEIN;
            break; // END case nameIn[1]==S
            case 'T' :
              if (nameIn[2] == 'H' && nameIn[3] == 'R') return PROTEIN;
              if (nameIn[2] == 'R' && nameIn[3] == 'P') return PROTEIN;
              if (nameIn[2] == 'Y' && nameIn[3] == 'R') return PROTEIN;
            break; // END case nameIn[1]==T
            case 'V' :
              if (nameIn[2] == 'A' && nameIn[3] == 'L') return PROTEIN;
            break; // END case nameIn[1]==V
          } // END switch nameIn[1]
        break; // END case nameIn[0]==N

        case 'T' :
          switch (nameIn[1]) {
            case 'I' :
              if (nameIn[2] == 'P' && nameIn[3] == '3') return SOLVENT; // TIP3P water
              if (nameIn[2] == 'P' && nameIn[3] == '4') return SOLVENT; // TIP4P water
              if (nameIn[2] == 'P' && nameIn[3] == '5') return SOLVENT; // TIP5P water
            break; // END case nameIn[1]==I
          } // END switch nameIn[1]
        break; // END case nameIn[0]==T

      } // END switch nameIn[0]
    break; // END case 4

    case 3 :
      // 3 letters. Could be lots of things
      switch (nameIn[0]) {
        case 'A' :
          if (nameIn[1] == 'C' && nameIn[2] == 'E') return PROTEIN; // Acetyl
          if (nameIn[1] == 'L' && nameIn[2] == 'A') return PROTEIN; // Alanine
          if (nameIn[1] == 'R' && nameIn[2] == 'G') return PROTEIN; // Arginine
          if (nameIn[1] == 'S' && nameIn[2] == 'H') return PROTEIN; // Protonated ASP
          if (nameIn[1] == 'S' && nameIn[2] == 'N') return PROTEIN; // Asparagine
          if (nameIn[1] == 'S' && nameIn[2] == 'P') return PROTEIN; // Aspartic acid
        break; // END case nameIn[0]==A
        case 'C' :
          if (nameIn[1] == 'Y' && nameIn[2] == 'M') return PROTEIN; // Deprotonated CYS
          if (nameIn[1] == 'Y' && nameIn[2] == 'S') return PROTEIN; // Protonated (normal) CYS
          if (nameIn[1] == 'Y' && nameIn[2] == 'X') return PROTEIN; // Disulphide CYS
        break; // END case nameIn[0]==C
        case 'G' :
          if (nameIn[1] == 'L' && nameIn[2] == 'H') return PROTEIN; // Protonated GLU
          if (nameIn[1] == 'L' && nameIn[2] == 'N') return PROTEIN; // Glutamine
          if (nameIn[1] == 'L' && nameIn[2] == 'U') return PROTEIN; // Glutamic acid
          if (nameIn[1] == 'L' && nameIn[2] == 'Y') return PROTEIN; // Glycine
        break; // END case nameIn[0]==G
        case 'H' :
          if (nameIn[1] == 'I' && nameIn[2] == 'D') return PROTEIN;
          if (nameIn[1] == 'I' && nameIn[2] == 'E') return PROTEIN;
          if (nameIn[1] == 'I' && nameIn[2] == 'P') return PROTEIN;
          if (nameIn[1] == 'I' && nameIn[2] == 'S') return PROTEIN;
          if (nameIn[1] == 'O' && nameIn[2] == 'H') return SOLVENT; // PDB water 
          if (nameIn[1] == 'Y' && nameIn[2] == 'P') return PROTEIN;
        break; // END case nameIn[0]==H
        case 'I' :
          if (nameIn[1] == 'L' && nameIn[2] == 'E') return PROTEIN;
        break; // END case nameIn[0]==I
        case 'L' :
          if (nameIn[1] == 'E' && nameIn[2] == 'U') return PROTEIN;
          if (nameIn[1] == 'Y' && nameIn[2] == 'N') return PROTEIN;
          if (nameIn[1] == 'Y' && nameIn[2] == 'S') return PROTEIN;
        break; // END case nameIn[0]==L
        case 'M' :
          if (nameIn[1] == 'E' && nameIn[2] == 'T') return PROTEIN;
        break; // END case nameIn[0]==M
        case 'N' :
          if (nameIn[1] == 'H' && nameIn[2] == 'E') return PROTEIN;
          if (nameIn[1] == 'M' && nameIn[2] == 'E') return PROTEIN;
        case 'P' :
          if (nameIn[1] == 'H' && nameIn[2] == 'E') return PROTEIN;
          if (nameIn[1] == 'R' && nameIn[2] == 'O') return PROTEIN;
        break; // END case nameIn[0]==P
        case 'S' :
          if (nameIn[1] == 'E' && nameIn[2] == 'R') return PROTEIN;
          if (nameIn[1] == 'O' && nameIn[2] == 'L') return SOLVENT; // CHARMM water
        break; // END case nameIn[0]==S
        case 'T' :
          if (nameIn[1] == 'H' && nameIn[2] == 'R') return PROTEIN;
          if (nameIn[1] == 'R' && nameIn[2] == 'P') return PROTEIN;
          if (nameIn[1] == 'Y' && nameIn[2] == 'R') return PROTEIN;
        break; // END case nameIn[0]==T
        case 'V' :
          if (nameIn[1] == 'A' && nameIn[2] == 'L') return PROTEIN;
        break; // END case nameIn[0]==V
        case 'W' :
          if (nameIn[1] == 'A' && nameIn[2] == 'T') return SOLVENT;
        break; // END case nameIn[0]==W

      } // END switch nameIn[0]
    break; // END case 3

    case 2 :
      // 2 letters.
      switch (nameIn[0]) {
        case 'A' :
          if (nameIn[1] == '3') return NUCLEIC;
          if (nameIn[1] == '5') return NUCLEIC;
        break; // END case nameIn[0]==A
        case 'D' :
          if (nameIn[1] == 'A') return NUCLEIC;
          if (nameIn[1] == 'C') return NUCLEIC;
          if (nameIn[1] == 'G') return NUCLEIC;
          if (nameIn[1] == 'T') return NUCLEIC;
        break; // END case nameIn[0]==D
      } // END switch nameIn[0]
    break; // END case 2

  } // END switch
  return UNKNOWN;
}

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
