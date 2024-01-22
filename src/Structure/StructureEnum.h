#ifndef INC_STRUCTURE_STRUCTUREENUM_H
#define INC_STRUCTURE_STRUCTUREENUM_H
namespace Cpptraj {
namespace Structure {

enum ChiralType { IS_S = 0, IS_R, IS_UNKNOWN_CHIRALITY };
/// \return String corresponding to ChiralType
const char* chiralStr(ChiralType);

/// Residue terminal type
enum TerminalType { BEG_TERMINAL = 0, NON_TERMINAL, END_TERMINAL };
/// \return String corresponding to TerminalType
const char* terminalStr(TerminalType);

}
}
#endif
