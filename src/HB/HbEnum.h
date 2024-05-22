#ifndef INC_HB_HBENUM_H
#define INC_HB_HBENUM_H
namespace Cpptraj {
namespace HB {
/// Heavy atom site types
enum SiteType { DONOR=0, ACCEPTOR, BOTH, UNKNOWN_SITE };
/// Strings corresponding to atom site types
const char* SiteTypeStr(SiteType);
/// Heavy atom hydrogen bond types
enum HbType { SOLUTE=0, SOLVENT, ION, UNKNOWN_HB };
/// Strings corresponding to hydrogen bond types
const char* HbTypeStr(HbType);
}
}
#endif
