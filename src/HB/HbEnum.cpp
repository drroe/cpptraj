#include "HbEnum.h"

const char* Cpptraj::HB::TypeStr(Type siteType) {
  switch (siteType) {
    case DONOR     : return "Solute Donor";
    case ACCEPTOR  : return "Solute Acceptor";
    case BOTH      : return "Solute Both";
    case VDONOR    : return "Solvent Donor";
    case VACCEPTOR : return "Solvent Acceptor";
    case VBOTH     : return "Solvent Both";
    case UNKNOWN   : return "Unknown";
  }
  return 0;
}

/*const char* Cpptraj::HB::SiteTypeStr(SiteType siteType) {
  switch (siteType) {
    case DONOR : return "Donor";
    case ACCEPTOR : return "Acceptor";
    case BOTH : return "Both";
    case UNKNOWN_SITE : break;
  }
  return 0;
}

const char* Cpptraj::HB::HbTypeStr(HbType hbtype) {
  switch (hbtype) {
    case SOLUTE : return "Solute";
    case SOLVENT : return "Solvent";
    case ION : return "Ion" ;
    case UNKNOWN_HB : break;
  }
  return 0;
}*/
