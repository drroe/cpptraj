#ifndef INC_MEAD_MEADERROR_H
#define INC_MEAD_MEADERROR_H
// MEAD fwd declares
class MEADexcept;
namespace Cpptraj {
namespace Mead {
/// Used to print MEAD error message
int ERR(const char*, MEADexcept&);

}
}
#endif
