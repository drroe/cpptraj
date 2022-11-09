#ifndef INC_PARAM_LES_ATOMTYPE_H
#define INC_PARAM_LES_ATOMTYPE_H
#include <vector>
namespace Cpptraj {
namespace Param {
/// Hold LES (locally enhanced sampling) atom parameters.
/** In LES, copies of certain regions are created with the
  * atoms duplicated. LES parameters describe how the duplicated
  * atoms interact with everything else.
  */
class LES_AtomType {
  public:
    /// CONSTRUCTOR
    LES_AtomType() : type_(0), cnum_(0), id_(0) {}
    /// CONSTRUCTOR - LES type, LES copy #, LES region ID
    LES_AtomType(int t, int c, int i) : type_(t), cnum_(c), id_(i) {}
    /// \return LES type for the atom
    int Type() const { return type_; }
    /// \return LES copy # for the atom. 0 == atom exists for all copies.
    int Copy() const { return cnum_; }
    /// \return LES region ID for the atom
    int ID()   const { return id_;   }
    /// Set the LES atom type
    void SetType(int t) { type_ = t; }
    /// Set the LES copy #
    void SetCopy(int c) { cnum_ = c; }
    /// Set the LES region ID
    void SetID(int i)   { id_ = i;   }
  private:
    int type_; ///< LES atom type
    int cnum_; ///< LES copy #
    int id_;   ///< LES region ID
};
/// Hold array of LES parameters
typedef std::vector<LES_AtomType> LES_Array;
}
}
#endif
