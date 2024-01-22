#ifndef INC_STRUCTURE_BUILDATOM_H
#define INC_STRUCTURE_BUILDATOM_H
#include "StructureEnum.h"
#include <vector>
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
/// Hold chirality/orientation information for an atom. Used when building/modelling new coordinates.
class BuildAtom {
  public:
    /// Blank constructor 
    BuildAtom() : tors_(0), ctype_(IS_UNKNOWN_CHIRALITY), orientation_(IS_UNKNOWN_CHIRALITY) {}
    /// CONSTRUCTOR - Take only chiral type. Usually used to indicate an error
    BuildAtom(ChiralType e) : tors_(0), ctype_(e), orientation_(e) {}
    /// CONSTRUCTOR - torsion, chirality, orientation, priority array
    BuildAtom(double t, ChiralType c, ChiralType o, std::vector<int> const& p) :
                  tors_(t), ctype_(c), orientation_(o), priority_(p) {}
    /// Set chirality, orientation, and priority around specified atom
    int DetermineChirality(int, Topology const&, Frame const&, int);
    /// Only set the priority around specified atom
    int SetPriority(int, Topology const&, Frame const&, int);

    /// \return Atom chirality
    ChiralType Chirality()             const { return ctype_; }
    /// \return Orientation around atom
    ChiralType Orientation()           const { return orientation_; }
    /// \return Priority array
    std::vector<int> const& Priority() const { return priority_; }
    /// \return Value of the orientation torsion around atom (in radians).
    double TorsionVal()                const { return tors_; }
    /// \return True if an error occurred determining chirality
    bool ChiralError()                 const { return (ctype_ == CHIRALITY_ERR);}
  private:
    /// Total priority (i.e. sum of atomic numbers) of atoms bonded to given atom.
    static int totalPriority(Topology const&, int, int, int, int, std::vector<bool>&);
    /// Determine priority and optionally chirality as well
    int determineChirality(int, Topology const&, Frame const&, int, bool);

    double tors_;               ///< Torsion around the atom in radians (via priority, 1-2-3-0, where 0 is this atom).
    ChiralType ctype_;          ///< Chirality around atom.
    ChiralType orientation_;    ///< Orientation when chirality cannot be determined.
    std::vector<int> priority_; ///< Indices of bonded atoms, sorted by priority.
};
}
}
#endif
