#ifndef INC_SOLVATE_H
#define INC_SOLVATE_H
class ArgList;
class Frame;
class Topology;
namespace Cpptraj {
namespace Parm {
class ParameterSet;
}
namespace Structure {
/// Used to add solvent to a frame/topology pair.
class Solvate {
  public:
    /// CONSTRUCTOR
    Solvate();
    /// Initialize
    int InitSolvate(ArgList&, int);
    /// Solvate with box
    int SolvateBox(Topology&, Frame&, Cpptraj::Parm::ParameterSet const&); 
  private:
    int debug_;
    static const double ATOM_DEFAULT_RADIUS_; ///< Atom default radius from LEaP
};
}
}
#endif
