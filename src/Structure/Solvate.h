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
    // Set vdW bounding box
    int setVdwBoundingBox(double&, double&, double&, Topology const&, Frame&,
                          Cpptraj::Parm::ParameterSet const&) const;

    int debug_;
    double bufferX_;
    double bufferY_;
    double bufferZ_;
    bool isotropic_;
    static const double ATOM_DEFAULT_RADIUS_; ///< Atom default radius from LEaP
};
}
}
#endif
