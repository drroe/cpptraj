#ifndef INC_SOLVATE_H
#define INC_SOLVATE_H
#include <string>
class ArgList;
class DataSet_Coords;
class DataSetList;
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
    int SolvateBox(Topology&, Frame&, Cpptraj::Parm::ParameterSet const&, DataSet_Coords&) const;

    /// \return Solvent unit selected from given DataSetList
    DataSet_Coords* GetSolventUnit(DataSetList const&) const;
    //std::string const& SolventBoxName() const { return solventBoxName_; }
  private:
    // Set vdW bounding box
    int setVdwBoundingBox(double&, double&, double&, Topology const&, Frame&,
                          Cpptraj::Parm::ParameterSet const&) const;
    // Add solvent unit boxes
    int addSolventUnits(int, int, int, double, double, double, double, double, double,
                        Frame&, Topology const&, Frame&, Topology&) const;

    int debug_;
    double bufferX_;
    double bufferY_;
    double bufferZ_;
    bool isotropic_;
    bool clip_;
    std::string solventBoxName_;
    static const double ATOM_DEFAULT_RADIUS_; ///< Atom default radius from LEaP
};
}
}
#endif
