#ifndef INC_SOLVATE_H
#define INC_SOLVATE_H
#include <string>
#include <vector>
class ArgList;
class DataSet_Coords;
class DataSetList;
class Frame;
class Topology;
class Vec3;
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
    int SolvateBox(Topology&, Frame&, Cpptraj::Parm::ParameterSet const&, DataSet_Coords&);

    /// \return Solvent unit selected from given DataSetList
    DataSet_Coords* GetSolventUnit(DataSetList const&) const;
    //std::string const& SolventBoxName() const { return solventBoxName_; }
  private:
    /// Set vdW bounding box
    int setVdwBoundingBox(double&, double&, double&, double&, std::vector<double>&,
                          Topology const&, Frame&,
                          Cpptraj::Parm::ParameterSet const&) const;
    /// Find solute atoms within a solvent box at given center
    int findCloseSoluteAtoms(std::vector<int>&, double, int, Frame const&, Vec3 const&, double, double, double) const;
    /// Determine which solvent residues do not clash with given solute atoms
    int determineValidSolventResidues(std::vector<int>&, std::vector<int> const&,
                                      Frame const&, Topology const&, Frame const&,
                                      std::vector<double> const&, std::vector<double> const&) const;

    // Add solvent unit boxes
    int addSolventUnits(int, int, int, double, double, double, double, double, double, double,
                        Frame&, Topology const&, Frame&, Topology&,
                        std::vector<double> const&, std::vector<double> const&) const;

    int debug_;
    double bufferX_;
    double bufferY_;
    double bufferZ_;
    double closeness_;
    double clipX_;
    double clipY_;
    double clipZ_;
    bool isotropic_;
    bool clip_;
    std::string solventBoxName_;
    static const double ATOM_DEFAULT_RADIUS_; ///< Atom default radius from LEaP
    static const double CLOSENESSMODIFIER_;   ///< Overlap closeness modifier from LEaP
};
}
}
#endif
