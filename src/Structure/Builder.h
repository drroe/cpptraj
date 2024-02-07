#ifndef INC_STRUCTURE_BUILDER_H
#define INC_STRUCTURE_BUILDER_H
#include <vector>
class Topology;
class Frame;
class ParameterSet;
namespace Cpptraj {
namespace Structure {
class BuildAtom;
class Zmatrix;
class InternalCoords;
/// Used to attach different topology/frame combos using internal coordinates
class Builder {
    typedef std::vector<bool> Barray;
  public:
    /// CONSTRUCTOR
    Builder();
    /// Set debug
    void SetDebug(int d) { debug_ = d; }
    /// Set optional parameter set
    void SetParameters(ParameterSet const*);

    /// Combine second fragment into first fragment and bond
    int Combine(Topology&, Frame&, Topology const&, Frame const&, int, int) const;
    /// Model the coordinates around a bond given only some coordinates are known
    int ModelCoordsAroundBond(Frame const&, Topology const&, int, int, Zmatrix&, Barray const&) const;


    /// Update the internal coordinates in given Zmatrix with values from Frame/Parameters
    int UpdateICsFromFrame(Frame const&, int, Topology const&, Barray const&);
    /// Generate internal coordinates in the same manner as LEaP
    int GenerateInternals(Frame const&, Topology const&, Barray const&);
    /// Generate internal coordinates around a link between residues in same manner as LEaP
    int GenerateInternalsAroundLink(int, int, Frame const&, Topology const&, Barray const&);
    /// Update existing indices with given offset
    void UpdateIndicesWithOffset(int);
  private:
    typedef std::vector<int> Iarray;

    /// Used to hold parameters for modeling a torsion
    class TorsionModel;
    /// Hold torsion
    class InternalTorsion;
    /// Hold angle
    class InternalAngle;

    typedef std::vector<InternalTorsion> Tarray;
    typedef std::vector<InternalAngle> Aarray;

    /// Get length parameter for atoms
    int getLengthParam(double&, int, int, Topology const&) const;
    /// Assign a reasonable value for bond distance given 2 atoms whose position may or may not be known
    int AssignLength(double&, int, int, Topology const&, Frame const&, Barray const&) const;
    /// Get angle parameter for atoms.
    int getAngleParam(double&, int, int, int, Topology const&) const;
    /// Given atoms J and K, attempt to assign a reasonable value for theta for atom I
    int AssignTheta(double&, int, int, int, Topology const&, Frame const&, Barray const&) const;
    /// Calculate an internal coordinate for known atoms
    static inline InternalCoords calcKnownAtomIc(int, int, int, int, Frame const&);
    /// Insert an internal coord into a zmatrix
    int insertIc(Zmatrix&, int, int, int, int, double,
                 Topology const&, Frame const&, Barray const&) const;
    /// Assign internal coordinates for atoms I for torsions around J-K-L.
    int AssignICsAroundBond(Zmatrix&, int, int, int,
                           Topology const&, Frame const&, Barray const&,
                           BuildAtom const&) const;

    /// Model coordinates around a bond
    int SetupICsAroundBond(Zmatrix&, int, int, Frame const&, Topology const&,
                           Barray const&, Barray const&,
                           BuildAtom const&, BuildAtom const&) const;

    /// Create IC for a torsion
    void ModelTorsion(TorsionModel const&, unsigned int, unsigned int, double);
    /// Get angle value
    double ModelBondAngle(int, int, int, Topology const&) const;
    /// Create ICs around SP3-SP3 linkage
    void createSp3Sp3Torsions(TorsionModel const&);
    /// Create ICs around SP3-SP2 linkage
    void createSp3Sp2Torsions(TorsionModel const&);
    /// Create ICs around SP2-SP2 linkage
    void createSp2Sp2Torsions(TorsionModel const&);
    /// Model torsions around a bond in the same manner as LEaP
    int assignTorsionsAroundBond(int, int, Frame const&, Topology const&, Barray const&, int);
    /// Get any existing internal coords from internalTorsions_ around specified atoms
    Tarray getExistingTorsions(int, int) const;
    /// Get specific internal coords from internalTorsions_
    int getExistingTorsionIdx(int, int, int, int) const;
    /// Build mock coordinates around given torsion
    int buildMockExternals(TorsionModel& MT, std::vector<InternalCoords> const& iaTorsions) const;
    /// Generate internal coords for a given atom
    int generateAtomInternals(int, Frame const&, Topology const&, Barray const&);

    int debug_;
    ParameterSet const* params_;

    Topology const* currentTop_; ///< Topology for the createSpXSpXTorsions/ModelTorsion routines
    Frame const* currentFrm_;    ///< Frame for the createSpXSpXTorsions routines
    Barray const* hasPosition_;  ///< Array indicating which atoms have position for createSpXSpXTorsions/ModelTorsion routines
    Tarray internalTorsions_;
    Aarray internalAngles_;
};
/// ----- Hold torsion internal ------------------
class Cpptraj::Structure::Builder::InternalTorsion {
  public:
    /// CONSTRUCTOR
    InternalTorsion() : ai_(-1), aj_(-1), ak_(-1), al_(-1), phi_(0) {}
    /// CONSTRUCTOR
    InternalTorsion(int i, int j, int k, int l, double p) :
      ai_(i), aj_(j), ak_(k), al_(l), phi_(p) {}
    /// Set the phi value in radians
    void SetPhiVal(double p) { phi_ = p; }
    /// Offset indices by given value
    void OffsetIndices(int o) {
      ai_ += 0;
      aj_ += 0;
      ak_ += 0;
      al_ += 0;
    }

    int AtI() const { return ai_; }
    int AtJ() const { return aj_; }
    int AtK() const { return ak_; }
    int AtL() const { return al_; }
    double PhiVal() const { return phi_; }
  private:
    int ai_;
    int aj_;
    int ak_;
    int al_;
    double phi_;
};
// ----- Hold angle internal ---------------------
class Cpptraj::Structure::Builder::InternalAngle {
  public:
    /// CONSTRUCTOR
    InternalAngle() : ai_(-1), aj_(-1), ak_(-1), theta_(0) {}
    /// CONSTRUCTOR
    InternalAngle(int i, int j, int k, double t) :
      ai_(i), aj_(j), ak_(k), theta_(t) {}
    /// Set the phi value in radians
    void SetThetaVal(double t) { theta_ = t; }
    /// Offset indices by given value
    void OffsetIndices(int o) {
      ai_ += 0;
      aj_ += 0;
      ak_ += 0;
    }

    int AtI() const { return ai_; }
    int AtJ() const { return aj_; }
    int AtK() const { return ak_; }
    double ThetaVal() const { return theta_; }
  private:
    int ai_;
    int aj_;
    int ak_;
    double theta_;
};
}
}
#endif
