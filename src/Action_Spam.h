#ifndef INC_ACTION_SPAM_H
#define INC_ACTION_SPAM_H
#include "Action.h"
#include "ImageOption.h"
#include "Vec3.h"
#include "Timer.h"
#include "PairList.h"
// Forward declares
class DataSet_Vector_Scalar;
/**
SPAM is a water profiling technique developed by Guanglei Cui at
GlaxoSmithKline (GSK). The original implementation involved a set of specialized
Python scripts interacting with VMD (via the VolMap tool), numpy, NAMD (for the
SPAM energy calculations) and R (for the free energy calculation using a
specialized kernel density estimate). While that implementation demonstrated
proof of principle, simply re-ordering the trajectory for use with NAMD proved
to be a performance bottleneck because it was written in Python. SPAM was
rewritten from the ground up in cpptraj, significantly improving efficiency and
providing a simpler interface.

The original C++ implementation of SPAM in cpptraj was done by Jason Swails
while interning at GSK. This code was built as a patch on top of cpptraj v.12
and was rewritten by Jason Swails for the current cpptraj version.

 (C) 2012 - 2013
*/
class Action_Spam: public Action {
  public:
    Action_Spam();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Spam(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif

    /// Hold info for a solvent type
    class SolventInfo;
    /// Hold info for a solvent within a peak site
    class SolventPeak;
    /// Hold information for a solvent peak site
    class PeakSite;

    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> Parray; ///< Peak array type
    typedef std::vector<double> Darray;
    typedef std::vector<Vec3> Varray;
    typedef std::vector<Residue> Rarray;
    typedef std::vector<DataSet*> DSarray;

    // ------------------- Functions -------------------
    int SetupParms(Topology const&);
    double Calculate_Energy(Frame const&, Residue const&);
    int Calc_G_Wat(DataSet*, unsigned int);
    // Custom Do- routines
    Action::RetType DoPureWater(int, Frame const&);
    Action::RetType DoSPAM(int, Frame&);

    DataSet_Vector_Scalar* GetPeaksData(std::string const&, DataSetList const&);
    typedef bool (Action_Spam::*FxnType)(Vec3, Vec3, double) const;
    bool inside_box(Vec3, Vec3, double) const;
    bool inside_sphere(Vec3, Vec3, double) const;
    inline double Ecalc(int, int, double) const;

    int debug_;
    FxnType Inside_;          ///< Function for determining if water is inside peak.
    ImageOption imageOpt_;    ///< Used to determine if imaging should be used.
    PairList pairList_;       ///< Atom pair list (purewater_ only)
    Iarray watidx_;           ///< Hold water index for each atom (starting from 0).
    std::string solvname_;    ///< Name of the solvent residues
    double DG_BULK_;          ///< SPAM free energy of the bulk solvent
    double DH_BULK_;          ///< SPAM enthalpy of the bulk solvent
    double temperature_;      ///< Temperature at which SPAM simulation was run
    bool purewater_;          ///< True if running a pure water simulation to derive bulk properties
    bool reorder_;            ///< True if solvent should be reordered
    bool calcEnergy_;         ///< True if energy needs to be calculated.
    double cut2_;             ///< Non-bonded cutoff in Angstroms (squared)
    double onecut2_;          ///< 1 / cut2_
    double doublecut_;        ///< twice the cutoff (to test if boxes are big enough)
    CpptrajFile* infofile_;   ///< SPAM info file
    AtomMask mask_;           ///< Mask for selecting individual solvent residues
    Iarray resPeakNum_;       ///< Peak that each solvent residue is assigned to; -1 is unassigned
    double site_size_;        ///< Size of the water site. This is a full edge length or diameter
    Topology* CurrentParm_;   ///< Current topology (for NB params).
    Darray atom_charge_;      ///< Charges that have been converted to Amber units
    bool sphere_;             ///< Is our site shape a sphere? If no, it's a box.
    DataSet* ds_dg_;          ///< Hold final delta G values for each peak
    DataSet* ds_dh_;          ///< Hold final delta H values for each peak
    DataSet* ds_ds_;          ///< Hold final -T*S values for each peak
    Parray peakFrameData_;    ///< A list of all omitted frames for each peak
    DSarray myDSL_;           ///< Hold energy data sets
    Varray comlist_;          ///< For given frame, each residue C.O.M. coords.
    Rarray solvent_residues_; ///< List of each solvent residue
    int Nframes_;             ///< Total number of frames
    bool overflow_;           ///< True if cutoff overflowed our box coordinates
    DataSetList peaksdsl_;    ///< Will allocate DataSet for peaks data if loading from a file.
    DataSet_Vector_Scalar* peaksData_; ///< Hold peaks DataSet

    std::vector<SolventInfo> solvents_; ///< Hold info for each solvent type
    std::vector<PeakSite> peakSites_;   ///< Hold info for every solvent peak
    // Timers
    Timer t_action_;
    Timer t_resCom_;
    Timer t_assign_;
    Timer t_occupy_;
    Timer t_energy_;
    Timer t_reordr_;
};

// ----- SolventInfo class -----------------------------------------------------
/** Hold information for a specific solvent type. */
class Action_Spam::SolventInfo {
  public:
    SolventInfo();
    /// Construct with peaks data, size size, name
    SolventInfo(DataSet_Vector_Scalar const*, double, std::string const&);
    /// Print info to stdout
    void PrintInfo() const;
  private:
    DataSet_Vector_Scalar const* peaksData_; ///< Hold peaks DataSet for this solvent.
    double site_size_;                       ///< Size of solvent site (Ang.). Full edge length or diameter
    std::string name_;                       ///< Solvent residue name.
    Iarray resIdxs_;                         ///< Solvent residue indices.
};

// ----- SolventPeak class -----------------------------------------------------
/** Hold information for specific solvent occupying a site. */
class Action_Spam::SolventPeak {
  public:
    SolventPeak();
  private:
    DataSet* energies_; ///< Hold solvent energies for this peak.
    Iarray ommitted_;   ///< Hold info on frames for which no solvent energies calcd.
};

// ----- PeakSite class --------------------------------------------------------
/** Hold all information related to a solvent peak site. */
class Action_Spam::PeakSite {
  public:
    PeakSite();
    /// Construct from given peak position.
    PeakSite(Vec3 const&);
    /// \return XYZ coords of peak location
    Vec3 const& XYZ() const { return xyz_; }
  private:
    typedef std::vector<SolventPeak> SolvPeakArray;
    Vec3 xyz_;                ///< Solvent peak location in Cartesian space.
    SolvPeakArray solvPeaks_; ///< Hold information for each solvent that might occupy this site.
};
#endif
