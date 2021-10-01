#ifndef INC_ACTION_SPAM_H
#define INC_ACTION_SPAM_H
#include "Action.h"
#include "ImageOption.h"
#include "Vec3.h"
#include "Timer.h"
#include "PairList.h"
// Forward declares
class DataSet_Vector_Scalar;
class DataSet_double;
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

    /// Hold info for solvent residue
    class SolventRes;
    /// Hold info for a solvent type
    class SolventInfo;
    /// Hold info for a solvent within a peak site
    class SolventPeak;
    /// Hold information for a solvent peak site
    class PeakSite;

    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;

    // ------------------- Functions -------------------
    DataSet_Vector_Scalar* GetPeaksData(std::string const&, DataSetList const&);
    int SetupParms(Topology const&);
    RetType DoPureWater(int, Frame const&);
    int Peaks_Ene_Calc(Iarray const&, Iarray const&, Frame const&, int);
    RetType SpamCalc(int, Frame&);
    int Calc_G(double&, int, double, double, double, DataSet_double const&) const;
    int Calc_Bulk() const;
    int Calc_G_Peak(unsigned int, PeakSite const&) const;

    typedef bool (Action_Spam::*FxnType)(Vec3 const&, Vec3 const&, double) const;
    bool inside_box(Vec3 const&, Vec3 const&, double) const;
    bool inside_sphere(Vec3 const&, Vec3 const&, double) const;

    inline double Ecalc(int, int, double) const;

    int debug_;
    FxnType Inside_;          ///< Function for determining if water is inside peak.
    ImageOption imageOpt_;    ///< Used to determine if imaging should be used.
    PairList pairList_;       ///< Atom pair list for energy calculations. 
    double DG_BULK_;          ///< SPAM free energy of the bulk solvent
    double DH_BULK_;          ///< SPAM enthalpy of the bulk solvent
    double temperature_;      ///< Temperature at which SPAM simulation was run
    bool purewater_;          ///< True if running a pure water simulation to derive bulk properties
    bool reorder_;            ///< True if solvent should be reordered
    bool calcEnergy_;         ///< True if energy needs to be calculated.
    bool printFrameInfo_;     ///< If true print details for omitted frames.
    double cut2_;             ///< Non-bonded cutoff in Angstroms (squared)
    double onecut2_;          ///< 1 / cut2_ (for truncation scheme)
    double doublecut_;        ///< twice the cutoff (to test if boxes are big enough)
    CpptrajFile* infofile_;   ///< SPAM info file
    AtomMask mask_;           ///< Mask for selecting atoms for pair list 
    Iarray resPeakNum_;       ///< Peak that each solvent residue is assigned to each frame; -1 is unassigned
    Iarray watidx_;           ///< For each atom in mask_, solvResArray index (purewater only)
    Topology* CurrentParm_;   ///< Current topology (for NB params).
    Darray atom_charge_;      ///< Charges that have been converted to Amber units
    bool sphere_;             ///< Is our site shape a sphere? If no, it's a box.
    DataSet* bulk_ene_set_;   ///< Hold bulk solvent energies (purewater_)

    int Nframes_;             ///< Total number of frames
    bool overflow_;           ///< True if cutoff overflowed our box coordinates
    DataSetList peaksdsl_;    ///< Will allocate DataSet for peaks data if loading from a file.

    std::vector<SolventInfo> solvents_;    ///< Hold info for each solvent type
    std::vector<PeakSite> peakSites_;      ///< Hold info for every solvent peak
    std::vector<SolventRes> solvResArray_; ///< Hold every solvent residue
    // Timers
    Timer t_action_;
    Timer t_resCom_;
    Timer t_assign_;
    Timer t_occupy_;
    Timer t_energy_;
    Timer t_reordr_;
#   ifdef _OPENMP
    std::vector<Darray> threadResEne_; ///< For each thread, hold 'purewater' ene for each res
#   endif
};

// ----- SolventRes class ------------------------------------------------------
/** Hold info for a solvent residue */
class Action_Spam::SolventRes {
  public:
    SolventRes();
    /// Construct with res first atom, last atom, and index into solvents_ array
    SolventRes(int, int, int);
    /// Print solvent res info to stdout
    void PrintInfo() const;
    /// \return First atom of the solvent residue.
    int At0() const { return at0_; }
    /// \return (one after the) Final atom of the solvent residue.
    int At1() const { return at1_; }
    /// \return Index into solvents_ array for this solvent residue.
    int Sidx() const { return sidx_; }
  private:
    int at0_;  ///< Residue first atom
    int at1_;  ///< Residue last atom
    int sidx_; ///< Index into the solvents_ array
};

// ----- SolventInfo class -----------------------------------------------------
/** Hold information for a specific solvent type. */
class Action_Spam::SolventInfo {
  public:
    SolventInfo();
    /// Construct with solvent name (bulk)
    SolventInfo(std::string const&);
    /// Construct with peaks data, size size, name
    SolventInfo(DataSet_Vector_Scalar const*, double, std::string const&);
    /// Create total delta E sets for solvent
    int CreateDeltaEneSets(std::string const&, int, DataSetList&, DataFileList&, DataFile*);
    /// Print info to stdout
    void PrintInfo() const;
    /// \return solvent residue name
    std::string const& Name() const { return name_; }
    /// \return Solvent site size
    double SiteSize() const { return site_size_; }
    /// \return Pointer to delta G of solvent for each peak set
    DataSet* DG() const { return ds_dg_; }
    /// \return Pointer to delta H of solvent for each peak set
    DataSet* DH() const { return ds_dh_; }
    /// \return Pointer to -T * delta S of solvent for each peak set
    DataSet* TDS() const { return ds_ds_; }
    /// \return Pointer to DataFile that solvent sets will be written to
    DataFile* SolventFile() const { return sfile_; }
  private:
    DataSet_Vector_Scalar const* peaksData_; ///< Hold peaks DataSet for this solvent.
    double site_size_;                       ///< Size of solvent site (Ang.). Full edge length or diameter
    std::string name_;                       ///< Solvent residue name.
    DataSet* ds_dg_;                         ///< Solvent Delta G for each peak
    DataSet* ds_dh_;                         ///< Solvent Delta H for each peak
    DataSet* ds_ds_;                         ///< Solvent -T * Delta S for each peak
    DataFile* sfile_;                        ///< File that solvent sets will be written to (for info only).
    //Iarray resIdxs_;                         ///< Solvent residue indices. TODO needed?
};

// ----- SolventPeak class -----------------------------------------------------
/** Hold information for specific solvent occupying a site. */
class Action_Spam::SolventPeak {
  public:
    SolventPeak();
    /// Construct with given energy DataSet
    SolventPeak(DataSet*);
    /// Add an omitted frame number to omitted_
    void AddOmittedFrameNum(int fn) { omitted_.push_back( fn ); }
    /// \return the energy DataSet pointer
    DataSet* DS() const { return energies_; }
    /// \return Omitted frames array
    Iarray const& Omitted() const { return omitted_; }
  private:
    DataSet* energies_; ///< Hold solvent energies for this peak.
    Iarray omitted_;   ///< Hold info on frames for which no solvent energies calcd.
};

// ----- PeakSite class --------------------------------------------------------
/** Hold all information related to a solvent peak site. */
class Action_Spam::PeakSite {
    typedef std::vector<SolventPeak> SolvPeakArray;
  public:
    PeakSite();
    /// Construct from given peak position.
    PeakSite(Vec3 const&);
    /// \return XYZ coords of peak location
    Vec3 const& XYZ() const { return xyz_; }
    /// Add an energy DataSet for each solvent in given array
    int AddEneDataSets(std::vector<SolventInfo> const&, std::string const&, DataSetList&, DataFile*, unsigned int);
    /// Add given frame number to omitted array for all solvents
    void AddOmittedFrame(int fn, unsigned int count) {
      for (SolvPeakArray::iterator it = solvPeaks_.begin(); it != solvPeaks_.end(); ++it) {
        if (count == 0)
          it->AddOmittedFrameNum( fn );
        else if (count > 0)
          it->AddOmittedFrameNum( -fn - 1 );
        if (it->DS() != 0) it->DS()->Add(fn, &ZERO_);
      }
    }
    /// Add an empty SolventPeak for each solvent in given array (when not calculating energy)
    int AddSolventPeaks(std::vector<SolventInfo> const&);
    /// Add given energy for specified solvent; all other solvents get omitted.
    void AddSolventEne(int fn, double ene, unsigned int tgtSidx) {
      for (unsigned int sidx = 0; sidx != solvPeaks_.size(); sidx++)
        if (sidx == tgtSidx)
          solvPeaks_[sidx].DS()->Add(fn, &ene);
        else {
          solvPeaks_[sidx].AddOmittedFrameNum( fn );
          solvPeaks_[sidx].DS()->Add(fn, &ZERO_);
        }
    }

    typedef SolvPeakArray::const_iterator const_iterator;
    /// \return iterator to beginning of solvent peak array
    const_iterator begin() const { return solvPeaks_.begin(); }
    /// \return iterator to end of solvent peak array
    const_iterator end()   const { return solvPeaks_.end(); }
  private:
    static const double ZERO_;

    Vec3 xyz_;                ///< Solvent peak location in Cartesian space.
    SolvPeakArray solvPeaks_; ///< Hold information for each solvent that might occupy this site.
};
#endif
