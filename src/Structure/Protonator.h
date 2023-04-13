#ifndef INC_STRUCTURE_PROTONATOR_H
#define INC_STRUCTURE_PROTONATOR_H
class ArgList;
class DataSet;
class CpptrajFile;
class CpptrajState;
namespace Cpptraj {
namespace Mead {
class MultiFlexResults;
}
namespace Structure {
/// Class for assigning protonation states to a structure
class Protonator {
  public:
    Protonator();
    /// Set up from arguments
    int SetupProtonator(CpptrajState&, ArgList&, Cpptraj::Mead::MultiFlexResults const&);
    /// Print options to stdout
    void PrintOptions() const;
  private:
    DataSet* site_intrinsic_pKas_; ///< DataSet containing calculated intrinsic pKas for each site.
    DataSet* site_site_matrix_;    ///< DataSet containing site-site interactions in e^2/ang
    int n_mc_steps_;               ///< Number of monte carlo steps
    double start_pH_;              ///< Starting pH
    double stop_pH_;               ///< final pH
    double pH_increment_;          ///< pH increment
    double min_wint_;              ///< Min inter (pK units) for 2-site transitions
    double fract_toler_;           ///< Tolerance for reduced sites
    int iseed_;                    ///< RNG seed
    enum MCModeType { MC_FULL = 0, MC_REDUCED, MC_CLUSTER };
    MCModeType mcmode_;            ///< MC mode for computing protonation vs pH
    CpptrajFile* logfile_;         ///< Output log file
};
}
}
#endif
