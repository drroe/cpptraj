#ifndef INC_STRUCTURE_PROTONATOR_H
#define INC_STRUCTURE_PROTONATOR_H
#include <vector>
class ArgList;
class DataSet;
class CpptrajFile;
class CpptrajState;
class DataSet_1D;
class DataSet_2D;
class Random_Number;
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
    /// Calculate titration curves using MC
    int CalcTitrationCurves() const;
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;

    int assign_random_state(Iarray&, Random_Number&) const; // FIXME Random_Number const&

    static double mc_deltae(Iarray const&, int, unsigned int, DataSet_2D const&, DataSet_1D const&, Darray const&);

    int perform_MC_at_pH(double, Iarray const&, DataSet_1D const&,
                         DataSet_2D const&, DataSet_1D const&,
                         Random_Number const&) const;

    DataSet* site_intrinsic_pKas_; ///< DataSet containing calculated intrinsic pKas for each site.
    DataSet* site_site_matrix_;    ///< DataSet containing site-site interactions in e^2/ang
    DataSet* site_qunprot_;        ///< DataSet containing charge of unprotonated state for each site.
    int n_mc_steps_;               ///< Number of monte carlo steps
    double start_pH_;              ///< Starting pH
    double stop_pH_;               ///< final pH
    double pH_increment_;          ///< pH increment
    double min_wint_;              ///< Min inter (pK units) for 2-site transitions
    double fract_toler_;           ///< Tolerance for reduced sites
    double beta_;                  ///< Value to convert from charge to kcal ((1/kT)*C*C)
    int iseed_;                    ///< RNG seed
    enum MCModeType { MC_FULL = 0, MC_REDUCED, MC_CLUSTER };
    MCModeType mcmode_;            ///< MC mode for computing protonation vs pH
    CpptrajFile* logfile_;         ///< Output log file
};
}
}
#endif
