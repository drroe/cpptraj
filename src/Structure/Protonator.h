#ifndef INC_STRUCTURE_PROTONATOR_H
#define INC_STRUCTURE_PROTONATOR_H
#include <vector>
#include <string>
#include <utility> // std::pair
class ArgList;
class DataSet;
class CpptrajFile;
class CpptrajState;
class DataSet_1D;
class DataSet_2D;
class DataSet_double;
class Random_Number;
namespace Cpptraj {
namespace Mead {
class MultiFlexResults;
}
namespace Structure {
/// Class for assigning protonation states to a structure
/** This class is largely based off of the MCTI program by Paul Beroza. */
class Protonator {
  public:
    Protonator();
    /// Print help to stdout
    static void HelpText();
    /// Set up from arguments
    int SetupProtonator(CpptrajState&, ArgList&, Cpptraj::Mead::MultiFlexResults const&);
    /// Print options to stdout
    void PrintOptions() const;
    /// Calculate titration curves using MC
    int CalcTitrationCurves() const;
  private:
    class StateArray;
    class MC_Corr;
    typedef std::vector<double> Darray;
    typedef std::vector<int> Iarray;
    typedef std::pair<int,int> StatePair;
    typedef std::vector<StatePair> PairArray;

    int read_files(CpptrajState&, std::string const&);

    static double mc_deltae(StateArray const&, int, unsigned int, DataSet_2D const&, DataSet_1D const&, Darray const&);

    void mc_pair_flip(double&, unsigned int, Darray const&, StatePair const&, StateArray&,
                      DataSet_2D const&, DataSet_1D const&, Random_Number const&) const;

    void mc_step(double&, unsigned int, Darray const&, PairArray const&, StateArray&,
                 DataSet_2D const&, DataSet_1D const&, Random_Number const&) const;

    int perform_MC_at_pH(double, StateArray&, MC_Corr&, int, DataSet_1D const&,
                         DataSet_2D const&, DataSet_1D const&,
                         Random_Number const&, PairArray const&) const;

    int reduce_sites(DataSet_double&, DataSet_1D&, DataSet_2D&, Iarray&,
                     DataSet_1D const&, DataSet_1D const&, DataSet_2D const&, MC_Corr const&) const;

    PairArray get_pairs(double, unsigned int, DataSet_2D const&) const;

    DataSet* site_intrinsic_pKas_; ///< Input DataSet containing calculated intrinsic pKas for each site.
    DataSet* site_site_matrix_;    ///< Input DataSet containing site-site interactions in e^2/ang
    DataSet* site_qunprot_;        ///< Input DataSet containing charge of unprotonated state for each site.
    DataSet* site_names_;          ///< Input DataSet containing name of each site.
    int n_mc_steps_;               ///< Number of monte carlo steps
    int n_reduced_mc_steps_;       ///< Number of monte carlo steps for reduced sites
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
    CpptrajFile* pkoutfile_;       ///< Pkout formatted file
    DataSet* pkhalf_;              ///< pkhalf for each site
    std::vector<DataSet*> TitrationCurves_; ///< Hold titration curve for each site
};
}
}
#endif
