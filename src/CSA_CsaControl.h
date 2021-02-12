#ifndef INC_CSA_CSACONTROL_H
#define INC_CSA_CSACONTROL_H
#include "CSA.h"
#include "CSA_Bank.h"
// Forward declares
class ArgList;
namespace Cpptraj {
namespace CSA {

/** Used to control a conformational space annealing run. */
class CsaControl {
  public:
    CsaControl();

    static void Help();

    int InitCsa(ArgList&);
  private:
    Bank currentBank_;
    Bank firstBank_;

    // User-defined variables
    int num_conf_;     ///< Number of conformations to store in bank/firstbank
    int num_iter_;     ///< Number of CSA iterations to do.
    int num_seed_;     ///< Number of seeds to use in a round.
    int rnd_per_iter_; ///< Rounds per iteration; if not set, use num_conf/num_seed, rounded up
    int nmove1_;       ///< Number of moves of type 1 for each seed
    int nmove2_;       ///< Number of moves of type 2 for each seed
    int nmove3_;       ///< Number of moves of type 3 for each seed
    int nmove4_;       ///< Number of moves of type 4 for each seed
    int move_maxnres_; ///< For moves of multiple residues, max number of residues to use
    int move_res0_;    ///< For moves of residue blocks, the min # residues to use.
    int move_res1_;    ///< For moves of residue blocks, the max # residues to use.
    DistType idist_;   ///< Type of distance calculation (parameters IDIST_)
    ScoreType ietyp_;  ///< Type of energy evaluation (parameters IETYP_)
    ManipType imanip_; ///< How to manipulate coordinates (parameters IMANIP_)
    DcutType idcut_;   ///< How to decrease Dcut after each bank update
    int iranseed_;     ///< Seed for random number generator
    int nolower_cut_;  ///< Exit when lower conf is not found this many times in a row.

    bool use_seeds_;   ///< If true, use seeds in trials.
    bool free_ene_;    ///< If true, energy will be modified by -T*S
    bool write_banks_; ///< If true, first bank and banks during iterations will be output to files.

    double min_ecut_;    ///< End calculation when energy of minimum conformation is below this
    double dcut0_fac_;   ///< Factor for determining initial dcut value.
    double dcut1_fac_;   ///< Factor for determining final dcut value.
    double temperature_; ///< Temperature : used when setting trial energy. Keyword 'TEMP'
};

}
}
#endif
