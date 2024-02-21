#ifndef INC_STRUCTURE_PDBCLEANER_H
#define INC_STRUCTURE_PDBCLEANER_H
#include <vector>
#include <string>
class ArgList;
namespace Cpptraj {
namespace Structure {
/// Used to clean up structures read from a PDB
class PdbCleaner {
    typedef std::vector<int> Iarray;
  public:
    PdbCleaner();
    /// Initialize
    int InitPdbCleaner(ArgList&, std::string const&, Iarray const&);
  private:
    bool remove_water_;      ///< If true, remove any water.
    bool remove_h_;          ///< If true, remove hydrogen atoms.
    std::string waterMask_;  ///< Mask expression for selecting water.
    std::string altLocArg_;  ///< Only keep atoms with this alternate location identifier.
    std::string stripMask_;  ///< General mask for removing atoms.
    Iarray resnumsToRemove_; ///< Other residue numbers to remove.
};
}
}
#endif
