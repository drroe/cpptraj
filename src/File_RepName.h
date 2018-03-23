#ifndef INC_FILE_REPNAME_H
#define INC_FILE_REPNAME_H
#include "File.h"
namespace File {

class RepName {
  public:
    RepName() : ExtWidth_(0), lowestRepnum_(-1), extChar_('.') {}
    RepName(Name const&, int);
    bool Error() const { return Prefix_.empty(); }
    /// \return Replica file name for given offset from lowest replica number.
    Name RepFilename(int) const;
  private:
    std::string Prefix_;      ///< File name up to the numerical extension.
    std::string ReplicaExt_;  ///< Numerical extension.
    std::string CompressExt_; ///< Optional compression extension after numerical extension.
    int ExtWidth_;            ///< Width of the numerical extension. TODO remove
    int lowestRepnum_;        ///< Integer value of numerical extension.
    char extChar_;            ///< Character preceding numerical extension
};

} /* END namespace File */
#endif
