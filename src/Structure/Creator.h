#ifndef INC_STRUCTURE_CREATOR_H
#define INC_STRUCTURE_CREATOR_H
class ArgList;
class DataSet_Parameters;
class DataSetList;
namespace Cpptraj {
namespace Structure {
/// Used to create a system from individual units
class Creator {
  public:
    /// CONSTRUCTOR
    Creator();
    /// DESTRUCTOR
    ~Creator();
    /// Associated keywords for InitCreator
    static const char* keywords_;
    /// Initialize the Creator
    int InitCreator(ArgList&, DataSetList const&, int);

    /// \return True if a parameter set is defined
    bool HasMainParmSet() const { return (mainParmSet_ != 0); }
  private:
    /// Get parameter sets
    int getParameterSets(ArgList&, DataSetList const&);

    DataSet_Parameters* mainParmSet_; ///< Hold optional parameter set.
    int debug_;                       ///< Debug level
    bool free_parmset_mem_;           ///< True if main parm set is combined and should be freed
};
}
}
#endif
