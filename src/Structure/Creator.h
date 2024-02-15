#ifndef INC_STRUCTURE_CREATOR_H
#define INC_STRUCTURE_CREATOR_H
#include <vector>
class ArgList;
class DataSet_Coords;
class DataSet_Parameters;
class DataSetList;
namespace Cpptraj {
namespace Structure {
/// Used to create a system from individual units
class Creator {
    typedef std::vector<DataSet_Coords*> Carray;
  public:
    /// CONSTRUCTOR
    Creator();
    /// DESTRUCTOR
    ~Creator();
    /// Associated parameter keywords for InitCreator
    static const char* parm_keywords_;
    /// Associated templtae keywords for InitCreator
    static const char* template_keywords_;
    /// Initialize the Creator
    int InitCreator(ArgList&, DataSetList const&, int);

    /// \return True if a parameter set is defined
    bool HasMainParmSet() const { return (mainParmSet_ != 0); }
    /// \return Main parm set
    DataSet_Parameters const& MainParmSet() const { return *mainParmSet_; }
    /// \return True if there are templates
    bool HasTemplates() const { return (!Templates_.empty()); }
  private:
    /// Get templates
    int getTemplates(ArgList&, DataSetList const&);
    /// Get parameter sets
    int getParameterSets(ArgList&, DataSetList const&);

    DataSet_Parameters* mainParmSet_; ///< Hold optional parameter set.
    Carray Templates_;                ///< Hold unit templates.
    int debug_;                       ///< Debug level
    bool free_parmset_mem_;           ///< True if main parm set is combined and should be freed
};
}
}
#endif
