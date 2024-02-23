#ifndef INC_STRUCTURE_CREATOR_H
#define INC_STRUCTURE_CREATOR_H
#include <vector>
#include <string>
#include "StructureEnum.h" // TerminalType
class ArgList;
class DataSet_Coords;
class DataSet_NameMap;
class DataSet_Parameters;
class DataSetList;
class NameType;
namespace Cpptraj {
namespace Structure {
/// Used to create a system from individual units
class Creator {
    typedef std::vector<DataSet_Coords*> Carray;
    typedef std::vector<DataSet_NameMap*> Narray;
  public:
    /// CONSTRUCTOR
    Creator();
    /// DESTRUCTOR
    ~Creator();
    /// Associated parameter keywords for InitCreator
    static const char* parm_keywords_;
    /// Associated template keywords for InitCreator
    static const char* template_keywords_;
    /// Other keywords for InitCreator
    static const char* other_keywords_;
    /// Initialize the Creator
    int InitCreator(ArgList&, DataSetList const&, int);

    /// \return True if a parameter set is defined
    bool HasMainParmSet() const { return (mainParmSet_ != 0); }
    /// \return Main parm set
    DataSet_Parameters const* MainParmSetPtr() const { return mainParmSet_; }
    /// \return True if there are templates
    bool HasTemplates() const { return (!Templates_.empty()); }
    /// Identify residue template from name
    DataSet_Coords* IdTemplateFromName(std::string const&) const;
    /// Identify residue template from residue name
    DataSet_Coords* IdTemplateFromResname(NameType const&, TerminalType) const;
    /// Get name map if its present
    bool GetAlias(NameType&, NameType const&) const;
  private:
    /// Get templates
    int getTemplates(ArgList&, DataSetList const&);
    /// Get parameter sets
    int getParameterSets(ArgList&, DataSetList const&);

    DataSet_Parameters* mainParmSet_; ///< Hold optional parameter set.
    Carray Templates_;                ///< Hold unit templates.
    Narray NameMaps_;                 ///< Hold atom name maps.
    int debug_;                       ///< Debug level
    bool free_parmset_mem_;           ///< True if main parm set is combined and should be freed
};
}
}
#endif
