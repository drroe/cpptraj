#ifndef INC_DATASET_NAMEMAP_H
#define INC_DATASET_NAMEMAP_H
#include "DataSet.h"
#include "NameType.h"
#include <map>
/// Used to map old atom names to new atom names 
class DataSet_NameMap : public DataSet {
    typedef std::pair<NameType,NameType> AmapPair;
    typedef std::map<NameType,NameType> AmapType;
  public:
    DataSet_NameMap();
    static DataSet* Alloc() { return (DataSet*)new DataSet_NameMap(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return nameMap_.size(); }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&)                         { return 1; }
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
    size_t MemUsageInBytes()                         const { return Size() * (2*NameType::max()); }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
    void AddNameMap(NameType const& oldName, NameType const& newName) {
      nameMap_.insert( AmapPair(oldName, newName) );
    }
    bool GetName(NameType& newName, NameType const& oldName) const {
      AmapType::const_iterator it = nameMap_.find( oldName );
      if (it == nameMap_.end()) return false;
      newName = it->second;
      return true;
    }
  private:
    AmapType nameMap_;
};
#endif
