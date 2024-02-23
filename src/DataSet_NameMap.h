#ifndef INC_DATASET_NAMEMAP_H
#define INC_DATASET_NAMEMAP_H
#include "DataSet.h"
#include "NameType.h"
#include <map>
/// Used to map old atom names to new atom names 
class DataSet_NameMap : public DataSet {
  public:
    DataSet_NameMap();
    static DataSet* Alloc() { return (DataSet*)new DataSet_NameMap(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return 0; }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&)                         { return 1; }
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
    size_t MemUsageInBytes()                         const { return 0; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
  private:
    typedef std::pair<NameType,NameType> AmapPair;
    typedef std::map<NameType,NameType> AmapType;
    AmapType nameMap_;
};
#endif
