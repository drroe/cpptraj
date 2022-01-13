#ifndef INC_DATASET_POTENTIALFXN_H
#define INC_DATASET_POTENTIALFXN_H
#include "DataSet.h"
class PotentialFunction;
/// Hold 1 or more Potential functions for energy calculations 
class DataSet_PotentialFxn : public DataSet {
  public:
    DataSet_PotentialFxn();
    ~DataSet_PotentialFxn();
    static DataSet* Alloc() { return (DataSet*)new DataSet_PotentialFxn(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return functions_.size(); }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&);
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
    size_t MemUsageInBytes()                         const;
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
    /// \return Pointer to new potential function added to the set
    PotentialFunction* AddNewFunction();
    /// \return Pointer to first function
    PotentialFunction* Pfunction(unsigned int i) const { return functions_[i]; }
  private:
    typedef std::vector<PotentialFunction*> FnArray;
    FnArray functions_;
};
#endif
