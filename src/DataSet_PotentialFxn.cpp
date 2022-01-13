#include "DataSet_PotentialFxn.h"
#include "CpptrajStdio.h"
#include "Potential/PotentialFunction.h"

/** CONSTRUCTOR */
DataSet_PotentialFxn::DataSet_PotentialFxn() :
  // 0 dim indicates DataSet-specific write
  DataSet(POTENTIALFXN, GENERIC, TextFormat(), 0)
{ }

/** DESTRUCTOR */
DataSet_PotentialFxn::~DataSet_PotentialFxn() {
  for (FnArray::iterator it = functions_.begin(); it != functions_.end(); ++it)
    delete *it;
}

/** Reserve space for a certain number of potential functions. */
int DataSet_PotentialFxn::Allocate(SizeArray const& sizes) {
  if (!sizes.empty())
    functions_.reserve( sizes[0] );
  return 0;
}

/** \return size in bytes used by potential function. */
size_t DataSet_PotentialFxn::MemUsageInBytes() const {
  // TODO each potential term should return usage in bytes
  return (functions_.size() * sizeof(PotentialFunction));
}

/** \return pointer to new potential function added to this set */
PotentialFunction* DataSet_PotentialFxn::AddNewFunction() {
  functions_.push_back( new PotentialFunction() );
  return (functions_.back());
}
