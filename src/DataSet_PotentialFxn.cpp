#include "DataSet_PotentialFxn.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
DataSet_PotentialFxn::DataSet_PotentialFxn() :
  // 0 dim indicates DataSet-specific write
  DataSet(POTENTIALFXN, GENERIC, TextFormat(), 0)
{

}
