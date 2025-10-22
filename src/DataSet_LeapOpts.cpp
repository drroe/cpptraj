#include "DataSet_LeapOpts.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_LeapOpts::DataSet_LeapOpts() :
  // 0 dim indicates DataSet-specific write
  DataSet(LEAPOPTS, GENERIC, TextFormat(TextFormat::STRING, 12, 0), 0)
{

}
