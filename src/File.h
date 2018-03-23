#ifndef INC_FILE_H
#define INC_FILE_H
#include <string>
#include <vector>
#include "FileName.h"
/// This namespace contains useful file-related routines.
namespace File {
  /// Array of file names
  typedef std::vector<Name> NameArray;
  /// Expand given expression to array of file names
  NameArray ExpandToFilenames(std::string const&);
  /// File access types
  enum AccessType { READ=0, WRITE, APPEND, UPDATE };
  /// File compression types
  enum CompressType { NO_COMPRESSION=0, GZIP, BZIP2, ZIP };
  /// Print error message with given file name corresponding to 'false' value from 'Exists()'
  void ErrorMsg(const char*);
  /// return last file error message
  const char* StrError();
  /// \return true if file exists and is accessible.
  bool Exists(Name const&);
  //bool Exists(std::string const&); // TODO remove?
  /// Given lowest replica name search for other replica names
  NameArray SearchForReplicas(Name const&, int);
# ifdef MPI
  /// Each rank searches for replica based on lowest replica number.
  NameArray SearchForReplicas(Name const&, bool, int, int, int);
# endif
} /* END namespace File */
#endif
