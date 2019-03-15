#ifndef FILE_TEMPNAME_H
#define FILE_TEMPNAME_H
namespace File {
  // Forward declaration
  class Name;
  /// \return A temporary file name.
  Name GenTempName();
  /// Indicate temporary file name no longer needed.
  void FreeTempName(Name const&);
} /* END namespace File */
#endif
