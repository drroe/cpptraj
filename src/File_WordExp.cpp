#ifndef _WIN32
#   include <wordexp.h>
#endif
#include "File_WordExp.h"
#include "CpptrajStdio.h"

using namespace File;

#ifndef _WIN32
/** Print error messages from the wordexp() function. */
static void WexpErr(int err) {
  switch ( err ) {
    case WRDE_BADCHAR :
      mprinterr("Error: Illegal occurrence of newline or one of |, &, ;, <, >, (, ), {, }.\n");
      break;
    //case WRDE_BADVAL
    case WRDE_CMDSUB :
      mprinterr("Error: Command substitution is not allowed in file names.\n");
      break;
    case WRDE_NOSPACE :
      mprinterr("Error: Out of memory.\n");
      break;
    case WRDE_SYNTAX :
      mprinterr("Error: Bad syntax (unbalanced parentheses, unmatched quotes.\n");
      break;
  }
}
#endif /* _WIN32 */

/** Given a name, expand it to a name or names with wordexp. */
Sarray WordExp(std::string const& nameIn) {
  if (nameIn.empty()) return Sarray();
# ifdef _WIN32
  // Expansion NOT supported.
  return Sarray(1, nameIn);
# else
  Sarray namesOut;
  wordexp_t expanded;
  int err = wordexp( nameIn.c_str(), &expanded, WRDE_NOCMD );
  WexpErr( err );
  if (err == 0) {
    for (unsigned int i = 0; i != expanded.we_wordc; i++) {
      if (expanded.we_wordv[i] == 0)
        mprinterr("Internal Error: Bad expansion at %i\n", i);
      else
        namesOut.push_back( expanded.we_wordv[i] );
    }
    wordfree( &expanded );
  }
  return namesOut;
# endif /* _WIN32 */
}
