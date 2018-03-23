#ifndef INC_FILE_WORDEXP_H
#define INC_FILE_WORDEXP_H
#include <vector>
#include <string>
namespace File {

typedef std::vector<std::string> Sarray;

Sarray WordExp(std::string const&);

} /* END namespace File */
#endif
