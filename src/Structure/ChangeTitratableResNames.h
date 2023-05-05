#ifndef INC_STRUCTURE_CHANGETITRATABLERESNAMES_H
#define INC_STRUCTURE_CHANGETITRATABLERESNAMES_H
#include <string>
class Topology;
namespace Cpptraj {
namespace Structure {
/// Change recognized titratable residue names to ones that correspond to sites
int ChangeTitratableResNames(Topology&, std::string const&);
}
}
#endif
