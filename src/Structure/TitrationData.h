#ifndef INC_TITRATIONDATA_H
#define INC_TITRATIONDATA_H
#include <vector>
namespace Cpptraj {
namespace Structure {
// Fwd declares
class TitratableSite;
/// Hold information for all titratable sites
class TitrationData {
  public:
    TitrationData();
  private:
    std::vector<TitratableSite> Sites_; ///< Hold data for potential sites
};
}
}
#endif
