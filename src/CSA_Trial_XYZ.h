#ifndef INC_CSA_TRIAL_XYZ_H
#define INC_CSA_TRIAL_XYZ_H
#include "CSA_Trial.h"
#include <vector>
namespace Cpptraj {
namespace CSA {

/** Trial stored by Cartesian coordinates */
class Trial_XYZ : public Trial {
  public:
    Trial_XYZ {}
  private:
    typedef std::vector<double> Darray;

    Darray xyz_;
    unsigned int natom_;
};

}
}
#endif
