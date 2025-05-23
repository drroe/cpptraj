#ifndef INC_MINIMIZE_STEEPESTDESCENT_H
#define INC_MINIMIZE_STEEPESTDESCENT_H
#include <string>
class PotentialFunction;
class Frame;
class CpptrajFile;
class DataSet;
class DataSet_1D;
class Minimize_SteepestDescent {
  public:
    Minimize_SteepestDescent();
    /** Set up minimization with optional output trajectory, RMS tolerance, step size, # steps, and optional energy data set. */
    int SetupMin(std::string const&, double, double,int,DataSet*);
    /** Set up minimization with optional output trajectory, RMS tolerance, step size, and # steps. */
    int SetupMin(std::string const&, double, double,int);
    /** Run minimization with given potential function and coordinates. */
    int RunMin(PotentialFunction&, Frame&, CpptrajFile&) const;
  private:
    std::string trajoutName_;
    double min_tol_; ///< Min RMS tolerance
    double dx0_;     ///< Initial step size
    int nMinSteps_;  ///< Number of minimization steps.
    DataSet_1D* eneSet_; ///< Optional energy output set
};
#endif
