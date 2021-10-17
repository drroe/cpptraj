#ifndef INC_COUNTER_REGULAR_H
#define INC_COUNTER_REGULAR_H
#include "Counter.h"
namespace Cpptraj {
/// Regular count with a start, stop, and offset
class Counter_Regular : public Counter {
  public:
    /// CONSTRUCTOR
    Counter_Regular();
    /// CONSTRUCTOR - start/stop/offset
    Counter_Regular(int,int,int);
    /// \return current number
    int CurrentNumber() const { return current_; }
    /// \return True if the count is finished
    bool IsFinished() const { return !(current_ < stop_ || stop_ == -1); }
  private:
    void update() { current_ += offset_; }

    int start_;   ///< Start number
    int stop_;    ///< Stop number
    int offset_;  ///< Increment
    int current_; ///< Current number
};

}
#endif
