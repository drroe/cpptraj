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
    /// \return Number corresponding to given index
    int NumberAtIdx(int idx) const { return (idx * offset_) + start_; }
    /// \return first number
    int FirstNumber() const { return start_; }
    /// \return current number
    int CurrentNumber() const { return current_; }
    /// \return previous number
    int PreviousNumber() const { return current_ - offset_; }
    /// \return True if the count is finished
    bool IsFinished() const { return !(current_ < stop_ || stop_ == -1); }
    /// \return string with counter info (<start>-<stop>, <offset>)
    std::string CounterInfo() const;
    /// \return Total number to be counted
    int CounterTotal() const { return total_read_frames_; }
  private:
    void update() { current_ += offset_; }
    void start()  { current_ = start_; }

    /// Determine total from arguments
    int determineTotal() const;

    int start_;   ///< Start number
    int stop_;    ///< Stop number
    int offset_;  ///< Increment
    int current_; ///< Current number
    int total_read_frames_; ///< Total number to be accessed based on args
};

}
#endif
