#ifndef INC_COUNTER_ARRAY_H
#define INC_COUNTER_ARRAY_H
#include "Counter.h"
#include <vector>
namespace Cpptraj {
/// Counter for numbers that may not be in a monotonic order
class Counter_Array : public Counter {
    typedef std::vector<int> Iarray;
  public:
    /** CONSTRUCTOR */
    Counter_Array();
    /** CONSTRUCTOR - take array of numbers. */
    Counter_Array(Iarray const&);
    /// \return Number at given index
    int NumberAtIdx(int idx) const { return numbers_[idx]; }
    /// \return current number
    int CurrentNumber() const { return numbers_[CurrentIdx()]; }
    /// \return previous number
    int PreviousNumber() const { return numbers_[CurrentIdx()-1]; }
    /// \return true if count is finished
    bool IsFinished() const { return CurrentIdx() >= numbers_.size(); }
    /// \return string containing counter info
    std::string CounterInfo() const;
    /// \return Total to be counted
    int CounterTotal() const { return (int)numbers_.size(); }
  private:
    /// update internal counter. Nothing needed since using CurrentIdx()
    void update() { return; }
    /// position at first number. Nothing needed since using CurrentIdx()
    void start() { return; }
    /// position at CurrentIdx();
    void assign() { return; }

    Iarray numbers_;
};
}
#endif
