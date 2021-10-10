#ifndef INC_COUNTER_H
#define INC_COUNTER_H
namespace Cpptraj {
/// Used to generate and keep track of a series of numbers to be iterated over.
class Counter {
  public:
    /// CONSTRUCTOR
    Counter();
    /// \return the current number.
    virtual int CurrentNumber() const = 0;
    /// Go to the next number.
    virtual void UpdateCounter() const = 0;
  private:
    int currentIdx_; ///< Current internal index.
};

}
#endif
