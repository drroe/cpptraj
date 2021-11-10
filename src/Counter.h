#ifndef INC_COUNTER_H
#define INC_COUNTER_H
#include <string>
//namespace Cpptraj {
/// Used to generate and keep track of a series of numbers to be iterated over.
class Counter {
  public:
    /// CONSTRUCTOR
    Counter();
    /// DESTRUCTOR - virtual since inherited
    virtual ~Counter() {}
    /// \return number at the specified index
    virtual int NumberAtIdx(int) const = 0;
    /// \return the first number.
    virtual int FirstNumber() const = 0;
    /// \return the current number.
    virtual int CurrentNumber() const = 0;
    /// \return previous frame number (last # before UpdateCounter() was called).
    virtual int PreviousNumber() const = 0;
    /// \return True if the count is finished
    virtual bool IsFinished() const = 0;
    /// \return String containing counter info
    virtual std::string CounterInfo() const = 0;
    /// \return String with more verbose counter info
    virtual std::string Verbose(int) const = 0;
    /// \return Total number of frames represented by the counter
    virtual int CounterTotal() const = 0;

    /// Start the counter
    void StartCounter() { currentIdx_ = 0; start(); }
    /// Go to the next number and update internal index
    void UpdateCounter() { ++currentIdx_; update(); }
    /// \return the current internal index
    unsigned int CurrentIdx() const { return currentIdx_; }
  protected:
    /// Go to the next number.
    virtual void update() = 0;
    /// Position at first number.
    virtual void start() = 0;
  private:
    unsigned int currentIdx_; ///< Current internal index. Also serves to count calls to UpdateCounter()
};

//}
#endif
