#ifndef INC_COUNTER_H
#define INC_COUNTER_H
namespace Cpptraj {
/// Used to generate and keep track of a series of numbers to be iterated over.
class Counter {
  public:
    /// CONSTRUCTOR
    Counter();
    /// DESTRUCTOR - virtual since inherited
    virtual ~Counter() {}
    /// \return the current number.
    virtual int CurrentNumber() const = 0;
    /// \return True if the count is finished
    virtual bool IsFinished() const = 0;

    /// Go to the next number and update internal index
    void UpdateCounter() { ++currentIdx_; update(); }
  protected:
    /// Go to the next number.
    virtual void update() = 0;
  private:
    int currentIdx_; ///< Current internal index.
};

}
#endif
