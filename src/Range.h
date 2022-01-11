#ifndef INC_RANGE_H
#define INC_RANGE_H
#include <list>
#include <string>
// Class: Range
/// The Range class is used to hold an ordered list of numbers. 
/** Range can either be set using a range expression, e.g. 
  * X-Y,C-D (X to Y and C to D) etc, or can be set using a beginning 
  * and end number.
*/
class Range {
  public:
    /// CONSTRUCTOR
    Range();
    /// CONSTRUCTOR - Takes argument string as input
    Range(std::string const&);
    /// CONSTRUCTOR - Single number
    Range(int);
    /// CONSTRUCTOR - Argument with offset
    Range(std::string const&,int);
    /// COPY CONSTRUCTOR
    Range(const Range&);
    /// ASSIGNMENT
    Range& operator=(const Range&);

    enum SortType { SORTED = 0, UNSORTED };

    typedef std::list<int>::const_iterator const_iterator;
    const_iterator begin() const { return rangeList_.begin();      }
    const_iterator end()   const { return rangeList_.end();        }
    bool Empty()           const { return rangeList_.empty();      }
    int Size()             const { return (int) rangeList_.size(); }
    int Back()             const { return rangeList_.back();       }
    int Front()            const { return rangeList_.front();      }
    void Clear() { rangeArg_.clear(); rangeList_.clear(); }

    /// Set range from argument
    int SetRange(std::string const&, SortType);
    /// Set range from argument, sorted
    int SetRange(std::string const& a) { return SetRange(a, SORTED); }
    /// Set range from start up to but not including end
    int SetRange(int,int);
    /// \return the range argument
    const char *RangeArg() const { return rangeArg_.c_str(); }
    void PrintRange(const char*,int) const;

    void ShiftBy(int);
    /// Add a number to the range. Range is NOT explicitly sorted in this case. 
    void AddToRange(int num) { rangeList_.push_back(num); }
    /// Remove number from range
    void RemoveFromRange(int);
    /// \return true if given number is in the Range list.
    bool InRange(int) const;
  private:
    std::string rangeArg_;
    std::list<int> rangeList_;
};
#endif
