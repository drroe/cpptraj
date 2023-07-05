#ifndef INC_RANGE_H
#define INC_RANGE_H
#include <vector>
#include <string>
// Class: Range
/// The Range class is used to hold an ordered list of numbers. 
/** Range can either be set using a range expression, e.g. 
  * X-Y,C-D (X to Y and C to D) etc, or can be set using a beginning 
  * and end number.
*/
class Range {
    typedef std::vector<int> RangeType;
  public:
    /// CONSTRUCTOR - empty range
    Range();
    /// CONSTRUCTOR - range expression
    Range(std::string const&);
    /// CONSTRUCTOR - single number
    Range(int);
    /// CONSTRUCTOR - range expression and offset
    Range(std::string const&,int);
    /// COPY CONSTRUCTOR
    Range(const Range&);
    /// ASSIGNMENT
    Range& operator=(const Range&);

    enum SortType { SORTED = 0, UNSORTED };

    typedef RangeType::const_iterator const_iterator;
    /// \return const iterator to beginning of range
    const_iterator begin() const { return rangeList_.begin();      }
    /// \return const iterator to end of range
    const_iterator end()   const { return rangeList_.end();        }
    /// \return true if range is empty
    bool Empty()           const { return rangeList_.empty();      }
    /// \return Total numbers in range
    unsigned int Size()    const { return rangeList_.size();       }
    /// \return Last number in range
    int Back()             const { return rangeList_.back();       }
    /// \return First number in range
    int Front()            const { return rangeList_.front();      }
    /// Clear the range
    void Clear() { rangeArg_.clear(); rangeList_.clear(); }

    /// Set range from argument, optionally sorted
    int SetRange(std::string const&, SortType);
    /// Set range from argument, sorted
    int SetRange(std::string const& a) { return SetRange(a, SORTED); }
    /// Set range from start up to but not including end
    int SetRange(int,int);

    /// Remove number from range
    void RemoveFromRange(int);
    /// \return true if given number is in the Range list.
    bool InRange(int) const;
    /// Add a number to the range. Range will still be sorted. FIXME save sort stat
    void AddToRange(int);
    /// Shift all numbers in the range by a constant
    void ShiftBy(int);

    /// \return the range expression used to set up range
    const char *RangeArg() const { return rangeArg_.c_str(); }
    /// \return string containing actual range, with optional offset.
    std::string PrintRange(int) const;
    /// Print entire range to stdout
    void PrintToStdout() const;
  private:
    /// Add numbers from start up to but not including end to range.
    int setRange(int, int);

    std::string rangeArg_; ///< Expression that describes the range
    RangeType rangeList_;  ///< Array containing numbers in the range
};
#endif
