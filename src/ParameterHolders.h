#ifndef INC_PARAMETERHOLDERS_H
#define INC_PARAMETERHOLDERS_H
#include <vector>
#include <utility> // std::pair
#include "TypeNameHolder.h"
#include "ParameterTypes.h"
//#inc lude "CpptrajStdio.h" // DEBUG

namespace ParameterHolders {
  enum RetType { ADDED = 0, SAME, UPDATED, ERR };
} /* END namespace ParameterHolders */

/// Used to associate atom type names with an object (parameter etc)
template <class T> class ParmHolder {
    // TODO may want to actually use a map one day for performance reasons.
    typedef std::pair<TypeNameHolder,T> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    ParmHolder() {}
    void clear()              { bpmap_.clear(); }
    size_t size()       const { return bpmap_.size(); }
    bool empty()        const { return bpmap_.empty(); }
    /// Set wildcard character
    void SetWildcard(char wc) { wc_ = NameType(std::string(1, wc)); }
    /// Add (or update if allowed) given parameter to holder.
    ParameterHolders::RetType AddParm(TypeNameHolder const& types, T const& bp, bool allowUpdate) {
      // Check if parm for these types exist
      typename Bmap::iterator it = bpmap_.begin();
      for (; it != bpmap_.end(); ++it)
        if (it->first.Match_NoWC( types )) break;
      if (it == bpmap_.end()) {
        // New parm
        //mprintf("DEBUG: New parameter:");
        //for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
        //  mprintf(" '%s'", *(*it));
        //mprintf("\n");
        bpmap_.push_back( Bpair(types, bp) );
      } else {
        if (bp < it->second || it->second < bp) {
          //mprintf("DEBUG: Potential update of existing parameter:");
          //for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
          //  mprintf(" '%s'", *(*it));
          //mprintf("\n");
          if (allowUpdate) {
            it->second = bp;
            return ParameterHolders::UPDATED;
          } else {
            return ParameterHolders::ERR;
          }
        } else {
          //mprintf("DEBUG: Existing parameter:");
          //for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
          //  mprintf(" '%s'", *(*it));
          //mprintf("\n");
          return ParameterHolders::SAME;
        }
      }
      return ParameterHolders::ADDED;
    }
    /// Constant iterator
    typedef typename Bmap::const_iterator const_iterator;
    /// \return constant iterator to beginning
    const_iterator begin() const { return bpmap_.begin(); }
    /// \return constant iterator to end.
    const_iterator end()   const { return bpmap_.end();   }
    /// Iterator
    typedef typename Bmap::iterator iterator;
    /// \return iterator to beginning
    iterator begin() { return bpmap_.begin(); }
    /// \return iterator to end
    iterator end()   { return bpmap_.end();   }
    /// \return Parameter matching given types, or empty parameter if not found.
    T FindParam(TypeNameHolder const& types, bool& found) const { // TODO only use GetParam()?
      found = true;
      for (const_iterator it = begin(); it != end(); ++it)
        if (it->first.Match_NoWC( types )) return it->second;
      if (wc_.len() > 0) {
        for (const_iterator it = begin(); it != end(); ++it)
          if (it->first.Match_WC( types, wc_)) return it->second;
      }
      found = false;
      return T();
    }
    /// \return iterator to parameter matching the given types.
    iterator GetParam(TypeNameHolder const& types) {
      for (iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first.Match_NoWC( types )) return it;
      if (wc_.len() > 0) {
        for (iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
          if (it->first.Match_WC( types, wc_)) return it;
      }
      return bpmap_.end();
    }
    /// \return const iterator to parameter matching the given types.
    const_iterator GetParam(TypeNameHolder const& types) const {
      for (const_iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first.Match_NoWC( types )) return it;
      if (wc_.len() > 0) {
        for (const_iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
          if (it->first.Match_WC( types, wc_)) return it;
      }
      return bpmap_.end();
    }
    /// \return size in memory in bytes
    size_t DataSize() const {
      if (bpmap_.empty()) return 0;
      const_iterator elt0 = begin();
      // Assume all TypeNameHolders are the same size
      return (bpmap_.size() * elt0->first.DataSize()) +
             (bpmap_.size() * sizeof(T)) +
             sizeof(Bmap);
    }
  private:
    Bmap bpmap_;
    NameType wc_; ///< Wildcard character
};

// -----------------------------------------------------------------------------
/// Specialized class for associating atom types with dihedral parameters.
/** NOTE: Instead of using a specialize template here I'm creating a new
  *       class because while I want AddParm() to accept DihedralParmType,
  *       I want FindParam to return an array of DihedralParmType, one for
  *       each unique multiplicity.
  */
class DihedralParmHolder {
    typedef std::pair<TypeNameHolder,DihedralParmArray> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    DihedralParmHolder() {}
    virtual ~DihedralParmHolder() {} ///< Virtual since inherited
    void clear()              { bpmap_.clear();        }
    size_t size()       const { return bpmap_.size();  }
    bool empty()        const { return bpmap_.empty(); }
    /// Set wildcard character
    void SetWildcard(char wc) { wc_ = NameType(std::string(1, wc)); }
    /** Add (or update) a single dihedral parameter for given atom types. */
    ParameterHolders::RetType
    AddParm(TypeNameHolder const& types, DihedralParmType const& dp, bool allowUpdate) {
      // Check if parm for these types exist
      Bmap::iterator it0 = bpmap_.begin();
      for (; it0 != bpmap_.end(); ++it0)
      {
        if (it0->first.Match_NoWC( types ))
          break;
      }
      if (it0 == bpmap_.end()) {
        // Brand new dihedral for these types.
        //mprintf("DEBUG: New dihedral parm: %s %s %s %s pk=%12.4f pn=%12.4f pp=%12.4f\n",
        //        *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase());
        bpmap_.push_back( Bpair(types, DihedralParmArray(1, dp)) );
      } else {
        // If we are here types match - check multiplicity.
        DihedralParmArray::iterator it1 = it0->second.begin();
        for (; it1 != it0->second.end(); ++it1)
        {
          if (FEQ(it1->Pn(), dp.Pn()))
            break;
        }
        if (it1 == it0->second.end()) {
          // Brand new multiplicity for this dihedral.
          //mprintf("DEBUG: Dihedral new mult: %s %s %s %s pk=%12.4f pn=%12.4f pp=%12.4f\n",
          //        *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase());
          if (it0->second.empty())
            it0->second.push_back( dp );
          else if (dp.Pn() > it0->second.back().Pn())
            it0->second.push_back( dp );
          else {
            // Try to keep multiplicities in order.
            DihedralParmArray sorted;
            bool isInserted = false;
            for (DihedralParmArray::const_iterator jt = it0->second.begin(); jt != it0->second.end(); ++jt) {
              if (!isInserted) {
                if (dp.Pn() < jt->Pn()) {
                  sorted.push_back( dp );
                  isInserted = true;
                }
              }
              sorted.push_back( *jt );
            }
            it0->second = sorted;
          }
        } else {
          if (dp < *it1 || *it1 < dp) {
            //mprintf("DEBUG: Attempt dihedral update mult (allow=%i): %s %s %s %s pk=%6.2f pn=%3.1f pp=%6.3f (orig pk=%6.2f pn=%3.1f pp=%6.3f )\n",
            //        (int)allowUpdate, *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase(), it1->Pk(), it1->Pn(), it1->Phase());
            if (allowUpdate) {
              *it1 = dp;
              return ParameterHolders::UPDATED;
            } else {
              return ParameterHolders::ERR;
            }
          } else
            return ParameterHolders::SAME;
        }
      }
      return ParameterHolders::ADDED;
    }

    /** This version takes an array of dihedral parameters. */
    ParameterHolders::RetType
    AddParm(TypeNameHolder const& types, DihedralParmArray const& dpa, bool allowUpdate) {
      // Check if parm for these types exist
      Bmap::iterator it0 = bpmap_.begin();
      for (; it0 != bpmap_.end(); ++it0)
      {
        if (it0->first.Match_NoWC( types ))
          break;
      }
      if (it0 == bpmap_.end()) {
        // Brand new dihedral for these types.
        bpmap_.push_back( Bpair(types, dpa) );
      } else {
        if (!allowUpdate) return ParameterHolders::ERR;
        // Check if sizes are the same.
        bool update = false;
        if (it0->second.size() != dpa.size())
          update = true;
        else {
          // Sizes are the same. See if parameters are the same.
          for (unsigned int i = 0; i != it0->second.size(); i++) {
            if (it0->second[i] < dpa[i] || dpa[i] < it0->second[i]) {
              update = true;
              break;
            }
          }
        }
        if (update) {
          it0->second = dpa;
          return ParameterHolders::UPDATED;
        } else
          return ParameterHolders::SAME;
      }
      return ParameterHolders::ADDED;
    }

    typedef typename Bmap::const_iterator const_iterator;
    const_iterator begin() const { return bpmap_.begin(); }
    const_iterator end()   const { return bpmap_.end();   }
    /// \return Array of dihedral parameters matching given atom types.
    DihedralParmArray FindParam(TypeNameHolder const& types, bool& found) const {
      found = true;
      for (const_iterator it = begin(); it != end(); ++it)
        if (it->first.Match_NoWC( types )) return it->second;
      if (wc_.len() > 0) {
        for (const_iterator it = begin(); it != end(); ++it)
          if (it->first.Match_WC( types, wc_)) return it->second;
      }
      found = false;
      return DihedralParmArray();
    }
    /// \return size in memory in bytes
    size_t DataSize() const {
      if (bpmap_.empty()) return 0;
      const_iterator elt0 = begin();
      // Assume all TypeNameHolders are the same size
      return (bpmap_.size() * elt0->first.DataSize()) +
             (bpmap_.size() * sizeof(DihedralParmArray)) +
             sizeof(Bmap);
    }
  protected:
    NameType wc_; ///< Wildcard character
  private:
    Bmap bpmap_;
};
// -----------------------------------------------------------------------------
/// Specialized class for associating atom types with improper parameters.
/** Impropers are a little tricky in that by convention the third atom is
  * the central atom, and all other atoms can be in any order.
  * The Amber convention is usually (but not always) to have the non-central
  * improper atom types sorted alphabetically, with wildcards given
  * precedence, but this is not always the case and does not always work.
  * For example, using straight up backwards/forwards matching, the wildcard
  * type X-X-CW-H4 will not match the given alphabetized type C*-H4-CW-NA.
  * All combinations of A1, A2, and A4 should be checked.
  */
class ImproperParmHolder : private DihedralParmHolder {
    /// Function for matching wildcards (WildCard Match)
    static inline bool wcm(NameType const& t0, NameType const& t1, NameType const& wc) {
      return (t0 == wc || t0 == t1);
    }
  public:
    /// Denote returned parameter atom type order
    enum OrderType { O_013,
                     O_031,
                     O_103,
                     O_130,
                     O_301,
                     O_310 };
    ImproperParmHolder() {}
    /// \return Number of improper parameter sets
    size_t size()       const { return DihedralParmHolder::size();  }
    /// \return True if no parameters
    bool empty()        const { return DihedralParmHolder::empty(); }
    /// \return ordering of last type
    //OrderType LastOrder() const { return lastOrder_; }
    typedef typename DihedralParmHolder::const_iterator const_iterator;
    const_iterator begin() const { return DihedralParmHolder::begin(); }
    const_iterator end()   const { return DihedralParmHolder::end();   }

    /** Set Wildcard char */
    void SetWildcard(char wc) { DihedralParmHolder::SetWildcard(wc); }
    /** Add (or update) a single improper parameter for given atom types. */
    ParameterHolders::RetType
    AddParm(TypeNameHolder const& types, DihedralParmType const& dp, bool allowUpdate) {
      return DihedralParmHolder::AddParm( types, dp, allowUpdate );
    }
    /** This version takes an array of dihedral parameters. */
    ParameterHolders::RetType
    AddParm(TypeNameHolder const& types, DihedralParmArray const& dpa, bool allowUpdate) {
      return DihedralParmHolder::AddParm( types, dpa, allowUpdate );
    }
    DihedralParmArray FindParam(TypeNameHolder const& types, bool& found) const {
      OrderType lastOrder;
      return FindParam(types, found, lastOrder);
    }
    /// \return Array of improper parameters matching given atom types.
    DihedralParmArray FindParam(TypeNameHolder const& types, bool& found, OrderType& lastOrder_) const {
      //mprintf("DEBUG: FindParam wc=%s Inco=%s-%s-%s-%s\n",*wc_, *(types[0]), *(types[1]),   *(types[2]),   *(types[3]));
      found = true;
      // First, no wildcard
      for (const_iterator it = begin(); it != end(); ++it) {
        TypeNameHolder const& myTypes = it->first;
        // Central (third) type must match
        if (myTypes[2] == types[2]) {
          //mprintf("DEBUG: FindParam (improper) central atom match %s", *(types[2]));
          //mprintf(" This=%s-%s-%s-%s", *(myTypes[0]), *(myTypes[1]), *(myTypes[2]), *(myTypes[3]));
          //mprintf(" Inco=%s-%s-%s-%s\n", *(types[0]), *(types[1]),   *(types[2]),   *(types[3]));
          // Try all permutations
          if (       myTypes[0] == types[0] && myTypes[1] == types[1] && myTypes[3] == types[3]) {
              // 0 1 2 3
              lastOrder_ = O_013;
              return it->second;
          } else if (myTypes[0] == types[0] && myTypes[1] == types[3] && myTypes[3] == types[1]) {
              // 0 3 2 1
              lastOrder_ = O_031;
              return it->second;
          } else if (myTypes[0] == types[1] && myTypes[1] == types[0] && myTypes[3] == types[3]) {
              // 1 0 2 3
              lastOrder_ = O_103;
              return it->second;
          } else if (myTypes[0] == types[1] && myTypes[1] == types[3] && myTypes[3] == types[0]) {
              // 1 3 2 0
              lastOrder_ = O_130;
              return it->second;
          } else if (myTypes[0] == types[3] && myTypes[1] == types[0] && myTypes[3] == types[1]) {
              // 3 0 2 1
              lastOrder_ = O_301;
              return it->second;
          } else if (myTypes[0] == types[3] && myTypes[1] == types[1] && myTypes[3] == types[0]) {
              // 3 1 2 0
              lastOrder_ = O_310;
              return it->second;
          }
        }
      } // END loop over parameters
      // Wildcard if present
      if (wc_.len() > 0) {
        for (const_iterator it = begin(); it != end(); ++it) {
          TypeNameHolder const& myTypes = it->first;
          // Central (third) type must match
          if (wcm(myTypes[2], types[2], wc_)) {
            // Try all permutations
            if (       wcm(myTypes[0], types[0], wc_) && wcm(myTypes[1], types[1], wc_) && wcm(myTypes[3], types[3], wc_)) {
                // 0 1 2 3
                lastOrder_ = O_013;
                return it->second;
            } else if (wcm(myTypes[0], types[0], wc_) && wcm(myTypes[1], types[3], wc_) && wcm(myTypes[3], types[1], wc_)) {
                // 0 3 2 1
                lastOrder_ = O_031;
                return it->second;
            } else if (wcm(myTypes[0], types[1], wc_) && wcm(myTypes[1], types[0], wc_) && wcm(myTypes[3], types[3], wc_)) {
                // 1 0 2 3
                lastOrder_ = O_103;
                return it->second;
            } else if (wcm(myTypes[0], types[1], wc_) && wcm(myTypes[1], types[3], wc_) && wcm(myTypes[3], types[0], wc_)) {
                // 1 3 2 0
                lastOrder_ = O_130;
                return it->second;
            } else if (wcm(myTypes[0], types[3], wc_) && wcm(myTypes[1], types[0], wc_) && wcm(myTypes[3], types[1], wc_)) {
                // 3 0 2 1
                lastOrder_ = O_301;
                return it->second;
            } else if (wcm(myTypes[0], types[3], wc_) && wcm(myTypes[1], types[1], wc_) && wcm(myTypes[3], types[0], wc_)) {
                // 3 1 2 0
                lastOrder_ = O_310;
                return it->second;
            }
          }
        } // END loop over parameters
      } // END wildcard matches
      found = false;
      return DihedralParmArray();
    } // END FindParam()
    /// \return size in memory in bytes
    size_t DataSize() const { return DihedralParmHolder::DataSize(); }
  private:
    //OrderType lastOrder_; ///< Atom ordering of the last returned parameter
};
#endif
