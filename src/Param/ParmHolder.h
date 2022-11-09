#ifndef INC_PARAM_PARMHOLDER_H
#define INC_PARAM_PARMHOLDER_H
#include <vector>
#include <utility> // std::pair
#include "Param.h"
#include "../TypeNameHolder.h"
namespace Cpptraj {
namespace Param {
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
    /// Add (or update if allowed) given parameter to holder.
    RetType AddParm(TypeNameHolder const& types, T const& bp, bool allowUpdate) {
      // Check if parm for these types exist
      typename Bmap::iterator it = bpmap_.begin();
      for (; it != bpmap_.end(); ++it)
        if (it->first == types) break;
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
            return UPDATED;
          } else {
            return ERR;
          }
        } else {
          //mprintf("DEBUG: Existing parameter:");
          //for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
          //  mprintf(" '%s'", *(*it));
          //mprintf("\n");
          return SAME;
        }
      }
      return ADDED;
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
        if (it->first == types) return it->second;
      found = false;
      return T();
    }
    /// \return iterator to parameter matching the given types.
    iterator GetParam(TypeNameHolder const& types) {
      for (iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first == types) return it;
      return bpmap_.end();
    }
    /// \return const iterator to parameter matching the given types.
    const_iterator GetParam(TypeNameHolder const& types) const {
      for (const_iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first == types) return it;
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
};
}
}
#endif
