#ifndef INC_EXEC_FIXATOMNAMES_H
#define INC_EXEC_FIXATOMNAMES_H
#include "Exec.h"
//#incl ude <map>
/// Fix topology so all atom names in residues are unique.
class Exec_FixAtomNames : public Exec {
  public:
    Exec_FixAtomNames() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_FixAtomNames(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
/*    typedef std::pair<NameType, NameType> NamePair;
    typedef std::map<NameType, NameType> NameMap;
    typedef std::pair<NameType, NameMap> ResPair;
    typedef std::map<NameType, NameMap> ResMap;

    ResMap resMap_; ///< Hold a map of residues that have had their atom names changed*/
    typedef std::vector<int> Iarray;
    /// Hold a part of an atom name
    class NameToken {
      public:
        enum Type { STRING = 0, ///< A letter (A-Z, a-z)
                    INT,        ///< A digit (0-9)
                    OTHER,      ///< Any other character. Will not change.
                    UNKNOWN };
        NameToken(std::string const& s) : ival_(-1), sval_(s), type_(STRING) {}
        NameToken(int i) : ival_(i), sval_(""), type_(INT) {}
        NameToken(Type t, std::string const& s) : ival_(-1), sval_(s), type_(t) {}

        int Ival() const { return ival_; }
        std::string const& Sval() const { return sval_; }
        Type TokenType() const { return type_; }

        void SetIval(int i) { ival_ = i; }
      private:
        int ival_;
        std::string sval_;
        Type type_;
    };
    /// Hold all parts that make up an atom name
    class AtomName {
      public:
        AtomName() {}
        //void AddToken(NameToken const& n) { tokens_.push_back( n ); }
        void PrintName() const;
        std::string NameString() const;
        int DissectName(NameType const&);
        int NameWidth() const;
        /// Naively increment the name
        int IncrementName(int);
      private:
        typedef std::vector<NameToken> NameTokenArray;
        NameTokenArray tokens_;
    };
    /// Hold info for a duplicated atom
    class DupAtom {
      public:
        DupAtom() {}
      private:
        int idx_; ///< Atom index in topology
        int sum_; ///< Used to prioritize duplicated atoms; sum of atomic #s of bonded atoms
    };
    /// Used for sorting the duplicated name map
    class DupName {
      public:
        DupName(NameType const& n, Iarray const& i, int w) : weight_(w), aname_(n), idxs_(i) {}
        int weight_; ///< 9999 for C, atomic number for others
        NameType aname_;
        Iarray idxs_; ///< Atom indices in topology

        bool operator<(const DupName& rhs) const {
          if (weight_ == rhs.weight_) {
            if (aname_ == rhs.aname_) {
              return (idxs_.front() < rhs.idxs_.front());
            } else {
              return (aname_ < rhs.aname_);
            }
          } else {
            return (weight_ > rhs.weight_);
          }
        }
    };
        

    static inline NameToken::Type get_type(char);

    int fixDuplicate(Topology&, int, int) const;

    int fixDuplicateSet(Topology&, Iarray const&) const;
};
#endif
