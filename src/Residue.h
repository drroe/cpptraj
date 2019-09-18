#ifndef INC_RESIDUE_H
#define INC_RESIDUE_H
#include <map> // for Residue name map
#include "NameType.h"
// Class: Residue
/// Hold information for a residue.
class Residue {
  public:
    /// Residue types.
    enum ResidueType { PROTEIN = 0, NUCLEIC, LIPID, SOLVENT, UNKNOWN };

    /// CONSTRUCTOR
    Residue() :
      resname_(""), type_(UNKNOWN), firstAtom_(0), lastAtom_(0), originalResNum_(0), segID_(-1),
      icode_(' '), chainID_(BLANK_CHAINID_), isTerminal_(false)
    {}
    /// CONSTRUCTOR - Copy given Residue, set first and last atom indices.
    Residue(Residue const& r, int first, int last) :
      resname_(r.resname_), type_(r.type_), firstAtom_(first), lastAtom_(last),
      originalResNum_(r.originalResNum_), segID_(r.segID_), icode_(r.icode_),
      chainID_(r.chainID_), isTerminal_(false)
    {}
    /// CONSTRUCTOR - Res name, original resnum, icode, chain ID
    Residue(NameType const& n, int r, char ic, char cid) :
      resname_(n), type_(GetTypeFromName(n)), firstAtom_(-1), lastAtom_(-1), originalResNum_(r), segID_(-1),
      icode_(ic), chainID_(cid), isTerminal_(false)
    {}
    /// CONSTRUCTOR - Res name, first atom, last atom, original resnum, icode, chain ID
    Residue(NameType const& n, int first, int last, int r, char ic, char cid) :
      resname_(n), type_(GetTypeFromName(n)), firstAtom_(first), lastAtom_(last),
      originalResNum_(r), segID_(-1), icode_(ic), chainID_(cid),
      isTerminal_(false)
    {}
    /// CONSTRUCTOR - Res name, original resnum, res icode, segment ID
    Residue(NameType const& n, int r, char i, int s) :
      resname_(n), type_(GetTypeFromName(n)), firstAtom_(-1), lastAtom_(-1), originalResNum_(r), segID_(s),
       icode_(i), chainID_(BLANK_CHAINID_), isTerminal_(false)
    {}
    void SetFirstAtom(int i)        { firstAtom_ = i;      }
    void SetLastAtom(int i)         { lastAtom_ = i;       }
    void SetOriginalNum(int i)      { originalResNum_ = i; }
    void SetSegID(int s)            { segID_ = s;          }
    void SetIcode(char c)           { icode_ = c;          }
    void SetChainID(char c)         { chainID_ = c;        }
    void SetName(NameType const& n) { resname_ = n; type_ = GetTypeFromName(n); }
    void SetTerminal(bool t)        { isTerminal_ = t;     }
    /// \return First atom in residue, indexing from 0
    int FirstAtom()        const { return firstAtom_;      }
    /// \return Atom _after_ the last in residue, indexing from 0
    int LastAtom()         const { return lastAtom_;       }
    int OriginalResNum()   const { return originalResNum_; }
    int SegID()            const { return segID_;          }
    char Icode()           const { return icode_;          }
    inline char ChainID()         const;
    bool HasChainID()      const { return (chainID_ != BLANK_CHAINID_); }
    const char *c_str()    const { return *resname_;       }
    NameType const& Name() const { return resname_;        }
    ResidueType Type()     const { return type_; }
    bool IsTerminal()      const { return isTerminal_;     }
    int NumAtoms()         const { return (lastAtom_ - firstAtom_); }
    /// Convert 3-letter residue code to single letter.
    static char ConvertResName(std::string const&);
    /// Convert 1-letter residue code to 3 letters.
    static const char* ConvertResName(char);
    /// \return Default chain ID
    static char DefaultChainID() { return DEFAULT_CHAINID_; }
    /// \return Blank (no) chain ID
    static char BlankChainID() { return BLANK_CHAINID_; }
    /// Convert this residue name to single letter.
    char SingleCharName() const { return ConvertResName( *resname_ ); }
    /// Initialize the residue-type name map. Called from Cpptraj constructor.
    static void InitResNameMap();
    /// \return string corresponding to ResidueType
    static const char* ResTypeStr(ResidueType);
  private:
    /// /return Residue type from given residue name
    static ResidueType GetTypeFromName(NameType const&);

    /// Used to map recognized residue names to types
    typedef std::map<NameType,ResidueType> ResNameMapType;
    /// Used to pair residue names to types
    typedef std::pair<NameType,ResidueType> ResNamePairType;
    /// Map residue names to their associated type
    static ResNameMapType resNameMap_;

    static const char* resTypeStr_[];
    static const char BLANK_CHAINID_;
    static const char DEFAULT_CHAINID_;
    NameType resname_;   ///< Residue name.
    ResidueType type_;   ///< Residue type based on name.
    int firstAtom_;      ///< Index of first atom (from 0).
    int lastAtom_;       ///< Atom index after last atom in residue.
    int originalResNum_; ///< Original residue number.
    int segID_;          ///< Segment ID index.
    char icode_;         ///< Residue insertion code.
    char chainID_;       ///< Residue chain ID
    bool isTerminal_;    ///< True if residue was originally a terminal residue
};
// ----- INLINE ROUTINES -------------------------
char Residue::ChainID() const {
  if (chainID_ == BLANK_CHAINID_)
    return ' ';
  else
    return chainID_;
}
#endif
