#ifndef INC_EXEC_SEQUENCEALIGN_H
#define INC_EXEC_SEQUENCEALIGN_H
#include "Exec.h"
#include <vector>
// EXPERIMENTAL ALPHA CODE
class Exec_SequenceAlign : public Exec {
  public:
    /// CONSTRUCTOR
    Exec_SequenceAlign();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_SequenceAlign(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<char> Carray;
    /// Used to hold corresponding residue character and mapped residue number/character
    class ResNumChar {
      public:
        ResNumChar() : resChar_(' '), mappedNum_(-1), mappedChar_(' ') {}
        ResNumChar(char rc, int mn, char mc) :
          resChar_(rc), mappedNum_(mn), mappedChar_(mc) {}
        char ResChar() const { return resChar_; }
        int MappedNum() const { return mappedNum_; }
        char MappedChar() const { return mappedChar_; }
      private:
        char resChar_; ///< This residues character
        int mappedNum_; ///< Residue # this residue is mapped to.
        char mappedChar_; ///< Character of residue this residue is mapped to.
    };

    typedef std::pair<int, ResNumChar> ResNumCharPair;
    typedef std::map<int, ResNumChar> ResNumCharMap;
/*    class ResMap {
      public:
        ResMap() : resNum_(-1), resChar_(' '), mappedNum_(-1), mappedChar_(' ') {}
        ResMap(int rn, char rc, int mn, char mc) :
          resNum_(rn), resChar_(rc), mappedNum_(mn), mappedChar_(mc) {}
        int ResNum() const { return resNum_; }
        char ResChar() const { return resChar_; }
        int MappedNum() const { return mappedNum_; }
        char MappedChar() const { return mappedChar_; }
      private:
        int resNum_;   ///< This residues number
        char resChar_; ///< This residues character
        int mappedNum_; ///< Residue # this residue is mapped to.
        char mappedChar_; ///< Character of residue this residue is mapped to.
    };
    typedef std::vector<ResMap> ResNumCharMap;*/

    /// Read BLAST alignment, generate query/sbjct residue strings and query->sbjct map
    int read_blast(std::string&, std::string&, ResNumCharMap&, std::string const&);
    /// \return any residue number after advancing to string of residues
    static int advance_past_rnum(std::string::const_iterator&, std::string const&);
    /// Build subject from query using mapping
    int buildFromReference(Topology&, Frame&, ReferenceFrame const&, ResNumCharMap const&) const;

    int queryOffset_; ///< Query residue offset
    int sbjctOffset_; ///< Sbjct residue offset
    Iarray SBJCT_NUMS_;
    Carray SBJCT_CHAR_; 
};
#endif
