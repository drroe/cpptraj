#ifndef INC_AMBERPARAMFILE_H
#define INC_AMBERPARAMFILE_H
#include <vector>
#include <string>
class ParameterSet;
//class BufferedLine;
class FileName;
/// Used to read in Amber parameters from Amber FF/FRCMOD file.
class AmberParamFile {
  public:
    AmberParamFile();
    /// Read main Amber FF file
    int ReadParams(ParameterSet&, FileName const&, std::string const&, int) const;
    /// Read Amber frcmod file
    int ReadFrcmod(ParameterSet&, FileName const&, int) const;
    /// Write main Amber FF file
    int WriteParams(ParameterSet&, FileName const&, int) const;
  private:
    static const int MAXSYMLEN;

    enum SectionType { ATYPE = 0, HYDROPHILIC, BOND, ANGLE, DIHEDRAL, IMPROPER, 
                       LJ1012, NB_EQUIV, NONBOND, LJEDIT, UNKNOWN };

    class NonbondSet;

    static int read_symbols(const char*, std::vector<std::string>&, int);
    /// Read atom type line
    int read_atype(ParameterSet&, const char*) const;
    /// Read bond line
    int read_bond(ParameterSet&, const char*) const;
    /// Read angle line
    int read_angle(ParameterSet&, const char*) const;
    /// Read dihedral line
    int read_dihedral(ParameterSet&, const char*) const;
    /// Read improper line
    int read_improper(ParameterSet&, const char*) const;
    /// Read LJ 10-12 hbond line
    int read_lj1012(ParameterSet&, const char*) const;
    /// Read LJ 6-12 R/depth line
    int read_nb_RE(NonbondSet&, const char*) const;
    //int ReadInput(std::string&, BufferedLine&) const;
};
#endif
