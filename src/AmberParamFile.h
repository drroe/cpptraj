#ifndef INC_AMBERPARAMFILE_H
#define INC_AMBERPARAMFILE_H
#include <vector>
#include <string>
class ParameterSet;
class CmapGridType;
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
    /// Set debug level
    void SetDebug(int);
  private:
    static const int MAXSYMLEN;

    enum SectionType { ATYPE = 0, HYDROPHILIC, BOND, ANGLE, DIHEDRAL, IMPROPER, 
                       LJ1012, NB_EQUIV, NONBOND, LJEDIT, CMAP, UNKNOWN };
    enum CmapType { CMAP_INITIAL, CMAP_TITLE, CMAP_RESLIST, CMAP_PARAMETER };

    class NonbondSet;
    class OffdiagNB;
    typedef std::vector<OffdiagNB> Oarray;

    static int read_symbols(const char*, std::vector<std::string>&, int);
    /// Read atom type line
    int read_atype(ParameterSet&, const char*) const;
    /// Read bond line
    int read_bond(ParameterSet&, const char*) const;
    /// Read angle line
    int read_angle(ParameterSet&, const char*) const;
    /// Read dihedral line
    int read_dihedral(ParameterSet&, const char*, std::vector<std::string>&, bool) const;
    /// Read improper line
    int read_improper(ParameterSet&, const char*) const;
    /// Read LJ 10-12 hbond line
    int read_lj1012(ParameterSet&, const char*) const;
    /// Read LJ 6-12 R/depth line
    int read_nb_RE(NonbondSet&, const char*) const;
    /// Read LJ 6-12 off-diagonal modifications
    int read_ljedit(Oarray&, const char*) const;
    /// Read CMAP section
    int read_cmap(CmapGridType&, ParameterSet&, CmapType&, std::string const&) const;
    /// Assign parameters from NonbondSet to ParameterSet
    int assign_nb(ParameterSet&, NonbondSet const&) const;
    /// Assign parameters from OffdiagNB array to ParameterSet
    int assign_offdiag(ParameterSet&, Oarray const&) const;
    //int ReadInput(std::string&, BufferedLine&) const;

    int debug_;
};
#endif