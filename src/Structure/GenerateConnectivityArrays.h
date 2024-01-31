#ifndef INC_STRUCTURE_GENERATECONNECTIVITYARRAYS_H
#define INC_STRUCTURE_GENERATECONNECTIVITYARRAYS_H
#include <vector>
class BondArray;
class AngleArray;
class DihedralArray;
class Residue;
class Atom;
namespace Cpptraj {
namespace Structure {
/// Specify direction in which atoms in residues should be scanned
enum AtomScanDirectionType { SCAN_ATOMS_BACKWARDS = 0, SCAN_ATOMS_FORWARDS };
/// Set default atom scan direction
void SetAtomScanDirection(AtomScanDirectionType);
/// Generate bond array in same order as leap
BondArray GenerateBondArray(std::vector<Residue> const&, std::vector<Atom> const&);
/// Generate angle array in same order as leap
AngleArray GenerateAngleArray(std::vector<Residue> const&, std::vector<Atom> const&);
/// Generate proper dihedral array in same order as leap
DihedralArray GenerateDihedralArray(std::vector<Residue> const&, std::vector<Atom> const&);
/// Generate improper dihedral array in same order as leap
DihedralArray GenerateImproperArray(std::vector<Residue> const&, std::vector<Atom> const&);
/// Generate angle and torsion arrays from bonds
void GenerateAngleAndTorsionArraysFromBonds(AngleArray&, DihedralArray&, std::vector<Atom> const&, BondArray const&);
/// Generate array of dihedrals to be used as internal coordinates
DihedralArray GenerateInternals(std::vector<Residue> const&, std::vector<Atom> const&);
}
}
#endif
