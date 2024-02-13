#ifndef INC_STRUCTURE_CHIRALITY_H
#define INC_STRUCTURE_CHIRALITY_H
class Atom;
class Vec3;
namespace Cpptraj {
namespace Structure {
namespace Chirality {

/// Calculate chirality in same manner as LEaP. All atom positions should be known.
double VectorAtomChirality(Vec3 const&, Vec3 const&, Vec3 const&, Vec3 const&);
/// Calculate chirality in same manner as LEaP. Some atom positions may not be known.
double VectorAtomNormalizedChirality(Vec3 const&,
                                     Vec3 const&, bool, Vec3 const&, bool,
                                     Vec3 const&, bool, Vec3 const&, bool);
/// Order atoms for chirality calculation like LEaP
void chiralityOrderNeighbors(Atom const&, int&, int&, int&, int&);
/// Transform given chirality to an orientation
double chiralityToOrientation(double, Atom const&, int, int, int, int);

}
}
}
#endif
