#include "Chirality.h"
#include "../Atom.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../Vec3.h"
#include <cmath> // fabs

/** LEaP routine for determining atom chirality.
  * This is done by crossing A to B and then dotting the
  * result with C. TODO use Chirality in BuildAtom? Use an Enum?
  * The chirality of the vectors is determined by the sign of
  * the result, which is determined by whether or not C has
  * a component in the direction AxB or in the opposite direction.
  */
double Cpptraj::Structure::Chirality::VectorAtomChirality(Vec3 const& Center,
                                                          Vec3 const& A, Vec3 const& B, Vec3 const& C)
{
  Vec3 vA = A - Center;
  Vec3 vB = B - Center;
  Vec3 vC = C - Center;
  Vec3 vCross = vA.Cross( vB );
  double dot = vCross * vC;
  if (dot > 0)
    return 1.0;
  else if (dot < 0)
    return -1.0;
  return 0.0;
}

/** LEaP routine for determining atom chirality when positions may or
  * may not be defined. Currently the criteria for chirality is 
  * absolute orientation of the vectors joining this atom to its neighbors.
  * The neighbors are passed as vPA, vPB, vPC, vPD and bA, bB, bC, bD
  * define whether or not the position is defined.
  *
  * This routine calculates the chirality using the defined vectors and
  * then flips the sign depending on which vectors were used to calculate
  * the chirality. If the ATOMs have (fKnown) set then their coordinate
  * is considered to be known.
  */
double Cpptraj::Structure::Chirality::VectorAtomNormalizedChirality(Vec3 const& Center,
                                                                    Vec3 const& vPA, bool bA,
                                                                    Vec3 const& vPB, bool bB,
                                                                    Vec3 const& vPC, bool bC,
                                                                    Vec3 const& vPD, bool bD)
{
  double dChi = 0;
  
  if (!bA) {
    // If A is not known then use B,C,D to calc chirality.
    // The chirality calculated will be negative w.r.t. the
    // correct chirality.
    if (!bB || !bC || !bD) return dChi;
    dChi = -VectorAtomChirality( Center, vPB, vPC, vPD );
    return dChi;
  }

  if (!bB) {
    // If B is not known then use A,C,D to calc chirality.
    // The chirality calculated will be correct.
    if (!bB || !bD) return dChi;
    dChi = VectorAtomChirality( Center, vPA, vPC, vPD );
    return dChi;
  }

  if (!bC) {
    // If C is not known then use A,B,D to calc chirality.
    // The chirality calculated will be negative w.r.t. the
    // correct chirality.
    if (!bD) return dChi;
    dChi = -VectorAtomChirality( Center, vPA, vPB, vPD );
    return dChi;
  }

  dChi = VectorAtomChirality( Center, vPA, vPB, vPC );

  return dChi;
}

/** \return index of atom less than all others but larger than aLast */
static inline int findLeastLargerThan(Atom const& aAtom, int aLast)
{
  int aSmall = -1;
  for (Atom::bond_iterator aCur = aAtom.bondbegin(); aCur != aAtom.bondend(); ++aCur)
  {
    if (aLast != -1) {
      if (aLast >= *aCur) continue;
    }
    if (aSmall == -1)
      aSmall = *aCur;
    else if ( *aCur < aSmall )
      aSmall = *aCur;
  }
  return aSmall;
}

/** Transform the orientation that has been measured with 
 *      respect to the ordering in aaOrig[4] to the ordering
 *      in aaNew[4].  Return the result.
 *
 *      The transformation is done by swapping ATOMs in aaOrig until
 *      the order matches that of aaNew, each time two ATOMs are swapped,
 *      flip the sign of the orientation.
 *
 *      SIDE EFFECT:   The order in aaOrig is changed.
  */
static inline void chiralityTransformOrientation(double dOrig, int* aaOrig, double& dPNew, const int* aaNew)
{
  dPNew = dOrig;
  for (int i=0; i<4; i++ ) {
    int j = i;
    for ( ; j<4; j++ ) {
      if ( aaOrig[j] == aaNew[i] )
        break;
    }
    if ( j >= 4 ) {
      mprinterr("Error: Comparing atoms %i %i %i and %i to atoms %i %i %i and %i.\n",
                aaOrig[0]+1, aaOrig[1]+1, aaOrig[2]+1, aaOrig[3]+1,
                aaNew[0]+1, aaNew[1]+1, aaNew[2]+1, aaNew[3]+1);
      mprinterr("Error: This error may be due to faulty Connection atoms.\n");
      // TODO fatal
    }
    // Swap elements and flip sign
    if ( j != i ) {
      std::swap( aaOrig[j], aaOrig[i] );
      dPNew = -dPNew;
    }
  }
}

/** Sort neighbors of given atom in the same manner as LEaP. */
void Cpptraj::Structure::Chirality::chiralityOrderNeighbors(Atom const& aAtom,
                                                            int& aPAtomA, int& aPAtomB,
                                                            int& aPAtomC, int& aPAtomD)
{
  aPAtomA = -1;
  aPAtomB = -1;
  aPAtomC = -1;
  aPAtomD = -1;

  if (aAtom.Nbonds() < 1) return;

  aPAtomA = findLeastLargerThan(aAtom, -1);
  if (aAtom.Nbonds() < 2) return;

  aPAtomB = findLeastLargerThan(aAtom, aPAtomA);
  if (aAtom.Nbonds() < 3) return;

  aPAtomC = findLeastLargerThan(aAtom, aPAtomB);
  if (aAtom.Nbonds() < 4) return;

  aPAtomD = findLeastLargerThan(aAtom, aPAtomC);
}

/** Transform the chirality which has been measured with
 *  respect to ATOM ID ordering to an arbitrary ordering.
 */
double Cpptraj::Structure::Chirality::chiralityToOrientation(double dChirality, Atom const& aCenter,
                                            int aAtomA, int aAtomB, int aAtomC, int aAtomD)
{
  if (fabs(dChirality) < Constants::SMALL) return 0.0;

  int aaOrig[4];
  chiralityOrderNeighbors( aCenter, aaOrig[0], aaOrig[1], aaOrig[2], aaOrig[3] );

  int aaNew[4];
  aaNew[0] = aAtomA;
  aaNew[1] = aAtomB;
  aaNew[2] = aAtomC;
  aaNew[3] = aAtomD;

  bool newNull = (aaNew[3] == -1);
  bool origNull = (aaOrig[3] == -1);
  if (newNull && !origNull) {
    for (int i = 0; i < 4; i++) {
      bool found = false;
      for (int j=0; j<3; j++) found |= (aaOrig[i] == aaNew[j]);
      if ( !found ) {
        aaNew[3] = aaOrig[i];
        break;
      }
    }
  } else if (!newNull && origNull) {
    mprinterr("Error: Only three neighbors around: aCenter, but orientation has 4\n");
  }

  double dOrient;
  chiralityTransformOrientation( dChirality, aaOrig, dOrient, aaNew );

  return dOrient;
}


