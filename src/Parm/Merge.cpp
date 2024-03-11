#include "Merge.h"
#include "../Atom.h"
#include "../CpptrajStdio.h"
#include "../ParameterHolders.h"
#include "../ParameterTypes.h"

/// Append bnd1 to bonds0 arrays along with parameters
static inline void append_bond(BondArray& bonds0,
                               BondArray& bondsh0,
                               BondParmArray& bp0,
                               unsigned int atomOffset,
                               BondType const& bnd1,
                               ParmHolder<int>& currentTypes0,
                               ParmHolder<int> const& currentTypes1,
                               BondParmArray const& bp1,
                               std::vector<Atom> const& atoms1)
{
  mprintf("DEBUG: Bond from top1 %i - %i will be %i - %i in top0\n",
          bnd1.A1()+1, bnd1.A2()+1, bnd1.A1()+1+atomOffset, bnd1.A2()+1+atomOffset);
  Atom const& A1 = atoms1[bnd1.A1()];
  Atom const& A2 = atoms1[bnd1.A2()];
  // Do we have an existing parameter in top0
  TypeNameHolder types(2);
  types.AddName( A1.Type() );
  types.AddName( A2.Type() );
  bool found;
  int idx = currentTypes0.FindParam(types, found);
  if (!found) {
    // No parameter yet.
    // Do we have an existing parameter in top1
    idx = currentTypes1.FindParam(types, found);
    if (!found) {
      // No parameter in either top.
      idx = -1;
    } else {
      // Found a parameter in top1, add it to top0.
      int newIdx = bp0.size();
      bp0.push_back( bp1[idx] );
      idx = newIdx;
    }
    // Add to existing parameters in top0.
    // Do this even if a parameter was not found so we dont keep looking.
    if (idx < 0)
      mprintf("Warning: No bond parameters for types %s - %s\n", *(A1.Type()), *(A2.Type()));
    currentTypes0.AddParm(types, idx, false);
  }
  // At this point we have either found a parameter or not.
  if (A1.Element() == Atom::HYDROGEN ||
      A2.Element() == Atom::HYDROGEN)
    bondsh0.push_back( BondType(bnd1.A1()+atomOffset, bnd1.A2()+atomOffset, idx) );
  else
    bonds0.push_back( BondType(bnd1.A1()+atomOffset, bnd1.A2()+atomOffset, idx) );
}

/// Index existing bond types in bond arrays
static inline void index_bond_types(ParmHolder<int>& currentTypes,
                                    BondArray const& bonds,
                                    std::vector<Atom> const& atoms)
{
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    if (bnd->Idx() > -1) {
      TypeNameHolder types(2);
      types.AddName( atoms[bnd->A1()].Type() );
      types.AddName( atoms[bnd->A2()].Type() );
      bool found;
      currentTypes.FindParam(types, found);
      if (!found) {
        currentTypes.AddParm(types, bnd->Idx(), false);
      }
    }
  }
}

/** Given bond/bond parameter arrays from top0 and bond/bond parameter
  * arrays from top1, merge the bond arrays and consolidate the
  * parameters.
  */
void Cpptraj::Parm::MergeBondArrays(BondArray& bonds0,
                                    BondArray& bondsh0,
                                    BondParmArray& bp0,
                                    AtArray const& atoms0,
                                    BondArray const& bonds1,
                                    BondArray const& bondsh1,
                                    BondParmArray const& bp1,
                                    AtArray const& atoms1)
{
  // First index existing parameters
  ParmHolder<int> currentTypes0, currentTypes1;
  index_bond_types(currentTypes0, bonds0, atoms0);
  index_bond_types(currentTypes0, bondsh0, atoms0);
  index_bond_types(currentTypes1, bonds1, atoms1);
  index_bond_types(currentTypes1, bondsh1, atoms1);

  // Loop over separate bond arrays from top1 in the correct order
  unsigned int atomOffset = atoms0.size();
  BondArray::const_iterator bx = bonds1.begin();
  BondArray::const_iterator by = bondsh1.begin();
  while (bx != bonds1.end() && by != bondsh1.end()) {
    // Which one goes next?
    Atom const& bx0 = atoms1[bx->A1()];
    Atom const& by0 = atoms1[by->A1()];
    //bool has_parm_indices = (bx->Idx() != -1 && by->Idx() != -1);
    if (bx0.ResNum() == by0.ResNum()) {
      if (bx->A1() == by->A1()) {
        // Same residue, same A1. Lower A2 goes first.
        //if (by->A2() < bx->A2()) {
        if (*by < *bx) {
          append_bond( bonds0, bondsh0, bp0, atomOffset, *by, currentTypes0, currentTypes1, bp1, atoms1 );
          ++by;
        } else {
          append_bond( bonds0, bondsh0, bp0, atomOffset, *bx, currentTypes0, currentTypes1, bp1, atoms1 );
          ++bx;
        }
      } else {
        // Both bonds in same residue, different A1.
        // Higher A1 goes first.
        // FIXME fix for scan direction forwards.
        if (by->A1() > bx->A1()) {
          append_bond( bonds0, bondsh0, bp0, atomOffset, *by, currentTypes0, currentTypes1, bp1, atoms1 );
          ++by;
        } else {
          append_bond( bonds0, bondsh0, bp0, atomOffset, *bx, currentTypes0, currentTypes1, bp1, atoms1 );
          ++bx;
        }
      }
    } else {
      // Lower residue goes first.
      if (by0.ResNum() < bx0.ResNum()) {
        append_bond( bonds0, bondsh0, bp0, atomOffset, *by, currentTypes0, currentTypes1, bp1, atoms1 );
        ++by;
      } else {
        append_bond( bonds0, bondsh0, bp0, atomOffset, *bx, currentTypes0, currentTypes1, bp1, atoms1 );
        ++bx;
      }
    }
  } // END loop over both bond arrays from top1
  if (bx != bonds1.end()) {
    for (; bx != bonds1.end(); ++bx)
      append_bond( bonds0, bondsh0, bp0, atomOffset, *bx, currentTypes0, currentTypes1, bp1, atoms1 );
  }
  if (by != bondsh1.end()) {
    for (; by != bondsh1.end(); ++by)
      append_bond( bonds0, bondsh0, bp0, atomOffset, *by, currentTypes0, currentTypes1, bp1, atoms1 );
  }
}

// -----------------------------------------------------------------------------
/// Append ang1 to angles0 arrays along with parameters
static inline void append_angle(AngleArray& angles0,
                               AngleArray& anglesh0,
                               AngleParmArray& ap0,
                               unsigned int atomOffset,
                               AngleType const& ang1,
                               ParmHolder<int>& currentTypes0,
                               ParmHolder<int> const& currentTypes1,
                               AngleParmArray const& ap1,
                               std::vector<Atom> const& atoms1)
{
  mprintf("DEBUG: Angle from top1 %i - %i  - %i will be %i - %i - %i in top0\n",
          ang1.A1()+1, ang1.A2()+1, ang1.A3()+1, ang1.A1()+1+atomOffset, ang1.A2()+1+atomOffset, ang1.A3()+1+atomOffset);
  Atom const& A1 = atoms1[ang1.A1()];
  Atom const& A2 = atoms1[ang1.A2()];
  Atom const& A3 = atoms1[ang1.A3()];
  // Do we have an existing parameter in top0
  TypeNameHolder types(3);
  types.AddName( A1.Type() );
  types.AddName( A2.Type() );
  types.AddName( A3.Type() );
  bool found;
  int idx = currentTypes0.FindParam(types, found);
  if (!found) {
    // No parameter yet.
    // Do we have an existing parameter in top1
    idx = currentTypes1.FindParam(types, found);
    if (!found) {
      // No parameter in either top.
      idx = -1;
    } else {
      // Found a parameter in top1, add it to top0.
      int newIdx = ap0.size();
      ap0.push_back( ap1[idx] );
      idx = newIdx;
    }
    // Add to existing parameters in top0.
    // Do this even if a parameter was not found so we dont keep looking.
    if (idx < 0)
      mprintf("Warning: No angle parameters for types %s - %s - %s\n", *(A1.Type()), *(A2.Type()), *(A3.Type()));
    currentTypes0.AddParm(types, idx, false);
  }
  // At this point we have either found a parameter or not.
  if (A1.Element() == Atom::HYDROGEN ||
      A2.Element() == Atom::HYDROGEN ||
      A3.Element() == Atom::HYDROGEN)
    anglesh0.push_back( AngleType(ang1.A1()+atomOffset, ang1.A2()+atomOffset, ang1.A3()+atomOffset, idx) );
  else
    angles0.push_back( AngleType(ang1.A1()+atomOffset, ang1.A2()+atomOffset, ang1.A3()+atomOffset, idx) );
}

/// Index existing angle types in angle arrays
static inline void index_angle_types(ParmHolder<int>& currentTypes,
                                    AngleArray const& angles,
                                    std::vector<Atom> const& atoms)
{
  for (AngleArray::const_iterator ang = angles.begin(); ang != angles.end(); ++ang)
  {
    if (ang->Idx() > -1) {
      TypeNameHolder types(3);
      types.AddName( atoms[ang->A1()].Type() );
      types.AddName( atoms[ang->A2()].Type() );
      types.AddName( atoms[ang->A3()].Type() );
      bool found;
      currentTypes.FindParam(types, found);
      if (!found) {
        currentTypes.AddParm(types, ang->Idx(), false);
      }
    }
  }
}

/** Given angle/angle parameter arrays from top0 and angle/angle parameter
  * arrays from top1, merge the angle arrays and consolidate the
  * parameters.
  */
void Cpptraj::Parm::MergeAngleArrays(AngleArray& angles0,
                                    AngleArray& anglesh0,
                                    AngleParmArray& bp0,
                                    AtArray const& atoms0,
                                    AngleArray const& angles1,
                                    AngleArray const& anglesh1,
                                    AngleParmArray const& bp1,
                                    AtArray const& atoms1)
{
  // First index existing parameters
  ParmHolder<int> currentTypes0, currentTypes1;
  index_angle_types(currentTypes0, angles0, atoms0);
  index_angle_types(currentTypes0, anglesh0, atoms0);
  index_angle_types(currentTypes1, angles1, atoms1);
  index_angle_types(currentTypes1, anglesh1, atoms1);

  // Loop over separate angle arrays from top1 in the correct order
  unsigned int atomOffset = atoms0.size();
  AngleArray::const_iterator bx = angles1.begin();
  AngleArray::const_iterator by = anglesh1.begin();
  while (bx != angles1.end() && by != anglesh1.end()) {
    // Which one goes next?
    Atom const& bx0 = atoms1[bx->A1()];
    Atom const& by0 = atoms1[by->A1()];
    //bool has_parm_indices = (bx->Idx() != -1 && by->Idx() != -1);
    if (bx0.ResNum() == by0.ResNum()) {
      if (bx->A1() == by->A1()) {
        // Same residue, same A1. Lower A2 goes first.
        //if (by->A2() < bx->A2()) {
        if (*by < *bx) {
          append_angle( angles0, anglesh0, bp0, atomOffset, *by, currentTypes0, currentTypes1, bp1, atoms1 );
          ++by;
        } else {
          append_angle( angles0, anglesh0, bp0, atomOffset, *bx, currentTypes0, currentTypes1, bp1, atoms1 );
          ++bx;
        }
      } else {
        // Both angles in same residue, different A1.
        // Higher A1 goes first.
        // FIXME fix for scan direction forwards.
        if (by->A1() > bx->A1()) {
          append_angle( angles0, anglesh0, bp0, atomOffset, *by, currentTypes0, currentTypes1, bp1, atoms1 );
          ++by;
        } else {
          append_angle( angles0, anglesh0, bp0, atomOffset, *bx, currentTypes0, currentTypes1, bp1, atoms1 );
          ++bx;
        }
      }
    } else {
      // Lower residue goes first.
      if (by0.ResNum() < bx0.ResNum()) {
        append_angle( angles0, anglesh0, bp0, atomOffset, *by, currentTypes0, currentTypes1, bp1, atoms1 );
        ++by;
      } else {
        append_angle( angles0, anglesh0, bp0, atomOffset, *bx, currentTypes0, currentTypes1, bp1, atoms1 );
        ++bx;
      }
    }
  } // END loop over both angle arrays from top1
  if (bx != angles1.end()) {
    for (; bx != angles1.end(); ++bx)
      append_angle( angles0, anglesh0, bp0, atomOffset, *bx, currentTypes0, currentTypes1, bp1, atoms1 );
  }
  if (by != anglesh1.end()) {
    for (; by != anglesh1.end(); ++by)
      append_angle( angles0, anglesh0, bp0, atomOffset, *by, currentTypes0, currentTypes1, bp1, atoms1 );
  }
}


