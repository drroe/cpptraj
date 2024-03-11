#include "Merge.h"
#include "../Atom.h"
#include "../CpptrajStdio.h"
#include "../ParameterTypes.h"

static inline void append_bond(BondArray& bonds0,
                               BondArray& bondsh0,
                               BondParmArray& bp0,
                               unsigned int atomOffset,
                               BondType const& bnd1)
{
  mprintf("DEBUG: Bond from top1 %i - %i will be %i - %i in top0\n",
          bnd1.A1()+1, bnd1.A2()+1, bnd1.A1()+1+atomOffset, bnd1.A2()+1+atomOffset);
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
        if (by->A2() < bx->A2()) {
          append_bond( bonds0, bondsh0, bp0, atomOffset, *by );
          ++by;
        } else {
          append_bond( bonds0, bondsh0, bp0, atomOffset, *bx );
          ++bx;
        }
      } else {
        // Both bonds in same residue, different A1.
        // Higher A1 goes first.
        // FIXME fix for scan direction forwards.
        if (by->A1() > bx->A1()) {
          append_bond( bonds0, bondsh0, bp0, atomOffset, *by );
          ++by;
        } else {
          append_bond( bonds0, bondsh0, bp0, atomOffset, *bx );
          ++bx;
        }
      }
    } else {
      // Lower residue goes first.
      if (by0.ResNum() < bx0.ResNum()) {
        append_bond( bonds0, bondsh0, bp0, atomOffset, *by );
        ++by;
      } else {
        append_bond( bonds0, bondsh0, bp0, atomOffset, *bx );
        ++bx;
      }
    }
  } // END loop over both bond arrays from top1

}
