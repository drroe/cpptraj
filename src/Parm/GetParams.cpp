#include "GetParams.h"
#include "../Atom.h"
#include "../AtomType.h"
#include "../CpptrajStdio.h"
#include "../ParameterHolders.h"
#include "../ParameterTypes.h"
#include "../TypeNameHolder.h"

/** \param atomTypesOut Output array of atom types and indivudual LJ parameters.
  * \param LJ612out Output array of LJ 6-12 pair parameters.
  * \param LJ14out Output array of LJ 6-12 1-4 pair parameters.
  * \param LJ1012out Output array of LJ 10-12 pair parameters.
  * \param atoms Current array of atoms.
  * \param NB0 Current nonbond parameters.
  */
void Cpptraj::Parm::GetParams::GetLJAtomTypes(ParmHolder<AtomType>& atomTypesOut,
                                  ParmHolder<NonbondType>& LJ612out,
                                  ParmHolder<NonbondType>& LJ14out,
                                  ParmHolder<double>& LJCout,
                                  ParmHolder<HB_ParmType>& LJ1012out,
                                  std::vector<Atom> const& atoms,
                                  NonbondParmType const& NB0)
const
{
  if (NB0.HasNonbond()) {
    mprintf("DEBUG: Topology has nonbond parameters.\n");
    bool hasLJ14 = !NB0.LJ14().empty();
    if (hasLJ14)
      mprintf("DEBUG: Topology has 1-4 nonbond parameters.\n");
    bool hasLJC = !NB0.LJC_Array().empty();
    if (hasLJC)
      mprintf("DEBUG: Topology has LJC nonbond paramters.\n");
    // Nonbonded parameters are present.
    for (std::vector<Atom>::const_iterator atm = atoms.begin(); atm != atoms.end(); ++atm)
    {
      if (!atm->HasType()) {
        mprintf("Warning: Topology has nonbond parameters but atom %s has no type.\n", *(atm->Name()));
        continue;
      }
      TypeNameHolder atype( atm->Type() );
      // Check for self parameters to back-calculate LJ depth/radius
      int idx = NB0.GetLJindex( atm->TypeIndex(), atm->TypeIndex() );
      AtomType thisType;
      if (idx > -1) {
        // Has LJ 6-12 parameters
        NonbondType const& LJ = NB0.NBarray( idx );
        thisType = AtomType(LJ.Radius(), LJ.Depth(), atm->Mass(), atm->Polar());
        if (hasLJ14) {
          NonbondType const& lj14 = NB0.LJ14( idx );
          thisType.SetLJ14( LJparmType(lj14.Radius(), lj14.Depth()) );
        }
        // FIXME do LJ C
      } else {
        // Has LJ 10-12 parameters
        thisType = AtomType(atm->Mass(), atm->Polar());
      }
      thisType.SetTypeIdx( atm->TypeIndex() );
      ParameterHolders::RetType ret = atomTypesOut.AddParm( atype, thisType, true );
      //if (debug_ > 0 && ret == ParameterHolders::ADDED) { // FIXME
      if (ret == ParameterHolders::ADDED) {
        mprintf("DEBUG: New atom type: %s R=%g D=%g M=%g P=%g\n", *(atype[0]),
                thisType.LJ().Radius(), thisType.LJ().Depth(), thisType.Mass(), thisType.Polarizability());
        if (hasLJ14)
          mprintf("DEBUG: New LJ14 atom type: %s R=%g D=%g\n", *(atype[0]), thisType.LJ14().Radius(), thisType.LJ14().Depth());
      }
    }
    // Do atom type pairs, check for off-diagonal elements.
    // Explicitly store pairs instead of regenerating to avoid round-off issues.
    //GetLJterms(atomTypesOut, LJ612out, &LJ1012out, NB0, false, debug_);
    //if (hasLJ14)
    //  GetLJterms(LJ14typesOut, LJ14out, 0, NB0, true, debug_);
    unsigned int nModifiedOffDiagonal = 0;
    unsigned int nModified14OffDiagonal = 0;
    for (ParmHolder<AtomType>::const_iterator i1 = atomTypesOut.begin(); i1 != atomTypesOut.end(); ++i1)
    {
      for (ParmHolder<AtomType>::const_iterator i2 = i1; i2 != atomTypesOut.end(); ++i2)
      {
        NameType const& name1 = i1->first[0];
        NameType const& name2 = i2->first[0];
        TypeNameHolder types(2);
        types.AddName( name1 );
        types.AddName( name2 );
        // Extract original nonbonded parameters for this type pair.
        AtomType const& type1 = i1->second;
        AtomType const& type2 = i2->second;
        int idx1 = type1.OriginalIdx();
        int idx2 = type2.OriginalIdx();
        int idx = NB0.GetLJindex( idx1, idx2 );
        if (idx < 0) {
          // This is LJ 10-12.
          //mprinterr("Error: No off-diagonal LJ for  %s %s (%i %i)\n",
          //          *name1, *name2, idx1, idx2);
          //return;
          mprintf("DEBUG: LJ 10-12 parameters detected for %s %s (%i %i)\n",
                  *name1, *name2, idx1, idx2);
          LJ1012out.AddParm( types, NB0.HBarray((-idx)-1), false );
        } else {
          // This is LJ 6-12.
          // Determine what A and B parameters would be.
          NonbondType lj0 = type1.LJ().Combine_LB( type2.LJ() );
        
          NonbondType lj1 = NB0.NBarray( idx );
          // Compare them
          if (lj0 != lj1) {
            nModifiedOffDiagonal++;
            //if (debug_ > 0) {
              double deltaA = fabs(lj0.A() - lj1.A());
              double deltaB = fabs(lj0.B() - lj1.B());
              mprintf("DEBUG: Potential off-diagonal LJ: %s %s expect A=%g B=%g, actual A=%g B=%g\n",
                      *name1, *name2, lj0.A(), lj0.B(), lj1.A(), lj1.B());
              mprintf("DEBUG:\tdeltaA= %g    deltaB= %g\n", deltaA, deltaB);
              double pe_a = (fabs(lj0.A() - lj1.A()) / lj0.A());
              double pe_b = (fabs(lj0.B() - lj1.B()) / lj0.B());
              mprintf("DEBUG:\tPEA= %g  PEB= %g\n", pe_a, pe_b);
            //}
          }
          LJ612out.AddParm( types, lj1, false );
          if (hasLJ14) {
            // This is LJ 6-12 1-4.
            // Determine what A and B parameters would be.
            lj0 = type1.LJ14().Combine_LB( type2.LJ14() );
        
            lj1 = NB0.LJ14( idx );
            // Compare them
            if (lj0 != lj1) {
              nModified14OffDiagonal++;
              //if (debug_ > 0) {
                double deltaA = fabs(lj0.A() - lj1.A());
                double deltaB = fabs(lj0.B() - lj1.B());
                mprintf("DEBUG: Potential off-diagonal LJ 1-4: %s %s expect A=%g B=%g, actual A=%g B=%g\n",
                        *name1, *name2, lj0.A(), lj0.B(), lj1.A(), lj1.B());
                mprintf("DEBUG:\tdeltaA= %g    deltaB= %g\n", deltaA, deltaB);
                double pe_a = (fabs(lj0.A() - lj1.A()) / lj0.A());
                double pe_b = (fabs(lj0.B() - lj1.B()) / lj0.B());
                mprintf("DEBUG:\tPEA= %g  PEB= %g\n", pe_a, pe_b);
              //}
            }
            LJ14out.AddParm( types, lj1, false );
          } // END hasLJ14
          if (hasLJC) {
            // This is LJC
            double ljc1 = NB0.LJC_Array( idx );
            LJCout.AddParm( types, ljc1, false );
          }
        }
      } // END inner loop over atom types
    } // END outer loop over atom types
    if (nModifiedOffDiagonal > 0)
      mprintf("Warning: %u modified off-diagonal LJ terms present.\n", nModifiedOffDiagonal);
    if (nModified14OffDiagonal > 0)
      mprintf("Warning: %u modified off-diagonal LJ 1-4 terms present.\n", nModified14OffDiagonal);
  } else {
    if (!atoms.empty())
      mprintf("DEBUG: Topology does not have nonbond parameters.\n");
    // No nonbonded parameters. Just save mass/polarizability.
    for (std::vector<Atom>::const_iterator atm = atoms.begin(); atm != atoms.end(); ++atm)
      if (atm->HasType() > 0)
        atomTypesOut.AddParm( TypeNameHolder(atm->Type()), AtomType(atm->Mass(), atm->Polar()), true );
  }
}

