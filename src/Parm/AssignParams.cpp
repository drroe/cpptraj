#include "AssignParams.h"
#include "../Atom.h"
#include "../AtomType.h"
#include "../CpptrajStdio.h"
#include "../ParameterHolders.h"
#include "../Topology.h"
#include "../TypeNameHolder.h"

using namespace Cpptraj::Parm;

/** CONSTRUCTOR */
AssignParams::AssignParams() :
  debug_(0)
{}

/** Set parameters for atoms via given atom type parameter holder. */
void AssignParams::AssignAtomTypeParm(AtArray& atoms, ParmHolder<AtomType> const& newAtomTypeParams)
const
{
  for (AtArray::iterator iat = atoms.begin(); iat != atoms.end(); ++iat)
  {
    bool found;
    TypeNameHolder atype( iat->Type() );
    AtomType AT = newAtomTypeParams.FindParam( atype, found );
    if (found) {
      // Update mass
      iat->SetMass( AT.Mass() );
      // Update polarizability
      iat->SetPolar( AT.Polarizability() );
    } else
      mprintf("Warning: Atom type parameter not found for %s\n", *atype[0]);
  }
}

/** Set parameters for bonds in given bond array. */
void AssignParams::AssignBondParm(Topology const& topOut,
                                  ParmHolder<BondParmType> const& newBondParams,
                                  BondArray& bonds, BondParmArray& bpa, const char* desc)
const
{
  ParmHolder<int> currentTypes;
  for (BondArray::iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd) {
    TypeNameHolder types(2);
    types.AddName( topOut[bnd->A1()].Type() );
    types.AddName( topOut[bnd->A2()].Type() );
    bool found;
    // See if parameter is present.
    int idx = currentTypes.FindParam(types, found);
    if (!found) {
      idx = -1;
      // Not yet present in current types.
      BondParmType bp = newBondParams.FindParam( types, found );
      if (!found) {
        mprintf("Warning: parameter not found for %s %s-%s (%s-%s)\n", desc,
                topOut.TruncResAtomNameNum(bnd->A1()).c_str(),
                topOut.TruncResAtomNameNum(bnd->A2()).c_str(),
                *types[0], *types[1]);
      } else {
        //idx = addBondParm( bpa, bp ); TODO handle array packing
        idx = (int)bpa.size();
        bpa.push_back( bp );
        currentTypes.AddParm(types, idx, false);
      }
    }
    bnd->SetIdx( idx );
  }
/*
  for (BondArray::iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd) {
    TypeNameHolder types(2);
    types.AddName( atoms_[bnd->A1()].Type() );
    types.AddName( atoms_[bnd->A2()].Type() );
    bool found;
    // See if parameter is present.
    int idx = -1;
    BondParmType bp = newBondParams.FindParam( types, found );
    if (!found) {
      mprintf("Warning: parameter not found for %s %s-%s (%s-%s)\n", desc,
              TruncResAtomNameNum(bnd->A1()).c_str(),
              TruncResAtomNameNum(bnd->A2()).c_str(),
              *types[0], *types[1]);
    } else {
      idx = addBondParm( bpa, bp );
    }
    bnd->SetIdx( idx );
  }*/
}

/** Replace any current bond parameters with given bond parameters. */
void AssignParams::AssignBondParams(Topology& topOut, ParmHolder<BondParmType> const& newBondParams) const {
  topOut.ModifyBondParm().clear();
  AssignBondParm( topOut, newBondParams, topOut.ModifyBonds(),  topOut.ModifyBondParm(), "bond" );
  AssignBondParm( topOut, newBondParams, topOut.ModifyBondsH(), topOut.ModifyBondParm(), "bond" );
}

/** Replace any current Urey-Bradley parameters with given UB parameters. */
void AssignParams::AssignUBParams(Topology& topOut, ParmHolder<BondParmType> const& newBondParams) const {
  topOut.ModifyUBparm().clear();
  AssignBondParm( topOut, newBondParams, topOut.ModifyUB(), topOut.ModifyUBparm(), "UB term" );
}

/** Set parameters for angles in given angle array. */
AngleArray AssignParams::AssignAngleParm(Topology const& topOut,
                                         ParmHolder<AngleParmType> const& newAngleParams,
                                         AngleArray const& angles, AngleParmArray& apa)
const
{
  AngleArray newAngles;
  ParmHolder<int> currentTypes;
  for (AngleArray::const_iterator ang = angles.begin(); ang != angles.end(); ++ang) {
    TypeNameHolder types(3);
    types.AddName( topOut[ang->A1()].Type() );
    types.AddName( topOut[ang->A2()].Type() );
    types.AddName( topOut[ang->A3()].Type() );
    // Skip extra points // FIXME make this an option
    if ( topOut[ang->A1()].Element() == Atom::EXTRAPT ||
         topOut[ang->A3()].Element() == Atom::EXTRAPT)
    {
      mprintf("DEBUG: Skipping angle with extra point: %4i %4i %4i (%2s %2s %2s)\n",
              ang->A1()+1, ang->A2()+1, ang->A3()+1,
              *types[0], *types[1], *types[2]);
      continue;
    }
    // Skip water angles // FIXME make this an option
    if (topOut[ang->A1()].Element() == Atom::HYDROGEN &&
        topOut[ang->A2()].Element() == Atom::OXYGEN &&
        topOut[ang->A3()].Element() == Atom::HYDROGEN)
    {
      // H-O-H Angle. If there is an H-H bond assume this is a rigid water model.
      if (topOut[ang->A1()].IsBondedTo( ang->A3() )) {
        mprintf("DEBUG: H-O-H angle, H-H bond detected. Assuming rigid water, skipping angle: "
                "%4i %4i %4i (%2s %2s %2s)\n",
                ang->A1()+1, ang->A2()+1, ang->A3()+1,
                *types[0], *types[1], *types[2]);
        continue;
      }
    }
    if (topOut[ang->A1()].Element() == Atom::OXYGEN &&
        topOut[ang->A2()].Element() == Atom::HYDROGEN &&
        topOut[ang->A3()].Element() == Atom::HYDROGEN)
    {
      // O-H-H Angle. Assume rigid water model
      mprintf("DEBUG: O-H-H angle. Assuming rigid water, skipping angle: "
              "%4i %4i %4i (%2s %2s %2s)\n",
              ang->A1()+1, ang->A2()+1, ang->A3()+1,
              *types[0], *types[1], *types[2]);
      continue;
    }

    bool found;
    // See if parameter is present.
    int idx = currentTypes.FindParam( types, found );
    if (!found) {
      idx = -1;
      // Not yet present in current types
      AngleParmType ap = newAngleParams.FindParam( types, found );
      if (!found) {
        mprintf("Warning: Angle parameter not found for angle %s-%s-%s (%s-%s-%s)\n",
                topOut.TruncResAtomNameNum(ang->A1()).c_str(),
                topOut.TruncResAtomNameNum(ang->A2()).c_str(),
                topOut.TruncResAtomNameNum(ang->A3()).c_str(),
                *types[0], *types[1], *types[2]);
      } else {
        //idx = addAngleParm( angleparm_, ap ); // TODO uncomment for array packing
        idx = (int)apa.size();
        apa.push_back( ap );
        //mprintf("DEBUG: New angle type for %s-%s-%s (Tk=%g Teq=%g) idx=%i\n", *types[0], *types[1], *types[2], ap.Tk(), ap.Teq(), idx);
        currentTypes.AddParm(types, idx, false);
      }
    }
    newAngles.push_back( *ang );
    newAngles.back().SetIdx( idx );
  }
  return newAngles;
/*
  for (AngleArray::iterator ang = angles.begin(); ang != angles.end(); ++ang) {
    TypeNameHolder types(3);
    types.AddName( atoms_[ang->A1()].Type() );
    types.AddName( atoms_[ang->A2()].Type() );
    types.AddName( atoms_[ang->A3()].Type() );
    bool found;
    // See if parameter is present.
    int idx = -1;
    AngleParmType ap = newAngleParams.FindParam( types, found );
    if (!found) {
      mprintf("Warning: Angle parameter not found for angle %s-%s-%s (%s-%s-%s)\n",
              TruncResAtomNameNum(ang->A1()).c_str(),
              TruncResAtomNameNum(ang->A2()).c_str(),
              TruncResAtomNameNum(ang->A3()).c_str(),
              *types[0], *types[1], *types[3]);
    } else {
      idx = addAngleParm( angleparm_, ap );
    }
    ang->SetIdx( idx );
  }*/
}

/** Replace any current angle parameters with given angle parameters. */
void AssignParams::AssignAngleParams(Topology& topOut, ParmHolder<AngleParmType> const& newAngleParams) const {
  topOut.ModifyAngleParm().clear();
  topOut.ModifyAngles() = AssignAngleParm( topOut, newAngleParams, topOut.Angles(), topOut.ModifyAngleParm() );
  topOut.ModifyAnglesH() = AssignAngleParm( topOut, newAngleParams, topOut.AnglesH(), topOut.ModifyAngleParm() );
}

/** Warn if improper atoms have been reordered so they match the parameter. */
void AssignParams::warn_improper_reorder(DihedralType const& imp0, DihedralType const& imp)
const
{
  if (debug_ < 1) return;
  mprintf("Warning: Improper types have been reordered from %4s %4s %4s %4s",
          *(atoms_[imp0.A1()].Type()),
          *(atoms_[imp0.A2()].Type()),
          *(atoms_[imp0.A3()].Type()),
          *(atoms_[imp0.A4()].Type()));
  mprintf(" to %4s %4s %4s %4s to match improper parameter.\n",
          *(atoms_[imp.A1()].Type()),
          *(atoms_[imp.A2()].Type()),
          *(atoms_[imp.A3()].Type()),
          *(atoms_[imp.A4()].Type()));
}

/** Set parameters for improper dihedrals in given improper dihedral array. */
void AssignParams::AssignImproperParm(ImproperParmHolder const& newImproperParams,
                                  DihedralArray& impropers,
                                  DihedralParmArray& dpa)
const
{
  for (DihedralArray::iterator imp = impropers.begin(); imp != impropers.end(); ++imp) {
    TypeNameHolder types(4);
    types.AddName( atoms_[imp->A1()].Type() );
    types.AddName( atoms_[imp->A2()].Type() );
    types.AddName( atoms_[imp->A3()].Type() );
    types.AddName( atoms_[imp->A4()].Type() );
    bool found;
    // See if parameter is present.
    int idx = -1;
    DihedralType imp0 = *imp;
    bool reordered;
    DihedralParmArray ipa = newImproperParams.FindParam( types, found, *imp, reordered );
    if (!found) {
      mprintf("Warning: Parameter not found for improper %s-%s-%s-%s (%s-%s-%s-%s)\n",
              TruncResAtomNameNum(imp->A1()).c_str(),
              TruncResAtomNameNum(imp->A2()).c_str(),
              TruncResAtomNameNum(imp->A3()).c_str(),
              TruncResAtomNameNum(imp->A4()).c_str(),
              *types[0], *types[1], *types[3], *types[4]);
    } else {
      if (ipa.size() > 1)
        mprintf("Warning: %zu improper parameters found for types %s - %s - %s - %s, expected only one."
                        "Warning: Only using first parameter.\n", ipa.size(), *(types[0]), *(types[1]), *(types[2]), *(types[3]));
      if (reordered) warn_improper_reorder( imp0, *imp );
      idx = addTorsionParm( dpa, ipa.front() );
    }
    imp->SetIdx( idx );
  }
}

/** Replace any current improper parameters with given improper parameters. */
void AssignParams::AssignImproperParams(ImproperParmHolder const& newImproperParams) const {
  improperparm_.clear();
  AssignImproperParm( newImproperParams, impropers_, improperparm_ );
}

/** Set parameters for dihedrals in given dihedral array.
  * Bond and angle information must be set up prior to calling
  * this function in order for improper and 1-4 detection to work.
  * \param newDihedralParams New proper dihedral parameters.
  * \param newImproperParams New improper dihedral parameters.
  * \param dihedrals Array containing only unique dihedrals.
  * \param sort_improper_cache If true, sort improper types in cache (to match current leap behavior)
  */
DihedralArray AssignParams::AssignDihedralParm(DihedralParmHolder const& newDihedralParams,
                                           ImproperParmHolder const& newImproperParams,
                                           ParmHolder<AtomType> const& AT,
                                           DihedralArray const& dihedrals,
                                           bool sort_improper_cache)
const
{ // TODO skip extra points
  DihedralArray dihedralsIn;
  // Improper cache
  ImproperParmHolder improperCache;
  improperCache.SetWildcard( newImproperParams.Wildcard() );
  improperCache.SetRequireExactMatch( newImproperParams.RequireExactMatch() );
  // Keep track of 1-4 interactions
  typedef std::pair<int,int> Ipair;
  typedef std::set<Ipair> Imap;
  Imap PairMap;
  // Loop over all dihedrals
  for (DihedralArray::const_iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih) {
    TypeNameHolder types(4);
    types.AddName( atoms_[dih->A1()].Type() );
    types.AddName( atoms_[dih->A2()].Type() );
    types.AddName( atoms_[dih->A3()].Type() );
    types.AddName( atoms_[dih->A4()].Type() );
    // Skip extra points // FIXME make this an option
    if ( atoms_[dih->A1()].Element() == Atom::EXTRAPT ||
         atoms_[dih->A4()].Element() == Atom::EXTRAPT)
    {
      mprintf("DEBUG: Skipping dihedral with extra point: %4i %4i %4i %4i (%2s %2s %2s %2s)\n",
              dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1,
              *types[0], *types[1], *types[2], *types[3]);
      continue;
    }
//    mprintf("DEBUG: Assigning dihedral %4i %4i %4i %4i (%2s %2s %2s %2s) isImproper=%i skip14=%i\n",
//            dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1,
//            *types[0], *types[1], *types[2], *types[3],
//            (int)dih->IsImproper(), (int)dih->Skip14());
    // Determine improper
    bool isImproper = (!atoms_[dih->A1()].IsBondedTo(dih->A2())) ||
                      (!atoms_[dih->A2()].IsBondedTo(dih->A3())) ||
                      (!atoms_[dih->A3()].IsBondedTo(dih->A4()));
    if (isImproper != dih->IsImproper()) {
      mprintf("Warning: dihedral %s-%s-%s-%s improper status %i does not match detected (%i)\n",
              TruncResAtomNameNum(dih->A1()).c_str(),
              TruncResAtomNameNum(dih->A2()).c_str(),
              TruncResAtomNameNum(dih->A3()).c_str(),
              TruncResAtomNameNum(dih->A4()).c_str(),
              (int)dih->IsImproper(), (int)isImproper);
    }
    bool found;
    if (dih->IsImproper()) {
      // ----- This is actually an improper dihedral. ----------------
      // Impropers are treated differently than other topology types. If
      // no parameter is found, do not add it to the list of dihedrals.
      // However, if no parameter is found and the central atom is
      // SP2, print a warning.
      found = false;
      DihedralType mydih = *dih;
      bool reordered;
      bool is_new_improper = false;
      DihedralParmArray ipa;
      TypeNameHolder paramTypes;
      ImproperParmHolder::const_iterator impit = improperCache.GetParam( types, mydih, reordered );
      if (impit == improperCache.end()) {
        impit = newImproperParams.GetParam( types, mydih, reordered );
        if (impit != newImproperParams.end()) {
          is_new_improper = true;
          paramTypes = impit->first;
          ipa = impit->second;
          found = true;
          if (debug_ > 1)
            mprintf("DEBUG: Found new value for improper %2s %2s %2s %2s (%2s %2s %2s %2s)\n",
                    *types[0], *types[1], *types[2], *types[3],
                    *paramTypes[0], *paramTypes[1], *paramTypes[2], *paramTypes[3]);
        }
      } else {
        // If the cached improper has wildcards, see if there is a more
        // specific one in the new improper parameters.
        int n_wildcards = 0;
        if (newImproperParams.Wildcard().len() > 0) {
          for (TypeNameHolder::const_iterator tt = impit->first.begin(); tt != impit->first.end(); ++tt) {
            if ( *tt == newImproperParams.Wildcard() ) {
              n_wildcards++;
            }
          }
        }
        if (n_wildcards > 0) {
          // Check if there is a more specific improper.
          bool newReordered;
          DihedralType newdih = *dih;
          ImproperParmHolder::const_iterator moreSpecificImpit = newImproperParams.GetParam( types, newdih, newReordered );
          if (moreSpecificImpit != newImproperParams.end()) {
            // See if the new improper has fewer wildcards than the cached one.
            int new_wildcards = 0;
            for (TypeNameHolder::const_iterator tt = moreSpecificImpit->first.begin();
                                                tt != moreSpecificImpit->first.end(); ++tt)
            {
              if ( *tt == newImproperParams.Wildcard() ) {
                new_wildcards++;
              }
            }
            if (new_wildcards < n_wildcards) {
              // More specific improper was found.
              mprintf("DEBUG: A more specific improper was found for %2s %2s %2s %2s (%2s %2s %2s %2s)\n",
                      *(impit->first[0]), *(impit->first[1]), *(impit->first[2]), *(impit->first[3]),
                      *(moreSpecificImpit->first[0]), *(moreSpecificImpit->first[1]), *(moreSpecificImpit->first[2]), *(moreSpecificImpit->first[3]));
              is_new_improper = true;
              reordered = newReordered;
              mydih = newdih;
              impit = moreSpecificImpit;
            }
          }
        } // END found improper had wildcards
        paramTypes = impit->first;
        ipa = impit->second;
        found = true;
        if (debug_ > 1)
          mprintf("DEBUG: Using cached value for improper %2s %2s %2s %2s (%2s %2s %2s %2s)\n",
                  *types[0], *types[1], *types[2], *types[3],
                  *paramTypes[0], *paramTypes[1], *paramTypes[2], *paramTypes[3]);
      }
      int idx = -1;
      if (!found) {
        if (debug_ > 0)
          mprintf("Warning: Improper parameters not found for improper dihedral %s-%s-%s-%s (%s-%s-%s-%s)\n",
                  TruncResAtomNameNum(dih->A1()).c_str(),
                  TruncResAtomNameNum(dih->A2()).c_str(),
                  TruncResAtomNameNum(dih->A3()).c_str(),
                  TruncResAtomNameNum(dih->A4()).c_str(),
                  *types[0], *types[1], *types[2], *types[3]);
        // Central atom
        Atom const& AJ = atoms_[dih->A3()];
        AtomType::HybridizationType hybrid = AtomType::UNKNOWN_HYBRIDIZATION;
        AtomType atype = AT.FindParam(TypeNameHolder(AJ.Type()), found);
        if (found)
          hybrid = atype.Hybridization();
        if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION) {
          mprintf("Warning: Guessing hybridization for improper central atom %s\n", AtomMaskName(dih->A3()).c_str());
          hybrid = Cpptraj::GuessAtomHybridization( AJ, atoms_ );
        }
        if (hybrid == AtomType::SP2) {
          mprintf("Warning: No improper parameters for SP2 hybridized atom %s\n", AtomMaskName(dih->A3()).c_str());
        }
      } else {
        if (ipa.size() > 1)
          mprintf("Warning: %zu improper parameters found for types %s - %s - %s - %s, expected only one."
                  "Warning: Only using first parameter.\n", ipa.size(), *(types[0]), *(types[1]), *(types[2]), *(types[3]));
        if (reordered) warn_improper_reorder( *dih, mydih );
        idx = addTorsionParm( dihedralparm_, ipa.front() );
        mydih.SetIdx( idx );
        mydih.SetImproper( true );
        // Always skip 1-4 for impropers
        mydih.SetSkip14( true );
        dihedralsIn.push_back( mydih );
        // Add to the cache
        if (is_new_improper) {
          // To match leap behavior, make sure paramTypes are sorted alphabetically.
          //mprintf("DEBUG: Improper wildcard: %s\n", *(newImproperParams.Wildcard()));
          if (sort_improper_cache) paramTypes.SortImproperByAlpha( newImproperParams.Wildcard() );
          improperCache.AddParm( paramTypes, ipa.front(), false );
        }
      }
    } else {
      // -----Regular dihedral. See if parameter already present. ----
      DihedralParmArray dpa = newDihedralParams.FindParam( types, found );
      if (!found) {
        mprintf("Warning: Dihedral parameters not found for dihedral %s-%s-%s-%s (%s-%s-%s-%s)\n",
                TruncResAtomNameNum(dih->A1()).c_str(),
                TruncResAtomNameNum(dih->A2()).c_str(),
                TruncResAtomNameNum(dih->A3()).c_str(),
                TruncResAtomNameNum(dih->A4()).c_str(),
                *types[0], *types[1], *types[2], *types[3]);
        DihedralType mydih = *dih;
        mydih.SetIdx( -1 );
        dihedralsIn.push_back( mydih );
      } else {
        // Actually add parameters for this dihedral.
        // Determine if this is actually a 1-4 interaction by making
        // sure that A4 isnt part of a bond or angle with A1.
        bool skip14 = false;
        for (Atom::bond_iterator bat1 = atoms_[dih->A1()].bondbegin();
                                 bat1 != atoms_[dih->A1()].bondend(); ++bat1)
        {
          if (*bat1 != dih->A2()) {
            if (*bat1 == dih->A4()) {
              skip14 = true;
              break;
            }
            // Loop over angles, dih->A1() - bat1 - bat2
            for (Atom::bond_iterator bat2 = atoms_[*bat1].bondbegin();
                                     bat2 != atoms_[*bat1].bondend(); ++bat2)
            {
              //if (dih->A1() == 442 && dih->A4() == 444) { // DEBUG
              //  mprintf("DEBUG: %4i %4i %4i %4i Checking angle %4i %4i %4i\n", dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1, dih->A1()+1, *bat1 + 1, *bat2 + 1);
              //}
              if (*bat2 != *bat1) {
                if (*bat2 == dih->A4()) {
                  skip14 = true;
                  break;
                }
              }
            } // END loop over bat1 bonded atoms
            if (skip14) break;
          }
        } // END loop over dih->A1() bonded atoms
        if (!skip14) {
          // Determine if 1-4 interaction already calculated by previous dihedral.
          // Make the lower end atom first. 
          Ipair pair14;
          if (dih->A1() < dih->A4()) {
            pair14.first = dih->A1();
            pair14.second = dih->A4();
          } else {
            pair14.first = dih->A4();
            pair14.second = dih->A1();
          }
          Imap::const_iterator it = PairMap.find( pair14 );
          if (it == PairMap.end()) {
            skip14 = false;
            PairMap.insert( pair14 );
          } else {
            skip14 = true;
            //mprintf("DEBUG: Prior 1-4 calc detected.\n");
          }
        }
        // Loop over multiplicities
        //if (dih->A1() == 442 && dih->A4() == 444) { // DEBUG
        //  mprintf("DEBUG: Skip14= %i\n", (int)skip14);
        //}
        for (DihedralParmArray::const_iterator it = dpa.begin(); it != dpa.end(); ++it) {
          DihedralType mydih = *dih;
          int idx = addTorsionParm( dihedralparm_, *it );
          // If there are multiple parameters for the same dihedral, all but
          // one of the 1-4 calcs need to be skipped.
          if (it == dpa.begin())
            mydih.SetSkip14(skip14);
          else
            mydih.SetSkip14(true);
          mydih.SetIdx( idx );
          dihedralsIn.push_back( mydih );
//          mprintf("DEBUG: Assigned %4i %4i %4i %4i (%2s %2s %2s %2s) isImproper=%i skip14=%i\n",
//            mydih.A1()+1, mydih.A2()+1, mydih.A3()+1, mydih.A4()+1,
//            *types[0], *types[1], *types[2], *types[3],
//            (int)mydih.IsImproper(), (int)mydih.Skip14());
        }
      }
    }
  } // END loop over all dihedrals
  return dihedralsIn;
}

/** \return Array of only unique dihedrals (get rid of multiplicities) */
DihedralArray AssignParams::get_unique_dihedrals(DihedralArray const& dihedralsIn) const {
  // Assume dihedrals with multiple terms are consecutive already
  DihedralArray dihedrals;
  for (DihedralArray::const_iterator dih = dihedralsIn.begin(); dih != dihedralsIn.end(); ++dih) {
    if (dihedrals.empty())
      dihedrals.push_back( *dih );
    else {
      if ( *dih != dihedrals.back() ) // TODO skip end dihedrals, impropers?
        dihedrals.push_back( *dih );
    }
  }
  if (debug_ > 0)
    mprintf("DEBUG: %zu incoming dihedrals, %zu unique dihedrals.\n",
            dihedralsIn.size(), dihedrals.size());
  return dihedrals;
}

/** Replace any current dihedral parameters with given dihedral parameters. */
void AssignParams::AssignDihedralParams(DihedralParmHolder const& newDihedralParams,
                                        ImproperParmHolder const& newImproperParams,
                                        ParmHolder<AtomType> const& AT)
const
{
  dihedralparm_.clear();
  // Dihedrals can be a bit of a pain since there can be multiple
  // multiplicities for a single dihedral type. In case multiplicities
  // change, start with a fresh dihedral array containing only unique
  // dihedrals.
  dihedrals_  = AssignDihedralParm( newDihedralParams, newImproperParams, AT, get_unique_dihedrals(dihedrals_), false  );
  dihedralsh_ = AssignDihedralParm( newDihedralParams, newImproperParams, AT, get_unique_dihedrals(dihedralsh_), false );
}

/** Replace current nonbond parameters with given nonbond parameters. */
//TODO Accept array of atom type names?
void AssignParams::AssignNonbondParams(ParmHolder<AtomType> const& newTypes,
                                   ParmHolder<NonbondType> const& newNB,
                                   ParmHolder<NonbondType> const& new14,
                                   ParmHolder<double> const& newLJC,
                                   ParmHolder<HB_ParmType> const& newHB,
                                   int verbose)
const
{
  bool hasLJ14 = !new14.empty();
  bool hasLJC = !newLJC.empty();
  // Generate array of only the types that are currently in Topology. TODO should this be a permanent part of Topology?
  ParmHolder<AtomType> currentAtomTypes;
  for (AtArray::const_iterator atm = atoms_.begin(); atm != atoms_.end(); ++atm)
  {
    if (atm->HasType()) {
      TypeNameHolder types(1);
      types.AddName(atm->Type());
      // Find in newTypes.
      bool found;
      AtomType newAT = newTypes.FindParam( types, found );
      if (!found) {
        mprintf("Warning: No atom type information for type '%s'\n", *types[0]);
        newAT = AtomType(atm->Mass());
      }
      // Check LJ14 types if needed
      if (hasLJ14) {
        if (!newAT.HasLJ14()) {
          mprintf("Warning: No 1-4 atom type information for type '%s'\n", *types[0]); // TODO error out here?
        }
      }
      currentAtomTypes.AddParm( types, newAT, true );
    }
  }
  if (currentAtomTypes.size() < 1) {
    mprintf("Warning: No atom type information in %s - cannot assign nonbond parameters.\n",
            c_str());
    return;
  }
  // Regenerate nonbond params for existing types
  nonbond_.Clear();
  // Set type indices in order.
  for (ParmHolder<AtomType>::iterator t1 = currentAtomTypes.begin(); t1 != currentAtomTypes.end(); ++t1)
    t1->second.SetTypeIdx( -1 );
  int n_unique_lj_types = 0;
  for (ParmHolder<AtomType>::iterator t1 = currentAtomTypes.begin(); t1 != currentAtomTypes.end(); ++t1)
  {
    if (t1->second.OriginalIdx() == -1) {
      t1->second.SetTypeIdx( n_unique_lj_types );
      // Look for equivalent nonbond types
      for (ParmHolder<AtomType>::iterator t2 = t1 + 1; t2 != currentAtomTypes.end(); ++t2) {
        // Everything must match
        if (!hasLJ14) {
          if (t2->second.OriginalIdx() == -1 && t1->second.LJ() == t2->second.LJ()) {
            mprintf("DEBUG: Type %s equivalent to type %s\n", *(t2->first[0]), *(t1->first[0]));
            t2->second.SetTypeIdx( n_unique_lj_types );
          }
        } else {
          if (t2->second.OriginalIdx() == -1 &&
              t1->second.LJ() == t2->second.LJ() &&
              t1->second.LJ14() == t2->second.LJ14())
          {
            mprintf("DEBUG: Type %s (with 1-4 params) equivalent to type %s\n", *(t2->first[0]), *(t1->first[0]));
            t2->second.SetTypeIdx( n_unique_lj_types );
          }
        }
      }
      n_unique_lj_types++;
    }
  }
  mprintf("DEBUG: Setting up nonbond array for %i unique LJ types.\n", n_unique_lj_types);
  nonbond_.SetupLJforNtypes( n_unique_lj_types );
  if (hasLJ14)
    nonbond_.SetNLJ14terms( nonbond_.NBarray().size() );
  if (hasLJC)
    nonbond_.SetNLJCterms( nonbond_.NBarray().size() );
  // Loop over all atom type pairs
  for (ParmHolder<AtomType>::const_iterator t1 = currentAtomTypes.begin(); t1 != currentAtomTypes.end(); ++t1)
  {
    NameType const& name1 = t1->first[0];
    //mprintf("DEBUG: Type1= %s (%i)\n", *name1, nidx1);
    AtomType const& type1 = t1->second;
    for (ParmHolder<AtomType>::const_iterator t2 = t1; t2 != currentAtomTypes.end(); ++t2)
    {
      NameType const& name2 = t2->first[0];
      //mprintf("DEBUG:\t\tType2= %s (%i)\n", *name2, nidx2);
      AtomType const& type2 = t2->second;
      TypeNameHolder types(2);
      types.AddName( name1 );
      types.AddName( name2 );
      int newnbidx;
      // Check for LJ10-12 first
      ParmHolder<HB_ParmType>::const_iterator hb = newHB.GetParam( types );
      if (hb != newHB.end()) {
        if (verbose > 0) mprintf("LJ 10-12 parameter found for %s %s\n", *name1, *name2);
        nonbond_.AddHBterm(t1->second.OriginalIdx(), t2->second.OriginalIdx(), hb->second);
        // Also add a blank 6-12 term since that seems to be the convention
        newnbidx = nonbond_.AddLJterm( t1->second.OriginalIdx(), t2->second.OriginalIdx(), NonbondType() );
      } else {
        // See if this parameter exists in the given nonbond array.
        NonbondType LJAB;
        ParmHolder<NonbondType>::const_iterator it = newNB.GetParam( types );
        if (it == newNB.end()) {
          if (verbose > 0) mprintf("NB parameter for %s %s not found. Generating.\n", *name1, *name2);
          LJAB = type1.LJ().Combine_LB( type2.LJ() );
        } else {
          if (verbose > 0) mprintf("Using existing NB parameter for %s %s\n", *name1, *name2);
          LJAB = it->second;
        }
        newnbidx = nonbond_.AddLJterm(t1->second.OriginalIdx(), t2->second.OriginalIdx(), LJAB);
      }
      // LJ 1-4
      if (hasLJ14) {
        // Get parameter if it exists
        ParmHolder<NonbondType>::const_iterator it = new14.GetParam( types );
        if (it == new14.end()) {
          if (verbose > 0) mprintf("NB 1-4 parameter for %s %s not found. Generating.\n", *name1, *name2);
          nonbond_.SetLJ14( newnbidx ) = type1.LJ14().Combine_LB( type2.LJ14() );
        } else {
          if (verbose > 0) mprintf("Using existing NB 1-4 parameter for %s %s\n", *name1, *name2);
          nonbond_.SetLJ14( newnbidx ) = it->second;
        }
      } // END hasLJ14
      // LJC
      if (hasLJC) {
        // Get parameter if it exists
        ParmHolder<double>::const_iterator it = newLJC.GetParam( types );
        if (it == newLJC.end()) {
          mprintf("LJC parameter for %s %s not found. Setting to 0.\n", *name1, *name2); // FIXME generate
          nonbond_.SetLJC( newnbidx, 0 );
        } else {
          mprintf("Using existing LJC parameter for %s %s\n", *name1, *name2);
          nonbond_.SetLJC( newnbidx, it->second );
        }
      } // END hasLJC
    } // END inner loop over current types
  } // END outer loop over current types
  // Reset the atom type indices.
  for (AtArray::iterator atm = atoms_.begin(); atm != atoms_.end(); ++atm)
  {
    int tidx = -1;
    ParmHolder<AtomType>::const_iterator it = currentAtomTypes.GetParam( TypeNameHolder(atm->Type()) );
    if (it == currentAtomTypes.end()) {
      mprintf("Warning: Atom type not found for %s (type %s)\n",
              TruncResAtomNameNum( atm - atoms_.begin() ).c_str(),
              *(atm->Type()));
    } else {
      if (it->second.OriginalIdx() < 0 || it->second.OriginalIdx() >= (int)currentAtomTypes.size()) {
        mprinterr("Internal Error: Type index for %s (type %s) out of range: %i\n",
                  TruncResAtomNameNum( atm - atoms_.begin() ).c_str(),
                  *(atm->Type()), it->second.OriginalIdx());
      } else
        tidx = it->second.OriginalIdx();
    }
    atm->SetTypeIndex( tidx );
  }
}

/// \return index of 5th cmap atom if atom names match those in dihedral, -1 otherwise
static inline int cmap_anames_match(DihedralType const& dih,
                                     AtArray const& atoms,
                                     std::vector<std::string> const& cmapAtomNames)
{
/*  mprintf("DEBUG: Check res %i %s %s %s %s against %s %s %s %s\n",
          atoms[dih.A2()].ResNum()+1,
            atoms[dih.A1()].Name().Truncated().c_str(),
            atoms[dih.A2()].Name().Truncated().c_str(),
            atoms[dih.A3()].Name().Truncated().c_str(),
            atoms[dih.A4()].Name().Truncated().c_str(),
            cmap.AtomNames()[0].c_str(),
            cmap.AtomNames()[1].c_str(),
            cmap.AtomNames()[2].c_str(),
            cmap.AtomNames()[3].c_str());*/
  // TODO without Truncated
  if (atoms[dih.A1()].Name().Truncated() == cmapAtomNames[0] &&
      atoms[dih.A2()].Name().Truncated() == cmapAtomNames[1] &&
      atoms[dih.A3()].Name().Truncated() == cmapAtomNames[2] &&
      atoms[dih.A4()].Name().Truncated() == cmapAtomNames[3])
  {
    /*mprintf("DEBUG: Partial match %s %s %s %s = %s %s %s %s\n",
            atoms[dih.A1()].Name().Truncated().c_str(),
            atoms[dih.A2()].Name().Truncated().c_str(),
            atoms[dih.A3()].Name().Truncated().c_str(),
            atoms[dih.A4()].Name().Truncated().c_str(),
            cmap.AtomNames()[0].c_str(),
            cmap.AtomNames()[1].c_str(),
            cmap.AtomNames()[2].c_str(),
            cmap.AtomNames()[3].c_str());*/
    for (Atom::bond_iterator bat = atoms[dih.A4()].bondbegin();
                             bat != atoms[dih.A4()].bondend(); ++bat)
    {
      if (*bat != dih.A3()) { // TODO can A4 ever be bonded to A1??
        if (atoms[*bat].Name().Truncated() == cmapAtomNames[4]) {
          return *bat;
        }
      }
    }
  }
  return -1;
}

/// \return true if given name matches a residue name
static inline bool MatchesResName(std::vector<std::string> const& resNames, std::string const& nameIn) {
  for (std::vector<std::string>::const_iterator it = resNames.begin(); it != resNames.end(); ++it)
    if (nameIn == *it) return true;
  return false;
}

/// Remap cmap indices
int AssignParams::remap_cmap_indices(std::vector<int>& originalCmapIndices,
                                     CmapGridArray& cmapGrids,
                                     CmapArray& cmapTerms,
                                 CmapParmHolder const& cmapIn)
const
{
  mprintf("DEBUG: Cmap parameter indices:\n");
  for (unsigned int idx = 0; idx < originalCmapIndices.size(); idx++)
    mprintf("DEBUG:\t\tCurrent idx=%i  Actual idx=%u\n", originalCmapIndices[idx], idx);
  if (originalCmapIndices.empty()) {
    mprintf("DEBUG: No CMAP indices in %s\n", c_str());
    return 0;
  }
  std::sort(originalCmapIndices.begin(), originalCmapIndices.end());
  std::vector<int> currentToNew(cmapIn.size(), -1);
  // Add the grid terms in original parameter file order
  // and create map from assigned index to new index.
  cmapGrids.reserve( originalCmapIndices.size() );
  for (unsigned int idx = 0; idx < originalCmapIndices.size(); idx++) {
    mprintf("DEBUG: Will change cmap parameter index %i to %u\n",
            originalCmapIndices[idx]+1, idx+1);
    currentToNew[ originalCmapIndices[idx] ] = idx;
    cmapGrids.push_back( cmapIn[ originalCmapIndices[idx] ] );
  }
  mprintf("DEBUG: Assigned CMAP parameters:\n");
  for (CmapGridArray::const_iterator it = cmapGrids.begin(); it != cmapGrids.end(); ++it)
    mprintf("DEBUG:\t\t%li : %s\n", it-cmapGrids.begin(), it->Title().c_str());
  // Renumber the parameter indices
  for (CmapArray::iterator it = cmapTerms.begin(); it != cmapTerms.end(); ++it) {
    int newIdx = currentToNew[ it->Idx() ];
    if (newIdx < 0) {
      mprinterr("Internal Error: CMAP term index is not mapped.\n");
      return 1;
    }
    it->SetIdx( newIdx );
  }
  return 0;
}

/** Assign CMAP parameters to existing CMAP terms. */
int AssignParams::AssignCmapParams(CmapArray& cmapTerms, CmapParmHolder const& cmapIn,
                               CmapGridArray& cmapGrids)
const
{
  // Keep track of residue names each cmap applies to
  typedef std::vector<std::string> Sarray;
  std::vector<Sarray> CmapResNames;
  std::vector<Sarray> CmapAtomNames;
  // The LEaP convention is to number the CMAP parameters in the same
  // order as they are listed in the original parameter file. This
  // variable keeps track of those indices.
  std::vector<int> originalCmapIndices;
  for (CmapArray::iterator cm = cmapTerms.begin(); cm != cmapTerms.end(); ++cm)
  {
    // Get residue name for A2
    Atom const& A2 = atoms_[cm->A2()];
    // TODO make sure A2-A4 in same residue?
    NameType const& rn2 = residues_[A2.ResNum()].Name();
    // Is this residue already in cmapGrids?
    int cidx = -1;
    for (unsigned int idx = 0; idx != CmapResNames.size(); idx++) {
      if ( MatchesResName(CmapResNames[idx], rn2.Truncated()) ) {
        if (atoms_[cm->A1()].Name() == CmapAtomNames[idx][0] &&
            atoms_[cm->A2()].Name() == CmapAtomNames[idx][1] &&
            atoms_[cm->A3()].Name() == CmapAtomNames[idx][2] &&
            atoms_[cm->A4()].Name() == CmapAtomNames[idx][3] &&
            atoms_[cm->A5()].Name() == CmapAtomNames[idx][4])
        {
          cidx = originalCmapIndices[idx];
          break;
        }
      }
    }
    if (cidx > -1) {
      mprintf("DEBUG: Potential existing cmap %i found for %s (%li)\n", cidx, TruncResNameNum(A2.ResNum()).c_str(), cm-cmapTerms.begin());
      mprintf("DEBUG:\t\t%i - %i - %i - %i - %i %i\n", cm->A1()+1, cm->A2()+1, cm->A3()+1, cm->A4()+1, cm->A5()+1, cidx+1);
      cm->SetIdx( cidx );
    }
    // If not already in cmapGrids, check cmapIn
    if (cidx == -1) {
      for (unsigned int idx = 0; idx != cmapIn.size(); idx++) {
        if ( MatchesResName( cmapIn[idx].ResNames(), rn2.Truncated() ) ) {
          if (atoms_[cm->A1()].Name() == cmapIn[idx].AtomNames()[0] &&
              atoms_[cm->A2()].Name() == cmapIn[idx].AtomNames()[1] &&
              atoms_[cm->A3()].Name() == cmapIn[idx].AtomNames()[2] &&
              atoms_[cm->A4()].Name() == cmapIn[idx].AtomNames()[3] &&
              atoms_[cm->A5()].Name() == cmapIn[idx].AtomNames()[4])
          {
            cidx = (int)idx;
            CmapResNames.push_back( cmapIn[idx].ResNames() );
            CmapAtomNames.push_back( cmapIn[idx].AtomNames() );
            originalCmapIndices.push_back( idx );
            break;
          }
        }
      }
      if (cidx > -1) {
        mprintf("DEBUG: Potential new cmap %i found for %s (%li)\n", cidx, TruncResNameNum(A2.ResNum()).c_str(), cm-cmapTerms.begin());
        mprintf("DEBUG:\t\t%i - %i - %i - %i - %i %i\n", cm->A1()+1, cm->A2()+1, cm->A3()+1, cm->A4()+1, cm->A5()+1, cidx+1);
        //mprintf("DEBUG:\t\t%s - %s - %s - %s - %s\n",
        //        AtomMaskName(dih->A1()).c_str(),
        //        AtomMaskName(dih->A2()).c_str(),
        //        AtomMaskName(dih->A3()).c_str(),
        //        AtomMaskName(dih->A4()).c_str(),
        //        AtomMaskName(a5).c_str());
        cm->SetIdx( cidx );
      }
    }
  } // END loop over cmap terms

  return remap_cmap_indices(originalCmapIndices, cmapGrids, cmapTerms, cmapIn);
}

/** Assign CMAP parameters. Also generates CMAP terms from dihedrals. */
int AssignParams::AssignCmapParams(DihedralArray const& allDih, CmapParmHolder const& cmapIn,
                               CmapGridArray& cmapGrids, CmapArray& cmapTerms)
const
{
  // TODO check for duplicates, look for reverse match?
  // Keep track of residue names each cmap applies to
  //typedef std::vector<NameType> Narray;
  typedef std::vector<std::string> Sarray;
  std::vector<Sarray> CmapResNames;
  std::vector<Sarray> CmapAtomNames;
  // The LEaP convention is to number the CMAP parameters in the same
  // order as they are listed in the original parameter file. This
  // variable keeps track of those indices.
  std::vector<int> originalCmapIndices;
  // TODO combine with AssignDihedrals?
  for (DihedralArray::const_iterator dih = allDih.begin(); dih != allDih.end(); ++dih)
  {
    // Ignore end (repeated) or improper dihedrals
    if (dih->IsImproper() || dih->Skip14()) continue;
    // Get residue name for A2
    Atom const& A2 = atoms_[dih->A2()];
    // TODO make sure A2-A4 in same residue?
    NameType const& rn2 = residues_[A2.ResNum()].Name();
    // Is this residue already in cmapGrids?
    int cidx = -1;
    int a5 = -1;
    //for (unsigned int idx = 0; idx != cmapGrids.size(); idx++) {
    //  if (cmapGrids[idx].MatchesResName( rn2.Truncated() )) {
    for (unsigned int idx = 0; idx != CmapResNames.size(); idx++) {
      if ( MatchesResName(CmapResNames[idx], rn2.Truncated()) ) {
        //a5 = cmap_anames_match(*dih, atoms_, cmapGrids[idx]);
        a5 = cmap_anames_match(*dih, atoms_, CmapAtomNames[idx]);
        if (a5 > -1) {
          //cidx = (int)idx;
          cidx = originalCmapIndices[idx];
          break;
        }
      }
    }
    if (cidx > -1) {
      mprintf("DEBUG: Potential existing cmap %i found for %s (%li)\n", cidx, TruncResNameNum(A2.ResNum()).c_str(), dih-allDih.begin());
      mprintf("DEBUG:\t\t%i - %i - %i - %i - %i %i\n", dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1, a5+1, cidx+1);
//      mprintf("DEBUG:\t\t%s - %s - %s - %s - %s\n",
//              AtomMaskName(dih->A1()).c_str(),
//              AtomMaskName(dih->A2()).c_str(),
//              AtomMaskName(dih->A3()).c_str(),
//              AtomMaskName(dih->A4()).c_str(),
//              AtomMaskName(a5).c_str());
      cmapTerms.push_back( CmapType(dih->A1(), dih->A2(), dih->A3(), dih->A4(), a5, cidx) );
    }
    // If not already in cmapGrids, check cmapIn
    if (cidx == -1) {
      //int nidx = -1;
      for (unsigned int idx = 0; idx != cmapIn.size(); idx++) {
        //if (cmapIn[idx].MatchesResName( rn2.Truncated() )) {
        if ( MatchesResName( cmapIn[idx].ResNames(), rn2.Truncated() ) ) {
          //a5 = cmap_anames_match(*dih, atoms_, cmapIn[idx]);
          a5 = cmap_anames_match(*dih, atoms_, cmapIn[idx].AtomNames());
          if (a5 > -1) {
            //cidx = (int)cmapGrids.size();
            cidx = (int)idx;
            //cmapGrids.push_back( cmapIn[idx] );
            CmapResNames.push_back( cmapIn[idx].ResNames() );
            CmapAtomNames.push_back( cmapIn[idx].AtomNames() );
            originalCmapIndices.push_back( idx );
            //nidx = (int)idx;
            break;
          }
        }
      }
      if (cidx > -1) {
        mprintf("DEBUG: Potential new cmap %i found for %s (%li)\n", cidx, TruncResNameNum(A2.ResNum()).c_str(), dih-allDih.begin());
        mprintf("DEBUG:\t\t%i - %i - %i - %i - %i %i\n", dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1, a5+1, cidx+1);
      //mprintf("DEBUG:\t\t%s - %s - %s - %s - %s\n",
      //        AtomMaskName(dih->A1()).c_str(),
      //        AtomMaskName(dih->A2()).c_str(),
      //        AtomMaskName(dih->A3()).c_str(),
      //        AtomMaskName(dih->A4()).c_str(),
      //        AtomMaskName(a5).c_str());
        cmapTerms.push_back( CmapType(dih->A1(), dih->A2(), dih->A3(), dih->A4(), a5, cidx) );
      }
    }
  } // END loop over dihedrals
  return remap_cmap_indices(originalCmapIndices, cmapGrids, cmapTerms, cmapIn);
/*  mprintf("DEBUG: Cmap parameter indices:\n");
  for (unsigned int idx = 0; idx < originalCmapIndices.size(); idx++)
    mprintf("DEBUG:\t\tCurrent idx=%i  Actual idx=%u\n", originalCmapIndices[idx], idx);
  if (originalCmapIndices.empty()) {
    mprintf("DEBUG: No CMAP indices in %s\n", c_str());
    return 0;
  }
  std::sort(originalCmapIndices.begin(), originalCmapIndices.end());
  std::vector<int> currentToNew(cmapIn.size(), -1);
  // Add the grid terms in original parameter file order
  // and create map from assigned index to new index.
  cmapGrids.reserve( originalCmapIndices.size() );
  for (unsigned int idx = 0; idx < originalCmapIndices.size(); idx++) {
    mprintf("DEBUG: Will change cmap parameter index %i to %u\n",
            originalCmapIndices[idx]+1, idx+1);
    currentToNew[ originalCmapIndices[idx] ] = idx;
    cmapGrids.push_back( cmapIn[ originalCmapIndices[idx] ] );
  }
  mprintf("DEBUG: Assigned CMAP parameters:\n");
  for (CmapGridArray::const_iterator it = cmapGrids.begin(); it != cmapGrids.end(); ++it)
    mprintf("DEBUG:\t\t%li : %s\n", it-cmapGrids.begin(), it->Title().c_str());
  // Renumber the parameter indices
  for (CmapArray::iterator it = cmapTerms.begin(); it != cmapTerms.end(); ++it) {
    int newIdx = currentToNew[ it->Idx() ];
    if (newIdx < 0) {
      mprinterr("Internal Error: CMAP term index is not mapped.\n");
      return 1;
    }
    it->SetIdx( newIdx );
  }
  return 0;*/
}

/** Replace existing parameters with the given parameter set. */
int AssignParams::AssignParameters(Topology& topOut, ParameterSet const& set0) const {

  // Bond parameters
  mprintf("\tAssigning bond parameters.\n");
  bondparm_.clear();
  // Regenerate bond array in LEaP order
  bonds_.clear();
  bondsh_.clear();
  BondArray allBonds = Cpptraj::Structure::GenerateBondArray(residues_, atoms_);
  AssignBondParm( topOut, set0.BP(), allBonds, bondparm_, "bond" );
  for (BondArray::const_iterator bnd = allBonds.begin(); bnd != allBonds.end(); ++bnd)
    AddToBondArrays( *bnd );
  // Angle parameters
  mprintf("\tAssigning angle parameters.\n");
  angleparm_.clear();
  // Regenerate angle array in LEaP order
  angles_.clear();
  anglesh_.clear();
  AngleArray allAngles = Cpptraj::Structure::GenerateAngleArray(residues_, atoms_);
  allAngles = AssignAngleParm( set0.AP(), allAngles );
  for (AngleArray::const_iterator ang = allAngles.begin(); ang != allAngles.end(); ++ang)
    AddToAngleArrays( *ang );
  // Dihedral parameters
  mprintf("\tAssigning dihedral parameters.\n");
  dihedralparm_.clear();
  // Regenerate dihedral array in LEaP order
  dihedrals_.clear();
  dihedralsh_.clear();
  DihedralArray allDihedrals = Cpptraj::Structure::GenerateDihedralArray(residues_, atoms_);
  // If we need CMAP terms, do it here before the dihedrals array is modified
  if (!set0.CMAP().empty()) {
    cmap_.clear();
    cmapGrid_.clear();
    mprintf("\tAssigning CMAP parameters.\n");
    AssignCmapParams( allDihedrals, set0.CMAP(), cmapGrid_, cmap_ );
  } 
  // Now modify the dihedrals for any multiplicities
  allDihedrals = AssignDihedralParm( set0.DP(), set0.IP(), set0.AT(), allDihedrals, false );
  for (DihedralArray::const_iterator dih = allDihedrals.begin(); dih != allDihedrals.end(); ++dih)
    AddToDihedralArrays( *dih );
  // Urey-Bradley
  mprintf("\tAssigning Urey-Bradley parameters.\n");
  AssignUBParams( set0.UB() );
  if (!impropers_.empty()) {
    // Charmm Improper parameters
    mprintf("\tAssigning CHARMM improper parameters.\n");
    AssignImproperParams( set0.IP() );
  } else {
    // Amber improper parameters
    mprintf("\tAssigning improper parameters.\n");
    DihedralArray allImpropers = Cpptraj::Structure::GenerateImproperArray(residues_, atoms_);
    allImpropers = AssignDihedralParm( set0.DP(), set0.IP(), set0.AT(), allImpropers, true );
    for (DihedralArray::const_iterator imp = allImpropers.begin(); imp != allImpropers.end(); ++imp)
      AddToDihedralArrays( *imp );
  }
  // Atom types
  mprintf("\tAssigning atom type parameters.\n");
  AssignAtomTypeParm( topOut.ModifyAtoms(), set0.AT() );
  // LJ 6-12
  mprintf("\tAssigning nonbond parameters.\n");
  AssignNonbondParams( set0.AT(), set0.NB(), set0.NB14(), set0.LJC(), set0.HB(), debug_ );
  mprintf("DEBUG: CMAP size %zu\n", set0.CMAP().size());

  return 0;
}

/** Update/add to parameters in this topology with those from given set.
  */
int AssignParams::UpdateParameters(Topology& topOut, ParameterSet const& set1) const {
  ParameterSet set0 = topOut.GetParameters();

  //set1.Summary(); // DEBUG
  // Check TODO is this necessary?
  if (set0.AT().size() < 1)
    mprintf("Warning: No atom type information in '%s'\n", c_str());
  if (debug_ > 0) {
    mprintf("DEBUG: Saving original parameters in originalp.dat, new parameters in newp.dat\n");
    set0.Debug("originalp.dat");
  }
  // Update existing parameters with new parameters
  ParameterSet::UpdateCount UC;
  if (set0.UpdateParamSet( set1, UC, debug_, debug_ )) { // FIXME verbose
    mprinterr("Error: Could not merge topology '%s' parameters with '%s' parameters.\n",
              c_str(), set1.ParamSetName().c_str());
    return 1;
  }
  mprintf("DEBUG: Updated parameters.\n");
  set0.Summary();
  //set0.Debug();

//  unsigned int updateCount;
  // Bond parameters
//  updateCount = UpdateParameters< ParmHolder<BondParmType> >(set0.BP(), set1.BP(), "bond");
  if (UC.nBondsUpdated_ > 0) {
    mprintf("\tRegenerating bond parameters.\n");
    AssignBondParams( set0.BP() );
  }
  // Angle parameters
//  updateCount = UpdateParameters< ParmHolder<AngleParmType> >(set0.AP(), set1.AP(), "angle");
  if (UC.nAnglesUpdated_ > 0) {
    mprintf("\tRegenerating angle parameters.\n");
    AssignAngleParams( set0.AP() );
  }
  // Dihedral parameters
//  updateCount = UpdateParameters< DihedralParmHolder >(set0.DP(), set1.DP(), "dihedral");
  if (UC.nDihedralsUpdated_ > 0) {
    mprintf("\tRegenerating dihedral parameters.\n");
    AssignDihedralParams( set0.DP(), set0.IP(), set0.AT() );
  }
  // Urey-Bradley
//  updateCount = UpdateParameters< ParmHolder<BondParmType> >(set0.UB(), set1.UB(), "Urey-Bradley");
  if (UC.nUreyBradleyUpdated_ > 0) {
    mprintf("\tRegenerating UB parameters.\n");
    AssignUBParams( set0.UB() );
  }
  // Improper parameters
//  updateCount = UpdateParameters< ParmHolder<DihedralParmType> >(set0.IP(), set1.IP(), "improper");
  if (UC.nImpropersUpdated_ > 0) {
    mprintf("\tRegenerating improper parameters.\n");
    AssignImproperParams( set0.IP() );
  }
  // Atom types
//  updateCount = UpdateParameters< ParmHolder<AtomType> >(set0.AT(), set1.AT(), "atom type");
  if (UC.nAtomTypeUpdated_ > 0) {
    mprintf("\tRegenerating atom type parameters.\n");
    AssignAtomTypeParm( topOut.ModifyAtoms(), set0.AT() );
  }
//  updateCount += UpdateParameters< ParmHolder<NonbondType> >(set0.NB(), set1.NB(), "LJ A-B");
  if (UC.nAtomTypeUpdated_ > 0 || UC.nLJCUpdated_ > 0 ||
      UC.nLJparamsUpdated_ > 0 || UC.nLJ14paramsUpdated_ > 0)
  {
    mprintf("\tRegenerating nonbond parameters.\n");
    AssignNonbondParams( set0.AT(), set0.NB(), set0.NB14(), set0.LJC(), set0.HB(), debug_ );
  }
  // CMAP
  if (UC.nCmapUpdated_ > 0) {
    cmapGrid_.clear();
    mprintf("\tRegenerating CMAP parameters.\n");
    AssignCmapParams(cmap_, set0.CMAP(), cmapGrid_);
  }
  // TODO LJ14

  if (debug_ > 0) set0.Debug("newp.dat");
  return 0;
}

