#include "GridAction.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "StringRoutines.h"

using namespace Cpptraj;

/** CONSTRUCTOR */
GridAction::GridAction() :
  gridOffsetType_(NO_OFFSET),
  increment_(1.0)
{}

// GridAction::HelpText
const char* GridAction::HelpText =
  "\t{ data <dsname> | boxref <ref name/tag> <nx> <ny> <nz> |\n"
  "\t  <nx> <dx> <ny> <dy> <nz> <dz>\n"
  "\t  [ { gridcenter <cx> <cy> <cz> |\n"
  "\t      boxcenter |\n"
  "\t      maskcenter <mask> |\n"
  "\t      rmsfit <mask> [noxalign]} ]\n"
  "\t[box|origin|center <mask>] [negative] [name <gridname>]";

static inline void CheckEven(int& N, char dir) {
  if (N % 2 == 1) {
    ++N;
    mprintf("Warning: number of grid points must be even. Incrementing N%c by 1 to %u\n", dir, N);
  }
}

/** Determine grid move type based on keywords. */
int GridAction::determineMoveType(ArgList& argIn, GridMover::MoveType& gridMoveType,
                                  bool& specifiedCenter)
{
  if (argIn.hasKey("boxcenter")) {
    specifiedCenter = true;
    gridMoveType = GridMover::TO_BOX_CTR;
  } else {
    std::string maskCenterArg = argIn.GetStringKey("maskcenter");
    std::string rmsFitArg = argIn.GetStringKey("rmsfit");
    if (!maskCenterArg.empty()) {
      specifiedCenter = true;
      gridMoveType = GridMover::TO_MASK_CTR;
      if (centerMask_.SetMaskString( maskCenterArg )) return 1;
    } else if (!rmsFitArg.empty()) {
      specifiedCenter = true;
      gridMoveType = GridMover::RMS_FIT;
      if (centerMask_.SetMaskString( rmsFitArg )) return 1;
    }
  }
  return 0;
}

// GridAction::GridInit()
DataSet_GridFlt* GridAction::GridInit(const char* callingRoutine, ArgList& argIn, DataSetList& DSL) 
{
  DataSet_GridFlt* Grid = 0; 
  bool specifiedCenter = false;
  GridMover::MoveType gridMoveType = GridMover::NO_MOVE;
  bool do_x_align = !argIn.hasKey("noxalign");
  std::string dsname = argIn.GetStringKey("data");
  std::string refname = argIn.GetStringKey("boxref");
  if (!dsname.empty()) { 
    // Get existing grid dataset
    Grid = (DataSet_GridFlt*)DSL.FindSetOfType( dsname, DataSet::GRID_FLT );
    if (Grid == 0) {
      mprinterr("Error: %s: Could not find grid data set with name %s\n",
                callingRoutine, dsname.c_str());
      return 0;
    }
#   ifdef MPI
    // DataSet will probably need syncing after action.
    Grid->SetNeedsSync( true );
#   endif
    if (determineMoveType(argIn, gridMoveType, specifiedCenter)) return 0;
  } else if (!refname.empty()) {
    // Get grid parameters from reference structure box.
    DataSet_Coords_REF* REF = (DataSet_Coords_REF*)DSL.FindSetOfType( refname, DataSet::REF_FRAME );
    if (REF == 0) {
      mprinterr("Error: Reference '%s' not found.\n", refname.c_str());
      return 0;
    }
    if (!(REF->CoordsInfo().TrajBox().HasBox())) {
      mprinterr("Error: Reference '%s' does not have box information.\n", refname.c_str());
      return 0;
    }
    int nx = argIn.getNextInteger(-1);
    int ny = argIn.getNextInteger(-1);
    int nz = argIn.getNextInteger(-1);
    if (nx < 1 || ny < 1 || nz < 1) {
      mprinterr("Error:  %s: Invalid grid sizes\n", callingRoutine);
      return 0;
    }
    Grid = (DataSet_GridFlt*)DSL.AddSet( DataSet::GRID_FLT, argIn.GetStringKey("name"), "GRID" );
    if (Grid == 0) return 0;
    if (Grid->Allocate_N_O_Box(nx,ny,nz,Vec3(0.0),REF->RefFrame().BoxCrd())) return 0;
  } else {
    // Create new data set.
    // Get nx, dx, ny, dy, nz, dz
    int nx = argIn.getNextInteger(-1);
    double dx = argIn.getNextDouble(-1);
    int ny = argIn.getNextInteger(-1);
    double dy = argIn.getNextDouble(-1);
    int nz = argIn.getNextInteger(-1);
    double dz = argIn.getNextDouble(-1);
    if (nx < 1 || ny < 1 || nz < 1 ||
        dx < 0 || dy < 0 || dz < 0)
    {
      mprinterr("Error: %s: Invalid grid size/spacing.\n", callingRoutine);
      mprinterr("       nx=%i ny=%i nz=%i | dx=%.3f dy=%.3f dz=%.3f\n",
                nx, ny, nz, dx, dy, dz);
      return 0;
    }
    // For backwards compat., enforce even grid spacing.
    CheckEven( nx, 'X' ); CheckEven( ny, 'Y' ); CheckEven( nz, 'Z' );
    Vec3 gridctr(0.0, 0.0, 0.0);
    // Check if we want to specifically re-center the grid during DoAction
    gridMoveType = GridMover::NO_MOVE;
    if (argIn.hasKey("gridcenter")) {
      double cx = argIn.getNextDouble(0.0);
      double cy = argIn.getNextDouble(0.0);
      double cz = argIn.getNextDouble(0.0);
      gridctr.SetVec(cx, cy, cz);
      specifiedCenter = true;
    } else {
      if (determineMoveType( argIn, gridMoveType, specifiedCenter )) return 0;
    }
    Grid = (DataSet_GridFlt*)DSL.AddSet( DataSet::GRID_FLT, argIn.GetStringKey("name"), "GRID" );
    if (Grid == 0) return 0;
    // Set up grid from dims, center, and spacing
    // NOTE: # of grid points in each direction with be forced to be even.
    if (Grid->Allocate_N_C_D(nx, ny, nz, gridctr, Vec3(dx,dy,dz))) return 0;
  }
  // Determine offset, Box/origin/center
  gridOffsetType_ = NO_OFFSET;
  if (argIn.hasKey("box"))
    gridOffsetType_ = BOX_CENTER;
  else if (argIn.hasKey("origin"))
    gridOffsetType_ = NO_OFFSET;
  else if (argIn.Contains("center")) {
    std::string maskexpr = argIn.GetStringKey("center");
    if (maskexpr.empty()) {
      mprinterr("Error: 'center' requires <mask>\n");
      return 0;
    }
    if (centerMask_.SetMaskString( maskexpr )) return 0;
    gridOffsetType_ = MASK_CENTER;
  }
  if (specifiedCenter) {
    // If center was specified, do not allow an offset
    if (gridOffsetType_ != NO_OFFSET)
      mprintf("Warning: Grid offset args (box/center) not allowed with 'gridcenter'.\n"
              "Warning: No offset will be used.\n");
    gridOffsetType_ = NO_OFFSET;
  }
  // Negative
  if (argIn.hasKey("negative"))
    increment_ = -1.0;
  else
    increment_ = 1.0;
  // Set up the grid mover
  if (mover_.MoverInit(gridMoveType, do_x_align)) return 0;

  return Grid;
}

#ifdef MPI
int GridAction::ParallelGridInit(Parallel::Comm const& commIn, DataSet_GridFlt* Grid) {
  if (commIn.Rank() > 0) {
    // Since this set may have been read in (and therefore is synced),
    // zero it out on all non-master processes to avoid overcounting.
    std::fill( Grid->begin(), Grid->end(), 0.0 );
  }
  trajComm_ = commIn;
  mover_.MoverSetComm( trajComm_ );
  return 0;
}
#endif

// GridAction::GridInfo() 
void GridAction::GridInfo(DataSet_GridFlt const& grid) const {
  if (gridOffsetType_ == NO_OFFSET)
    mprintf("\tNo offset will be applied to points.\n");
  else if (gridOffsetType_ == BOX_CENTER)
    mprintf("\tOffset for points is box center.\n");
  else if (gridOffsetType_ == MASK_CENTER)
    mprintf("\tOffset for points is center of atoms in mask [%s]\n",
            centerMask_.MaskString());
  mover_.MoverInfo(centerMask_);
  if (increment_ > 0)
    mprintf("\tCalculating positive density.\n");
  else
    mprintf("\tCalculating negative density.\n");
  grid.GridInfo();
}

// GridAction::GridSetup()
int GridAction::GridSetup(Topology const& currentParm, CoordinateInfo const& cInfo) {
  // Check box
  if (gridOffsetType_ == BOX_CENTER) {
    if (!cInfo.TrajBox().Is_X_Aligned_Ortho()) {
      mprintf("Warning: Code to shift to the box center is not yet\n");
      mprintf("Warning: implemented for non-orthorhombic unit cells.\n");
      mprintf("Warning: No offset will be used.\n");
      gridOffsetType_ = NO_OFFSET; // TODO Error?
    }
  }
  // Set up mask
  if (centerMask_.MaskStringSet()) {
    if ( currentParm.SetupIntegerMask( centerMask_ ) ) return 1;
    centerMask_.MaskInfo();
    if ( centerMask_.None() ) {
      mprinterr("Error: No atoms selected for grid mask [%s]\n", centerMask_.MaskString());
      return 1;
    }
  }
  // Set up mover if needed
  if (mover_.MoverSetup(currentParm, centerMask_)) return 1;

  return 0;
}
