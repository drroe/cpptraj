#ifndef INC_MASK_TO_XYZ_H
#define INC_MASK_TO_XYZ_H
#include "AtomMask.h"
#include "Frame.h"
template <typename T>
std::vector<T> mask_to_xyz(AtomMask const& Mask, Frame const& frm)
{
  std::vector<T> outerxyz;
  outerxyz.reserve( Mask.Nselected()*3 );
  for (AtomMask::const_iterator at = Mask.begin(); at != Mask.end(); ++at) {
    const double* xyz = frm.XYZ( *at );
    outerxyz.push_back( (T)xyz[0] );
    outerxyz.push_back( (T)xyz[1] );
    outerxyz.push_back( (T)xyz[2] );
  }
  return outerxyz;
}
#endif
