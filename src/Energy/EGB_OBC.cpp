#include "EGB_OBC.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EGB_OBC::EGB_OBC() {}

/** Initialize alpha beta gamma arrays */
int EGB_OBC::InitABG(GBtype typeIn, int natomIn)
{
  gbalpha_.clear();
  gbbeta_.clear();
  gbgamma_.clear();

  /* Set up GB/OBC parameters: */
  if (typeIn == OBC2 || typeIn == OBC5 || typeIn == NECK || typeIn == NECK2) {
    gbalpha_.assign(natomIn, 0);
    gbbeta_.assign(natomIn, 0);
    gbgamma_.assign(natomIn, 0);
  }
  if (typeIn == OBC2) { 
      for (int i = 0; i < natomIn; i++) {
         gbalpha_[i] = 0.8; 
         gbbeta_[i] = 0.0;
         gbgamma_[i] = 2.909125;
      }
  } else if (typeIn == OBC5) {
      for (int i = 0; i < natomIn; i++) {
           gbalpha_[i] = 1.0;
           gbbeta_[i] = 0.8;
           gbgamma_[i] = 4.85;
       }
  }// else if (gb == 7) {

  return 0;
}
