#ifndef INC_CSA_H
#define INC_CSA_H
namespace Cpptraj {
namespace CSA {

enum DistType {
  NO_DIST = 0,  // No distance specified
  DIST_ANGLE,   // Cumulative angle distance
  DIST_RMSD,    // RMSD
  DIST_TMSCORE, // TM score
  DIST_FRECHET  // Frechet distance
};

enum ScoreType {
  SCORE_NONE = 0,
  SCORE_PE,
  SCORE_OM
};

enum ManipType {
  MANIP_NONE = 0,
  MANIP_IC,
  MANIP_CART
};

enum DcutType {
  DCUT_NONE = 0,
  DCUT_LINEAR,
  DCUT_GEOM
};

}
}
#endif
