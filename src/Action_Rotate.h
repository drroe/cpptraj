#ifndef INC_ACTION_ROTATE_H
#define INC_ACTION_ROTATE_H
#include "Action.h"
#include "ParallelSetFrameNum.h"
class DataSet_Mat3x3;
/// Rotate coordinates or calculate rotations from rotation matrices
class Action_Rotate : public Action {
  public:
    Action_Rotate();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Rotate(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    /// Get 3x3 matrix DataSet from name
    int Get3x3Set(DataSetList const&, std::string const&);
    /// Setup output DataSets
    int SetupOutputSets(DataSetList&, std::string const&, DataFile*);

    enum ModeType { ROTATE = 0, DATASET, AXIS, CALC };
    Matrix_3x3 RotMatrix_;      ///< Rotation matrix.
    AtomMask mask_;             ///< Mask of atoms to rotate.
    AtomMask axis0_;            ///< Mask of atoms defining 1 end of rotation axis.
    AtomMask axis1_;            ///< Mask of atoms defining other end of rotation axis.
    DataSet_Mat3x3* rmatrices_; ///< DataSet containing rotation matrices to use.
    double delta_;              ///< Degrees to rotate if mode is AXIS
    ModeType mode_;             ///< Mode to use.
    bool inverse_;              ///< If true perform an inverse rotation.
    bool all_atoms_selected_;   ///< If true all atoms selected for rotation
    DataSet* dsout_tx_;         ///< Hold output theta X (calc)
    DataSet* dsout_ty_;         ///< Hold output theta Y (calc)
    DataSet* dsout_tz_;         ///< Hold output theta Z (calc)
    DataSet* dsout_t_;          ///< Hold output theta (calc)
    Cpptraj::ParallelSetFrameNum parallelNum_;
};
#endif
