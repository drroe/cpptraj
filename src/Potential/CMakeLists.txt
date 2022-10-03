#CMake buildfile for CPPTRAJ Potential subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/Constraints.cpp
  ${CMAKE_CURRENT_LIST_DIR}/EnergyArray.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MdOpts.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Minimize_SteepestDescent.cpp
  ${CMAKE_CURRENT_LIST_DIR}/OMM_helpers.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PotentialFunction.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PotentialTerm_Angle.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PotentialTerm_Bond.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PotentialTerm.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PotentialTerm_Dihedral.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PotentialTerm_LJ_Coulomb.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PotentialTerm_OpenMM.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PotentialTerm_Replicate.cpp
)
