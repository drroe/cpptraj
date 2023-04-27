#CMake buildfile for CPPTRAJ Mead subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/MeadCalc.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadCalc_Multiflex.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadCalc_Potential.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadCalc_Solvate.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadError.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadGrid.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadOpts.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MultiFlexResults.cpp
)
