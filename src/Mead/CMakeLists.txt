#CMake buildfile for CPPTRAJ Mead subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/MeadAtoms.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadGrid.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadInterface.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MeadOpts.cpp
)
