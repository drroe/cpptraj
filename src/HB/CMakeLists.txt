#CMake buildfile for CPPTRAJ HB subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/Bridge.cpp
  ${CMAKE_CURRENT_LIST_DIR}/HbCalc.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Hbond.cpp
)
