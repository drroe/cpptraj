#CMake buildfile for CPPTRAJ Cluster subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/Algorithm.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Algorithm_DBscan.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Algorithm_DPeaks.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Algorithm_HierAgglo.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Algorithm_Kmeans.cpp
  ${CMAKE_CURRENT_LIST_DIR}/BestReps.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CentroidArray.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Centroid_Coord.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Cframes.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Cmatrix_Binary.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Cmatrix_NC.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Control.cpp
  ${CMAKE_CURRENT_LIST_DIR}/DBI.cpp
  ${CMAKE_CURRENT_LIST_DIR}/DrawGraph.cpp
  ${CMAKE_CURRENT_LIST_DIR}/DynamicMatrix.cpp
  ${CMAKE_CURRENT_LIST_DIR}/List.cpp
  ${CMAKE_CURRENT_LIST_DIR}/MetricArray.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Metric_DME.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Metric_Matrix2D.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Metric_QuatRMSD.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Metric_RMS.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Metric_SRMSD.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Metric_Scalar.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Metric_Torsion.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Node.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Output.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PseudoF.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Results_Coords.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Sieve.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Silhouette.cpp
)
