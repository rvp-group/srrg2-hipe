add_executable(partition_covariance_projection_test partition_covariance_projection_test.cpp)
target_link_libraries(partition_covariance_projection_test
  srrg2_partition_library
  ${catkin_LIBRARIES}
  )

add_executable(partition_manager_test partition_manager_test.cpp)
target_link_libraries(partition_manager_test
  srrg2_partition_library
  ${catkin_LIBRARIES}
  )

add_executable(chordal_initialization_example chordal_initialization_example.cpp)
target_link_libraries(chordal_initialization_example
  srrg2_chordal_library
  ${catkin_LIBRARIES}
  )

add_executable(hierarchical_initialization_example hierarchical_initialization_example.cpp)
target_link_libraries(hierarchical_initialization_example
  srrg2_hierarchical_optimizer_library
  ${catkin_LIBRARIES}
  )

add_executable(bipartition_example bipartition_example.cpp)
target_link_libraries(bipartition_example
  srrg2_partition_library
  srrg2_splitting_library
  ${catkin_LIBRARIES}
  )

add_executable(multi_partition_example multi_partition_example.cpp)
target_link_libraries(multi_partition_example
  srrg2_partition_library
  srrg2_splitting_library
  ${catkin_LIBRARIES}
  )

add_executable(reoptimize_partitions_test reoptimize_partitions_test.cpp)
target_link_libraries(reoptimize_partitions_test
  srrg2_partition_library
  srrg2_splitting_library
  ${catkin_LIBRARIES}
  )
