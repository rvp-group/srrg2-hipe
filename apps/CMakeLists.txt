add_executable(graph_optimizer graph_optimizer.cpp)
target_link_libraries(graph_optimizer 
  srrg2_partition_library
  srrg2_splitting_library
  srrg2_hierarchical_optimizer_library
  ${catkin_LIBRARIES})

add_executable(geodesic_to_chordal_converter geodesic_to_chordal_converter.cpp)
target_link_libraries(geodesic_to_chordal_converter 
  ${catkin_LIBRARIES})

add_executable(chordal_to_geodesic_converter chordal_to_geodesic_converter.cpp)
target_link_libraries(chordal_to_geodesic_converter 
  ${catkin_LIBRARIES})

add_executable(compute_metrics compute_metrics.cpp)
target_link_libraries(compute_metrics 
  ${catkin_LIBRARIES})

add_executable(align_pose_graph align_pose_graph.cpp)
target_link_libraries(align_pose_graph 
  ${catkin_LIBRARIES})
