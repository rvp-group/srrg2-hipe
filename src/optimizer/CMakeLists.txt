add_library(srrg2_hierarchical_optimizer_library 
  hierarchical_optimizer.cpp
  hierarchical_optimizer.h
  instances.cpp
  instances.h)

target_link_libraries(srrg2_hierarchical_optimizer_library 
  srrg2_partition_library 
  srrg2_splitting_library
  srrg2_chordal_library
)
