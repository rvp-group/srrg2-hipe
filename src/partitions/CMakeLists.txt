add_library(srrg2_partition_library
  partition.h
  partition.cpp
  partition_manager.h
  partition_manager.cpp
  partition_creation_utils.h
  partition_creation_utils.cpp
  instances.h
  instances.cpp
  )

target_link_libraries(srrg2_partition_library
  srrg2_chordal_library
  srrg2_solver_core_library
  srrg2_data_structures_library)
