add_library(srrg2_hipe_factors_library
  se3_chordal_factor.cpp
  se3_chordal_factor.h
  se3_langevin_factor.cpp
  se3_langevin_factor.h
  instances.cpp
  instances.h
  )

target_link_libraries(srrg2_hipe_factors_library
  srrg2_solver_core_library
  srrg2_solver_types3d_library
  )
