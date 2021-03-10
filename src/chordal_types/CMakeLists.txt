add_library(srrg2_chordal_library 
  variables_chordal_initialization.cpp
  variables_chordal_initialization.h
  factors_chordal_initialization.cpp
  factors_chordal_initialization.h
  chordal_initializer.cpp
  chordal_initializer.h
  )

target_link_libraries(srrg2_chordal_library
  srrg2_solver_core_library
  srrg2_solver_types3d_library
  )
