cmake_minimum_required(VERSION 2.8.3)
project(srrg2_hipe)

find_package(catkin REQUIRED COMPONENTS
  srrg2_core
  srrg2_solver
  srrg_cmake_modules
)


include(${srrg_cmake_modules_INCLUDE_DIRS}/CMakeCompileOptions.txt)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -pedantic")
message("${PROJECT_NAME}|compiling with these CXX flags: ${CMAKE_CXX_FLAGS}")
add_definitions(-DTEST_DATA_FOLDER="${CMAKE_CURRENT_SOURCE_DIR}/data/")
message("${PROJECT_NAME}|data dir    = ${CMAKE_CURRENT_SOURCE_DIR}/data/")

find_package(Eigen3 REQUIRED)

catkin_package(
  INCLUDE_DIRS
  src
)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${EIGEN3_INCLUDE_DIR}
  ${catkin_INCLUDE_DIRS}
)

link_directories(${catkin_LIBRARY_DIRS})

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(apps)
