#include <partitions/partition_creation_utils.h>
#include <partitions/partition_manager.h>
#include <random>
#include <splitting/isoperimetrical_bipartitioner.h>
#include <srrg_geometry/geometry3d.h>
#include <srrg_solver/solver_core/instances.h>
#include <srrg_solver/variables_and_factors/types_2d/instances.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_system_utils/shell_colors.h>
#include <srrg_system_utils/system_utils.h>

const std::string executable_name("bipartition_example");
#define LOG std::cerr << FG_YELLOW(executable_name) << "|"
using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;
using namespace srrg2_hipe;

const std::string data_folder(TEST_DATA_FOLDER);

int main(int argc, char** argv) {
  variables_and_factors_2d_registerTypes();
  variables_and_factors_3d_registerTypes();
  FactorGraphPtr graph = FactorGraph::read(data_folder + "boss/sphere_bignoise_vertex3.boss");
  IsoperimetricBipartitioner splitter;
  splitter.setGraph(graph);
  splitter.computeEdgeWeights();
  SystemUsageCounter::tic();
  auto partition = splitter.compute();
  double time             = SystemUsageCounter::toc();
  int variable_left_side  = partition.first->variables().size();
  int variable_right_side = partition.second->variables().size();
  LOG << FG_RED("Total graph variables : ") << graph->variables().size()
      << " Boundary variables : " << partition.first->boundaryVariableIds().size()
      << " Left partition size : " << variable_left_side << " "
      << "Right partition size : " << variable_right_side << std::endl;
  int factor_left_side  = partition.first->factors().size();
  int factor_right_side = partition.second->factors().size();
  LOG << FG_RED("Total graph factors : ") << graph->factors().size()
      << " Left partition size : " << factor_left_side << " "
      << "Right partition size : " << factor_right_side << std::endl;
  LOG << FG_YELLOW("Execution time : ") << time << std::endl;

  PartitionManager manager;
  manager.setGraph(graph);
  manager.param_covariance_projectors.pushBack(
    VirtualFactorCreatorBasePtr(new SE3PosePoseVirtualFactorCreator));
  manager.param_solver.value()->param_max_iterations.pushBack(5);
  manager.addBoundaryVariablesToPartition(partition.first, partition.first->boundaryVariableIds());
  manager.addBoundaryVariablesToPartition(partition.second,
                                          partition.second->boundaryVariableIds());
  return 0;
}

