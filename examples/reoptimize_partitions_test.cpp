#include <partitions/partition_creation_utils.h>
#include <partitions/partition_manager.h>
#include <random>
#include <splitting/breadth_first_splitter.h>
#include <splitting/nested_dissection_splitter.h>
#include <srrg_geometry/geometry3d.h>
#include <srrg_solver/solver_core/instances.h>
#include <srrg_solver/solver_core/iteration_algorithm_lm.h>
#include <srrg_solver/variables_and_factors/types_2d/instances.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_system_utils/shell_colors.h>
#include <srrg_system_utils/system_utils.h>

const std::string executable_name("reoptimize partitions test");
#define LOG std::cerr << FG_YELLOW(executable_name) << "|"
using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;
using namespace srrg2_hipe;

const std::string data_folder(TEST_DATA_FOLDER);

int main(int argc, char** argv) {
  variables_and_factors_2d_registerTypes();
  variables_and_factors_3d_registerTypes();
  FactorGraphPtr graph = FactorGraph::read(data_folder + "boss/grid3D.boss");
  PartitionManager manager;
  manager.setGraph(graph);

  manager.param_solver.value()->param_max_iterations.pushBack(50);
  manager.param_solver.value()->param_algorithm.setValue(
    IterationAlgorithmBasePtr(new IterationAlgorithmLM));
  manager.param_covariance_projectors.pushBack(
    VirtualFactorCreatorBasePtr(new SE3PosePoseVirtualFactorCreator));

  BreadthFirstSplitter splitter;
  splitter.param_min_partition_variables.setValue(150);
  splitter.setPartitionManager(&manager);
  SystemUsageCounter::tic();
  splitter.compute();
  double split_and_optimize_time = SystemUsageCounter::toc();

  Solver solver;
  solver.param_termination_criteria.setValue(nullptr);
  solver.param_max_iterations.pushBack(0);
  solver.param_max_iterations.pushBack(15);
  std::shared_ptr<IterationAlgorithmGN> gn(new IterationAlgorithmGN);
  gn->param_damping.setValue(0.1);
  solver.param_algorithm.setValue(gn);
  solver.setGraph(graph);
  SystemUsageCounter::tic();
  solver.compute();
  LOG << solver.iterationStats() << std::endl;
  solver.param_max_iterations.setValue(0, 15);
  solver.param_max_iterations.setValue(1, 0);
  solver.compute();
  double global_optimization_time = SystemUsageCounter::toc();
  LOG << solver.iterationStats() << std::endl;
  LOG << FG_YELLOW("Partitions splitting and optimization : ") << split_and_optimize_time
      << std::endl;
  LOG << FG_YELLOW("Num partitions : ") << manager.partitions().size() << std::endl;
  LOG << FG_YELLOW("Optimization time : ") << global_optimization_time << std::endl;
  return 0;
}

