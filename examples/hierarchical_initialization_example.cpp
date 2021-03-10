#include <optimizer/hierarchical_optimizer.h>
#include <partitions/partition_creation_utils.h>
#include <partitions/partition_manager.h>
#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/solver_core/iteration_algorithm_ddl.h>
#include <srrg_solver/solver_core/solver.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_system_utils/parse_command_line.h>
#include <srrg_system_utils/system_utils.h>

const std::string executable_name("hierarchical_initialization_example");
#define LOG std::cerr << FG_YELLOW(executable_name) << "|"
using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;
using namespace srrg2_hipe;

int main(int argc, char** argv) {
  variables_and_factors_3d_registerTypes();
  ParseCommandLine cmd_line(argv);
  ArgumentString param_input_file(
    &cmd_line, "i", "input-file", "file containing graph to optimize (*.boss)", "");
  cmd_line.parse();
  if (!param_input_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(executable_name + "|ERROR, no input file specified");
  }

  const std::string& input_file = param_input_file.value();
  FactorGraphPtr graph          = FactorGraph::read(input_file);
  if (!graph) {
    throw std::runtime_error(executable_name + "| no graph");
  }
  HierarchicalOptimizer init;
  init.param_verbose.setValue(false);
  auto& global_solver = init.param_solver.value();
  IterationAlgorithmBasePtr dl(new IterationAlgorithmDL);
  global_solver->setName("global_solver");
  global_solver->param_max_iterations.pushBack(10);
  global_solver->param_max_iterations.pushBack(10);
  global_solver->param_algorithm.setValue(dl);
  init.param_low_level_alg.setValue(dl);

  auto& partition_manager = init.param_partition_manager.value();
  partition_manager->setName("partition_manager");
  partition_manager->param_verbose.setValue(false);
  // tg condensed factors pose pose 3D
  VirtualFactorCreatorBasePtr condensed_factor_se3_posepose(new SE3PosePoseVirtualFactorCreator);
  partition_manager->param_covariance_projectors.pushBack(condensed_factor_se3_posepose);
  auto& partition_solver = partition_manager->param_solver.value();
  partition_solver->param_max_iterations.pushBack(50);
  partition_solver->param_algorithm.setValue(dl);
  Solver solver;
  solver.param_max_iterations.pushBack(1);
  solver.setGraph(graph);
  solver.compute();
  LOG << solver.iterationStats() << std::endl;
  init.setGraph(graph);
  double start = getTime();
  init.initializeGraph();
  double end = getTime();
  solver.setGraph(graph);
  solver.compute();
  LOG << solver.iterationStats() << std::endl;
  LOG << "Done in " << FG_YELLOW(end - start) << " sec" << std::endl;
  return 0;
}
