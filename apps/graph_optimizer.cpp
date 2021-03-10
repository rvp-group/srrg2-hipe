#include <optimizer/hierarchical_optimizer.h>
#include <optimizer/instances.h>
#include <partitions/instances.h>
#include <partitions/partition_creation_utils.h>
#include <partitions/partition_manager.h>
#include <splitting/instances.h>
#include <srrg_config/configurable_manager.h>
#include <srrg_solver/solver_core/instances.h>
#include <srrg_solver/solver_core/internals/linear_solvers/instances.h>
#include <srrg_solver/solver_core/iteration_algorithm_ddl.h>
#include <srrg_solver/solver_core/iteration_algorithm_lm.h>
#include <srrg_solver/variables_and_factors/types_2d/instances.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_system_utils/parse_command_line.h>
#include <srrg_system_utils/shell_colors.h>
#include <srrg_system_utils/system_utils.h>

const std::string exe_name("graph optimizer");
#define LOG std::cerr << FG_YELLOW(exe_name) << "|"
using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;
using namespace srrg2_hipe;
// tg register types
void initTypes() {
  variables_and_factors_2d_registerTypes();
  variables_and_factors_3d_registerTypes();
  solver_registerTypes();
  linear_solver_registerTypes();
  partition_registerTypes();
  splitting_registerTypes();
  hierarchicalOptimizer_registerTypes();
}

void generateConfig(const std::string& config_file_);

int main(int argc, char** argv) {
  initTypes();

  ParseCommandLine cmd_line(argv);
  ArgumentString param_input_file(
    &cmd_line, "i", "input-file", "file containing graph to optimize (*.boss)", "");
  ArgumentFlag param_gen_config(
    &cmd_line, "j", "generate-config", "generates a config file and quits");
  ArgumentString param_stats_file(
    &cmd_line, "s", "stats-file", "file where to save the solver statistics", "");
  ArgumentString param_config_file(
    &cmd_line, "c", "config-file", "config file to read/write", "hierarchical_solver.config");
  ArgumentString output_file(
    &cmd_line,
    "o",
    "output-file",
    "file where to save the output (hgraph_ means the higher level of the hierarchy)",
    "");
  cmd_line.parse();

  if (param_gen_config.isSet()) {
    generateConfig(param_config_file.value());
    LOG << "exit\n";
    return 0;
  }

  if (!param_input_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, no input file specified");
  }

  const std::string& input_file = param_input_file.value();
  FactorGraphPtr graph          = FactorGraph::read(input_file);
  if (!graph) {
    throw std::runtime_error(exe_name + "|ERROR, invalid graph file [ " + input_file + " ]");
  }
  LOG << "loaded [ " << FG_YELLOW(input_file) << " ] ----- "
      << "variables [ " << FG_YELLOW(graph->variables().size()) << " ] | "
      << "factors [ " << FG_YELLOW(graph->factors().size()) << " ]" << std::endl;
  ConfigurableManager manager;
  manager.read(param_config_file.value());

  HierarchicalOptimizerPtr optimizer = manager.getByName<HierarchicalOptimizer>("optimizer");
  if (!optimizer) {
    throw std::runtime_error(exe_name + "| ERROR, no CondensedGraphOptimizer was found in " +
                             param_config_file.value());
  }
  LOG << "Configuration in [ " << FG_YELLOW(param_config_file.value()) " ] loaded " << std::endl;
  optimizer->setGraph(graph);
  double start = getTime();
  optimizer->compute();
  double end = getTime();

  if (output_file.isSet() && !output_file.value().empty()) {
    const string hfile("hgraph_" + output_file.value());
    LOG << "saving output in [ " << FG_YELLOW(output_file.value()) << " ]  and "
        << "[ " << FG_YELLOW(hfile) " ]" << std::endl;
    optimizer->writeLevel(output_file.value(), 0);
    optimizer->writeLevel(hfile, 1);
  }

  LOG << "Done in " << FG_YELLOW(end - start) << " sec" << std::endl;
  return 0;
}

void generateConfig(const std::string& config_file_) {
  ConfigurableManager manager;

  auto robustifier = manager.create<RobustifierCauchy>();
  robustifier->param_chi_threshold.setValue(1.0);

  auto policy = manager.create<RobustifierPolicyByType>();
  policy->param_factor_class_name.setValue("SE2PosePoseGeodesicErrorFactor");
  policy->param_robustifier.setValue(robustifier);

  auto optimizer = manager.create<HierarchicalOptimizer>("optimizer");
  auto lm        = optimizer->param_low_level_alg.value();
  lm->setName("LM");
  auto dl = manager.create<IterationAlgorithmDL>("dogleg");
  optimizer->param_low_level_alg.setValue(dl);
  auto term_criteria =
    manager.create<RelativeGradientChiTerminationCriteria>("relative_term_criteria");
  auto& global_solver = optimizer->param_solver.value();
  global_solver->setName("global_solver");
  global_solver->param_max_iterations.pushBack(10);
  global_solver->param_max_iterations.pushBack(10);
  global_solver->param_robustifier_policies.pushBack(policy);
  global_solver->param_algorithm.setValue(dl);

  auto& partition_manager = optimizer->param_partition_manager.value();
  partition_manager->setName("partition_manager");
  // tg virtual factors pose pose 3D
  auto virtual_factor_se3_posepose =
    manager.create<SE3PosePoseVirtualFactorCreator>("se3_posepose_condended");
  partition_manager->param_covariance_projectors.pushBack(virtual_factor_se3_posepose);
  // tg virtual factors pose pose chordal
  auto virtual_factor_se3_chordal =
    manager.create<SE3PosePoseChordalVirtualFactorCreator>("se3_chordal_virtual");
  partition_manager->param_covariance_projectors.pushBack(virtual_factor_se3_chordal);

  auto& partition_solver = partition_manager->param_solver.value();
  partition_solver->param_max_iterations.pushBack(50);
  partition_solver->setName("partition_solver");
  partition_solver->param_algorithm.setValue(dl);

  manager.write(config_file_);
  LOG << "created default configuration in [" << FG_YELLOW(config_file_) << "]\n";
}
