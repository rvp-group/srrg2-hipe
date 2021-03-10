#include <chordal_types/chordal_initializer.h>
#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/solver_core/solver.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_system_utils/parse_command_line.h>
#include <srrg_system_utils/system_utils.h>

const std::string executable_name("chordal_initialization_example");
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
  Solver solver;
  solver.param_max_iterations.pushBack(1);
  solver.setGraph(graph);
  solver.compute();
  LOG << solver.iterationStats() << std::endl;
  std::set<VariableBase::Id> pose_variables;
  double start = getTime();
  initializeChordal(graph, &pose_variables);
  double end = getTime();
  solver.setGraph(graph);
  solver.compute();
  LOG << solver.iterationStats() << std::endl;
  LOG << "Done in " << FG_YELLOW(end - start) << " sec" << std::endl;
  return 0;
}
