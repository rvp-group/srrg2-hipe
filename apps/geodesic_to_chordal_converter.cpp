#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_system_utils/parse_command_line.h>
#include <srrg_system_utils/shell_colors.h>

const std::string exe_name("geodesic_to_chordal_converter");
#define LOG std::cerr << FG_YELLOW(exe_name) << "|"
using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;

using VariableInType  = VariableSE3QuaternionRight;
using VariableOutType = VariableSE3EulerLeft;
using FactorInType    = SE3PosePoseGeodesicErrorFactor;
using FactorOutType   = SE3PosePoseChordalHessianFactor;

void convertGraph(FactorGraphPtr& out, const FactorGraphPtr& input);

int main(int argc, char** argv) {
  variables_and_factors_3d_registerTypes();
  ParseCommandLine cmd_line(argv);
  ArgumentString input_file(&cmd_line, "i", "input-file", "file containing graph to convert", "");
  ArgumentString output_file(&cmd_line, "o", "output-file", "file where to save the output ", "");
  cmd_line.parse();
  if (!input_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, no input file specified");
  }
  if (!output_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, no output file specified");
  }
  const std::string& file = input_file.value();
  FactorGraphPtr graph    = FactorGraph::read(file);
  if (!graph) {
    throw std::runtime_error(exe_name + "|ERROR, invalid graph file [ " + file + " ]");
  }
  LOG << "loaded [ " << FG_YELLOW(file) << " ] ----- "
      << "variables [ " << FG_YELLOW(graph->variables().size()) << " ] | "
      << "factors [ " << FG_YELLOW(graph->factors().size()) << " ]" << std::endl;
  FactorGraphPtr chordal_graph(new FactorGraph);
  convertGraph(chordal_graph, graph);
  chordal_graph->write(output_file.value());
  LOG << "output [ " FG_YELLOW(output_file.value()) " ] ----- "
      << "variables [ " << FG_YELLOW(graph->variables().size()) << " ] | "
      << "factors [ " << FG_YELLOW(graph->factors().size()) << " ]" << std::endl;
  return 0;
}

void convertGraph(FactorGraphPtr& out, const FactorGraphPtr& input) {
  out->clear();
  const auto& factors = input->factors();
  for (const auto& id_factor : factors) {
    const FactorBase* f = id_factor.second;
    const FactorInType* f_in = dynamic_cast<const FactorInType*>(f);
    if (f_in) {
      FactorOutType* f_out     = new FactorOutType;
      VariableBase::Id id_from = f_in->variableId(0);
      VariableBase::Id id_to   = f_in->variableId(1);
      f_out->setVariableId(0, id_from);
      f_out->setVariableId(1, id_to);
      f_out->setMeasurement(f_in->measurement());
      f_out->setInformationMatrix(f_in->informationMatrix());
      f_out->setEnabled(f_in->enabled());
      f_out->setGraphId(f_in->graphId());
      out->addFactor(FactorBasePtr(f_out));

      const VariableInType* v_from = dynamic_cast<VariableInType*>(input->variable(id_from));
      if (!v_from) {
        throw std::runtime_error("convertGraph| type misamatch in variables inside factor");
      }
      VariableOutType* v_from_out = new VariableOutType;
      v_from_out->setGraphId(v_from->graphId());
      v_from_out->setEstimate(v_from->estimate());
      v_from_out->setStatus(v_from->status());
      out->addVariable(VariableBasePtr(v_from_out));

      const VariableInType* v_to = dynamic_cast<VariableInType*>(input->variable(id_to));
      if (!v_to) {
        throw std::runtime_error("convertGraph| type misamatch in variables inside factor");
      }
      VariableOutType* v_to_out = new VariableOutType;
      v_to_out->setGraphId(v_to->graphId());
      v_to_out->setEstimate(v_to->estimate());
      v_to_out->setStatus(v_to->status());
      out->addVariable(VariableBasePtr(v_to_out));
    }
  }
}

