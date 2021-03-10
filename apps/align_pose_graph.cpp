#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/solver_core/solver.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_system_utils/parse_command_line.h>
#include <srrg_system_utils/shell_colors.h>

const std::string exe_name("align_pose_graph");
#define LOG std::cerr << FG_YELLOW(exe_name) << "|"
using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;

using FactorType = SE3PriorErrorFactorAD;
using VariableType = VariableSE3QuaternionRightAD;
using PointVector  = std::vector<Vector3f, Eigen::aligned_allocator<Vector3f>>;

void computeAlignmentTransform(const FactorGraphPtr& graph,
                               const FactorGraphPtr& target_graph,
                               Isometry3f& delta);
void align(FactorGraphPtr& graph, const Isometry3f& T);

int main(int argc, char** argv) {
  variables_and_factors_3d_registerTypes();
  ParseCommandLine cmd_line(argv);
  ArgumentString input_file(
    &cmd_line, "i", "input-file", "file containing the query graph (.boss)", "");
  ArgumentString out_file(
    &cmd_line, "o", "out-file", "file containing the output graph (.boss)", "");
  ArgumentString gt_file(
    &cmd_line, "gt", "ground-truth", "file containing the ground truth (.boss)", "");
  ArgumentInt id_gauge(
    &cmd_line, "id", "id-gauge", "id of the fixed variable in the ground truth", -1);
  cmd_line.parse();
  if (!input_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, no input file specified");
  }
  if (!out_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, no output file specified");
  }
  if (!gt_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, no output file specified");
  }
  if (!id_gauge.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, id gauge not specified");
  }
  const std::string& file = input_file.value();
  FactorGraphPtr graph    = FactorGraph::read(file);
  if (!graph) {
    throw std::runtime_error(exe_name + "|ERROR, invalid graph file [ " + file + " ]");
  }
  const std::string& gt   = gt_file.value();
  FactorGraphPtr graph_gt = FactorGraph::read(gt);
  if (!graph_gt) {
    throw std::runtime_error(exe_name + "|ERROR, invalid graph file [ " + gt + " ]");
  }

  if (graph->variables().size() != graph_gt->variables().size()) {
    throw std::runtime_error(exe_name +
                             "|number of variables mismatch between query and ground truth");
  }

  VariableSE3Base* gauge = dynamic_cast<VariableSE3Base*>(graph_gt->variable(id_gauge.value()));
  VariableSE3Base* origin_query = dynamic_cast<VariableSE3Base*>(graph->variable(id_gauge.value()));
  if (!gauge || !origin_query) {
    throw std::runtime_error(exe_name + "|ERROR, gauge id is not present");
  }
  LOG << "Compute delta" << std::endl;
  Isometry3f delta = gauge->estimate() * origin_query->estimate().inverse();
  Matrix3f R       = delta.linear();
  fixRotation(R);
  delta.linear() = R;
  computeAlignmentTransform(graph, graph_gt, delta);
  LOG << "Align ...." << std::endl;
  align(graph, delta);
  graph->write(out_file.value());
  LOG << "Done" << std::endl;
  return 0;
}

void computeAlignmentTransform(const FactorGraphPtr& graph,
                               const FactorGraphPtr& target_graph,
                               Isometry3f& delta) {
  const auto& variables = target_graph->variables();
  PointVector target, query;
  target.reserve(variables.size());
  query.reserve(variables.size());
  for (const auto& id_var : variables) {
    const VariableBase::Id& id = id_var.first;
    VariableSE3Base* v_target  = dynamic_cast<VariableSE3Base*>(id_var.second);
    VariableSE3Base* v         = dynamic_cast<VariableSE3Base*>(graph->variable(id));
    query.emplace_back(v->estimate().translation());
    target.emplace_back(v_target->estimate().translation());
  }
  Vector3f mean_query = Vector3f::Zero(), mean_target = Vector3f::Zero();
  for (size_t i = 0; i < query.size(); ++i) {
    mean_query += query.at(i);
    mean_target += target.at(i);
  }
  const float inv_size = 1.f / (float) query.size();
  mean_query *= inv_size;
  mean_target *= inv_size;
  Matrix3f sigma = Matrix3f::Zero();
  for (size_t i = 0; i < query.size(); ++i) {
    sigma += (target.at(i) - mean_target) * (query.at(i) - mean_query).transpose();
  }
  sigma *= inv_size;
  Eigen::JacobiSVD<Matrix3f> svd(sigma, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Matrix3f W = Matrix3f::Identity();
  if (svd.matrixU().determinant() * svd.matrixV().determinant() < 0) {
    W(2, 2) = -1;
  }
  Matrix3f R          = svd.matrixU() * W * svd.matrixV().transpose();
  fixRotation(R);
  delta.translation() = mean_target - R * mean_query;
  delta.linear()      = R;
}

void align(FactorGraphPtr& graph, const Isometry3f& T) {
  for (const auto& id_var : graph->variables()) {
    VariableSE3Base* v     = dynamic_cast<VariableSE3Base*>(id_var.second);
    Isometry3f estimate    = v->estimate();
    Matrix3f R_es          = estimate.linear();
    estimate.linear()           = R_es;
    Isometry3f aligned_estimate = T * estimate;
    v->setEstimate(aligned_estimate);
  }
}
