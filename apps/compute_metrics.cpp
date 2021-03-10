#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/solver_core/solver.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_system_utils/parse_command_line.h>
#include <srrg_system_utils/shell_colors.h>

const std::string exe_name("compute_metrics");
#define LOG std::cerr << FG_YELLOW(exe_name) << "|"
#define NLOG std::cerr << "\n" << FG_YELLOW(exe_name) << "|"
using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;
using PointVector = std::vector<Vector3f, Eigen::aligned_allocator<Vector3f>>;

void computeAlignmentTransform(const FactorGraphPtr& graph,
                               const FactorGraphPtr& target_graph,
                               Isometry3f& delta);
void align(FactorGraphPtr& graph, const Isometry3f& T);

void evaluateAte(FactorGraphPtr& graph,
                 const FactorGraphPtr& gt,
                 float& ate_rotation,
                 float& ate_translation,
                 float& rmse_translation);

void evaluateRE(const FactorGraphPtr& graph,
                const FactorGraphPtr& gt,
                float& re_rotation,
                float& re_translation,
                float& avg_translation,
                float& avg_absolute_rotation);

float getAngle(const Matrix3f& R);

int main(int argc, char** argv) {
  variables_and_factors_3d_registerTypes();
  ParseCommandLine cmd_line(argv);
  ArgumentString input_file(
    &cmd_line, "i", "input-file", "file containing the query graph (.boss)", "");
  ArgumentString gt_file(
    &cmd_line, "gt", "ground-truth", "file containing the ground truth (.boss)", "");
  cmd_line.parse();
  if (!input_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, no input file specified");
  }
  if (!gt_file.isSet()) {
    std::cerr << cmd_line.options() << std::endl;
    throw std::runtime_error(exe_name + "|ERROR, no output file specified");
  }
  const std::string& file = input_file.value();
  FactorGraphPtr graph    = FactorGraph::read(file);
  if (!graph) {
    throw std::runtime_error(exe_name + "|ERROR, invalid graph file [ " + file + " ]");
  }
  const std::string& gt = gt_file.value();
  FactorGraphPtr graph_gt = FactorGraph::read(gt);
  if (!graph_gt) {
    throw std::runtime_error(exe_name + "|ERROR, invalid graph file [ " + gt + " ]");
  }

  Solver solver;
  solver.param_max_iterations.pushBack(1);
  // tg compute normalized chi2 for the query
  solver.setGraph(graph);
  float chi_normalized_query = solver.normalizedChi2();
  // tg compute normalized chi2 for the ground truth
  solver.setGraph(graph_gt);
  float chi_normalized_gt = solver.normalizedChi2();

  float ate_rotation    = 0.f;
  float ate_translation = 0.f;
  float rmse_translation = 0.f;
  float re_rotation      = 0.f;
  float re_translation   = 0.f;
  float avg_translation       = 0.f;
  float avg_absolute_rotation = 0.f;
  evaluateAte(graph, graph_gt, ate_rotation, ate_translation, rmse_translation);
  evaluateRE(graph, graph_gt, re_rotation, re_translation, avg_translation, avg_absolute_rotation);

  NLOG << "chi normalized query graph : " << FG_YELLOW(chi_normalized_query) << std::endl;
  LOG << "chi normalized ground truth graph : " << FG_YELLOW(chi_normalized_gt) << std::endl;
  LOG << "ratio chi_gt/chi_query : " << FG_YELLOW(chi_normalized_gt / chi_normalized_query)
      << std::endl;
  NLOG << "ATE rotation [RMSE] : " << FG_GREEN(ate_rotation) << std::endl;
  LOG << "ATE translation [RMSE] : " << FG_GREEN(ate_translation) << std::endl;
  NLOG << "Absolute translation [RMSE]: " << FG_BWHITE(rmse_translation) << std::endl;
  NLOG << "RE rotation [RMSE] : " << FG_GREEN(re_rotation) << std::endl;
  LOG << "RE translation [RMSE] : " << FG_GREEN(re_translation) << std::endl;
  NLOG << "Avg distance travelled : " << FG_YELLOW(avg_translation) << std::endl;
  LOG << "Avg absolute orientation : " << FG_YELLOW(avg_absolute_rotation) << std::endl;
  return 0;
}

float getAngle(const Matrix3f& R) {
  // to compute the angular error we refer to the axis-angle representation of the error
  // rotation matrix from that we can extract the sin and cos of the angle and then use the
  // atan2 (see Alessandro De Luca slides for the Robotics 1 course)
  const Matrix3f skew_R = R - R.transpose();
  // sin of the error angle
  float sin_angle = 0.5f * std::sqrt(skew_R(0, 1) * skew_R(0, 1) + skew_R(0, 2) * skew_R(0, 2) +
                                     skew_R(1, 2) * skew_R(1, 2));
  // cos of the error angle
  float cos_angle = 0.5f * (R.trace() - 1.f);

  float angle = std::atan2(sin_angle, cos_angle);
  return angle;
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
  Matrix3f R = svd.matrixU() * W * svd.matrixV().transpose();
  fixRotation(R);
  delta.translation() = mean_target - R * mean_query;
  delta.linear()      = R;
}

void align(FactorGraphPtr& graph, const Isometry3f& T) {
  for (const auto& id_var : graph->variables()) {
    VariableSE3Base* v          = dynamic_cast<VariableSE3Base*>(id_var.second);
    Isometry3f estimate         = v->estimate();
    Isometry3f aligned_estimate = T * estimate;
    v->setEstimate(aligned_estimate);
  }
}

void evaluateAte(FactorGraphPtr& graph,
                 const FactorGraphPtr& gt,
                 float& ate_rotation,
                 float& ate_translation,
                 float& rmse_translation) {
  Isometry3f delta_alignment = Isometry3f::Identity();
  computeAlignmentTransform(graph, gt, delta_alignment);
  align(graph, delta_alignment);
  auto& gt_variables    = gt->variables();
  auto& query_variables = graph->variables();
  // tg check it the two graph are same number of variables
  if (gt_variables.size() != query_variables.size()) {
    throw std::runtime_error(exe_name +
                             "|number of variables mismatch between query and ground truth");
  }
  // tg counter for skipped variables (due to nans) and processed
  size_t skipped   = 0;
  size_t processed = 0;
  // tg for each ground truth variable
  for (const auto& id_var : gt_variables) {
    VariableSE3Base* gt_v = dynamic_cast<VariableSE3Base*>(id_var.second);
    auto it_query            = query_variables.find(id_var.first);
    VariableSE3Base* query_v = dynamic_cast<VariableSE3Base*>(it_query.value());

    if (gt_v == nullptr) {
      ++skipped;
      continue;
    }

    if (query_v == nullptr) {
      ++skipped;
      continue;
    }
    // extract rotation and translation for both variables
    const Matrix3f Ri          = gt_v->estimate().linear();
    const Matrix3f Rj          = query_v->estimate().linear();
    const Vector3f ti          = gt_v->estimate().translation();
    const Vector3f tj          = query_v->estimate().translation();
    // compute rotation error
    Matrix3f R_error = Ri * Rj.transpose();
    // fix round off due to the multiplication
    fixRotation(R_error);
    // compute rotation error
    float angle_error    = getAngle(R_error);
    Vector3f t_error_ATE = ti - R_error * tj;
    Vector3f t_error     = ti - tj;
    // check for nans
    if (std::isnan(t_error_ATE.norm()) || std::isnan(angle_error)) {
      LOG << "WARNING, detected nan, skipping vertex [ " << id_var.first << " ]\n";
      ++skipped;
      continue;
    }
    // accumulate squared errors
    ate_translation += t_error_ATE.transpose() * t_error_ATE;
    ate_rotation += angle_error * angle_error;
    rmse_translation += t_error.transpose() * t_error;
    ++processed;
  }

  LOG << "processed variables [ " << processed << "/" << gt_variables.size() << " ]\n";
  LOG << "skipped variables [ " << skipped << "/" << gt_variables.size() << " ]\n";
  if (!processed) {
    throw std::runtime_error(exe_name + "|invalid types");
  }
  // get number of variables processed
  const float inverse_num_variables = 1.f / (float) (processed);
  // compute final ATE
  ate_translation *= inverse_num_variables;
  ate_rotation *= inverse_num_variables;
  rmse_translation *= inverse_num_variables;

  ate_translation    = std::sqrt(ate_translation);
  ate_rotation       = std::sqrt(ate_rotation);
  rmse_translation   = std::sqrt(rmse_translation);
}

void evaluateRE(const FactorGraphPtr& graph,
                const FactorGraphPtr& gt,
                float& re_rotation,
                float& re_translation,
                float& avg_translation,
                float& avg_absolute_rotation) {
  const auto& factors_gt = gt->factors();
  size_t skipped         = 0;
  size_t processed       = 0;
  // tg for each ground truth variable
  for (const auto& id_factor : factors_gt) {
    FactorBase* f = id_factor.second;
    VariableBase::Id id_from = f->variableId(0);
    VariableBase::Id id_to   = f->variableId(1);

    // tg compute relative pose for ground truth
    VariableSE3Base* from = dynamic_cast<VariableSE3Base*>(gt->variable(id_from));
    VariableSE3Base* to   = dynamic_cast<VariableSE3Base*>(gt->variable(id_to));
    Isometry3f T_gt = from->estimate().inverse() * to->estimate();
    // tg compute relative pose estimated
    from                  = dynamic_cast<VariableSE3Base*>(graph->variable(id_from));
    to                    = dynamic_cast<VariableSE3Base*>(graph->variable(id_to));
    Isometry3f T_estimate = from->estimate().inverse() * to->estimate();

    const Matrix3f& R_gt = T_gt.linear();
    const Vector3f& t_gt       = T_gt.translation();
    const float abs_angle_gt   = std::abs(getAngle(R_gt));
    const Matrix3f& R_estimate = T_estimate.linear();
    const Vector3f& t_estimate = T_estimate.translation();
    Matrix3f R_error           = R_gt * R_estimate.transpose();
    // fix round off due to the multiplication
    fixRotation(R_error);
    // compute errors
    float angle_error = getAngle(R_error);
    Vector3f t_error = t_gt - R_error * t_estimate;
    // check for nans
    if (std::isnan(t_error.norm()) || std::isnan(angle_error)) {
      LOG << "WARNING, detected nan, skipping edge [ " << id_factor.first << " ]\n";
      ++skipped;
      continue;
    }
    re_rotation += angle_error * angle_error;
    re_translation += t_error.transpose() * t_error;
    avg_translation += t_gt.norm();
    avg_absolute_rotation += abs_angle_gt;
    ++processed;
  }
  LOG << "processed factors [ " << processed << "/" << factors_gt.size() << " ]\n";
  LOG << "skipped factors [ " << skipped << "/" << factors_gt.size() << " ]\n";
  if (!processed) {
    throw std::runtime_error(exe_name + "|invalid types");
  }
  // get number of variables processed
  const float inverse_num_factors = 1.f / (float) (processed);
  // compute final RPE
  re_translation *= inverse_num_factors;
  re_rotation *= inverse_num_factors;
  avg_translation *= inverse_num_factors;
  avg_absolute_rotation *= inverse_num_factors;

  re_translation = std::sqrt(re_translation);
  re_rotation    = std::sqrt(re_rotation);
}
