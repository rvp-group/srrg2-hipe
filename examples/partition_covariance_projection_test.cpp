#include <partitions/partition.h>
#include <partitions/partition_manager.h>
#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/solver_core/solver.h>
#include <srrg_system_utils/shell_colors.h>

using namespace srrg2_core;
using namespace srrg2_solver;
using namespace srrg2_hipe;

const std::string data(TEST_DATA_FOLDER);

int main(int argc, char** argv) {
  FactorGraphPtr graph = FactorGraph::read(data + "boss/sphere_bignoise_vertex3.boss");

  auto& variables = graph->variables();
  auto& factors   = graph->factors();

  VariablePtrVector partition_variables;
  FactorPtrVector partition_factors;
  VariableIdSet shared_variables;
  // First variable is a shared variable
  shared_variables.insert(0);
  // Gauge is the last variable
  VariableBase::Id anchor_id = 59;
  for (size_t v_id = 0; v_id < variables.size(); ++v_id) {
    partition_variables.push_back(graph->variable(v_id));
  }
  for (size_t f_id = 0; f_id < factors.size(); ++f_id) {
    partition_factors.push_back(graph->factor(f_id));
  }

  PartitionManager partition_maker;
  partition_maker.param_covariance_projectors.pushBack(
    VirtualFactorCreatorBasePtr(new SE3PosePoseVirtualFactorCreator));

  partition_maker.param_solver.value()->param_max_iterations.pushBack(10);
  partition_maker.param_solver.value()->param_termination_criteria.setValue(nullptr);
  partition_maker.setGraph(graph);
  partition_maker.createPartition(
    partition_variables, partition_factors, anchor_id, shared_variables);
  PartitionPtr partition     = partition_maker.partition(anchor_id);
  FactorBase* virtual_factor = partition->virtualFactors().front();
  const SE3PosePoseGeodesicErrorFactor* factor =
    dynamic_cast<const SE3PosePoseGeodesicErrorFactor*>(virtual_factor);

  std::cerr << FG_GREEN("Virtual factor error: \n") << factor->error() << std::endl;
  std::cerr << FG_GREEN("\nProjected information matrix : \n") << factor->informationMatrix()
            << std::endl;

  Eigen::JacobiSVD<Matrix6f> svd(factor->informationMatrix());
  std::cerr << FG_GREEN("Singular value to verify that is positive semi-definite : ")
            << svd.singularValues().transpose() << std::endl;
  if ((svd.singularValues().array() < 0).any()) {
    throw std::runtime_error("Possibly non semi-positive definite matrix!");
  }

  VariableBase* anchor = partition->anchor();
  anchor->setStatus(VariableBase::Status::Fixed);
  Solver solver;
  solver.param_max_iterations.pushBack(10);
  solver.setGraph(graph);
  solver.compute();

  auto& active_variable = solver.activeVariables();
  VariablePairVector v_pair;
  for (VariableBase* v : active_variable) {
    if (shared_variables.count(v->graphId())) {
      v_pair.push_back(std::make_pair(v, v));
    }
  }
  MatrixBlockVector cov;
  solver.computeMarginalCovariance(cov, v_pair);
  Matrix6f covariance_xx = cov[0]->eigenType<Matrix6f>();

  solver.setGraph(partition);
  solver.compute();

  auto& active_variable_partition = solver.activeVariables();
  v_pair.clear();
  for (VariableBase* v : active_variable_partition) {
    if (shared_variables.count(v->graphId())) {
      v_pair.push_back(std::make_pair(v, v));
    }
  }
  solver.computeMarginalCovariance(cov, v_pair);
  Matrix6f cov_virtual_xx = cov[0]->eigenType<Matrix6f>();
  std::cerr << FG_GREEN("Virtual covariance : \n") << cov_virtual_xx << std::endl;
  std::cerr << FG_GREEN("\nCovariance graph : \n") << covariance_xx << std::endl;
  float det_sigma_1   = covariance_xx.determinant();
  float det_sigma_0   = cov_virtual_xx.determinant();
  float kl_divergence = 0.5f * ((covariance_xx.inverse() * cov_virtual_xx).trace() - 6.f +
                                std::log(det_sigma_1 / det_sigma_0));
  std::cerr << FG_GREEN("\n KL divergence : ") << kl_divergence << std::endl;
  return 0;
}
