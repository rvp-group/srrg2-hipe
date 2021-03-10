#pragma once
#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>

namespace srrg2_hipe {
  using namespace srrg2_solver;
  using VariableIdSet = std::set<VariableBase::Id>;

  class OdometryPropagation {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    bool compute(const VariableBase* origin, const FactorBase* edge, VariableIdSet& variable_set_);

  protected:
    Isometry3f _odom_measurement = Isometry3f::Identity();

    bool propagateTo(const VariableBase* from, FactorBase* f, VariableIdSet& variable_set_);
    bool propagateFrom(const VariableBase* to, FactorBase* f, VariableIdSet& variable_set_);
    bool getMeasurement(FactorBase* factor);
  };
  using OdometryPropagationPtr = std::unique_ptr<OdometryPropagation>;
}
