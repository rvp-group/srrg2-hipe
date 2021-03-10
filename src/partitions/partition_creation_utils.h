#pragma once
#include <srrg_config/configurable.h>
#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  class VirtualFactorCreatorBase : public Configurable {
  public:
    virtual FactorBasePtr
    compute(VariableBase* variable_, VariableBase* anchor, MatrixBlockBase* covariance_) = 0;
  };

  using VirtualFactorCreatorBasePtr = std::shared_ptr<VirtualFactorCreatorBase>;

  class SE3PosePoseVirtualFactorCreator : public VirtualFactorCreatorBase {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    FactorBasePtr compute(VariableBase* variable_,
                          VariableBase* anchor_,
                          MatrixBlockBase* covariance_matrix_) override;
  };

  class SE3PosePoseChordalVirtualFactorCreator : public VirtualFactorCreatorBase {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    FactorBasePtr compute(VariableBase* variable_,
                          VariableBase* anchor_,
                          MatrixBlockBase* covariance_matrix_) override;
  };
} // namespace srrg2_hipe
