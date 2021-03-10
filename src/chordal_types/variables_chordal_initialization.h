#pragma once
#include <srrg_geometry/geometry_defs.h>
#include <srrg_solver/solver_core/variable.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  using Vector9f = Eigen::Matrix<float, 9, 1>;

  class VariableSE3ChordalInitializationRotation : public Variable_<9, Isometry3_> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void setZero() final {
      _estimate.setIdentity();
    }

    void applyPerturbation(const Vector9f& perturbation) final;
  };

  class VariableSE3ChordalInitializationTranslation : public Variable_<3, Isometry3_> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void setZero() final {
      _estimate.setIdentity();
    }

    void applyPerturbation(const Vector3f& perturbation) final;
  };
} // namespace srrg2_hipe
