#pragma once
#include "variables_chordal_initialization.h"
#include <srrg_solver/solver_core/error_factor.h>
#include <srrg_solver/solver_core/measurement_owner.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  class SE3ChordalInitializationRotationErrorFactor
    : public ErrorFactor_<9,
                          VariableSE3ChordalInitializationRotation,
                          VariableSE3ChordalInitializationRotation>,
      public MeasurementOwnerEigen_<Isometry3f> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void errorAndJacobian(bool error_only) final;
    void setInformationMatrix(const Matrix6f& omega);
  };

  class SE3ChordalInitializationTranslationErrorFactor
    : public ErrorFactor_<3,
                          VariableSE3ChordalInitializationTranslation,
                          VariableSE3ChordalInitializationTranslation>,
      public MeasurementOwnerEigen_<Isometry3f> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void errorAndJacobian(bool error_only) final;
    void setInformationMatrix(const Matrix6f& omega);
  };
} // namespace srrg2_hipe
