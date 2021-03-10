#include "partition_creation_utils.h"

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  FactorBasePtr SE3PosePoseVirtualFactorCreator::compute(VariableBase* variable_,
                                                         VariableBase* anchor_,
                                                         MatrixBlockBase* covariance_matrix_) {
    std::shared_ptr<SE3PosePoseGeodesicErrorFactor> f(new SE3PosePoseGeodesicErrorFactor);
    VariableSE3QuaternionRight* anchor = dynamic_cast<VariableSE3QuaternionRight*>(anchor_);
    if (!anchor) {
      return nullptr;
    }
    VariableSE3QuaternionRight* variable = dynamic_cast<VariableSE3QuaternionRight*>(variable_);
    if (!variable) {
      return nullptr;
    }
    // Anchor is always assumed to be the first variable
    f->setVariableId(0, anchor->graphId());
    f->setVariableId(1, variable->graphId());
    // required to connect the variables inside the factor
    Isometry3f z = anchor->estimate().inverse() * variable->estimate();
    f->setMeasurement(z);
    // The covariance of the virtual factor equals the covariance of the distribution of the
    // boundary variable conditioned on the anchor ---> p(X_b|X_a)
    Matrix6f sigma_xx = covariance_matrix_->eigenType<Matrix6f>();
    Eigen::CompleteOrthogonalDecomposition<Matrix6f> dec(sigma_xx);
    Matrix6f omega = dec.pseudoInverse();
    f->setInformationMatrix(omega);
    return f;
  }

  FactorBasePtr
  SE3PosePoseChordalVirtualFactorCreator::compute(VariableBase* variable_,
                                                  VariableBase* anchor_,
                                                  MatrixBlockBase* covariance_matrix_) {
    std::shared_ptr<SE3PosePoseChordalHessianFactor> f(new SE3PosePoseChordalHessianFactor);
    VariableSE3EulerLeft* anchor = dynamic_cast<VariableSE3EulerLeft*>(anchor_);
    if (!anchor) {
      return nullptr;
    }
    VariableSE3EulerLeft* variable = dynamic_cast<VariableSE3EulerLeft*>(variable_);
    if (!variable) {
      return nullptr;
    }
    // Anchor is always assumed to be the first variable
    f->setVariableId(0, anchor->graphId());
    f->setVariableId(1, variable->graphId());
    // required to connect the variables inside the factor
    Isometry3f z = anchor->estimate().inverse() * variable->estimate();
    f->setMeasurement(z);
    // The covariance of the virtual factor equals the covariance of the distribution of the
    // boundary variable conditioned on the anchor ---> p(X_b|X_a)
    Matrix6f sigma_xx = covariance_matrix_->eigenType<Matrix6f>();
    Eigen::CompleteOrthogonalDecomposition<Matrix6f> dec(sigma_xx);
    Matrix6f omega = dec.pseudoInverse();
    f->setInformationMatrix(omega);
    return f;
  }
} // namespace srrg2_hipe
