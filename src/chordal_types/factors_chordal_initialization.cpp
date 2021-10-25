#include "factors_chordal_initialization.h"
#include <srrg_solver/solver_core/error_factor_impl.cpp>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;
  using Matrix9f = Eigen::Matrix<float, 9, 9>;
  using JacobianType = Eigen::Matrix<float, 9, 3>;

  void SE3ChordalInitializationRotationErrorFactor::errorAndJacobian(bool error_only) {
    const Matrix3f M = _measurement.linear().transpose();

    const auto& from = _variables.at<0>()->estimate();
    const auto& to   = _variables.at<1>()->estimate();
    _e.setZero();
    _e.block<3, 1>(0, 0) = M * from.linear().row(0).transpose() - to.linear().row(0).transpose();
    _e.block<3, 1>(3, 0) = M * from.linear().row(1).transpose() - to.linear().row(1).transpose();
    _e.block<3, 1>(6, 0) = M * from.linear().row(2).transpose() - to.linear().row(2).transpose();

    Matrix9f J_from = Matrix9f::Zero();
    J_from.block<3, 3>(0, 0) = M;
    J_from.block<3, 3>(3, 3) = M;
    J_from.block<3, 3>(6, 6) = M;
    jacobian<0>()            = J_from;
    jacobian<1>() = -1.f * Matrix9f::Identity();
  }

  void SE3ChordalInitializationRotationErrorFactor::setInformationMatrix(const Matrix6f& omega) {
    Matrix3f sigma        = omega.block<3, 3>(3, 3).inverse();
    float sum_eigenvalues = sigma.trace();
    float kappa           = 3 / (2 * sum_eigenvalues);
    _information_matrix   = kappa * Matrix9f::Identity();
  }

  void SE3ChordalInitializationTranslationErrorFactor::errorAndJacobian(bool error_only) {
    const Vector3f& t_z = _measurement.translation();

    const auto& from = _variables.at<0>()->estimate();
    const auto& to   = _variables.at<1>()->estimate();
    const Matrix3f& Ri = from.linear();
    const Vector3f& tj = to.translation();
    const Vector3f& ti = from.translation();
    _e.setZero();
    _e = Ri.transpose() * (tj - ti) - t_z;

    jacobian<0>() = -1.f * Ri.transpose();
    jacobian<1>() = Ri.transpose();
  }
  void SE3ChordalInitializationTranslationErrorFactor::setInformationMatrix(const Matrix6f& omega) {
    _information_matrix = omega.block<3, 3>(0, 0);
  }
} // namespace srrg2_hipe
