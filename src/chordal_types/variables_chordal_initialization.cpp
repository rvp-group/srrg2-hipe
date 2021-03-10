#include "variables_chordal_initialization.h"
#include "srrg_solver/solver_core/variable_impl.cpp"
#include <Eigen/SVD>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  void VariableSE3ChordalInitializationRotation::applyPerturbation(const Vector9f& perturbation) {
    Matrix3f M = Matrix3f::Identity();
    M.row(0)   = perturbation.block<3, 1>(0, 0).transpose();
    M.row(1)   = perturbation.block<3, 1>(3, 0).transpose();
    M.row(2)   = perturbation.block<3, 1>(6, 0).transpose();
    _estimate.linear() += M;
    const Matrix3f& R = _estimate.linear();
    Eigen::JacobiSVD<Matrix3f> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    float det_uv = (svd.matrixU() * svd.matrixV().transpose()).determinant();
    Vector3f singular_values(1.f, 1.f, det_uv);
    _estimate.linear() = svd.matrixU() * singular_values.asDiagonal() * svd.matrixV().transpose();
  }

  void
  VariableSE3ChordalInitializationTranslation::applyPerturbation(const Vector3f& perturbation) {
    _estimate.translation() += perturbation;
  }
} // namespace srrg2_hipe
