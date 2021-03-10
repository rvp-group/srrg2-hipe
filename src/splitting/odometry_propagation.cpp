#include "odometry_propagation.h"

namespace srrg2_hipe {
  using namespace srrg2_solver;

  bool OdometryPropagation::compute(const VariableBase* origin,
                                    const FactorBase* edge,
                                    VariableIdSet& variable_set_) {
    FactorBase* factor = const_cast<FactorBase*>(edge);
    assert(factor && "OdometryPropagation::compute| Cannot const_cast factor");
    if (!getMeasurement(factor)) {
      return false;
    }
    bool is_init = false;
    if (origin->graphId() == factor->variableId(0)) {
      is_init = propagateTo(origin, factor, variable_set_);
    } else {
      is_init = propagateFrom(origin, factor, variable_set_);
    }
    return is_init;
  }

  bool OdometryPropagation::propagateTo(const VariableBase* from,
                                        FactorBase* factor,
                                        VariableIdSet& variable_set_) {
    bool is_init  = false;
    const VariableSE3Base* v = dynamic_cast<const VariableSE3Base*>(from);
    if (!v) {
      return is_init;
    }
    VariableSE3Base* to = dynamic_cast<VariableSE3Base*>(factor->variable(1));
    if (to && !variable_set_.count(to->graphId())) {
      to->setEstimate(v->estimate() * _odom_measurement);
      is_init = true;
    }
    return is_init;
  }

  bool OdometryPropagation::propagateFrom(const VariableBase* to,
                                          FactorBase* factor,
                                          VariableIdSet& variable_set_) {
    bool is_init  = false;
    const VariableSE3Base* v = dynamic_cast<const VariableSE3Base*>(to);
    if (!v) {
      return is_init;
    }
    VariableSE3Base* from = dynamic_cast<VariableSE3Base*>(factor->variable(0));
    if (from && !variable_set_.count(from->graphId())) {
      from->setEstimate(v->estimate() * _odom_measurement.inverse());
      is_init = true;
    }
    return is_init;
  }

  bool OdometryPropagation::getMeasurement(FactorBase* factor) {
    SE3PosePoseGeodesicErrorFactor* geodesic =
      dynamic_cast<SE3PosePoseGeodesicErrorFactor*>(factor);
    if (geodesic) {
      _odom_measurement = geodesic->measurement();
      return true;
    }
    SE3PosePoseChordalHessianFactor* chordal = dynamic_cast<SE3PosePoseChordalHessianFactor*>(factor);
    if (chordal) {
      _odom_measurement = chordal->measurement();
      return true;
    }
    return false;
  }
} // namespace srrg2_hipe
