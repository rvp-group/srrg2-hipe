#include "partition.h"

namespace srrg2_hipe {

  Partition::Partition(const size_t& level_) {
    _level = level_;
  }

  Partition::Partition(const VariablePtrVector& variables_,
                       const FactorPtrVector& low_level_factors_,
                       const size_t& level_) {
    for (VariableBase* v : variables_) {
      addVariable(v);
    }
    _level          = level_;
    int lower_level = level_ - 1;
    if (!std::all_of(low_level_factors_.begin(),
                     low_level_factors_.end(),
                     [lower_level](FactorBase* f) -> bool { return f->level() == lower_level; })) {
      throw std::runtime_error(
        "Partition::Partition| Not all the factors belong to the lower level of the partition");
    }
    for (FactorBase* f : low_level_factors_) {
      addGraphFactor(f);
    }
  }

  bool Partition::addVariable(VariableBase* variable_) {
    assert(variable_ && "Partition::addVariable| variable pointer is null");
    auto id_variable_pair = std::make_pair(variable_->graphId(), variable_);
    return _variables.insert(id_variable_pair);
  }

  bool Partition::addGraphFactor(FactorBase* factor_) {
    assert(factor_ && "Partition::addLowLevelFactor| factor_ pointer is null");
    auto id_factor_pair = std::make_pair(factor_->graphId(), factor_);
    return _factors.insert(id_factor_pair);
  }

  bool Partition::addVirtualFactor(FactorBase* factor_) {
    assert(factor_ && "Partition::addVirtualFactor| factor_ pointer is null");
    bool correctly_added = addGraphFactor(factor_);
    if (correctly_added) {
      FactorBase* f = const_cast<FactorBase*>(factor_);
      assert(f && "Partition::addVirtualFactor| cannot const_cast factor pointer");
      f->setLevel(_level);
      _virtual_factors.push_back(f);
    }
    return correctly_added;
  }

  IdVariablePtrContainer& Partition::variables() {
    return _variables;
  }

  void Partition::addBoundaryVariableIds(const VariableIdSet& ids_) {
    for (const VariableBase::Id& v_id : ids_) {
      assert(variable(v_id) &&
             "Partition::addBoundaryVariableIds| variable is not in the partition");
      _boundary_variable_ids.insert(v_id);
    }
  }

  VariableIdSet& Partition::boundaryVariableIds() {
    return _boundary_variable_ids;
  }

  const VariableIdSet& Partition::boundaryVariableIds() const {
    return _boundary_variable_ids;
  }

  FactorPtrVector& Partition::virtualFactors() {
    return _virtual_factors;
  }

  const FactorPtrVector& Partition::virtualFactors() const {
    return _virtual_factors;
  }

  IdFactorPtrContainer& Partition::factors() {
    return _factors;
  }

} // namespace srrg2_hipe
