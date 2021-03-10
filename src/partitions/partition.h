#pragma once
#include <Eigen/Core>
#include <srrg_solver/solver_core/factor_graph.h>

namespace srrg2_hipe {
  using namespace srrg2_solver;
  using namespace srrg2_core;

  using VariablePtrVector = std::vector<VariableBase*>;
  using FactorPtrVector   = std::vector<FactorBase*>;
  using VariableIdSet     = std::set<VariableBase::Id>;

  class Partition : public FactorGraphView {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Partition(const size_t& level_ = 1);

    Partition(const VariablePtrVector& variables_,
              const FactorPtrVector& low_level_factors_,
              const size_t& level_);

    bool addVariable(VariableBase* variable_);
    // add factor that directly belong to the factor graph
    bool addGraphFactor(FactorBase* factor_);
    // add virtual factor at the level of which the partition is constructed
    bool addVirtualFactor(FactorBase* factor_);

    inline Id anchorId() const {
      return _anchor->graphId();
    }

    inline void setAnchor(VariableBase* anchor_) {
      _anchor = anchor_;
    }

    inline VariableBase* anchor() {
      return _anchor;
    }

    inline const VariableBase* anchor() const {
      return _anchor;
    }

    inline size_t level() {
      return _level;
    }

    void addBoundaryVariableIds(const VariableIdSet& ids_);

    VariableIdSet& boundaryVariableIds();
    const VariableIdSet& boundaryVariableIds() const;

    IdVariablePtrContainer& variables() override;

    FactorPtrVector& virtualFactors();
    const FactorPtrVector& virtualFactors() const;

    IdFactorPtrContainer& factors() override;

  protected:
    VariableBase* _anchor = nullptr;
    FactorPtrVector _virtual_factors;
    VariableIdSet _boundary_variable_ids;
    size_t _level;
  };

  using PartitionPtr       = std::shared_ptr<Partition>;
  using PartitionPtrVector = std::vector<PartitionPtr>;
} // namespace srrg2_hipe
