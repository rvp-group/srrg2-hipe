#pragma once
#include <partitions/partition_manager.h>
#include <srrg_config/configurable.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  class SplitterBase : public Configurable {
  public:
    PARAM(PropertyInt,
          min_partition_variables,
          "Minimum number of variables in a partition",
          100,
          nullptr);
    inline void setPartitionManager(PartitionManager* manager_) {
      _manager = manager_;
      assert(manager_->graph() && "SplitterBase::setPartitionManager| No graph in manager");
      _graph = manager_->graph();
    }

    inline VariableBase::Id conditioner() const {
      return _conditioner;
    }

    virtual void compute() = 0;

  protected:
    PartitionManager* _manager     = nullptr;
    FactorGraphInterfacePtr _graph = nullptr;
    VariableBase::Id _conditioner  = -1;
  };

  using SplitterBasePtr = std::shared_ptr<SplitterBase>;
} // namespace srrg2_hipe
