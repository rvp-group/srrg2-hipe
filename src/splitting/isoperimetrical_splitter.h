#pragma once
#include "isoperimetrical_bipartitioner.h"
#include <deque>
#include <partitions/partition_manager.h>
#include <srrg_config/configurable.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  class IsoperimetricSplitter : public Configurable {
  public:
    PARAM(PropertyInt,
          max_partition_size,
          "Maximum number of variables in a partition",
          50,
          nullptr);

    IsoperimetricSplitter();
    inline void setPartitionManager(PartitionManager* manager_) {
      _manager = manager_;
      assert(manager_->graph() &&
             "IsoperimetricSplitter::setPartitionManager| No graph in manager");
      _graph = manager_->graph();
    }
    
    void compute();

  protected:
    using PartitionPtrDeque = std::deque<PartitionPtr>;
    PartitionPtrDeque _partitions_deque;
    FactorGraphInterfacePtr _graph = nullptr;
    PartitionManager* _manager     = nullptr;
    std::unique_ptr<IsoperimetricBipartitioner> _bipartitioner;
    // @brief divide a partition in two
    void _split(PartitionPtr& partition_);
  };
}
