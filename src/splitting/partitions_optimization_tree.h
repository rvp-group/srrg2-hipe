#pragma once
#include "odometry_propagation.h"
#include <partitions/partition_manager.h>

namespace srrg2_hipe {
  using namespace srrg2_core;

  class PartitionsOptimizationTree {
  public:
    PartitionsOptimizationTree();
    inline void setPartitionManager(PartitionManager* manager_) {
      _manager = manager_;
    }
    void compute();

  protected:
    void _optimizeAndExpand(PartitionPtr& partition_);
    void _initPartition(PartitionPtr& partition_, const VariableBase::Id& origin);
    VariableIdSet _visited_partitions;
    VariableIdSet _visited_variables;
    OdometryPropagationPtr _init_action;
    std::deque<PartitionPtr> _partitions_to_be_expanded;
    PartitionManager* _manager = nullptr;
  };
} // namespace srrg2_hipe
