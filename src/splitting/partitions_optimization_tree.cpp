#include "partitions_optimization_tree.h"

namespace srrg2_hipe {
  using namespace srrg2_core;

  PartitionsOptimizationTree::PartitionsOptimizationTree() : _init_action(new OdometryPropagation) {
  }

  void PartitionsOptimizationTree::compute() {
    assert(_manager && "PartitionsOptimizationTree::compute| no partition manager");
    _visited_partitions.clear();
    _partitions_to_be_expanded.clear();
    auto partitions                       = _manager->partitions();
    PartitionPtr most_connected_partition = nullptr;
    size_t max_boundary_variables         = 0;
    for (const auto& id_anchor_partition : partitions) {
      PartitionPtr partition    = id_anchor_partition.second;
      size_t boundary_variables = partition->boundaryVariableIds().size();
      if (boundary_variables > max_boundary_variables) {
        most_connected_partition = partition;
        max_boundary_variables   = boundary_variables;
      }
    }
    _visited_partitions.insert(most_connected_partition->anchorId());
    _partitions_to_be_expanded.push_back(most_connected_partition);
    while (!_partitions_to_be_expanded.empty()) {
      PartitionPtr partition = _partitions_to_be_expanded.front();
      _optimizeAndExpand(partition);
      _partitions_to_be_expanded.pop_front();
    }
  }

  void PartitionsOptimizationTree::_optimizeAndExpand(PartitionPtr& partition_) {
    const VariableIdSet& boundary_variables = partition_->boundaryVariableIds();
    _manager->addBoundaryVariablesToPartition(partition_, boundary_variables);
    for (const VariableBase::Id& id : boundary_variables) {
      if (_visited_variables.count(id)) {
        continue;
      }
      const VariableBase* v               = partition_->variable(id);
      const VariableSE3Base* pose_type_3d = dynamic_cast<const VariableSE3Base*>(v);
      if (!pose_type_3d) {
        continue;
      }
      PartitionPtrVector connected_partitions;
      if (!_manager->getConnectedPartitions(connected_partitions, id)) {
        continue;
      }
      for (PartitionPtr s : connected_partitions) {
        if (_visited_partitions.count(s->anchorId())) {
          continue;
        }
        _initPartition(s, id);
        if (!_visited_partitions.count(s->anchorId())) {
          _partitions_to_be_expanded.push_back(s);
          _visited_partitions.insert(s->anchorId());
        }
      }
    }
  }

  void PartitionsOptimizationTree::_initPartition(PartitionPtr& partition_,
                                                  const VariableBase::Id& origin) {
    partition_->bindFactors();
    VariableBase* v = partition_->variable(origin);
    std::deque<VariableBase*> variable_deque;
    variable_deque.push_back(v);
    _visited_variables.insert(v->graphId());
    while (!variable_deque.empty()) {
      v               = variable_deque.front();
      auto lower      = partition_->lowerFactor(v);
      auto upper      = partition_->upperFactor(v);
      for (auto it = lower; it != upper; ++it) {
        FactorBase* f = it->second;
        if (f->numVariables() < 2) {
          continue;
        }
        VariableBase::Id other_id =
          f->variableId(0) == v->graphId() ? f->variableId(1) : f->variableId(0);
        if (_init_action->compute(v, f, _visited_variables)) {
          variable_deque.push_back(partition_->variable(other_id));
          _visited_variables.insert(other_id);
        }
      }
      variable_deque.pop_front();
    }
  }
} // namespace srrg2_hipe
