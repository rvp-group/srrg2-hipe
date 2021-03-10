#include "isoperimetrical_splitter.h"

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  IsoperimetricSplitter::IsoperimetricSplitter() : _bipartitioner(new IsoperimetricBipartitioner) {
  }

  void IsoperimetricSplitter::compute() {
    // tg partition the initial graph
    _bipartitioner->setGraph(_graph);
    _bipartitioner->computeEdgeWeights();
    auto first_split = _bipartitioner->compute();
    // tg add splitting in the deque
    _partitions_deque.emplace_back(first_split.first);
    _partitions_deque.emplace_back(first_split.second);
    // tg while the partition deque is not empty
    while (!_partitions_deque.empty()) {
      // tg split
      PartitionPtr current_partition = _partitions_deque.front();
      _split(current_partition);
      // tg remove the front partition
      _partitions_deque.pop_front();
    }
  }

  void IsoperimetricSplitter::_split(PartitionPtr& partition_) {
    // tg check if required max size is reached
    int max_size       = param_max_partition_size.value();
    int partition_size = partition_->variables().size();
    std::cerr << "Partition " << partition_->anchorId() << " " << partition_size << " Max size "
              << max_size << std::endl;
    if (partition_size <= max_size) {
      _manager->addBoundaryVariablesToPartition(partition_, partition_->boundaryVariableIds());
      _manager->addPartition(partition_);
      return;
    }
    partition_->bindFactors();
    // tg split
    _bipartitioner->setGraph(partition_);
    auto split = _bipartitioner->compute();
    // tg divide the boundary variable
    auto& original_boundary_variables = partition_->boundaryVariableIds();
    for (const VariableBase::Id& id : original_boundary_variables) {
      if (_bipartitioner->isVariableOnLeftSideOfThePartition(id)) {
        split.first->addBoundaryVariableIds(VariableIdSet{id});
      } else {
        split.second->addBoundaryVariableIds(VariableIdSet{id});
      }
    }
    std::cerr << "Total graph variables : " << partition_->variables().size()
              << " Boundary variables left : " << split.first->boundaryVariableIds().size()
              << " Boundary variables right : " << split.second->boundaryVariableIds().size()
              << " Left partition size : " << split.first->variables().size()
              << " Right partition size : " << split.second->variables().size() << std::endl;
    // tg add the partitions to deque
    _partitions_deque.emplace_back(split.first);
    _partitions_deque.emplace_back(split.second);
  }
} // namespace srrg2_hipe
