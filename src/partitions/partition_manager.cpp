#include "partition_manager.h"
#include <chordal_types/chordal_initializer.h>

namespace srrg2_hipe {

  PartitionPtr PartitionManager::partition(const VariableBase::Id& partition_anchor_id_) {
    auto element_iterator = _partitions.find(partition_anchor_id_);
    if (element_iterator != _partitions.end()) {
      return element_iterator->second;
    }
    return nullptr;
  }

  size_t PartitionManager::variableConnectivity(const VariableBase::Id& variable_id_) {
    return _variables_to_partitions.count(variable_id_);
  }

  void PartitionManager::addPartition(PartitionPtr partition_) {
    _partitions.insert(std::make_pair(partition_->anchorId(), partition_));
    // add entities in the internal maps
    std::for_each(partition_->variables().begin(),
                  partition_->variables().end(),
                  [this, partition_](const auto& id_v) {
                    this->_variables_to_partitions.insert(std::make_pair(id_v.first, partition_));
                  });
  }

  bool PartitionManager::getConnectedPartitions(PartitionPtrVector& partitions_,
                                                const VariableBase::Id& variable_id_) {
    partitions_.clear();
    auto start_and_end = _variables_to_partitions.equal_range(variable_id_);
    if (start_and_end.first == start_and_end.second) {
      return false;
    }
    int num_elements = _variables_to_partitions.count(variable_id_);
    partitions_.reserve(num_elements);
    for (auto it = start_and_end.first; it != start_and_end.second; ++it) {
      partitions_.emplace_back(it->second);
    }
    return true;
  }

  void PartitionManager::optimizePartition(PartitionPtr partition_) {
    SolverPtr solver = this->param_solver.value();
    assert(solver && "PartitionManager::optimizePartition| solver not present");
    assert(partition_ && "PartitionManager::optimizePartition| partition ptr is void");
    // tg set graph and fix anchor
    solver->setGraph(partition_);
    VariableBase* anchor = partition_->anchor();
    anchor->setStatus(VariableBase::Status::Fixed);
    initializeChordal(partition_);
    // tg optimize partition
    solver->compute();
    // tg unlock hierarchical structure
    anchor->setStatus(VariableBase::Status::Active);
    if (param_verbose.value()) {
      std::cerr << "PartitionManager | " << solver->iterationStats() << std::endl;
    }
  }

  void PartitionManager::optimizePartition(const VariableBase::Id& anchor_id_) {
    auto it = _partitions.find(anchor_id_);
    if (it == _partitions.end()) {
      std::cerr << "PartitionManager::optimizePartition| Partition is not present in the manager"
                << std::endl;
      return;
    }
    optimizePartition(it->second);
  }

  void
  PartitionManager::addBoundaryVariablesToPartition(PartitionPtr partition_,
                                                    const VariableIdSet& boundary_variables_ids_) {
    assert(partition_ &&
           "PartitionManager::addBoundaryVariablesToPartition| partition ptr is null");
    if (boundary_variables_ids_.empty()) {
      return;
    }
    optimizePartition(partition_);
    // variables pairs pointes used to compute partition edges after the optimization of the
    // partition
    int num_boundary_variables = boundary_variables_ids_.size();
    VariablePairVector variable_pairs;
    variable_pairs.reserve(num_boundary_variables);
    SolverPtr solver = this->param_solver.value();
    // need to access the variables in solver otherwise i will miss the hessian index
    // required for marginal covariance computation
    auto& active_variables_in_partition = solver->activeVariables();
    for (VariableBase* v : active_variables_in_partition) {
      if (boundary_variables_ids_.count(v->graphId())) {
        variable_pairs.emplace_back(std::make_pair(v, v));
      }
    }

    if (variable_pairs.empty()) {
      throw std::runtime_error(
        "PartitionManager::addBoundaryVariablesToPartition| none of the requested "
        "boundary variables is contained in the partition");
    }
    // compute marginal covariance
    MatrixBlockVector marginal_covariances;
    solver->computeMarginalCovariance(marginal_covariances, variable_pairs);
    // compute virtual factors and free the anchor
    _computeVirtualFactors(partition_, variable_pairs, marginal_covariances);
    partition_->addBoundaryVariableIds(boundary_variables_ids_);
  }

  bool PartitionManager::createPartition(const VariablePtrVector& variables,
                                         const FactorPtrVector& factors,
                                         const VariableBase::Id& anchor_id,
                                         const VariableIdSet& boundary_variables_ids) {
    auto found = _partitions.find(anchor_id);
    if (found != _partitions.end()) {
      std::cerr << "PartitionManage::createPartition| Partition with the same anchor already exist"
                << std::endl;
      return false;
    }
    size_t level = factors.front()->level() + 1;
    PartitionPtr partition(new Partition(variables, factors, level));
    VariableBase* anchor = _graph->variable(anchor_id);
    assert(anchor && "PartitionManager::createPartition| anchor id is not present in the graph");
    // set anchor and add boundary variables with corresponding virtual
    // factors
    partition->setAnchor(anchor);
    addBoundaryVariablesToPartition(partition, boundary_variables_ids);
    addPartition(partition);
    return true;
  }

  bool PartitionManager::destroyPartition(const VariableBase::Id& anchor_id_) {
    IdAnchorPartitionMap::iterator partition_it = _partitions.find(anchor_id_);
    if (partition_it == _partitions.end()) {
      return false;
    }
    PartitionPtr& partition = partition_it->second;
    for (const auto& id_variable : partition->variables()) {
      const VariableBase::Id& id = id_variable.first;
      assert(_variables_to_partitions.count(id) &&
             "PartitionManager::destroyPartition| variable is not connected to any partition");
      auto begin_end = _variables_to_partitions.equal_range(id);
      auto element   = std::find_if(begin_end.first,
                                  begin_end.second,
                                  [partition](VariablePartitionMultiMap::value_type& v) -> bool {
                                    return v.second->anchorId() == partition->anchorId();
                                  });
      _variables_to_partitions.erase(element);
    }
    for (FactorBase* f : partition->virtualFactors()) {
      _graph->removeFactor(f);
    }
    _partitions.erase(partition_it);
    return true;
  }

  void PartitionManager::_computeVirtualFactors(PartitionPtr& partition_,
                                                VariablePairVector& variable_pairs_,
                                                MatrixBlockVector& covariances_) {
    size_t num_virtual_factors = variable_pairs_.size();
    // project covariance of the boundary variables in the error space of the virtual factors
    // anchor is assumed to be a pose
    const int num_projectors       = this->param_covariance_projectors.size();
    FactorBasePtr virtual_factor   = nullptr;
    VariableBase* anchor           = partition_->anchor();
    // for each virtual factors
    for (size_t v_idx = 0; v_idx < num_virtual_factors; ++v_idx) {
      // extract variable and covariance block
      VariableBase* variable      = variable_pairs_[v_idx].first;
      MatrixBlockBase* covariance = covariances_[v_idx];
      // try all the projectors until one succed
      for (int proj_idx = 0; proj_idx < num_projectors; ++proj_idx) {
        VirtualFactorCreatorBasePtr proj   = this->param_covariance_projectors.value(proj_idx);
        virtual_factor                     = proj->compute(variable, anchor, covariance);
        if (virtual_factor) {
          break;
        }
      }
      assert(virtual_factor &&
             "PartitionManager::_computeVirtualFactors| none of the projector was "
             "compatible with the type of variables");
      // add virtual factor to graph and the partition
      virtual_factor->setLevel(partition_->level());
      _graph->addFactor(virtual_factor);
      partition_->addVirtualFactor(virtual_factor.get());
    }
  }
} // namespace srrg2_hipe
