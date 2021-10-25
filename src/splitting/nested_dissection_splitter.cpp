#include "nested_dissection_splitter.h"
#include "partitions_optimization_tree.h"
#include <srrg_solver/variables_and_factors/types_3d/instances.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;
  using BlockIndices = std::pair<int, int>;
  struct CompareBlockIndices {
    bool operator()(const BlockIndices& a, const BlockIndices& b) const {
      return a.second < b.second || (a.second == b.second && a.first < b.first);
    }
  };
  using BlockIndicesSet = std::set<BlockIndices, CompareBlockIndices>;

  NestedDissectionSplitter::NestedDissectionSplitter() {
    cholmod_start(&_cholmodCommon);
    // tg set up cholmod parameters
    _cholmodCommon.nmethods           = 0;
    _cholmodCommon.method[0].ordering = CHOLMOD_AMD;
    _cholmodCommon.supernodal         = CHOLMOD_AUTO;
  }

  NestedDissectionSplitter::~NestedDissectionSplitter() {
    cholmod_free_sparse(&_A, &_cholmodCommon);
    cholmod_finish(&_cholmodCommon);
  }

  void NestedDissectionSplitter::_fillVariableIndeces() {
    assert(_graph && "NestedDissectionSplitter::_fillVariableIndeces| graph not set");
    int hessian_index = -1;
    for (const auto& id_variable : _graph->variables()) {
      const VariableBase* v       = id_variable.second;
      const VariableBase::Id v_id = v->graphId();
      // tg map from id to index
      _variable_to_index[v_id] = ++hessian_index;
      _index_to_variable[hessian_index] = v_id;
      if (v->status() == VariableBase::Status::Fixed) {
        _conditioner = v_id;
      }
    }
  }

  void NestedDissectionSplitter::_buildApproximateHessianPattern() {
    assert(!_variable_to_index.empty() &&
           "NestedDissectionSplitter::_buildApproximateHessianPattern| miss "
           "variable indices...something is wrong in the graph");
    cholmod_free_sparse(&_A, &_cholmodCommon);
    // tg fill block structure
    int num_variables = _variable_to_index.size();
    BlockIndicesSet block_layout;
    // tg for each factor extract indices for the variables involved
    for (const auto& id_factor : _graph->factors()) {
      const FactorBase* f = id_factor.second;
      int num_variables_in_factor = f->numVariables();
      for (int r = 0; r < num_variables_in_factor; ++r) {
        VariableBase::Id id_row = f->variableId(r);
        for (int c = r; c < num_variables_in_factor; ++c) {
          VariableBase::Id id_col = f->variableId(c);
          // upper triangular part of H
          int row_index           = _variable_to_index[id_row];
          int col_index           = _variable_to_index[id_col];
          if (row_index > col_index) {
            std::swap(row_index, col_index);
          }
          block_layout.insert(std::make_pair(row_index, col_index));
        }
      }
    }

    // tg construct the pattern of H in cholmod sparse format partitionting from triplets
    const size_t nnz = block_layout.size();
    cholmod_triplet triplet;

    triplet.itype = CHOLMOD_INT;
    triplet.xtype = CHOLMOD_REAL;
    triplet.dtype = CHOLMOD_DOUBLE;

    triplet.nrow  = num_variables;
    triplet.ncol  = num_variables;
    triplet.nzmax = nnz;
    triplet.nnz   = nnz;
    triplet.stype = 1;

    std::unique_ptr<int[]> row_indices(new int[nnz]);
    std::unique_ptr<int[]> col_indices(new int[nnz]);
    std::unique_ptr<double[]> values(new double[nnz]);

    triplet.i = row_indices.get();
    triplet.j = col_indices.get();
    triplet.x = values.get();
    int idx   = 0;
    for (auto it = block_layout.begin(); it != block_layout.end(); ++it, ++idx) {
      const BlockIndices& element = *it;
      row_indices[idx] = element.first;
      col_indices[idx] = element.second;
      values[idx]      = 1.;
    }
    _A = cholmod_triplet_to_sparse(&triplet, triplet.nnz, &_cholmodCommon);
  }

  void NestedDissectionSplitter::_collapseSeparatorTreeNodes(ComponentIdToPartition& partitions,
                                                             int* parents) {
    int min_variables  = param_min_partition_variables.value();
    int num_components = partitions.size();
    for (int c = 0; c < num_components; ++c) {
      // tg for each partition check all the children partitions in the separator tree
      if (!partitions.count(c)) {
        continue;
      }
      PartitionPtr& parent = partitions[c];
      for (int child_idx = 0; child_idx < num_components; ++child_idx) {
        // tg if a child is found check if it has enough variables
        if (parents[child_idx] == c && partitions.count(child_idx)) {
          PartitionPtr& child    = partitions[child_idx];
          int variables_in_child = child->variables().size();
          if (variables_in_child < min_variables) {
            // tg if not merge child into parent
            for (const auto& e : child->variables()) {
              parent->addVariable(e.second);
            };
            partitions.erase(child_idx);
            std::replace_if(parents,
                            parents + num_components,
                            [child_idx](const int& a) -> bool { return a == child_idx; },
                            c);
          }
        }
      }
    }
  }

  void NestedDissectionSplitter::compute() {
    assert(_manager && "NestedDissectionSplitter::compute| Partition manager not found");
    _fillVariableIndeces();
    _buildApproximateHessianPattern();
    // tg extract min variables and initialize pointers for cholmod
    int num_variables = _variable_to_index.size();
    std::unique_ptr<int[]> fset(new int[num_variables]), perm(new int[num_variables]);
    std::iota(fset.get(), fset.get() + num_variables, 0);
    std::unique_ptr<int[]> parents(new int[num_variables]), membership(new int[num_variables]);
    // tg nested dissection by cholmod
    int num_components = cholmod_nested_dissection(
      _A, fset.get(), _A->nrow, perm.get(), parents.get(), membership.get(), &_cholmodCommon);
    // tg fill component to partition map with corresponding variables
    ComponentIdToPartition partitions;
    for (int v = 0; v < num_variables; ++v) {
      const int& component_index = membership[v];
      bool found                 = partitions.count(component_index);
      if (!found) {
        partitions[component_index] = PartitionPtr(new Partition(1));
      }
      PartitionPtr& s      = partitions[component_index];
      VariableBase* vertex = _graph->variable(_index_to_variable[v]);
      s->addVariable(vertex);
    }
    // prune separator tree (partitions)
    // _collapseSeparatorTreeNodes(partitions, parents.get());
    // tg bookeeping variable partition
    std::map<VariableBase::Id, int> variable_to_partition;
    for (const auto& comp_index_partition : partitions) {
      const int& c                  = comp_index_partition.first;
      const PartitionPtr& partition = comp_index_partition.second;
      for (const auto& id_variable : partition->variables()) {
        variable_to_partition.insert(std::make_pair(id_variable.first, c));
        // assert(result.second && "NestedDissectionSplitter::compute| variable appear in multiple "
        //                         "partitions after nested dissection algorithm");
      }
    }
    num_components = partitions.size();
    std::cerr << "NestedDissectionSplitter::compute| Graph partitioned in " << num_components
              << " partitions" << std::endl;
    std::cerr << "NestedDissectionSplitter::compute| Num parents " << num_components
              << " partitions" << std::endl;
    std::set<FactorBase::Id> processed_factors;
    // tg add factors and boundary variables to each partition
    for (auto& id_partition : partitions) {
      PartitionPtr& partition = id_partition.second;
      // tg nested loop to check if any of the factor in the partition
      // has an end-point in the parent
      for (const auto& id_variable : partition->variables()) {
        VariableBase* v = id_variable.second;
        auto lower      = _graph->lowerFactor(v);
        auto upper      = _graph->upperFactor(v);
        for (auto it = lower; it != upper; ++it) {
          FactorBase* f = it->second;
          if (processed_factors.count(f->graphId())) {
            continue;
          }
          int num_variables   = f->numVariables();
          switch (num_variables) {
            case 1:
              partition->addGraphFactor(f);
              processed_factors.insert(f->graphId());
              break;
            case 2:
              VariableBase::Id v_id = v->graphId();
              // tg get the id of the other variable
              VariableBase::Id other_id =
                v_id == f->variableId(0) ? f->variableId(1) : f->variableId(0);
              // tg if both variable are in the same partition add the factor to it
              if (partition->variable(other_id)) {
                partition->addGraphFactor(f);
              } else {
                PartitionPtr other_partition = partitions.at(variable_to_partition.at(other_id));
                other_partition->addVariable(v);
                other_partition->addGraphFactor(f);
                partition->addBoundaryVariableIds(VariableIdSet{v_id});
                other_partition->addBoundaryVariableIds(VariableIdSet{v_id});
              }
              processed_factors.insert(f->graphId());
              break;
          }
        }
      }
    }
    // tg determine anchor as the pose variable with the highest degree,
    // add partition to manager and compute condensed factors
    for (auto& partition : partitions) {
      _determineAnchorAndAddToManager(partition.second);
    }
    PartitionsOptimizationTree opt;
    opt.setPartitionManager(_manager);
    opt.compute();
  }

  void NestedDissectionSplitter::_determineAnchorAndAddToManager(PartitionPtr& partition_) {
    int max_degree = -1, degree = -1;
    VariableBase::Id anchor_id              = -1;
    const VariableIdSet& boundary_variables = partition_->boundaryVariableIds();
    for (const auto& v_id : partition_->variables()) {
      if (boundary_variables.count(v_id.first)) {
        continue;
      }

      VariableBase* v = v_id.second;
      // double cast to avoid template
      VariableSE3Base* anchor_type_3d = dynamic_cast<VariableSE3Base*>(v);
      if (!anchor_type_3d) {
        continue;
      }
      auto lower       = _graph->lowerFactor(v);
      auto upper       = _graph->upperFactor(v);
      degree           = std::distance(lower, upper);
      if (degree > max_degree) {
        anchor_id  = v_id.first;
        max_degree = degree;
      }
    }
    VariableBase* v_anchor = _graph->variable(anchor_id);
    assert(v_anchor &&
           "NestedDissectionSplitter::_determineAnchorAndAddToManager| anchor not found");
    partition_->setAnchor(v_anchor);
    partition_->bindFactors();
    _manager->addPartition(partition_);
  }
} // namespace srrg2_hipe
