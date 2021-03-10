#include "isoperimetrical_bipartitioner.h"
#include <srrg_solver/variables_and_factors/types_2d/instances.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  IsoperimetricBipartitioner::IsoperimetricBipartitioner() {
    cholmod_start(&_cholmodCommon);
    // tg set up cholmod parameters
    _cholmodCommon.nmethods           = 0;
    _cholmodCommon.method[0].ordering = CHOLMOD_AMD;
    _cholmodCommon.supernodal         = CHOLMOD_AUTO;
  }

  IsoperimetricBipartitioner::~IsoperimetricBipartitioner() {
    // tg free workspace and terminate cholmod
    cholmod_free_sparse(&_A, &_cholmodCommon);
    cholmod_free_dense(&_b, &_cholmodCommon);
    cholmod_free_dense(&_x, &_cholmodCommon);
    cholmod_free_factor(&_L, &_cholmodCommon);
    cholmod_finish(&_cholmodCommon);
  }

  void IsoperimetricBipartitioner::computeEdgeWeights() {
    _factor_to_weight.clear();
    for (const auto& id_factor : _graph->factors()) {
      FactorBase* f = const_cast<FactorBase*>(id_factor.second);
      f->compute(true, false);
      auto stats                         = f->stats();
      double negative_log_likelihood     = static_cast<double>(stats.chi);
      double w_ij = 1. / (std::max(0.01, std::exp(-negative_log_likelihood)));
      // double w_ij                        = std::exp(-negative_log_likelihood);
      int num_variables                  = f->numVariables();
      std::vector<VariableBase::Id> variables;
      variables.reserve(num_variables);
      for (int i = 0; i < num_variables; ++i) {
        variables.emplace_back(f->variableId(i));
      }
      _factor_to_weight[f->graphId()] = std::make_pair(variables, w_ij);
    }
  }

  void IsoperimetricBipartitioner::_constructLaplacian() {
    _variable_id_to_index.clear();
    _variable_id_degree.clear();
    int hessian_index = -1;
    // select a anchor for the partition
    VariableBase::Id anchor_variable_id = _graph->variables().begin().key();
    // tg for each variable
    for (const auto& id_variable : _graph->variables()) {
      const VariableBase* v       = id_variable.second;
      const VariableBase::Id v_id = v->graphId();
      // tg compute the degree of the vertex
      const auto lower = _graph->lowerFactor(v);
      const auto upper = _graph->upperFactor(v);
      double degree    = 0;
      for (auto it = lower; it != upper; ++it) {
        double weight = _factor_to_weight[it->second->graphId()].second;
        degree += weight;
      }
      _variable_id_degree[v_id] = degree;
      if (v_id != anchor_variable_id) {
        _variable_id_to_index[v_id] = ++hessian_index;
      }
    }
    // tg initialize laplacian triplets
    _laplacian.clear();
    size_t num_variables = _variable_id_degree.size();
    _laplacian.reserve(num_variables * (num_variables + 1) / 2);
    // tg for each factor in the graph
    for (const auto& id_factor : _graph->factors()) {
      const FactorBase* f           = id_factor.second;
      VariableVectorWeightPair edge = _factor_to_weight[f->graphId()];
      const std::vector<VariableBase::Id>& variables = edge.first;
      const double& weight                           = edge.second;
      // tg inizialize diagonal and off diagonal terms
      LaplacianEntry l_ij;
      int num_variables = variables.size();
      for (int r = 0; r < num_variables; ++r) {
        VariableBase::Id id_row = variables[r];
        if (id_row == anchor_variable_id) {
          continue;
        }
        for (int c = r; c < num_variables; ++c) {
          VariableBase::Id id_col = variables[c];
          if (id_col == anchor_variable_id) {
            continue;
          }
          // tg add just the upper triangular part of the laplacian
          int row_index = _variable_id_to_index[id_row];
          int col_index = _variable_id_to_index[id_col];
          if (row_index > col_index) {
            std::swap(row_index, col_index);
          }
          std::get<0>(l_ij) = row_index;
          std::get<1>(l_ij) = col_index;
          if (id_row == id_col) {
            std::get<2>(l_ij) = _variable_id_degree.at(id_row);
          } else {
            std::get<2>(l_ij) = -weight;
          }
          _laplacian.emplace_back(l_ij);
        }
      }
    }
    // tg sort element column major
    std::sort(_laplacian.begin(),
              _laplacian.end(),
              [](const LaplacianEntry& a, const LaplacianEntry& b) -> bool {
                return std::get<1>(a) < std::get<1>(b) ||
                       (std::get<1>(a) == std::get<1>(b) && std::get<0>(a) < std::get<0>(b));
              });
    // tg remove duplicates
    auto last =
      std::unique(_laplacian.begin(),
                  _laplacian.end(),
                  [](const LaplacianEntry& a, const LaplacianEntry& b) -> bool {
                    return std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b);
                  });
    _laplacian.erase(last, _laplacian.end());
  }

  void IsoperimetricBipartitioner::_buildLinearSystem() {
    // clean internal structure
    cholmod_free_sparse(&_A, &_cholmodCommon);
    cholmod_free_dense(&_b, &_cholmodCommon);
    cholmod_free_factor(&_L, &_cholmodCommon);

    const size_t num_variables = _variable_id_to_index.size();
    const size_t nnz           = _laplacian.size();

    cholmod_triplet laplacian_triplet;

    laplacian_triplet.itype = CHOLMOD_INT;
    laplacian_triplet.xtype = CHOLMOD_REAL;
    laplacian_triplet.dtype = CHOLMOD_DOUBLE;

    laplacian_triplet.nrow  = num_variables;
    laplacian_triplet.ncol  = num_variables;
    laplacian_triplet.nzmax = nnz;
    laplacian_triplet.nnz   = nnz;
    laplacian_triplet.stype = 1;

    std::unique_ptr<int[]> row_indices(new int[nnz]);
    std::unique_ptr<int[]> col_indices(new int[nnz]);
    std::unique_ptr<double[]> values(new double[nnz]);

    laplacian_triplet.i = row_indices.get();
    laplacian_triplet.j = col_indices.get();
    laplacian_triplet.x = values.get();
    for (size_t idx = 0; idx < nnz; ++idx) {
      LaplacianEntry l = _laplacian.at(idx);
      row_indices[idx] = std::get<0>(l);
      col_indices[idx] = std::get<1>(l);
      values[idx]      = std::get<2>(l);
    }
    _A = cholmod_triplet_to_sparse(&laplacian_triplet, laplacian_triplet.nnz, &_cholmodCommon);
    _b = cholmod_ones(_A->nrow, 1, _A->xtype, &_cholmodCommon);
    _L = cholmod_analyze(_A, &_cholmodCommon);
    cholmod_factorize(_A, _L, &_cholmodCommon);
  }

  void IsoperimetricBipartitioner::_solveLinearSystem() {
    cholmod_free_dense(&_x, &_cholmodCommon);
    _x                         = cholmod_solve(CHOLMOD_A, _L, _b, &_cholmodCommon);
    const double* values       = (const double*) _x->x;
    const size_t num_variables = _A->nrow;
    _partition_vector.clear();
    _partition_vector.reserve(num_variables);
    for (size_t i = 0; i < num_variables; ++i) {
      assert(values[i] == values[i] &&
             "IsoperimetricBipartitioner::_solveLinearSystem| Nan detected");
      // std::cerr << values[i] << std::endl;
      _partition_vector.emplace_back(values[i]);
    }
  }

  void IsoperimetricBipartitioner::_computeCutThreshold() {
    assert(!_partition_vector.empty() &&
           "IsoperimetricBipartitioner::_computeCutThreshold())| Partion vector is void");

    // tg avoid copy before sorting (MOST PROBABLY THIS IS MORE EXPENSIVE)
    size_t num_variables = _partition_vector.size();
    std::vector<size_t> partition_indices(_partition_vector.size(), 0);
    std::iota(partition_indices.begin(), partition_indices.end(), 0);
    std::sort(partition_indices.begin(),
              partition_indices.end(),
              [this](const size_t& a, const size_t& b) -> bool {
                return this->_partition_vector.at(a) < this->_partition_vector.at(b);
              });
    // tg compute median
    size_t half_size        = num_variables / 2;
    double median_threshold = 0;
    if (half_size % 2 == 0) {
      median_threshold += _partition_vector[partition_indices.at(half_size - 1)];
      median_threshold += _partition_vector[partition_indices.at(half_size)];
      median_threshold /= 2;
    } else {
      median_threshold = _partition_vector[partition_indices.at(half_size)];
    }
    _threshold = median_threshold;
  }

  void
  IsoperimetricBipartitioner::_addVariables(PartitionPtr& star,
                                            std::set<VariableBase::Id>& variable_set_,
                                            const std::set<VariableBase::Id>& boundary_variables_) {
    int max_degree            = -1;
    VariableBase::Id anchor_id = -1;
    variable_set_.insert(boundary_variables_.begin(), boundary_variables_.end());
    for (const VariableBase::Id& v_id : variable_set_) {
      VariableBase* v = _graph->variable(v_id);
      star->addVariable(v);
      if (boundary_variables_.count(v_id)) {
        continue;
      }
      // double cast to avoid template
      const VariableSE3Base* anchor_type_3d = dynamic_cast<const VariableSE3Base*>(v);
      const VariableSE2Base* anchor_type_2d = dynamic_cast<const VariableSE2Base*>(v);
      if (!anchor_type_2d && !anchor_type_3d) {
        continue;
      }
      auto it_degree = _variable_id_degree.find(v_id);
      if (it_degree != _variable_id_degree.end() && it_degree->second > max_degree) {
        anchor_id  = v_id;
        max_degree = it_degree->second;
      }
    }
    VariableBase* v_anchor = _graph->variable(anchor_id);
    if (!v_anchor) {
      throw std::runtime_error(
        "IsoperimetricBipartitioner::_addVariables| Cannot const_cast variable");
    }
    star->setAnchor(v_anchor);
    star->addBoundaryVariableIds(boundary_variables_);
  }

  IsoperimetricBipartitioner::GraphBipartition IsoperimetricBipartitioner::compute() {
    _constructLaplacian();
    _buildLinearSystem();
    _solveLinearSystem();
    _computeCutThreshold();
    // tg temporary structures for variables and factors
    std::set<VariableBase::Id> left_variables, right_variables, boundary_variables;
    GraphBipartition returned(std::make_shared<Partition>(), std::make_shared<Partition>());
    // tg compute left and right variables for the partition
    for (const auto& id_variable : _graph->variables()) {
      const VariableBase::Id& v_id = id_variable.first;
      auto it_hessian_index        = _variable_id_to_index.find(v_id);
      if (it_hessian_index == _variable_id_to_index.end()) {
        left_variables.insert(v_id);
        continue;
      }
      const int& hessian_index = it_hessian_index->second;
      if (_partition_vector[hessian_index] < _threshold) {
        left_variables.insert(v_id);
      } else {
        right_variables.insert(v_id);
      }
    }

    // tg move the factor from the border to the inside of a partition
    for (const auto& id_factor : _graph->factors()) {
      FactorBase* f               = id_factor.second;
      int num_variables_in_factor = f->numVariables();
      // tg check where the first variable is and store the degree
      VariableBase::Id first_variable = f->variableId(0);
      bool is_left                    = left_variables.count(first_variable);
      const double& degree_first      = _variable_id_degree[first_variable];
      double factor_weight            = _factor_to_weight[f->graphId()].second;
      double gain_first               = degree_first - factor_weight;
      // tg for all other variables
      for (int i = 1; i < num_variables_in_factor; ++i) {
        // check where the variable is and store the degree
        VariableBase::Id other_variable = f->variableId(i);
        bool other_is_left              = left_variables.count(other_variable);
        const double& degree_other      = _variable_id_degree[other_variable];
        double gain_other               = degree_other - factor_weight;
        // tg if both are left put the factor in the left partition
        if (is_left && other_is_left) {
          returned.first->addGraphFactor(f);
          // tg same for right
        } else if (!is_left && !other_is_left) {
          returned.second->addGraphFactor(f);
          // tg if one is left and the other is right
        } else if (is_left && !other_is_left) {
          // tg move the variable with smallest degree to the other partition
          // and add the factor
          if (gain_first > gain_other) {
            boundary_variables.insert(first_variable);
            returned.second->addGraphFactor(f);
          } else {
            boundary_variables.insert(other_variable);
            returned.first->addGraphFactor(f);
          }
          // tg same with opposite position
        } else if (gain_first > gain_other) {
          boundary_variables.insert(first_variable);
          returned.first->addGraphFactor(f);
        } else {
          boundary_variables.insert(other_variable);
          returned.second->addGraphFactor(f);
        }
      }
    }
    _addVariables(returned.first, left_variables, boundary_variables);
    _addVariables(returned.second, right_variables, boundary_variables);
    return returned;
  }
} // namespace srrg2_hipe
