#include "breadth_first_splitter.h"

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  BreadthFirstSplitter::BreadthFirstSplitter() : _init_action(new OdometryPropagation){
  }

  const VariableBase* BreadthFirstSplitter::_computeOuterDegreeMap() {
    const auto& variables           = _graph->variables();
    int max_degree_pose = -1;
    const VariableBase* conditioner = nullptr;
    for (const auto& id_variable : variables) {
      VariableBase::Id id                           = id_variable.first;
      const VariableBase* v                         = id_variable.second;
      auto lower                                    = _graph->lowerFactor(v);
      auto upper                                    = _graph->upperFactor(v);
      int degree = std::distance(lower, upper);
      _variable_outer_degree_map[id]                = degree;
      const VariableSE3Base* anchor_type_3d         = dynamic_cast<const VariableSE3Base*>(v);
      if (!anchor_type_3d) {
        continue;
      }
      if (v->status() == VariableBase::Status::Fixed) {
        conditioner     = v;
        max_degree_pose = std::numeric_limits<int>::max();
        continue;
      }
      if (degree > max_degree_pose) {
        max_degree_pose = degree;
        conditioner     = v;
      }
    }
    return conditioner;
  }
  void BreadthFirstSplitter::_determineAnchorAndCreatePartition(
    const VariablePtrVector& variables_,
    const FactorPtrVector& factors_,
    const VariableIdSet& boundary_variables_) {
    if (variables_.empty() || factors_.empty() || variables_.size() <= boundary_variables_.size()) {
      return;
    }
    int max_degree = -1;
    VariableBase::Id anchor_id = -1;
    for (const VariableBase* v : variables_) {
      if (boundary_variables_.count(v->graphId())) {
        continue;
      }
      auto lower                            = _graph->lowerFactor(v);
      auto upper = _graph->upperFactor(v);
      int degree = std::distance(lower, upper);
      const VariableSE3Base* anchor_type_3d = dynamic_cast<const VariableSE3Base*>(v);
      if (!anchor_type_3d) {
        continue;
      }

      if (v->status() == VariableBase::Status::Fixed) {
        anchor_id = v->graphId();
        break;
      }
      if (degree > max_degree) {
        max_degree  = degree;
        anchor_id   = v->graphId();
      }
    }

    if (anchor_id < 0) {
      return;
    }

    _manager->createPartition(variables_, factors_, anchor_id, boundary_variables_);
  }

  void BreadthFirstSplitter::_constructPartition(const VariableBase* root_, bool root_is_boundary) {
    // tg add root to variable deque
    std::deque<const VariableBase*> variable_deque;
    variable_deque.push_back(root_);
    // tg initialize quantity for partition creation
    VariablePtrVector variables;
    FactorPtrVector factors;
    VariableIdSet boundary_variables;
    if (root_is_boundary) {
      boundary_variables.insert(root_->graphId());
    }
    // tg extract break conditions
    const int& min_variables  = param_min_partition_variables.value();
    const int& min_num_levels = param_min_partition_diameter.value();
    int num_levels            = 0;
    // tg while deque is not empty
    while (!variable_deque.empty()) {
      // tg take front of the deque and extract all the factors
      // where the variable is involved
      const VariableBase* top = variable_deque.front();
      VariableBase* v         = const_cast<VariableBase*>(top);
      VariableBase::Id v_id   = v->graphId();
      auto lower              = _graph->lowerFactor(top);
      auto upper              = _graph->upperFactor(top);
      // add variable to partition and initialize degree
      variables.push_back(v);
      // tg add a new level in the diameter counter
      ++num_levels;
      // tg for all the factors
      for (auto it = lower; it != upper; ++it) {
        auto it_var = _variable_outer_degree_map.find(v_id);
        assert(it_var != _variable_outer_degree_map.end() &&
               "BreadthFirstSplitter::_constructPartition| bookeeping variables is broken");
        it_var->second--;
        // if its already processed skip
        const FactorBase* f = it->second;
        if (_processed_factors.count(f->graphId())) {
          continue;
        }
        // tg add factor to partition and to processed factors
        FactorBase* factor = const_cast<FactorBase*>(f);
        factors.push_back(factor);
        _processed_factors.insert(f->graphId());
        if (f->numVariables() < 2) {
          continue;
        }
        // tg determine id of the other pose variable
        VariableBase::Id other_id =
          f->variableId(0) == v->graphId() ? f->variableId(1) : f->variableId(0);
        auto it_other = _variable_outer_degree_map.find(other_id);
        assert(it_other != _variable_outer_degree_map.end() &&
               "BreadthFirstSplitter::_constructPartition| bookeeping variables is broken");
        it_other->second--;
        // tg check if the variable is an open variable in the graph (variable to be expanded)
        auto is_open = std::find(_open_variables.begin(), _open_variables.end(), other_id);
        if (is_open != _open_variables.end()) {
          boundary_variables.insert(*is_open);
          variable_deque.push_back(_graph->variable(other_id));
          _open_variables.erase(is_open);
          continue;
        }
        // tg initialize other variable from current vertex and factor
        if (_init_action->compute(v, f, _visited_variables)) {
          variable_deque.push_back(_graph->variable(other_id));
          _visited_variables.insert(other_id);
        }
      }
      variable_deque.pop_front();
      // tg check current size and diameter, if sufficient create a partition
      int current_variables_size = variables.size() + variable_deque.size();
      if (current_variables_size > min_variables && num_levels > min_num_levels) {
        break;
      }
    }
    // tg check which non-expanded variables are eligible to be boundary variables
    for (const VariableBase* vs : variable_deque) {
      VariableBase::Id vs_id = vs->graphId();
      VariableBase* v        = const_cast<VariableBase*>(vs);
      variables.push_back(v);
      if (_variable_outer_degree_map.at(vs_id) > 0) {
        _open_variables.push_back(vs_id);
        boundary_variables.insert(vs_id);
      }
    }

    this->_determineAnchorAndCreatePartition(variables, factors, boundary_variables);
  }

  void BreadthFirstSplitter::compute() {
    assert(_graph && "BreadthFirstSplitter::compute| no graph, call setPartitionManager");
    assert(_manager && "BreadthFirstSplitter::compute| no manager, call setPartitionManager");
    // tg connect factors and variables
    _graph->bindFactors();
    // tg determine root of the expansion as the maximum degree variable in the graph
    const VariableBase* root = _computeOuterDegreeMap();
    // tg compute first partition
    _constructPartition(root, false);
    _conditioner = root->graphId();
    _visited_variables.insert(root->graphId());
    // tg while there are variables to be expanded
    while (!_open_variables.empty()) {
      // tg extract front variable
      VariableBase::Id id   = _open_variables.front();
      const VariableBase* v = _graph->variable(id);
      _open_variables.pop_front();
      // tg expand next partition
      _constructPartition(v, true);
    }
  }
} // namespace srrg2_hipe
