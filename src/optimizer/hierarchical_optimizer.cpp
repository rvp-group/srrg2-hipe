#include "hierarchical_optimizer.h"
#include <chordal_types/chordal_initializer.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  void HierarchicalOptimizer::computePartitions() {
    assert(_graph && "HierarchicalOptimizer::computePartitions| no graph found, call setGraph");
    PartitionManagerPtr manager = param_partition_manager.value();
    manager->setGraph(_graph);

    SplitterBasePtr splitter = param_splitter.value();
    splitter->setPartitionManager(manager.get());

    splitter->compute();
    VariableBase::Id conditioner = param_splitter.value()->conditioner();
    _anchor                      = _graph->variable(conditioner);
    // tg store the skeleton of the graph
    _coarse_grained_variables.clear();
  }

  void HierarchicalOptimizer::_initializeSkeleton() {
    assert(_graph && "HierarchicalOptimizer::_initializeSkeleton| no graph found, call setGraph");
    // tg chordal initialization for high layer
    _anchor->setStatus(VariableBase::Status::Fixed);
    initializeChordal(_graph, &_coarse_grained_variables, 1);
    _anchor->setStatus(VariableBase::Status::Active);
  }

  void HierarchicalOptimizer::optimizeSkeleton() {
    assert(_graph && "HierarchicalOptimizer::optimizeSkeleton| no graph found, call setGraph");
    _initializeSkeleton();
    _anchor->setStatus(VariableBase::Status::Fixed);
    // tg optimize structure using condensed factors
    SolverPtr solver         = param_solver.value();
    int low_level_iterations = solver->param_max_iterations.value(0);
    solver->param_max_iterations.setValue(0, 0);
    solver->setGraph(_graph);
    solver->compute();
    _anchor->setStatus(VariableBase::Status::Active);
    solver->param_max_iterations.setValue(0, low_level_iterations);
    if (param_verbose.value()) {
      std::cerr << "HierarchicalOptimizer::optimizeSkeleton|" << param_solver->iterationStats()
                << std::endl;
    }
  }

  void HierarchicalOptimizer::optimizeLowLevel() {
    assert(_graph && "HierarchicalOptimizer::optimizeLowLevel| no graph found, call setGraph");
    _anchor->setStatus(VariableBase::Status::Fixed);
    SolverPtr solver              = param_solver.value();
    int coarse_grained_iterations = solver->param_max_iterations.value(1);
    auto alg                      = solver->param_algorithm.value();
    solver->param_algorithm.setValue(param_low_level_alg.value());
    solver->param_max_iterations.setValue(1, 0);
    solver->setGraph(_graph);
    solver->compute();
    _anchor->setStatus(VariableBase::Status::Active);
    solver->param_max_iterations.setValue(1, coarse_grained_iterations);
    solver->param_algorithm.setValue(alg);
  }

  void HierarchicalOptimizer::optimizeRemainingVariables() {
    // tg fix hierarchical structure
    for (const VariableBase::Id& v_id : _coarse_grained_variables) {
      _graph->variable(v_id)->setStatus(VariableBase::Status::Fixed);
    }
    initializeChordal(_graph);
    // tg optimize free variables
    SolverPtr solver              = param_solver.value();
    int coarse_grained_iterations = solver->param_max_iterations.value(1);
    solver->param_max_iterations.setValue(1, 0);
    solver->setGraph(_graph);
    solver->compute();
    solver->param_max_iterations.setValue(1, coarse_grained_iterations);
    // tg unlock skeleton variables
    for (const VariableBase::Id& v_id : _coarse_grained_variables) {
      _graph->variable(v_id)->setStatus(VariableBase::Status::Active);
    }

    if (param_verbose.value()) {
      std::cerr << "HierarchicalOptimizer::optimizeRemainingVariables|"
                << param_solver->iterationStats() << std::endl;
    }
  }

  void HierarchicalOptimizer::initializeGraph() {
    // initialize graph using the HiPe strategy
    computePartitions();
    optimizeSkeleton();
    optimizeRemainingVariables();
  }

  void HierarchicalOptimizer::compute() {
    initializeGraph();
    // optimize fine-grained graph
    optimizeLowLevel();
    if (param_verbose.value()) {
      std::cerr << "HierarchicalOptimizer::compute|" << param_solver->iterationStats() << std::endl;
    }
  }

  void HierarchicalOptimizer::writeLevel(const std::string& filename, const int& level_) {
    _graph->setSerializationLevel(level_);
    _graph->write(filename);
    _graph->setSerializationLevel(0);
  }

} // namespace srrg2_hipe
