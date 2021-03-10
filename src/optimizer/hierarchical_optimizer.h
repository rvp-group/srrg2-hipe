#pragma once
#include <partitions/partition_manager.h>
#include <splitting/breadth_first_splitter.h>
#include <splitting/nested_dissection_splitter.h>
#include <srrg_solver/solver_core/iteration_algorithm_lm.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  class HierarchicalOptimizer : public Configurable {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PARAM(PropertyConfigurable_<PartitionManager>,
          partition_manager,
          "Create and store the partitions. Perform the optimization and add the condensed factors",
          PartitionManagerPtr(new PartitionManager),
          nullptr);

    PARAM(PropertyConfigurable_<Solver>,
          solver,
          "Optimize the graph at both levels",
          SolverPtr(new Solver),
          nullptr);

    PARAM(PropertyConfigurable_<IterationAlgorithmBase>,
          low_level_alg,
          "Iteration algorithm for the low level (original) graph",
          IterationAlgorithmBasePtr(new IterationAlgorithmLM),
          nullptr);

    PARAM(PropertyConfigurable_<SplitterBase>,
          splitter,
          "Split the original graph in variable sets from which the partitions are created",
          SplitterBasePtr(new BreadthFirstSplitter),
          nullptr);

    PARAM(PropertyBool, verbose, "Show statistics of the optimization", true, nullptr);

    inline void setGraph(FactorGraphPtr& graph_) {
      _graph = graph_;
    }

    void computePartitions();
    void optimizeSkeleton();
    void optimizeRemainingVariables();
    void optimizeLowLevel();

    void initializeGraph();
    void compute();

    void writeLevel(const std::string& filename, const int& level_);

  protected:
    void _initializeSkeleton();
    FactorGraphPtr _graph = nullptr;
    VariableIdSet _coarse_grained_variables;
    // tg store the anchor of the optimization
    VariableBase* _anchor = nullptr;
  };

  using HierarchicalOptimizerPtr = std::shared_ptr<HierarchicalOptimizer>;
} // namespace srrg2_hipe
