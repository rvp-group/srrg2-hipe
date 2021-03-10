#pragma once
#include "splitter_base.h"
#include "odometry_propagation.h"

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  class BreadthFirstSplitter : public SplitterBase {
  public:
    BreadthFirstSplitter();
    PARAM(PropertyInt,
          min_partition_diameter,
          "Minimim number of expansions to construct a partition",
          50,
          nullptr);

    void compute() final;

  protected:
    void _constructPartition(const VariableBase* root_, bool root_is_boundary);
    const VariableBase* _computeOuterDegreeMap();
    void _determineAnchorAndCreatePartition(const VariablePtrVector& variables_,
                                            const FactorPtrVector& factors,
                                            const VariableIdSet& boundary_variables_);
    std::list<VariableBase::Id> _open_variables;
    VariableIdSet _visited_variables;
    std::set<FactorBase::Id> _processed_factors;
    std::unordered_map<VariableBase::Id, int> _variable_outer_degree_map;
    OdometryPropagationPtr _init_action;
  };
} // namespace srrg2_hipe
