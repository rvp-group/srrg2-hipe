#pragma once
#include "partition.h"
#include "partition_creation_utils.h"
#include <srrg_config/configurable.h>
#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/solver_core/solver.h>
#include <unordered_map>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  class PartitionManager : public Configurable {
  public:
    // partitions are identified by the variable graph id of their anchor
    using IdAnchorPartitionMap = std::unordered_map<VariableBase::Id, PartitionPtr>;
    // get all the partitions that include a certain variable
    using VariablePartitionMultiMap = std::multimap<VariableBase::Id, PartitionPtr>;

    PARAM(PropertyConfigurable_<Solver>,
          solver,
          "Non linear solver to optimize the partitions",
          SolverPtr(new Solver),
          nullptr);

    PARAM(PropertyBool, verbose, "Show statistics of the partition optimization", true, nullptr);

    PARAM_VECTOR(PropertyConfigurableVector_<VirtualFactorCreatorBase>,
                 covariance_projectors,
                 "Projectors for the covariance of the virtual measurements",
                 nullptr);

    // Access partition using its anchor id
    PartitionPtr partition(const VariableBase::Id& partition_anchor_id_);
    // Get variable connectivity in the current partitions configuration
    size_t variableConnectivity(const VariableBase::Id& variable_id_);
    // Given a variable id takes all the partitions that contain that variable
    bool getConnectedPartitions(PartitionPtrVector& partitions_,
                                const VariableBase::Id& variable_id_);
    // add an already created partition to this manager
    void addPartition(PartitionPtr partition);
    // create partition out of variables and factors
    // the anchor variable MUST be contained in variables
    bool createPartition(const VariablePtrVector& variables,
                         const FactorPtrVector& factors,
                         const VariableBase::Id& anchor_id,
                         const VariableIdSet& boundary_variables_ids);
    // add new boundary variables in an existing partition
    void addBoundaryVariablesToPartition(PartitionPtr partition_,
                                         const VariableIdSet& boundary_variables_ids_);
    // optimize partition_ (might be external to this manager)
    void optimizePartition(PartitionPtr partition_);
    // optimize partition_ identified in this manager by an anchor_id_
    void optimizePartition(const VariableBase::Id& anchor_id_);
    // cancel partition_ from the internal structure, the virtual factors will be canceled
    // in the original graph
    bool destroyPartition(const VariableBase::Id& anchor_id_);
    // set factor graph, the partitions will refer to this graph
    // and the virtual factors will be added to it
    inline void setGraph(const FactorGraphPtr& graph_) {
      _graph = graph_;
    }
    // get underlying graph (cannot add remove variable if you access the graph in this way)
    inline FactorGraphInterfacePtr graph() {
      return _graph;
    }
    // get partitions created until now
    inline IdAnchorPartitionMap& partitions() {
      return _partitions;
    }
    // same but read only
    inline const IdAnchorPartitionMap& partitions() const {
      return _partitions;
    }

  protected:
    // Compute virtual factors for the variables in the pairs
    // using the covariance blocks
    void _computeVirtualFactors(PartitionPtr& partition_,
                                VariablePairVector& variable_pairs_,
                                MatrixBlockVector& covariances_);
    FactorGraphPtr _graph = nullptr;
    IdAnchorPartitionMap _partitions;
    VariablePartitionMultiMap _variables_to_partitions;
  };

  using PartitionManagerPtr = std::shared_ptr<PartitionManager>;
} // namespace srrg2_hipe
