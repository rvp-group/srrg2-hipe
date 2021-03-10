#pragma once
#include <partitions/partition.h>
#include <srrg_config/configurable.h>
#include <srrg_solver/solver_core/factor_graph.h>
#include <suitesparse/cholmod.h>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace srrg2_hipe {

  using namespace srrg2_core;
  using namespace srrg2_solver;

  class IsoperimetricBipartitioner : public Configurable {
  public:
    using GraphBipartition = std::pair<PartitionPtr, PartitionPtr>;
    // @brief ctor
    IsoperimetricBipartitioner();
    // @brief dtor
    ~IsoperimetricBipartitioner();
    // @brief set graph to split
    void setGraph(FactorGraphInterfacePtr graph_) {
      _graph = graph_;
    }
    // @brief check if a variable is on the left or right side of the partition
    // if is not an active variable it assume that is the gauge of the partition
    // does return true
    inline bool isVariableOnLeftSideOfThePartition(const VariableBase::Id& id_) const {
      bool flag = true;
      auto it   = _variable_id_to_index.find(id_);
      if (it != _variable_id_to_index.end()) {
        flag = _partition_vector[it->second] < _threshold;
      }
      return flag;
    }
    // @brief compute edge weights, in case of multiple partitions this have to be called
    // just for the original graph
    void computeEdgeWeights();

    // @brief get the partition vector
    inline const std::vector<double>& partitionVector() const {
      return _partition_vector;
    }
    // @brief compute partition
    GraphBipartition compute();

  protected:
    // @brief custom types
    using LaplacianEntry       = std::tuple<int, int, double>;
    using LaplacianEntryVector = std::vector<LaplacianEntry>;
    using VariableIdIntMap     = std::unordered_map<VariableBase::Id, int>;
    using VariableIdDoubleMap      = std::unordered_map<VariableBase::Id, double>;
    using VariableVectorWeightPair = std::pair<std::vector<VariableBase::Id>, double>;
    using FactorWeightMap          = std::unordered_map<FactorBase::Id, VariableVectorWeightPair>;

    void _constructLaplacian();
    void _buildLinearSystem();
    void _solveLinearSystem();
    void _computeCutThreshold();
    void _addVariables(PartitionPtr& partition_,
                       std::set<VariableBase::Id>& variable_set_,
                       const std::set<VariableBase::Id>& shared_variables_);
    // @brief graph to be splitted
    FactorGraphInterfacePtr _graph = nullptr;
    // @brief partition threshold
    double _threshold = 0;
    // @brief partition vector
    std::vector<double> _partition_vector;
    // @brief laplacian entries
    LaplacianEntryVector _laplacian;
    // @brief variable id to matrix index (required since i don't want to make this class
    // friend with VariableBase)
    VariableIdIntMap _variable_id_to_index;
    // @brief variable to degree of a vertex, used for construction the laplacian and
    // select the shared variables from the edge separetor
    VariableIdDoubleMap _variable_id_degree;
    // @brief store factor weight, as the likelihood
    FactorWeightMap _factor_to_weight;
    // @brief cholmod members
    cholmod_common _cholmodCommon;
    cholmod_sparse* _A = nullptr;
    cholmod_dense* _b  = nullptr;
    cholmod_dense* _x  = nullptr;
    cholmod_factor* _L = nullptr;
  };
} // namespace srrg2_hipe
