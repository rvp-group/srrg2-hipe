#pragma once
#include "splitter_base.h"

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  class NestedDissectionSplitter : public SplitterBase {
  public:
    using VariableIdIntMap       = std::unordered_map<VariableBase::Id, int>;
    using IntVariableIdMap       = std::unordered_map<int, VariableBase::Id>;
    using ComponentIdToPartition = std::map<int, PartitionPtr>;
    NestedDissectionSplitter();
    ~NestedDissectionSplitter();
    void compute() final;

  protected:
    void _fillVariableIndeces();
    void _buildApproximateHessianPattern();
    void _collapseSeparatorTreeNodes(ComponentIdToPartition& stars, int* parents);
    void _determineAnchorAndAddToManager(PartitionPtr& star_);
    VariableIdIntMap _variable_to_index;
    IntVariableIdMap _index_to_variable;
    cholmod_common _cholmodCommon;
    cholmod_sparse* _A = nullptr;
  };

  using NestedDissectionSplitterPtr = std::shared_ptr<NestedDissectionSplitter>;
} // namespace srrg2_hipe
