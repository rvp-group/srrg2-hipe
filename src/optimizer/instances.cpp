#include "instances.h"
#include "hierarchical_optimizer.h"

namespace srrg2_hipe {
  using namespace srrg2_core;
  void hierarchicalOptimizer_registerTypes() {
    BOSS_REGISTER_CLASS(HierarchicalOptimizer)
  }
} // namespace srrg2_hipe
