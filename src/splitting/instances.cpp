#include "instances.h"
#include "breadth_first_splitter.h"
#include "nested_dissection_splitter.h"

namespace srrg2_hipe {
  void splitting_registerTypes() {
    BOSS_REGISTER_CLASS(NestedDissectionSplitter)
    BOSS_REGISTER_CLASS(BreadthFirstSplitter)
  }
} // namespace srrg2_hipe
