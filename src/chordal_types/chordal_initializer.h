#pragma once
#include "factors_chordal_initialization.h"
#include "variables_chordal_initialization.h"
#include <srrg_solver/solver_core/factor_graph.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  void initializeChordal(FactorGraphInterfacePtr graph_,
                         std::set<VariableBase::Id>* selected_variables = nullptr,
                         const int& level_                              = 0,
                         const bool do_rotation                         = true,
                         const bool do_translation                      = true);
} // namespace srrg2_hipe
