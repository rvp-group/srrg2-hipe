#include "chordal_initializer.h"
#include <srrg_solver/solver_core/internals/sparse_block_matrix/matrix_block_factory.h>
#include <srrg_solver/solver_core/iteration_algorithm_lm.h>
#include <srrg_solver/solver_core/solver.h>
#include <srrg_solver/variables_and_factors/types_3d/se3_pose_pose_chordal_hessian_factor.h>
#include <srrg_solver/variables_and_factors/types_3d/se3_pose_pose_geodesic_error_factor.h>
#include <srrg_solver/variables_and_factors/types_3d/variable_se3.h>

namespace srrg2_hipe {
  using namespace srrg2_core;
  using namespace srrg2_solver;

  void initializeChordal(FactorGraphInterfacePtr graph,
                         std::set<VariableBase::Id>* selected_variables,
                         const int& level_,
                         const bool do_rotation,
                         const bool do_translation) {
    FactorGraphPtr rotation_graph    = std::make_shared<FactorGraph>();
    FactorGraphPtr translation_graph = std::make_shared<FactorGraph>();

    MatrixBlockFactory* factory = MatrixBlockFactory::instance();
    factory->addAllocator<9, 9>();
    factory->addAllocator<9, 1>();
    factory->addAllocator<3, 3>();
    factory->addAllocator<3, 1>();
    factory->addAllocator<4, 4>();
    factory->addAllocator<4, 1>();
    factory->addAllocator<2, 2>();
    factory->addAllocator<2, 1>();

    const auto& factors = graph->factors();
    std::set<VariableBase::Id> pose_variables;
    for (const auto& id_factor : factors) {
      const FactorBase* f      = id_factor.second;
      if (f->level() != level_) {
        continue;
      }

      const SE3PosePoseChordalHessianFactor* se3 =
        dynamic_cast<const SE3PosePoseChordalHessianFactor*>(f);
      if (se3) {
        SE3ChordalInitializationRotationErrorFactor* chord_rotation_factor =
          new SE3ChordalInitializationRotationErrorFactor;
        chord_rotation_factor->setMeasurement(se3->measurement());
        chord_rotation_factor->setVariableId(0, f->variableId(0));
        chord_rotation_factor->setVariableId(1, f->variableId(1));
        chord_rotation_factor->setEnabled(f->enabled());
        rotation_graph->addFactor(FactorBasePtr(chord_rotation_factor));

        SE3ChordalInitializationTranslationErrorFactor* chord_translation_factor =
          new SE3ChordalInitializationTranslationErrorFactor;
        chord_translation_factor->setMeasurement(se3->measurement());
        chord_translation_factor->setVariableId(0, f->variableId(0));
        chord_translation_factor->setVariableId(1, f->variableId(1));
        chord_translation_factor->setEnabled(f->enabled());
        translation_graph->addFactor(FactorBasePtr(chord_translation_factor));

        pose_variables.insert(f->variableId(0));
        pose_variables.insert(f->variableId(1));
        continue;
      }

      const SE3PosePoseGeodesicErrorFactor* se3_geo =
        dynamic_cast<const SE3PosePoseGeodesicErrorFactor*>(f);
      if (se3_geo) {
        SE3ChordalInitializationRotationErrorFactor* chord_rotation_factor =
          new SE3ChordalInitializationRotationErrorFactor;
        chord_rotation_factor->setMeasurement(se3_geo->measurement());
        chord_rotation_factor->setVariableId(0, f->variableId(0));
        chord_rotation_factor->setVariableId(1, f->variableId(1));
        chord_rotation_factor->setEnabled(f->enabled());
        rotation_graph->addFactor(FactorBasePtr(chord_rotation_factor));

        SE3ChordalInitializationTranslationErrorFactor* chord_translation_factor =
          new SE3ChordalInitializationTranslationErrorFactor;
        chord_translation_factor->setMeasurement(se3_geo->measurement());
        chord_translation_factor->setVariableId(0, f->variableId(0));
        chord_translation_factor->setVariableId(1, f->variableId(1));
        chord_translation_factor->setEnabled(f->enabled());
        translation_graph->addFactor(FactorBasePtr(chord_translation_factor));

        pose_variables.insert(f->variableId(0));
        pose_variables.insert(f->variableId(1));
        continue;
      }
    }

    for (const VariableBase::Id& id : pose_variables) {
      const VariableBase* v = graph->variable(id);
      const VariableSE3Base* se3 = dynamic_cast<const VariableSE3Base*>(v);
      if (se3) {
        VariableSE3ChordalInitializationRotation* chord_rotation_variable =
          new VariableSE3ChordalInitializationRotation;
        chord_rotation_variable->setGraphId(v->graphId());
        chord_rotation_variable->setEstimate(se3->estimate());
        chord_rotation_variable->setStatus(v->status());
        rotation_graph->addVariable(VariableBasePtr(chord_rotation_variable));

        VariableSE3ChordalInitializationTranslation* chord_translation_variable =
          new VariableSE3ChordalInitializationTranslation;
        chord_translation_variable->setGraphId(v->graphId());
        chord_translation_variable->setEstimate(se3->estimate());
        chord_translation_variable->setStatus(v->status());
        translation_graph->addVariable(VariableBasePtr(chord_translation_variable));

        continue;
      }
    }

    Solver solver;
    solver.param_max_iterations.pushBack(1);
    solver.param_algorithm.setValue(IterationAlgorithmBasePtr(new IterationAlgorithmLM));
    // std::shared_ptr<SimpleTerminationCriteria> term(new SimpleTerminationCriteria);
    // term->param_epsilon.setValue(1e-1);
    // solver.param_termination_criteria.setValue(term);
    solver.setGraph(rotation_graph);
    if (do_rotation) {
      solver.compute();
      std::cerr << "chordalInitializer|Rotation" << solver.iterationStats() << std::endl;
    }
    for (const VariableBase::Id& id : pose_variables) {
      VariableBase* v      = graph->variable(id);
      VariableSE3Base* se3 = dynamic_cast<VariableSE3Base*>(v);
      if (se3) {
        const VariableSE3ChordalInitializationRotation* chord_rotation_variable =
          dynamic_cast<const VariableSE3ChordalInitializationRotation*>(
            rotation_graph->variable(id));
        VariableSE3ChordalInitializationTranslation* chord_translation_variable =
          dynamic_cast<VariableSE3ChordalInitializationTranslation*>(
            translation_graph->variable(id));
        chord_translation_variable->setEstimate(chord_rotation_variable->estimate());
        continue;
      }
    }
    solver.setGraph(translation_graph);
    if (do_translation) {
      solver.compute();
      std::cerr << "chordalInitializer|Translation" << solver.iterationStats() << std::endl;
    }
    for (const VariableBase::Id& id : pose_variables) {
      VariableBase* v      = graph->variable(id);
      VariableSE3Base* se3 = dynamic_cast<VariableSE3Base*>(v);
      if (se3) {
        VariableSE3ChordalInitializationTranslation* chord_variable =
          dynamic_cast<VariableSE3ChordalInitializationTranslation*>(
            translation_graph->variable(id));
        se3->setEstimate(chord_variable->estimate());
        continue;
      }
    }
    // tg copy to output set
    if (selected_variables) {
      selected_variables->clear();
      selected_variables->insert(pose_variables.begin(), pose_variables.end());
    }
  }
} // namespace srrg2_hipe
