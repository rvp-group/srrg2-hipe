#include <partitions/partition_manager.h>
#include <srrg_geometry/geometry3d.h>
#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>

using namespace srrg2_core;
using namespace srrg2_solver;
using namespace srrg2_hipe;

using VariableType = VariableSE3QuaternionRight;
using FactorType   = SE3PosePoseGeodesicErrorFactor;

void printPartition(PartitionPtr partition);

int main(int argc, char** argv) {
  int pose_for_partition = 3;

  PartitionManager manager;
  FactorGraphPtr graph(new FactorGraph);
  manager.setGraph(graph);
  manager.param_solver.value()->param_max_iterations.pushBack(1);
  manager.param_covariance_projectors.pushBack(
    VirtualFactorCreatorBasePtr(new SE3PosePoseVirtualFactorCreator));
  VariablePtrVector variables;
  FactorPtrVector factors;
  // first variable
  VariableType* v = new VariableType;
  v->setEstimate(Isometry3f::Identity());
  graph->addVariable(VariableBasePtr(v));
  variables.push_back(v);
  FactorType* f  = nullptr;
  Vector6f delta = Vector6f::Zero();
  // create data for first partition
  for (int i = 1; i < pose_for_partition; ++i) {
    delta.setRandom();
    f            = new FactorType;
    Isometry3f z = geometry3d::ta2t(delta);
    f->setMeasurement(z);
    f->setVariableId(0, v->graphId());
    Isometry3f from = v->estimate();
    Isometry3f to   = from * z;
    v               = new VariableType;
    v->setEstimate(to);
    graph->addVariable(VariableBasePtr(v));
    f->setVariableId(1, v->graphId());
    graph->addFactor(FactorBasePtr(f));
    variables.push_back(v);
    factors.push_back(f);
  }
  // create partition
  VariableBase::Id anchor_id_1 = variables.at(1)->graphId();
  VariableIdSet boundary;
  boundary.insert(variables.at(2)->graphId());
  manager.createPartition(variables, factors, anchor_id_1, boundary);
  printPartition(manager.partition(anchor_id_1));
  // data for the second partition
  variables.clear();
  factors.clear();
  variables.push_back(v);
  for (int i = 1; i < pose_for_partition; ++i) {
    delta.setRandom();
    f            = new FactorType;
    Isometry3f z = geometry3d::ta2t(delta);
    f->setMeasurement(z);
    f->setVariableId(0, v->graphId());
    Isometry3f from = v->estimate();
    Isometry3f to   = from * z;
    v               = new VariableType;
    v->setEstimate(to);
    graph->addVariable(VariableBasePtr(v));
    f->setVariableId(1, v->graphId());
    graph->addFactor(FactorBasePtr(f));
    variables.push_back(v);
    factors.push_back(f);
  }
  // create partition
  VariableBase::Id anchor_id_2 = variables.at(1)->graphId();
  manager.createPartition(variables, factors, anchor_id_2, boundary);
  printPartition(manager.partition(anchor_id_2));
  // manager.mergePartitions(anchor_id_1, anchor_id_2, anchor_id_1);
  printPartition(manager.partition(anchor_id_1));
  manager.destroyPartition(anchor_id_1);
  for (const auto& var : graph->variables()) {
    if (manager.variableConnectivity(var.second->graphId())) {
      std::cerr << "There is still some variables in the partitions" << std::endl;
    }
  }
  return 0;
}

void printPartition(PartitionPtr partition) {
  std::cerr << "Partition " << partition->anchorId() << std::endl;
  std::cerr << "Variables : " << std::endl;
  for (const auto& v : partition->variables()) {
    std::cerr << v.first << " " << v.second << std::endl;
  }
  std::cerr << "Factors : " << std::endl;
  for (const auto& f : partition->factors()) {
    std::cerr << f.first << " " << f.second << std::endl;
  }
  std::cerr << "Boundary variables id : " << std::endl;
  for (const auto& id : partition->boundaryVariableIds()) {
    std::cerr << id << std::endl;
  }
  std::cerr << "Virtual factors : " << std::endl;
  for (const auto& f : partition->virtualFactors()) {
    std::cerr << f->graphId() << " " << f << std::endl;
  }
}
