#include "instances.h"
#include "partition.h"
#include "partition_creation_utils.h"
#include "partition_manager.h"

namespace srrg2_hipe {
  void partition_registerTypes() {
    BOSS_REGISTER_CLASS(PartitionManager)
    BOSS_REGISTER_CLASS(SE3PosePoseVirtualFactorCreator)
    BOSS_REGISTER_CLASS(SE3PosePoseChordalVirtualFactorCreator)
  }
} // namespace srrg2_hipe
