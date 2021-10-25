#include "instances.h"
#include "se3_chordal_factor.h"
#include "se3_langevin_factor.h"
using namespace srrg2_solver;
void se3FactorsType() {
  BOSS_REGISTER_CLASS(SE3ChordalFactor);
  BOSS_REGISTER_CLASS(SE3LangevinFactor);
}
