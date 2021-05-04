
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

bool
SprayParticleContainer::injectParticles(
  Real time,
  Real dt,
  int nstep,
  int lev,
  int finest_level,
  ProbParmHost const& prob_parm,
  ProbParmDevice const& prob_parm_d)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(
  ProbParmHost const& prob_parm, ProbParmDevice const& prob_parm_d)
{
}
