
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

using namespace amrex;

bool
SprayParticleContainer::insertParticles (Real time, Real dt, int nstep, int lev)
{
  return false;
}

bool
SprayParticleContainer::injectParticles (Real time, Real dt, int nstep, int lev)
{
  return false;
} 

void
SprayParticleContainer::InitSprayParticles(AmrLevel* pelec, const int& lev, const int& num_ppc)
{
}
