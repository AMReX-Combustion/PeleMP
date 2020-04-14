
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

using namespace amrex;

bool
SprayParticleContainer::insertParticles (Real time, int nstep, int lev)
{
  return false;
}

bool
SprayParticleContainer::injectParticles (Real time, int nstep, int lev)
{
  return false;
} 

void
SprayParticleContainer::InitParticlesUniform(AmrLevel* pelec, const int& lev, const int& num_ppc)
{
}
