
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

using namespace amrex;

bool
SprayParticleContainer::insertParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level,
                                        ProbParmHost const& prob_parm)
{
  return false;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level,
                                        ProbParmHost const& prob_parm)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(ProbParmHost const& prob_parm)
{
}
