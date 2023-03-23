
#include "SprayParticles.H"
#include "SprayInjection.H"
#include "pelelm_prob.H"

bool
SprayParticleContainer::injectParticles(
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  int /*nstep*/,
  int /*lev*/,
  int /*finest_level*/,
  ProbParm const& /*prob_parm*/)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts, ProbParm const& prob_parm)
{
  // Check if we are simply restarting
  if (!init_parts) {
    return;
  }
  const int level = 0;
  amrex::ParmParse pp("prob");
  amrex::IntVect partNum(AMREX_D_DECL(100, 100, 100));
  pp.query("num_particles", partNum);
  amrex::RealVect partVel = amrex::RealVect::TheZeroVector();
  amrex::Real partDia = PeleLM::prob_parm->partDia;
  amrex::Real partTemp = PeleLM::prob_parm->partTemp;
  std::array<amrex::Real, SPRAY_FUEL_NUM> partY = {0.0};
  int numRedist = 1;
  partY[0] = 1.;
  uniformSprayInit(
    partNum, partVel, partDia, partTemp, partY.begin(), level, numRedist);
}
