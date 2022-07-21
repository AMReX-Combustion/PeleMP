
#include "SprayParticles.H"
#include "SprayInjectionTemplate.H"
#include <PeleC.H>
#include "prob.H"

bool
SprayParticleContainer::injectParticles(
  amrex::Real time,
  amrex::Real dt,
  int nstep,
  int lev,
  int finest_level,
  ProbParmHost const& prob_parm,
  ProbParmDevice const& prob_parm_d)
{
  if (lev != 0) {
    return false;
  }
  if (time < prob_parm.jet_start_time || time > prob_parm.jet_end_time) {
    return false;
  }
  amrex::ignore_unused(nstep, finest_level, prob_parm_d);
  amrex::Real mass_flow_rate = prob_parm.mass_flow_rate;
  amrex::Real jet_vel = prob_parm.jet_vel;
  const amrex::Geometry& geom = this->m_gdb->Geom(lev);
  const auto dx = geom.CellSize();
  amrex::Real jet_dia = prob_parm.jet_dia;
  // This absolutely must be included with any injection or insertion
  // function or significant issues will arise
  if (jet_vel * dt / dx[0] > 0.5) {
    amrex::Real max_vel = dx[0] * 0.5 / dt;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::string warn_msg =
        "Injection velocity of " + std::to_string(jet_vel) +
        " is reduced to maximum " + std::to_string(max_vel);
      amrex::Warning(warn_msg);
    }
    m_injectVel = jet_vel;
    jet_vel = max_vel;
  }
  int curProc = amrex::ParallelDescripter::MyProc();
  amrex::Real part_dia = prob_parm.part_mean_dia;
  amrex::Real part_stdev = prob_parm.part_stdev_dia;
  LogNormDist log_dist(part_dia, part_stdev);
  const int norm_dir = 1;
  amrex::RealVect jet_norm(amrex::RealVect::TheZeroVector());
  jet_norm[norm_dir] = 1.;
  for (int jindx = 0; jindx < prob_parm.num_jets; ++jindx) {
    if (curProc == jindx) {
      amrex::RealVect cur_jet_cent = prob_parm.jet_cents[jindx];
      sprayInjection<LogNormDist>(
        log_dist, cur_jet_cent, jet_norm, jet_dia, prob_parm.part_temp, mass_flow_rate,
        jet_vel, prob_parm.spray_angle, dt, prob_parm.Y_jet.data(), lev, curProc);
    }
  }
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  ProbParmHost const& prob_parm, ProbParmDevice const& prob_parm_d)
{
  // This ensures the initial time step size stays reasonable
  m_injectVel = prob_parm.jet_vel;
  // Start without any particles
  return;
}
