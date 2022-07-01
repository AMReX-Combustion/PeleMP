
#include "SprayParticles.H"
#include "SprayInjectTemplate.H"
#include "pelelm_prob.H"

bool
SprayParticleContainer::injectParticles(
  amrex::Real time,
  amrex::Real dt,
  int nstep,
  int lev,
  int finest_level,
  ProbParm const& prob_parm)
{
  amrex::ignore_unused(nstep, finest_level);
  if (lev != 0) {
    return false;
  }
  if (time < prob_parm.jet_start_time || time > prob_parm.jet_end_time) {
    return false;
  }
  const amrex::Geometry& geom = this->m_gdb->Geom(lev);
  const auto dx = geom.CellSize();
  amrex::Real mass_flow_rate = prob_parm.mass_flow_rate;
  amrex::Real jet_vel = prob_parm.jet_vel;
  amrex::Real jet_dia = prob_parm.jet_dia;
  amrex::Real part_temp = prob_parm.part_temp;
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
  amrex::Real part_dia = prob_parm.part_mean_dia;
  amrex::Real part_stdev = prob_parm.part_stdev_dia;
  LogNormDist log_dist(part_dia, part_stdev);
  sprayInjection<LogNormDist>(
    log_dist, prob_parm.jet_cent, jet_dia, part_temp, mass_flow_rate, jet_vel,
    prob_parm.spray_angle, dt, prob_parm.Y_jet.data(), lev);
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(ProbParm const& prob_parm)
{
  // Start without any particles
  m_injectVel = prob_parm.jet_vel;
  return;
}
