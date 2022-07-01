
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
  amrex::Real jet_vel = prob_parm.jet_vel;
  const auto dx = this->m_gdb->Geom(lev).CellSize();
  // This absolutely must be included with any injection or insertion
  // function or significant issues will arise
  if (jet_vel * dt / dx[0] > m_partCFL) {
    amrex::Real max_vel = dx[0] * m_partCFL / dt;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::string warn_msg =
        "Injection velocity of " + std::to_string(jet_vel) +
        " is reduced to maximum " + std::to_string(max_vel);
      amrex::Warning(warn_msg);
    }
    m_injectVel = jet_vel;
    jet_vel = max_vel;
  }
  if (prob_parm.part_stdev_dia >= 0.) {
    LogNormDist log_dist(prob_parm.part_mean_dia, prob_parm.part_stdev_dia);
    sprayInjection<LogNormDist>(
      log_dist, prob_parm.jet_cent, prob_parm.jet_dia, prob_parm.part_temp,
      prob_parm.mass_flow_rate, jet_vel, prob_parm.spread_angle, dt,
      prob_parm.Y_jet.data(), lev, prob_parm.jet_angle, prob_parm.jet_azi_lo,
      prob_parm.jet_azi_hi);
  } else {
    WeibullDist w_dist(prob_parm.part_mean_dia, prob_parm.part_weibull_k);
    sprayInjection<WeibullDist>(
      w_dist, prob_parm.jet_cent, prob_parm.jet_dia, prob_parm.part_temp,
      prob_parm.mass_flow_rate, jet_vel, prob_parm.spread_angle, dt,
      prob_parm.Y_jet.data(), lev, prob_parm.jet_angle, prob_parm.jet_azi_lo,
      prob_parm.jet_azi_hi);
  }
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
