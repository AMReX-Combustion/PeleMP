
#include "SprayParticles.H"
#include "SprayInjectTemplate.H"
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

int
interpolateInjectTime(
  const amrex::Real& time, const int nvals, const amrex::Real* inject_time)
{
  int i = 0;
  amrex::Real ctime = inject_time[i];
  while (ctime < time) {
    ctime = inject_time[++i];
  }
  return i;
}

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
  if (prob_parm.inject_N > 0) {
    const int time_indx = interpolateInjectTime(
      time, prob_parm.inject_N, prob_parm.inject_time.dataPtr());
    const amrex::Real time1 = prob_parm.inject_time[time_indx - 1];
    const amrex::Real time2 = prob_parm.inject_time[time_indx];
    const amrex::Real mf1 = prob_parm.inject_mass[time_indx - 1];
    const amrex::Real mf2 = prob_parm.inject_mass[time_indx];
    const amrex::Real invt = (time - time1) / (time2 - time1);
    mass_flow_rate = mf1 + (mf2 - mf1) * invt;
    const amrex::Real jv1 = prob_parm.inject_vel[time_indx - 1];
    const amrex::Real jv2 = prob_parm.inject_vel[time_indx];
    jet_vel = jv1 + (jv2 - jv1) * invt;
  }
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
  amrex::Real part_dia = prob_parm.part_mean_dia;
  amrex::Real part_stdev = prob_parm.part_stdev_dia;
  LogNormDist log_dist(part_dia, part_stdev);
  amrex::Real jet_angle = 0.;
  const int norm_dir = 1;
  sprayInjection<LogNormDist>(
    log_dist, prob_parm.jet_cent, jet_dia, prob_parm.part_temp, mass_flow_rate,
    jet_vel, prob_parm.spray_angle, dt, prob_parm.Y_jet.data(), lev, jet_angle,
    0., 360., norm_dir);
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
