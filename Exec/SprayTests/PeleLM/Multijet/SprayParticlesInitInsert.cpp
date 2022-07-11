
#include "SprayParticles.H"
#include "SprayInjectionTemplate.H"
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
  if (lev != 0) {
    return false;
  }
  if (time < prob_parm.jet_start_time || time > prob_parm.jet_end_time) {
    return false;
  }
  const amrex::Geometry& geom = this->m_gdb->Geom(lev);
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const auto dx = geom.CellSize();
  amrex::Vector<amrex::RealVect> jet_cents(prob_parm.num_jets);
  amrex::Real div_lenx =
    (phi[0] - plo[0]) / (amrex::Real(prob_parm.jets_per_dir[0]));
  int jetz = 1;
  amrex::Real div_lenz = 0.;
  amrex::Real zlo = 0.;
#if AMREX_SPACEDIM == 3
  div_lenz = (phi[2] - plo[2]) / (amrex::Real(prob_parm.jets_per_dir[2]));
  zlo = plo[2];
  jetz = prob_parm.jets_per_dir[2];
#endif
  amrex::Real yloc = plo[1];
  int jindx = 0;
  for (int i = 0; i < prob_parm.jets_per_dir[0]; ++i) {
    amrex::Real xloc = plo[0] + div_lenx * (amrex::Real(i) + 0.5);
    for (int k = 0; k < jetz; ++k) {
      amrex::Real zloc = zlo + div_lenz * (amrex::Real(k) + 0.5);
      jet_cents[jindx] = amrex::RealVect(AMREX_D_DECL(xloc, yloc, zloc));
      jindx++;
    }
  }
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
  int norm_dir = 1;
  for (int jindx = 0; jindx < prob_parm.num_jets; ++jindx) {
    amrex::RealVect cur_jet_cent = jet_cents[jindx];
    sprayInjection<LogNormDist>(
      log_dist, jet_dia, part_temp, mass_flow_rate, jet_vel,
      prob_parm.spray_angle, dt, prob_parm.Y_jet.data(), lev, 0., 0., 360.,
      norm_dir);
  }
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(ProbParm const& prob_parm)
{
  // This ensures the initial time step size stays reasonable
  m_injectVel = prob_parm.jet_vel;
  // Start without any particles
  return;
}
