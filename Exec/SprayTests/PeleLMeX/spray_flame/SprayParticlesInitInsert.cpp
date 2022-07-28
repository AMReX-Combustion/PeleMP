
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
  SprayJet* js = m_sprayJets[0].get();
  if (lev != 0) {
    return false;
  }
  if (!js->jet_active(time)) {
    return false;
  }
  amrex::Real orig_vel = js->jet_vel();
  amrex::Real jet_vel = orig_vel;
  const auto dx = this->m_gdb->Geom(lev).CellSize();
  // This absolutely must be included with any injection or insertion
  // function or significant issues will arise
  bool vel_changed = false;
  if (jet_vel * dt / dx[0] > m_partCFL) {
    amrex::Real max_vel = dx[0] * m_partCFL / dt;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::string warn_msg =
        "Injection velocity of " + std::to_string(jet_vel) +
        " is reduced to maximum " + std::to_string(max_vel);
      amrex::Warning(warn_msg);
    }
    m_injectVel = jet_vel;
    js->set_jet_vel(jet_vel);
    vel_changed = true;
  }

  sprayInjection(js, dt, lev, 0, 1, 5);

  // Be sure to reset jet velocity if it was changed
  if (vel_changed) {
    js->set_jet_vel(orig_vel);
  }
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  bool init_parts, ProbParm const& prob_parm)
{
  m_sprayJets.resize(1);
  std::string jet_name = "jet1";
  m_sprayJets[0] = std::make_unique<SprayJet>(jet_name);
  // Start without any particles
  m_injectVel = m_sprayJets[0]->jet_vel();
  return;
}
