
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

  sprayInjection(js, dt, lev, 0, 1, 5);
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts, ProbParm const& prob_parm)
{
  amrex::ignore_unused(init_parts, prob_parm);
  m_sprayJets.resize(1);
  std::string jet_name = "jet1";
  m_sprayJets[0] = std::make_unique<SprayJet>(jet_name);
  // Start without any particles
  m_injectVel = m_sprayJets[0]->jet_vel();
  if (init_parts) {
    const auto dx = this->m_gdb->Geom(0).CellSize();
    amrex::Real fakedt = dx[0] * m_partCFL / m_injectVel;
    sprayInjection(m_sprayJets[0].get(), fakedt, 0);
  }
  return;
}
