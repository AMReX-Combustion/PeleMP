
#include "SprayParticles.H"
#include "SprayInjection.H"
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
  amrex::ignore_unused(nstep, finest_level, prob_parm);
  if (lev != 0) {
    return false;
  }
  SprayJet* js = m_sprayJets[0].get();
  if (!js->jet_active(time)) {
    return false;
  }
  sprayInjection(time, js, dt, lev);
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts, ProbParm const& prob_parm)
{
  m_sprayJets.resize(1);
  std::string jet_name = "jet1";
  m_sprayJets[0] = std::make_unique<SprayJet>(jet_name, Geom(0));
  // Start without any particles
  m_injectVel = m_sprayJets[0]->jet_vel();
  if (init_parts) {
    const auto dx = this->m_gdb->Geom(0).CellSize();
    amrex::Real fakedt = dx[0] * m_partCFL / m_injectVel;
    amrex::Real time = 0.;
    sprayInjection(time, m_sprayJets[0].get(), fakedt, 0);
    Redistribute();
  }
  return;
}
