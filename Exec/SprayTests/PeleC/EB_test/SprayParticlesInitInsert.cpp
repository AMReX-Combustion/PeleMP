
#include "SprayParticles.H"
#include "SprayInjectTemplate.H"
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
  amrex::ignore_unused(nstep, finest_level, prob_parm_d);
  if (lev != 0) {
    return false;
  }
  SprayJets* js = m_sprayJets[0].get();

  sprayInjection(js, dt, lev);

  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts,
  ProbParmHost const& prob_parm,
  ProbParmDevice const& prob_parm_d)
{
  amrex::ignore_unused(prob_parm, init_parts);
  m_sprayJets.resize(1);
  std::string jet_name = "jet1";
  bool do_inject = false;
  ParmParse ps("spray.jet1");
  pp.query("do_inject", do_inject);
  if (do_inject) {
    m_sprayJets[0] = std::make_unique<SprayJet>(jet_name, Geom(0));
  }
  return;
}
