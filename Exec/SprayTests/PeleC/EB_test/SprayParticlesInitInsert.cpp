
#include "SprayParticles.H"
#include "SprayInjectTemplate.H"
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
  amrex::ignore_unused(nstep, finest_level, prob_parm_d);
  if (lev != 0) {
    return false;
  }
  SprayJets* js = m_sprayJets[0].get();
  if (prob_parm.inject_N > 0) {
    const int time_indx = interpolateInjectTime(
      time, prob_parm.inject_N, prob_parm.inject_time.dataPtr());
    const amrex::Real time1 = prob_parm.inject_time[time_indx - 1];
    const amrex::Real time2 = prob_parm.inject_time[time_indx];
    const amrex::Real mf1 = prob_parm.inject_mass[time_indx - 1];
    const amrex::Real mf2 = prob_parm.inject_mass[time_indx];
    const amrex::Real invt = (time - time1) / (time2 - time1);
    const amrex::Real mass_flow_rate = mf1 + (mf2 - mf1) * invt;
    js->set_mass_rate(mass_flow_rate);
    const amrex::Real jv1 = prob_parm.inject_vel[time_indx - 1];
    const amrex::Real jv2 = prob_parm.inject_vel[time_indx];
    const amrex::Real jet_vel = jv1 + (jv2 - jv1) * invt;
    js->set_jet_vel(jet_vel);
  }

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
    m_sprayJets[0] = std::make_unique<SprayJet>(jet_name);
  }
  return;
}
