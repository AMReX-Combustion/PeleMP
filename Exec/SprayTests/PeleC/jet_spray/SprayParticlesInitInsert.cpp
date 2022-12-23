
#include "SprayParticles.H"
#include "SprayInjection.H"
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
  SprayJet* js = m_sprayJets[0].get();
  if (prob_parm.inject_N > 0) {
    const int time_indx = interpolateInjectTime(
      time, prob_parm.inject_N, prob_parm.inject_time.dataPtr());
    const amrex::Real time1 = prob_parm.inject_time[time_indx - 1];
    const amrex::Real time2 = prob_parm.inject_time[time_indx];
    const amrex::Real mf1 = prob_parm.inject_mass[time_indx - 1];
    const amrex::Real mf2 = prob_parm.inject_mass[time_indx];
    const amrex::Real invt = (time - time1) / (time2 - time1);
    const amrex::Real mass_flow_rate = mf1 + (mf2 - mf1) * invt;
    js->set_mass_flow(mass_flow_rate);
    const amrex::Real jv1 = prob_parm.inject_vel[time_indx - 1];
    const amrex::Real jv2 = prob_parm.inject_vel[time_indx];
    const amrex::Real jet_vel = jv1 + (jv2 - jv1) * invt;
    js->set_jet_vel(jet_vel);
  }

  sprayInjection(time, js, dt, lev);

  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts,
  ProbParmHost const& prob_parm,
  ProbParmDevice const& prob_parm_d)
{
  amrex::ignore_unused(prob_parm_d, init_parts);
  m_sprayJets.resize(1);
  std::string jet_name = "jet1";
  m_sprayJets[0] = std::make_unique<SprayJet>(jet_name, Geom(0));
  if (prob_parm.inject_N <= 0) {
    m_injectVel = m_sprayJets[0]->jet_vel();
  } else {
    amrex::Real start_time = prob_parm.inject_time[0];
    amrex::Real end_time = prob_parm.inject_time[prob_parm.inject_N - 1];
    m_sprayJets[0]->set_start_time(start_time);
    m_sprayJets[0]->set_end_time(end_time);
  }
  return;
}
