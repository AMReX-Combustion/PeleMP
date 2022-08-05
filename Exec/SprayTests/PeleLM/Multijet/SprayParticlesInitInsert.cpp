
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
  amrex::ignore_unused(nstep, finest_level);
  if (lev != 0) {
    return false;
  }
  for (int jindx = 0; jindx < m_sprayJets.size(); ++jindx) {
    SprayJets* js = m_sprayJets[jindx].get();
    if (js->jet_active(time)) {
      sprayInjection(js, dt, 0);
    }
  }
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts, ProbParm const& prob_parm)
{
  amrex::ProbParm ps("spray");
  amrex::Real jet_vel = 50.;
  amrex::Real jet_dia = 1.E-4;
  amrex::Real part_mean_dia = 1.E-5;
  amrex::Real part_stdev_dia = 0.;
  amrex::Real mass_flow_rate = 2.3E-3;
  amrex::Real part_temp = 300.;
  amrex::Real jet_start_time = 0.;
  amrex::Real jet_end_time = 10000.;
  amrex::Real spray_angle = 20.;
  amrex::Real[SPRAY_FUEL_NUM] Y_jet = {0.0};
  ps.query("jet_vel", jet_vel);
  ps.query("jet_start_time", jet_start_time);
  ps.query("jet_end_time", jet_end_time);
  // The cells are divided by this value when prescribing the jet inlet
  ps.get("jet_dia", jet_dia);
  ps.get("part_mean_dia", part_mean_dia);
  ps.query("part_stdev_dia", part_stdev_dia);
  ps.get("part_temp", part_temp);
  ps.query("mass_flow_rate", mass_flow_rate);
  ps.get("spray_angle_deg", spray_angle);
  std::vector<int> jets_per_dir(AMREX_SPACEDIM);
  ps.getarr("jets_per_dir", jets_per_dir);
  std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
  in_Y_jet[0] = 1.;
  ps.queryarr("jet_mass_fracs", in_Y_jet);
  amrex::Real sumY = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_jet[spf] = in_Y_jet[spf];
    sumY += in_Y_jet[spf];
  }
  if (std::abs(sumY - 1.) > 1.E-8) {
    amrex::Abort("'jet_mass_fracs' must sum to 1");
  }
  // Total number of jets
  int num_jets = AMREX_D_TERM(jets_per_dir[0], *1, *jets_per_dir[2]);
  m_sprayJets.resize(num_jets);
  const amrex::Geometry& geom = this->m_gdb->Geom(lev);
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const auto dx = geom.CellSize();
  amrex::Vector<amrex::RealVect> jet_cents(num_jets);
  amrex::Real div_lenx = (phi[0] - plo[0]) / (amrex::Real(jets_per_dir[0]));
  int jetz = 1;
  amrex::Real div_lenz = 0.;
  amrex::Real zlo = 0.;
#if AMREX_SPACEDIM == 3
  div_lenz = (phi[2] - plo[2]) / (amrex::Real(jets_per_dir[2]));
  zlo = plo[2];
  jetz = jets_per_dir[2];
#endif
  amrex::Real yloc = plo[1];
  int jindx = 0;
  amrex::RealVect jet_norm(AMREX_D_DECL(0., 1., 0.));
  std::string dist_type = "Uniform";
  for (int i = 0; i < prob_parm.jets_per_dir[0]; ++i) {
    amrex::Real xloc = plo[0] + div_lenx * (amrex::Real(i) + 0.5);
    for (int k = 0; k < jetz; ++k) {
      amrex::Real zloc = zlo + div_lenz * (amrex::Real(k) + 0.5);
      amrex::RealVect jet_cent(AMREX_D_DECL(xloc, yloc, zloc));
      m_sprayJets[jindx] = std::make_unique<SprayJet>(
        jet_cent, jet_norm, spray_angle, jet_dia, jet_vel, mass_flow_rate,
        part_temp, Y_jet, dist_type, jet_start_time, jet_end_time);
      jindx++;
    }
  }
  if (init_parts) {
    const auto dx = this->m_gdb->Geom(0).CellSize();
    amrex::Real fakedt = dx[0] * m_partCFL / jet_vel;
    for (int jindx = 0; jindx < prob_parm.num_jets; ++jindx) {
      SprayJets* js = m_sprayJets[jiindx].get();
      if (js->jet_active(0.)) {
        sprayInjection(js, fakedt, 0);
      }
    }
  }
  return;
}
