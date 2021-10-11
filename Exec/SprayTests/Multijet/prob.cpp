#include "prob.H"

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("init_v", PeleC::h_prob_parm_device->v0);
  pp.get("ref_p", PeleC::h_prob_parm_device->p0);
  pp.get("ref_T", PeleC::h_prob_parm_device->T0);
  pp.query("init_N2", PeleC::h_prob_parm_device->Y_N2);
  pp.query("init_O2", PeleC::h_prob_parm_device->Y_O2);
  pp.query("jet_vel", PeleC::prob_parm_host->jet_vel);
  pp.query("jet_start_time", PeleC::prob_parm_host->jet_start_time);
  pp.query("jet_end_time", PeleC::prob_parm_host->jet_end_time);
  // The cells are divided by this value when prescribing the jet inlet
  pp.query("jet_dx_mod", PeleC::prob_parm_host->jet_dx_mod);
  pp.get("jet_dia", PeleC::prob_parm_host->jet_dia);
  pp.get("part_mean_dia", PeleC::prob_parm_host->part_mean_dia);
  pp.query("part_stdev_dia", PeleC::prob_parm_host->part_stdev_dia);
  pp.get("part_temp", PeleC::prob_parm_host->part_temp);
  pp.query("mass_flow_rate", PeleC::prob_parm_host->mass_flow_rate);
  pp.get("spray_angle_deg", PeleC::prob_parm_host->spray_angle);
  std::vector<int> jets_per_dir(AMREX_SPACEDIM);
  pp.getarr("jets_per_dir", jets_per_dir);
  std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
  in_Y_jet[0] = 1.;
  pp.queryarr("jet_mass_fracs", in_Y_jet);
  amrex::Real sumY = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    PeleC::prob_parm_host->Y_jet[spf] = in_Y_jet[spf];
    sumY += in_Y_jet[spf];
  }
  if (std::abs(sumY - 1.) > 1.E-8)
    amrex::Abort("'jet_mass_fracs' must sum to 1");
  // Convert to radians
  PeleC::prob_parm_host->spray_angle *= M_PI / 180.;
  // Total number of jets
  unsigned int total_jets = AMREX_DTERM(jets_per_dir[0], , *jets_per_dir[2]);
  PeleC::prob_parm_host->num_jets = total_jets;
  PeleC::prob_parm_host->jet_cents.resize(total_jets);
  amrex::Real div_lenx =
    (probhi[0] - problo[0]) / (amrex::Real(jets_per_dir[0]));
  int jetz = 1;
  amrex::Real div_lenz = 0.;
#if AMREX_SPACEDIM == 3
  div_lenz = (probhi[2] - problo[2]) / (amrex::Real(jets_per_dir[2]));
  jetz = jets_per_dir[2];
#endif
  amrex::Real yloc = problo[1];
  int jindx = 0;
  for (int i = 0; i < jets_per_dir[0]; ++i) {
    amrex::Real xloc = div_lenx * (amrex::Real(i) + 0.5);
    for (int k = 0; k < jetz; ++k) {
      amrex::Real zloc = div_lenz * (amrex::Real(k) + 0.5);
      PeleC::prob_parm_host->jet_cents[jindx] =
        amrex::RealVect(AMREX_D_DECL(xloc, yloc, zloc));
      jindx++;
    }
  }
  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = PeleC::h_prob_parm_device->Y_N2;
  massfrac[O2_ID] = PeleC::h_prob_parm_device->Y_O2;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p0, massfrac, PeleC::h_prob_parm_device->T0,
    PeleC::h_prob_parm_device->rho0, eint);
  eos.RTY2Cs(
    PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac,
    cs);
  eos.TY2Cp(PeleC::h_prob_parm_device->T0, massfrac, cp);
}
}

void
pc_prob_close()
{
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
