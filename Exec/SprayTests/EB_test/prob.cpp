#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("init_v", PeleC::h_prob_parm_device->v0);
  pp.query("init_p", PeleC::h_prob_parm_device->p0);
  pp.query("init_T", PeleC::h_prob_parm_device->T0);
  pp.query("init_N2", PeleC::h_prob_parm_device->Y_N2);
  pp.query("init_O2", PeleC::h_prob_parm_device->Y_O2);
  pp.query("jet_vel", PeleC::prob_parm_host->jet_vel);
  pp.get("jet_dia", PeleC::prob_parm_host->jet_dia);
  // The cells are divided by this value when prescribing the jet inlet
  pp.query("jet_dx_mod", PeleC::prob_parm_host->jet_dx_mod);
  pp.get("part_mean_dia", PeleC::prob_parm_host->part_mean_dia);
  pp.query("part_stdev_dia", PeleC::prob_parm_host->part_stdev_dia);
  pp.get("part_temp", PeleC::prob_parm_host->part_temp);
  pp.query("mass_flow_rate", PeleC::prob_parm_host->mass_flow_rate);
  pp.get("spray_angle_deg", PeleC::prob_parm_host->spray_angle);
  pp.query("jet_start_time", PeleC::prob_parm_host->jet_start_time);
  pp.query("jet_end_time", PeleC::prob_parm_host->jet_end_time);
  // Convert to radians
  PeleC::prob_parm_host->spray_angle *= M_PI / 180.;

  AMREX_D_TERM(
    PeleC::prob_parm_host->jet_cent[0] = 0.5 * (probhi[0] + problo[0]);
    , PeleC::prob_parm_host->jet_cent[1] = problo[1];
    , PeleC::prob_parm_host->jet_cent[2] = 0.5 * (probhi[2] + problo[2]););

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
    PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac, cs);
  eos.TY2Cp(PeleC::h_prob_parm_device->T0, massfrac, cp);
}
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
