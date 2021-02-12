#include "prob.H"

void
pc_prob_close()
{
}

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
  pp.query("init_p", PeleC::prob_parm_device->p0);
  pp.query("init_T", PeleC::prob_parm_device->T0);
  pp.query("init_N2", PeleC::prob_parm_device->Y_N2);
  pp.query("init_O2", PeleC::prob_parm_device->Y_O2);
  pp.query("init_fuel", PeleC::prob_parm_device->Y_Fuel);
  std::string fuel_name = "C2H4";
  pp.query("fuel_name", fuel_name);
  amrex::Vector<std::string> spec_names;
  EOS::speciesNames(spec_names);
  for (int sp = 0; sp != NUM_SPECIES; ++sp) {
    std::string spec_name = spec_names[sp];
    if (spec_name == fuel_name)
      PeleC::prob_parm_device->fuelIndx = sp;
  }
  if (PeleC::prob_parm_device->fuelIndx < 0)
    amrex::Abort("Fuel not found in chemistry mechanism");

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = PeleC::prob_parm_device->Y_N2;
  massfrac[O2_ID] = PeleC::prob_parm_device->Y_O2;
  massfrac[PeleC::prob_parm_device->fuelIndx] = PeleC::prob_parm_device->Y_Fuel;
  EOS::PYT2RE(
    PeleC::prob_parm_device->p0, massfrac, PeleC::prob_parm_device->T0,
    PeleC::prob_parm_device->rho0, eint);
  EOS::RTY2Cs(
    PeleC::prob_parm_device->rho0, PeleC::prob_parm_device->T0, massfrac, cs);
  EOS::TY2Cp(PeleC::prob_parm_device->T0, massfrac, cp);
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
