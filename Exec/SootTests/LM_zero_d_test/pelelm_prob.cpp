#include <PeleLM.H>
#include <pelelm_prob.H>

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
  pp.query("ref_p", PeleLM::prob_parm->P_mean);
  pp.query("init_T", PeleLM::prob_parm->T0);
  pp.query("init_vel", PeleLM::prob_parm->init_vel);
  pp.query("init_N2", PeleLM::prob_parm->Y_N2);
  pp.query("init_O2", PeleLM::prob_parm->Y_O2);
  pp.query("init_fuel", PeleLM::prob_parm->Y_Fuel);
  std::string fuel_name = "C2H4";
  pp.query("fuel_name", fuel_name);
  amrex::Vector<std::string> spec_names;
  EOS::speciesNames(spec_names);
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    std::string spec_name = spec_names[sp];
    if (spec_name == fuel_name)
      PeleLM::prob_parm->fuelIndx = sp;
  }
  if (PeleLM::prob_parm->fuelIndx < 0)
    amrex::Abort("Fuel not found in chemistry mechanism");
}
}
