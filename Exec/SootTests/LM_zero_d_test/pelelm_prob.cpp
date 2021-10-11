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

  amrex::Real moments[NUM_SOOT_MOMENTS + 1];
  SootData* const sd = PeleLM::soot_model->getSootData();
  sd->initialSmallMomVals(moments);
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; ++n) {
    PeleLM::prob_parm->soot_vals[n] = moments[n];
  }
  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::EosType>(spec_names);
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    std::string spec_name = spec_names[sp];
    if (spec_name == fuel_name)
      PeleLM::prob_parm->fuelIndx = sp;
  }
  if (PeleLM::prob_parm->fuelIndx < 0)
    amrex::Abort("Fuel not found in chemistry mechanism");
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real air_O2 = 0.233;
  amrex::Real air_N2 = 0.767;
  amrex::Real S_equil = 0.;
  {
    int ecompCHON[NUM_SPECIES * 4];
    pele::physics::eos::element_compositionCHON<pele::physics::EosType>(
      ecompCHON);
    int numc = ecompCHON[PeleLM::prob_parm->fuelIndx * 4];
    int numh = ecompCHON[PeleLM::prob_parm->fuelIndx * 4 + 1];
    amrex::Real m = (Real)numc;
    amrex::Real n = (Real)numh;
    amrex::Real s = 32. * (m + n / 4.) / (12. * m + n);
    S_equil = s / air_O2;
  }
  if (pp.contains("mixture_fraction")) {
    amrex::Real Z = 0.;
    pp.get("mixture_fraction", Z);
    PeleLM::prob_parm->Y_Fuel = Z;
    PeleLM::prob_parm->Y_O2 = air_O2 * (1. - Z);
    PeleLM::prob_parm->Y_N2 = (1. - air_O2) * (1. - Z);
  }
}
}
