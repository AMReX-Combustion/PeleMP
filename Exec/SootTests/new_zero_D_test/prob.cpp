#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_N2 = 0.6693174609128006;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_O2 = 0.20333443771796805;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_Fuel = 0.12734810136923122;
AMREX_GPU_DEVICE_MANAGED int fuelIndx = -1;
} // namespace ProbParm

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
  pp.query("init_p", ProbParm::p0);
  pp.query("init_T", ProbParm::T0);
  pp.query("init_N2", ProbParm::Y_N2);
  pp.query("init_O2", ProbParm::Y_O2);
  pp.query("init_fuel", ProbParm::Y_Fuel);
  std::string fuel_name = "C2H4";
  pp.query("fuel_name", fuel_name);
  amrex::Vector<std::string> spec_names;
  EOS::speciesNames(spec_names);
  for (int sp = 0; sp != NUM_SPECIES; ++sp) {
    std::string spec_name = spec_names[sp];
    if (spec_name == fuel_name) ProbParm::fuelIndx = sp;
  }
  if (ProbParm::fuelIndx < 0) amrex::Abort("Fuel not found in chemistry mechanism");

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = ProbParm::Y_N2;
  massfrac[O2_ID] = ProbParm::Y_O2;
  massfrac[ProbParm::fuelIndx] = ProbParm::Y_Fuel;
  EOS::PYT2RE(ProbParm::p0, massfrac, ProbParm::T0, ProbParm::rho0, eint);
  EOS::RTY2Cs(ProbParm::rho0, ProbParm::T0, massfrac, cs);
  EOS::TY2Cp(ProbParm::T0, massfrac, cp);
}
}

void pc_prob_close()
{}
