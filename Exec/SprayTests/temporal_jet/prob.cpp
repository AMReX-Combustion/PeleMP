#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 900.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v0 = 1.;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_O2 = 0.233;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_N2 = 0.767;
AMREX_GPU_DEVICE_MANAGED int partNum = 8000;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_jet = 0.01;
AMREX_GPU_DEVICE_MANAGED amrex::Real partTemp = 300.;
AMREX_GPU_DEVICE_MANAGED amrex::Real partRho = 0.64;
AMREX_GPU_DEVICE_MANAGED amrex::Real partDia = 0.002;
AMREX_GPU_DEVICE_MANAGED amrex::Real partVel = 0.;
AMREX_GPU_DEVICE_MANAGED amrex::Real velFluct = 1.E-2;
} // namespace ProbParm

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
  pp.query("ref_p", ProbParm::p0);
  pp.query("ref_T", ProbParm::T0);
  pp.query("init_O2", ProbParm::Y_O2);
  pp.query("init_N2", ProbParm::Y_N2);
  pp.query("jet_len", ProbParm::L_jet);
  pp.query("part_num", ProbParm::partNum);
  pp.query("part_vel", ProbParm::partVel);
  pp.query("vel_fluct", ProbParm::velFluct);
  pp.get("part_rho", ProbParm::partRho);
  pp.get("part_dia", ProbParm::partDia);
  pp.get("part_temp", ProbParm::partTemp);

  // Initial density, velocity, and material properties
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = ProbParm::Y_O2;
  massfrac[N2_ID] = ProbParm::Y_N2;

}
}

