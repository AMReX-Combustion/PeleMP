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
  pp.query("init_v", PeleC::prob_parm_device->v0);
  pp.query("ref_p", PeleC::prob_parm_device->p0);
  pp.query("ref_T", PeleC::prob_parm_device->T0);
  pp.query("init_N2", PeleC::prob_parm_device->Y_N2);
  pp.query("init_O2", PeleC::prob_parm_device->Y_O2);

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = PeleC::prob_parm_device->Y_N2;
  massfrac[O2_ID] = PeleC::prob_parm_device->Y_O2;
  EOS::PYT2RE(PeleC::prob_parm_device->p0, massfrac, PeleC::prob_parm_device->T0,
              PeleC::prob_parm_device->rho0, eint);
  EOS::RTY2Cs(PeleC::prob_parm_device->rho0, PeleC::prob_parm_device->T0, massfrac, cs);
  EOS::TY2Cp(PeleC::prob_parm_device->T0, massfrac, cp);
}
}

void pc_prob_close()
{}


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
