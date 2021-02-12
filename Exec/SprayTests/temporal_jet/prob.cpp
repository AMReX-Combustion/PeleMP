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
  pp.query("ref_p", PeleC::prob_parm_device->p0);
  pp.query("ref_T", PeleC::prob_parm_device->T0);
  pp.query("init_O2", PeleC::prob_parm_device->Y_O2);
  pp.query("init_N2", PeleC::prob_parm_device->Y_N2);
  pp.query("jet_len", PeleC::prob_parm_device->L_jet);

  // Determine how smooth the velocity profile for the gas phase
  pp.query("vel_smoothing", PeleC::prob_parm_host->velSmooth);

  // Initial density, velocity, and material properties
  PeleC::prob_parm_device->v0 = PeleC::prob_parm_host->partVel;
  pp.query("part_num", PeleC::prob_parm_host->partNum);
  pp.query("part_vel", PeleC::prob_parm_host->partVel);
  pp.query("vel_fluct", PeleC::prob_parm_host->velFluct);
  pp.get("part_dia", PeleC::prob_parm_host->partDia);
  pp.get("part_temp", PeleC::prob_parm_host->partTemp);
}
}
