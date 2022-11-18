#include <PeleLM.H>
#include <pelelm_prob.H>

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("init_T", PeleLM::prob_parm->T0);
  pp.query("init_vel", PeleLM::prob_parm->vel);
  pp.query("init_N2", PeleLM::prob_parm->Y_N2);
  pp.query("init_O2", PeleLM::prob_parm->Y_O2);
}
}
