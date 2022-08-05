#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");
  std::string type;
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("gas_T", PeleLM::prob_parm->T_ox);
  pp.query("gas_jet_dia", PeleLM::prob_parm->gas_jet_dia);
  pp.query("gas_jet_vel", PeleLM::prob_parm->gas_jet_vel);
}
