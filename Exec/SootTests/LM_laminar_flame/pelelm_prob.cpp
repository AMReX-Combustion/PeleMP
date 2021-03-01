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
  pp.query("standoff", PeleLM::prob_parm->standoff);
  std::string pmf_datafile;
  pp.query("pmf_datafile", pmf_datafile);
  int pmf_do_average = 1;
  PMF::read_pmf(pmf_datafile, pmf_do_average);
}
}
