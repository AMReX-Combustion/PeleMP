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
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("standoff", PeleLM::prob_parm->standoff);
  std::string pmf_datafile;
  pp.query("pmf_datafile", pmf_datafile);
  int pmf_do_average = 0;
  pp.query("pmf_average", pmf_do_average);
  PMF::read_pmf(pmf_datafile, pmf_do_average);
  amrex::Real moments[NUM_SOOT_MOMENTS + 1] = {0.0};
  if (PeleLM::do_soot_solve) {
    SootData* const sd = PeleLM::soot_model->getSootData();
    sd->initialSmallMomVals(moments);
  }
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; ++n) {
    PeleLM::prob_parm->soot_vals[n] = moments[n];
  }
  PmfData* pmf_data = pmf_data_g;
  const int pmfN = pmf_data->pmf_N;
  for (int i = 0; i < pmfN; ++i) {
    amrex::Real sumY = 0.;
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    for (int n = 0; n < NUM_SPECIES; ++n) {
      const int j = pmfN * (n + 3) + i;
      massfrac[n] = amrex::max(0., amrex::min(1., pmf_data->pmf_Y[j]));
      sumY += massfrac[n];
    }
    for (int n = 0; n < NUM_SPECIES; ++n) {
      const int j = pmfN * (n + 3) + i;
      pmf_data->pmf_Y[j] = massfrac[n] / sumY;
    }
  }
}
}
