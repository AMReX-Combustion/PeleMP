#include <PeleLM.H>
#include <pelelm_prob.H>
#include <pmf.H>

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo /*problo*/,
  const amrex_real* probhi /*probhi*/)
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("standoff", PeleLM::prob_parm->standoff);
  pp.query("pertmag", PeleLM::prob_parm->pertmag);

  std::string pmf_datafile;
  pp.query("pmf_datafile", pmf_datafile);
  int pmf_do_average = 1;
  PMF::read_pmf(pmf_datafile, pmf_do_average);

  pp.query("jet_vel", PeleLM::prob_parm->jet_vel);
  // The cells are divided by this value when prescribing the jet inlet
  pp.query("jet_dx_mod", PeleLM::prob_parm->jet_dx_mod);
  pp.get("jet_dia", PeleLM::prob_parm->jet_dia);
  pp.get("part_mean_dia", PeleLM::prob_parm->part_mean_dia);
  pp.query("part_stdev_dia", PeleLM::prob_parm->part_stdev_dia);
  pp.get("part_temp", PeleLM::prob_parm->part_temp);
  pp.query("mass_flow_rate", PeleLM::prob_parm->mass_flow_rate);
  pp.get("spray_angle_deg", PeleLM::prob_parm->spray_angle);
  int jets_per_dir = 0;
  pp.get("jets_per_dir", jets_per_dir);
  std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
  in_Y_jet[0] = 1.;
  pp.queryarr("jet_mass_fracs", in_Y_jet);
  amrex::Real sumY = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    PeleLM::prob_parm->Y_jet[spf] = in_Y_jet[spf];
    sumY += in_Y_jet[spf];
  }
  if (std::abs(sumY - 1.) > 1.E-8)
    amrex::Abort("'jet_mass_fracs' must sum to 1");
  // Convert to radians
  PeleLM::prob_parm->spray_angle *= M_PI / 180.;
  // Total number of jets
  // unsigned int total_jets = std::pow(jets_per_dir, AMREX_SPACEDIM - 1);
  unsigned int total_jets = 1;
  PeleLM::prob_parm->num_jets = total_jets;
  PeleLM::prob_parm->jet_cents.resize(total_jets);
  amrex::Real dom_len = probhi[0] - problo[0];
  amrex::Real div_len = dom_len / (amrex::Real(jets_per_dir));
  amrex::Real yloc = (probhi[1] - problo[1]) * .5;
  amrex::Real xloc = div_len * .5;
  amrex::Real zloc = problo[2];
  PeleLM::prob_parm->jet_cents[0] =
    amrex::RealVect(AMREX_D_DECL(xloc, yloc, zloc));
}
}
