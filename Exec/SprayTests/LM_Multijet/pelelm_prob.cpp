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
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("ref_T", PeleLM::prob_parm->T0);
  pp.query("init_vel", PeleLM::prob_parm->vel);
  pp.query("init_N2", PeleLM::prob_parm->Y_N2);
  pp.query("init_O2", PeleLM::prob_parm->Y_O2);
  pp.query("jet_vel", PeleLM::prob_parm->jet_vel);
  pp.query("jet_start_time", PeleLM::prob_parm->jet_start_time);
  pp.query("jet_end_time", PeleLM::prob_parm->jet_end_time);
  // The cells are divided by this value when prescribing the jet inlet
  pp.query("jet_dx_mod", PeleLM::prob_parm->jet_dx_mod);
  pp.get("jet_dia", PeleLM::prob_parm->jet_dia);
  pp.get("part_mean_dia", PeleLM::prob_parm->part_mean_dia);
  pp.query("part_stdev_dia", PeleLM::prob_parm->part_stdev_dia);
  pp.get("part_temp", PeleLM::prob_parm->part_temp);
  pp.query("mass_flow_rate", PeleLM::prob_parm->mass_flow_rate);
  pp.get("spray_angle_deg", PeleLM::prob_parm->spray_angle);
  std::vector<int> jets_per_dir(AMREX_SPACEDIM);
  pp.getarr("jets_per_dir", jets_per_dir);
  std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
  in_Y_jet[0] = 1.;
  pp.queryarr("jet_mass_fracs", in_Y_jet);
  amrex::Real sumY = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    PeleLM::prob_parm->Y_jet[spf] = in_Y_jet[spf];
    sumY += in_Y_jet[spf];
  }
  if (std::abs(sumY - 1.) > 1.E-8) {
    amrex::Abort("'jet_mass_fracs' must sum to 1");
  }
  // Convert to radians
  PeleLM::prob_parm->spray_angle *= M_PI / 180.;
  // Total number of jets
  unsigned int total_jets = AMREX_D_TERM(jets_per_dir[0], *1, *jets_per_dir[2]);
  PeleLM::prob_parm->num_jets = total_jets;
  PeleLM::prob_parm->jet_cents.resize(total_jets);
  amrex::Real div_lenx =
    (probhi[0] - problo[0]) / (amrex::Real(jets_per_dir[0]));
  int jetz = 1;
  amrex::Real div_lenz = 0.;
#if AMREX_SPACEDIM == 3
  div_lenz = (probhi[2] - problo[2]) / (amrex::Real(jets_per_dir[2]));
  jetz = jets_per_dir[2];
#endif
  amrex::Real yloc = problo[1];
  int jindx = 0;
  for (int i = 0; i < jets_per_dir[0]; ++i) {
    amrex::Real xloc = div_lenx * (amrex::Real(i) + 0.5);
    for (int k = 0; k < jetz; ++k) {
      amrex::Real zloc = div_lenz * (amrex::Real(k) + 0.5);
      PeleLM::prob_parm->jet_cents[jindx] =
        amrex::RealVect(AMREX_D_DECL(xloc, yloc, zloc));
      jindx++;
    }
  }
}
}
