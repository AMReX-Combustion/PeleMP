#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("T_init", PeleLM::prob_parm->T0);

  pp.get("jet_vel", PeleLM::prob_parm->jet_vel);
  std::vector<amrex::Real> jcent(AMREX_SPACEDIM);
  pp.getarr("jet_cent", jcent);
  std::vector<amrex::Real> jnorm(AMREX_SPACEDIM);
  pp.getarr("jet_norm", jnorm);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    PeleLM::prob_parm->jet_cent[dir] = jcent[dir];
    PeleLM::prob_parm->jet_norm[dir] = jnorm[dir];
  }
  PeleLM::prob_parm->num_inj_procs = 1;
  // Number of processors to do the injection
  pp.query("num_inj_procs", PeleLM::prob_parm->num_inj_procs);
  // The cells are divided by this value when prescribing the jet inlet
  pp.get("jet_dia", PeleLM::prob_parm->jet_dia);
  pp.get("part_mean_dia", PeleLM::prob_parm->part_mean_dia);
  pp.query("part_stdev_dia", PeleLM::prob_parm->part_stdev_dia);
  if (PeleLM::prob_parm->part_stdev_dia < 0.) {
    // If no standard deviation is specified, assume we are using Weibull distribution
    if (!pp.contains("part_weibull_k")) {
      amrex::Abort("Must specify either standard deviation or Weibull k value");
    } else {
      pp.get("part_weibull_k", PeleLM::prob_parm->part_weibull_k);
    }
  }
  pp.get("part_temp", PeleLM::prob_parm->part_temp);
  pp.query("mass_flow_rate", PeleLM::prob_parm->mass_flow_rate);
  // All angles must be in degrees
  // This is the spreading angle
  pp.get("spread_angle", PeleLM::prob_parm->spread_angle);
  std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
  in_Y_jet[0] = 1.;
  pp.queryarr("Y_jet", in_Y_jet);
  amrex::Real sumY = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    PeleLM::prob_parm->Y_jet[spf] = in_Y_jet[spf];
    sumY += in_Y_jet[spf];
  }
  if (std::abs(sumY - 1.) > 1.E-8) {
    amrex::Abort("'jet_mass_fracs' must sum to 1");
  }
  // Convert to radians
  PeleLM::prob_parm->spread_angle *= M_PI / 180.;
}
