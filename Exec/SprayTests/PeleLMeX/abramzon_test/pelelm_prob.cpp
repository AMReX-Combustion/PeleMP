#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");
  auto eos = pele::physics::PhysicsType::eos();
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("T_gas", PeleLM::prob_parm->T0);
  pp.query("vel_gas", PeleLM::prob_parm->vel_gas);
  pp.query("init_N2", PeleLM::prob_parm->Y_N2);
  pp.query("init_O2", PeleLM::prob_parm->Y_O2);
  // Droplet initial values
  amrex::Real Reyn = -1.; // Reynolds number based on diameter
  pp.query("re_d", Reyn);
  pp.get("dia_drop", PeleLM::prob_parm->dia_drop);
  pp.get("T_drop", PeleLM::prob_parm->T_drop);
  std::vector<amrex::Real> vel_drop(AMREX_SPACEDIM, 0.);
  pp.queryarr("vel_drop", vel_drop);
  if (Reyn > 0. && pp.contains("vel_drop")) {
    amrex::Abort("Cannot specify droplet velocity and Reynolds number");
  }
  std::vector<amrex::Real> loc_drop(AMREX_SPACEDIM, 0.);
  pp.queryarr("loc_drop", loc_drop);
  if (pp.contains("loc_drop")) {
    PeleLM::prob_parm->set_drop_loc = true;
  }
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    PeleLM::prob_parm->vel_drop[dir] = vel_drop[dir];
    PeleLM::prob_parm->loc_drop[dir] = loc_drop[dir];
  }
  if (SPRAY_FUEL_NUM > 1) {
    std::vector<amrex::Real> Y_drop(SPRAY_FUEL_NUM, 0.);
    pp.getarr("Y_drop", Y_drop);
    amrex::Real sumtest = 0.;
    for (int n = 0; n < SPRAY_FUEL_NUM; ++n) {
      PeleLM::prob_parm->Y_drop[n] = Y_drop[n];
      sumtest += Y_drop[n];
    }
    if (std::abs(1. - sumtest) > 1.E-6) {
      amrex::Abort("Liquid mass fractions must sum to 1!");
    }
  } else {
    PeleLM::prob_parm->Y_drop[0] = 1.;
  }
  if (Reyn > 0.) {
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[N2_ID] = PeleLM::prob_parm->Y_N2;
    massfrac[O2_ID] = PeleLM::prob_parm->Y_O2;
    amrex::Real Ddiag[NUM_SPECIES] = {0.0};
    amrex::Real temp = PeleLM::prob_parm->T0;
    amrex::Real pres = PeleLM::prob_parm->P_mean * 10.;
    amrex::Real rho_cgs = 0.;
    eos.PYT2R(pres, massfrac, temp, rho_cgs);
    auto const* ltransparm = PeleLM::trans_parms.device_trans_parm();
    auto trans = pele::physics::PhysicsType::transport();
    amrex::Real mu, lambda, xi;
    trans.transport(
      false, true, false, false, false, temp, rho_cgs, massfrac, Ddiag, nullptr, mu, xi, lambda,
      ltransparm);
    amrex::Real dia_cgs = PeleLM::prob_parm->dia_drop * 100.;
    amrex::Real umax = mu * Reyn / (rho_cgs * dia_cgs) * 0.01;
    PeleLM::prob_parm->vel_gas = umax;
  }
}
