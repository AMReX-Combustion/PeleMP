#include <PeleLM.H>
#include <pelelm_prob.H>

void
PeleLM::readProbParm()
{
  // Parse params
  PeleLM::prob_parm->norm_dir = AMREX_SPACEDIM - 1;
  amrex::ParmParse pp("prob");
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("fuel_temp", PeleLM::prob_parm->fuel_T);
  pp.query("oxid_temp", PeleLM::prob_parm->oxid_T);
  pp.query("wall_temp", PeleLM::prob_parm->wall_T);
  pp.query("ref_temp", PeleLM::prob_parm->T_ref);
  pp.query("init_vel", PeleLM::prob_parm->init_vel);
  pp.query("norm_dir", PeleLM::prob_parm->norm_dir);
  std::string fuel_name = "C2H4";
  std::string ign_name = "N-C7H16";
  pp.query("fuel_name", fuel_name);
  pp.query("ign_name", ign_name);
  pp.query("ign_Y", PeleLM::prob_parm->ign_Y);
  pp.query("fuel_dia", PeleLM::prob_parm->fuel_dia);
  pp.query("wall_thick", PeleLM::prob_parm->wall_thick);
  pp.query("oxid_dia", PeleLM::prob_parm->oxid_dia);
  pp.query("fuel_vel", PeleLM::prob_parm->fuel_vel);
  pp.query("oxid_vel", PeleLM::prob_parm->oxid_vel);
  pp.query("ext_vel", PeleLM::prob_parm->ext_vel);
  pp.query("oxid_vel_flat", PeleLM::prob_parm->oxid_flat);
  pp.query("fuel_vel_flat", PeleLM::prob_parm->fuel_flat);
  pp.query("vel_ramp_time", PeleLM::prob_parm->vel_time);
  pp.query("jet_center", PeleLM::prob_parm->jet_center);
  amrex::Vector<std::string> spec_names(NUM_SPECIES);
  pele::physics::eos::speciesNames<pele::physics::EosType>(spec_names);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    if (spec_names[n] == fuel_name) {
      PeleLM::prob_parm->fuel_indx = n;
    }
    if (spec_names[n] == ign_name) {
      PeleLM::prob_parm->ign_indx = n;
    }
  }
  if (PeleLM::prob_parm->fuel_indx < 0) {
    amrex::Abort("Fuel species not found from reaction data");
  }
  if (PeleLM::prob_parm->ign_indx < 0) {
    amrex::Abort("Ignition species not found from reaction data");
  }
  {
    int ecompCHON[NUM_SPECIES * 4];
    pele::physics::eos::element_compositionCHON<pele::physics::EosType>(ecompCHON);
    int numc = ecompCHON[PeleLM::prob_parm->fuel_indx * 4];
    int numh = ecompCHON[PeleLM::prob_parm->fuel_indx * 4 + 1];
    amrex::Real m = (amrex::Real)numc;
    amrex::Real n = (amrex::Real)numh;
    amrex::Real s = 32.*(m + n / 4.)/(12. * m + n);
    PeleLM::prob_parm->S_equil = s / 0.233;
  }
#ifdef PELELM_USE_SOOT
  if (PeleLM::do_soot_solve) {
    amrex::Real moments[NUM_SOOT_MOMENTS + 1];
    SootData* const sd = PeleLM::soot_model->getSootData();
    sd->initialSmallMomVals(moments);
    for (int n = 0; n < NUM_SOOT_MOMENTS + 1; ++n) {
      PeleLM::prob_parm->soot_vals[n] = moments[n];
    }
  } else {
    PeleLM::prob_parm->solve_soot = false;
  }
#endif
}
