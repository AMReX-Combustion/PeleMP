#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   std::string type;
   pp.query("P_mean",   PeleLM::prob_parm->P_mean);
   pp.query("jet_vel", PeleLM::prob_parm->jet_vel);
   // The cells are divided by this value when prescribing the jet inlet
   pp.query("jet_dx_mod", PeleLM::prob_parm->jet_dx_mod);
   pp.get("jet_dia", PeleLM::prob_parm->jet_dia);
   pp.get("part_mean_dia", PeleLM::prob_parm->part_mean_dia);
   pp.query("part_stdev_dia", PeleLM::prob_parm->part_stdev_dia);
   pp.get("part_temp", PeleLM::prob_parm->part_temp);
   pp.query("mass_flow_rate", PeleLM::prob_parm->mass_flow_rate);
   pp.get("spray_angle_deg", PeleLM::prob_parm->spray_angle);
   pp.query("gas_jet_dia", PeleLM::prob_parm->gas_jet_dia);
   pp.query("gas_jet_vel", PeleLM::prob_parm->gas_jet_vel);
   std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
   in_Y_jet[0] = 1.;
   pp.queryarr("jet_mass_fracs", in_Y_jet);
   amrex::Real sumY = 0.;
   for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
     PeleLM::prob_parm->Y_jet[spf] = in_Y_jet[spf];
     sumY += in_Y_jet[spf];
   }
   // Convert to radians
   PeleLM::prob_parm->spray_angle *= M_PI / 180.;
   amrex::Real dom_len = geom[0].ProbHi(0) - geom[0].ProbLo(0);
   amrex::Real yloc = (geom[0].ProbHi(1) - geom[0].ProbLo(1)) * .5;
   amrex::Real xloc = dom_len * .5;
   amrex::Real zloc = geom[0].ProbLo(2);
   AMREX_D_TERM(PeleLM::prob_parm->jet_cents[0] = xloc;,
                PeleLM::prob_parm->jet_cents[1] = yloc;,
                PeleLM::prob_parm->jet_cents[2] = zloc;)
}
