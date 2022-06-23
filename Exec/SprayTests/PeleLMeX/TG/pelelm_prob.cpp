#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");

   std::string type;
   amrex::Real Stmod = 5.;
   amrex::Real rhoRatio = 1000.;
   amrex::Real mach = 0.05;
   pp.query("P_mean",   PeleLM::prob_parm->P_mean);
   pp.query("init_T", PeleLM::prob_parm->T0);
   pp.query("reynolds", PeleLM::prob_parm->reynolds);
   pp.query("mach", mach);
   pp.query("num_particles", PeleLM::prob_parm->partNum);
   pp.query("density_ratio", rhoRatio);
   pp.query("st_mod", Stmod);
   if (mach > 0.3) {
     amrex::Abort("Mach number should be less than 0.3");
   }
   const amrex::Real Pr = 0.71;
   const amrex::Real refL = geom[0].ProbHi(0) - geom[0].ProbLo(0);
   amrex::Real pres = PeleLM::prob_parm->P_mean;
   amrex::Real pres_cgs = pres * 10.;
   amrex::Real T = PeleLM::prob_parm->T0;
   amrex::Real reyn = PeleLM::prob_parm->reynolds;
   amrex::Real massfrac[NUM_SPECIES] = {0.0};
   massfrac[N2_ID] = 1.;
   amrex::Real rho, cs, cp;
   auto eos = pele::physics::PhysicsType::eos();
   eos.PYT2R(pres_cgs, massfrac, T, rho);
   eos.RTY2Cs(rho, T, massfrac, cs);
   eos.TY2Cp(T, massfrac, cp);
   rho *= 1.E3;
   cp *= 1.E-4;
   cs *= 0.01;
   amrex::Real Ugas = mach * cs;
   PeleLM::prob_parm->Ugas = Ugas;
   auto& ltransparm = PeleLM::trans_parms.host_trans_parm();
   amrex::Real mu = rho * Ugas * refL / reyn;
   amrex::Real lambda = mu * cp / Pr;
   ltransparm.const_bulk_viscosity = 0.;
   ltransparm.const_diffusivity = 0.;
   ltransparm.const_viscosity = mu * 10.;
   ltransparm.const_conductivity = lambda * 1.E5;
   PeleLM::trans_parms.sync_to_device();
   amrex::Real St_num = Stmod / (8. * M_PI);
   amrex::Real partRho;
   amrex::ParmParse ppp("particles");
   ppp.get("fuel_rho", partRho);
   bool mass_trans = false;
   ppp.get("mass_transfer", mass_trans);
   if (mass_trans) {
     amrex::Abort("particles.mass_transfer must be turned off");
   }
   if (std::abs(rhoRatio - partRho / rho) > 10.) {
     amrex::Real desrho = rho * rhoRatio;
     std::string errorstr = "Error: restart solution with particles.fuel_rho = " + std::to_string(rho * rhoRatio) +
       " to achieve desired density ratio";
     amrex::Abort(errorstr);
   }
   amrex::Real tau_g = refL / Ugas;
   amrex::Real tau_d = St_num * tau_g;
   amrex::Real dia = std::sqrt(18. * mu * tau_d / partRho);
   amrex::Real Re_d = dia * Ugas * rho / mu;
   amrex::Real error = 1000.;
   const amrex::Real tol = 1.E-6;
     int k = 0;
  const int maxIter = 500;
  // Iteratively solve for the diameter of the spray droplets
  while (error > tol) {
    amrex::Real oldDia = dia;
    Re_d = dia * Ugas * rho / mu;
    amrex::Real C_D = 24. / Re_d;
    if (Re_d > 1.) {
      C_D *= (1. + std::pow(Re_d, 2. / 3.) / 6.);
    }
    dia = 0.75 * rho * C_D * Ugas * tau_d / partRho;
    error = std::abs(oldDia - dia) / dia;
    k++;
    if (k > maxIter) {
      amrex::Abort("Failed to converge a particle diameter");
    }
  }
  PeleLM::prob_parm->partDia = dia;
  PeleLM::prob_parm->partTemp = T;
  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "number of particles: "
                      << PeleLM::prob_parm->partNum[0];
    for (int dir = 1; dir != AMREX_SPACEDIM; ++dir) {
      amrex::Print() << ", " << PeleLM::prob_parm->partNum[dir];
    }
    amrex::Print() << std::endl;
    amrex::Print() << "rho0: " << rho << std::endl;
    amrex::Print() << "cs: " << cs << std::endl;
    amrex::Print() << "U: " << Ugas << std::endl;
    amrex::Print() << "mu: " << mu << std::endl;
    amrex::Print() << "Re: " << reyn << std::endl;
    amrex::Print() << "Stokes number: " << Stmod << "*Stc " << std::endl;
    amrex::Print() << "particle diameter: " << dia << std::endl;
    amrex::Print() << "tau_d: " << tau_d << std::endl;
    amrex::Print() << "Re_d: " << Re_d << std::endl;
    amrex::Print() << "final time (72 tau_g): " << 72. * tau_g << std::endl;
  }
}
