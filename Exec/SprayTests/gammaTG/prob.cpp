#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds = 1600.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real mach = 0.1;
AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl = 0.71;
AMREX_GPU_DEVICE_MANAGED amrex::Real Stmod = 5.;
AMREX_GPU_DEVICE_MANAGED bool convecting = false;
AMREX_GPU_DEVICE_MANAGED amrex::Real L = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::IntVect partNum = amrex::IntVect(AMREX_D_DECL(100, 100, 100));
AMREX_GPU_DEVICE_MANAGED amrex::Real rhoRatio = 1000.;
AMREX_GPU_DEVICE_MANAGED amrex::Real partRho = 1.;
AMREX_GPU_DEVICE_MANAGED amrex::Real partDia = 1.E-3;
} // namespace ProbParm

void
pc_prob_close()
{
}

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
  pp.query("reynolds", ProbParm::reynolds);
  pp.query("mach", ProbParm::mach);
  pp.query("convecting", ProbParm::convecting);
  pp.query("ref_p", ProbParm::p0);
  pp.query("ref_T", ProbParm::T0);
  pp.query("st_mod", ProbParm::Stmod);
  pp.query("num_particles", ProbParm::partNum);
  pp.query("density_ratio", ProbParm::rhoRatio);

  // Define the length scale
  ProbParm::L = probhi[0] - problo[0];

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(ProbParm::p0, massfrac, ProbParm::T0, ProbParm::rho0, eint);
  EOS::RTY2Cs(ProbParm::rho0, ProbParm::T0, massfrac, cs);
  EOS::TY2Cp(ProbParm::T0, massfrac, cp);

  amrex::Real refL = ProbParm::L;
  ProbParm::v0 = ProbParm::mach * cs;
  transport_params::const_bulk_viscosity = 0.0;
  transport_params::const_diffusivity = 0.0;
  const amrex::Real mu =
    ProbParm::rho0 * ProbParm::v0 * refL / ProbParm::reynolds;
  transport_params::const_viscosity = mu;
  transport_params::const_conductivity = mu * cp / ProbParm::prandtl;

  const amrex::Real St_num = ProbParm::Stmod/(8.*M_PI);
  amrex::Real refU = ProbParm::v0;
  ProbParm::partRho = ProbParm::rho0*ProbParm::rhoRatio;
  amrex::Real refRho = ProbParm::rho0;
  // Time scale for Eulerian phase
  amrex::Real tau_g = refL/refU;
  amrex::Real tau_d = St_num*tau_g;
  amrex::Real dia = std::sqrt(18.*mu*tau_d/ProbParm::partRho);
  amrex::Real Re_d = dia*refU*ProbParm::rho0/mu;
  amrex::Real error = 1000.;
  const amrex::Real tol = 1.E-6;
  int k = 0;
  const int maxIter = 500;
  // Iteratively solve for the diameter of the spray droplets
  while (error > tol) {
    amrex::Real oldDia = dia;
    Re_d = dia*refU*ProbParm::rho0/mu;
    amrex::Real C_D = 24./Re_d;
    if (Re_d > 1.) C_D *= (1. + std::pow(Re_d, 2./3.)/6.);
    dia = 0.75*refRho*C_D*refU*tau_d/ProbParm::partRho;
    error = std::abs(oldDia - dia)/dia;
    k++;
    if (k > maxIter) {
      amrex::Abort("Failed to converge a particle diameter");
    }
  }
  ProbParm::partDia = dia;
  // Output IC
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "number of particles: " << ProbParm::partNum[0];
    for (int dir = 1; dir != AMREX_SPACEDIM; ++dir)
      amrex::Print(ofs) << ", " << ProbParm::partNum[dir];
    amrex::Print(ofs) << std::endl;
    amrex::Print(ofs) << "particle mass: " << 6./M_PI*ProbParm::partRho*std::pow(ProbParm::partDia, 3) << std::endl;
    amrex::Print(ofs) << "rho_d/rho_f: " << ProbParm::rhoRatio << std::endl;
    amrex::Print(ofs) << "rho0: " << ProbParm::rho0 << std::endl;
    amrex::Print(ofs) << "cs: " << cs << std::endl;
    amrex::Print(ofs) << "U: " << ProbParm::v0 << std::endl;
    amrex::Print(ofs) << "mu: " << transport_params::const_viscosity << std::endl;
    amrex::Print(ofs) << "Re: " << ProbParm::reynolds << std::endl;
    amrex::Print(ofs) << "Stokes number: " << ProbParm::Stmod << "*Stc " << std::endl;
    amrex::Print(ofs) << "particle diameter: " << ProbParm::partDia << std::endl;
    amrex::Print(ofs) << "particle density: " << ProbParm::partRho << std::endl;
    amrex::Print(ofs) << "tau_d: " << tau_d << std::endl;
    amrex::Print(ofs) << "Re_d: " << Re_d << std::endl;
    amrex::Print(ofs) << "final time (72 tau_g): " << 72.*tau_g << std::endl;
    ofs.close();
  }
}
}

