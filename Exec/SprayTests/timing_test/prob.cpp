#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real mach = 0.1;
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 1000.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_O2 = 0.233;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_N2 = 0.767;
AMREX_GPU_DEVICE_MANAGED amrex::IntVect partNum = amrex::IntVect(AMREX_D_DECL(100, 100, 100));
AMREX_GPU_DEVICE_MANAGED amrex::Real partTemp = 300.;
AMREX_GPU_DEVICE_MANAGED amrex::Real partRho = 1.;
AMREX_GPU_DEVICE_MANAGED amrex::Real partDia = 1.E-3;
AMREX_GPU_DEVICE_MANAGED amrex::RealVect partVel = amrex::RealVect(AMREX_D_DECL(0., 0., 0.));
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
  pp.query("mach", ProbParm::mach);
  pp.query("ref_p", ProbParm::p0);
  pp.query("ref_T", ProbParm::T0);
  pp.query("init_O2", ProbParm::Y_O2);
  pp.query("init_N2", ProbParm::Y_N2);
  pp.query("num_particles", ProbParm::partNum);
  std::array<amrex::Real, AMREX_SPACEDIM> pvel;
  pp.query<amrex::Real>("part_vel", pvel);
  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) ProbParm::partVel[dir] = pvel[dir];
  pp.get("part_rho", ProbParm::partRho);
  pp.get("part_dia", ProbParm::partDia);
  pp.get("part_temp", ProbParm::partTemp);

  // Initial density, velocity, and material properties
  amrex::Real cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = ProbParm::Y_O2;
  massfrac[N2_ID] = ProbParm::Y_N2;
  amrex::Real wbar, gamma0;
  EOS::TY2G(ProbParm::T0, massfrac, gamma0);
  EOS::Y2WBAR(massfrac, wbar);
  // Compute the speed of sound
  cs = std::sqrt(ProbParm::T0 * gamma0 * EOS::RU / wbar);
  ProbParm::v0 = ProbParm::mach * cs;

  // Output IC
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "number of particles: " << ProbParm::partNum[0];
    for (int dir = 1; dir != AMREX_SPACEDIM; ++dir)
      amrex::Print(ofs) << ", " << ProbParm::partNum[dir];
    amrex::Print(ofs) << std::endl;
    amrex::Print(ofs) << "particle velocity: " << ProbParm::partVel[0];
    for (int dir = 1; dir != AMREX_SPACEDIM; ++dir)
      amrex::Print(ofs) << ", " << ProbParm::partVel[dir];
    amrex::Print(ofs) << std::endl;
    amrex::Print(ofs) << "p0: " << ProbParm::p0 << std::endl;
    amrex::Print(ofs) << "cs: " << cs << std::endl;
    amrex::Print(ofs) << "U: " << ProbParm::v0 << std::endl;
    amrex::Print(ofs) << "particle diameter: " << ProbParm::partDia << std::endl;
    amrex::Print(ofs) << "particle density: " << ProbParm::partRho << std::endl;
    ofs.close();
  }
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
