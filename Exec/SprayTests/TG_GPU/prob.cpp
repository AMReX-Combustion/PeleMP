#include "prob.H"

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
  pp.query("reynolds", PeleC::h_prob_parm_device->reynolds);
  pp.query("mach", PeleC::h_prob_parm_device->mach);
  pp.query("convecting", PeleC::h_prob_parm_device->convecting);
  pp.query("init_redist", PeleC::prob_parm_host->numRedist);
  pp.query("ref_T", PeleC::h_prob_parm_device->T0);
  pp.query("init_O2", PeleC::h_prob_parm_device->Y_O2);
  pp.query("init_N2", PeleC::h_prob_parm_device->Y_N2);
  pp.query("num_particles", PeleC::prob_parm_host->partNum);
  pp.get("part_dia", PeleC::prob_parm_host->partDia);
  pp.get("part_temp", PeleC::prob_parm_host->partTemp);

  // Define the length scale
  PeleC::h_prob_parm_device->L = probhi[0] - problo[0];

  // Initial density, velocity, and material properties
  amrex::Real cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = PeleC::h_prob_parm_device->Y_O2;
  massfrac[N2_ID] = PeleC::h_prob_parm_device->Y_N2;
  bool get_xi = false;
  bool get_Ddiag = false;
  bool get_lambda = false;
  bool get_mu = true;
  amrex::Real wbar, gamma0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2G(PeleC::h_prob_parm_device->T0, massfrac, gamma0);
  eos.Y2WBAR(massfrac, wbar);
  // Compute the speed of sound
  cs = std::sqrt(
    PeleC::h_prob_parm_device->T0 * gamma0 * pele::physics::Constants::RU / wbar);
  amrex::Real refL = PeleC::h_prob_parm_device->L;
  PeleC::h_prob_parm_device->v0 = PeleC::h_prob_parm_device->mach * cs;
  pele::physics::transport::TransParm const* ltransparm = pele::physics::transport::trans_parm_g;
  amrex::Real rho, mu, xi, lambda;
  amrex::Real Ddiag[NUM_SPECIES]; // Should be unused
  auto trans = pele::physics::PhysicsType::transport();
  trans.transport(
    get_xi, get_mu, get_lambda, get_Ddiag, PeleC::h_prob_parm_device->T0,
    rho, // Should be unused
    massfrac, Ddiag, mu, xi, lambda, pele::physics::transport::trans_parm_g);
  // Compute the density from the Reynolds number
  PeleC::h_prob_parm_device->rho0 = PeleC::h_prob_parm_device->reynolds * mu /
                                  (refL * PeleC::h_prob_parm_device->v0);
  eos.RTY2P(
    PeleC::h_prob_parm_device->T0, PeleC::h_prob_parm_device->rho0, massfrac,
    PeleC::h_prob_parm_device->p0);

  // Output IC
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "number of particles: "
                      << PeleC::prob_parm_host->partNum[0];
    for (int dir = 1; dir != AMREX_SPACEDIM; ++dir)
      amrex::Print(ofs) << ", " << PeleC::prob_parm_host->partNum[dir];
    amrex::Print(ofs) << std::endl;
    amrex::Print(ofs) << "rho0: " << PeleC::h_prob_parm_device->rho0 << std::endl;
    amrex::Print(ofs) << "p0: " << PeleC::h_prob_parm_device->p0 << std::endl;
    amrex::Print(ofs) << "cs: " << cs << std::endl;
    amrex::Print(ofs) << "U: " << PeleC::h_prob_parm_device->v0 << std::endl;
    amrex::Print(ofs) << "Re: " << PeleC::h_prob_parm_device->reynolds
                      << std::endl;
    amrex::Print(ofs) << "particle diameter: " << PeleC::prob_parm_host->partDia
                      << std::endl;
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
