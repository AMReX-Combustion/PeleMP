
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

using namespace amrex;

bool
SprayParticleContainer::insertParticles (Real time, Real dt, int nstep, int lev)
{
  return false;
}

bool
SprayParticleContainer::injectParticles (Real time, Real dt, int nstep, int lev)
{
  // Move these to prob_params
  Real start_time = 0.;
  if (time < start_time) return false;
  Real inj_dt = 6.E-7;
  int num_per_old = int((time - dt)/inj_dt);
  int num_per_new = int(time/inj_dt);
  Real next_inj_time = (num_per_old + 1)*inj_dt;
  const Real eps =
    std::numeric_limits<Real>::epsilon()*10.*time;
  if ((num_per_old == num_per_new)
      && std::abs(time - next_inj_time) <= eps)
    num_per_new += 1;
  if ((num_per_old != num_per_new)
      && std::abs((time - dt) - next_inj_time) <= eps)
    num_per_old += 1;
  if (num_per_old == num_per_new) return false; 
  const Geometry& geom = this->m_gdb->Geom(lev);
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  Real dom_lenx = (phi[0] - plo[0]);
  Real jet_dia = ProbParm::jet_dia;
  Real jet_vel = ProbParm::jet_vel;
  Real part_temp = ProbParm::part_temp;
  Real part_rho = ProbParm::part_rho;
  Real part_dia = ProbParm::part_dia;
  Real part_y = 0.;
  Real jet_lo = dom_lenx/2. - jet_dia/2.;
  Real jet_hi = dom_lenx/2. + jet_dia/2.;
  Real part_x_loc[10];
  Real part_dx = (jet_hi - jet_lo)/10.;
  for (int np = 0; np != 10; ++np) {
    part_x_loc[np] = jet_lo + part_dx*(np + 0.5);
  }
  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const RealBox& temp = RealBox(bx,geom.CellSize(),geom.ProbLo());
    const Real* xlo = temp.lo();
    const Real* xhi = temp.hi();
    Gpu::HostVector<ParticleType> host_particles;
    if (xlo[1] == plo[1]) {
      for (int np = 0; np != 10; ++np) {
        if (part_x_loc[np] >= xlo[0] && 
            part_x_loc[np] < xhi[0]) {
          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
            p.rdata(PeleC::pstateVel+dir) = 0.;
          p.rdata(PeleC::pstateVel+1) = jet_vel;
          AMREX_D_TERM(p.pos(0) = part_x_loc[np];
                       , p.pos(1) = part_y;, p.pos(2) = 0.;);
          p.rdata(PeleC::pstateT) = part_temp;
          p.rdata(PeleC::pstateDia) = part_dia;
          p.rdata(PeleC::pstateRho) = part_rho;
          for (int sp = 0; sp != SPRAY_FUEL_NUM; ++sp)
            p.rdata(PeleC::pstateY + sp) = 0.;
          p.rdata(PeleC::pstateY) = 1.;
          host_particles.push_back(p);
        }
      }
      auto& particle_tile = GetParticles(lev)[std::make_pair(mfi.index(),
                                                         mfi.LocalTileIndex())];
      auto old_size = particle_tile.GetArrayOfStructs().size();
      auto new_size = old_size + host_particles.size();
      particle_tile.resize(new_size);

      Gpu::copy(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
                particle_tile.GetArrayOfStructs().begin() + old_size);
    }
  }
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(AmrLevel* pelec, const int& lev, const int& num_ppc)
{
  // Start without any particles
  return;
}
