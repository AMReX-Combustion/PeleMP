
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

using namespace amrex;

bool
SprayParticleContainer::insertParticles (Real time, int nstep, int lev)
{
  return false;
}

bool
SprayParticleContainer::injectParticles (Real time, int nstep, int lev)
{
  return false;
} 

void
SprayParticleContainer::InitParticlesUniform(AmrLevel* pelec, const int& lev, const int& num_ppc)
{
  const int numGrids = pelec->numGrids();
  Real part_rho = ProbParm::partRho;
  Real part_dia = ProbParm::partDia;
  Real T_ref = ProbParm::T0;
  const IntVect num_part = ProbParm::partNum;
  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();
  const auto phi = Geom(lev).ProbHiArray();
  const Box& boxDom = Geom(lev).Domain();
  const Real len = phi[0] - plo[0];
  const RealVect dx_part(AMREX_D_DECL(len/Real(num_part[0]),
					     len/Real(num_part[1]),
					     len/Real(num_part[2])));
  // Number of particles per cell per direction
  // Makes sure there is at least 1 per cell
  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    Box tile_box = mfi.tilebox();
    tile_box &= boxDom;
    const RealBox tile_realbox(tile_box, Geom(lev).CellSize(), Geom(lev).ProbLo());
    Gpu::HostVector<ParticleType> host_particles;
    RealVect lo_end;
    IntVect ppc_pdir;
    for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
      Real box_length = tile_realbox.length(dir);
      lo_end[dir] = tile_realbox.lo(dir);
      ppc_pdir[dir] = ProbParm::partNum[dir]*box_length/len;
    }
    RealVect part_loc = RealVect::Zero;
    for (int x_part = 0; x_part != ppc_pdir[0]; ++x_part) {
      part_loc[0] = lo_end[0] + (0.5 + Real(x_part))*dx_part[0];
      for (int y_part = 0; y_part != ppc_pdir[1]; ++y_part) {
	part_loc[1]= lo_end[1] + (0.5 + Real(y_part))*dx_part[1];
#if AMREX_SPACEDIM == 3
	for (int z_part = 0; z_part != ppc_pdir[2]; ++z_part) {
	  part_loc[2] = lo_end[2] + (0.5 + Real(z_part))*dx_part[2];
#endif
	  ParticleType p;
	  p.id() = ParticleType::NextID();
	  p.cpu() = ParallelDescriptor::MyProc();
	  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
	    p.pos(dir) = part_loc[dir];
	    p.rdata(PeleC::pstateVel+dir) = 0.;
	  }
	  p.rdata(PeleC::pstateT) = T_ref; // temperature
	  p.rdata(PeleC::pstateDia) = part_dia; // diameter
	  p.rdata(PeleC::pstateRho) = part_rho; // liquid fuel density
	  for (int sp = 0; sp != SPRAY_FUEL_NUM; ++sp) {
	    p.rdata(PeleC::pstateY + sp) = 0.;
	  }
	  p.rdata(PeleC::pstateY) = 1.; // Only use the first fuel species
	  host_particles.push_back(p);
#if AMREX_SPACEDIM == 3
	}
#endif
      }
    }
    auto& particles = GetParticles(lev);
    auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
    auto old_size = particle_tile.GetArrayOfStructs().size();
    auto new_size = old_size + host_particles.size();
    particle_tile.resize(new_size);

    // Copy the AoS part of the host particles to the GPU
    Gpu::copy(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
	      particle_tile.GetArrayOfStructs().begin() + old_size);
  }
  Redistribute();
}
