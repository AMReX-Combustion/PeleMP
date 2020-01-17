
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include <Transport_F.H>
//#include <Drag2d.H>
//#include <drag_F.H>

using namespace amrex;


void
SprayParticleContainer::insertParticles (Real time, int nstep, int lev)
{

}


void
SprayParticleContainer::injectParticles (Real time, int nstep, int lev)
{

} 

void
SprayParticleContainer::InitParticlesUniform(const int& lev, const int& num_ppc)
{

  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {   
      const Box& tile_box  = mfi.tilebox();

      auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
							 mfi.LocalTileIndex())];
    
      for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
	for (int i_part=0; i_part<num_ppc;i_part++) {

	  ParticleType p;
	  p.id()  = ParticleType::NextID();
	  p.cpu() = ParallelDescriptor::MyProc();
	  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
	    {
	      Real r = (rand()%100)/99.;
	      Real loc = plo[dir] + (iv[dir] + r)*dx[dir];
	      p.pos(dir) = loc;
	      p.rdata(PeleC::pstate_vel+dir) = 0.;
	    }

	  p.rdata(PeleC::pstate_T) = 293.; // temperature
	  p.rdata(PeleC::pstate_dia) = 0.0200; // diameter
	  p.rdata(PeleC::pstate_rho) = 0.68141; // fuel density

	  particles.push_back(p);
	}   
      }

    }
    
}
