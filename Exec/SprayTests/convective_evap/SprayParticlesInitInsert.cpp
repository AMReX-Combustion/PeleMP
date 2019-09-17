
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <Transport_F.H>
#include <drag_F.H>

using namespace amrex;


void
SprayParticleContainer::insertParticles (Real time, int nstep, int lev)
{
  const Geometry& geom = this->m_gdb->Geom(lev);

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {

    auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                                       mfi.LocalTileIndex())];

    const Box& bx = mfi.tilebox();
    const RealBox& temp = RealBox(bx,geom.CellSize(),geom.ProbLo());
    const Real* xlo = temp.lo();
    const Real* xhi = temp.hi();

// Take random draw from inside the jet of x width 0.02 and length
// geometry.prob_lo     =   -0.12 -0.06 0.0
// geometry.prob_hi     =  0.12 0.06 0.24

    for (int iter = 0; iter < 10; iter++)
    {
 // Real x = 0.04*(rand()%100)/99.-0.02;
    Real x = 0.10*(rand()%100)/99.-0.05;
    Real y = 0.12*(rand()%100)/99.-0.06;
    Real z = 0.24*(rand()%100)/99.;

// Check if injection point belongs to box
   if(x > xlo[0] && x < xhi[0] && y > xlo[1] && y < xhi[1] && z > xlo[2] && z < xhi[2]) {

     std::cout << " new particle at " <<x<< " " <<y<< " " <<z<< '\n' << std::endl;

     ParticleType p;
     p.id()   = ParticleType::NextID();
     p.cpu()  = ParallelDescriptor::MyProc();

     p.pos(0) = x;
     p.pos(1) = y;
     p.pos(2) = z;

     // fill in NSR_SPR items
     p.rdata(0) = 0.;
     p.rdata(1) = 0.;
     if (AMREX_SPACEDIM>2) {
       p.rdata(2) = 0.;
     }
//   p.rdata(AMREX_SPACEDIM) = 293.; // temperature
     p.rdata(AMREX_SPACEDIM) = 400.; // temperature
     p.rdata(AMREX_SPACEDIM+1) = 0.0004; // diameter
//   p.rdata(AMREX_SPACEDIM+1) = 0.0010; // diameter
     p.rdata(AMREX_SPACEDIM+2) = 0.68141; // fuel density
    
     particles.push_back(p);
   }
  } // for iter
 }
  Redistribute();
}


static bool my_first = true; // Hack to avoid adding new particle on every call
static Real t_next = 0.;

void
SprayParticleContainer::injectParticles (Real time, int nstep, int lev)
{
// Specialized for injection trough a 2D orifice - many values are hardcoded
  if (my_first) {
    srand(15);
    my_first = false;
    t_next = time;
  }

  std::cout << "TIME " << time << " at parent step "<< nstep <<"; next injection time at " << t_next << '\n' << std::endl;
  if(time < t_next) return;
  t_next+= 1e-5;

  const Real r0 = 1.5;
  const Real y = -30;
  const Real y0 = r0*3.7320; // virtual origin for 15 deg angle

  const Geometry& geom = this->m_gdb->Geom(lev);

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {

    auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                                       mfi.LocalTileIndex())];

    const Box& bx = mfi.tilebox();
    const RealBox& temp = RealBox(bx,geom.CellSize(),geom.ProbLo());
    const Real* xlo = temp.lo();
    const Real* xhi = temp.hi();

// Check if injection rake belongs to box
    if(xlo[0] <= r0 && xhi[0] >= -r0 && xlo[1] <= -30 && xhi[1] >= -30) {

      Real x = (rand()%10)/3.-r0;
      if (x > xlo[0] && x < xhi[0]) {

        ParticleType p;
        p.id()   = ParticleType::NextID();
        p.cpu()  = ParallelDescriptor::MyProc();

        p.pos(0) = x;
        p.pos(1) = y;
        if (AMREX_SPACEDIM>2) {
          p.pos(2) = 0.;
        }

        Real angle = atan(x/y0); // seen from virtual origin (0,-y0)
        // fill in NSR_SPR items
        p.rdata(0) = 1000.*sin(angle); // ux
        p.rdata(1) = 1000.*cos(angle); // uy
        if (AMREX_SPACEDIM>2) {
          p.rdata(2) = 0.0;
        }
        p.rdata(AMREX_SPACEDIM) = 293.; // temperature
        p.rdata(AMREX_SPACEDIM+1) = 0.0200; // diameter
        p.rdata(AMREX_SPACEDIM+2) = 0.68141; // fuel density
    
        particles.push_back(p);
      }
    }
  }

  Redistribute();
} 

void
SprayParticleContainer::InitParticlesUniform(const int& lev, const int& num_ppc)
{

  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {   
      const Box& tile_box  = mfi.tilebox();

      auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                                       mfi.LocalTileIndex())];
    
      for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
	for (int i_part=0; i_part<num_ppc;i_part++) {

	  Real r = (rand()%100)/99.;
	  Real x = plo[0] + (iv[0] + r)*dx[0];

	  r = (rand()%100)/99.;
	  Real y = plo[1] + (iv[1] + r)*dx[1];

	  ParticleType p;
	  p.id()  = ParticleType::NextID();
	  p.cpu() = ParallelDescriptor::MyProc();
	  p.pos(0) = x;
	  p.pos(1) = y;
	  if (AMREX_SPACEDIM>2) {
	    r = (rand()%100)/99.;
	    Real z = plo[2] + (iv[2] + r)*dx[2];
	    p.pos(2) = z;
	  }

	  //Now assign the real data carried by each particle
	  p.rdata(0) = 0.;//u-vel
	  p.rdata(1) = 0.;//v-vel
	  if (AMREX_SPACEDIM>2) {
	    p.rdata(2) = 0.;//w-vel
	  }    
	  p.rdata(AMREX_SPACEDIM) = 293.; // temperature
	  p.rdata(AMREX_SPACEDIM+1) = 0.0200; // diameter
	  p.rdata(AMREX_SPACEDIM+2) = 0.68141; // fuel density

	  particles.push_back(p);

	}   
      }

    }
    
}
