
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

using namespace amrex;

int
interpolateInjectTime(const Real& time)
{
  const int nvals = ProbParm::inject_N;
  int i = 0;
  while (i < nvals) {
    if (ProbParm::d_inject_time[i] > time) return i;
    ++i;
  }
  return -1;
}

bool
SprayParticleContainer::insertParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level)
{
  return false;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level)
{
  if (lev != 0) return false;
  if (time < ProbParm::jet_start_time
      || time > ProbParm::jet_end_time) return false;
  const Geometry& geom = this->m_gdb->Geom(lev);
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const auto dx = geom.CellSize();
  RealVect dom_len(AMREX_D_DECL(geom.ProbLength(0),
                                geom.ProbLength(1),
                                geom.ProbLength(2)));
  Real mass_flow_rate = ProbParm::mass_flow_rate;
  Real jet_vel = ProbParm::jet_vel;
  Real jet_dia = ProbParm::jet_dia;
  Real jr2 = jet_dia*jet_dia/4.; // Jet radius squared
#if AMREX_SPACEDIM == 3
  Real jet_area = M_PI*jr2;
#else
  Real jet_area = jet_dia;
#endif
  const Real Cd = 0.89;
  Real part_temp = ProbParm::part_temp;
  Real part_rho = ProbParm::part_rho;
  int ctime = 0;
  if (ProbParm::inject_N > 0) {
    const int time_indx = interpolateInjectTime(time);
    const Real time1 = ProbParm::d_inject_time[time_indx];
    const Real time2 = ProbParm::d_inject_time[time_indx - 1];
    const Real mf1 = ProbParm::d_inject_mass[time_indx];
    const Real mf2 = ProbParm::d_inject_mass[time_indx - 1];
    const Real invt = 1./(time2 - time1);
    mass_flow_rate = (mf1*(time2 - time) + mf1*(time - time1))*invt;
    const Real jv1 = ProbParm::d_inject_vel[time_indx];
    const Real jv2 = ProbParm::d_inject_vel[time_indx - 1];
    jet_vel = (jv1*(time2 - time) + jv1*(time - time1))*invt;
    ctime = time_indx;
  }
  if (jet_vel*dt/dx[0] > 0.5) {
    Real max_vel = dx[0]*0.5/dt;
    std::string warn_msg = "Injection velocity of " + std::to_string(jet_vel) 
      + " is reduced to maximum " + std::to_string(max_vel);
    jet_vel = max_vel;
    amrex::Warning(warn_msg);
  }
  Real part_dia = ProbParm::part_mean_dia;
  Real part_stdev = ProbParm::part_stdev_dia;
  Real stdsq = part_stdev*part_stdev;
  Real meansq = part_dia*part_dia;
  Real log_mean = 2.*std::log(part_dia) - 0.5*std::log(stdsq + meansq);
  Real log_stdev = std::sqrt(-2.*std::log(part_dia)
                             + std::log(stdsq + meansq));
  Real part_y = plo[1];
  RealVect jet_center(AMREX_D_DECL(dom_len[0]/2., plo[1],
                                   dom_len[2]/2.));
  Real Pi_six = M_PI/6.;
  Real spray_angle = ProbParm::spray_angle;
  Real lo_angle = -0.5*spray_angle;
  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const RealBox& temp = RealBox(bx,geom.CellSize(),geom.ProbLo());
    const Real* xlo = temp.lo();
    const Real* xhi = temp.hi();
    RealVect box_len(AMREX_D_DECL(temp.length(0), 0., temp.length(2)));
    Gpu::HostVector<ParticleType> host_particles;
    if (xlo[1] == plo[1]) {
      // Box locations relative to jet center
      const RealVect xloJ(AMREX_D_DECL(xlo[0] - jet_center[0], plo[1],
                                       xlo[1] - jet_center[2]));
      const RealVect xhiJ(AMREX_D_DECL(xhi[0] - jet_center[0], plo[1],
                                       xhi[2] - jet_center[2]));
      Real cur_jet_area = 0.;
      Real testdx = dx[0]/100.;
      Real testdx2 = testdx*testdx;
      Real curx = xloJ[0];
      // Loop over each cell and check how much overlap there is with the jet
#if AMREX_SPACEDIM == 3
      while (curx < xhiJ[0]) {
        Real curz = xloJ[2];
        while (curz < xhiJ[2]) {
          Real r2 = curx*curx + curz*curz;
          if (r2 <= jr2) cur_jet_area += testdx2;
          curz += testdx;
        }
        curx += testdx;
      }
#else
      while (curx < xhiJ[0]) {
        Real r2 = curx*curx;
        if (r2 <= jr2) cur_jet_area += testdx;
        curx += testdx;
      }
#endif
      Real jet_perc = cur_jet_area/jet_area;
      Real perc_mass = jet_perc*mass_flow_rate*dt;
      Real total_mass = 0.;
      int np = 0;
      while (total_mass < perc_mass) {
        RealVect part_loc(AMREX_D_DECL(xlo[0] + amrex::Random()*box_len[0],
                                       plo[1],
                                       xlo[2] + amrex::Random()*box_len[2]));
        Real r2 = AMREX_D_TERM(std::pow(part_loc[0] - jet_center[0], 2),,
                               + std::pow(part_loc[2] - jet_center[2], 2));
        if (r2 <= jr2) {
          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          Real theta = lo_angle + spray_angle*amrex::Random();
#if AMREX_SPACEDIM == 3
          Real phi = lo_angle + spray_angle*amrex::Random();
#else
          Real phi = 0.;
#endif
          AMREX_D_TERM(p.rdata(PeleC::pstateVel) = jet_vel*std::sin(theta)*std::cos(phi);,
                       p.rdata(PeleC::pstateVel+1) = jet_vel*std::cos(theta);,
                       p.rdata(PeleC::pstateVel+2) = jet_vel*std::sin(theta)*std::sin(phi););
          Real cur_dia = amrex::RandomNormal(log_mean, log_stdev);
          // Use a log normal distribution
          cur_dia = std::exp(cur_dia);
          for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
            p.pos(dir) = part_loc[dir];
          p.rdata(PeleC::pstateT) = part_temp;
          p.rdata(PeleC::pstateDia) = cur_dia;
          p.rdata(PeleC::pstateRho) = part_rho;
          for (int sp = 0; sp != SPRAY_FUEL_NUM; ++sp)
            p.rdata(PeleC::pstateY + sp) = 0.;
          p.rdata(PeleC::pstateY) = 1.;
          host_particles.push_back(p);
          Real pmass = Pi_six*part_rho*std::pow(cur_dia, 3);
          total_mass += pmass;
        }
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
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(AmrLevel* pelec, const int& lev, const int& num_ppc)
{
  // Start without any particles
  return;
}
