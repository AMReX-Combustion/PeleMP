
#include "SprayParticles.H"
#include <AMReX_Particles.H>
#include "pelelm_prob.H"

int
interpolateInjectTime(
  const amrex::Real& time, const int nvals, const amrex::Real* inject_time)
{
  int i = 0;
  amrex::Real ctime = inject_time[i];
  while (ctime < time) {
    ctime = inject_time[++i];
  }
  return i;
}

amrex::Real
jetOverlapArea(
  const amrex::Real testdx,
  const amrex::RealVect xloJ,
  const amrex::RealVect xhiJ,
  const amrex::RealVect jet_cent,
  const amrex::Real jr2,
#if AMREX_SPACEDIM == 3
  amrex::Real& loy,
  amrex::Real& hiy,
#endif
  amrex::Real& lox,
  amrex::Real& hix)
{
  // This is how much we reduce x to test for the area of the jet
  // This helps if the jet size is smaller than cell size
  amrex::Real cur_jet_area = 0.;
  amrex::Real testdx2 = testdx * testdx;
  amrex::Real curx = xloJ[0];
  hix = xloJ[0];
  lox = xhiJ[0];
  // Loop over each cell and check how much overlap there is with the jet
#if AMREX_SPACEDIM == 3
  hiy = xloJ[1];
  loy = xhiJ[1];
  while (curx < xhiJ[0]) {
    amrex::Real cury = xloJ[1];
    while (cury < xhiJ[1]) {
      amrex::Real r2 = curx * curx + cury * cury;
      if (r2 <= jr2) {
        cur_jet_area += testdx2;
        lox = amrex::min(curx, lox);
        hix = amrex::max(curx, hix);
        loy = amrex::min(cury, loy);
        hiy = amrex::max(cury, hiy);
      }
      cury += testdx;
    }
    curx += testdx;
  }
#else
  while (curx < xhiJ[0]) {
    amrex::Real r2 = curx * curx;
    if (r2 <= jr2) {
      cur_jet_area += testdx;
      lox = amrex::min(curx, lox);
      hix = amrex::max(curx, hix);
    }
    curx += testdx;
  }
#endif
  return cur_jet_area;
}

amrex::RealVect
jet_vel_comp(
  amrex::Real theta_1,
  amrex::Real phi_1,
  amrex::Real theta_2,
  amrex::Real phi_2)
{
  amrex::Real st1 = std::sin(theta_1);
  amrex::Real ct1 = std::cos(theta_1);
  amrex::Real st2 = std::sin(theta_2);
  amrex::Real ct2 = std::cos(theta_2);
  amrex::Real sp1 = std::sin(phi_1);
  amrex::Real cp1 = std::cos(phi_1);
  amrex::Real sp2 = std::sin(phi_2);
  amrex::Real cp2 = std::cos(phi_2);
  amrex::Real v1 = st1 * ct2 + st2 * cp2 * ct1;
  amrex::RealVect vel;
  AMREX_D_TERM(vel[0] = sp1 * v1 + sp2 * st2 * cp1;
               , vel[1] = cp1 * v1 - sp2 * st2 * sp1;
               , vel[2] = -st1 * st2 * cp2 + ct1 * ct2;);
  return vel;
}

bool
SprayParticleContainer::injectParticles(
  amrex::Real time,
  amrex::Real dt,
  int nstep,
  int lev,
  int finest_level,
  ProbParm const& prob_parm)
{
  amrex::ignore_unused(nstep, finest_level);
  int lowD = AMREX_SPACEDIM - 1;
  if (lev != 0) {
    return false;
  }
  if (time < prob_parm.jet_start_time || time > prob_parm.jet_end_time) {
    return false;
  }
  const int pstateVel = m_sprayIndx.pstateVel;
  const int pstateT = m_sprayIndx.pstateT;
  const int pstateDia = m_sprayIndx.pstateDia;
  const int pstateY = m_sprayIndx.pstateY;
  const SprayData* fdat = m_sprayData;
  amrex::Real rho_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    rho_part += prob_parm.Y_jet[spf] / fdat->rho[spf];
  }
  rho_part = 1. / rho_part;
  // Number of particles per parcel
  const amrex::Real num_ppp = fdat->num_ppp;
  const amrex::Geometry& geom = this->m_gdb->Geom(lev);
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const auto dx = geom.CellSize();
  amrex::RealVect dom_len(
    AMREX_D_DECL(geom.ProbLength(0), geom.ProbLength(1), geom.ProbLength(2)));
  amrex::Real mass_flow_rate = prob_parm.mass_flow_rate;
  amrex::Real jet_vel = prob_parm.jet_vel;
  amrex::Real jet_dia = prob_parm.jet_dia;
  amrex::Real newdxmod = dx[0] / jet_dia * 10.;
  amrex::Real dx_mod = amrex::max(50., newdxmod);
  amrex::Real jr2 = jet_dia * jet_dia / 4.; // Jet radius squared
#if AMREX_SPACEDIM == 3
  amrex::Real jet_area = M_PI * jr2;
#else
  amrex::Real jet_area = jet_dia;
#endif
  amrex::Real part_temp = prob_parm.part_temp;
  // This absolutely must be included with any injection or insertion
  // function or significant issues will arise
  if (jet_vel * dt / dx[0] > m_partCFL) {
    amrex::Real max_vel = dx[0] * m_partCFL / dt;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::string warn_msg =
        "Injection velocity of " + std::to_string(jet_vel) +
        " is reduced to maximum " + std::to_string(max_vel);
      amrex::Warning(warn_msg);
    }
    m_injectVel = jet_vel;
    jet_vel = max_vel;
  }
  amrex::Real part_dia = prob_parm.part_mean_dia;
  amrex::Real part_stdev = prob_parm.part_stdev_dia;
  amrex::Real stdsq = part_stdev * part_stdev;
  amrex::Real meansq = part_dia * part_dia;
  amrex::Real log_mean =
    2. * std::log(part_dia) - 0.5 * std::log(stdsq + meansq);
  amrex::Real log_stdev = std::sqrt(
    amrex::max(-2. * std::log(part_dia) + std::log(stdsq + meansq), 0.));
  amrex::Real Pi_six = M_PI / 6.;
  amrex::Real spray_angle = prob_parm.spread_angle;
  // We have 2 spherical coordinate systems
  // 1: Jet spherical coordinates for the center of the jet relative to
  // Cartesian 2: Spread spherical coordinates for spread of particle relative
  // center of jet Within these coordinates, theta is the polar angle and phi is
  // the azimuthal
  amrex::Real theta_1 = prob_parm.jet_angle;
  amrex::Real theta_2_lo = -0.5 * spray_angle;
  // Range for the azimuthal coordinate
  amrex::Real phi_range = prob_parm.jet_azi_hi - prob_parm.jet_azi_lo;
  for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.tilebox();
    const amrex::RealBox& temp =
      amrex::RealBox(bx, geom.CellSize(), geom.ProbLo());
    const amrex::Real* xloB = temp.lo();
    const amrex::Real* xhiB = temp.hi();
    amrex::RealVect cur_jet_cent;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      cur_jet_cent[dir] = prob_parm.jet_cent[dir];
    }
    const amrex::Real jet_norm = cur_jet_cent[lowD];
    if (xloB[lowD] <= jet_norm && xhiB[lowD] > jet_norm) {
      amrex::Gpu::HostVector<ParticleType> host_particles;
      amrex::RealVect xlo;
      amrex::RealVect xhi;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        xlo[dir] = amrex::max(xloB[dir], cur_jet_cent[dir] - 3. * jet_dia);
        xhi[dir] = amrex::min(xhiB[dir], cur_jet_cent[dir] + 3. * jet_dia);
      }
      amrex::RealVect xloJ;
      amrex::RealVect xhiJ;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        xloJ[dir] = xlo[dir] - cur_jet_cent[dir];
        xhiJ[dir] = xhi[dir] - cur_jet_cent[dir];
      }
      // Set jet to lower Y domain for 2D and lower Z domain for 3D
      xloJ[lowD] = cur_jet_cent[lowD];
      xhiJ[lowD] = cur_jet_cent[lowD];
      amrex::Real lox, hix;
#if AMREX_SPACEDIM == 3
      amrex::Real loy, hiy;
      amrex::Real cur_jet_area = jetOverlapArea(
        dx[0] / dx_mod, xloJ, xhiJ, cur_jet_cent, jr2, loy, hiy, lox, hix);
      amrex::Real ylen = hiy - loy;
      loy += cur_jet_cent[1];
#else
      amrex::Real cur_jet_area =
        jetOverlapArea(dx[0] / dx_mod, xloJ, xhiJ, cur_jet_cent, jr2, lox, hix);
#endif
      amrex::Real xlen = hix - lox;
      lox += cur_jet_cent[0];
      if (cur_jet_area > 0.) {
        amrex::Real jet_perc = cur_jet_area / jet_area;
        amrex::Real perc_mass = jet_perc * mass_flow_rate * dt;
        amrex::Real total_mass = 0.;
        while (total_mass < perc_mass) {
          amrex::RealVect part_loc;
          part_loc[0] = lox + amrex::Random() * xlen;
          AMREX_D_PICK(, part_loc[1] = jet_norm;
                       , part_loc[1] = loy + amrex::Random() * ylen;
                       part_loc[2] = jet_norm;)
          amrex::Real r2 = AMREX_D_TERM(
            std::pow(part_loc[0] - cur_jet_cent[0], 2), ,
            +std::pow(part_loc[1] - cur_jet_cent[1], 2));
          if (r2 <= jr2) {
            ParticleType p;
            p.id() = ParticleType::NextID();
            p.cpu() = amrex::ParallelDescriptor::MyProc();
            // Determine the polar angle for the spread where 0 is perpendicular
            // to the theta_2 axis
            amrex::Real theta_2 = theta_2_lo + amrex::Random() * spray_angle;
#if AMREX_SPACEDIM == 3
            // Determine the azimuthal angle for the jet center
            amrex::Real phi_1 =
              prob_parm.jet_azi_lo + amrex::Random() * phi_range;
            // Determine the azimuthal angle for the spread
            amrex::Real phi_2 = 2. * M_PI * amrex::Random();
#else
            amrex::Real phi_1 = 0.;
            amrex::Real phi_2 = 0.;
#endif
            amrex::RealVect part_vel =
              jet_vel_comp(theta_1, phi_1, theta_2, phi_2);
            AMREX_D_TERM(p.rdata(pstateVel) = jet_vel * part_vel[0];
                         , p.rdata(pstateVel + 1) = jet_vel * part_vel[1];
                         , p.rdata(pstateVel + 2) = jet_vel * part_vel[2];);
            amrex::Real cur_dia = amrex::RandomNormal(log_mean, log_stdev);
            // Use a log normal distribution
            cur_dia = std::exp(cur_dia);
            // Add particles as if they have advanced some random portion of
            // dt
            amrex::Real pmov = amrex::Random();
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
              p.pos(dir) = part_loc[dir] + pmov * dt * part_vel[dir];
            }
            p.rdata(pstateT) = part_temp;
            p.rdata(pstateDia) = cur_dia;
            for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
              p.rdata(pstateY + sp) = prob_parm.Y_jet[sp];
            }
            host_particles.push_back(p);
            amrex::Real pmass = Pi_six * rho_part * std::pow(cur_dia, 3);
            total_mass += num_ppp * pmass;
          }
        }
      }
      if (host_particles.size() > 0) {
        auto& particle_tile =
          GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, host_particles.begin(),
          host_particles.end(),
          particle_tile.GetArrayOfStructs().begin() + old_size);
      }
    }
  }
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(ProbParm const& prob_parm)
{
  // Start without any particles
  m_injectVel = prob_parm.jet_vel;
  return;
}
