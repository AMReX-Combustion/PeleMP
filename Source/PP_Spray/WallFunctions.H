
#ifndef WALLFUNCTIONS_H
#define WALLFUNCTIONS_H

#include "AhamedSplash.H"

// Check if tile is adjacent to non-periodic boundaries
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
bool
tile_at_bndry(
  const amrex::Box& in_box,
  const amrex::IntVect& bndry_lo,
  const amrex::IntVect& bndry_hi,
  const amrex::Box& domain)
{
  amrex::Box testBox(in_box);
  testBox.grow(1);
  if (domain.contains(testBox)) {
    return false;
  }
  bool at_bndry = false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    amrex::Box boxlo(in_box);
    boxlo.growLo(dir, 1);
    amrex::Box boxhi(in_box);
    boxhi.growHi(dir, 1);
    if (
      (!domain.contains(boxlo) && bndry_lo[dir] != 0) ||
      (!domain.contains(boxhi) && bndry_hi[dir] != 0)) {
      at_bndry = true;
    }
  }
  return at_bndry;
}

// This is the same as the check_bounds function but pos and phi have been
// normalized by plo and dx
AMREX_GPU_DEVICE AMREX_INLINE bool
check_bounds_norm(
  const amrex::RealVect pos,
  const amrex::RealVect phi,
  const amrex::IntVect bndry_lo,
  const amrex::IntVect bndry_hi,
  amrex::IntVect& bflags)
{
  for (int hilo = 0; hilo < 2; ++hilo) {
    int fact = (hilo == 0) ? 1 : -1;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::Real dom_loc = (hilo == 0) ? 0. : phi[dir];
      int bflag = (hilo == 0) ? bndry_lo[dir] : bndry_hi[dir];
      if (bflag != 0) { // Not periodic
        amrex::Real diff =
          static_cast<amrex::Real>(fact) * (pos[dir] - dom_loc);
        if (diff < 0.) {
          // Outside domain
          if (bflag == 1) { // Outside reflective BC
            bflags[dir] = -fact;
          } else { // Outside non-reflective BC
            return true;
          }
          // Adjacent to non-periodic boundary
        } else if (diff < 0.5) {
          bflags[dir] = -2 * fact;
        }
      }
    }
  }
  return false;
}

/**
Determines if particle must be reflected off BC or EB
@param[in] pos Position of the particle after advecting, normalized by dx
@param[in] bflags Flags if particle is outside reflective boundaries
@param[in] ijkc Grid cell index containing particle
@param[in] ijkc_prev Previous grid cell index containing particle
@param[in] vel_part Velocity of the particle
@param[in] use_EB Flag if EB is used in the current box
@param[in] flags Array of flags denoting if a cell has EB in it
@param[in] bcent Array of EB centroids for each cell
@param[in] bnorm Array of EB normal vectors for each cell
@param[in] vfrac Array of EB volume fractions for each cell
@param[in] min_eb_vfrac Minimum volume fraction to allow particles in
@param[out] par_dot Magnitude of vector from wall plane to particle, negative
means outside the domain
@param[out] normal Boundary normal vector
@return Flag if particle must be reflected
*/
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
bool
check_wall(
  const amrex::RealVect& pos,
  const amrex::IntVect& bflags,
  const amrex::IntVect& ijkc,
#ifdef AMREX_USE_EB
  const amrex::IntVect& ijkc_prev,
  const amrex::RealVect& vel_part,
  const bool use_EB,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<amrex::Real const> const& bcent,
  amrex::Array4<amrex::Real const> const& bnorm,
  amrex::Array4<amrex::Real const> const& vfrac,
  const amrex::Real min_eb_vfrac,
#endif
  amrex::Real& par_dot,
  amrex::RealVect& normal)
{
  bool wall_check = false; // Check if particle is reflected
  amrex::RealVect bcentv;  // Absolute location of point on boundary face
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    // -1 - Outside lower reflective boundary
    // 1 - Outside upper reflective boundary
    // Only accounts for 1  Cartesian direction
    if ((bflags[dir] == -1 || bflags[dir] == 1) && !wall_check) {
      normal[dir] = -static_cast<amrex::Real>(bflags[dir]);
      bcentv[dir] = static_cast<amrex::Real>(ijkc[dir]) +
                    0.5 * static_cast<amrex::Real>(1 - bflags[dir]);
      wall_check = true;
    } else {
      normal[dir] = 0.;
      bcentv[dir] = 0.;
    }
  }
#ifdef AMREX_USE_EB
  if (use_EB && !wall_check) {
    if (flags(ijkc).isSingleValued()) {
      wall_check = true;
      normal = {
        AMREX_D_DECL(-bnorm(ijkc, 0), -bnorm(ijkc, 1), -bnorm(ijkc, 2))};
      bcentv = {AMREX_D_DECL(
        ijkc[0] + 0.5 + bcent(ijkc, 0), ijkc[1] + 0.5 + bcent(ijkc, 1),
        ijkc[2] + 0.5 + bcent(ijkc, 2))};
    } else if (flags(ijkc).isCovered()) {
      wall_check = true;
      if (flags(ijkc_prev).isSingleValued()) {
        normal = {AMREX_D_DECL(
          -bnorm(ijkc_prev, 0), -bnorm(ijkc_prev, 1), -bnorm(ijkc_prev, 2))};
        bcentv = {AMREX_D_DECL(
          ijkc_prev[0] + 0.5 + bcent(ijkc_prev, 0),
          ijkc_prev[1] + 0.5 + bcent(ijkc_prev, 1),
          ijkc_prev[2] + 0.5 + bcent(ijkc_prev, 2))};
      } else {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          if (ijkc[dir] != ijkc_prev[dir]) {
            amrex::Real face_sgn = (ijkc[dir] > ijkc_prev[dir]) ? -1. : 1.;
            normal[dir] = face_sgn;
            bcentv[dir] = static_cast<amrex::Real>(ijkc[dir]) + 0.5 * face_sgn;
          } else {
            normal[dir] = 0.;
            bcentv[dir] = static_cast<amrex::Real>(ijkc[dir]);
          }
        }
        normal /= normal.vectorLength();
      }
    }
  }
#endif
  bool refl_check = false;
  if (wall_check) {
    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
    // Projection of vector pointing from EB centroid to particle onto EB normal
    par_dot = AMREX_D_TERM(
      (pos[0] - bcentv[0]) * normal[0], +(pos[1] - bcentv[1]) * normal[1],
      +(pos[2] - bcentv[2]) * normal[2]);
    if (par_dot < tolerance) {
      refl_check = true;
    }
  }
  return refl_check;
}

/*
Imposes BC or EB surfaces. Splash model is applied if it is turned on.
If the particle reflects multiple times, splash model is not applied to avoid
creating particles outside of boundary
*/
AMREX_GPU_DEVICE
AMREX_INLINE
void
impose_wall(
  bool isActive,
  int pid,
  SprayParticleContainer::ParticleType& p,
  const SprayData& fdat,
  const SprayComps& SPI,
  const amrex::RealVect& dx,
  const amrex::RealVect& plo,
  const amrex::RealVect& phi,
  const amrex::IntVect bndry_lo,
  const amrex::IntVect bndry_hi,
  amrex::IntVect& bflags,
  const amrex::Real* cBoilT,
  const amrex::Real p_fluid,
  const bool use_EB,
#ifdef AMREX_USE_EB
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<amrex::Real const> const& bcent,
  amrex::Array4<amrex::Real const> const& bnorm,
  amrex::Array4<amrex::Real const> const& vfrac,
  const amrex::Real min_eb_vfrac,
#endif
  const amrex::IntVect& ijkc_prev,
  splash_breakup* N_SB,
  const SBPtrs& rf,
  const amrex::Real rem_dt)
{
  amrex::Real par_dot = 2.;
  amrex::RealVect normal;
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SPI.pstateVel), p.rdata(SPI.pstateVel + 1),
    p.rdata(SPI.pstateVel + 2)));
  // Normalize positions
  amrex::RealVect normphi;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    normphi[dir] = (phi[dir] - plo[dir]) / dx[dir];
    p.pos(dir) = (p.pos(dir) - plo[dir]) / dx[dir];
  }
  amrex::IntVect ijkc = p.pos().floor();
  // Check if particle is behind a wall
  bool do_refl = check_wall(
    p.pos(), bflags, ijkc,
#ifdef AMREX_USE_EB
    ijkc_prev, vel_part, use_EB, flags, bcent, bnorm, vfrac, min_eb_vfrac,
#endif
    par_dot, normal);
  int bounce = 0;
  int max_bounce = 5;
  amrex::IntVect refl_prev = ijkc;
  bool do_deposit = false;
  while (do_refl) {
    amrex::Real Nw_Vp = normal.dotProduct(vel_part);
    // Reflect the particle to see if it interacts with another surface or
    // leaves the domain
    amrex::RealVect refl_pos = p.pos() - 1.8 * par_dot * normal;
    amrex::RealVect refl_vel = vel_part - 1.8 * Nw_Vp * normal;
    amrex::IntVect refl_ijkc = refl_pos.floor();
    bool left_dom =
      check_bounds_norm(refl_pos, normphi, bndry_lo, bndry_hi, bflags);
    if (left_dom) {
      p.id() = -1;
      return;
    }
    amrex::Real refl_par_dot = 2.;
    amrex::RealVect refl_normal;
    do_refl = check_wall(
      refl_pos, bflags, refl_ijkc,
#ifdef AMREX_USE_EB
      refl_prev, refl_vel, use_EB, flags, bcent, bnorm, vfrac, min_eb_vfrac,
#endif
      refl_par_dot, refl_normal);
    if (fdat.wall_T > 0.) {
      // If it reflects again, deposit particle at wall
      if (do_refl) {
        do_deposit = true;
        do_refl = false;
      } else {
        // If not bouncing and splash model is on
        if (isActive) {
          droplet_splashing(
            p, pid, dx, plo, SPI, fdat, p_fluid, vel_part, par_dot, normal,
            cBoilT, N_SB, rf, true, rem_dt);
        } else {
          // Otherwise, just reflect particle
          p.pos() = refl_pos;
          vel_part = refl_vel;
        }
      }
    } else {
      if (!use_EB) {
        refl_prev = p.pos().floor();
      }
#ifdef AMREX_USE_EB
      else if (!flags(ijkc).isCovered()) {
        refl_prev = p.pos().floor();
      }
#endif
      p.pos() = refl_pos;
      vel_part = refl_vel;
      normal = refl_normal;
      par_dot = refl_par_dot;
      if (do_refl) {
        bounce++;
        if (bounce >= max_bounce) {
          do_deposit = true;
          do_refl = false;
        }
      }
    }
  } // while (do_refl)
  // If deposited, place at wall
  // TODO: Implement wall film model for this scenario
  if (do_deposit) {
    amrex::Real Nw_Vp = normal.dotProduct(vel_part);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.pos(dir) -= 1.1 * par_dot * vel_part[dir] / Nw_Vp;
    }
    vel_part = amrex::RealVect::TheZeroVector();
  }
  // Undo normalization and set new velocity
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    p.pos(dir) = p.pos(dir) * dx[dir] + plo[dir];
    p.rdata(SPI.pstateVel + dir) = vel_part[dir];
  }
}

#if 0
/****************************************************************
 Functions computing wall film source terms
 WIP
 ***************************************************************/

AMREX_GPU_DEVICE AMREX_INLINE void
calculateWallFilmSource(
  const amrex::Real flow_dt,
  GasPhaseVals& gpv,
  const SprayComps& SPI,
  const SprayData& fdat,
  SprayParticleContainer::ParticleType& p,
  const amrex::Real /*face_area*/,
  const amrex::Real diff_cent,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* trans_parm)
{
  auto eos = pele::physics::PhysicsType::eos();
  SprayUnits SPU;
  const amrex::Real C_eps = 1.E-15;
  const amrex::Real ht_tol = 2.E-6;
  const amrex::Real part_dt = 0.5 * flow_dt;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_film;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_vapor;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> L_fuel;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mi_dot;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cBoilT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Psat;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_skin;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> h_skin;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Ddiag;
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    Y_skin[sp] = 0.;
  }
  // Get particle variables
  // Note: the diameter is set assuming a sphere
  // of the same volume as the film
  // and temperature is set to a negative value
  amrex::Real T_film = -p.rdata(SPI.pstateT);
  amrex::Real vol = p.rdata(SPI.pstateVol);
  amrex::Real ht_film = p.rdata(SPI.pstateHt);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_film[spf] = p.rdata(SPI.pstateY + spf);
  }
  calcBoilT(fdat, gpv, cBoilT.data());
  // TODO: Add model for when T > T_boil
  amrex::Real T_i = gpv.T_fluid;
  amrex::Real area_film = vol / ht_film;
  amrex::Real dia_film = std::sqrt(4. * area_film / M_PI);
  eos.T2Hi(T_film, h_skin.data());
  for (int n = 0; n < NUM_SPECIES; ++n) {
    h_skin[n] *= SPU.eng_conv;
  }
  calcVaporY(
    fdat, gpv, T_film, C_eps, Y_film.data(), h_skin.data(), cBoilT.data(),
    Y_vapor.data(), Psat.data(), L_fuel.data());
  amrex::Real T_vapor = 0.5 * (T_film + T_i);
  amrex::Real sumYSkin = 0.;
  amrex::Real sumYFuel = 0.;
  amrex::Real rho_film = 0.;
  amrex::Real lambda_film = 0.;
  amrex::Real cp_film = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    rho_film += Y_film[spf] / fdat.rho[spf];
    lambda_film += Y_film[spf] * fdat.lambda[spf];
    cp_film += Y_film[spf] * fdat.cp[spf];
    Y_skin[fspec] = 0.5 * (Y_vapor[spf] + gpv.Y_fluid[fspec]);
    sumYSkin += Y_skin[fspec];
    sumYFuel += gpv.Y_fluid[fspec];
  }
  rho_film = 1. / rho_film;
  const amrex::Real restYSkin = 1. - sumYSkin;
  const amrex::Real restYfluid = 1. - sumYFuel;
  amrex::Real renorm = restYSkin / restYfluid;
  amrex::Real mw_vap = 0.;
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    if (Y_skin[sp] == 0.) {
      Y_skin[sp] = gpv.Y_fluid[sp] * renorm;
    }
    mw_vap += Y_skin[sp] * gpv.invmw[sp];
  }
  mw_vap = 1. / mw_vap;
  amrex::Real pmass = vol * rho_film;
  amrex::Real lambda_skin = 0.;
  amrex::Real mu_skin = 0.;
  amrex::Real xi_skin = 0.;
  auto trans = pele::physics::PhysicsType::transport();
  trans.transport(
    get_xi, get_mu, get_lambda, get_Ddiag, T_vapor, gpv.rho_fluid,
    Y_skin.data(), Ddiag.data(), mu_skin, xi_skin, lambda_skin, trans_parm);
  mu_skin *= SPU.mu_conv;
  lambda_skin *= SPU.lambda_conv;
  amrex::Real dy_i = diff_cent - ht_film; // Distance from film surface to cell center
  // Determine the mass evaporation values
  amrex::Real m_dot = 0.;
  amrex::Real qvap = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    Ddiag[fspec] *= mw_vap * gpv.invmw[fspec] * SPU.rhod_conv;
    // TODO: Ensure condensed mass does not exceed gas phase mass
    mi_dot[spf] = -Ddiag[fspec] / (1. - Y_vapor[spf]) *
                  (gpv.Y_fluid[fspec] - Y_vapor[spf]) / dy_i;
    m_dot += mi_dot[spf];
    qvap += mi_dot[spf] * L_fuel[spf];
    gpv.fluid_Y_dot[spf] += mi_dot[spf];
    gpv.fluid_eng_src += mi_dot[spf] * h_skin[fspec];
  }
  gpv.fluid_mass_src = m_dot;
  // Determine the temperature at wall film surface using energy balance
  amrex::Real fs1 = lambda_skin * ht_film;
  amrex::Real fs2 = lambda_film * dy_i;
  amrex::Real T_s = (fs1 * T_i + fs2 * fdat.wall_T - ht_film * dy_i * qvap) / (fs1 + fs2);
  amrex::Real qconv = lambda_skin * (T_i - T_s) / dy_i;
  gpv.fluid_eng_src += qconv;
  // Update film height
  amrex::Real new_mass = pmass + part_dt * m_dot;
  amrex::Real new_rho = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    amrex::Real newY = (Y_film[spf] * pmass + mi_dot[spf] * part_dt) / new_mass;
    p.rdata(SPI.pstateY + spf) = newY;
    new_rho += newY / fdat.rho[spf];
  }
  new_rho = 1. / new_rho;
  amrex::Real new_vol = new_mass / new_rho;
  p.rdata(SPI.pstateVol) = new_vol;
  // Adjust the height by assuming radius remains unchanged
  amrex::Real new_ht = 4. * new_vol / (M_PI * dia_film * dia_film);
  p.rdata(SPI.pstateHt) = new_ht;
  p.rdata(SPI.pstateT) = -0.5 * (T_s + fdat.wall_T);
  if (new_ht < ht_tol) {
    p.id() = -1;
  }
}
#endif // #if 0
#endif // _WALLFUNCTIONS_H_
