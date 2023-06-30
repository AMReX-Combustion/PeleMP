
#ifndef WALLFUNCTIONS_H
#define WALLFUNCTIONS_H

#include "Drag.H"
#include "SprayFuelData.H"
#include "SprayInterpolation.H"

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
  const amrex::RealVect pos,
  const amrex::IntVect bflags,
  const amrex::IntVect ijkc,
#ifdef AMREX_USE_EB
  const amrex::IntVect ijkc_prev,
  const amrex::RealVect vel_part,
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
    // If particle has moved into a cut-cell
    if (flags(ijkc).isSingleValued()) {
      wall_check = true;
      // Use the normal and boundary centroid from the new cell
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        normal[dir] = -bnorm(ijkc, dir);
        bcentv[dir] =
          static_cast<amrex::Real>(ijkc[dir]) + 0.5 + bcent(ijkc, dir);
      }
      // If volume fraction is too small, treat cell face as the wall to prevent
      // particle from residing in cell
      if (vfrac(ijkc) < min_eb_vfrac) {
        // Direction of face to treat as wall
        int fdir = normal.maxDir(true);
        amrex::Real face_sgn = (bcent(ijkc, fdir) < 0.) ? 1. : -1.;
        normal = amrex::RealVect::TheZeroVector();
        normal[fdir] = -face_sgn;
        bcentv[fdir] =
          static_cast<amrex::Real>(ijkc[fdir]) + 0.5 * (1. - face_sgn);
      }
      // If the particle has moved into a covered cell
    } else if (flags(ijkc).isCovered()) {
      wall_check = true;
      // If the particle came from a cut-cell, check if it is behind that cells
      // EB
      amrex::Real test_par = 2.;
      if (flags(ijkc_prev).isSingleValued()) {
        test_par = 0.;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          normal[dir] = -bnorm(ijkc_prev, dir);
          bcentv[dir] = static_cast<amrex::Real>(ijkc_prev[dir]) + 0.5 +
                        bcent(ijkc_prev, dir);
          test_par += (pos[dir] - bcentv[dir]) * normal[dir];
        }
      }
      // If it didn't come from a cut-cell or isn't behind the EB of the
      // previous cell, more tests are required
      if (test_par > 0.) {
        amrex::Real max_diff = 0.;
        amrex::RealVect cur_bcent;
        amrex::RealVect cur_norm;
        // Determine the distance from the particle to the nearest wall
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          if (ijkc_prev[dir] != ijkc[dir]) {
            amrex::Real cur_diff = 0.;
            amrex::IntVect tsti = ijkc_prev;
            tsti[dir] = ijkc[dir];
            // If droplet passes through adjacent cell with EB, find distance to
            // EB
            if (flags(tsti).isSingleValued()) {
              amrex::Real diff1 = 0.;
              amrex::Real diff2 = 0.;
              for (int okdir = 0; okdir < AMREX_SPACEDIM; ++okdir) {
                cur_norm[okdir] = -bnorm(tsti, okdir);
                cur_bcent[okdir] = static_cast<amrex::Real>(tsti[okdir]) + 0.5 +
                                   bcent(tsti, okdir);
                diff1 += (pos[okdir] - cur_bcent[okdir]) * cur_norm[okdir];
                diff2 += vel_part[okdir] * cur_norm[okdir];
              }
              if (diff2 != 0.) {
                cur_diff = diff1 / diff2;
              }
              // If droplet passes through another covered cell, find distance
              // to normal cell face
            } else if (flags(tsti).isCovered()) {
              for (int okdir = 0; okdir < AMREX_SPACEDIM; ++okdir) {
                cur_bcent[okdir] = static_cast<amrex::Real>(tsti[okdir]) + 0.5;
                cur_norm[okdir] = 0.;
              }
              amrex::Real sgn = (ijkc_prev[dir] > ijkc[dir]) ? 1. : -1.;
              cur_bcent[dir] += sgn * 0.5;
              cur_norm[dir] = sgn;
              cur_diff = (pos[dir] - cur_bcent[dir]) / vel_part[dir];
            }
            if (cur_diff > max_diff) {
              max_diff = cur_diff;
              normal = cur_norm;
              bcentv = cur_bcent;
            }
          }
        } // for (int dir = 0;...
      }   // if (test_par > 0.)
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
Checks and reflects particles off BC or EB. If particle reflects too many times,
it places the particle near the wall sets the normal velocity component to zero.
*/
AMREX_GPU_DEVICE
AMREX_INLINE
void
impose_wall(
  SprayParticleContainer::ParticleType& p,
  const amrex::RealVect& dx,
  const amrex::RealVect& plo,
  const amrex::RealVect& phi,
  const amrex::IntVect bndry_lo,
  const amrex::IntVect bndry_hi,
  amrex::IntVect& bflags,
  const bool use_EB,
#ifdef AMREX_USE_EB
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<amrex::Real const> const& bcent,
  amrex::Array4<amrex::Real const> const& bnorm,
  amrex::Array4<amrex::Real const> const& vfrac,
  const amrex::Real min_eb_vfrac,
#endif
  const amrex::IntVect ijkc,
  const amrex::IntVect ijkc_prev)
{
  amrex::Real par_dot = 2.;
  amrex::RealVect normal;
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SprayComps::pstateVel), p.rdata(SprayComps::pstateVel + 1),
    p.rdata(SprayComps::pstateVel + 2)));
  amrex::RealVect normphi;
  amrex::RealVect normpos;
  // Normalize pos
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    normpos[dir] = (p.pos(dir) - plo[dir]) / dx[dir];
    normphi[dir] = (phi[dir] - plo[dir]) / dx[dir];
  }
  // Check if particle is behind a wall
  bool do_refl = check_wall(
    normpos, bflags, ijkc,
#ifdef AMREX_USE_EB
    ijkc_prev, vel_part, use_EB, flags, bcent, bnorm, vfrac, min_eb_vfrac,
#endif
    par_dot, normal);
  if (do_refl) {
    amrex::Real Nw_Vp = normal.dotProduct(vel_part);
    vel_part -= 2. * Nw_Vp * normal;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.pos(dir) -= 2. * par_dot * normal[dir];
    }
  }
}
#endif // _WALLFUNCTIONS_H
