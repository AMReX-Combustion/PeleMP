#ifndef NFLDMUPDATE_H
#define NFLDMUPDATE_H

#include "SelfSimFuncs.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool
updateNFLDM(
  SprayParticleContainer::ParticleType& p,
  const amrex::Real U_jet,
  const amrex::Real dia_jet,
  const amrex::Real dt,
  const amrex::Real rho_part,
  const amrex::Real mu_fluid)
{
  amrex::Real radial_loc0 =
    std::sqrt(p.pos(0) * p.pos(0) + p.pos(1) * p.pos(1));
  if (p.pos(2) / dia_jet > 60. || radial_loc0 / dia_jet > 3.8) {
    return false;
  }
  amrex::Real rhoAC = SelfSimRho(dt);
  amrex::Real dia_part = p.rdata(SprayComps::pstateDia);
  // Azimuthal angle relative to the jet center
  amrex::Real theta_azi = std::atan2(p.pos(1), p.pos(0));
  amrex::Real Nran = p.rdata(SprayComps::PstateNran);
  amrex::Real tau_p = rho_part * dia_part * dia_part / (18. * mu_fluid);
  // Covariance of UpRadial
  amrex::Real UURadial = p.rdata(SprayComps::pstateUURadial);
  // Update fluctuating radial velocity of the gas
  p.rdata(SprayComps::pstateUpRadial) *= rhoAC;
  p.rdata(SprayComps::pstateUpRadial) +=
    std::sqrt(UURadial) * (1. - rhoAC) * Nran;
  amrex::Real UpRadial = p.rdata(SprayComps::pstateUpRadial);
  // Update fluctuating radial velocity of the particle
  p.rdata(SprayComps::pstateVpRadial) =
    UpRadial -
    (UpRadial - p.rdata(SprayComps::pstateVpRadial)) * std::exp(-dt / tau_p);
  amrex::Real VpRadial = p.rdata(SprayComps::pstateVpRadial);

  // Mean axial velocity of the particle
  amrex::Real VAxial = p.rdata(SprayComps::pstateVAxial);

  // Update radial and axial locations, not sure why some things are done out of
  // order
  amrex::Real axial_loc = p.pos(2) + dt * VAxial;
  amrex::Real radial_loc = radial_loc0 + dt * VpRadial;
  amrex::Real z_norm = axial_loc / dia_jet;
  amrex::Real r_norm = radial_loc / dia_jet;
  amrex::Real rh_norm = SelfSimNormRadius(z_norm, radial_loc);
  amrex::Real U0 = SelfSimU0(z_norm, U_jet);
  amrex::Real Uz = SelfSimUz(z_norm, rh_norm, r_norm, U0, dia_jet);
  p.rdata(SprayComps::pstateVAxial) = VAxial + dt * (Uz - VAxial) / tau_p;
  p.rdata(SprayComps::pstateUURadial) = SelfSimUrpUrp(z_norm, rh_norm, U0);
  // Update positions in Cartesian coordinates
  p.pos(0) = radial_loc * std::cos(theta);
  p.pos(1) = radial_loc * std::sin(theta);
  p.pos(2) = axial_loc;
  // Update velocities in Cartesian coordinates
  p.rdata(SprayComps::pstateVel) = VpRadial * std::cos(theta);
  p.rdata(SprayComps::pstateVel + 1) = VpRadial * std::sin(theta);
  p.rdata(SprayComps::pstateVel + 2) = p.rdata(SprayComps::pstateVAxial);
  return true;
}

#endif
