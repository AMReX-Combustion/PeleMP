#ifndef NFLDMUPDATE_H
#define NFLDMUPDATE_H

#include "SelfSimFuncs.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool
useNFLDM(const amrex::RealVect pos, const amrex::Real dia_jet)
{
  amrex::Real radial_loc0 = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  if (pos[2] / dia_jet > 60. || radial_loc0 / dia_jet > 3.8) {
    return false;
  }
  return true;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
updateNFLDM(
  SprayParticleContainer::ParticleType& p,
  const amrex::Real dia_part,
  const amrex::Real U_jet,
  const amrex::Real dia_jet,
  const amrex::Real dt,
  const amrex::Real rho_part,
  const amrex::Real mu_fluid,
  const amrex::Real pmass,
  amrex::RealVect& pos,
  amrex::RealVect& vel_part,
  amrex::RealVect& fluid_mom_src)
{
  amrex::Real radial_loc0 = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  amrex::Real rhoAC = SelfSimRho(dt);
  // Azimuthal angle relative to the jet center
  amrex::Real theta = std::atan2(pos[1], pos[0]);
  const amrex::Real Nran = p.rdata(SprayComps::pstateNran);
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
  amrex::Real axial_loc = pos[2] + dt * VAxial;
  amrex::Real radial_loc = radial_loc0 + dt * VpRadial;
  amrex::Real z_norm = axial_loc / dia_jet;
  amrex::Real r_norm = radial_loc / dia_jet;
  amrex::Real rh_norm = SelfSimNormRadius(z_norm, r_norm);
  amrex::Real U0 = SelfSimU0(z_norm, U_jet);
  amrex::Real Uz = SelfSimUz(z_norm, rh_norm, r_norm, U0);
  amrex::Real dvpz_dt = (Uz - VAxial) / tau_p;
  p.rdata(SprayComps::pstateVAxial) = VAxial + dt * dvpz_dt;
  p.rdata(SprayComps::pstateUURadial) = SelfSimUrpUrp(z_norm, rh_norm, U0);
  amrex::RealVect new_vel = {
    VpRadial * std::cos(theta), VpRadial * std::sin(theta),
    p.rdata(SprayComps::pstateVAxial)};
  fluid_mom_src += pmass * dvpz_dt;
  vel_part = new_vel;
  pos[0] = radial_loc * std::cos(theta);
  pos[1] = radial_loc * std::sin(theta);
  pos[2] = axial_loc;
}

#endif
