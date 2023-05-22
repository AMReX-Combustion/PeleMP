#ifndef SELFSIMFUNCS_H
#define SELFSIMFUNCS_H

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
SelfSimRho(const amrex::Real dt)
{
  if (dt < 6.E-8) {
    return 0.9998 +
           dt * (6.315E4 + dt * (-5.851E13 + dt * (3.858E20 - 8.613E26 * dt)));
  } else {
    return 1.04 +
           dt *
             (-2.885E6 +
              dt * (-7.414E12 +
                    dt * (5.362E19 +
                          dt * (-1.108E26 + dt * (1.001E32 - 3.371E37 * dt)))));
  }
}

// Return local normalized half radius
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
SelfSimNormRadius(const amrex::Real z_norm, const amrex::Real radial_loc)
{
  if (z_norm <= 10.) {
    return radial_loc / (1.23E-6 * z_norm + 3.875E-5);
  } else {
    return radial_loc / (2.652E-6 * z_norm + 2.453E-5);
  }
}

// Return local centerline (r=0) axial velocity
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
SelfSimU0(const amrex::Real z_norm, const amrex::Real U_jet)
{
  if (z_norm >= 15.) {
    return 441.3 + z_norm * (2.47 + z_norm * (-0.09747 + z_norm * 6.661E-4));
  } else {
    return U_jet;
  }
}

// Return local mean axial velocity
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
SelfSimUz(
  const amrex::Real z_norm,
  const amrex::Real rh_norm,
  const amrex::Real r_norm,
  const amrex::Real U0,
  const amrex::Real dia_jet)
{
  if (rh_norm < 3.4) {
    // Flat hat weighting
    amrex::Real weight =
      1. +
      rh_norm *
        (-0.04867 +
         rh_norm *
           (-1.026 +
            rh_norm *
              (0.8134 + rh_norm * (-0.2779 +
                                   rh_norm * (0.04556 - 2.918E-3 * rh_norm)))));
    if (z_norm >= 25.) {
      return weight * U0;
    } else {
      amrex::Real weight1 = 1.;
      if (r_norm > 0.5) {
        weight1 = 0.;
      }
      return U0 * (z_norm / 25. * weight + (1. - z_norm / 25.) * weight1);
    }
  } else {
    return 0.;
  }
}

// Return local gas phase covariance fluctuating radial velocity
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
SelfSimUrpUrp(
  const amrex::Real z_norm, const amrex::Real rh_norm, const amrex::Real U0)
{
  const amrex::Real entryLength = 18.;
  if (z_norm > entryLength || rh_norm > 3.35) {
    amrex::Real vvNorm =
      2.501E-3 +
      rh_norm *
        (-6.764E-3 +
         rh_norm *
           (0.038 +
            rh_norm *
              (-0.04557 +
               rh_norm *
                 (0.02381 +
                  rh_norm *
                    (-6.358E-3 + rh_norm * (8.549E-4 - 6.606E-5 * rh_norm))))));
    vvNorm = amrex::max(0., vvNorm);
    if (rh_norm > 1.) {
      return vvNorm * U0 * U0;
    } else {
      amrex::Real vvNormFit0 =
        2.804E-5 +
        z_norm *
          (-1.828E-6 +
           z_norm * (-9.18E-7 + z_norm * (7.928E-8 - 7.537E-10 * z_norm)));
      vvNormFit0 = amrex::max(0., vvNormFit0);
      amrex::Real rh_norm2 = std::pow(rh_norm, 2);
      vvNorm = (1. - rh_norm2) * vvNormFit0 + vvNorm * rh_norm2;
      return U0 * U0 * vvNorm;
    }
  } else {
    return 0.;
  }
}

#endif
