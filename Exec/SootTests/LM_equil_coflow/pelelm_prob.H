#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>
#include <pmf.H>
#include <pmf_data.H>
#include <PeleLM_parm.H>

#include <IndexDefines.H>
#include <PelePhysics.H>
#ifdef SOOT_MODEL
#include "SootModel.H"
#endif

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real mixture_fraction (const amrex::Real S_equil,
                              const amrex::Real Y_ff,
                              const amrex::Real Y_oo)
{
  return (S_equil * Y_ff - Y_oo + 1.) / (S_equil + 1.);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
flame(amrex::Real Z,
      amrex::Real& T,
      Real* massfrac,
      ProbParm const& prob_parm)
{
  const int N = prob_parm.N_equil;
  if (Z > 1. - 1.E-12) {
    massfrac[prob_parm.fuel_indx] = 1.;
    T = prob_parm.fuel_T;
  } else if (Z < 1.E-12) {
    massfrac[N2_ID] = prob_parm.oxid_N2;
    massfrac[O2_ID] = prob_parm.oxid_O2;
    T = prob_parm.oxid_T;
  } else {
    int loside = -1;
    int hiside = -1;
    for (unsigned int i = 0; i < N - 1; i++) {
      if (
        (Z >= prob_parm.Z_equil[i]) &&
        (Z <= prob_parm.Z_equil[i + 1])) {
        loside = i;
        hiside = i + 1;
      }
    }
    amrex::Real x1 = prob_parm.Z_equil[loside];
    amrex::Real x2 = prob_parm.Z_equil[hiside];
    amrex::Real invdx = 1. / (x2 - x1);
    amrex::Real T1 = prob_parm.T_equil[loside];
    amrex::Real T2 = prob_parm.T_equil[hiside];
    T = T1 + (T2 - T1) * invdx * (Z - x1);
    for (unsigned int j = 0; j < NUM_SPECIES; j++) {
      const amrex::Real y1 = prob_parm.Y_equil[N * j + loside];
      const amrex::Real y2 = prob_parm.Y_equil[N * j + hiside];
      amrex::Real dydx = (y2 - y1) * invdx;
      massfrac[j] = y1 + dydx * (Z - x1);
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
smooth_function(const amrex::Real b,
                const amrex::Real L,
                const amrex::Real x)
{
  return 0.5 - 0.5*std::tanh(x/(b*L));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
weight_vals(amrex::Real* weights,
            const amrex::Real b,
            const amrex::Real L,
            const amrex::Real x,
            const amrex::Real xc,
            const amrex::Real oxid_dia,
            const amrex::Real fuel_dia)
{
  amrex::Real x1 = x - xc + oxid_dia/2.;
  amrex::Real x2 = x - xc + fuel_dia/2.;
  amrex::Real x3 = x - xc - fuel_dia/2.;
  amrex::Real x4 = x - xc - oxid_dia/2.;
  weights[0] = smooth_function(b, L, x1);
  weights[1] = smooth_function(b, L, x2);
  weights[2] = smooth_function(b, L, x3);
  weights[3] = smooth_function(b, L, x4);
  return;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
profile_vals(const amrex::Real* weights,
             const amrex::Real oxid_y,
             const amrex::Real fuel_y,
             const amrex::Real ext_y)
{
  amrex::Real y1 = (ext_y - oxid_y)*weights[0] + oxid_y;
  amrex::Real y2 = (y1 - fuel_y)*weights[1] + fuel_y;
  amrex::Real y3 = (y2 - oxid_y)*weights[2] + oxid_y;
  return (y3 - ext_y)*weights[3] + ext_y;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  PmfData const* /*pmf_data*/)
{

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real len = prob_hi[0] - prob_lo[0];
  const amrex::Real leny = prob_hi[1] - prob_lo[1];

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);
  amrex::Real o2in = prob_parm.oxid_O2;
  amrex::Real n2in = prob_parm.oxid_N2;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  amrex::Real temp = prob_parm.T_ref;
  amrex::Real top_flame = 0.7 * leny;
  // Adjust the fuel diameter linearly in the Y direction to form a triangle
  amrex::Real cur_dia = prob_parm.fuel_dia * (-(y - prob_lo[1]) / top_flame + 1.);
  amrex::Real xc = prob_parm.jet_center;
  amrex::GpuArray<amrex::Real, 4> weights;
  const amrex::Real smooth_b = prob_parm.smooth_b;
  weight_vals(weights.data(), smooth_b, len, x, xc, prob_parm.oxid_dia, cur_dia);
  amrex::Real fuel_mf = profile_vals(weights.data(), 0., 1., 0.);
  amrex::Real o2_mf = profile_vals(weights.data(), o2in, 0., o2in);
  amrex::Real n2_mf = profile_vals(weights.data(), n2in, 0., n2in);
  amrex::Real Yff = fuel_mf;
  amrex::Real Yoo = o2_mf / o2in;
  amrex::Real Z = mixture_fraction(prob_parm.S_equil, Yff, Yoo);
  flame(Z, temp, massfrac.data(), prob_parm);
  state(i, j, k, DEF_Temp) = temp;
  amrex::Real init_vel = prob_parm.init_vel;
  AMREX_D_TERM(state(i, j, k, Xvel) = 0.;
               , state(i, j, k, Yvel) = init_vel;, state(i, j, k, Zvel) = 0.;);
  amrex::Real rho_cgs, P_cgs;
  P_cgs = prob_parm.P_mean * 10.;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2R(P_cgs, massfrac.data(), state(i, j, k, DEF_Temp), rho_cgs);
  state(i, j, k, Density) = rho_cgs * 1.0e3; // CGS -> MKS conversion

  eos.TY2H(
    state(i, j, k, DEF_Temp), massfrac.data(), state(i, j, k, DEF_RhoH));
  state(i, j, k, DEF_RhoH) = state(i, j, k, DEF_RhoH) * 1.0e-4 *
                             state(i, j, k, Density); // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, DEF_first_spec + n) = massfrac[n] * state(i, j, k, Density);
  }
  for (int is = 0; is < DEF_NUM_SOOT_VARS; ++is) {
    state(i, j, k, DEF_first_soot + is) = prob_parm.soot_vals[is];
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  ACParm const& ac_parm,
  PmfData const* /*pmf_data*/)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  if (sgn == 1) {
    amrex::Real o2in = prob_parm.oxid_O2;
    amrex::Real n2in = prob_parm.oxid_N2;
    amrex::Real len = prob_hi[0] - prob_lo[0];
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
    amrex::Real xc = prob_parm.jet_center;
    amrex::GpuArray<amrex::Real, 4> weights;
    const amrex::Real smooth_b = prob_parm.smooth_b;
    weight_vals(weights.data(), smooth_b, len, x[0], xc, prob_parm.oxid_dia, prob_parm.fuel_dia);
    amrex::Real fuel_mf = profile_vals(weights.data(), 0., 1., 0.);
    amrex::Real o2_mf = profile_vals(weights.data(), o2in, 0., o2in);
    amrex::Real n2_mf = profile_vals(weights.data(), n2in, 0., n2in);
    amrex::Real jet_vel = profile_vals(weights.data(), prob_parm.oxid_vel, prob_parm.fuel_vel, 0.);
    amrex::Real temp = profile_vals(weights.data(), prob_parm.oxid_T, prob_parm.fuel_T, prob_parm.T_ref);
    amrex::Real Yff = fuel_mf;
    amrex::Real Yoo = o2_mf / o2in;
    amrex::Real Z = mixture_fraction(prob_parm.S_equil, Yff, Yoo);
    flame(Z, temp, massfrac.data(), prob_parm);
    s_ext[DEF_Temp] = temp;
    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;
    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(P_cgs, massfrac.data(), s_ext[DEF_Temp], rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;
    eos.TY2H(s_ext[DEF_Temp], massfrac.data(), RhoH_temp);
    s_ext[DEF_RhoH] =
      RhoH_temp * 1.0e-4 * s_ext[Density]; // CGS -> MKS conversion
    AMREX_D_TERM(s_ext[Xvel] = 0.;,
                 s_ext[Yvel] = jet_vel;,
                 s_ext[Zvel] = 0.;);
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[DEF_first_spec + n] = massfrac[n] * s_ext[Density];
    }
    for (int is = 0; is < DEF_NUM_SOOT_VARS; ++is) {
      s_ext[DEF_first_soot + is] = prob_parm.soot_vals[is];
    }
  }
}

AMREX_GPU_DEVICE
inline void
zero_visc(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& beta,
  amrex::GeometryData const& geomdata,
  amrex::Box const& dBox,
  const int dir,
  const int state_comp,
  const int nComp)
{

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
  // const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
  // const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_add_ext_source(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& ext_src,
  const amrex::Real time,
  const amrex::Real dt,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  PmfData const* /*pmf_data*/)
{
}
#endif
