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

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  amrex::Real fuel_mf = 0.;
  amrex::Real o2_mf = 0.233;
  amrex::Real n2_mf = 0.767;
  amrex::Real T_init = prob_parm.T_ref;
  amrex::Real ign_mf = 0.;
  if (prob_parm.hot_spot_temp > 0.) {
    amrex::Real smb = prob_parm.hot_spot_smooth;
    amrex::Real gaussR = std::abs(y - prob_parm.hot_spot_center);
    amrex::Real modv = std::exp(-gaussR * gaussR / (smb * smb));
    T_init += (prob_parm.hot_spot_temp - T_init) * modv;
    ign_mf = prob_parm.ign_Y_start * modv;
    o2_mf *= (1. - ign_mf);
    n2_mf = 1. - (ign_mf + o2_mf);
  }
  amrex::Real jet_vel = prob_parm.init_vel;
  massfrac[N2_ID] = n2_mf;
  massfrac[O2_ID] = o2_mf;
  massfrac[prob_parm.fuel_indx] = fuel_mf;
  massfrac[prob_parm.ign_indx] = ign_mf;

  state(i, j, k, DEF_Temp) = T_init;
  AMREX_D_TERM(state(i, j, k, Xvel) = 0.;
               , state(i, j, k, Yvel) = jet_vel;, state(i, j, k, Zvel) = 0.;);
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

AMREX_GPU_DEVICE
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
    const amrex::Real& fuel_dia = prob_parm.fuel_dia;
    const amrex::Real& oxid_dia = prob_parm.oxid_dia;
    const amrex::Real& wall_thick = prob_parm.wall_thick;
    amrex::Real len = prob_hi[0] - prob_lo[0];
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
    amrex::Real xc = prob_parm.jet_center;
    amrex::GpuArray<amrex::Real, 3> xcl;
    xcl[0] = oxid_dia / 2.;
    xcl[1] = fuel_dia / 2. + wall_thick;
    xcl[2] = fuel_dia / 2.;
    amrex::Real xcx = x[0] - xc;
    amrex::Real crad = std::abs(xcx);
#if AMREX_SPACEDIM == 3
    amrex::Real xcz = x[2] - xc;
    crad = std::sqrt(xcx * xcx + xcz * xcz);
#endif
    amrex::Real fuel_Y = 1. - prob_parm.ign_Y;
    amrex::Real ign_mf = 0.;
    amrex::Real fuel_mf = 0.;
    amrex::Real o2_mf = 0.233;
    amrex::Real n2_mf = 0.767;
    amrex::Real jet_vel = prob_parm.ext_vel;
    amrex::Real temp = prob_parm.T_ref;
    if (crad <= prob_parm.fuel_dia / 2.) {
      o2_mf = 0.;
      n2_mf = 0.;
      ign_mf = prob_parm.ign_Y;
      fuel_mf = fuel_Y;
      jet_vel = prob_parm.fuel_vel;
      temp = prob_parm.fuel_T;
    } else if (crad <= prob_parm.fuel_dia / 2. + prob_parm.wall_thick) {
      jet_vel = 0.;
      temp = prob_parm.fuel_T;
    } else if (crad <= prob_parm.oxid_dia / 2.) {
      jet_vel = prob_parm.oxid_vel;
      temp = prob_parm.oxid_T;
    }
    massfrac[N2_ID] = n2_mf;
    massfrac[O2_ID] = o2_mf;
    massfrac[prob_parm.fuel_indx] = fuel_mf;
    massfrac[prob_parm.ign_indx] = ign_mf;
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
#ifdef SOOT_MODEL
    for (int is = 0; is < DEF_NUM_SOOT_VARS; ++is) {
      s_ext[DEF_first_soot + is] = prob_parm.soot_vals[is];
    }
#endif
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

#endif
