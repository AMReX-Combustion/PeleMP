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
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& prob_parm,
  PmfData const* /*pmf_data*/)
{

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};

  AMREX_D_TERM(state(i, j, k, Xvel) = prob_parm.init_vel;
               , state(i, j, k, Yvel) = 0.;, state(i, j, k, Zvel) = 0.;);

  massfrac[N2_ID] = prob_parm.Y_N2;
  massfrac[O2_ID] = prob_parm.Y_O2;
  massfrac[prob_parm.fuelIndx] = prob_parm.Y_Fuel;
  state(i, j, k, DEF_Temp) = prob_parm.T0;
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
  ACParm const& /*ac_parm*/,
  PmfData const* /*pmf_data*/)
{
  AMREX_D_TERM(s_ext[Xvel] = 0.;, s_ext[Yvel] = 0.;, s_ext[Zvel] = 0.;);

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  massfrac[N2_ID] = prob_parm.Y_N2;
  massfrac[O2_ID] = prob_parm.Y_O2;
  massfrac[prob_parm.fuelIndx] = prob_parm.Y_Fuel;
  s_ext[DEF_Temp] = prob_parm.T0;
  amrex::Real rho_cgs, P_cgs;
  P_cgs = prob_parm.P_mean * 10.;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2R(P_cgs, massfrac.data(), prob_parm.T0, rho_cgs);
  amrex::Real rho_mks = rho_cgs * 1.0e3; // CGS -> MKS conversion
  s_ext[Density] = rho_mks;
  amrex::Real rhoH;
  eos.TY2H(prob_parm.T0, massfrac.data(), rhoH);
  s_ext[DEF_RhoH] = rhoH * 1.0e-4 * rho_mks;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[DEF_first_spec + n] = massfrac[n] * rho_mks;
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
