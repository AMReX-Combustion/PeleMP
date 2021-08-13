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
  PmfData const* pmf_data)
{

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  amrex::Real x1 = (x - prob_parm.standoff - 0.5 * dx[0]) * 100.;
  amrex::Real x2 = (x - prob_parm.standoff + 0.5 * dx[0]) * 100.;
  PMF::pmf(pmf_data, x1, x2, pmf_vals);
  state(i, j, k, DEF_Temp) = pmf_vals[1];
  amrex::Real norm = 0.;
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = amrex::max(0., amrex::min(1., pmf_vals[3 + n]));
    norm += massfrac[n];
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    massfrac[n] = massfrac[n] / norm;
  }
  AMREX_D_TERM(state(i, j, k, Xvel) = pmf_vals[0] * 1.E-2;
               , state(i, j, k, Yvel) = 0.;, state(i, j, k, Zvel) = 0.;);
  amrex::Real rho_cgs, P_cgs;
  P_cgs = prob_parm.P_mean * 10.;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2R(P_cgs, massfrac.data(), state(i, j, k, DEF_Temp), rho_cgs);
  state(i, j, k, Density) = rho_cgs * 1.0e3; // CGS -> MKS conversion

  eos.TY2H(state(i, j, k, DEF_Temp), massfrac.data(), state(i, j, k, DEF_RhoH));
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
  PmfData const* pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  if (sgn == 1) {
    PMF::pmf(pmf_data, prob_lo[idir], prob_lo[idir], pmf_vals);

    if (!ac_parm.ctrl_active || time <= 0.0) {
      AMREX_D_TERM(s_ext[Xvel] = pmf_vals[0] * 1.E-2;, s_ext[Yvel] = 0.0;
                   , s_ext[Zvel] = 0.0;);
    } else {
      AMREX_D_TERM(s_ext[Xvel] = ac_parm.ctrl_V_in +
                                 (time - ac_parm.ctrl_tBase) * ac_parm.ctrl_dV;
                   , s_ext[Yvel] = 0.0;, s_ext[Zvel] = 0.0;);
    }

    s_ext[DEF_Temp] = pmf_vals[1];

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = pmf_vals[3 + n];
    }

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(P_cgs, massfrac, s_ext[DEF_Temp], rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[DEF_Temp], massfrac, RhoH_temp);
    s_ext[DEF_RhoH] =
      RhoH_temp * 1.0e-4 * s_ext[Density]; // CGS -> MKS conversion

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
#endif
