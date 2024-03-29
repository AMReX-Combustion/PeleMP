#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include "pelelm_prob_parm.H"
#include "PMF.H"
#include "PMFData.H"
#include "PeleLM_parm.H"

#include "NavierStokesBase.H"
#include "IndexDefines.H"
#include "PelePhysics.H"

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
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{

  auto eos = pele::physics::PhysicsType::eos();
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};

  state(i, j, k, DEF_Temp) = prob_parm.T0;

  molefrac[O2_ID] = 0.21;
  molefrac[N2_ID] = 0.79;

  eos.X2Y(molefrac.data(), massfrac.data());

  state(i, j, k, Xvel) = 0;
  state(i, j, k, Yvel) = 0.;
#if (AMREX_SPACEDIM == 3)
  state(i, j, k, Zvel) = 0.;
#endif

  amrex::Real rho_cgs, P_cgs;
  P_cgs = prob_parm.P_mean * 10.0;

  eos.PYT2R(P_cgs, massfrac.data(), state(i, j, k, DEF_Temp), rho_cgs);
  state(i, j, k, Density) = rho_cgs * 1.0e3; // CGS -> MKS conversion

  eos.TY2H(state(i, j, k, DEF_Temp), massfrac.data(), state(i, j, k, DEF_RhoH));
  state(i, j, k, DEF_RhoH) = state(i, j, k, DEF_RhoH) * 1.0e-4 *
                             state(i, j, k, Density); // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, DEF_first_spec + n) = massfrac[n] * state(i, j, k, Density);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  amrex::Real s_ext[DEF_NUM_STATE],
  const int /*idir*/,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& prob_parm,
  ACParm const& /*ac_parm*/,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  auto eos = pele::physics::PhysicsType::eos();

  if (sgn == 1) {

    s_ext[DEF_Temp] = prob_parm.T0;

    molefrac[O2_ID] = 0.21;
    molefrac[N2_ID] = 0.79;
    eos.X2Y(molefrac.data(), massfrac.data());

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac.data(), s_ext[DEF_Temp], rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[DEF_Temp], massfrac.data(), RhoH_temp);
    s_ext[DEF_RhoH] =
      RhoH_temp * 1.0e-4 * s_ext[Density]; // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[DEF_first_spec + n] = massfrac[n] * s_ext[Density];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc(
  int /*i*/,
  int /*j*/,
  int /*k*/,
  amrex::Array4<amrex::Real> const& /*beta*/,
  amrex::GeometryData const& /*geomdata*/,
  amrex::Box const& /*dBox*/,
  const int /*dir*/,
  const int /*state_comp*/,
  const int /*nComp*/)
{
}
#endif
