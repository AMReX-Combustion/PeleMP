#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata(
  int i,
  int j,
  int k,
  int /*is_incompressible*/,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& /*aux*/,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

  auto eos = pele::physics::PhysicsType::eos();
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac = {{0.0}};
  molefrac[N2_ID] = 0.79;
  molefrac[O2_ID] = 0.21;
  eos.X2Y(molefrac.data(), massfrac.data());
  state(i, j, k, TEMP) = prob_parm.T_ox;

  state(i, j, k, VELX) = 0.0;
#if (AMREX_SPACEDIM == 2)
  state(i, j, k, VELY) = 0.0;
#elif (AMREX_SPACEDIM == 3)
  state(i, j, k, VELY) = 0.0;
  state(i, j, k, VELZ) = 0.0;
#endif

  amrex::Real P_cgs = prob_parm.P_mean * 10.0;

  // Density
  amrex::Real rho_cgs = 0.0;
  eos.PYT2R(P_cgs, massfrac.data(), state(i, j, k, TEMP), rho_cgs);
  state(i, j, k, DENSITY) = rho_cgs * 1.0e3;

  // Enthalpy
  amrex::Real h_cgs = 0.0;
  eos.TY2H(state(i, j, k, TEMP), massfrac.data(), h_cgs);
  state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);

  // Species mass
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int /*m_nAux*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  auto eos = pele::physics::PhysicsType::eos();
  if (sgn == 1) {
    amrex::Real diffx = x[0] - prob_parm.jet_cents[0];
    amrex::Real diffy = x[1] - prob_parm.jet_cents[1];
    amrex::Real r2 = diffx * diffx + diffy * diffy;
    amrex::Real jr2 = prob_parm.gas_jet_dia * prob_parm.gas_jet_dia / 4.;
    if (r2 < jr2) {
      amrex::Real v_loc = prob_parm.gas_jet_vel * std::exp(-r2 * 1e8);
      amrex::Real cos_th = diffx / (sqrt(r2) + 1e-16);
      amrex::Real sin_th = diffy / (sqrt(r2) + 1e-16);
      s_ext[VELX] = 0.0 * cos_th * v_loc;
      s_ext[VELY] = 0.2 * sin_th * v_loc;
      s_ext[VELZ] = v_loc;
    }
    s_ext[TEMP] = prob_parm.T_ox;
    molefrac[O2_ID] = 0.21;
    molefrac[N2_ID] = 0.79;
    eos.X2Y(molefrac.data(), massfrac.data());

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac.data(), s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac.data(), RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& beta,
  amrex::GeometryData const& geomdata,
  amrex::Box const& domainBox,
  const int dir,
  const int beta_comp,
  const int nComp)
{
  amrex::ignore_unused(
    i, j, k, beta, geomdata, domainBox, dir, beta_comp, nComp);
  // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
  // otherwise this routine could be called for other face diffusivity (Temp,
  // velocity, ...)
}
#endif
