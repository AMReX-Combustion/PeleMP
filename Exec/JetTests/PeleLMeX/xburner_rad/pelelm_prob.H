#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>

#include <PeleLM_Index.H>
#include <PelePhysics.H>
#include "SootModel.H"

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
pelelm_initdata(
  int i,
  int j,
  int k,
  int /*is_incompressible*/,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& aux,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real len = prob_hi[0] - prob_lo[0];

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);
  int ndir = prob_parm.norm_dir;
  int nvelComp = VELY;
  amrex::Real normlen = y;
#if AMREX_SPACEDIM == 3
  if (ndir == 2) {
    nvelComp = VELZ;
    normlen = z;
  }
#endif
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  amrex::Real fuel_mf = 0.;
  amrex::Real o2_mf = 0.233;
  amrex::Real n2_mf = 0.767;
  amrex::Real T_init = prob_parm.T_ref;
  amrex::Real ign_mf = 0.;
  amrex::Real jet_vel = prob_parm.init_vel;
  massfrac[N2_ID] = n2_mf;
  massfrac[O2_ID] = o2_mf;
  massfrac[prob_parm.fuel_indx] = fuel_mf;
  massfrac[prob_parm.ign_indx] = ign_mf;
  state(i, j, k, TEMP) = T_init;
  AMREX_D_TERM(state(i, j, k, VELX) = 0.;
               , state(i, j, k, VELY) = 0.;, state(i, j, k, VELZ) = 0.;);
  state(i, j, k, nvelComp) = jet_vel;
  amrex::Real rho_cgs, P_cgs;
  P_cgs = prob_parm.P_mean * 10.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2R(P_cgs, massfrac.data(), state(i, j, k, TEMP), rho_cgs);
  state(i, j, k, DENSITY) = rho_cgs * 1.0e3; // CGS -> MKS conversion

  eos.TY2H(state(i, j, k, TEMP), massfrac.data(), state(i, j, k, RHOH));
  state(i, j, k, RHOH) = state(i, j, k, RHOH) * 1.0e-4 *
                             state(i, j, k, DENSITY); // CGS -> MKS conversion
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
  }
  if (prob_parm.solve_soot) {
    for (int is = 0; is < NUM_SOOT_MOMENTS + 1; ++is) {
      state(i, j, k, FIRSTSOOT + is) = prob_parm.soot_vals[is];
    }
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
  pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  if (sgn == 1 && idir == 2) {
    const amrex::Real& fuel_dia = prob_parm.fuel_dia;
    const amrex::Real& oxid_dia = prob_parm.oxid_dia;
    const amrex::Real& wall_thick = prob_parm.wall_thick;
    amrex::Real len = prob_hi[0] - prob_lo[0];
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
    amrex::Real xc = prob_parm.jet_center;
    amrex::Real xcx = x[0] - xc;
    amrex::Real crad = std::abs(xcx);
    int nvelComp = VELY;
#if AMREX_SPACEDIM == 3
    amrex::Real xcz = x[2] - xc;
    if (prob_parm.norm_dir == 2) {
      xcz = x[1] - xc;
      nvelComp = VELZ;
    }
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
      if (time < prob_parm.vel_time) {
        jet_vel = prob_parm.oxid_vel +
          time / prob_parm.vel_time * (prob_parm.fuel_vel - prob_parm.oxid_vel);
      }
      if (!prob_parm.fuel_flat) {
        jet_vel *= (1. - std::pow(2. * crad  / prob_parm.fuel_dia, 2));
      }
      temp = prob_parm.fuel_T;
    } else if (crad < prob_parm.fuel_dia / 2. + prob_parm.wall_thick) {
      jet_vel = 0.;
      temp = prob_parm.wall_T;
    } else if (crad <= prob_parm.oxid_dia / 2.) {
      jet_vel = prob_parm.oxid_vel;
      temp = prob_parm.oxid_T;
    }
    massfrac[N2_ID] = n2_mf;
    massfrac[O2_ID] = o2_mf;
    massfrac[prob_parm.fuel_indx] = fuel_mf;
    massfrac[prob_parm.ign_indx] = ign_mf;
    s_ext[TEMP] = temp;
    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;
    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(P_cgs, massfrac.data(), s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.E3;
    eos.TY2H(s_ext[TEMP], massfrac.data(), RhoH_temp);
    s_ext[RHOH] =
      RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion
    AMREX_D_TERM(s_ext[VELX] = 0.;,
                 s_ext[VELY] = 0.;,
                 s_ext[VELZ] = 0.;);
    s_ext[RHORT] = prob_parm.P_mean;
    s_ext[nvelComp] = jet_vel;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
    }
    if (prob_parm.solve_soot) {
      for (int is = 0; is < NUM_SOOT_MOMENTS + 1; ++is) {
        s_ext[FIRSTSOOT + is] = prob_parm.soot_vals[is];
      }
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif
