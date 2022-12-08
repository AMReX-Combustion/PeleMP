#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 1013250.0_rt;
  amrex::Real T0 = 1500.;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
  amrex::Real vel_gas = 0.;
  amrex::Real dia_drop = -1.;
  amrex::Real T_drop = 300.;
  bool set_drop_loc = false;
  amrex::RealVect vel_drop = amrex::RealVect(AMREX_D_DECL(0., 0., 0.));
  amrex::RealVect loc_drop = amrex::RealVect(AMREX_D_DECL(0., 0., 0.));
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_drop;
};
#endif
