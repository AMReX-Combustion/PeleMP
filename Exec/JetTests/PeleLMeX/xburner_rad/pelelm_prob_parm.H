#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
#ifdef PELELM_USE_SOOT
  bool solve_soot = true;
  amrex::GpuArray<amrex::Real, NUM_SOOT_MOMENTS + 1> soot_vals;
#endif
  amrex::Real init_vel = 0.;
  amrex::Real P_mean = 101325.0;
  amrex::Real T_ref = 300.;
  amrex::Real fuel_T = 300.;
  amrex::Real oxid_T = 300.;
  amrex::Real wall_T = 300.;
  amrex::Real jet_center = 0.;
  amrex::Real wall_thick = 0.002;
  amrex::Real fuel_dia = 0.012;
  amrex::Real oxid_dia = 0.055;
  amrex::Real fuel_vel = 0.034;
  amrex::Real oxid_vel = 0.309;
  amrex::Real ext_vel = 0.;
  amrex::Real vel_time = 0.;
  amrex::Real ign_Y = 0.;
  amrex::Real S_equil = 1.;
  // Whether to use a developed velocity profile or flat
  int oxid_flat = 1;
  int fuel_flat = 1;
  int fuel_indx = -1;
  int ign_indx = -1;
  int norm_dir = 1;
};

#endif
