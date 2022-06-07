#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::GpuArray<amrex::Real, NUM_SOOT_MOMENTS + 1> soot_vals;
  amrex::Real init_vel = 0.;
  amrex::Real P_mean = 101325.0;
  amrex::Real T_ref = 300.;
  amrex::Real fuel_T = 300.;
  amrex::Real oxid_T = 300.;
  amrex::Real jet_center = 0.;
  amrex::Real wall_thick = 0.002;
  amrex::Real fuel_dia = 0.012;
  amrex::Real oxid_dia = 0.055;
  amrex::Real fuel_vel = 0.034;
  amrex::Real oxid_vel = 0.309;
  amrex::Real ext_vel = 0.;
  amrex::Real vel_time = 0.;
  amrex::Real ign_Y = 0.;
  // Whether to use a developed velocity profile or flat
  int oxid_flat = 1;
  int fuel_flat = 1;
  int fuel_indx = -1;
  int ign_indx = -1;
  int norm_dir = 1;

  // These are for the hot spot used to ignite the mixture
  amrex::Real hot_spot_center = -1.;
  amrex::Real hot_spot_temp = -1.;
  amrex::Real hot_spot_smooth = 0.0055;
};

#endif
