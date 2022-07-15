#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::GpuArray<amrex::Real, NUM_SOOT_MOMENTS + 1> soot_vals;
  amrex::Real P_mean = 101325.0;
  amrex::Real init_vel = 0.;
  amrex::Real T0 = 300.0;
  amrex::Real rho0 = 0.0;
  amrex::Real Y_N2 = 0.6693174609128006;
  amrex::Real Y_O2 = 0.20333443771796805;
  amrex::Real Y_Fuel = 0.12734810136923122;
  int fuelIndx = -1;
};

#endif
