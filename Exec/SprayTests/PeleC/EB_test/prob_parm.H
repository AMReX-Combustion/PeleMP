#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  amrex::Real rho0 = 0.0;
  amrex::Real v0 = 0.0;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
};

struct ProbParmHost
{
  amrex::Real jet_vel = 5000.;
  amrex::Real jet_dia = 1.E-2;
  amrex::Real part_mean_dia = 1.E-3;
  amrex::Real part_stdev_dia = 0.;
  amrex::Real mass_flow_rate = 2.3;
  amrex::Real part_temp = 300.;
  amrex::Real jet_start_time = 0.;
  amrex::Real jet_end_time = 10000.;
  amrex::Real spray_angle = 20.;
  amrex::Real jet_dx_mod = 50.;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> jet_cent = {{0.0}};
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_jet = {{0.0}};
};

#endif
