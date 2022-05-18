#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  amrex::Real rho0 = 0.0;
  amrex::Real v0 = 0.;
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
  amrex::RealVect jet_cent = amrex::RealVect::TheZeroVector();
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_jet = {{0.0}};
  unsigned int inject_N = 0;
  amrex::Vector<amrex::Real> inject_time;
  amrex::Vector<amrex::Real> inject_mass;
  amrex::Vector<amrex::Real> inject_vel;
  std::string input_file = "";
};

#endif
