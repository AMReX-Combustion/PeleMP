#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real T0 = 900.;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
  // Jet parameters
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_jet = {{0.0}};
  amrex::Real jet_start_time = 0.;
  amrex::Real jet_end_time = 10000.;
  amrex::Real jet_dia = 1.E-5;
  amrex::Real jet_vel = 50.;
  // This is the spread of the jet
  amrex::Real spread_angle = 20.;
  int start_inj_proc = 0;
  int num_inj_procs = 1;

  amrex::Real part_mean_dia = 1.E-6;
  amrex::Real part_stdev_dia = -1.;
  amrex::Real part_weibull_k = 3.;
  amrex::Real mass_flow_rate = 2.3E-3;
  amrex::Real part_temp = 300.;
  // Center of the jet
  amrex::RealVect jet_cent;
  amrex::RealVect jet_norm;
  bool hollow_jet = false;
};
#endif
