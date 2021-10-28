#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0;
  amrex::Real T0 = 500.;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
  amrex::Real vel = 0.;
  amrex::Real jet_vel = 50.;
  amrex::Real jet_dia = 1.E-4;
  amrex::Real part_mean_dia = 1.E-5;
  amrex::Real part_stdev_dia = 0.;
  amrex::Real mass_flow_rate = 2.3E-3;
  amrex::Real part_temp = 300.;
  amrex::Real jet_start_time = 0.;
  amrex::Real jet_end_time = 10000.;
  amrex::Real spray_angle = 20.;
  amrex::Real jet_dx_mod = 50.;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_jet = {{0.0}};
  unsigned int num_jets = 0;
  amrex::Real* jet_cents;
};

#endif
