#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real equil_start_time = 0.;
  amrex::Real equil_end_time = 0.;
  amrex::Real init_vel = 0.;
  amrex::Real P_mean = 101325.0;
  amrex::Real T_ref = 300.;
  amrex::Real fuel_T = 300.;
  amrex::Real oxid_T = 300.;
  amrex::Real jet_center = 0.;
  amrex::Real fuel_dia = 0.012;
  amrex::Real oxid_dia = 0.055;
  amrex::Real fuel_vel = 0.034;
  amrex::Real oxid_vel = 0.309;
  amrex::Real oxid_N2 = 0.767;
  amrex::Real oxid_O2 = 0.233;
  int fuel_indx = -1;
  amrex::Real S_equil = 0.;
  amrex::Real YF_equil = 1.;
  amrex::Real YO_equil = 1.;
  amrex::Real smooth_b = 0.004;
  unsigned int N_equil = 0;
  amrex::Real* Z_equil = nullptr;
  amrex::Real* T_equil = nullptr;
  amrex::Real* Y_equil = nullptr;
  amrex::Real* soot_equil = nullptr;
};

#endif
