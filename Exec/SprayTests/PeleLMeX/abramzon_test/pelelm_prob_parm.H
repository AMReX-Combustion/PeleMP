#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean   = 1013250.0_rt;
   amrex::Real T0 = 1500.;
   amrex::Real Y_O2 = 0.233;
   amrex::Real Y_N2 = 0.767;
   amrex::Real vel = 0.;
};
#endif