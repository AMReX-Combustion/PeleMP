#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm : amrex::Gpu::Managed
{
  amrex::Real P_mean = 101325.0;
  amrex::Real standoff = 0.;
};

#endif
