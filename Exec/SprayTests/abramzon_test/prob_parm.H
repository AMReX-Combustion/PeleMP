#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  amrex::Real rho0 = 0.0;
  amrex::Real v0 = 0.0;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
};

struct ProbParmHost : amrex::Gpu::Managed
{
};

#endif
