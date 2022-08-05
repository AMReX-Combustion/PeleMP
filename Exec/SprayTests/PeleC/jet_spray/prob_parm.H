#ifndef PROB_PARM_H
#define PROB_PARM_H

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
  unsigned int inject_N = 0;
  amrex::Vector<amrex::Real> inject_time;
  amrex::Vector<amrex::Real> inject_mass;
  amrex::Vector<amrex::Real> inject_vel;
  std::string input_file = "";
};

#endif
