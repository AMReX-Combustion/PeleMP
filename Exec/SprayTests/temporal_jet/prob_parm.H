#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 900.0;
  amrex::Real v0 = 1.;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
  amrex::Real L_jet = 0.01;
};

struct ProbParmHost : amrex::Gpu::Managed
{
  int partNum = 8000;
  amrex::Real partTemp = 300.;
  amrex::Real partDia = 0.002;
  amrex::Real partVel = 0.;
  amrex::Real velFluct = 1.E-2;
  amrex::Real velSmooth = 7.5E-4;
};

#endif
