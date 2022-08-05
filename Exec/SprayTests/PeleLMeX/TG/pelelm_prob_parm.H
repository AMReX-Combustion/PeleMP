#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real T0 = 500.;
  amrex::Real reynolds = 1600.;
  amrex::Real Ugas = 0.;
  amrex::IntVect partNum = amrex::IntVect(AMREX_D_DECL(100, 100, 100));
  amrex::Real partTemp = 300.;
  amrex::Real partDia = 1.E-3;
};
#endif
