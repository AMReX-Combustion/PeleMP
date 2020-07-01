#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real inlet_vel;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real jet_vel;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real jet_dia;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real part_dia;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real part_rho;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real part_temp;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Y_O2;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Y_N2;
} // namespace ProbParm

#endif
