#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Y_O2;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Y_N2;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Y_Fuel;
extern AMREX_GPU_DEVICE_MANAGED int fuelIndx;
} // namespace ProbParm

#endif
