#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real mach;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Stmod;
extern AMREX_GPU_DEVICE_MANAGED bool convecting;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v0;
extern AMREX_GPU_DEVICE_MANAGED amrex::IntVect partNum;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhoRatio;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real partRho;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real partDia;
} // namespace ProbParm

#endif
