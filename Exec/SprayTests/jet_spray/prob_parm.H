#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real jet_vel;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real jet_dia;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real part_mean_dia;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real part_stdev_dia;
extern AMREX_GPU_DEVICE_MANAGED unsigned int inject_N;
extern amrex::Gpu::ManagedVector<amrex::Real>* inject_time;
extern amrex::Gpu::ManagedVector<amrex::Real>* inject_mass;
extern amrex::Gpu::ManagedVector<amrex::Real>* inject_vel;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_time;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_mass;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_vel;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real part_rho;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real part_temp;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real mass_flow_rate;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Y_O2;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Y_N2;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real jet_start_time;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real jet_end_time;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real spray_angle;

extern std::string input_file;
} // namespace ProbParm

#endif
