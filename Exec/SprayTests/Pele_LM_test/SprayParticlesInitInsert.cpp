
#include "SprayParticles.H"
#include <AMReX_Particles.H>
#include <pelelm_prob.H>

using namespace amrex;

IntVect unflatten_particles (const Long idx, const IntVect& max_parts)
{
  IntVect indx;
  Long cidx = idx;
#if AMREX_SPACEDIM > 1
#if AMREX_SPACEDIM > 2
  indx[2] = cidx/(max_parts[0]*max_parts[1]);
  cidx -= indx[2]*max_parts[0]*max_parts[1];
#endif
  indx[1] = cidx/max_parts[0];
#endif
  indx[0] = cidx % max_parts[0];
  return indx;
}

IntVect unflatten_boxes (const Long     idx,
                         const IntVect& box_strides,
                         const IntVect& cell_strides)
{
  int cpb = AMREX_D_TERM(cell_strides[0],*cell_strides[1],*cell_strides[2]);
  Long cidx = idx;
  Long cblock = cidx/cpb;
  cidx -= cblock*cpb;
  IntVect curblocks = unflatten_particles(cblock, box_strides);
  IntVect local = unflatten_particles(cidx, cell_strides);
  local += curblocks*cell_strides;
  return local;
}

bool
SprayParticleContainer::insertParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level,
                                        ProbParm const& prob_parm)
{
  return false;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level,
                                        ProbParm const& prob_parm)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(ProbParm const& prob_parm)
{
  const int lev = 0;
  // Reset the particle ID in case this is during the initial iterations
  ParticleType::NextID(1);
  const int MyProc = ParallelDescriptor::MyProc();
  const int NProcs = ParallelDescriptor::NProcs();
  const int IOProc = ParallelDescriptor::IOProcessorNumber();
  Real part_rho = prob_parm.partRho;
  Real part_dia = prob_parm.partDia;
  Real T_ref = prob_parm.partTemp;
  const int pstateVel = m_sprayIndx.pstateVel;
  const int pstateDia = m_sprayIndx.pstateDia;
  const int pstateT = m_sprayIndx.pstateT;
  const int pstateRho = m_sprayIndx.pstateRho;
  const int pstateY = m_sprayIndx.pstateY;
  const IntVect num_part = prob_parm.partNum;
  const RealVect part_vel = prob_parm.partVel;
  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();
  const auto phi = Geom(lev).ProbHiArray();
  IntVect box_stride(AMREX_D_DECL(32,32,32));
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    while ((num_part[dir] % box_stride[dir]) != 0) {
      box_stride[dir] -= 2;
    }
  }
  const Long total_part_num = AMREX_D_TERM(num_part[0],*num_part[1],*num_part[2]);
  const IntVect cell_stride = num_part/box_stride;
  const Box& boxDom = Geom(lev).Domain();
  const Real len = phi[0] - plo[0];
  const RealVect dx_part(AMREX_D_DECL(len/Real(num_part[0]),
				      len/Real(num_part[1]),
				      len/Real(num_part[2])));
  Long parts_pp = total_part_num / NProcs;
  Long cur_parts_pp = parts_pp;
  if (MyProc == 0) cur_parts_pp += (total_part_num % NProcs);
  const Long first_part = (NProcs - MyProc - 1)*parts_pp;
  ParticleLocData pld;
  std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > host_particles;
#ifdef USE_SPRAY_SOA
  std::map<std::pair<int, int>, std::array<Gpu::HostVector<Real>, NAR_SPR > > host_real_attribs;
#endif
  for (int prc = first_part; prc != first_part + cur_parts_pp; ++prc) {
    //IntVect indx = unflatten_particles(prc, num_part);
    IntVect indx = unflatten_boxes(prc, box_stride, cell_stride);
    ParticleType p;
    p.id() = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();
    for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) p.pos(dir) = (Real(indx[dir]) + 0.5)*dx_part[dir];
    std::pair<int, int> ind(pld.m_grid, pld.m_tile);
#ifdef USE_SPRAY_SOA
    for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
      host_real_attribs[ind][pstateVel+dir].push_back(part_vel[dir]);
    host_real_attribs[ind][pstateT].push_back(T_ref);
    host_real_attribs[ind][pstateDia].push_back(part_dia);
    host_real_attribs[ind][pstateRho].push_back(part_rho);
    host_real_attribs[ind][pstateY].push_back(1.);
    for (int spf = 1; spf != SPRAY_FUEL_NUM; ++spf)
      host_real_attribs[ind][pstateY+spf].push_back(0.);
#else
    for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
      p.rdata(pstateVel+dir) = part_vel[dir];
    p.rdata(pstateT) = T_ref; // temperature
    p.rdata(pstateDia) = part_dia; // diameter
    p.rdata(pstateRho) = part_rho; // liquid fuel density
    for (int sp = 0; sp != SPRAY_FUEL_NUM; ++sp)
      p.rdata(pstateY+sp) = 0.;
    p.rdata(pstateY) = 1.; // Only use the first fuel species
#endif
    host_particles[ind].push_back(p);
  }
  for (auto& kv : host_particles) {
    auto grid = kv.first.first;
    auto tile = kv.first.second;
    const auto& src_tile = kv.second;
    auto& dst_tile = GetParticles(lev)[std::make_pair(grid,tile)];
    auto old_size = dst_tile.GetArrayOfStructs().size();
    auto new_size = old_size + src_tile.size();
    dst_tile.resize(new_size);

    // Copy the AoS part of the host particles to the GPU
    Gpu::copy(Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
              dst_tile.GetArrayOfStructs().begin() + old_size);
#ifdef USE_SPRAY_SOA
    for (int i = 0; i != NAR_SPR; ++i) {
      Gpu::copy(Gpu::hostToDevice,
                host_real_attribs[std::make_pair(grid,tile)][i].begin(),
                host_real_attribs[std::make_pair(grid,tile)][i].end(),
                dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
    }
#endif
  }
  Redistribute();
  Gpu::streamSynchronize();
}
