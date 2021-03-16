
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

using namespace amrex;

bool
SprayParticleContainer::insertParticles(
                                        Real time, Real dt, int nstep, int lev, int finest_level, 
                                        ProbParmHost const& prob_parm),
  ProbParmDevice const& prob_parm_d)
{
  return false;
}

bool
SprayParticleContainer::injectParticles(
  Real time,
  Real dt,
  int nstep,
  int lev,
  int finest_level,
  ProbParmHost const& prob_parm,
  ProbParmDevice const& prob_parm_d)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(
  AmrLevel* pelec,
  const int& lev,
  const int& num_ppc,
  ProbParmHost const& prob_parm,
  ProbParmDevice const& prob_parm_d)
{
  const int numGrids = pelec->numGrids();
  const int MyProc = ParallelDescriptor::MyProc();
  const int NProcs = ParallelDescriptor::NProcs();
  const int IOProc = ParallelDescriptor::IOProcessorNumber();
  Real part_dia = prob_parm.partDia;
  const int pstateVel = m_sprayIndx.pstateVel;
  const int pstateT = m_sprayIndx.pstateT;
  const int pstateDia = m_sprayIndx.pstateDia;
  const int pstateY = m_sprayIndx.pstateY;
  const SprayData* fdat = m_sprayData;
  Real rho_part = fdat->rho[0];
  Real T_ref = prob_parm.partTemp;
  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();
  const auto phi = Geom(lev).ProbHiArray();
  const Real jet_len = prob_parm_d.L_jet;
  RealVect jet_loc(AMREX_D_DECL(
    plo[0] + 0.5 * (phi[0] - plo[0]), plo[1],
    plo[2] + 0.5 * (phi[2] - plo[2])));
  const Long total_part_num = prob_parm.partNum;
  // Length in the jet normal direction
  const Real len = phi[1] - plo[1];
  // Number of particles per processor
  Long parts_pp = total_part_num / NProcs;
  // Number of particles for this processor
  Long cur_parts_pp = parts_pp;
  // The first processor gets the left-overs
  if (MyProc == 0)
    cur_parts_pp += (total_part_num % NProcs);
  // First particle index for the current processor
  const int first_part = (NProcs - MyProc - 1) * parts_pp;
  const Real len_pp = Real(parts_pp) / Real(total_part_num) * len;
  const Real cur_len_pp = Real(cur_parts_pp) / Real(total_part_num) * len;
  const Real start_loc = plo[1] + (NProcs - MyProc - 1) * len_pp;
  ParticleLocData pld;
  std::map<std::pair<int, int>, Gpu::HostVector<ParticleType>> host_particles;
#ifdef USE_SPRAY_SOA
  std::map<std::pair<int, int>, std::array<Gpu::HostVector<Real>, NAR_SPR>>
    host_real_attribs;
#endif
  for (int prc = first_part; prc != first_part + cur_parts_pp; ++prc) {
    ParticleType p;
    p.id() = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();
    AMREX_D_TERM(p.pos(0) = jet_loc[0] + (amrex::Random() - 0.5) * jet_len;
                 , p.pos(1) = start_loc + amrex::Random() * cur_len_pp;
                 , p.pos(2) = jet_loc[2] + (amrex::Random() - 0.5) * jet_len;);
    std::pair<int, int> ind(pld.m_grid, pld.m_tile);
    RealVect cpartVel(AMREX_D_DECL(
      prob_parm.velFluct * (amrex::Random() - 0.5),
      prob_parm.partVel + prob_parm.velFluct * (amrex::Random() - 0.5),
      prob_parm.velFluct * (amrex::Random() - 0.5)));
#ifdef USE_SPRAY_SOA
    for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
      host_real_attribs[ind][pstateVel + dir].push_back(cpartVel[dir]);
    host_real_attribs[ind][pstateT].push_back(T_ref);
    host_real_attribs[ind][pstateDia].push_back(part_dia);
    host_real_attribs[ind][pstateY].push_back(1.);
    for (int spf = 1; spf != SPRAY_FUEL_NUM; ++spf)
      host_real_attribs[ind][pstateY + spf].push_back(0.);
#else
    for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
      p.rdata(pstateVel + dir) = cpartVel[dir];
    p.rdata(pstateT) = T_ref;      // temperature
    p.rdata(pstateDia) = part_dia; // diameter
    for (int sp = 0; sp != SPRAY_FUEL_NUM; ++sp)
      p.rdata(pstateY + sp) = 0.;
    p.rdata(pstateY) = 1.; // Only use the first fuel species
#endif
    host_particles[ind].push_back(p);
  }
  for (auto& kv : host_particles) {
    auto grid = kv.first.first;
    auto tile = kv.first.second;
    const auto& src_tile = kv.second;
    auto& dst_tile = GetParticles(lev)[std::make_pair(grid, tile)];
    auto old_size = dst_tile.GetArrayOfStructs().size();
    auto new_size = old_size + src_tile.size();
    dst_tile.resize(new_size);

    // Copy the AoS part of the host particles to the GPU
    Gpu::copy(
      Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
      dst_tile.GetArrayOfStructs().begin() + old_size);
#ifdef USE_SPRAY_SOA
    for (int i = 0; i != NAR_SPR; ++i) {
      Gpu::copy(
        Gpu::hostToDevice,
        host_real_attribs[std::make_pair(grid, tile)][i].begin(),
        host_real_attribs[std::make_pair(grid, tile)][i].end(),
        dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
    }
#endif
  }
  Redistribute();
  Gpu::streamSynchronize();
}
