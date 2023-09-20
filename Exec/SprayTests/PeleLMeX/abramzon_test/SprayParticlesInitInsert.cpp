
#include "SprayParticles.H"
#include <pelelmex_prob.H>

bool
SprayParticleContainer::injectParticles(
  amrex::Real time,
  amrex::Real dt,
  int nstep,
  int lev,
  int finest_level,
  ProbParm const& prob_parm)
{
  amrex::ignore_unused(time, dt, nstep, lev, finest_level, prob_parm);
  return false;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts, ProbParm const& prob_parm)
{
  if (!init_parts) {
    return;
  }
  const amrex::Geometry& geom = this->m_gdb->Geom(0);
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  if (amrex::ParallelDescriptor::MyProc() == 0) {
    ParticleType p;
    p.id() = ParticleType::NextID();
    p.cpu() = amrex::ParallelDescriptor::MyProc();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (prob_parm.set_drop_loc) {
        p.pos(dir) = prob_parm.loc_drop[dir];
      } else {
        p.pos(dir) = plo[dir] + (phi[dir] - plo[dir]) / 2.;
      }
      p.rdata(SprayComps::pstateVel + dir) = prob_parm.vel_drop[dir];
    }
    p.rdata(SprayComps::pstateT) = prob_parm.T_drop;
    p.rdata(SprayComps::pstateDia) = prob_parm.dia_drop;
    for (int n = 0; n < SPRAY_FUEL_NUM; ++n) {
      p.rdata(SprayComps::pstateY + n) = prob_parm.Y_drop[n];
    }
    amrex::ParticleLocData pld;
    std::map<std::pair<int, int>, amrex::Gpu::HostVector<ParticleType>>
      host_particles;
    bool where = Where(p, pld);
    if (!where) {
      amrex::Abort("Bad particle");
    }
    std::pair<int, int> ind(pld.m_grid, pld.m_tile);
    host_particles[ind].push_back(p);
    for (auto& kv : host_particles) {
      auto grid = kv.first.first;
      auto tile = kv.first.second;
      const auto& src_tile = kv.second;
      auto& dst_tile = GetParticles(0)[std::make_pair(grid, tile)];
      auto old_size = dst_tile.GetArrayOfStructs().size();
      auto new_size = old_size + src_tile.size();
      dst_tile.resize(new_size);

      // Copy the AoS part of the host particles to the GPU
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
        dst_tile.GetArrayOfStructs().begin() + old_size);
    }
  }
}
