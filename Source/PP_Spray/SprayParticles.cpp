
#include "SprayParticles.H"
#include <AMReX_Particles.H>
#include <AMReX_ParticleReduce.H>
#ifdef SPRAY_PELE_LM
#include "PeleLM.H"
#endif
#include "Transport.H"
#include "Drag.H"

using namespace amrex;

void
SprayParticleContainer::init_bcs()
{
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (phys_bc->lo(dir) == Symmetry   ||
        phys_bc->lo(dir) == SlipWall   ||
        phys_bc->lo(dir) == NoSlipWall) {
      reflect_lo[dir] = true;
    } else {
      reflect_lo[dir] = false;
    }
    if (phys_bc->hi(dir) == Symmetry   ||
        phys_bc->hi(dir) == SlipWall   ||
        phys_bc->hi(dir) == NoSlipWall) {
      reflect_hi[dir] = true;
    } else {
      reflect_hi[dir] = false;
    }
  }
}

void
SprayParticleContainer::moveKick (MultiFab&   state,
                                  MultiFab&   source,
                                  const int   level,
                                  const Real& dt,
                                  const Real  time,
                                  const bool  isVirtualPart,
                                  const bool  isGhostPart,
                                  const int   state_ghosts,
                                  const int   source_ghosts,
                                  MultiFab*   u_mac)
{
  bool do_move = false;
  int width = 0;
  moveKickDrift(state, source, level, dt, time, isVirtualPart, isGhostPart,
                state_ghosts, source_ghosts, do_move, width, u_mac);
}

void
SprayParticleContainer::moveKickDrift (MultiFab&   state,
                                       MultiFab&   source,
                                       const int   level,
                                       const Real& dt,
                                       const Real  time,
                                       const bool  isVirtualPart,
                                       const bool  isGhostPart,
                                       const int   state_ghosts,
                                       const int   source_ghosts,
                                       const bool  do_move,
                                       const int   where_width,
                                       MultiFab*   u_mac)
{
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  AMREX_ASSERT(u_mac == nullptr || u_mac[0].nGrow() >= 1);
  AMREX_ASSERT(level >= 0);
  AMREX_ASSERT(state.nGrow() >= 2);

  //If there are no particles at this level
  if (level >= this->GetParticles().size())
    return;

  const Real strttime = ParallelDescriptor::second();

  BL_PROFILE_VAR("SprayParticles::updateParticles()", UPD_PART);
  updateParticles(level, state, source, dt, time, state_ghosts, source_ghosts, do_move, u_mac);
  BL_PROFILE_VAR_STOP(UPD_PART);

  // Fill ghost cells after we've synced up ..
  // TODO: Check to see if this is needed at all
  // if (level > 0)
  //   source.FillBoundary(Geom(level).periodicity());

  // ********************************************************************************

  // ********************************************************************************

  if (this->m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;
    ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor()) {
      if (do_move) {
        Print() << "SprayParticleContainer::moveKickDrift() time: "
                << stoptime << '\n';
      } else {
        Print() << "SprayParticleContainer::moveKick() time: "
                << stoptime << '\n';
      }
    }
  }
}

Real
SprayParticleContainer::estTimestep (int level, Real cfl) const
{
  BL_PROFILE("ParticleContainer::estTimestep()");
  AMREX_ASSERT(m_setFuelData);
  // TODO: Clean up this mess and bring the num particle functionality back
  Real dt = std::numeric_limits<Real>::max();
  if (level >= this->GetParticles().size() ||
      m_sprayIndx[SprayComps::mom_tran] == 0)
    return -1.;

  const Real strttime = ParallelDescriptor::second();
  const Geometry& geom = this->m_gdb->Geom(level);
  const auto dx = Geom(level).CellSizeArray();
  const auto dxi = Geom(level).InvCellSizeArray();
  {
    amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MyParConstIter pti(*this, level); pti.isValid(); ++pti) {
      const AoS& pbox = pti.GetArrayOfStructs();
      const ParticleType* pstruct = pbox().data();
      const Long n = pbox.numParticles();
#ifdef USE_SPRAY_SOA
      auto& attribs = pti.GetAttribs();
      AMREX_D_TERM(const Real * up = attribs[0].data();,
                   const Real * vp = attribs[1].data();,
                   const Real * wp = attribs[2].data(););
#endif
      reduce_op.eval(n, reduce_data,
                     [=] AMREX_GPU_DEVICE (const Long i) -> ReduceTuple
      {
        const ParticleType& p = pstruct[i];
        // TODO: This assumes that pstateVel = 0 and dxi[0] = dxi[1] = dxi[2]
        if (p.id() > 0) {
#ifdef USE_SPRAY_SOA
          const Real max_mag_vdx = amrex::max(AMREX_D_DECL(std::abs(up[i]),
                                                           std::abs(vp[i]),
                                                           std::abs(wp[i])))*dxi[0];
#else
          const Real max_mag_vdx = amrex::max(AMREX_D_DECL(std::abs(p.rdata(0)),
                                                           std::abs(p.rdata(1)),
                                                           std::abs(p.rdata(2))))*dxi[0];
#endif
          Real dt_part = (max_mag_vdx > 0.) ? (cfl/max_mag_vdx) : 1.E50;
#ifdef SPRAY_PELE_LM
          // Conversion since particle velocities are in cm and dx is in m for PeleLM
          dt_part *= 100.;
#endif
          return dt_part;
        }
        return 1.E50;
      });
    }
    ReduceTuple hv = reduce_data.value();
    Real ldt_cpu = amrex::get<0>(hv);
    dt = amrex::min(dt,ldt_cpu);
  }
  ParallelDescriptor::ReduceRealMin(dt);
  // Check if the velocity of particles being injected
  // is greater existing particle velocities
  if (m_injectVel > 0.)
    dt = amrex::min(dt, cfl*dx[0]/m_injectVel);

  if (this->m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;
    ParallelDescriptor::ReduceRealMax(stoptime,
                                      ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor())
      std::cout << "SprayParticleContainer::estTimestep() time: "
                << stoptime << '\n';
  }
  return dt;
}

void
SprayParticleContainer::updateParticles(const int&  level,
                                        MultiFab&   state,
                                        MultiFab&   source,
                                        const Real& flow_dt,
                                        const Real& time,
                                        const int   state_ghosts,
                                        const int   source_ghosts,
                                        const bool  do_move,
                                        MultiFab*   u_mac)
{
  AMREX_ASSERT(m_setFuelData);
  AMREX_ASSERT(OnSameGrids(level, state));
  AMREX_ASSERT(OnSameGrids(level, source));
  const auto dxi = this->Geom(level).InvCellSizeArray();
  const auto dx = this->Geom(level).CellSizeArray();
  const auto plo = this->Geom(level).ProbLoArray();
  const auto phi = this->Geom(level).ProbHiArray();
  const auto domain = this->Geom(level).Domain();
  IntVect dom_lo = domain.smallEnd();
  IntVect dom_hi = domain.bigEnd();
  IntVect lo_bound;
  IntVect hi_bound;
  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
    if (!this->Geom(level).isPeriodic(dir)) {
      if (reflect_lo[dir]) lo_bound[dir] = 1;
      else lo_bound[dir] = -1;
      if (reflect_hi[dir]) hi_bound[dir] = 1;
      else hi_bound[dir] = -1;
    } else {
      lo_bound[dir] = 0;
      hi_bound[dir] = 0;
    }
    // Only concerned with reflective boundaries
    if (!reflect_lo[dir]) dom_lo[dir] -= 100;
    if (!reflect_hi[dir]) dom_hi[dir] += 100;
  }
  const Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
  const Real inv_vol = 1./vol;
  // Set all constants
  const Real num_ppp = m_parcelSize;
  const Real C_eps = 1.E-15;
  const Real B_eps = 1.E-7;
  const Real dia_eps = 2.E-6;
  const int nSubMax = 100;
  const Real third = 1./3.;
  const Real rule = third; // This determines how average values for the vapor are approximated
  const Real Pi_six = M_PI/6.;
  Real mw_fluid[NUM_SPECIES];
  Real invmw[NUM_SPECIES];
  EOS::molecular_weight(mw_fluid);
  EOS::inv_molecular_weight(invmw);
  // Extract control parameters for mass, heat, and momentum transfer
  const int heat_trans = m_sprayIndx[SprayComps::heat_tran];
  const int mass_trans = m_sprayIndx[SprayComps::mass_tran];
  const int mom_trans = m_sprayIndx[SprayComps::mom_tran];
  const Real inv_Ru = 1./EOS::RU;
  const Real ref_T = m_sprayRefT;
  // Particle components indices
  const int pstateVel = m_sprayIndx[SprayComps::pstateVel];
  const int pstateT   = m_sprayIndx[SprayComps::pstateT];
  const int pstateRho = m_sprayIndx[SprayComps::pstateRho];
  const int pstateDia = m_sprayIndx[SprayComps::pstateDia];
  const int pstateY   = m_sprayIndx[SprayComps::pstateY];
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  if (!mass_trans && !heat_trans) {
    get_Ddiag = false;
    get_lambda = false;
  }
  Real vel_conv = 1.;
  Real pos_conv = 1.;
  Real rho_conv = 1.;
  Real eng_conv = 1.;
  Real mom_src_conv = 1.;
  Real mass_src_conv = 1.;
  Real eng_src_conv = 1.;
#ifdef SPRAY_PELE_LM
  vel_conv = 100.;  // Turn m/s to cm/s
  rho_conv = 0.001; // Turn kg/m^3 to g/cm^3
  pos_conv = 0.01;  // Turn cm to m for updating position
  eng_conv = 1.E4; // For converting enthalpy to CGS
  // This makes no sense, conversions should be independent
  // of dimensions but numerical tests show this isn't the case
#if AMREX_SPACEDIM == 2
  mom_src_conv = 1.E-3;
  mass_src_conv = 1.E-1;
  eng_src_conv = 1.E-5;
#elif AMREX_SPACEDIM == 3
  mom_src_conv = 1.E-5;
  mass_src_conv = 1.E-3;
  eng_src_conv = 1.E-7;
#endif
#endif

  // Component indices for conservative state
  const int rhoIndx = m_sprayIndx[SprayComps::rhoIndx];
  const int momIndx = m_sprayIndx[SprayComps::momIndx];
  const int engIndx = m_sprayIndx[SprayComps::engIndx];
  const int utempIndx = m_sprayIndx[SprayComps::utempIndx];
  const int specIndx = m_sprayIndx[SprayComps::specIndx];
  // Start the ParIter, which loops over separate sets of particles in different boxes
  for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
    const Box& tile_box = pti.tilebox();
    const Box& state_box = pti.growntilebox(state_ghosts);
    const Box& src_box = pti.growntilebox(source_ghosts);
    const Long Np = pti.numParticles();
    ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);
    // Get particle attributes if StructOfArrays are used
#ifdef USE_SPRAY_SOA
    auto& attribs = pti.GetAttribs();
    Real * velp[AMREX_SPACEDIM];
    AMREX_D_TERM(velp[0] = attribs[pstateVel].data();,
                 velp[1] = attribs[pstateVel+1].data();,
                 velp[2] = attribs[pstateVel+2].data(););
    Real * Tp = attribs[pstateT].data();
    Real * diap = attribs[pstateDia].data();
    Real * rhop = attribs[pstateRho].data();
    std::array<Real *, SPRAY_FUEL_NUM> Yp;
    for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
      Yp[spf] = attribs[pstateY+spf].data();
    }
#endif
    SprayData fdat = m_fuelData.getSprayData();
    Array4<const Real> const& statearr = state.array(pti);
    Array4<Real> const& sourcearr = source.array(pti);
// #ifdef SPRAY_PELE_LM
//     GpuArray<
//       Array4<const Real>, AMREX_SPACEDIM> const
//       umac{AMREX_D_DECL(u_mac[0].array(pti), u_mac[1].array(pti), u_mac[2].array(pti))};
// #endif
    AMREX_FOR_1D ( Np, i,
    {
      ParticleType& p = pstruct[i];
      if (p.id() > 0) {
        Real dt = flow_dt;
        Real sub_source = inv_vol;
        // TODO: I was hoping to not have to instantiate this everytime
        Real Y_fluid[NUM_SPECIES];
        Real Y_skin[NUM_SPECIES];
        Real h_skin[NUM_SPECIES];
        Real cp_n[NUM_SPECIES];
        Real mass_frac[NUM_SPECIES];
        Real Ddiag[NUM_SPECIES];
        Real B_M_num[SPRAY_FUEL_NUM];
        Real Sh_num[SPRAY_FUEL_NUM];
        Real Y_dot[SPRAY_FUEL_NUM];
        Real L_fuel[SPRAY_FUEL_NUM];
        // Weights for interpolation
        Real coef[AMREX_D_PICK(2, 4, 8)];
        // Indices of adjacent cells
        IntVect indx_array[AMREX_D_PICK(2, 4, 8)];
        // Cell-based length of particle in domain
        RealVect len(AMREX_D_DECL((p.pos(0) - plo[0])*dxi[0] + 0.5,
                                  (p.pos(1) - plo[1])*dxi[1] + 0.5,
                                  (p.pos(2) - plo[2])*dxi[2] + 0.5));
        RealVect vel_fluid(RealVect::TheZeroVector());
// #ifdef SPRAY_PELE_LM
//         InterpolateFaceVelocity(len, dom_lo, dom_hi, umac, vel_fluid);
// #endif
        Real T_fluid = 0.;
        Real rho_fluid = 0.;
        for (int sp = 0; sp != NUM_SPECIES; ++sp) Y_fluid[sp] = 0.;
        // Do initial interpolation
        AdjIndexWeights(len, indx_array, coef, dom_lo, dom_hi);
        // Extract adjacent values and interpolate fluid at particle location
        for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
          IntVect cur_indx = indx_array[aindx];
#ifdef AMREX_DEBUG
          if (!state_box.contains(cur_indx))
            Abort("SprayParticleContainer::updateParticles() -- state box too small");
#endif
          Real cur_coef = coef[aindx];
          Real cur_rho = statearr(cur_indx, rhoIndx);
          rho_fluid += cur_coef*cur_rho*rho_conv;
          Real inv_rho = 1./cur_rho;
          for (int sp = 0; sp != NUM_SPECIES; ++sp) {
            int mf_indx = sp + specIndx;
            Real cur_mf = statearr(cur_indx, mf_indx)*inv_rho;
            Y_fluid[sp] += cur_coef*cur_mf;
            mass_frac[sp] = cur_mf;
          }
#ifdef SPRAY_PELE_LM
          inv_rho = 1.; // Since velocity is provided instead of momentum
#endif
          AMREX_D_TERM(Real velx = statearr(cur_indx, momIndx)*inv_rho*vel_conv;
                       Real ke = 0.5*velx*velx;
                       vel_fluid[0] += cur_coef*velx;,
                       Real vely = statearr(cur_indx, momIndx+1)*inv_rho*vel_conv;
                       ke += 0.5*vely*vely;
                       vel_fluid[1] += cur_coef*vely;,
                       Real velz = statearr(cur_indx, momIndx+2)*inv_rho*vel_conv;
                       ke += 0.5*velz*velz;
                       vel_fluid[2] += cur_coef*velz;);
          Real T_val = statearr(cur_indx, utempIndx);
#ifdef SPRAY_PELE_LM
          // Not sure if this is needed
//           Real intH = statearr(cur_indx, engIndx)/cur_rho*eng_conv;
//           EOS::HY2T(intH, mass_frac, T_val);
#else
          Real intEng = statearr(cur_indx, engIndx)*inv_rho - ke;
          EOS::EY2T(intEng, mass_frac, T_val);
#endif
          T_fluid += cur_coef*T_val;
        }
        int isub = 1; // Initialize the number of sub-cycles
        int nsub = 1; // This is set in the first run through the loop
        while (isub <= nsub) {
#ifdef USE_SPRAY_SOA
          RealVect vel_part(AMREX_D_DECL(velp[0][i], velp[1][i], velp[2][i]));
          Real T_part = Tp[i];
          Real dia_part = diap[i];
          Real rho_part = rhop[i];
#else
          RealVect vel_part(AMREX_D_DECL(p.rdata(pstateVel),
                                         p.rdata(pstateVel+1),
                                         p.rdata(pstateVel+2)));
          Real T_part = p.rdata(pstateT);
          Real dia_part = p.rdata(pstateDia);
          Real rho_part = p.rdata(pstateRho);
#endif
          Real dia2_part = dia_part*dia_part;
          Real pmass = Pi_six*rho_part*dia_part*dia2_part;
          Real part_ke = 0.5*vel_part.radSquared();
          // If multiple sub-cycle iterations are needed, we might need to
          // re-interpolate values at the particle
          // However, since we don't allow the particle to move more than
          // 5% of a cell (particle cfl = 0.05), the adjacent fluid values
          // should remain the same
          // Model the fuel vapor using the one-third rule
          Real delT = amrex::max(T_fluid - T_part, 0.);
          Real T_skin = T_part + rule*delT;
          // Calculate the C_p at the skin temperature for each species
          EOS::T2Cpi(T_skin, cp_n);
          EOS::T2Hi(T_part, h_skin);
          Real mw_mix = 0.;  // Average molar mass of gas mixture
          if (heat_trans || mass_trans) {
            for (int sp = 0; sp != NUM_SPECIES; ++sp) {
              mw_mix += Y_fluid[sp]*invmw[sp];
              Y_skin[sp] = 0.;
            }
          } else {
            for (int sp = 0; sp != NUM_SPECIES; ++sp) {
              mw_mix += Y_fluid[sp]*invmw[sp];
              Y_skin[sp] = Y_fluid[sp];
            }
          }
          Real p_fluid = rho_fluid*EOS::RU*mw_mix*T_fluid;
          mw_mix = 1./mw_mix;

          // Solve for state of the vapor and mass transfer coefficient B_M
          Real sumYSkin = 0.; // Mass fraction of the fuel in skin film, uses one-thirds rule
          Real sumYFuel = 0.; // Mass fraction of the fuel in the gas phase
          Real cp_skin = 0.; // Averaged C_p at particle surface
          Real cp_L_av = 0.; // Cp of the liquid state
          Real mw_vap = 0.; // Average molar mass of vapor mixture
          if (heat_trans || mass_trans) {
            for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
              const int fspec = fdat.indx(spf);
              const Real mw_fuel = mw_fluid[fspec];
              // Compute latent heat
              Real part_latent = h_skin[fspec] + fdat.latent(spf)
                - fdat.cp(spf)*(T_part - ref_T);
              L_fuel[spf] = part_latent;
              // Compute the mass fraction of the fuel vapor at droplet surface
              Real pres_sat = EOS::PATM*std::exp(part_latent*inv_Ru*mw_fuel*
                                                 (1./fdat.boilT(spf) - 1./T_part)) + C_eps;
              Real Yfv = mw_fuel*pres_sat/(mw_mix*p_fluid + (mw_fuel - mw_mix)*pres_sat);
              Yfv = amrex::max(0., amrex::min(1. - C_eps, Yfv));
              B_M_num[spf] = amrex::max(C_eps, (Yfv - Y_fluid[fspec])/(1. - Yfv));
              Y_skin[fspec] = Yfv + rule*(Y_fluid[fspec] - Yfv);
              sumYSkin += Y_skin[fspec];
#ifdef USE_SPRAY_SOA
              cp_L_av += Yp[spf][i]*fdat.cp(spf);
#else
              cp_L_av += p.rdata(pstateY+spf)*fdat.cp(spf);
#endif
              sumYFuel += Y_fluid[fspec];
            }
            const Real restYSkin = 1. - sumYSkin;
            for (int sp = 0; sp != NUM_SPECIES; ++sp) {
              Y_skin[sp] += restYSkin*Y_fluid[sp];
              cp_skin += Y_skin[sp]*cp_n[sp];
              mw_vap += Y_skin[sp]*invmw[sp];
            }
            mw_vap = 1./mw_vap;
          }
          Real lambda_skin = 0.;
          Real mu_skin = 0.;
          Real xi_skin = 0.;
          transport(get_xi, get_mu, get_lambda, get_Ddiag,
                    T_skin, rho_fluid, Y_skin, Ddiag,
                    mu_skin, xi_skin, lambda_skin);
          // Ensure gas is not all fuel to allow evaporation
          bool evap_fuel = (sumYFuel >= 1.) ? false : true;
          RealVect diff_vel = vel_fluid - vel_part;
          Real diff_vel_mag = diff_vel.vectorLength();
          // Local Reynolds number
          Real Reyn = rho_fluid*diff_vel_mag*dia_part/mu_skin;
          Real Nu_0 = 1.;

          // Solve mass transfer source terms
          Real m_dot = 0.;
          Real d_dot = 0.;
          if ((mass_trans || heat_trans) && evap_fuel) {
            Real Pr_skin = mu_skin*cp_skin/lambda_skin;
            Real powR = amrex::max(std::pow(Reyn, 0.077), 1.);
            Nu_0 = 1. + powR*std::cbrt(1. + Reyn*Pr_skin);
            for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
              const int fspec = fdat.indx(spf);
              // Convert mass diffusion coefficient from mixture average
              // to binary for fuel only, not concerned with other species
              Ddiag[fspec] *= mw_vap*invmw[fspec];
              const Real rhoD = Ddiag[fspec];
              const Real Sc_skin = mu_skin/rhoD;
              const Real B_M = B_M_num[spf];
              Real logB = std::log(1. + B_M);
              // Calculate Sherwood number and evaporation rate
              Real invFM = B_M/(logB*std::pow(1. + B_M, 0.7));
              Real Sh_0 = 1. + powR*std::cbrt(1. + Reyn*Sc_skin);
              Sh_num[spf] = 2. + (Sh_0 - 2.)*invFM;
              if (mass_trans) {
                Y_dot[spf] = -amrex::max(M_PI*rhoD*dia_part*Sh_num[spf]*logB, 0.);
                m_dot += Y_dot[spf];
              }
            }
            d_dot = m_dot/(0.5*M_PI*rho_part*dia2_part);
            Real inv_tau_d = -m_dot/(3.*pmass);
            if (isub == 1)
              nsub = amrex::min(int(flow_dt*inv_tau_d) + 1, nSubMax);
          }

          // Solve for momentum source terms
          const Real inv_pmass = 1./pmass;
          RealVect fluid_mom_src(RealVect::TheZeroVector());
          RealVect part_mom_src(RealVect::TheZeroVector());
          Real fluid_eng_src = 0.;
          if (mom_trans) {
            Real drag_coef = 0.;
            if (Reyn > 0.)
	      drag_coef =
	        (Reyn > 1.) ? 24./Reyn*(1. + std::cbrt(Reyn*Reyn)/6.) : 24./Reyn;
            Real drag_force = 0.125*rho_fluid*drag_coef*M_PI*dia2_part*diff_vel_mag;
            part_mom_src = drag_force*diff_vel;
            fluid_mom_src = part_mom_src + vel_part*m_dot;
            // s_d,mu dot u_d
            Real S_dmu_dot_u = part_mom_src.dotProduct(vel_part);
#ifndef SPRAY_PELE_LM
            fluid_eng_src += S_dmu_dot_u + m_dot*part_ke;
#endif
            Real inv_tau_var = drag_force*inv_pmass;
            if (isub == 1)
              nsub = amrex::min(amrex::max(nsub, int(flow_dt*inv_tau_var) + 1), nSubMax);
          }

          // Solve for energy source terms
          Real part_temp_src = 0.;
          if (heat_trans && evap_fuel) {
            const Real inv_pm_cp = inv_pmass/cp_L_av;
            Real coeff_heat = 0.;
            for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
              const int fspec = fdat.indx(spf);
              Real ratio = cp_n[fspec]*Sh_num[spf]*Ddiag[fspec]/lambda_skin;
              Real heatC = calcHeatCoeff(ratio, B_M_num[spf], B_eps, C_eps, Nu_0);
              // Convection term
              coeff_heat += heatC;
              fluid_eng_src += Y_dot[spf]*h_skin[fspec];
              part_temp_src += Y_dot[spf]*L_fuel[spf];
            }
            Real conv_src = M_PI*lambda_skin*dia_part*delT*coeff_heat;
            fluid_eng_src += conv_src;
            part_temp_src += conv_src;
            part_temp_src *= inv_pm_cp;
            if (isub == 1 && delT > C_eps) {
              Real inv_tau_T = conv_src*inv_pm_cp/delT;
              nsub = amrex::min(amrex::max(nsub, int(flow_dt*inv_tau_T) + 1), nSubMax);
            }
          }
          if (isub == 1) {
            sub_source /= Real(nsub);
            dt = flow_dt/Real(nsub);
          }
          const Real part_dt = 0.5*dt;
          if (mom_trans || mass_trans || heat_trans) {
            bool remove_particle = false;
            // Modify particle velocity by half time step
#ifdef USE_SPRAY_SOA
            if (mom_trans) {
              AMREX_D_TERM(Gpu::Atomic::Add(&velp[0][i], part_dt*part_mom_src[0]*inv_pmass);,
                           Gpu::Atomic::Add(&velp[1][i], part_dt*part_mom_src[1]*inv_pmass);,
                           Gpu::Atomic::Add(&velp[2][i], part_dt*part_mom_src[2]*inv_pmass););
              if (do_move) {
                for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
                  Gpu::Atomic::Add(&p.pos(dir), dt*velp[dir][i]*pos_conv);
                  remove_particle = checkWall(p.pos(dir), velp[dir][i],
                                              phi[dir], plo[dir], hi_bound[dir], lo_bound[dir]);
                }
              }
            }
            if (heat_trans) {
              // Modify particle temperature
              Gpu::Atomic::Add(&Tp[i], part_dt*part_temp_src);
            }
            if (mass_trans) {
              // Compute new particle diameter
              Real new_dia = dia_part + part_dt*d_dot;
              if (new_dia < dia_eps) {
                remove_particle = true;
              } else {
                diap[i] = new_dia;
              }
            }
#else
            if (mom_trans) {
              AMREX_D_TERM(Gpu::Atomic::Add(&p.rdata(pstateVel), part_dt*part_mom_src[0]*inv_pmass);,
                           Gpu::Atomic::Add(&p.rdata(pstateVel+1), part_dt*part_mom_src[1]*inv_pmass);,
                           Gpu::Atomic::Add(&p.rdata(pstateVel+2), part_dt*part_mom_src[2]*inv_pmass););
              // Modify particle position by whole time step
              if (do_move) {
                for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
                  Gpu::Atomic::Add(&p.pos(dir), dt*p.rdata(pstateVel+dir)*pos_conv);
                  remove_particle = checkWall(p.pos(dir), p.rdata(pstateVel+dir),
                                              phi[dir], plo[dir], hi_bound[dir], lo_bound[dir]);
                }
              }
            }
            if (heat_trans) {
              // Modify particle temperature
              Gpu::Atomic::Add(&p.rdata(pstateT), part_dt*part_temp_src);
            }
            if (mass_trans) {
              // Compute new particle diameter
              Real new_dia = dia_part + part_dt*d_dot;
              if (new_dia < dia_eps) {
                remove_particle = true;
              } else {
                p.rdata(pstateDia) = new_dia;
              }
            }
#endif // End of check for SOA
            if (remove_particle) {
              p.id() = -1;
              isub = nsub + 1;
              continue;
            }
            for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
              Real cur_coef = -coef[aindx]*sub_source*num_ppp;
              IntVect cur_indx = indx_array[aindx];
#ifdef AMREX_DEBUG
              if (!src_box.contains(cur_indx))
                Abort("SprayParticleContainer::updateParticles() -- source box too small");
#endif
              if (mom_trans) {
                for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
                  const int nf = momIndx + dir;
                  Gpu::Atomic::Add(&sourcearr(cur_indx, nf), cur_coef*fluid_mom_src[dir]*mom_src_conv);
                }
              }
              if (mass_trans) {
                Gpu::Atomic::Add(&sourcearr(cur_indx, rhoIndx), cur_coef*m_dot*mass_src_conv);
                for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
                  const int nf = specIndx + fdat.indx(spf);
                  Gpu::Atomic::Add(&sourcearr(cur_indx, nf), cur_coef*Y_dot[spf]*mass_src_conv);
                }
              }
              Gpu::Atomic::Add(&sourcearr(cur_indx, engIndx), cur_coef*fluid_eng_src*eng_src_conv);
            }
          }
          // Increment sub-iteration
          ++isub;
        } // End of isub while loop
      } // End of p.id() > 0 check
    }); // End of loop over particles
  }
}
