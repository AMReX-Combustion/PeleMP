#include <PeleLM.H>
#include <pelelm_prob.H>


extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");

        pp.query("ref_p", PeleLM::prob_parm->P_mean);
        pp.query("ref_T", PeleLM::prob_parm->T0);
        pp.query("init_vel", PeleLM::prob_parm->vel);
        pp.query("init_N2", PeleLM::prob_parm->Y_N2);
        pp.query("init_O2", PeleLM::prob_parm->Y_O2);
        pp.query("num_particles", PeleLM::prob_parm->partNum);
        std::array<amrex::Real, AMREX_SPACEDIM> pvel;
        pp.query<amrex::Real>("part_vel", pvel);
        for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
          PeleLM::prob_parm->partVel[dir] = pvel[dir];
        pp.get("part_rho", PeleLM::prob_parm->partRho);
        pp.get("part_dia", PeleLM::prob_parm->partDia);
        pp.get("part_temp", PeleLM::prob_parm->partTemp);
    }
}
