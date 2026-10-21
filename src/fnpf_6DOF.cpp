/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_6DOF.h"
#include"6DOF_obj_fnpf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"solver.h"
#include"fnpf_laplace.h"
#include"fnpf_fsf.h"
#include"fnpf_bed_update.h"
#include"fnpf_fsf_update.h"
#include"fnpf_amr.h"
#include"reefamr.h"
#include"solver_print.h"
#include<functional>
#include<mpi.h>

// Laplace solver seen by the time stepping when bodies are present: body geometry on the
// new sigma grid before the phi solve, body-band extrapolation after it.
class fnpf_laplace_6DOF : public fnpf_laplace
{
public:
    fnpf_laplace_6DOF(fnpf_6DOF *fb, fnpf_laplace *plap) : fb(fb), plap(plap) {}
    
    void start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, double *f, slice &Fifsf) override
    {
        fb->pre_solve(p,c,pgc);
        
        plap->start(p,c,pgc,psolv,pf,f,Fifsf);
        
        fb->post_solve(p,c,pgc,f);
    }
    
private:
    fnpf_6DOF *fb;
    fnpf_laplace *plap;
};

// X 30: the psi solves of the loads (the hull pressure from psi_0, the added mass from the unit modes)
// to X 30 x N 44 instead of the tolerance of the phi solves; the solvers (REEFMG, the composite solve
// of fnpf_amr) read N 44 at every solve, so it is set for the solve and restored after it (the
// stagnation stop of the composite solve scales with it)
fnpf_psi_tol::fnpf_psi_tol(lexer *p) : p(p), n44(p->N44)
{
    p->N44 = n44*p->X30;
}

fnpf_psi_tol::~fnpf_psi_tol()
{
    p->N44 = n44;
}

// X 35: the initial guess of a psi solve.  A psi solve started from the last solution of its kind
// (psi_0 of the previous stage, a unit mode of the previous added mass).  With X 35 1 it starts
// from the linear extrapolation in time of the last two solutions of the same kind: psi_0 of the
// same stage of the previous two steps (the stages of the SSP schemes are not equally spaced in
// time, the same stage of consecutive steps is), the unit mode of the previous two added masses,
// with subcycling psi_0 of the same stage of the previous two finest steps.  The lagged sigma
// cross-derivatives of the right-hand side (fnpf_laplace_cds2) take the guess as well, so the
// results change within the solver tolerance (as a tighter N 44 would).  The extrapolation
// factor is limited to 2 (a time step that grows).  The Dirichlet rows, the bed ghost nodes and
// the body band are set after the guess by the caller.
double fnpf_6DOF::psi_time(lexer *p, int slot) const
{
    // the time the solution of the slot stands for: the stage time of psi_0 (RK3: 0, 1, 1/2;
    // RK4: 0, 1/2, 1/2, 1), the step start for the unit modes
    const int s = (slot>=10) ? slot-10 : (slot<4 ? slot : -1);
    double c = 0.0;
    if(s>=0)
    {
        if(p->A310==3)
        c = (s==1) ? 1.0 : (s==2 ? 0.5 : 0.0);
        else
        c = (s==3) ? 1.0 : (s>=1 ? 0.5 : 0.0);
    }
    return p->simtime + c*p->dt;
}

void fnpf_6DOF::psi_guess(fnpf_6DOF_grid &G, int slot, double *f, double t)
{
    if(G.p->X35!=1 || slot>=(int)G.h1.size() || G.h1[slot]==nullptr)
    return;

    const double *p1 = G.h1[slot];
    const double *p2 = G.h2[slot];
    const int size = G.p->imax*G.p->jmax*(G.p->kmax+2);

    if(p2==nullptr || !(G.ht1[slot]>G.ht2[slot]))
    {
        for(int n=0; n<size; ++n)
        f[n] = p1[n];
        return;
    }

    const double a = MIN((t-G.ht1[slot])/(G.ht1[slot]-G.ht2[slot]), 2.0);

    for(int n=0; n<size; ++n)
    f[n] = p1[n] + a*(p1[n]-p2[n]);
}

void fnpf_6DOF::psi_record(fnpf_6DOF_grid &G, int slot, const double *f, double t)
{
    if(G.p->X35!=1)
    return;

    if(slot>=(int)G.h1.size())
    {
        G.h1.resize(slot+1,nullptr);
        G.h2.resize(slot+1,nullptr);
        G.ht1.resize(slot+1,0.0);
        G.ht2.resize(slot+1,0.0);
    }

    const int size = G.p->imax*G.p->jmax*(G.p->kmax+2);

    // the older one takes the new solution (no copy of the newer one)
    double *q = G.h2[slot];
    if(q==nullptr)
    q = new double[size];

    for(int n=0; n<size; ++n)
    q[n] = f[n];

    G.h2[slot] = G.h1[slot];
    G.ht2[slot] = G.ht1[slot];
    G.h1[slot] = q;
    G.ht1[slot] = t;
}

fnpf_6DOF::fnpf_6DOF(lexer *p, fdm_fnpf *c, ghostcell *pgc) : initialized(false), foot(p), psiD(p), zeroslice(p),
                                                               eta_ext(p), fi_ext(p), amr(nullptr), amr_layout(-1)
{
    if(p->mpirank==0)
    cout<<"6DOF FNPF startup ..."<<endl;
    
    // X 20 bodies (geometry: floating.stl, floating-1.stl, ...), coupled through the added mass
    nbody = MAX(1, p->X20);
    
    for(int nb=0; nb<nbody; ++nb)
    fb_obj.push_back(new sixdof_obj_fnpf(p,pgc,nb));
    
    Across.assign(nbody*nbody, Eigen::Matrix<double,6,6>::Zero());

    // hull pressure with the body-following time derivative of phi (exact rigid-body Neumann
    // data, no m-terms neglected): needs the velocity field of the one body on the free surface;
    // with several bodies psi_0 = phi_t with the m-terms neglected
    sixdof_obj_fnpf::chi_on = (nbody==1);
    
    if(p->mpirank==0 && nbody>1)
    cout<<"6DOF FNPF: "<<nbody<<" bodies, coupled added mass, psi_0 = phi_t (m-terms neglected)"<<endl;

    if(p->mpirank==0 && p->X28>1)
    {
    cout<<"6DOF FNPF: the added mass every "<<p->X28<<" time steps (X 28)";
    if(p->X29_d>0.0 || p->X29_a>0.0)
    cout<<", earlier after a motion by more than "<<p->X29_d<<" m or "<<p->X29_a<<" deg (X 29)";
    cout<<endl;
    }

    if(p->mpirank==0 && p->X30!=1.0)
    cout<<"6DOF FNPF: the psi solves of the loads to "<<p->X30<<" x N 44 = "<<p->X30*p->N44<<" (X 30)"<<endl;

    if(p->mpirank==0 && p->X35==1)
    cout<<"6DOF FNPF: the psi solves start from the extrapolation of the last two of their kind (X 35)"<<endl;

    gcval = (p->j_dir==0) ? 150 : 250;
    
    const int size = p->imax*p->jmax*(p->kmax+2);
    
    p->Darray(psi0,size);
    p->Darray(zero,size);
    
    for(int m=0; m<6; ++m)
    {
    psi[m] = nullptr;
    
    // 2D: sway, roll and yaw do not exist
    if(p->j_dir==0 && (m==1 || m==3 || m==5))
    continue;
    
    p->Darray(psi[m],size);
    }
    
    mark = new int[size]();
    
    pbed = new fnpf_bed_update(p);
    pvel = new fnpf_fsf_update(p,c,pgc);
    plap = nullptr;
    
    // level 0
    g0.p = p;
    g0.c = c;
    g0.id = -1;
    g0.l0 = true;
    g0.psi0 = psi0;
    for(int m=0; m<6; ++m)
    g0.psi[m] = psi[m];
    g0.zero = zero;
    g0.mark = mark;
    g0.foot = &foot;
    g0.psiD = &psiD;
    g0.zeroslice = &zeroslice;
    g0.eta_ext = &eta_ext;
    g0.fi_ext = &fi_ext;
    g0.pbed = pbed;
    g0.pvel = pvel;
}

fnpf_laplace* fnpf_6DOF::laplace(fnpf_laplace *pl)
{
    // keep the plain solver for the psi solves, hand the decorated one to the RK scheme
    plap = pl;
    
    return new fnpf_laplace_6DOF(this,pl);
}

void fnpf_6DOF::initialize(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(!initialized)
    ini(p,c,pgc);
}

void fnpf_6DOF::stage(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, slice &Keta, slice &Kfi, int iter)
{
    // loads from the current state (Fi, eta, sigma grid of the previous stage) and its
    // tendencies, then the body RK stage; the geometry follows with the next phi solve
    if(!initialized)
    ini(p,c,pgc);
    
    // external dynamics: connect once at the first step (the external solver is set up after the
    // fluid initialisation)
    if(!ext_checked)
    ext_connect(p,pgc);
    
    pvel->velcalc_sig(p,c,pgc,c->Fi);
    pgc->gcparax7(p,c->U,7);
    pgc->gcparax7(p,c->V,7);
    pgc->gcparax7(p,c->W,7);
    body_velocities(g0,pgc);
    
    // subcycling (G 7 1): level 0 steps first with a predicted copy of the body; the loads and
    // the added mass of all grids at the start of the step (all grids at the same time), the
    // finest level advances the body later (fnpf_6DOF_sub.cpp)
    if(amr_sub())
    {
        g0.pf = pf;
        
        if(iter==0)
        {
            amr_grids(p,pgc);
            
            {
            reefamr_comms_off guard(pgc);
            for(auto &G : gp)
            {
            G.pvel->velcalc_sig(G.p,G.c,pgc,G.c->Fi);
            body_velocities(G,pgc);
            }
            }
            
            forces_amr(p,c,pgc,psolv,pf,Keta,Kfi,iter);
        }
        
        // external dynamics: the copy on the poses the external solver predicts
        if(pext!=nullptr)
        ext_predict_copy(p,pgc,iter,p->simtime,p->dt);
        
        else
        for(int nb=0; nb<nbody; ++nb)
        fb_obj[nb]->amr_predict_fnpf(p,pgc,iter);
        
        return;
    }
    
    if(amr_on())
    {
        g0.pf = pf;
        
        // the patches may have changed at the end of the last step (zone following the body)
        amr_grids(p,pgc);
        
        reefamr_comms_off guard(pgc);
        for(auto &G : gp)
        {
        G.pvel->velcalc_sig(G.p,G.c,pgc,G.c->Fi);
        body_velocities(G,pgc);
        }
        
        forces_amr(p,c,pgc,psolv,pf,Keta,Kfi,iter);
    }
    else
    forces(p,c,pgc,psolv,pf,Keta,Kfi,iter);
    
    motion(p,c,pgc,iter);
    
    if(iter == ((p->A310==3) ? 2 : 3))
    psi_report(p,pgc);
}

// the psi solves of the step (all stages, all grids): their number and Laplace iterations, and the
// time of the loads since the start (the maximum over the ranks: rank 0 often holds no patch), one line
void fnpf_6DOF::psi_report(lexer *p, ghostcell *pgc)
{
    if(p->count%p->P12==0)
    {
        double t[3] = {t_loads,t_psi,t_band};
        pgc->globalmax(t,3);

        if(p->mpirank==0 && psi_solves>0)
        cout<<"6DOF: "<<nbody<<" bod"<<(nbody>1 ? "ies" : "y")<<", "<<psi_solves<<" psi solves"
            <<(am_count==p->count ? " (with the added mass)" : "")
            <<", iterations "<<psi_iter<<" (max "<<psi_itmax<<"); loads "<<setprecision(4)<<t[0]
            <<" s in total (psi solves "<<t[1]<<" s, body band "<<t[2]<<" s)"<<endl;
    }

    psi_solves = 0;
    psi_iter = 0;
    psi_itmax = 0;
}

void fnpf_6DOF::pre_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(initialized)
    geometry(g0,pgc);
}

void fnpf_6DOF::post_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *f)
{
    if(!initialized)
    return;
    
    pgc->start7V(p,f,c->bc,gcval);
    extrapolate(g0,pgc,f);
}

void fnpf_6DOF::surface(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf, int gcval_eta, int gcval_fifsf)
{
    if(initialized)
    footprint(g0,pgc,eta,Fifsf,gcval_eta,gcval_fifsf);
}

void fnpf_6DOF::velocity(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(initialized)
    body_velocities(g0,pgc);
}

void fnpf_6DOF::body_velocities(fnpf_6DOF_grid &G, ghostcell *pgc)
{
    // Velocities at the body nodes, layer by layer from the fluid inwards: the average of the
    // fluid (or already filled) neighbours.  velcalc_sig differentiates the body-band values
    // of Fi there with 4th-order stencils that reach two nodes into the body, where the
    // extrapolation is only a continuation; at the free-surface node of a footprint column
    // W is Fz from the extended Fifsf over the body interior.  Both gave |U|,|W| of several
    // m/s that alternate from column to column (moored box, X 311, from about 7 s on: 4-5
    // m/s with a body velocity of 0.2-0.5 m/s).  The hull loads sample U,V,W half a cell off
    // the hull with trilinear weights that include these nodes, and umax/wmax (time step,
    // N 61 stop, output) took them up.  The fluid velocity continued into the body keeps the
    // slip velocity at the hull and the normal velocity of the fluid next to it.
    // On a patch (mesh refinement) the fill is local to the patch.
    lexer *p = G.p;
    fdm_fnpf *c = G.c;
    int *mark = G.mark;
    double *U = c->U, *V = c->V, *W = c->W;
    const double *FBF = c->FBF;
    const int size = p->imax*p->jmax*(p->kmax+2);
    const int sI = p->jmax*p->kmaxF;
    const int sJ = p->kmaxF;
    
    for(int n=0; n<size; ++n)
    mark[n]=0;
    
    const int maxpass = 4*(p->gknox + p->gknoy + p->gknoz) + 2;
    
    for(int pass=0; pass<maxpass; ++pass)
    {
        int filled=0;
        
        ILOOP
        JLOOP
        FKLOOP
        {
            const int q = FIJK;
            
            if(p->flag7[q]>0 && FBF[q]>0.5 && mark[q]==0)
            {
                int nbr[6] = {q-sI, q+sI, q-sJ, q+sJ, q-1, q+1};
                int nn = 6;
                
                if(p->j_dir==0)
                {
                nbr[2] = q-1;
                nbr[3] = q+1;
                nn = 4;
                }
                
                double su=0.0, sv=0.0, sw=0.0;
                int cnt=0;
                
                for(int m=0; m<nn; ++m)
                {
                    const int r = nbr[m];
                    
                    if(p->flag7[r]>0 && (FBF[r]<0.5 || (mark[r]>0 && mark[r]<=pass)))
                    {
                    su += U[r];
                    sv += V[r];
                    sw += W[r];
                    ++cnt;
                    }
                }
                
                if(cnt>0)
                {
                U[q] = su/double(cnt);
                V[q] = sv/double(cnt);
                W[q] = sw/double(cnt);
                mark[q] = pass+1;
                ++filled;
                }
            }
        }
        
        if(G.l0)
        {
        pgc->gcparax7(p,U,7);
        pgc->gcparax7(p,V,7);
        pgc->gcparax7(p,W,7);
        pgc->gcparax7int(p,mark,7);
        
        filled = pgc->globalisum(filled);
        }
        
        if(filled==0)
        break;
    }
}

fnpf_6DOF::~fnpf_6DOF()
{
    for(auto &G : gp)
    free_grid(G);
}

// X 28: the added mass (the unit-mode psi solves of the free DOFs, with several bodies the cross terms)
// at the first stage of the first step and of every X 28-th step after the last one; in between the
// bodies keep the last one (the loads of psi_0 are new at every stage). X 29: also as soon as a body
// has moved by more than X29_d [m] or turned by more than X29_a [deg] since. The same decision on all
// ranks: the poses are those of the body kernel, identical on every rank
bool fnpf_6DOF::am_due(lexer *p, int iter)
{
    if(iter!=0)
    return false;

    bool due = (p->X28<=1 || am_count<0 || (int)am_c.size()!=nbody || p->count-am_count>=p->X28);

    if(!due && (p->X29_d>0.0 || p->X29_a>0.0))
    for(int nb=0; nb<nbody; ++nb)
    {
        if(p->X29_d>0.0 && (fb_obj[nb]->fnpf_cog()-am_c[nb]).norm()>p->X29_d)
        due = true;

        if(p->X29_a>0.0 && (fb_obj[nb]->fnpf_angles()-am_ang[nb]).cwiseAbs().maxCoeff()*(180.0/PI)>p->X29_a)
        due = true;
    }

    if(due)
    {
        am_count = p->count;
        am_c.resize(nbody);
        am_ang.resize(nbody);
        for(int nb=0; nb<nbody; ++nb)
        {
            am_c[nb] = fb_obj[nb]->fnpf_cog();
            am_ang[nb] = fb_obj[nb]->fnpf_angles();
        }
    }

    return due;
}

void fnpf_6DOF::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(nbody>1 && amr_sub())
    {
        if(p->mpirank==0)
        cout<<"6DOF FNPF: several bodies (X 20 > 1) with mesh refinement need G 7 0 (no subcycling) for now"<<endl;

        pgc->final();
        exit(1);
    }

    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->initialize_fnpf(p,c,pgc);
    
    geometry(g0,pgc);
    extrapolate(g0,pgc,c->Fi);
    
    initialized = true;
}

void fnpf_6DOF::forces(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf,
                       slice &Keta, slice &Kfi, int iter)
{
    const double t0 = MPI_Wtime();

    // psi_0: phi_t at fixed z on the free surface, phi_t|z = dFifsf/dt - Fz*deta/dt
    // (one body: + V.grad(phi), the body-following derivative chi, see face_data_fnpf)
    SLICELOOP4
    psiD(i,j) = Kfi(i,j) - c->Fz(i,j)*Keta(i,j);
    
    if(nbody==1)
    fb_obj[0]->chi_fsf(p,c,psiD,foot);
    
    pgc->gcsl_start4(p,psiD,50);
    
    // X 18: running mean of chi at body-fixed points removed next to the body
    if(nbody==1)
    fb_obj[0]->chi_mean(p,pgc,psiD,foot);
    
    zero_face(g0);
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->face_data_fnpf(p,c,pgc,-2,c->FBF,c->FBu,c->FBv,c->FBw);
    exchange_face(g0,pgc);
    
    solve_psi(p,c,pgc,psolv,plap,pf,psi0,psiD,psi_slot(-1,iter,false));
    
    // added mass: unit-mode solves of the free DOFs, once per time step, with X 28 every X 28 steps
    // (am_due; the rows and columns of fixed or prescribed DOFs are not used in the motion); with
    // several bodies the modes of body nb also give its cross added mass on all other bodies
    const bool refresh = am_due(p,iter);
    
    if(refresh && nbody>1)
    for(auto &A : Across)
    A.setZero();

    for(int nb=0; nb<nbody; ++nb)
    {
        const bool computeA = refresh && !fb_obj[nb]->fnpf_fixed(p);
        
        double *pm[6] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

        for(int m=0; m<6; ++m)
        if(computeA && psi[m]!=nullptr && !fb_obj[nb]->dof_fixed_fnpf(p,m))
        pm[m] = psi[m];

        for(int m=0; m<6; ++m)
        if(pm[m]!=nullptr)
        {
            zero_face(g0);
            fb_obj[nb]->face_data_fnpf(p,c,pgc,m,c->FBF,c->FBu,c->FBv,c->FBw);
            exchange_face(g0,pgc);
            
            solve_psi(p,c,pgc,psolv,plap,pf,psi[m],zeroslice,psi_slot(m,iter,false));
        }
        
        if(computeA && nbody>1)
        for(int nb2=0; nb2<nbody; ++nb2)
        if(nb2!=nb)
        fb_obj[nb2]->added_mass_block_fnpf(p,c,pgc,zero,pm,Across[nb2*nbody+nb]);

        fb_obj[nb]->forces_fnpf(p,c,pgc,psi0,pm,computeA);

        // force log: first stage = state at the current time level
        if(iter==0)
        fb_obj[nb]->print_force_fnpf(p);
    }

    t_loads += MPI_Wtime()-t0;
}

void fnpf_6DOF::motion(lexer *p, fdm_fnpf *c, ghostcell *pgc, int iter)
{
    if(pext!=nullptr)
    {
    motion_external(p,c,pgc,iter);
    return;
    }
    
    // last stage: 2 for RK3 (A 310 3), 3 for RK4 (A 310 4)
    const bool finalize = (iter == ((p->A310==3) ? 2 : 3));
    
    if(nbody>1)
    motion_coupled(p,c,pgc,iter,finalize);

    for(int nb=0; nb<nbody; ++nb)
    {
        if(nbody==1)
        fb_obj[nb]->solve_eqmotion_fnpf(p,c,pgc,iter,finalize);

        fb_obj[nb]->update_position_fnpf(p,pgc,finalize);
        
        if(finalize)
        fb_obj[nb]->print_fnpf(p,pgc,iter);
    }
}

void fnpf_6DOF::motion_coupled(lexer *p, fdm_fnpf *c, ghostcell *pgc, int iter, bool finalize)
{
    // Several bodies: each body assembles its block of the instantaneous added-mass system
    //   (M_i + A_ii) a_i + sum_j A_ij a_j = F_i
    // (own added mass, gyroscopic and PTO terms, fixed DOFs), the cross blocks A_ij come from
    // the unit modes of body j on the hull of body i, symmetrised as (A_ij + A_ji^T)/2.
    // The coupled system is solved on every rank (identical data), then each body takes its
    // accelerations and does its RK stage.
    const int n6 = 6*nbody;

    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->eqmotion_assemble_fnpf(p,c,pgc,iter,finalize);

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n6,n6);
    Eigen::VectorXd r = Eigen::VectorXd::Zero(n6);
    std::vector<int> has(nbody,0);
    std::vector<int> fixed(n6,1);

    for(int nb=0; nb<nbody; ++nb)
    {
        Eigen::Matrix<double,6,6> L;
        Eigen::Matrix<double,6,1> rb;
        bool fx[6];

        has[nb] = fb_obj[nb]->eqmotion_block_fnpf(L,rb,fx) ? 1 : 0;

        if(!has[nb])
        {
        L.setIdentity();
        rb.setZero();
        for(int n=0; n<6; ++n)
        fx[n] = true;
        }

        K.block<6,6>(6*nb,6*nb) = L;
        r.segment<6>(6*nb) = rb;

        for(int n=0; n<6; ++n)
        fixed[6*nb+n] = fx[n] ? 1 : 0;
    }

    for(int i=0; i<nbody; ++i)
    for(int j=0; j<nbody; ++j)
    if(i!=j && has[i] && has[j])
    {
        const Eigen::Matrix<double,6,6> Aij = 0.5*(Across[i*nbody+j] + Across[j*nbody+i].transpose());

        for(int a=0; a<6; ++a)
        for(int b=0; b<6; ++b)
        if(!fixed[6*i+a] && !fixed[6*j+b])
        K(6*i+a,6*j+b) = Aij(a,b);
    }

    const Eigen::VectorXd acc = K.partialPivLu().solve(r);

    for(int nb=0; nb<nbody; ++nb)
    {
        const Eigen::Matrix<double,6,1> a = acc.segment<6>(6*nb);
        fb_obj[nb]->eqmotion_advance_fnpf(p,pgc,iter, has[nb] ? &a : nullptr);
    }
}

void fnpf_6DOF::footprint(fnpf_6DOF_grid &G, ghostcell *pgc, slice &eta, slice &Fifsf, int gcval_eta, int gcval_fifsf)
{
    // Harmonic extension of eta and Fifsf over the footprint, the surrounding free surface
    // acts as Dirichlet data.
    // The stage value in the footprint contains the RK tendency of columns whose free-surface
    // node lies inside the body (Fz there is meaningless and can be large), so it is discarded:
    // the relaxation restarts from the last extension (eta_ext, fi_ext) and is iterated to
    // convergence. Columns entering the footprint start from their last free-surface value.
    // On a patch (mesh refinement) the iteration is local to the patch.
    lexer *p = G.p;
    slice4 &foot = *G.foot;
    slice4 &eta_ext = *G.eta_ext;
    slice4 &fi_ext = *G.fi_ext;
    
    if(!G.ext_ini)
    {
        SLICELOOP4
        {
        eta_ext(i,j) = eta(i,j);
        fi_ext(i,j)  = Fifsf(i,j);
        }
        
        pgc->gcsl_start4(p,eta_ext,gcval_eta);
        pgc->gcsl_start4(p,fi_ext,gcval_fifsf);
        
        G.ext_ini = true;
    }
    
    if(G.footcount>0)
    {
        SLICELOOP4
        if(foot(i,j)>0.5)
        {
        eta(i,j)   = eta_ext(i,j);
        Fifsf(i,j) = fi_ext(i,j);
        }
        
        pgc->gcsl_start4(p,eta,gcval_eta);
        pgc->gcsl_start4(p,Fifsf,gcval_fifsf);
        
        const double tol = 1.0e-9*MAX(p->DXM,1.0e-12);
        
        for(int it=0; it<1000; ++it)
        {
            double dmax=0.0;
            
            SLICELOOP4
            if(foot(i,j)>0.5)
            {
                double se = eta(i-1,j) + eta(i+1,j);
                double sf = Fifsf(i-1,j) + Fifsf(i+1,j);
                double nn = 2.0;
                
                if(p->j_dir==1)
                {
                se += eta(i,j-1) + eta(i,j+1);
                sf += Fifsf(i,j-1) + Fifsf(i,j+1);
                nn += 2.0;
                }
                
                // over-relaxed Gauss-Seidel
                const double de = 1.6*(se/nn - eta(i,j));
                const double df = 1.6*(sf/nn - Fifsf(i,j));
                
                eta(i,j)   += de;
                Fifsf(i,j) += df;
                
                dmax = MAX(dmax,fabs(de));
            }
            
            pgc->gcsl_start4(p,eta,gcval_eta);
            pgc->gcsl_start4(p,Fifsf,gcval_fifsf);
            
            if(G.l0)
            dmax = pgc->globalmax(dmax);
            
            if(dmax<tol)
            break;
        }
        
        // X 17 w: relaxation of eta in the free-surface columns next to the footprint,
        //     eta = (1 - w) eta + w (mean of the 4 neighbours),
        // once per stage. A hull moving at speed uncovers columns at its waterline (stern
        // quarters) every step; there a 2-cell disturbance of eta grows until the emergency stop
        // (smooth hull, Fr 0.14, from about 10 s on). The relaxation of eta alone stabilises it
        // without changing the hull loads noticeably; Fifsf is left alone (relaxing the
        // potential next to the hull changes phi_t and the hull pressure).
        // On by default (0.5): on fine grids (0.05 m patch or uniform grid at a 1 m float, LUPA)
        // the columns at the hull also blow up in the first deep wave trough, without forward
        // speed, and a smaller time step does not help. X 17 0 switches it off.
        if(p->X17>0.0)
        {
            const double w = MIN(p->X17,1.0);
            
            SLICELOOP4
            eta_ext(i,j) = eta(i,j);
            
            pgc->gcsl_start4(p,eta_ext,gcval_eta);
            
            SLICELOOP4
            if(foot(i,j)<0.5 && (foot(i-1,j)>0.5 || foot(i+1,j)>0.5 || (p->j_dir==1 && (foot(i,j-1)>0.5 || foot(i,j+1)>0.5))))
            {
                double se = eta_ext(i-1,j) + eta_ext(i+1,j);
                double nn = 2.0;
                
                if(p->j_dir==1)
                {
                se += eta_ext(i,j-1) + eta_ext(i,j+1);
                nn += 2.0;
                }
                
                eta(i,j) = (1.0 - w)*eta_ext(i,j) + w*se/nn;
            }
            
            pgc->gcsl_start4(p,eta,gcval_eta);
        }
    }
    
    // remember the clean state: extension in the footprint, free surface elsewhere
    SLICELOOP4
    {
    eta_ext(i,j) = eta(i,j);
    fi_ext(i,j)  = Fifsf(i,j);
    }
    
    pgc->gcsl_start4(p,eta_ext,gcval_eta);
    pgc->gcsl_start4(p,fi_ext,gcval_fifsf);
}

void fnpf_6DOF::geometry(fnpf_6DOF_grid &G, ghostcell *pgc)
{
    lexer *p = G.p;
    fdm_fnpf *c = G.c;
    slice4 &foot = *G.foot;
    
    // fnpf_sigma_update sets ZSN in the interior columns only; the sampling of the hull
    // loads (ccipol7V) reaches into the ghost columns at subdomain boundaries
    if(G.l0)
    pgc->gcparax7(p,p->ZSN,7);
    
    const int size = p->imax*p->jmax*(p->kmax+2);
    
    for(int n=0; n<size; ++n)
    c->FBF[n] = 0.0;
    
    SLICELOOP4
    foot(i,j) = 0.0;
    
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->ray_cast_fnpf(p,c,pgc,c->FBF,foot);
    
    if(G.l0)
    pgc->gcparax7(p,c->FBF,7);
    pgc->gcsl_start4(p,foot,50);
    
    int count=0;
    SLICELOOP4
    if(foot(i,j)>0.5)
    ++count;
    
    G.footcount = G.l0 ? pgc->globalisum(count) : count;

        // Neumann data for phi: rigid-body velocity
    zero_face(G);
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->face_data_fnpf(p,c,pgc,-1,c->FBF,c->FBu,c->FBv,c->FBw);
    exchange_face(G,pgc);
}

void fnpf_6DOF::extrapolate(fnpf_6DOF_grid &G, ghostcell *pgc, double *f)
{
    // All body nodes, layer by layer from the fluid inwards. A body node takes the
    // average over its fluid (or already filled) neighbours r of the value continued
    // with the Neumann data of the face between them,
    //     f_q = f_r + V_e(q)*(x_q - x_r),   V_e = (FBu,FBv,FBw).e,
    // i.e. the same ghost value the Laplace assembly eliminates on that face
    // (fnpf_body_bc), and deeper inside a linear continuation with the body velocity.
    // The body nodes act as consistent ghost nodes for every stencil that reaches into
    // the body: the free-surface Fz of columns over a submerged part of the hull (tilted
    // walls, corners of a body that is not aligned with the grid), velcalc_sig next to
    // the hull, the lagged cross terms of the Laplace rhs and the sampling of the hull
    // loads. A plain average gave a zero normal gradient there instead of V_e: wrong Fz
    // and kinematic FSBC in those columns, a local dip or spike of eta at the hull that
    // grew until the emergency stop.
    // Face data: the mode of the solve that produced f (geometry() sets -1 before the
    // phi solve, solve_psi the mode of each psi solve).
    // The whole interior has to be filled, not only two layers: the body nodes are
    // identity rows, so deeper nodes would keep the value from the time they were
    // covered. phi drifts in time (Bernoulli constant, set-down), the stale interior
    // then gives large spurious gradients, which velcalc_sig turns into body-band and
    // footprint velocities (|W| of O(10) m/s after a few wave periods). These leak into
    // the hull loads through the trilinear sampling and blow the case up.
    // On a patch (mesh refinement) the fill is local to the patch.
    lexer *p = G.p;
    fdm_fnpf *c = G.c;
    int *mark = G.mark;
    const int size = p->imax*p->jmax*(p->kmax+2);
    const int sI = p->jmax*p->kmaxF;
    const int sJ = p->kmaxF;
    const double *FBF = c->FBF;
    const double *FBu = c->FBu;
    const double *FBv = c->FBv;
    const double *FBw = c->FBw;
    const double *ZSN = p->ZSN;
    
    for(int n=0; n<size; ++n)
    mark[n]=0;
    
    // fill until no body node is left that touches the fluid or a filled node;
    // nodes without any connection to the fluid keep their value
    const int maxpass = 4*(p->gknox + p->gknoy + p->gknoz) + 2;
    
    for(int pass=0; pass<maxpass; ++pass)
    {
        int filled=0;
        
        ILOOP
        JLOOP
        FKLOOP
        {
            const int q = FIJK;
            
            if(p->flag7[q]>0 && FBF[q]>0.5 && mark[q]==0)
            {
                // neighbour and continuation f_q - f_r = V_e*(x_q - x_r) towards it
                int nbr[6] = {q-sI, q+sI, q-sJ, q+sJ, q-1, q+1};
                double dv[6] = { FBu[q]*(p->XP[IP]-p->XP[IM1]),
                                -FBu[q]*(p->XP[IP1]-p->XP[IP]),
                                 FBv[q]*(p->YP[JP]-p->YP[JM1]),
                                -FBv[q]*(p->YP[JP1]-p->YP[JP]),
                                 FBw[q]*(ZSN[q]-ZSN[q-1]),
                                -FBw[q]*(ZSN[q+1]-ZSN[q])};
                const int nn = (p->j_dir==1) ? 6 : 4;
                
                if(p->j_dir==0)
                {
                nbr[2] = q-1;
                nbr[3] = q+1;
                dv[2] = dv[4];
                dv[3] = dv[5];
                }
                
                double sum=0.0;
                int cnt=0;
                
                for(int m=0; m<nn; ++m)
                {
                    const int r = nbr[m];
                    
                    if(p->flag7[r]>0 && (FBF[r]<0.5 || (mark[r]>0 && mark[r]<=pass)))
                    {
                    sum += f[r] + dv[m];
                    ++cnt;
                    }
                }
                
                if(cnt>0)
                {
                f[q] = sum/double(cnt);
                mark[q] = pass+1;
                ++filled;
                }
            }
        }
        
        // values and fill state of the ghost nodes, so the result does not depend on the
        // domain decomposition
        if(G.l0)
        {
        pgc->gcparax7(p,f,7);
        pgc->gcparax7int(p,mark,7);
        
        filled = pgc->globalisum(filled);
        }
        
        if(filled==0)
        break;
    }
}

void fnpf_6DOF::solve_psi(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_laplace *plap, fnpf_fsf *pf,
                          double *f, slice &D, int slot)
{
    const double t0 = MPI_Wtime();
    const double ts = psi_time(p,slot);
    
    psi_guess(g0,slot,f,ts);

    // Dirichlet on the free surface outside the footprint
    FFILOOP4
    if(foot(i,j)<0.5)
    {
        f[FIJK]   = D(i,j);
        f[FIJKp1] = D(i,j);
        f[FIJKp2] = D(i,j);  
        f[FIJKp3] = D(i,j);
    }
    
    pbed->bedbc_sig(p,c,pgc,f,pf);
    pgc->start7V(p,f,c->bc,gcval);
    
    // the A329 inflow flux belongs to phi only
    double *uin = c->Uin;
    c->Uin = zero;
    
    solver_print_quiet() = true;
    {
    fnpf_psi_tol tol(p);
    plap->start(p,c,pgc,psolv,pf,f,D);
    }
    solver_print_quiet() = false;
    psi_count(p->poissoniter);
    
    c->Uin = uin;

    pgc->start7V(p,f,c->bc,gcval);

    const double t1 = MPI_Wtime();
    extrapolate(g0,pgc,f);
    t_band += MPI_Wtime()-t1;
    
    psi_record(g0,slot,f,ts);
    t_psi += MPI_Wtime()-t0;
}

void fnpf_6DOF::zero_face(fnpf_6DOF_grid &G)
{
    lexer *p = G.p;
    fdm_fnpf *c = G.c;
    const int size = p->imax*p->jmax*(p->kmax+2);
    
    for(int n=0; n<size; ++n)
    {
    c->FBu[n] = 0.0;
    c->FBv[n] = 0.0;
    c->FBw[n] = 0.0;
    }
}

void fnpf_6DOF::exchange_face(fnpf_6DOF_grid &G, ghostcell *pgc)
{
    if(!G.l0)
    return;
    
    lexer *p = G.p;
    fdm_fnpf *c = G.c;
    pgc->gcparax7(p,c->FBu,7);
    pgc->gcparax7(p,c->FBv,7);
    pgc->gcparax7(p,c->FBw,7);
}

// --------------------------------------------------------------------- mesh refinement
//  With fnpf_amr every grid of the hierarchy carries the body at its own resolution: ray
//  cast, footprint extension, body-band extrapolation and the Neumann data of the body faces
//  are done per grid (local on the patches), the phi and psi solves are composite solves
//  over all grids (fnpf_amr::lap_solve, lap_solve_psi), and every hull triangle is integrated
//  once, on the finest grid that holds its centroid.

bool fnpf_6DOF::amr_on() const
{
    return amr!=nullptr && amr->active();
}

bool fnpf_6DOF::amr_sub() const
{
    return amr!=nullptr && amr->sub_on();
}

void fnpf_6DOF::amr_attach(fnpf_amr *a)
{
    amr = a;
}

void fnpf_6DOF::amr_bodies(vector<sixdof_obj*> &obj)
{
    for(auto o : fb_obj)
    obj.push_back(o);
}

void fnpf_6DOF::free_grid(fnpf_6DOF_grid &G)
{
    // the patch lexer may already be gone (regrid): plain deletes only
    delete [] G.psi0;
    delete [] G.zero;
    for(int m=0; m<6; ++m)
    delete [] G.psi[m];
    delete [] G.mark;
    for(auto q : G.h1) delete [] q;
    for(auto q : G.h2) delete [] q;
    G.h1.clear(); G.h2.clear();
    delete G.foot; delete G.psiD; delete G.zeroslice; delete G.eta_ext; delete G.fi_ext;
    delete G.pbed;
    delete G.pvel;
}

// body grids of the patches, after the (re)gridding of fnpf_amr: grids of patches that are
// still there are kept (with their footprint state), the others are freed, new patches get a
// grid with the body geometry of the current position
void fnpf_6DOF::amr_grids(lexer *p, ghostcell *pgc)
{
    if(amr==nullptr || amr_layout==amr->layout())
    return;
    
    vector<fnpf_6DOF_grid> old;
    old.swap(gp);
    vector<char> used(old.size(),0);
    
    reefamr_comms_off guard(pgc);
    
    for(int n=0; n<amr->patches(); ++n)
    {
        const int serial = amr->patch_serial(n);
        
        bool found=false;
        for(size_t k=0; k<old.size(); ++k)
        if(!used[k] && old[k].serial==serial)
        {
            used[k] = 1;
            old[k].id = n;
            old[k].Keta = &amr->patch_tendency(n,0);
            old[k].Kfi = &amr->patch_tendency(n,1);
            gp.push_back(old[k]);
            found = true;
            break;
        }
        if(found)
        continue;
        
        fnpf_6DOF_grid G;
        lexer *pp = amr->patch_lexer(n);
        fdm_fnpf *cc = amr->patch_fdm(n);
        const int size = pp->imax*pp->jmax*(pp->kmax+2);
        
        G.p = pp;
        G.c = cc;
        G.pf = amr->patch_fsf(n);
        G.id = n;
        G.serial = serial;
        G.l0 = false;
        G.del = 0.5*fb_obj[0]->fnpf_dsm()/double(1<<amr->patch_level(n));
        
        pp->Darray(G.psi0,size);
        pp->Darray(G.zero,size);
        for(int m=0; m<6; ++m)
        if(psi[m]!=nullptr)
        pp->Darray(G.psi[m],size);
        G.mark = new int[size]();
        G.foot = new slice4(pp);
        G.psiD = new slice4(pp);
        G.zeroslice = new slice4(pp);
        G.eta_ext = new slice4(pp);
        G.fi_ext = new slice4(pp);
        G.pbed = new fnpf_bed_update(pp);
        G.pvel = new fnpf_fsf_update(pp,cc,pgc);
        G.Keta = &amr->patch_tendency(n,0);
        G.Kfi = &amr->patch_tendency(n,1);
        G.fresh = true;
        
        gp.push_back(G);
    }
    
    for(size_t k=0; k<old.size(); ++k)
    if(!used[k])
    free_grid(old[k]);
    
    amr_layout = amr->layout();
    
    // restart: the grids of the patches of the checkpoint take its arrays, the others are set up
    // after the bodies (fnpf_6DOF_ckpt.cpp)
    if(!ckpt_grids)
    amr_fresh(pgc);
}

// psi0 and the unit modes of a fresh grid: prolonged from the coarser grid as the initial
// guess of the next psi solves (gp is in patch order, coarse levels first).  Left at zero
// they started the composite solves with a jump at the patch edge
// (all ranks together, level by level: the parent may be on another rank, G 40)
void fnpf_6DOF::amr_fresh(ghostcell *pgc)
{
    vector<double*> fg(gp.size()+1);
    
    reefamr_comms_off guard(pgc);
    
    if(initialized)
    {
        vector<char> need(gp.size(),0);
        for(auto &G : gp)
        if(G.fresh)
        {
            geometry(G,pgc);
            extrapolate(G,pgc,G.c->Fi);
            amr->patch_walls_fi(G.id,G.c->Fi);
            need[G.id] = 1;
        }
        
        for(int m=-1; m<6; ++m)
        {
            fg[0] = (m<0) ? g0.psi0 : g0.psi[m];
            vector<char> nd = need;
            for(size_t k=0; k<gp.size(); ++k)
            {
                fg[k+1] = (m<0) ? gp[k].psi0 : gp[k].psi[m];
                if(fg[k+1]==nullptr)
                nd[k] = 0;
            }
            if(fg[0]!=nullptr)
            amr->prolong_cols(&fg[0],nd);
        }
        
        for(auto &G : gp)
        G.fresh = false;
    }
}

// before the composite phi solve of a stage: the body on the new sigma grid of every grid
void fnpf_6DOF::amr_geometry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(!initialized)
    return;
    
    amr_grids(p,pgc);
    
    geometry(g0,pgc);
    
    reefamr_comms_off guard(pgc);
    for(auto &G : gp)
    geometry(G,pgc);
}

// after it: ghost nodes and body-band extrapolation on every grid
void fnpf_6DOF::amr_post_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *f)
{
    if(!initialized)
    return;
    
    pgc->start7V(p,f,c->bc,gcval);
    extrapolate(g0,pgc,f);
    
    reefamr_comms_off guard(pgc);
    for(auto &G : gp)
    {
        amr->patch_walls_fi(G.id,G.c->Fi);
        extrapolate(G,pgc,G.c->Fi);
    }
}

// patch stage values: footprint extension on the patch
void fnpf_6DOF::amr_surface(lexer *p, ghostcell *pgc, int id, slice &eta, slice &Fifsf)
{
    if(!initialized)
    return;
    
    for(auto &G : gp)
    if(G.id==id)
    {
        reefamr_comms_off guard(pgc);
        const int ge = (G.p->j_dir==0) ? 155 : 55;
        const int gf = (G.p->j_dir==0) ? 160 : 60;
        footprint(G,pgc,eta,Fifsf,ge,gf);
    }
}

// psi solve on all grids: m -1 psi_0 with the free-surface data psiD, m >= 0 the unit mode m
// with homogeneous free-surface data
void fnpf_6DOF::solve_psi_amr(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, int m, bool withD, int slot)
{
    const double t0 = MPI_Wtime();
    const double ts = psi_time(p,slot);
    const int ng = 1+(int)gp.size();
    vector<double*> fs(ng);
    vector<slice*> Ds(ng);
    
    auto prep = [&](fnpf_6DOF_grid &G, double *f, slice &D)
    {
        lexer *p = G.p;
        slice4 &foot = *G.foot;
        
        // Dirichlet on the free surface outside the footprint
        FFILOOP4
        if(foot(i,j)<0.5)
        {
            f[FIJK]   = D(i,j);
            f[FIJKp1] = D(i,j);
            f[FIJKp2] = D(i,j);  
            f[FIJKp3] = D(i,j);
        }
        
        G.pbed->bedbc_sig(p,G.c,pgc,f,G.l0 ? pf : G.pf);
    };
    
    for(int g=0; g<ng; ++g)
    {
        fnpf_6DOF_grid &G = (g==0) ? g0 : gp[g-1];
        fs[g] = (m<0) ? G.psi0 : G.psi[m];
        Ds[g] = withD ? static_cast<slice*>(G.psiD) : static_cast<slice*>(G.zeroslice);
        psi_guess(G,slot,fs[g],ts);
    }
    
    prep(g0,fs[0],*Ds[0]);
    pgc->start7V(p,fs[0],c->bc,gcval);
    
    {
    reefamr_comms_off guard(pgc);
    for(size_t n=0; n<gp.size(); ++n)
    {
        prep(gp[n],fs[n+1],*Ds[n+1]);
        amr->patch_walls_fi(gp[n].id,fs[n+1]);
    }
    }
    
    // the A329 inflow flux belongs to phi only
    double *uin = c->Uin;
    c->Uin = zero;
    
    {
    fnpf_psi_tol tol(p);
    amr->lap_solve_psi(p,c,pgc,psolv,pf,&fs[0],&Ds[0]);
    }
    psi_count(p->solveriter);
    
    c->Uin = uin;

    pgc->start7V(p,fs[0],c->bc,gcval);

    const double t1 = MPI_Wtime();
    extrapolate(g0,pgc,fs[0]);

    {
    reefamr_comms_off guard(pgc);
    for(size_t n=0; n<gp.size(); ++n)
    {
        amr->patch_walls_fi(gp[n].id,fs[n+1]);
        extrapolate(gp[n],pgc,fs[n+1]);
    }
    }
    t_band += MPI_Wtime()-t1;
    
    for(int g=0; g<ng; ++g)
    psi_record((g==0) ? g0 : gp[g-1],slot,fs[g],ts);
    
    t_psi += MPI_Wtime()-t0;
}

void fnpf_6DOF::forces_amr(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf,
                           slice &Keta, slice &Kfi, int iter)
{
    const double t0 = MPI_Wtime();

    amr_grids(p,pgc);
    
    // psi_0: phi_t at fixed z on the free surface, phi_t|z = dFifsf/dt - Fz*deta/dt
    SLICELOOP4
    psiD(i,j) = Kfi(i,j) - c->Fz(i,j)*Keta(i,j);
    
    if(nbody==1)
    fb_obj[0]->chi_fsf(p,c,psiD,foot);
    
    pgc->gcsl_start4(p,psiD,50);
    
    // X 18: running mean of chi at body-fixed points removed next to the body
    if(nbody==1)
    fb_obj[0]->chi_mean(p,pgc,psiD,foot);
    
    {
    reefamr_comms_off guard(pgc);
    for(auto &G : gp)
    {
        lexer *p = G.p;
        slice4 &D = *G.psiD;
        slice &Ke = *G.Keta, &Kf = *G.Kfi;
        
        SLICELOOP4
        D(i,j) = Kf(i,j) - G.c->Fz(i,j)*Ke(i,j);
        
        if(nbody==1)
        fb_obj[0]->chi_fsf(p,G.c,D,*G.foot);
        
        pgc->gcsl_start4(p,D,50);
        
        // X 18 on the patches
        if(nbody==1)
        fb_obj[0]->chi_mean_patch(p,pgc,D,*G.foot,amr->patch_level(G.id));
    }
    }
    
    // Neumann data on all grids: mode -2 (psi_0) for every body, a unit mode (0-5) only for the
    // body sel, the others held (its added mass and its cross added mass on the others)
    auto face = [&](int mode, int sel)
    {
        zero_face(g0);
        for(int nb=0; nb<nbody; ++nb)
        if(sel<0 || nb==sel)
        fb_obj[nb]->face_data_fnpf(p,c,pgc,mode,c->FBF,c->FBu,c->FBv,c->FBw);
        exchange_face(g0,pgc);
        
        reefamr_comms_off guard(pgc);
        for(auto &G : gp)
        {
            zero_face(G);
            for(int nb=0; nb<nbody; ++nb)
            if(sel<0 || nb==sel)
            fb_obj[nb]->face_data_fnpf(G.p,G.c,pgc,mode,G.c->FBF,G.c->FBu,G.c->FBv,G.c->FBw);
        }
    };
    
    face(-2,-1);
    solve_psi_amr(p,c,pgc,psolv,pf,-1,true,psi_slot(-1,iter,false));
    
    // added mass: unit-mode solves of the free DOFs, once per time step, with X 28 every X 28 steps
    // (am_due; rows and columns of fixed or prescribed DOFs are not used in the motion); with several
    // bodies the modes of body nb also give its cross added mass on all other bodies (one patch per
    // body, G 12)
    const bool refresh = am_due(p,iter);
    
    if(refresh && nbody>1)
    for(auto &A : Across)
    A.setZero();

    // psi arrays of grid G restricted to the solved modes
    auto modes = [&](fnpf_6DOF_grid &G, const bool *on, double **out)
    {
        for(int m=0; m<6; ++m)
        out[m] = on[m] ? G.psi[m] : nullptr;
    };

    // hull integral of body nbi over all grids, every triangle on the finest grid that holds its centroid
    auto hull_sum = [&](int nbi, const bool *on, bool computeA, sixdof_obj_fnpf::fnpf_force_sum &S, bool zero_psi0)
    {
        fb_obj[nbi]->forces_fnpf_zero(p,S);
        
        for(int g=-1; g<(int)gp.size(); ++g)
        {
            fnpf_6DOF_grid &G = (g<0) ? g0 : gp[g];
            const int id = G.id;
            std::function<bool(double,double)> own = [&](double x, double y)
            {
                return amr->owns_point(x,y,id);
            };
            const double del = (g<0) ? 0.5*fb_obj[nbi]->fnpf_dsm() : G.del;
            double *pm[6];
            modes(G,on,pm);
            fb_obj[nbi]->forces_fnpf_sum(G.p,G.c,zero_psi0 ? G.zero : G.psi0,pm,computeA,del,&own,S);
        }
    };

    for(int nb=0; nb<nbody; ++nb)
    {
        const bool computeA = refresh && !fb_obj[nb]->fnpf_fixed(p);

        bool on[6];
        for(int m=0; m<6; ++m)
        on[m] = computeA && psi[m]!=nullptr && !fb_obj[nb]->dof_fixed_fnpf(p,m);

        for(int m=0; m<6; ++m)
        if(on[m])
        {
            face(m,nb);
            solve_psi_amr(p,c,pgc,psolv,pf,m,false,psi_slot(m,iter,false));
        }

        if(computeA && nbody>1)
        for(int nb2=0; nb2<nbody; ++nb2)
        if(nb2!=nb)
        {
            sixdof_obj_fnpf::fnpf_force_sum S2;
            hull_sum(nb2,on,true,S2,true);
            fb_obj[nb2]->added_mass_from_sum(pgc,S2,Across[nb2*nbody+nb]);
        }

        sixdof_obj_fnpf::fnpf_force_sum S;
        hull_sum(nb,on,computeA,S,false);
        
        fb_obj[nb]->forces_fnpf_set(p,pgc,S,computeA);

        // force log: first stage = state at the current time level
        if(iter==0)
        fb_obj[nb]->print_force_fnpf(p);
    }

    t_loads += MPI_Wtime()-t0;
}
