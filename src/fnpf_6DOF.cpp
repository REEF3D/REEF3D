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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_6DOF.h"
#include"6DOF_obj.h"
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
#include<functional>

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

fnpf_6DOF::fnpf_6DOF(lexer *p, fdm_fnpf *c, ghostcell *pgc) : initialized(false), foot(p), psiD(p), zeroslice(p),
                                                               eta_ext(p), fi_ext(p), amr(nullptr), amr_layout(-1)
{
    if(p->mpirank==0)
    cout<<"6DOF FNPF startup ..."<<endl;
    
    // one body, as in NHFLOW
    nbody = 1;
    
    for(int nb=0; nb<nbody; ++nb)
    fb_obj.push_back(new sixdof_obj(p,pgc,nb));
    
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
    
    pvel->velcalc_sig(p,c,pgc,c->Fi);
    pgc->gcparax7(p,c->U,7);
    pgc->gcparax7(p,c->V,7);
    pgc->gcparax7(p,c->W,7);
    
    if(amr_on())
    {
        g0.pf = pf;
        
        reefamr_comms_off guard(pgc);
        for(auto &G : gp)
        G.pvel->velcalc_sig(G.p,G.c,pgc,G.c->Fi);
        
        forces_amr(p,c,pgc,psolv,pf,Keta,Kfi,iter);
    }
    else
    forces(p,c,pgc,psolv,pf,Keta,Kfi,iter);
    
    motion(p,c,pgc,iter);
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

fnpf_6DOF::~fnpf_6DOF()
{
    for(auto &G : gp)
    free_grid(G);
}

void fnpf_6DOF::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->initialize_fnpf(p,c,pgc);
    
    geometry(g0,pgc);
    extrapolate(g0,pgc,c->Fi);
    
    initialized = true;
}

void fnpf_6DOF::forces(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, 
                       slice &Keta, slice &Kfi, int iter)
{
    // psi_0: phi_t at fixed z on the free surface, phi_t|z = dFifsf/dt - Fz*deta/dt
    SLICELOOP4
    psiD(i,j) = Kfi(i,j) - c->Fz(i,j)*Keta(i,j);
    
    pgc->gcsl_start4(p,psiD,50);
    
    zero_face(g0);
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->face_data_fnpf(p,c,pgc,-2,c->FBF,c->FBu,c->FBv,c->FBw);
    exchange_face(g0,pgc);
    
    solve_psi(p,c,pgc,psolv,plap,pf,psi0,psiD);
    
    // added mass: six unit-mode solves, once per time step
    const bool refresh = (iter==0);
    
    for(int nb=0; nb<nbody; ++nb)
    {
        const bool computeA = refresh && !fb_obj[nb]->fnpf_fixed(p);
        
        if(computeA)
        for(int m=0; m<6; ++m)
        if(psi[m]!=nullptr)
        {
            zero_face(g0);
            fb_obj[nb]->face_data_fnpf(p,c,pgc,m,c->FBF,c->FBu,c->FBv,c->FBw);
            exchange_face(g0,pgc);
            
            solve_psi(p,c,pgc,psolv,plap,pf,psi[m],zeroslice);
        }
        
        fb_obj[nb]->forces_fnpf(p,c,pgc,psi0,psi,computeA);
        
        // force log: first stage = state at the current time level
        if(iter==0)
        fb_obj[nb]->print_force_fnpf(p);
    }
}

void fnpf_6DOF::motion(lexer *p, fdm_fnpf *c, ghostcell *pgc, int iter)
{
    // last stage: 2 for RK3 (A 310 3), 3 for RK4 (A 310 4)
    const bool finalize = (iter == ((p->A310==3) ? 2 : 3));
    
    for(int nb=0; nb<nbody; ++nb)
    {
        fb_obj[nb]->solve_eqmotion_fnpf(p,pgc,iter,finalize);
        fb_obj[nb]->update_position_fnpf(p,pgc,finalize);
        
        if(finalize)
        fb_obj[nb]->print_fnpf(p,pgc,iter);
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
    // All body nodes, layer by layer from the fluid inwards: average of the fluid (or
    // already filled) neighbours. The first layers feed the lagged cross terms of the
    // Laplace rhs and the sampling of the hull pressure.
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
                int nbr[6] = {q-sI, q+sI, q-sJ, q+sJ, q-1, q+1};
                const int nn = (p->j_dir==1) ? 6 : 4;
                
                if(p->j_dir==0)
                {
                nbr[2] = q-1;
                nbr[3] = q+1;
                }
                
                double sum=0.0;
                int cnt=0;
                
                for(int m=0; m<nn; ++m)
                {
                    const int r = nbr[m];
                    
                    if(p->flag7[r]>0 && (FBF[r]<0.5 || (mark[r]>0 && mark[r]<=pass)))
                    {
                    sum += f[r];
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
                          double *f, slice &D)
{
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
    
    plap->start(p,c,pgc,psolv,pf,f,D);
    
    c->Uin = uin;
    
    pgc->start7V(p,f,c->bc,gcval);
    
    extrapolate(g0,pgc,f);
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
    lexer *p = G.p;
    const int size = p->imax*p->jmax*(p->kmax+2);
    
    p->del_Darray(G.psi0,size);
    p->del_Darray(G.zero,size);
    for(int m=0; m<6; ++m)
    if(G.psi[m]!=nullptr)
    p->del_Darray(G.psi[m],size);
    delete [] G.mark;
    delete G.foot; delete G.psiD; delete G.zeroslice; delete G.eta_ext; delete G.fi_ext;
    delete G.pbed;
    delete G.pvel;
}

// body grids of the patches, after the (re)gridding of fnpf_amr
void fnpf_6DOF::amr_grids(lexer *p, ghostcell *pgc)
{
    if(amr==nullptr || amr_layout==amr->layout())
    return;
    
    for(auto &G : gp)
    free_grid(G);
    gp.clear();
    
    reefamr_comms_off guard(pgc);
    
    for(int n=0; n<amr->patches(); ++n)
    {
        fnpf_6DOF_grid G;
        lexer *pp = amr->patch_lexer(n);
        fdm_fnpf *cc = amr->patch_fdm(n);
        const int size = pp->imax*pp->jmax*(pp->kmax+2);
        
        G.p = pp;
        G.c = cc;
        G.pf = amr->patch_fsf(n);
        G.id = n;
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
        
        gp.push_back(G);
    }
    
    amr_layout = amr->layout();
    
    if(initialized)
    for(auto &G : gp)
    {
        geometry(G,pgc);
        extrapolate(G,pgc,G.c->Fi);
        amr->patch_walls_fi(G.id,G.c->Fi);
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
void fnpf_6DOF::solve_psi_amr(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, int m, bool withD)
{
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
    
    amr->lap_solve_psi(p,c,pgc,psolv,pf,&fs[0],&Ds[0]);
    
    c->Uin = uin;
    
    pgc->start7V(p,fs[0],c->bc,gcval);
    extrapolate(g0,pgc,fs[0]);
    
    reefamr_comms_off guard(pgc);
    for(size_t n=0; n<gp.size(); ++n)
    {
        amr->patch_walls_fi(gp[n].id,fs[n+1]);
        extrapolate(gp[n],pgc,fs[n+1]);
    }
}

void fnpf_6DOF::forces_amr(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, 
                           slice &Keta, slice &Kfi, int iter)
{
    amr_grids(p,pgc);
    
    // psi_0: phi_t at fixed z on the free surface, phi_t|z = dFifsf/dt - Fz*deta/dt
    SLICELOOP4
    psiD(i,j) = Kfi(i,j) - c->Fz(i,j)*Keta(i,j);
    
    pgc->gcsl_start4(p,psiD,50);
    
    {
    reefamr_comms_off guard(pgc);
    for(auto &G : gp)
    {
        lexer *p = G.p;
        slice4 &D = *G.psiD;
        slice &Ke = *G.Keta, &Kf = *G.Kfi;
        
        SLICELOOP4
        D(i,j) = Kf(i,j) - G.c->Fz(i,j)*Ke(i,j);
        
        pgc->gcsl_start4(p,D,50);
    }
    }
    
    auto face = [&](int mode)
    {
        zero_face(g0);
        for(int nb=0; nb<nbody; ++nb)
        fb_obj[nb]->face_data_fnpf(p,c,pgc,mode,c->FBF,c->FBu,c->FBv,c->FBw);
        exchange_face(g0,pgc);
        
        reefamr_comms_off guard(pgc);
        for(auto &G : gp)
        {
            zero_face(G);
            for(int nb=0; nb<nbody; ++nb)
            fb_obj[nb]->face_data_fnpf(G.p,G.c,pgc,mode,G.c->FBF,G.c->FBu,G.c->FBv,G.c->FBw);
        }
    };
    
    face(-2);
    solve_psi_amr(p,c,pgc,psolv,pf,-1,true);
    
    // added mass: six unit-mode solves, once per time step
    const bool refresh = (iter==0);
    
    for(int nb=0; nb<nbody; ++nb)
    {
        const bool computeA = refresh && !fb_obj[nb]->fnpf_fixed(p);
        
        if(computeA)
        for(int m=0; m<6; ++m)
        if(psi[m]!=nullptr)
        {
            face(m);
            solve_psi_amr(p,c,pgc,psolv,pf,m,false);
        }
        
        // every hull triangle on the finest grid that holds its centroid
        sixdof_obj::fnpf_force_sum S;
        fb_obj[nb]->forces_fnpf_zero(p,S);
        
        for(int g=-1; g<(int)gp.size(); ++g)
        {
            fnpf_6DOF_grid &G = (g<0) ? g0 : gp[g];
            const int id = G.id;
            std::function<bool(double,double)> own = [&](double x, double y)
            {
                if(!(x >= p->originx && x < p->endx))
                return false;
                if(p->j_dir==1 && !(y >= p->originy && y < p->endy))
                return false;
                return amr->finest_at(x,y)==id;
            };
            const double del = (g<0) ? 0.5*fb_obj[nb]->fnpf_dsm() : G.del;
            fb_obj[nb]->forces_fnpf_sum(G.p,G.c,G.psi0,G.psi,computeA,del,&own,S);
        }
        
        fb_obj[nb]->forces_fnpf_set(p,pgc,S,computeA);
        
        // force log: first stage = state at the current time level
        if(iter==0)
        fb_obj[nb]->print_force_fnpf(p);
    }
}
