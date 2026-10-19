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

#include"cfd_amr.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heaviside.h"
#include"interface_width.h"
#include"momentum_rk.h"
#include"heat_print.h"
#include"concentration_io.h"
#include"fluid_update_fsf_heat.h"
#include"weno_flux_nug.h"
#include"diffusion.h"
#include"solver.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include<cmath>
#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<algorithm>

//  Heat (H 10) and concentration (C 10) with the mesh refinement: the refinement transports the
//  scalar on all grids, level 0 included (the driver leaves its heat / concentration step out),
//  in flux form, so that the fluxes through the faces between a patch and the coarser level can
//  be matched:
//
//    SSP-RK3 with the velocities of the start of the step; the flux through a face u_f s_f, the
//    face value first-order upwind (H 15 / C 15 1) or WENO5 of weno_flux_nug (4) by the sign of
//    the face velocity; the diffusion (D 20 1) explicit in the stages with the face flux
//    0.5 (k_i + k_n + (nu_t,i + nu_t,n)/sig) (s_n - s_i)/dx_p, or (D 20 2) implicit after the
//    stages (the diffusion object of the grid, backward Euler, the same face coefficients), the
//    cells around a patch Dirichlet values;
//    the flux registers: per face between a patch cell and a coarse leaf cell (the faces of the
//    composite projection, cf) the time integral of the coarse flux and of the fine fluxes (with
//    G 7 1 over the two steps of the finer level); after the steps of the finer level the coarse
//    cell takes the difference (refluxing, Berger-Colella), then the patch is restricted (the
//    mean of the fine cells: the volume mean of the halved coarse cells) onto the covered cells.
//
//  The sum of s V over the leaf cells of all grids changes only by the fluxes through the domain
//  boundary (and, for C, the concentration set to zero in the air as concentration_RK3 does).
//  The diffusivity: heat the phase values H 1, H 2 (H 9), concentration the viscosity of the
//  mixture; the eddy viscosity divided by the Prandtl number (heat sigT, concentration 1).  The
//  density and viscosity from the level set and the scalar (as fluid_update_fsf_heat with H 3 1,
//  fluid_update_fsf_concentration) also at the cells around the patches and at the patch edges
//  of the composite projection (rmix, vmix).

namespace
{
inline int cix(const lexer *q, int i, int j, int k)
{
    return (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
}

struct sscope
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    sscope(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~sscope()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};

// SSP-RK3: the weight of the fluxes of stage s in the step, the input time of the stage
const double sc_w[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
const double sc_tin[3] = {0.0, 1.0, 0.5};
}

void cfd_amr::scalar_objects(heat *ph, concentration *pc)
{
    pheat0 = ph;
    pconc0 = pc;
    if(scal==0)
    return;

    if((scal==1 && dynamic_cast<heat_print*>(ph)==nullptr) || (scal==2 && dynamic_cast<concentration_io*>(pc)==nullptr))
    {
        cout<<"CFD AMR: the heat / concentration object of level 0 holds no field"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3130);
    }
    if(scal==1)
    smat = new fluid_update_fsf_heat(p0,a0,pgc0,pheat0);
}

field& cfd_amr::sfield(int g)
{
    if(scal==1)
    return static_cast<heat_print*>(g<0 ? pheat0 : CP(g)->pheat)->T;
    return static_cast<concentration_io*>(g<0 ? pconc0 : CP(g)->pconc)->C;
}

double cfd_amr::sval(int g, int i, int j, int k)
{
    if(cav)
    return cfield(g)(i,j,k);
    return scal>0 ? sfield(g)(i,j,k) : 0.0;
}

double cfd_amr::sval(int g, int q)
{
    if(cav)
    return cfield(g).V[q];
    return scal>0 ? sfield(g).V[q] : 0.0;
}

// the density and viscosity of the mixture; with cavitation s is the vapour fraction of the liquid
// (as cavitation::rho_mix, cavitation::properties)
double cfd_amr::rmix(const lexer *q, double H, double s)
{
    if(cav)
    return ((1.0-s)*q->W1 + s*q->W62)*H + q->W3*(1.0-H);
    if(scal==0)
    return q->W1*H + q->W3*(1.0-H);
    if(scal==2)
    return (q->W1+s*q->C1)*H + (q->W3+s*q->C3)*(1.0-H);
    double r1, r2, v1, v2;
    smat->props(const_cast<lexer*>(q),s,r1,r2,v1,v2);
    return r1*H + r2*(1.0-H);
}

double cfd_amr::vmix(const lexer *q, double H, double s)
{
    if(cav)
    {
        const double r = (1.0-s)*q->W1 + s*q->W62;
        return (((1.0-s)*q->W1*q->W2 + s*q->W62*q->W63)/r)*H + q->W4*(1.0-H);
    }
    if(scal==0)
    return q->W2*H + q->W4*(1.0-H);
    if(scal==2)
    return (q->W2+s*q->C2)*H + (q->W4+s*q->C4)*(1.0-H);
    double r1, r2, v1, v2;
    smat->props(const_cast<lexer*>(q),s,r1,r2,v1,v2);
    return v1*H + v2*(1.0-H);
}

// the work fields of grid g
void cfd_amr::sc_make(int g)
{
    lexer *q = glex(g);
    cfd_amr_scgrid &G = SG(g);
    G.s0 = new field4(q);
    G.L = new field4(q);
    G.fx = new field1(q);
    G.fy = new field2(q);
    G.fz = new field3(q);
    if(scal==1)
    G.kap = new field4(q);
    if(sc_conv==4)
    G.weno = new weno_flux_nug(q);
}

void cfd_amr::sc_release(cfd_amr_scgrid &G)
{
    delete G.s0;
    delete G.L;
    delete G.fx;
    delete G.fy;
    delete G.fz;
    delete G.kap;
    delete G.weno;
    G = cfd_amr_scgrid();
}

// after the setup or a regrid: the work fields of the new grids, the flux registers (the faces of
// the composite projection)
void cfd_amr::sc_alloc()
{
    if(sc0.s0==nullptr)
    sc_make(-1);
    for(int id=0; id<(int)P.size(); ++id)
    if(CP(id)->sc==nullptr)
    {
        CP(id)->sc = new cfd_amr_scgrid;
        sc_make(id);
    }
    fr_c.assign(cf.size(),0.0);
    fr_f.assign(cf.size(),0.0);
}

void cfd_amr::sc_free()
{
    sc_release(sc0);
    for(int id=0; id<(int)P.size(); ++id)
    if(CP(id)->sc!=nullptr)
    {
        sc_release(*CP(id)->sc);
        delete CP(id)->sc;
        CP(id)->sc = nullptr;
    }
}

// heat: the diffusivity of the phases (heat_RK3::diff_update) on the whole array of grid g
void cfd_amr::sc_kappa(int g)
{
    if(scal!=1)
    return;
    lexer *q = glex(g);
    fdm *a = gfd(g);
    field &kap = *SG(g).kap;
    const double a1 = (q->H9==2) ? q->H2 : q->H1;
    const double a2 = (q->H9==2) ? q->H1 : q->H2;
    const int m = q->margin;
    for(int i=-m; i<q->knox+m; ++i)
    for(int j=(q->j_dir==1 ? -m : 0); j<(q->j_dir==1 ? q->knoy+m : q->knoy); ++j)
    for(int k=-m; k<q->knoz+m; ++k)
    {
        const int iw = std::min(std::max(i,1-m),q->knox+m-2);
        const int jw = q->j_dir==1 ? std::min(std::max(j,1-m),q->knoy+m-2) : j;
        const int kw = std::min(std::max(k,1-m),q->knoz+m-2);
        const double H = heaviside(a->phi(i,j,k),interface_width(q,a->phi,iw,jw,kw));
        kap(i,j,k) = a1*H + a2*(1.0-H);
    }
}

// the fluxes through the faces of the cells of grid g (per unit area, along +d) from its scalar,
// and L = -div
void cfd_amr::sc_flux(int g)
{
    lexer *q = glex(g);
    fdm *a = gfd(g);
    cfd_amr_scgrid &G = SG(g);
    field &S = sfield(g);
    field &kap = (scal==1) ? *G.kap : a->visc;
    field &ev = a->eddyv;
    field *F[3] = {G.fx,G.fy,G.fz};
    field *U[3] = {&a->u,&a->v,&a->w};
    const double *dp[3] = {q->DXP,q->DYP,q->DZP};
    const double *dn[3] = {q->DXN,q->DYN,q->DZN};
    const int n[3] = {q->knox,q->knoy,q->knoz};

    auto act = [&](int i, int j, int k) -> bool
    {
        return i>=0 && i<n[0] && j>=0 && j<n[1] && k>=0 && k<n[2] && q->flag4[cix(q,i,j,k)]>0;
    };

    for(int d=0; d<3; ++d)
    {
        if(d==1 && q->j_dir==0)
        continue;
        int e[3] = {0,0,0};
        e[d] = 1;
        field &Fd = *F[d];
        field &Ud = *U[d];
        int lo[3] = {0,0,0};
        lo[d] = -1;
        for(int i=lo[0]; i<n[0]; ++i)
        for(int j=lo[1]; j<n[1]; ++j)
        for(int k=lo[2]; k<n[2]; ++k)
        {
            const int i2=i+e[0], j2=j+e[1], k2=k+e[2];
            if(!act(i,j,k) && !act(i2,j2,k2))
            {
                Fd(i,j,k) = 0.0;
                continue;
            }
            const double vel = Ud(i,j,k);
            double f = 0.0;
            if(sc_conv==1)
            f = vel*(vel>=0.0 ? S(i,j,k) : S(i2,j2,k2));
            if(sc_conv==4)
            f = vel*G.weno->face_value(q,a,S,Ud,d,i,j,k,vel);
            if(sc_diff==1)
            {
                const int L = (d==0 ? i : (d==1 ? j : k)) + marge;
                const double kf = 0.5*(kap(i,j,k)+kap(i2,j2,k2) + (ev(i,j,k)+ev(i2,j2,k2))/sc_sig);
                f -= kf*(S(i2,j2,k2)-S(i,j,k))/dp[d][L];
            }
            Fd(i,j,k) = f;
        }
    }

    field &Lf = *G.L;
    for(int i=0; i<n[0]; ++i)
    for(int j=0; j<n[1]; ++j)
    for(int k=0; k<n[2]; ++k)
    {
        if(!act(i,j,k))
        continue;
        double v = -((*G.fx)(i,j,k)-(*G.fx)(i-1,j,k))/dn[0][i+marge]
                   -((*G.fz)(i,j,k)-(*G.fz)(i,j,k-1))/dn[2][k+marge];
        if(q->j_dir==1)
        v -= ((*G.fy)(i,j,k)-(*G.fy)(i,j-1,k))/dn[1][j+marge];
        Lf(i,j,k) = v;
    }
}

namespace
{
// the area of the face of a patch-edge entry (the fine face)
inline double cf_area(const lexer *q, const cfd_amr_cf &F, int mg)
{
    const double *dn[3] = {q->DXN,q->DYN,q->DZN};
    double A = 1.0;
    for(int t=0; t<3; ++t)
    if(t!=F.d)
    A *= dn[t][F.fi[t]+mg];
    return A;
}
}

// the fluxes of grid g through the faces of the flux registers: coarse (g the coarse cell of the
// face) or fine (g the patch), for the faces of the patches of level lfine; weight w (time)
void cfd_amr::sc_record(int g, double w, bool coarse, int lfine)
{
    if(lfine<1 || lfine>maxlev)
    return;
    cfd_amr_scgrid &G = SG(g);
    field *Fs[3] = {G.fx,G.fy,G.fz};
    for(size_t m=0; m<cf.size(); ++m)
    {
        const cfd_amr_cf &F = cf[m];
        if(glevel(F.pid)!=lfine)
        continue;
        const double A = cf_area(CP(F.pid)->pp,F,marge);
        if(coarse && F.gc==g)
        {
            int c[3] = {F.ci[0],F.ci[1],F.ci[2]};
            if(F.dir>0)
            c[F.d] -= 1;
            fr_c[m] += w*(*Fs[F.d])(c[0],c[1],c[2])*A;
        }
        if(!coarse && F.pid==g)
        fr_f[m] += w*(*Fs[F.d])(F.ff[0],F.ff[1],F.ff[2])*A;
    }
}

// the implicit diffusion: the diffusive fluxes of the solution through the faces of the registers
void cfd_amr::sc_record_diff(int g, double w, bool coarse, int lfine)
{
    if(lfine<1 || lfine>maxlev)
    return;
    lexer *q = glex(g);
    fdm *a = gfd(g);
    field &S = sfield(g);
    field &kap = (scal==1) ? *SG(g).kap : a->visc;
    field &ev = a->eddyv;
    const double *dp[3] = {q->DXP,q->DYP,q->DZP};

    auto flux = [&](int i, int j, int k, int d) -> double
    {
        int e[3] = {0,0,0};
        e[d] = 1;
        const int i2=i+e[0], j2=j+e[1], k2=k+e[2];
        const int L = (d==0 ? i : (d==1 ? j : k)) + marge;
        const double kf = 0.5*(kap(i,j,k)+kap(i2,j2,k2) + (ev(i,j,k)+ev(i2,j2,k2))/sc_sig);
        return -kf*(S(i2,j2,k2)-S(i,j,k))/dp[d][L];
    };

    for(size_t m=0; m<cf.size(); ++m)
    {
        const cfd_amr_cf &F = cf[m];
        if(glevel(F.pid)!=lfine)
        continue;
        const double A = cf_area(CP(F.pid)->pp,F,marge);
        if(coarse && F.gc==g)
        {
            int c[3] = {F.ci[0],F.ci[1],F.ci[2]};
            if(F.dir>0)
            c[F.d] -= 1;
            fr_c[m] += w*flux(c[0],c[1],c[2],F.d)*A;
        }
        if(!coarse && F.pid==g)
        fr_f[m] += w*flux(F.ff[0],F.ff[1],F.ff[2],F.d)*A;
    }
}

// the scalar of grid g at the start of the step
void cfd_amr::sc_save(int g)
{
    lexer *q = glex(g);
    const size_t n = (size_t)q->imax*q->jmax*q->kmax;
    const double *src = sfield(g).V;
    double *dst = SG(g).s0->V;
    for(size_t c=0; c<n; ++c)
    dst[c] = src[c];
}

// stage s of the SSP-RK3 on grid g
void cfd_amr::sc_update(int g, int s, double dt)
{
    lexer *q = glex(g);
    field &S = sfield(g), &S0 = *SG(g).s0, &L = *SG(g).L;
    const double a0 = (s==0) ? 0.0 : (s==1 ? 0.75 : 1.0/3.0);
    const double b0 = 1.0-a0;
    for(int i=0; i<q->knox; ++i)
    for(int j=0; j<q->knoy; ++j)
    for(int k=0; k<q->knoz; ++k)
    {
        const int c = cix(q,i,j,k);
        if(q->flag4[c]<=0)
        continue;
        S.V[c] = a0*S0.V[c] + b0*(S.V[c] + dt*L.V[c]);
    }
}

// boundary values of the scalar of grid g (the code of heat_RK3 / concentration_RK3, 80; the
// concentration zero in the air as bc_concentration)
void cfd_amr::sc_bc(int g)
{
    lexer *q = glex(g);
    fdm *a = gfd(g);
    field &S = sfield(g);
    if(scal==2)
    for(int i=0; i<q->knox; ++i)
    for(int j=0; j<q->knoy; ++j)
    for(int k=0; k<q->knoz; ++k)
    {
        const int c = cix(q,i,j,k);
        if(q->flag4[c]>0 && a->phi.V[c]<0.0)
        S.V[c] = 0.0;
    }
    if(g<0)
    pgc0->start4(p0,S,80);
    else
    {
        sscope ps(pgc0,a,a0);
        pgc0->start4(q,S,80);
    }
}

// the implicit diffusion of grid g (D 20 2): the diffusion object of the grid, backward Euler over
// the step (the cells around a patch Dirichlet values, amr_dsolver)
void cfd_amr::sc_diffuse(int g)
{
    lexer *q = glex(g);
    fdm *a = gfd(g);
    field &S = sfield(g);
    field &kap = (scal==1) ? *SG(g).kap : a->visc;
    diffusion *pd = (g<0) ? mom0->amr_diffusion() : CP(g)->pdiff;
    solver *ps = (g<0) ? mom0->amr_solver() : CP(g)->psolv;

    int n = 0;
    for(int i=0; i<q->knox; ++i)
    for(int j=0; j<q->knoy; ++j)
    for(int k=0; k<q->knoz; ++k)
    if(q->flag4[cix(q,i,j,k)]>0)
    a->rhsvec.V[n++] = 0.0;

    if(g<0)
    pd->diff_scalar(q,a,pgc0,ps,S,kap,a->eddyv,sc_sig,1.0);
    else
    {
        sscope sp(pgc0,a,a0);
        pd->diff_scalar(q,a,pgc0,ps,S,kap,a->eddyv,sc_sig,1.0);
    }
}

// the cells around the patches of level l: mode 0 from the scalar of the coarser grids (and the
// patches of the level), mode 1 from the time level of the coarser grids (G 7 1, after interp)
void cfd_amr::sc_fill(int l, int mode)
{
    vector<fspec> fs;
    if(mode==0)
    fs.push_back({CELL,1,[this](int g) -> field& { return sfield(g); }});
    else
    {
        const int m = tl_s;
        fs.push_back({CELL,1,[this,l,m](int g) -> field& { return glevel(g)==l ? sfield(g) : tsf(g,m); }});
    }
    fill(l,fs);
}

void cfd_amr::sc_reset(int lfine)
{
    for(size_t m=0; m<cf.size(); ++m)
    if(lfine<0 || glevel(cf[m].pid)==lfine)
    fr_c[m] = fr_f[m] = 0.0;
}

// the flux registers of the patch edges of level l: the coarse cell takes the difference of the
// fine and the coarse flux through the face
void cfd_amr::sc_reflux(int l)
{    // REEF3D_AMR_REFLUX=0 in the environment: no refluxing (a diagnostic for the comparison)
    static const bool off = (getenv("REEF3D_AMR_REFLUX")!=nullptr && atoi(getenv("REEF3D_AMR_REFLUX"))==0);
    if(off)
    return;

    // G 40: the coarse cells of the shadows collect their difference, added on the ranks that hold
    // the cells
    const bool up = (place && l==1);
    const std::function<field&(int)> fl = [this](int g) -> field& { return *SG(g).L; };
    if(up)
    sh_zero(pu_cc,fl);

    for(size_t m=0; m<cf.size(); ++m)
    {
        const cfd_amr_cf &F = cf[m];
        if(glevel(F.pid)!=l)
        continue;
        lexer *qc = glex(F.gc);
        const double V = qc->DXN[F.ci[0]+marge]*qc->DYN[F.ci[1]+marge]*qc->DZN[F.ci[2]+marge];
        field &S = (up && is_shadow(F.gc)) ? *SG(F.gc).L : sfield(F.gc);
        S.V[F.qc] += double(F.dir)*(fr_f[m]-fr_c[m])/V;
    }

    if(up)
    {
        vector<std::function<field&(int)>> f(1,[this](int g) -> field& { return g<0 ? sfield(g) : *SG(g).L; });
        sh_run(pu_cc,f,true);
    }
}

// G 40: the coarse fluxes of level 0 through the faces at the edges of the level-1 patches, recorded
// in the registers from the shadows (diff: the diffusive fluxes of the implicit diffusion)
void cfd_amr::sc_record_shadows(double w, bool diff)
{
    if(!place)
    return;
    vector<std::function<field&(int)>> f;
    if(!diff)
    {
        f.push_back([this](int g) -> field& { return *SG(g).fx; });
        f.push_back([this](int g) -> field& { return *SG(g).fy; });
        f.push_back([this](int g) -> field& { return *SG(g).fz; });
    }
    else
    {
        f.push_back([this](int g) -> field& { return sfield(g); });
        f.push_back([this](int g) -> field& { return gfd(g)->eddyv; });
        if(scal==1)
        f.push_back([this](int g) -> field& { return *SG(g).kap; });
        else
        f.push_back([this](int g) -> field& { return gfd(g)->visc; });
    }
    sh_down_f(f);
    for(int g=npat; g<(int)P.size(); ++g)
    if(diff)
    sc_record_diff(g,w,true,1);
    else
    sc_record(g,w,true,1);
}

// G 7 0: the step of the scalar on all grids (the step of the finest level)
void cfd_amr::sc_step_comp()
{
    const double dt = p0->dt;
    const int ng = npat;

    sc_reset(-1);
    for(int g=-1; g<ng; ++g)
    {
        sc_kappa(g);
        sc_save(g);
    }

    for(int s=0; s<3; ++s)
    {
        for(int l=1; l<=maxlev; ++l)
        sc_fill(l,0);

        for(int g=-1; g<ng; ++g)
        {
            sc_flux(g);
            const int lg = glevel(g);
            sc_record(g,sc_w[s]*dt,true,lg+1);
            sc_record(g,sc_w[s]*dt,false,lg);
        }
        sc_record_shadows(sc_w[s]*dt,false);
        for(int g=-1; g<ng; ++g)
        {
            sc_update(g,s,dt);
            sc_bc(g);
        }
    }

    if(sc_diff==2)
    for(int l=0; l<=maxlev; ++l)
    {
        if(l>=1)
        sc_fill(l,0);
        const vector<int> grids = (l==0) ? vector<int>{-1} : lev[l];
        for(int g : grids)
        {
            sc_diffuse(g);
            sc_bc(g);
            sc_record_diff(g,dt,true,l+1);
            sc_record_diff(g,dt,false,l);
        }
        if(l==0)
        sc_record_shadows(dt,true);
    }

    vector<std::function<field&(int)>> cells, none;
    cells.push_back([this](int g) -> field& { return sfield(g); });
    for(int l=maxlev; l>=1; --l)
    {
        sc_reflux(l);
        restrict_level(l,cells,none);
    }

    for(int g=-1; g<ng; ++g)
    sc_bc(g);
    for(int l=1; l<=maxlev; ++l)
    sc_fill(l,0);
}

// G 7 1: the step of the scalar on the grids of level l (step dt, from th0 in the step of the
// coarser level); the flux registers of the faces of its patches start with the coarse fluxes
void cfd_amr::sc_step_level(int l, double dt, double th0)
{
    const vector<int> grids = (l==0) ? vector<int>{-1} : lev[l];

    sc_reset(l+1);
    for(int g : grids)
    {
        sc_kappa(g);
        sc_save(g);
    }

    for(int s=0; s<3; ++s)
    {
        if(l>=1)
        {
            interp(l-1,th0+0.5*sc_tin[s]);
            sc_fill(l,1);
        }
        for(int g : grids)
        {
            sc_flux(g);
            sc_record(g,sc_w[s]*dt,true,l+1);
            sc_record(g,sc_w[s]*dt,false,l);
        }
        if(l==0)
        sc_record_shadows(sc_w[s]*dt,false);
        for(int g : grids)
        {
            sc_update(g,s,dt);
            sc_bc(g);
        }
    }

    if(sc_diff==2)
    {
        if(l>=1)
        {
            interp(l-1,th0+0.5);
            sc_fill(l,1);
        }
        for(int g : grids)
        {
            sc_diffuse(g);
            sc_bc(g);
            sc_record_diff(g,dt,true,l+1);
            sc_record_diff(g,dt,false,l);
        }
        if(l==0)
        sc_record_shadows(dt,true);
    }
}

// the sum of s V over the leaf cells of all grids (smin, smax: the range of s there)
double cfd_amr::sc_total(double &smin, double &smax)
{
    double sum = 0.0;
    smin = 1.0e300;
    smax = -1.0e300;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        field &S = sfield(g);
        for(int c : lq[g+1])
        {
            const int jk = q->jmax*q->kmax;
            const int i = c/jk + q->imin;
            const int r = c%jk;
            const int j = r/q->kmax + q->jmin;
            const int k = r%q->kmax + q->kmin;
            sum += S.V[c]*q->DXN[i+marge]*q->DYN[j+marge]*q->DZN[k+marge];
            smin = std::min(smin,S.V[c]);
            smax = std::max(smax,S.V[c]);
        }
    }
    smin = pgc0->globalmin(smin);
    smax = pgc0->globalmax(smax);
    return pgc0->globalsum(sum);
}
