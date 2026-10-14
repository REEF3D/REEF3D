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
#include"turbulence.h"
#include"6DOF_cfd.h"
#include"cfd_amr_6dof.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include<cmath>
#include<iostream>
#include<iomanip>
#include<algorithm>

//  Subcycling in time (G 7 1) and the composite volume correction of the level set (F 46 2, 3).
//
//  Every level takes its own step: level l+1 two steps of half the size per step of level l
//  (Berger-Oliger).  A step of level l: its SSP stages (level set, momentum, the projection of
//  the level alone), then two steps of level l+1, then the restriction of level l+1 onto level l.
//  The cells around a patch get the values of the coarser level linear in time between the start
//  and the end of its step, at the time of the stage (SSP-RK3: input times 0, 1, 1/2, output times
//  1, 1/2, 1 of the step); the pressure linear between the middles of the steps.  The projection
//  of a level (cfd_amr_press.cpp, pr_level) takes the pressure of the coarser level at the patch
//  edges (Dirichlet, its hydrostatic part as in G 7 0): the patch makes its own flux there, so its
//  pressure answers what happens inside it (a floating body).  At the end of the level-0 step the
//  synchronisation projection (pr_sync_project), one composite correction solve over all levels,
//  makes the coarse faces at the patch edges the mean of the fine ones and the composite velocity
//  field divergence free.  The fine grid acts on the coarse one through the restriction
//  (velocities, level set, density, viscosity) and the synchronisation.  Floating bodies (X 10 1):
//  the finest level advances them, the coarser levels step with a predicted copy (the loads of the
//  last finest stage).

namespace
{
struct pscope3
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    pscope3(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~pscope3()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};

inline void cdecode(const lexer *q, int qq, int &i, int &j, int &k)
{
    const int jk = q->jmax*q->kmax;
    i = qq/jk + q->imin;
    const int r = qq%jk;
    j = r/q->kmax + q->jmin;
    k = r%q->kmax + q->kmin;
}
}

cfd_amr_tlev::cfd_amr_tlev(lexer *q, int nn) : n(nn)
{
    o[0] = new field1(q); o[1] = new field2(q); o[2] = new field3(q);
    t[0] = new field1(q); t[1] = new field2(q); t[2] = new field3(q);
    for(int m=3; m<8; ++m)
    {
        o[m] = (m<n) ? new field4(q) : nullptr;
        t[m] = (m<n) ? new field4(q) : nullptr;
    }
}

cfd_amr_tlev::~cfd_amr_tlev()
{
    for(int m=0; m<8; ++m)
    {
        delete o[m];
        delete t[m];
    }
}

// the fields of grid g that the finer level interpolates: u, v, w, phi, press, eddyv; kin, eps
field& cfd_amr::stf(int g, int m)
{
    fdm *a = gfd(g);
    if(m==0) return a->u;
    if(m==1) return a->v;
    if(m==2) return a->w;
    if(m==3) return a->phi;
    if(m==4) return a->press;
    if(m==5) return a->eddyv;
    return kfield(g,m-6);
}

namespace
{
inline size_t asize(const lexer *q) { return (size_t)q->imax*q->jmax*q->kmax; }
}

// the state of the grids of level l at the start of their step
void cfd_amr::save_old(int l)
{
    vector<int> gs;
    if(l==0)
    gs.push_back(-1);
    else
    gs = lev[l];

    for(int g : gs)
    {
        if(gts(g)==nullptr)
        continue;
        const size_t n = asize(glex(g));
        for(int m=0; m<gts(g)->n; ++m)
        std::copy(stf(g,m).V,stf(g,m).V+n,gts(g)->o[m]->V);
    }
}

// the values of the grids of level l at the time th of their step (0 start, 1 end)
void cfd_amr::interp(int l, double th)
{
    vector<int> gs;
    if(l==0)
    gs.push_back(-1);
    else
    gs = lev[l];

    for(int g : gs)
    {
        if(gts(g)==nullptr)
        continue;
        const size_t n = asize(glex(g));
        for(int m=0; m<gts(g)->n; ++m)
        {
            const double *o = gts(g)->o[m]->V, *c = stf(g,m).V;
            double *t = gts(g)->t[m]->V;
            // the pressure of a step belongs to its middle: linear between the middles of the
            // steps (the start value at -1/2, the end value at 1/2), beyond the middle extrapolated
            // from level 0, held at the end value from a patch; in the first step the end value
            // (the initial pressure is not that of a step).  The pressure of a patch follows the
            // edge pressure of its own coarser level and the synchronisation, it is not smooth from
            // step to step: extrapolated it fed its jumps back through the edges of the finer level
            // (two levels, the free surface across the edge of level 2: unstable after 1.3 s)
            double w = (m!=4) ? th : (press_step ? th+0.5 : 1.0);
            if(m==4 && l>=1)
            w = MIN(w,1.0);
            for(size_t q=0; q<n; ++q)
            t[q] = o[q] + w*(c[q]-o[q]);
        }
    }
}

// the cells around the patches of level l at the time th of the step of level l-1
void cfd_amr::fill_level(int l, double th)
{
    interp(l-1,th);

    vector<fspec> fs;
    const int mode[8] = {1,1,1,0,0,1,1,1};
    for(int m=0; m<nts; ++m)
    fs.push_back({m<3 ? m : CELL,mode[m],[this,m,l](int g) -> field& { return glevel(g)==l ? stf(g,m) : tsf(g,m); }});
    fill(l,fs);

    for(int id : lev[l])
    margin_rovisc(id);
}

// restriction of level l onto the grids of level l-1 at the end of their step: velocities (faces),
// level set, density, viscosity
void cfd_amr::restrict_sub(int l)
{
    vector<std::function<field&(int)>> cells, faces;
    cells.push_back([this](int g) -> field& { return gfd(g)->phi; });
    cells.push_back([this](int g) -> field& { return gfd(g)->ro; });
    cells.push_back([this](int g) -> field& { return gfd(g)->visc; });
    if(turb)
    {
        cells.push_back([this](int g) -> field& { return gfd(g)->eddyv; });
        cells.push_back([this](int g) -> field& { return kfield(g,0); });
        cells.push_back([this](int g) -> field& { return kfield(g,1); });
    }
    faces.push_back([this](int g) -> field& { return gfd(g)->u; });
    faces.push_back([this](int g) -> field& { return gfd(g)->v; });
    faces.push_back([this](int g) -> field& { return gfd(g)->w; });

    restrict_level(l,cells,faces);

    if(l==1)
    {
        pgc0->start1(p0,a0->u,10);
        pgc0->start2(p0,a0->v,11);
        pgc0->start3(p0,a0->w,12);
        pgc0->start4(p0,a0->phi,gcval_phi);
        pgc0->start4(p0,a0->ro,1);
        pgc0->start4(p0,a0->visc,1);
        if(turb)
        pturb0->gcupdate(p0,a0,pgc0);
    }
    else
    fill_level(l-1,th_end[l-1]);
}

// one step of size dt of the patches of level l, starting at th0 of the step of level l-1
void cfd_amr::patch_step(int l, double dt, double th0, double t)
{
    ghostcell *pgc = pgc0;

    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        pp->dt_old = (c->dt_l>0.0) ? c->dt_l : dt;
        pp->dt = dt;
        c->dt_l = dt;
        pp->simtime = t;
        pp->count = p0->count;
    }

    // the turbulence model (the cells around the patches filled at the start of the step)
    if(turb)
    turb_step(l);

    for(int id : lev[l])
    {
        pscope3 ps(pgc,CP(id)->a,a0);
        CP(id)->pmom->amr_step_begin(CP(id)->pp,CP(id)->a,pgc);
    }

    const int stages = mom0->amr_stages();

    // floating bodies: a level coarser than the finest steps with a predicted copy
    if(bodies && l<maxlev)
    for(int nb=0; nb<b6->objects(); ++nb)
    b6->object(nb)->amr_save();

    for(int s=0; s<stages; ++s)
    {
        const double thi = th0 + 0.5*tin[s];
        const double tho = th0 + 0.5*tout[s];

        // stage input around the patches: the coarser level at the input time of the stage
        interp(l-1,thi);
        {
            vector<fspec> fs;
            for(int c=0; c<3; ++c)
            fs.push_back({c,1,[this,s,c,l](int g) -> field& { return glevel(g)==l ? gmom(g)->amr_vel(gfd(g),c,s) : tsf(g,c); }});
            fs.push_back({CELL,0,[this,l](int g) -> field& { return glevel(g)==l ? gfd(g)->phi : tsf(g,3); }});
            fs.push_back({CELL,0,[this,l](int g) -> field& { return glevel(g)==l ? gfd(g)->press : tsf(g,4); }});
            fs.push_back({CELL,1,[this,l](int g) -> field& { return glevel(g)==l ? gfd(g)->eddyv : tsf(g,5); }});
            if(s==0)
            fs.push_back({CELL,0,[this,l](int g) -> field& { return glevel(g)==l ? gfd(g)->phi : tsf(g,3); },
                                 [this](int g) -> field& { return gmom(g)->amr_phi_in(0); }});
            else
            fs.push_back({CELL,0,[this,s,l](int g) -> field& { return glevel(g)==l ? gmom(g)->amr_phi_in(s) : tsf(g,3); }});
            fill(l,fs);
            for(int id : lev[l])
            margin_rovisc(id);
        }

        for(int id : lev[l])
        {
            pscope3 ps(pgc,CP(id)->a,a0);
            CP(id)->pmom->amr_ls_transport(CP(id)->pp,CP(id)->a,pgc,s);
        }

        // the coarser level at the output time of the stage
        interp(l-1,tho);

        vector<fspec> fo;
        fo.push_back({CELL,0,[this,s,l](int g) -> field& { return glevel(g)==l ? gmom(g)->amr_phi_out(s) : tsf(g,3); }});
        fill(l,fo);

        const int iters = mom0->amr_reini_iters(p0,s);
        for(int it=0; it<iters; ++it)
        {
            if(it>0)
            fill(l,fo);

            for(int id : lev[l])
            {
                pscope3 ps(pgc,CP(id)->a,a0);
                if(it==0)
                CP(id)->pmom->amr_ls_finish(CP(id)->pp,CP(id)->a,pgc,s,false,1);
                else
                CP(id)->pmom->amr_ls_reini(CP(id)->pp,CP(id)->a,pgc,s,1);
            }
        }

        // momentum: the Dirichlet values of the implicit diffusion around the patch
        {
            vector<fspec> fs;
            for(int c=0; c<3; ++c)
            fs.push_back({c,1,[this,s,c,l](int g) -> field& { return glevel(g)==l ? gmom(g)->amr_vel(gfd(g),c,s) : tsf(g,c); },
                              [this,s,c](int g) -> field& { return gmom(g)->amr_velout(gfd(g),c,s); }});
            fill(l,fs);
        }
        edge_faces(l,s);

        // the bodies at the output time of the stage
        if(bodies)
        body_stage(l,s,t,dt);

        for(int id : lev[l])
        {
            pscope3 ps(pgc,CP(id)->a,a0);
            CP(id)->pmom->amr_momentum(CP(id)->pp,CP(id)->a,pgc,CP(id)->pvrans,(sixdof*)CP(id)->p6dof,s,patch_forcing());
        }

        // the projection of the level with the coarse flux at the patch edges
        pr_level(p0,pgc,l,s);

        for(int id : lev[l])
        {
            pscope3 ps(pgc,CP(id)->a,a0);
            CP(id)->pmom->amr_stage_end(CP(id)->pp,CP(id)->a,pgc,s);
        }

        // the finest level: the loads after its stage, on the finest grid at every hull triangle;
        // the force output once per level-0 step
        if(bodies && l==maxlev)
        body_loads(s,t,s==stages-1 && bfin==(1<<maxlev)-1);
    }

    if(bodies && l==maxlev)
    {
        const double dts = p0->dt;
        p0->dt = dt;
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->saveTimeStep(p0,stages-1);
        p0->dt = dts;
        ++bfin;
    }

    if(bodies && l<maxlev)
    for(int nb=0; nb<b6->objects(); ++nb)
    b6->object(nb)->amr_restore(p0,pgc0);

    if(turb)
    for(int id : lev[l])
    {
        pscope3 ps(pgc,CP(id)->a,a0);
        CP(id)->pturb->ktimesave(CP(id)->pp,CP(id)->a,pgc);
        CP(id)->pturb->etimesave(CP(id)->pp,CP(id)->a,pgc);
    }
}

// two steps of level l per step (of size dtc, from the time tc) of level l-1, each followed by the
// steps of the finer levels; then level l restricted onto level l-1
void cfd_amr::advance(int l, double dtc, double tc)
{
    for(int m=0; m<2; ++m)
    {
        const double th0 = 0.5*m;

        fill_level(l,th0);
        if(l<maxlev)
        save_old(l);

        patch_step(l,0.5*dtc,th0,tc+th0*dtc);

        th_end[l] = th0+0.5;
        fill_level(l,th_end[l]);

        if(l<maxlev)
        advance(l+1,0.5*dtc,tc+th0*dtc);
    }

    restrict_sub(l);
}

void cfd_amr::step_sub(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{
    const double t0 = MPI_Wtime();

    for(auto q : P)
    q->pp->count = p->count;

    if(p->F46==2 || p->F46==3)
    {
        double V, dV;
        volume(V,dV);
        vol_step = V;
        if(!vol_ini)
        {
            vol_net = V;
            vol_ini = true;
        }
    }

    // level 0: its step as step_ssp (the volume correction over all levels at the end)
    save_old(0);

    // the bodies: level 0 steps with a predicted copy (the finest level advances them)
    if(bodies)
    {
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->amr_save();
        b6->amr_predict = true;
    }
    bfin = 0;

    mom0->amr_step_begin(p,a,pgc);
    for(int s=0; s<mom0->amr_stages(); ++s)
    {
        mom0->amr_ls_transport(p,a,pgc,s);
        mom0->amr_ls_finish(p,a,pgc,s,false,-1);
        mom0->amr_momentum(p,a,pgc,pvrans,p6dof,s,1);
        mom0->amr_projection(p,a,pgc,s);
        mom0->amr_stage_end(p,a,pgc,s);
    }

    const double t1 = MPI_Wtime();
    tm[4] += t1-t0;

    if(bodies)
    {
        b6->amr_predict = false;
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->amr_restore(p0,pgc0);
    }

    advance(1,p->dt,p->simtime);

    // the coarser levels take the fluxes of the patches at their edges: one composite projection
    pr_sync_project(p,pgc);

    press_step = true;

    if(bodies)
    body_end();

    volume_correct(p,pgc);

    if(adapt && p->G2>0 && p->count%p->G2==0)
    regrid(p,pgc);

    tm[7] += MPI_Wtime()-t0;

    if(p->count%p->P12==0)
    {
        double V, dV;
        volume(V,dV);
        if(p->mpirank==0)
        {
            const std::streamsize pr = cout.precision();
            cout<<"CFD AMR: volume of the level set over the leaf cells "<<setprecision(8)<<V<<endl;
            cout.precision(pr);
        }
    }
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
        cout<<"CFD AMR: subcycling, level-0 step "<<setprecision(3)<<tm[4]<<" s, the patches "<<tm[7]-tm[4]<<" s; pressure iterations (last level solve) "<<pr_it_last
            <<"  res "<<pr_res_last<<", synchronisation "<<pr_it_sync<<" ("<<setprecision(3)<<t_sync<<" s)"<<endl;
    }
}

// --------------------------------------------------------------------- floating bodies (G 7 1)
// stage s of the bodies with the step of level l (from t with dt): the finest level advances the
// bodies (the loads of its last stage), a coarser one its predicted copy (the same loads)
void cfd_amr::body_stage(int l, int s, double t, double dt)
{
    const double dts = p0->dt, ts = p0->simtime;
    p0->dt = dt;
    p0->simtime = t;
    for(int nb=0; nb<b6->objects(); ++nb)
    b6->object(nb)->amr_stage_cfd(p0,a0,pgc0,s);
    p0->dt = dts;
    p0->simtime = ts;
}

// the loads after stage s of the finest level (hull triangles on the finest grid there, else on
// level 0, sixdof_obj_cfd::forces_stl); the force file labels them with the start of the finest step
void cfd_amr::body_loads(int s, double t, bool fin)
{
    // the cells around the patches that copy a patch of the same level: up to date (the loads
    // take the patch values there, finest_grid)
    for(int l=1; l<=maxlev; ++l)
    {
        vector<fspec> fs;
        for(int c=0; c<3; ++c)
        fs.push_back({c,1,[this,c](int g) -> field& { if(c==0) return gfd(g)->u; if(c==1) return gfd(g)->v; return gfd(g)->w; }});
        fs.push_back({CELL,0,[this](int g) -> field& { return gfd(g)->press; }});
        fill(l,fs,true);
    }

    const double ts = p0->simtime;
    p0->simtime = t;
    for(int nb=0; nb<b6->objects(); ++nb)
    b6->object(nb)->hydrodynamic_forces_cfd(p0,a0,pgc0,a0->u,a0->v,a0->w,s,fin);
    p0->simtime = ts;
}

// end of the level-0 step: the hull on level 0 at the final pose, the output of the bodies once per
// level-0 step
void cfd_amr::body_end()
{
    for(int nb=0; nb<b6->objects(); ++nb)
    {
        sixdof_obj_cfd *o = b6->object(nb);
        o->update_position_3D(p0,a0,pgc0,true);

        if(p0->X50==1)
        o->print_vtp(p0,pgc0);

        if(p0->X50==2)
        o->print_stl(p0,pgc0);

        o->print_parameter(p0,pgc0);
    }
    pgc0->gcdf_update(p0,a0);
}

// --------------------------------------------------------------------- composite volume (F 46)
// the volume of the level set over the leaf cells of all grids and its derivative with respect to
// a uniform shift (as picard_lsm::volcalc2)
void cfd_amr::volume(double &V, double &dV)
{
    V = dV = 0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        field &f = gfd(g)->phi;
        for(int qq : lq[g+1])
        {
            int i, j, k;
            cdecode(q,qq,i,j,k);
            const double dv = q->DXN[i+marge]*q->DYN[j+marge]*q->DZN[k+marge];
            const double e = interface_width(q,f,i,j,k);
            V += dv*heaviside(f.V[qq],e);
            dV += dv*heaviside_delta(f.V[qq],e);
        }
    }
    double s[2] = {V,dV};
    MPI_Allreduce(MPI_IN_PLACE,s,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    V = s[0];
    dV = s[1];
}

// the uniform shift of the level set of all grids that restores the volume (F 46 2: the volume at
// the start of the step, F 46 3: the initial volume plus inflow minus outflow), F 47 Newton
// iterations; then density and viscosity of all grids
void cfd_amr::volume_correct(lexer *p, ghostcell *pgc)
{
    if(p->F46!=2 && p->F46!=3)
    return;

    double target = vol_step;
    if(p->F46==3)
    {
        vol_net += p->Qi*p->dt - p->Qo*p->dt;
        target = vol_net;
    }

    for(int n=0; n<p->F47; ++n)
    {
        double V, dV;
        volume(V,dV);
        if(dV<=1.0e-20)
        break;
        const double w = (target-V)/dV;
        for(int g=-1; g<(int)P.size(); ++g)
        {
            field &f = gfd(g)->phi;
            const size_t na = asize(glex(g));
            for(size_t q=0; q<na; ++q)
            f.V[q] += w;
        }
    }

    rovisc_all();
}

// density and viscosity of all grids from the level set (as fluid_update_fsf)
void cfd_amr::rovisc_all()
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        fdm *a = gfd(g);
        for(int i=0; i<q->knox; ++i)
        for(int j=0; j<q->knoy; ++j)
        for(int k=0; k<q->knoz; ++k)
        {
            const double H = heaviside(a->phi(i,j,k),interface_width(q,a->phi,i,j,k));
            a->ro(i,j,k) = q->W1*H + q->W3*(1.0-H);
            a->visc(i,j,k) = q->W2*H + q->W4*(1.0-H);
        }
        if(g>=0)
        margin_rovisc(g);
    }

    pgc0->start4(p0,a0->phi,gcval_phi);
    pgc0->start4(p0,a0->ro,1);
    pgc0->start4(p0,a0->visc,1);
}
