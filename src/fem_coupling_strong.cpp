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


#include"fem_coupling.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"sliceint.h"
#include<cmath>
#include<cstring>
#include<algorithm>
#include<iostream>

typedef fem_solid::Vec3 Vec3;

// Strong coupling of rigid bodies with NHFLOW (coupling option 'strong on').
//
// The fluid loads of a step are probed from the flow that was forced with the body
// velocity of the previous step; the added-mass stabilisation (M + A) a = F + A a_prev
// removes the instability of light bodies, but it delays their response when A is far from
// the true added mass, and NHFLOW needs a large A (5 times the estimate, at most 30 body
// masses, patch 0020). With strong coupling the fluid step is repeated: every pass forces
// the flow with the estimate of the body motion over the step and uses the acceleration of
// the estimate in the stabilisation; the passes stop when the end-of-step velocity of every
// rigid body in the water (translation, and rotation times the radius of gyration) changes by
// less than strong_tol relative to its speed (at least 0.05 m/s). Then the A terms cancel and
// the body motion is that of the converged loads, whatever A; A is the plain estimate. The
// estimates are relaxed with Aitken's factor per body (between 0.01 and 20: the passes of a
// body that the stabilisation holds back converge slowly by themselves, the factor finds
// the fixed point of the nearly linear map in a few passes). The fluid fields, p->DF and the FEM
// state (solid, filters, marks) are saved before the momentum step and restored before every
// further pass; the outputs of the step are written once, after the last pass.
// 'strong test' repeats every step once with unchanged inputs: the results must be those of
// 'strong off' (check that the saved state is complete).

void fem_coupling::sc_fluid(lexer *p, fdm_nhf *d, bool save)
{
    const size_t n3 = size_t(p->imax)*p->jmax*(p->kmax+2);
    const size_t n2 = size_t(p->imax)*p->jmax;
    const size_t nd = size_t(p->imax)*p->jmax*p->kmax;

    double *a3[] = {d->U,d->V,d->W,d->omegaF,d->UH,d->VH,d->WH,d->P,d->RO,d->VISC,d->EV,
                    p->sig,p->sigx,p->sigy,p->sigt,p->sigxx};
    slice *s2[] = {&d->eta,&d->eta_n,&d->WL,&d->detadt,&d->detadt_n,&d->un,&d->vn,&d->dudt,
                   &d->bed,&d->depth,&d->K,&d->fs,&d->Hs,&d->breaking_print,&d->hx,&d->hy};
    sliceint *i2[] = {&d->etaloc,&d->wet_n,&d->breaking,&d->breaklog,&d->bc,&d->nodeval2D};

    const size_t na = sizeof(a3)/sizeof(a3[0]), ns = sizeof(s2)/sizeof(s2[0]), ni = sizeof(i2)/sizeof(i2[0]);
    if(save)
    {
        sc->farr.resize(na+ns+2);
        sc->iarr.resize(ni+6);
    }
    std::vector<std::vector<double>>& F = sc->farr;
    std::vector<std::vector<int>>& I = sc->iarr;

    auto cpd = [&](std::vector<double>& b, double *x, size_t n)
    {
        if(save) b.assign(x,x+n);
        else std::memcpy(x,b.data(),n*sizeof(double));
    };
    auto cpi = [&](std::vector<int>& b, int *x, size_t n)
    {
        if(save) b.assign(x,x+n);
        else std::memcpy(x,b.data(),n*sizeof(int));
    };
    for(size_t q=0; q<na; ++q)
    cpd(F[q],a3[q],n3);
    for(size_t q=0; q<ns; ++q)
    cpd(F[na+q],s2[q]->data(),n2);
    cpd(F[na+ns],p->sigz,n2);
    cpd(F[na+ns+1],p->WL,n2);
    for(size_t q=0; q<ni; ++q)
    cpi(I[q],i2[q]->data(),n2);
    cpi(I[ni],p->DF,nd);
    cpi(I[ni+1],p->wet,n2);
    cpi(I[ni+2],p->wet_n,n2);
    cpi(I[ni+3],p->deep,n2);
    cpi(I[ni+4],p->flag4,nd);
    cpi(I[ni+5],p->flagslice4,n2);
}

void fem_coupling::sc_save()
{
    sc->fs = fs;
    sc->pts = pts;
    sc->surf_version = surf_version;
    sc->fdeb = fdeb; sc->fprb = fprb; sc->fp_bar = fp_bar; sc->fr_bar = fr_bar;
    sc->lim_u = lim_u; sc->lim_vt = lim_vt; sc->rb_lo = rb_lo; sc->rb_hi = rb_hi;
    sc->pnh_bar = pnh_bar; sc->pnh_bar_b = pnh_bar_b; sc->sp_w = sp_w; sc->sp_u = sp_u;
    sc->df_cells = df_cells; sc->df_save = df_save;
    sc->hybrid_ini = hybrid_ini; sc->warned_overtop = warned_overtop;
    sc->printtime = printtime; sc->printcount = printcount; sc->nstep = nstep;
}

void fem_coupling::sc_restore()
{
    fs = sc->fs;
    pts = sc->pts;
    surf_version = sc->surf_version;
    fdeb = sc->fdeb; fprb = sc->fprb; fp_bar = sc->fp_bar; fr_bar = sc->fr_bar;
    lim_u = sc->lim_u; lim_vt = sc->lim_vt; rb_lo = sc->rb_lo; rb_hi = sc->rb_hi;
    pnh_bar = sc->pnh_bar; pnh_bar_b = sc->pnh_bar_b; sp_w = sc->sp_w; sp_u = sc->sp_u;
    df_cells = sc->df_cells; df_save = sc->df_save;
    hybrid_ini = sc->hybrid_ini; warned_overtop = sc->warned_overtop;
    printtime = sc->printtime; printcount = sc->printcount; nstep = sc->nstep;
}

void fem_coupling::nhf_step_begin(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    (void)pgc;
    sc_on = false;
    sc_defer = false;
    sc_pass = 0;
    const int mode = fs.coupling().strong;
    if(mode==0 || !nhflow || !initialised || fs.n_rigid()==0)
    return;

    // rigid bodies in the water (wet at the end of the last step); test mode: every step
    const int nr = fs.n_rigid();
    sc_wet.assign(nr,0);
    bool any = false;
    for(int k=0; k<nr; ++k)
    if(mode==2 || ((size_t)k<lim_vt.size() && lim_vt[k]>0.0))
    {
        sc_wet[k] = 1;
        any = true;
    }
    if(!any)
    return;

    if(sc==nullptr)
    sc = new sc_snapshot;
    sc_fluid(p,d,true);
    sc_save();
    sc_on = true;
    sc_defer = true;
    sc_Vin.assign(nr,Vec3::Zero()); sc_win.assign(nr,Vec3::Zero());
    sc_rV.assign(nr,Vec3::Zero()); sc_rw.assign(nr,Vec3::Zero());
    sc_V0.resize(nr); sc_w0.resize(nr);
    sc_om.assign(nr,1.0);
    for(int k=0; k<nr; ++k)
    {
        sc_V0[k] = fs.rigid(k).V;
        sc_w0[k] = fs.rigid(k).w;
    }
}

bool fem_coupling::nhf_step_redo(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(!sc_on)
    return false;

    const int mode = fs.coupling().strong;
    const int nr = fs.n_rigid();
    ++sc_pass;

    bool done;
    if(mode==2)
    done = (sc_pass>=2);
    else
    {
        // end-of-step velocities of this pass against the estimate the pass was forced with
        done = true;
        if(fs.n_rigid()!=(int)sc_wet.size())
        done = true;        // bodies changed in the step (fragments): accept
        else
        for(int k=0; k<nr; ++k)
        if(sc_wet[k])
        {
            const fem_solid::rigid_body& rb = fs.rigid(k);
            const double rg = rb.M>0.0 ? std::sqrt(std::max(0.0,rb.I0.trace())/(2.0*rb.M)) : 0.0;
            const Vec3 rV = rb.V - (sc_pass==1 ? sc_V0[k] : sc_Vin[k]);
            const Vec3 rw = rb.w - (sc_pass==1 ? sc_w0[k] : sc_win[k]);
            if(sc_pass==1)
            {
                // the first pass is the explicit step: its result is the first estimate
                sc_Vin[k] = rb.V;
                sc_win[k] = rb.w;
                done = false;
                continue;
            }
            const double res = std::sqrt(rV.squaredNorm() + rg*rg*rw.squaredNorm());
            const double ref = std::max(0.05, rb.V.norm() + rg*rb.w.norm());
            if(res > fs.coupling().strong_tol*ref)
            done = false;

            // Aitken relaxation of the estimate
            if(sc_pass>2)
            {
                const Vec3 dV = rV - sc_rV[k], dw = rw - sc_rw[k];
                const double den = dV.squaredNorm() + rg*rg*dw.squaredNorm();
                if(den>0.0)
                sc_om[k] = -sc_om[k]*(sc_rV[k].dot(dV) + rg*rg*sc_rw[k].dot(dw))/den;
                sc_om[k] = std::max(0.01,std::min(20.0,sc_om[k]));
            }
            sc_rV[k] = rV;
            sc_rw[k] = rw;
            sc_Vin[k] += sc_om[k]*rV;
            sc_win[k] += sc_om[k]*rw;
        }
        if(!done && sc_pass>=fs.coupling().strong_maxit)
        {
            done = true;
            ++sc_unconv;
        }
    }

    if(done)
    {
        ++sc_steps;
        sc_passes += sc_pass;
        sc_on = false;
        sc_defer = false;
        step_output(p,pgc);
        return false;
    }

    // next pass: fluid and solid back to the start of the step, the bodies in the water
    // forced with the mean velocity of the estimate over the step and stabilised with its
    // acceleration
    sc_fluid(p,d,false);
    sc_restore();
    if(mode==1)
    for(int k=0; k<nr; ++k)
    if(sc_wet[k])
    {
        const fem_solid::rigid_body& rb = fs.rigid(k);
        const Vec3 Vm = 0.5*(sc_V0[k] + sc_Vin[k]), wm = 0.5*(sc_w0[k] + sc_win[k]);
        for(int i : rb.nodes)
        {
            Vec3 v = Vm + wm.cross(fs.pos(i) - rb.c);
            if(p->j_dir==0) v(1) = 0.0;
            fs.set_vel_mean(i,v);
        }
        fs.set_rigid_accel(k,(sc_Vin[k]-sc_V0[k])/p->dt,(sc_win[k]-sc_w0[k])/p->dt);
    }
    return true;
}
