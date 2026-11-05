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

#include"nhflow_tracer_SI.h"
#include"nhflow_momentum_SI.h"
#include"reef3d_mpmd.h"
#include"reefckpt.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<cmath>
#include<cstring>
#include<iomanip>
#include<iostream>
#include<mpi.h>
#include<sys/stat.h>
#include<cstdlib>

nhflow_tracer_SI::nhflow_tracer_SI(lexer *p, fdm_nhf *d, nhflow_momentum_SI *mom)
                : lx(p),dnx(p),ly(p),dny(p),ar(p),wl0(p),wl1(p),trini(p)
{
    m = mom;
    nt = p->A600>=2 ? 2 : 1;
    K = p->knoz;
    N = p->imax*p->jmax*(p->kmax+2);
    scheme = p->A609;
    dt = 0.0;
    ini = false;

    riv[0] = p->A606_s;
    riv[1] = p->A606_t;
    sea[0] = p->A607_s;
    sea[1] = p->A607_t;
    tth = p->A608;
    kh = p->A604;
    sigt = p->A605_sig>0.0 ? p->A605_sig : 1.0;
    k0 = p->A605_k0;

    // geometry of the faces and cells (static)
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
    lx(ii,jj) = dnx(ii,jj) = ly(ii,jj) = dny(ii,jj) = 0.0;
    ar(ii,jj) = m->carea(p,ii,jj);
    trini(ii,jj) = 0.0;
    wl0(ii,jj) = wl1(ii,jj) = 0.0;
    }

    for(int ii=p->imin; ii<p->imin+p->imax-1; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
    lx(ii,jj) = m->flen(p,1,ii,jj);
    dnx(ii,jj) = m->fdn(p,1,ii,jj);
    }

    if(p->j_dir==1)
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax-1; ++jj)
    {
    ly(ii,jj) = m->flen(p,2,ii,jj);
    dny(ii,jj) = m->fdn(p,2,ii,jj);
    }

    om.assign((size_t)p->imax*p->jmax*(K+1),0.0);
    ca.assign(K,0.0); cb.assign(K,0.0); cc.assign(K,0.0); cr.assign(K,0.0); cx.assign(K,0.0);
    bfx.assign(p->imax*p->jmax,-1);
    bfy.assign(p->imax*p->jmax,-1);
    violc.assign(p->imax*p->jmax,0.0);

    for(int n=0; n<2; ++n)
    {
    Co[n] = nullptr;
    thc[n] = nullptr;
    Cx[n] = nullptr;
    fvx[n] = fvy[n] = nullptr;
    cmn[n].assign(p->imax*p->jmax,0.0);
    cmx[n].assign(p->imax*p->jmax,0.0);
    clipc[n].assign(p->imax*p->jmax,0.0);
    binc[n].assign(p->imax*p->jmax,0.0);
    }

    for(int n=0; n<nt; ++n)
    {
    p->Darray(Co[n],N);
    p->Darray(thc[n],N);
    p->Darray(Cx[n],N);
    p->Darray(fvx[n],N);
    p->Darray(fvy[n],N);

        for(int q=0; q<N; ++q)
        Cx[n][q] = -1.0e30;

        for(int q=0; q<N; ++q)
        thc[n][q] = sea[n];
    }

    p->Darray(tht,N);
    for(int q=0; q<N; ++q)
    tht[q] = 1.0e30;

    p->Darray(mask,N);
    mask_own(p);
}

nhflow_tracer_SI::~nhflow_tracer_SI()
{
    for(int n=0; n<nt; ++n)
    {
    delete [] Co[n];
    delete [] thc[n];
    delete [] Cx[n];
    delete [] fvx[n];
    delete [] fvy[n];
    }
    delete [] tht;
    delete [] mask;
}

double* nhflow_tracer_SI::field(fdm_nhf *d, int mm)
{
    return mm==0 ? d->SALT : d->TEMP;
}

int nhflow_tracer_SI::cidx(lexer *p, int ii, int jj, int kk) const
{
    return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
}

int nhflow_tracer_SI::sidx(lexer *p, int ii, int jj) const
{
    return (ii-p->imin)*p->jmax + (jj-p->jmin);
}

bool nhflow_tracer_SI::fluid(lexer *p, int ii, int jj) const
{
    // the fluid cells of the grid and of its neighbours (partition halo, cells around a block): the mask,
    // exchanged and filled as the tracers (the flags of the lexer are not set deep into the halo)
    if(ii<p->imin || ii>=p->imin+p->imax || jj<p->jmin || jj>=p->jmin+p->jmax)
    return false;

    return mask[cidx(p,ii,jj,0)]>0.5;
}

void nhflow_tracer_SI::mask_own(lexer *p)
{
    for(int q=0; q<N; ++q)
    mask[q] = 0.0;

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    if(p->flagslice4[sidx(p,ii,jj)]>0)
    for(int kk=0; kk<K; ++kk)
    mask[cidx(p,ii,jj,kk)] = 1.0;
}

// ---- initial values (A 601, A 602, A 603), all cells of the grid
void nhflow_tracer_SI::initial(lexer *p, fdm_nhf *d)
{
    const int mg = increment::marge;

    for(int n=0; n<nt; ++n)
    {
        double *C = field(d,n);
        const bool lin = (n==0) ? p->A602==1 : p->A603==1;
        const double c0 = (n==0) ? p->A601_s : p->A601_t;
        const double x1 = n==0 ? p->A602_x1 : p->A603_x1, y1 = n==0 ? p->A602_y1 : p->A603_y1;
        const double x2 = n==0 ? p->A602_x2 : p->A603_x2, y2 = n==0 ? p->A602_y2 : p->A603_y2;
        const double c1 = n==0 ? p->A602_c1 : p->A603_c1, c2 = n==0 ? p->A602_c2 : p->A603_c2;

        for(int q=0; q<N; ++q)
        C[q] = c0;

        if(lin)
        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
        {
            double x,y;

            if(p->geo_curv==1)
            {
            x = p->XC2D[sidx(p,ii,jj)];
            y = p->YC2D[sidx(p,ii,jj)];
            }

            else
            {
            x = p->XP[ii+mg];
            y = p->j_dir==1 ? p->YP[jj+mg] : y1;
            }

            const double dx = x2-x1, dy = y2-y1, L2 = dx*dx + dy*dy;
            double t = L2>0.0 ? ((x-x1)*dx + (y-y1)*dy)/L2 : 0.0;
            t = t<0.0 ? 0.0 : (t>1.0 ? 1.0 : t);

            for(int kk=0; kk<K; ++kk)
            C[cidx(p,ii,jj,kk)] = c1 + t*(c2-c1);
        }
    }

    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    trini(ii,jj) = 1.0;

    ini = true;
}

// boundary value of tracer mm at an inflow face: discharge boundaries (rivers) A 606; water level and
// radiating boundaries (sea) A 607, with the Thatcher-Harleman return from the last outflow value thcv
// over A 608 (thtv: the time of inflow so far)
double nhflow_tracer_SI::bvalue(int mm, int type, double thcv, double thtv) const
{
    if(type==nhflow_momentum_SI::bc_flux)
    return riv[mm];

    if(tth<=0.0 || thtv>=tth)
    return sea[mm];

    return thcv + (sea[mm]-thcv)*0.5*(1.0 - cos(M_PI*thtv/tth));
}

double nhflow_tracer_SI::lim(double a, double b) const
{
    // limited slope from the upwind difference a and the downwind difference b (both along the flow):
    // van Leer (A 609 2) or monotonised central (A 609 3); 0 at an extremum
    if(a*b<=0.0)
    return 0.0;

    if(scheme==3)
    {
    const double s = a>0.0 ? 1.0 : -1.0;
    return s*fmin(fmin(2.0*fabs(a),2.0*fabs(b)),0.5*fabs(a+b));
    }

    return 2.0*a*b/(a+b);
}

// the value of tracer mm on face (ii,jj) (dir 1: i+1/2, dir 2: j+1/2) in layer kk with the transport T
// (positive towards +i / +j), from the values at the start of the substep
double nhflow_tracer_SI::face_value(lexer *p, int mm, int dir, int ii, int jj, int kk, double T)
{
    const double *C = Co[mm];

    // open boundary face: inflow the boundary value, outflow the cell value
    const int b = (dir==1 ? bfx : bfy)[sidx(p,ii,jj)];

    if(b>=0)
    {
        const auto &f = m->bcf[b];

        if(f.s*T>0.0)
        return bval[mm][(size_t)b*K+kk];

        const int ci = cidx(p,f.i,f.j,kk);
        return (cx_on && Cx[mm][ci]>-1.0e29) ? Cx[mm][ci] : C[ci];
    }

    const int di = dir==1 ? 1 : 0, dj = dir==2 ? 1 : 0;
    int ui,uj,dni,dnj,uui,uuj;

    if(T>=0.0)
    {
    ui = ii; uj = jj;
    dni = ii+di; dnj = jj+dj;
    uui = ii-di; uuj = jj-dj;
    }

    else
    {
    ui = ii+di; uj = jj+dj;
    dni = ii; dnj = jj;
    uui = ii+2*di; uuj = jj+2*dj;
    }

    const double cu = C[cidx(p,ui,uj,kk)];

    // a flushed upwind cell: its mixed content
    if(cx_on && Cx[mm][cidx(p,ui,uj,kk)]>-1.0e29)
    return Cx[mm][cidx(p,ui,uj,kk)];

    if(scheme==1)
    return cu;

    const double cd = fluid(p,dni,dnj) ? C[cidx(p,dni,dnj,kk)] : cu;
    const double cuu = fluid(p,uui,uuj) ? C[cidx(p,uui,uuj,kk)] : cu;

    // Courant number of the face for the substep, relative to the content of the upwind cell: the
    // limited correction (1 - nu)/2 (Lax-Wendroff in the smooth part, Sweby 1984), second order in time
    // and space with the forward step (without (1 - nu) the step is anti-diffusive, -u^2 dt/2)
    const double wu = wl0(ui,uj) + sfa*(wl1(ui,uj) - wl0(ui,uj));
    const double L = dir==1 ? lx(ii,jj) : ly(ii,jj);
    double nu = wu>1.0e-12 ? fabs(T)*L*sdt/(ar(ui,uj)*wu) : 1.0;
    nu = nu>1.0 ? 1.0 : nu;

    return cu + 0.5*(1.0-nu)*lim(cu-cuu,cd-cu);
}

int nhflow_tracer_SI::substeps(lexer *p, double rmax)
{
    // the ratio of prepare below 0.9 per substep, at most A 610 substeps
    const double rt = 0.9;
    const int nmax = p->A610>0 ? p->A610 : 1;
    int n = int(ceil(rmax/rt - 1.0e-12));
    return n<1 ? 1 : (n>nmax ? nmax : n);
}

// ---- the step: initial values, omega of the transport, boundary values, the ratio of the substeps
double nhflow_tracer_SI::prepare(lexer *p, fdm_nhf *d)
{
    const int mg = increment::marge;
    dt = p->dt;

    // the first step: the initial values, or the values of the checkpoint (trini set)
    if(!ini)
    {
        bool rs = false;

        for(int ii=0; ii<p->knox && !rs; ++ii)
        for(int jj=0; jj<p->knoy && !rs; ++jj)
        if(trini(ii,jj)>0.5)
        rs = true;

        if(rs)
        ini = true;

        else
        {
        initial(p,d);
        ghosts(p,d);
        }
    }

    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
    wl0(ii,jj) = m->WL0(ii,jj);
    wl1(ii,jj) = d->WL(ii,jj);
    }

    // omega from the layer continuity with the transport of the step
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        const size_t o = (size_t)sidx(p,ii,jj)*(K+1);

        for(int kk=0; kk<=K; ++kk)
        om[o+kk] = 0.0;

        if(p->flagslice4[sidx(p,ii,jj)]<0)
        continue;

        const double dw = (wl1(ii,jj) - wl0(ii,jj))/dt;
        double w = 0.0;

        for(int kk=0; kk<K; ++kk)
        {
            double dv = m->Tx[cidx(p,ii,jj,kk)]*lx(ii,jj) - m->Tx[cidx(p,ii-1,jj,kk)]*lx(ii-1,jj);

            if(p->j_dir==1)
            dv += m->Ty[cidx(p,ii,jj,kk)]*ly(ii,jj) - m->Ty[cidx(p,ii,jj-1,kk)]*ly(ii,jj-1);

            w -= p->DZN[kk+mg]*(dw + dv/ar(ii,jj));
            om[o+kk+1] = w;
        }
    }

    // open faces of the step, their boundary values per layer; the Thatcher-Harleman state of the cell
    for(size_t q=0; q<bfx.size(); ++q)
    bfx[q] = bfy[q] = -1;

    const auto &bcf = m->bcf;

    for(int n=0; n<nt; ++n)
    bval[n].assign(bcf.size()*K,0.0);

    for(size_t b=0; b<bcf.size(); ++b)
    {
        const auto &f = bcf[b];
        const bool xf = (f.cs==1 || f.cs==4);

        (xf ? bfx : bfy)[sidx(p,f.fi,f.fj)] = int(b);

        const double *T = xf ? m->Tx : m->Ty;

        for(int kk=0; kk<K; ++kk)
        {
            const int c = cidx(p,f.i,f.j,kk);
            const bool in = f.s*T[cidx(p,f.fi,f.fj,kk)]>0.0;

            if(in)
            tht[c] += dt;

            else
            tht[c] = 0.0;

            for(int n=0; n<nt; ++n)
            {
                if(!in)
                thc[n][c] = field(d,n)[c];

                bval[n][b*K+kk] = bvalue(n,f.type,thc[n][c],tht[c]);
            }
        }
    }

    // the largest ratio of outflow plus diffusion to content of a layer over the step (wet cells)
    double rmax = 0.0;

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        const int c2 = sidx(p,ii,jj);

        if(p->flagslice4[c2]<0)
        continue;

        const double h = fmin(wl0(ii,jj),wl1(ii,jj));

        // the shoreline columns do not set the substeps (flushed instead)
        if(h<=10.0*p->A544)
        continue;

        double dif = 0.0;

        if(kh>0.0)
        {
        if(m->Hx(ii,jj)>0.0)
        dif += kh*m->Hx(ii,jj)*lx(ii,jj)/dnx(ii,jj);
        if(m->Hx(ii-1,jj)>0.0)
        dif += kh*m->Hx(ii-1,jj)*lx(ii-1,jj)/dnx(ii-1,jj);

            if(p->j_dir==1)
            {
            if(m->Hy(ii,jj)>0.0)
            dif += kh*m->Hy(ii,jj)*ly(ii,jj)/dny(ii,jj);
            if(m->Hy(ii,jj-1)>0.0)
            dif += kh*m->Hy(ii,jj-1)*ly(ii,jj-1)/dny(ii,jj-1);
            }
        }

        for(int kk=0; kk<K; ++kk)
        {
            const double te = m->Tx[cidx(p,ii,jj,kk)]*lx(ii,jj), tw = m->Tx[cidx(p,ii-1,jj,kk)]*lx(ii-1,jj);
            double out = fmax(te,0.0) + fmax(-tw,0.0);
            double in = fmax(-te,0.0) + fmax(tw,0.0);

            if(p->j_dir==1)
            {
            const double tn = m->Ty[cidx(p,ii,jj,kk)]*ly(ii,jj), ts = m->Ty[cidx(p,ii,jj-1,kk)]*ly(ii,jj-1);
            out += fmax(tn,0.0) + fmax(-ts,0.0);
            in += fmax(-tn,0.0) + fmax(ts,0.0);
            }

            // upwind: positive when the outflow of a substep stays below the content; limited (TVD): bounded
            // when inflow plus outflow stay below the content after the horizontal step (>= h - outflow)
            const double f = scheme==1 ? out + dif : in + 2.0*out + dif;
            rmax = fmax(rmax,dt*f/(ar(ii,jj)*h));
        }
    }

    for(int n=0; n<nt; ++n)
    for(size_t q=0; q<clipc[n].size(); ++q)
    clipc[n][q] = binc[n][q] = 0.0;

    for(size_t q=0; q<violc.size(); ++q)
    violc[q] = 0.0;

    return rmax;
}

// ---- pass 1 of substep s: the layers beyond the bound of the substep are flushed, their outflow value Cx
void nhflow_tracer_SI::flushed(lexer *p, fdm_nhf *d, int s, int nsub)
{
    const double dts = dt/double(nsub);
    const double fa = double(s)/double(nsub);
    sdt = dts;
    sfa = fa;
    cx_on = false;

    for(int n=0; n<nt; ++n)
    {
    memcpy(Co[n],field(d,n),N*sizeof(double));

        for(int q=0; q<N; ++q)
        Cx[n][q] = -1.0e30;
    }

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        const int c2 = sidx(p,ii,jj);

        if(p->flagslice4[c2]<0)
        continue;

        const double wa = wl0(ii,jj) + fa*(wl1(ii,jj) - wl0(ii,jj));
        const double V = wa*ar(ii,jj);

        double dif = 0.0;

        if(kh>0.0)
        {
        if(m->Hx(ii,jj)>0.0)
        dif += kh*m->Hx(ii,jj)*lx(ii,jj)/dnx(ii,jj);
        if(m->Hx(ii-1,jj)>0.0)
        dif += kh*m->Hx(ii-1,jj)*lx(ii-1,jj)/dnx(ii-1,jj);

            if(p->j_dir==1)
            {
            if(m->Hy(ii,jj)>0.0)
            dif += kh*m->Hy(ii,jj)*ly(ii,jj)/dny(ii,jj);
            if(m->Hy(ii,jj-1)>0.0)
            dif += kh*m->Hy(ii,jj-1)*ly(ii,jj-1)/dny(ii,jj-1);
            }
        }

        for(int kk=0; kk<K; ++kk)
        {
            // the faces of the layer: transport (outward positive) and length
            const int nf = p->j_dir==1 ? 4 : 2;
            const int fd[4] = {1,1,2,2};
            const int fi[4] = {ii,ii-1,ii,ii};
            const int fj[4] = {jj,jj,jj,jj-1};
            const double sg[4] = {1.0,-1.0,1.0,-1.0};

            double in=0.0, out=0.0, tq[4];

            for(int q=0; q<nf; ++q)
            {
            const double T = (fd[q]==1 ? m->Tx : m->Ty)[cidx(p,fi[q],fj[q],kk)];
            const double L = fd[q]==1 ? lx(fi[q],fj[q]) : ly(fi[q],fj[q]);
            tq[q] = sg[q]*T*L;
            out += fmax(tq[q],0.0);
            in += fmax(-tq[q],0.0);
            }

            const double f = scheme==1 ? out + dif : in + 2.0*out + dif;

            if(V>1.0e-12 && dts*f<=0.9*V)
            continue;

            // flushed: the content mixed with the inflow of the substep
            for(int n=0; n<nt; ++n)
            {
                const int c = cidx(p,ii,jj,kk);
                double num = V*Co[n][c];

                for(int q=0; q<nf; ++q)
                if(tq[q]<0.0)
                {
                const double T = (fd[q]==1 ? m->Tx : m->Ty)[cidx(p,fi[q],fj[q],kk)];
                num -= dts*tq[q]*face_value(p,n,fd[q],fi[q],fj[q],kk,T);
                }

                Cx[n][c] = (V + dts*in)>1.0e-30 ? num/(V + dts*in) : Co[n][c];
            }
        }
    }
}

// ---- pass 2 of substep s of nsub: horizontal explicit, vertical implicit
void nhflow_tracer_SI::substep(lexer *p, fdm_nhf *d, int s, int nsub)
{
    const int mg = increment::marge;
    const double dts = dt/double(nsub);
    const double fa = double(s)/double(nsub), fb = double(s+1)/double(nsub);
    sdt = dts;
    sfa = fa;
    cx_on = true;

    // the range of every fluid column (bounds) and the face values of the faces of the own cells, once
    auto own = [&](int a, int b) { return a>=0 && a<p->knox && b>=0 && b<p->knoy && p->flagslice4[sidx(p,a,b)]>0; };

    for(int n=0; n<nt; ++n)
    {
        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
        if(fluid(p,ii,jj))
        {
            double lo = Co[n][cidx(p,ii,jj,0)], hi = lo;
            for(int kk=1; kk<K; ++kk)
            {
            lo = fmin(lo,Co[n][cidx(p,ii,jj,kk)]);
            hi = fmax(hi,Co[n][cidx(p,ii,jj,kk)]);
            }
            cmn[n][sidx(p,ii,jj)] = lo;
            cmx[n][sidx(p,ii,jj)] = hi;
        }

        for(int ii=-1; ii<p->knox; ++ii)
        for(int jj=0; jj<p->knoy; ++jj)
        if(own(ii,jj) || own(ii+1,jj))
        for(int kk=0; kk<K; ++kk)
        fvx[n][cidx(p,ii,jj,kk)] = face_value(p,n,1,ii,jj,kk,m->Tx[cidx(p,ii,jj,kk)]);

        if(p->j_dir==1)
        for(int ii=0; ii<p->knox; ++ii)
        for(int jj=-1; jj<p->knoy; ++jj)
        if(own(ii,jj) || own(ii,jj+1))
        for(int kk=0; kk<K; ++kk)
        fvy[n][cidx(p,ii,jj,kk)] = face_value(p,n,2,ii,jj,kk,m->Ty[cidx(p,ii,jj,kk)]);
    }

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        const int c2 = sidx(p,ii,jj);

        if(p->flagslice4[c2]<0)
        continue;

        const double wa = wl0(ii,jj) + fa*(wl1(ii,jj) - wl0(ii,jj));
        const double wb = wl0(ii,jj) + fb*(wl1(ii,jj) - wl0(ii,jj));

        // a column (nearly) empty at the end of the substep (below 1e-6 m: the round-off of the transport
        // of a through-flow film exceeds its content): the flushed value (its inflow), else its value;
        // the content is negligible
        if(wb<1.0e-6)
        {
            for(int n=0; n<nt; ++n)
            for(int kk=0; kk<K; ++kk)
            {
            const int c = cidx(p,ii,jj,kk);
            if(Cx[n][c]>-1.0e29)
            field(d,n)[c] = Cx[n][c];
            }
        continue;
        }

        const double A = ar(ii,jj);
        const size_t o = (size_t)c2*(K+1);

        // the open faces of the cell (boundary inflow)
        const int bx[2] = {bfx[c2], bfx[sidx(p,ii-1,jj)]};
        const int by[2] = {p->j_dir==1 ? bfy[c2] : -1, p->j_dir==1 ? bfy[sidx(p,ii,jj-1)] : -1};

        // vertical operator of the column: implicit upwind advection with omega, diffusion
        for(int kk=0; kk<K; ++kk)
        {
        ca[kk] = 0.0;
        cc[kk] = 0.0;
        cb[kk] = p->DZN[kk+mg]*wb;
        }

        for(int kk=0; kk<K-1; ++kk)
        {
            const double w = om[o+kk+1];
            const double wp = fmax(w,0.0), wm = fmax(-w,0.0);
            const double ev = 0.5*(d->EV[cidx(p,ii,jj,kk)] + d->EV[cidx(p,ii,jj,kk+1)]);
            const double kv = ev/sigt + k0;
            const double dz = wb*(p->ZP[kk+1+mg] - p->ZP[kk+mg]);
            const double D = dz>0.0 ? kv/dz : 0.0;

            cb[kk] += dts*(wp + D);
            cc[kk] = -dts*(wm + D);
            cb[kk+1] += dts*(wm + D);
            ca[kk+1] = -dts*(wp + D);
        }

        // the remainder of omega at the surface carries the value of the top layer (not below half the
        // volume of the layer: round-off of a thin column)
        cb[K-1] = fmax(cb[K-1] + dts*om[o+K],0.5*p->DZN[K-1+mg]*wb);

        for(int n=0; n<nt; ++n)
        {
            const double *C = Co[n];
            double *Cn = field(d,n);

            // the inflow values of the column (bounds)
            double ilo = 1.0e30, ihi = -1.0e30;

            for(int kk=0; kk<K; ++kk)
            {
                const int c = cidx(p,ii,jj,kk);
                const double dz = p->DZN[kk+mg];

                // horizontal flux out of the cell (per unit sigma)
                double F=0.0;

                auto face = [&](int dir, int fi, int fj, double sg)
                {
                    const double T = (dir==1 ? m->Tx : m->Ty)[cidx(p,fi,fj,kk)];
                    const double v = (dir==1 ? fvx : fvy)[n][cidx(p,fi,fj,kk)];
                    F += sg*T*v*(dir==1 ? lx(fi,fj) : ly(fi,fj));

                    if(sg*T<0.0)
                    {
                    ilo = fmin(ilo,v);
                    ihi = fmax(ihi,v);
                    }
                };

                face(1,ii,jj,1.0);
                face(1,ii-1,jj,-1.0);

                if(p->j_dir==1)
                {
                face(2,ii,jj,1.0);
                face(2,ii,jj-1,-1.0);
                }

                // horizontal diffusion through the open interior faces
                if(kh>0.0)
                {
                    if(m->Hx(ii,jj)>0.0 && bx[0]<0)
                    F -= kh*m->Hx(ii,jj)*(C[cidx(p,ii+1,jj,kk)] - C[c])/dnx(ii,jj)*lx(ii,jj);
                    if(m->Hx(ii-1,jj)>0.0 && bx[1]<0)
                    F -= kh*m->Hx(ii-1,jj)*(C[cidx(p,ii-1,jj,kk)] - C[c])/dnx(ii-1,jj)*lx(ii-1,jj);

                    if(p->j_dir==1)
                    {
                    if(m->Hy(ii,jj)>0.0 && by[0]<0)
                    F -= kh*m->Hy(ii,jj)*(C[cidx(p,ii,jj+1,kk)] - C[c])/dny(ii,jj)*ly(ii,jj);
                    if(m->Hy(ii,jj-1)>0.0 && by[1]<0)
                    F -= kh*m->Hy(ii,jj-1)*(C[cidx(p,ii,jj-1,kk)] - C[c])/dny(ii,jj-1)*ly(ii,jj-1);
                    }
                }

                // tracer over the open boundary faces of the cell, inward positive
                for(int q=0; q<4; ++q)
                {
                    const int b = q<2 ? bx[q] : by[q-2];

                    if(b<0)
                    continue;

                    const auto &f = m->bcf[b];
                    const bool xf = (f.cs==1 || f.cs==4);
                    const double Tf = (xf ? m->Tx : m->Ty)[cidx(p,f.fi,f.fj,kk)];
                    const double L = xf ? lx(f.fi,f.fj) : ly(f.fi,f.fj);

                    binc[n][c2] += dts*f.s*Tf*(xf ? fvx : fvy)[n][cidx(p,f.fi,f.fj,kk)]*L*dz;
                }

                cr[kk] = dz*(wa*C[c] - dts*F/A);
            }

            // tridiagonal solve (Thomas)
            for(int kk=1; kk<K; ++kk)
            {
                const double w = ca[kk]/cb[kk-1];
                cb[kk] -= w*cc[kk-1];
                cr[kk] -= w*cr[kk-1];
            }

            cx[K-1] = cr[K-1]/cb[K-1];

            for(int kk=K-2; kk>=0; --kk)
            cx[kk] = (cr[kk] - cc[kk]*cx[kk+1])/cb[kk];

            // bounds: the values of the column and of its fluid neighbour columns (the vertical step mixes the
            // layers) and the inflow values of its faces at the start of the substep; the new values are a
            // convex combination of them (the horizontal step is bounded by the substeps or flushed, the column
            // operator is an M-matrix with the layer volumes of the horizontal step as row sums)
            double lo = ilo, hi = ihi;
            {
                const int nb[5][2] = {{0,0},{1,0},{-1,0},{0,1},{0,-1}};

                for(int q=0; q<(p->j_dir==1 ? 5 : 3); ++q)
                if(q==0 || fluid(p,ii+nb[q][0],jj+nb[q][1]))
                {
                const int s2 = sidx(p,ii+nb[q][0],jj+nb[q][1]);
                lo = fmin(lo,cmn[n][s2]);
                hi = fmax(hi,cmx[n][s2]);
                }
            }

            for(int kk=0; kk<K; ++kk)
            {
                const int c = cidx(p,ii,jj,kk);
                const double tol = 1.0e-9*(1.0 + fabs(hi) + fabs(lo));
                double v = cx[kk];

                if(v>hi+tol || v<lo-tol)
                {
                violc[c2] += 1.0;
                const double vc = fmin(fmax(v,lo),hi);
                clipc[n][c2] += (vc - v)*p->DZN[kk+mg]*wb*A;
                v = vc;
                }

                Cn[c] = v;
            }

            // the column operator again for the next tracer (the elimination changed cb)
            if(n+1<nt)
            {
                for(int kk=0; kk<K; ++kk)
                cb[kk] = p->DZN[kk+mg]*wb;

                for(int kk=0; kk<K-1; ++kk)
                {
                    const double w = om[o+kk+1];
                    const double wp = fmax(w,0.0), wm = fmax(-w,0.0);
                    const double ev = 0.5*(d->EV[cidx(p,ii,jj,kk)] + d->EV[cidx(p,ii,jj,kk+1)]);
                    const double kv = ev/sigt + k0;
                    const double dz = wb*(p->ZP[kk+1+mg] - p->ZP[kk+mg]);
                    const double D = dz>0.0 ? kv/dz : 0.0;

                    cb[kk] += dts*(wp + D);
                    cb[kk+1] += dts*(wm + D);
                }

                cb[K-1] = fmax(cb[K-1] + dts*om[o+K],0.5*p->DZN[K-1+mg]*wb);
            }
        }
    }
}

// ---- ghost values for the output
void nhflow_tracer_SI::ghosts(lexer *p, fdm_nhf *d)
{
    const int off[8][2] = {{1,0},{-1,0},{0,1},{0,-1},{1,1},{1,-1},{-1,1},{-1,-1}};

    for(int n=0; n<nt; ++n)
    {
        double *C = field(d,n);

        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
        if(!fluid(p,ii,jj))
        for(int q=0; q<8; ++q)
        {
            const int a = ii+off[q][0], b = jj+off[q][1];

            if(p->j_dir==0 && off[q][1]!=0)
            continue;

            if(fluid(p,a,b))
            {
            for(int kk=0; kk<K; ++kk)
            C[cidx(p,ii,jj,kk)] = C[cidx(p,a,b,kk)];
            break;
            }
        }

        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
        {
        C[cidx(p,ii,jj,-1)] = C[cidx(p,ii,jj,0)];
        C[cidx(p,ii,jj,K)] = C[cidx(p,ii,jj,K-1)];
        }
    }
}

// ---- sums over the own cells
void nhflow_tracer_SI::stats_zero(stats &s)
{
    for(int n=0; n<2; ++n)
    {
    s.mass[n] = s.bin[n] = s.clip[n] = 0.0;
    s.cmin[n] = 1.0e30;
    s.cmax[n] = -1.0e30;
    }
    s.nviol = 0.0;
}

void nhflow_tracer_SI::sums(lexer *p, fdm_nhf *d, const std::function<bool(int,int)> &own, stats &s)
{
    const int mg = increment::marge;

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        const int c2 = sidx(p,ii,jj);

        if(p->flagslice4[c2]<0 || !own(ii,jj))
        continue;

        s.nviol += violc[c2];

        for(int n=0; n<nt; ++n)
        {
            const double *C = field(d,n);

            for(int kk=0; kk<K; ++kk)
            s.mass[n] += p->DZN[kk+mg]*d->WL(ii,jj)*C[cidx(p,ii,jj,kk)]*ar(ii,jj);

            s.bin[n] += binc[n][c2];
            s.clip[n] += clipc[n][c2];

            if(d->WL(ii,jj)>p->A544)
            for(int kk=0; kk<K; ++kk)
            {
            s.cmin[n] = fmin(s.cmin[n],C[cidx(p,ii,jj,kk)]);
            s.cmax[n] = fmax(s.cmax[n],C[cidx(p,ii,jj,kk)]);
            }
        }
    }
}

void nhflow_tracer_SI::reduce(lexer *p, stats &s)
{
    double a[7] = {s.mass[0],s.mass[1],s.bin[0],s.bin[1],s.clip[0],s.clip[1],s.nviol}, b[7];
    MPI_Allreduce(a,b,7,MPI_DOUBLE,MPI_SUM,reef3d_comm());
    s.mass[0] = b[0]; s.mass[1] = b[1]; s.bin[0] = b[2]; s.bin[1] = b[3]; s.clip[0] = b[4]; s.clip[1] = b[5]; s.nviol = b[6];

    double lo[2] = {s.cmin[0],s.cmin[1]}, hi[2] = {s.cmax[0],s.cmax[1]}, r[2];
    MPI_Allreduce(lo,r,2,MPI_DOUBLE,MPI_MIN,reef3d_comm());
    s.cmin[0] = r[0]; s.cmin[1] = r[1];
    MPI_Allreduce(hi,r,2,MPI_DOUBLE,MPI_MAX,reef3d_comm());
    s.cmax[0] = r[0]; s.cmax[1] = r[1];
}

void nhflow_tracer_SI::write_log(lexer *p, const stats &s, int nsub, double rmax)
{
    if(p->mpirank!=0)
    return;

    if(!tlog.is_open())
    {
        mkdir("./REEF3D_NHFLOW_Log",0777);

        if(reefckpt_open(tlog,"./REEF3D_NHFLOW_Log/REEF3D_NHFLOW_tracer.dat",p->I45))   // I 45: continued
        {
            tlog<<"# NHFLOW semi-implicit tracers (A 600): per tracer (salinity, temperature) the mass sum(C V) [unit m3],"<<endl;
            tlog<<"# the net inflow over the open boundaries and the clipped mass of the step, min and max of the wet cells"<<endl;
            tlog<<"# count  simtime  dt  substeps  ratio  bound_violations";
            for(int n=0; n<nt; ++n)
            {
            const char *nm = n==0 ? "S" : "T";
            tlog<<"  mass_"<<nm<<"  inflow_"<<nm<<"  clipped_"<<nm<<"  min_"<<nm<<"  max_"<<nm;
            }
            tlog<<endl;
        }
    }

    tlog<<p->count<<" "<<setprecision(10)<<p->simtime+p->dt<<" "<<p->dt<<" "<<nsub<<" "<<setprecision(6)<<rmax<<" "<<s.nviol;

    for(int n=0; n<nt; ++n)
    tlog<<" "<<setprecision(16)<<s.mass[n]<<" "<<s.bin[n]<<" "<<s.clip[n]<<" "<<setprecision(10)<<s.cmin[n]<<" "<<s.cmax[n];

    tlog<<"\n"<<flush;

    if(p->count%p->P12==0)
    {
    cout<<"tracers: substeps "<<nsub;
    for(int n=0; n<nt; ++n)
    cout<<(n==0?"  S ":"  T ")<<setprecision(5)<<s.cmin[n]<<" .. "<<s.cmax[n];
    if(s.nviol>0.0)
    cout<<"  clipped cells "<<s.nviol;
    cout<<endl;
    }
}

// ---- gauges P 51: the nearest own cell centre
void nhflow_tracer_SI::gauge_find(lexer *p, const std::function<bool(int,int)> &own, vector<double> &dist, vector<int> &cell)
{
    const int mg = increment::marge;
    const int ng = p->P51;

    dist.assign(ng,1.0e30);
    cell.assign(ng,-1);

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        const int c2 = sidx(p,ii,jj);

        if(p->flagslice4[c2]<0 || !own(ii,jj))
        continue;

        const double xc = p->geo_curv==1 ? p->XC2D[c2] : p->XP[ii+mg];

        for(int g=0; g<ng; ++g)
        {
            const double yc = p->geo_curv==1 ? p->YC2D[c2] : (p->j_dir==1 ? p->YP[jj+mg] : p->P51_y[g]);
            const double r = sqrt((xc-p->P51_x[g])*(xc-p->P51_x[g]) + (yc-p->P51_y[g])*(yc-p->P51_y[g]));

            if(r<dist[g])
            {
            dist[g] = r;
            cell[g] = c2;
            }
        }
    }
}

void nhflow_tracer_SI::gauge_values(lexer *p, fdm_nhf *d, const vector<int> &cell, vector<double> &val)
{
    const int mg = increment::marge;
    const int ng = p->P51;

    val.assign((size_t)ng*3*nt,0.0);

    for(int g=0; g<ng; ++g)
    if(cell[g]>=0)
    {
        const int ii = cell[g]/p->jmax + p->imin;
        const int jj = cell[g]%p->jmax + p->jmin;

        for(int n=0; n<nt; ++n)
        {
            const double *C = field(d,n);
            double cm = 0.0;

            for(int kk=0; kk<K; ++kk)
            cm += p->DZN[kk+mg]*C[cidx(p,ii,jj,kk)];

            val[(size_t)g*3*nt + 3*n] = cm;
            val[(size_t)g*3*nt + 3*n + 1] = C[cidx(p,ii,jj,0)];
            val[(size_t)g*3*nt + 3*n + 2] = C[cidx(p,ii,jj,K-1)];
        }
    }
}

// the values of the nearest cell over all ranks (equal distances: their mean)
void nhflow_tracer_SI::gauge_combine(lexer *p, vector<double> &dist, vector<double> &val)
{
    const int ng = p->P51;
    const size_t nv = ng>0 ? val.size()/ng : 0;

    vector<double> gmin(ng);
    MPI_Allreduce(dist.data(),gmin.data(),ng,MPI_DOUBLE,MPI_MIN,reef3d_comm());

    vector<double> loc(val.size()+ng,0.0), all(val.size()+ng);

    for(int g=0; g<ng; ++g)
    if(dist[g]<=gmin[g] && gmin[g]<1.0e29)
    {
        for(size_t q=0; q<nv; ++q)
        loc[g*nv+q] = val[g*nv+q];
        loc[val.size()+g] = 1.0;
    }

    MPI_Allreduce(loc.data(),all.data(),int(loc.size()),MPI_DOUBLE,MPI_SUM,reef3d_comm());

    for(int g=0; g<ng; ++g)
    for(size_t q=0; q<nv; ++q)
    val[g*nv+q] = all[val.size()+g]>0.0 ? all[g*nv+q]/all[val.size()+g] : 0.0;

    dist = gmin;
}

void nhflow_tracer_SI::gauge_write(lexer *p, const vector<double> &val)
{
    static ofstream out;

    if(p->mpirank!=0 || p->P51<=0)
    return;

    const int nt = p->A600>=2 ? 2 : 1;

    if(!out.is_open())
    {
        mkdir("./REEF3D_NHFLOW_Log",0777);

        if(reefckpt_open(out,"./REEF3D_NHFLOW_Log/REEF3D_NHFLOW_tracer_gauges.dat",p->I45))
        {
            out<<"# NHFLOW tracers at the gauges P 51 (nearest cell centre): per gauge and tracer the depth mean, bottom and top layer"<<endl;
            out<<"# simtime";
            for(int g=0; g<p->P51; ++g)
            for(int n=0; n<nt; ++n)
            {
            const char *nm = n==0 ? "S" : "T";
            out<<" \t "<<nm<<"_mean("<<p->P51_x[g]<<","<<p->P51_y[g]<<") \t "<<nm<<"_bottom \t "<<nm<<"_top";
            }
            out<<endl;
        }
    }

    out<<setprecision(10)<<p->simtime;
    for(size_t q=0; q<val.size(); ++q)
    out<<" \t "<<setprecision(10)<<val[q];
    out<<"\n"<<flush;
}

// ---- checkpoints
void nhflow_tracer_SI::ckpt_fields(const std::function<void(const char*, char, double*)> &fd)
{
    fd("trc_thc0",'C',thc[0]);
    if(nt>1)
    fd("trc_thc1",'C',thc[1]);
    fd("trc_tht",'C',tht);
    fd("trc_ini",'S',trini.V);
}
