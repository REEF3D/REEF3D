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

#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_curv.h"
#include<cmath>
#include<unordered_map>

// ---------------------------------------------------------------------------
// A 549 1: semi-Lagrangian (Eulerian-Lagrangian) advection of the momentum along the sigma layers in the predictor of
// the hydrostatic semi-implicit NHFLOW, in the velocity form of Casulli (1990) and Casulli and Cheng (1992):
//
//     h^n u^{n+1} = h^n u(x_d) + dt F^n - dt g h^n grad(eta)^{n+theta} + (column: vertical diffusion, advection, bed)
//     x_d = x - int_t^{t+dt} u dt   (the foot of the characteristic in the layer, backwards along the flow)
//
// The predictor (si_stage_rhs, one stage) is UH^n + h^n (u(x_d) - u^n) + dt F^n (= h^n u(x_d) + dt F^n except in the
// desingularised shallow cells), F^n the other explicit terms (sources, horizontal diffusion). The face discharges
// of the eta system take it as h^n u (as Casulli: H_f^n times the face velocity), so the continuity is that of the
// flux form. The column operator takes the vertical advection in the advective form omega du/dsigma (col_coef, upwind,
// implicit), and cell_update the momentum of the new depth, UH^{n+1} = h^{n+1} u^{n+1}. The momentum is not conserved
// (as every semi-Lagrangian scheme), the volume and the tracers are.
//
// The foot is found in index space (i, j) of the layer by the semi-analytical tracking of Pollock (1988): in every
// cell the index velocity of each direction varies linearly between its faces, the time to the next face is exact;
// the faces take the mean normal velocity of their two cells (u.n L/A per cell area, the Courant number of the face;
// on the curvilinear grid with the face normal and length).
// A face to a dry or solid cell has zero velocity: the track never crosses it (wetting and drying, walls). An open
// boundary face takes the velocity of its interior cell; a track that leaves through it takes the inflow velocities
// of the ghost cell (bc_face ug, vg, wg). The track stays within two cells of the halo (multi-block grids, MPI
// partitions: their filled cells), the time step limits its length (Courant number A549_c, nhflow_momentum_SI::adv_cmin).
//
// The velocity at the foot is the bilinear interpolation (index space) of the cell values of t^n of the layer over
// the wet cells around it (a convex combination: bounded, no new extrema).
// ---------------------------------------------------------------------------

namespace
{
// time to leave [0,1] from x with the velocity a + b x (inf: no exit), the exit side in s (+1 / -1)
inline double pollock_time(double x, double a, double b, int &s)
{
    const double v = a + b*x;
    s = 0;

    if(fabs(v)<1.0e-14)
    return 1.0e30;

    const double xe = v>0.0 ? 1.0 : 0.0;
    s = v>0.0 ? 1 : -1;
    const double ve = a + b*xe;

    if(fabs(b)<1.0e-10*(fabs(a)+fabs(b)+1.0e-30))
    return MAX((xe - x)/v, 0.0);

    if(ve*v<=0.0)
    {
    s = 0;
    return 1.0e30;
    }

    return MAX(log(ve/v)/b, 0.0);
}

inline double pollock_move(double x, double a, double b, double t)
{
    const double v = a + b*x;

    if(fabs(b*t)<1.0e-8)
    return x + v*t*(1.0 + 0.5*b*t);

    return x + v*(exp(b*t) - 1.0)/b;
}
}

bool nhflow_momentum_SI::elm_fluid(lexer *p, int ii, int jj)
{
    if(ii<p->imin+1 || ii>p->imin+p->imax-2)
    return false;

    if(p->j_dir==1 && (jj<p->jmin+1 || jj>p->jmin+p->jmax-2))
    return false;

    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);

    return p->flagslice4[q]>0 && p->wet[q]==1 && WL0(ii,jj)>p->A544;
}

void nhflow_momentum_SI::elm_departure(lexer *p, fdm_nhf *d)
{
    const double dt = p->dt;
    const int ii0=i, jj0=j, kk0=k;

    auto c3 = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    auto c2 = [&](int a, int b) { return (a-p->imin)*p->jmax + (b-p->jmin); };

    // the open boundary faces: (dir, face i, face j) -> bcf
    std::unordered_map<long long,int> obf;
    auto fkey = [&](int dir, int a, int b) { return ((long long)dir*1000003LL + (long long)(a+16))*1000003LL + (long long)(b+16); };
    for(int n=0; n<(int)bcf.size(); ++n)
    {
        const bc_face &f = bcf[n];
        const int dir = (f.cs==1 || f.cs==4) ? 1 : 2;
        obf[fkey(dir,f.fi,f.fj)] = n;
    }

    // normal velocity of face a+1/2 (dir 1) or b+1/2 (dir 2) of cell (a,b), layer c, times its length over the area of
    // cell (ca,cb): the index velocity of the face in that cell; 0 to a dry or solid cell, an open boundary face the
    // velocity of its interior cell
    auto face_iv = [&](int dir, int a, int b, int c, int ca, int cb) -> double
    {
        const int a2 = dir==1 ? a+1 : a;
        const int b2 = dir==2 ? b+1 : b;
        const bool f1 = elm_fluid(p,a,b);
        const bool f2 = elm_fluid(p,a2,b2);
        double u, v;

        if(f1 && f2)
        {
        u = 0.5*(d->U[c3(a,b,c)] + d->U[c3(a2,b2,c)]);
        v = 0.5*(d->V[c3(a,b,c)] + d->V[c3(a2,b2,c)]);
        }

        else if((f1 || f2) && obf.count(fkey(dir,a,b))>0)
        {
        const int q = f1 ? c3(a,b,c) : c3(a2,b2,c);
        u = d->U[q];
        v = d->V[q];
        }

        else
        return 0.0;

        if(p->geo_curv==1)
        {
        double nx, ny, len;
        nhflow_curv::normal(p,dir,a,b,nx,ny,len);
        return (u*nx + v*ny)*len/p->AREA2D[c2(ca,cb)];
        }

        return dir==1 ? u/p->DXN[ca+marge] : v*p->y_dir/p->DYN[cb+marge];
    };

    LOOP
    {
        double &ud = ELU[IJK], &vd = ELV[IJK], &wd = ELW[IJK];
        ud = d->U[IJK];
        vd = d->V[IJK];
        wd = d->W[IJK];

        if(!(p->flagslice4[IJ]>0 && p->wet[IJ]==1 && WL0(i,j)>p->A544))
        continue;

        // ---- the track backwards from the cell centre
        int ci=i, cj=j;
        const int ck=k;
        double x=0.5, y=0.5, tau=dt;
        int inflow=-1;

        for(int it=0; it<64 && tau>0.0; ++it)
        {
            // the index velocities of the faces of the cell, backwards
            const double vw = -face_iv(1,ci-1,cj,ck,ci,cj);
            const double ve = -face_iv(1,ci,cj,ck,ci,cj);
            double vs=0.0, vn=0.0;
            if(p->j_dir==1)
            {
            vs = -face_iv(2,ci,cj-1,ck,ci,cj);
            vn = -face_iv(2,ci,cj,ck,ci,cj);
            }

            int sx, sy=0;
            const double tx = pollock_time(x,vw,ve-vw,sx);
            const double ty = p->j_dir==1 ? pollock_time(y,vs,vn-vs,sy) : 1.0e30;
            const double tm = MIN(tx,ty);

            if(tm>=tau)
            {
                x = pollock_move(x,vw,ve-vw,tau);
                if(p->j_dir==1)
                y = pollock_move(y,vs,vn-vs,tau);
                tau = 0.0;
                break;
            }

            x = pollock_move(x,vw,ve-vw,tm);
            if(p->j_dir==1)
            y = pollock_move(y,vs,vn-vs,tm);
            tau -= tm;

            // the face reached
            if(tm==tx)
            {
                x = sx>0 ? 1.0 : 0.0;
                const int ni = ci + sx;

                if(elm_fluid(p,ni,cj))
                {
                ci = ni;
                x = sx>0 ? 0.0 : 1.0;
                continue;
                }

                auto it2 = obf.find(fkey(1,sx>0?ci:ci-1,cj));
                if(it2!=obf.end())
                inflow = it2->second;
                break;
            }

            if(tm==ty && tm<1.0e29)
            {
                y = sy>0 ? 1.0 : 0.0;
                const int nj = cj + sy;

                if(elm_fluid(p,ci,nj))
                {
                cj = nj;
                y = sy>0 ? 0.0 : 1.0;
                continue;
                }

                auto it2 = obf.find(fkey(2,ci,sy>0?cj:cj-1));
                if(it2!=obf.end())
                inflow = it2->second;
                break;
            }

            break;
        }

        x = MAX(0.0,MIN(1.0,x));
        y = MAX(0.0,MIN(1.0,y));

        // ---- the velocity at the foot
        if(inflow>=0)
        {
            const bc_face &f = bcf[inflow];
            if((int)f.ug.size()>ck)
            {
            ud = f.ug[ck];
            vd = f.vg[ck];
            wd = f.wg[ck];
            }
            continue;
        }

        const int ia = x>=0.5 ? ci : ci-1;
        const double wx = x>=0.5 ? x-0.5 : x+0.5;
        const int ja = p->j_dir==1 ? (y>=0.5 ? cj : cj-1) : cj;
        const double wy = p->j_dir==1 ? (y>=0.5 ? y-0.5 : y+0.5) : 0.0;

        double su=0.0, sv=0.0, sw=0.0, sm=0.0;

        for(int a=0; a<2; ++a)
        for(int b=0; b<(p->j_dir==1?2:1); ++b)
        {
            const int ni = ia+a, nj = ja+b;
            const double wh = (a==1 ? wx : 1.0-wx)*(p->j_dir==1 ? (b==1 ? wy : 1.0-wy) : 1.0);

            if(wh<=0.0 || !elm_fluid(p,ni,nj))
            continue;

            const int q = c3(ni,nj,ck);
            su += wh*d->U[q];
            sv += wh*d->V[q];
            sw += wh*d->W[q];
            sm += wh;
        }

        if(sm>1.0e-12)
        {
        ud = su/sm;
        vd = sv/sm;
        wd = sw/sm;
        }

        else
        {
        const int q = c3(ci,cj,ck);
        ud = d->U[q];
        vd = d->V[q];
        wd = d->W[q];
        }
    }

    i=ii0; j=jj0; k=kk0;
}

// the advective time step limit of cell (i,j) as cmin (dt = 2 cfl cmin, nhflow_timestep_SI, nhflow_amr::dt_si_cmin);
// A 549 0: the Courant number of the flow plus the speed of the first internal mode ci (A 620); A 549 1
// (semi-Lagrangian advection): the flow alone up to the Courant number A549_c, at most 1 with a turbulence model (its
// convection is explicit upwind), and the internal speed (the baroclinic term is explicit) with 2 cfl
double nhflow_momentum_SI::adv_cmin(lexer *p, double u, double v, double ci, int ii, int jj, double cfl)
{
    double sigma;

    if(elm_mode(p)==0)
    {
        if(p->geo_curv==1)
        {
        const double sp = sqrt(u*u + v*v) + ci;
        sigma = sp/nhflow_curv::width1(p,ii,jj) + sp/nhflow_curv::width2(p,ii,jj);
        }

        else
        {
        sigma = (fabs(u) + ci)/p->DXN[ii+marge];

        if(p->j_dir==1)
        sigma += (fabs(v) + ci)/p->DYN[jj+marge];
        }

        return sigma>1.0e-20 ? 1.0/sigma : 1.0e20;
    }

    double sa, sc;

    if(p->geo_curv==1)
    {
    const double sp = sqrt(u*u + v*v);
    const double r = 1.0/nhflow_curv::width1(p,ii,jj) + 1.0/nhflow_curv::width2(p,ii,jj);
    sa = sp*r;
    sc = ci*r;
    }

    else
    {
    const double r = 1.0/p->DXN[ii+marge] + (p->j_dir==1 ? 1.0/p->DYN[jj+marge] : 0.0);
    sa = (fabs(u)/p->DXN[ii+marge] + (p->j_dir==1 ? fabs(v)/p->DYN[jj+marge] : 0.0));
    sc = ci*r;
    }

    const double cmax = p->A560>0 ? MIN(elm_c(p),1.0) : elm_c(p);
    double cm = 1.0e20;

    if(sa>1.0e-20)
    cm = MIN(cm, cmax/(2.0*cfl*sa));

    if(sc>1.0e-20)
    cm = MIN(cm, 1.0/sc);

    return cm;
}
