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

#include"ioflow_amr_block.h"
#include"fdm2D.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"wave_interface.h"
#include<cmath>
#include<iostream>

// Relaxation and beach zones on the faces of a block 1.. of a multi-block grid (B 520 method 1 or 2 with
// B 522 b, NHFLOW), as the zones of iowave on the native grid (iowave_nhflow_relax.cpp, bc_zone_set::relax_weight,
// beach_weight).  The zone of the -i (+i) face reaches its B 521 width w into the block: the distance d of a cell
// is measured from the centre of its face along the grid line through the cell centres (on a straight block the
// distance to the face, as iowave measures it to the edge of the native grid).  The weight of a relaxation zone
// is r = 1 - (exp(x^3.5) - 1)/(e - 1) with x = 1 - d/w, of a beach zone with x = (w - d)/(w fac) (fac 2 with
// B 99 1 on level 0), r = 1 outside, and a field relaxes as f = (1 - r) f_target + r f.
//
// Targets: a relaxation zone the still water depth + its background (B 523) + the ramped waves (B 101, B 102) of
// its B 524 sources (all sources without B 524), the waves only with B 90 1 and B 98 2 on level 0 (the native
// generation zones need them as well); a beach zone its background (still water without B 523), the non-hydrostatic pressure
// to zero.  With B 95 the targets are those at the time of the stage output and the weights r^(dt/dt_ref).
// The waves are evaluated at the cell centres in the frame of B 105 (iowave xgen, ygen), the velocities at the
// sigma levels of the cells.

namespace
{
inline double zone_r(double x)
{
    return 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);
}
}

void ioflow_amr_block::open_zones(lexer *p, const bc_zone *zm, const bc_zone *zp, const background_state &bg, bool waves,
                                  int i0, int ni, int j0, int nj, const bc_zone *zjm, const bc_zone *zjp)
{
    zi0 = i0; zni = ni; zj0 = j0; znj = nj;
    rzone.clear();
    rz[0] = rz[1] = rz[2] = rz[3] = -1;
    const bc_zone *zf[4] = {zm,zp,zjm,zjp};
    for(int f=0; f<4; ++f)
    if(zf[f]!=nullptr)
    {
        rz[f] = int(rzone.size());
        rzone.push_back(*zf[f]);
    }
    gen_waves = waves;
    zbuilt = false;

    if(!rzone.empty())
    {
        bgs = bg;
        bgs.update(p,p->simtime);
    }
}

// the wave members of the lexer (RUN members of lexer_amr.cpp)
#define IOAMR_WAVE_D(X) X(phiin) X(phiout) X(wavetime) X(wts) X(wte) X(wT) X(wV) X(wH) X(wA) X(wL) X(wd) X(ww) X(wk) \
    X(wC) X(wHs) X(wAs) X(wwp) X(ww_s) X(ww_e) X(wTp) X(wLp)

void ioflow_amr_block::wave_get(const lexer *p, wave_run &w)
{
    int q=0;
#define G(m) w.v[q++] = p->m;
    IOAMR_WAVE_D(G)
#undef G
    w.n = p->wN;
}

void ioflow_amr_block::wave_set(lexer *p, const wave_run &w)
{
    int q=0;
#define S(m) p->m = w.v[q++];
    IOAMR_WAVE_D(S)
#undef S
    p->wN = w.n;
}

// the block's wave members in, the ones of level 0 kept in w
void ioflow_amr_block::wave_in(lexer *p, wave_run &w) const
{
    wave_get(p,w);
    wave_set(p,wown);
}

// the wave members of level 0 back (the block keeps its own, the time aside)
void ioflow_amr_block::wave_out(lexer *p, const wave_run &w) const
{
    wave_set(p,w);
}

// the waves (B 92, the sources of B 500-504) on the grid of the block, built at the first use
bool ioflow_amr_block::wave_on(lexer *p, ghostcell *pgc)
{
    if(pwave==nullptr)
    {
        wave_run l0;
        wave_get(p,l0);
        pwave = new wave_interface(p,pgc);
        wave_get(p,wown);
        wave_set(p,l0);
    }
    return true;
}

// iowave::ramp (B 101, B 102)
double ioflow_amr_block::wramp(lexer *p) const
{
    double f=1.0;
    const double tr = (p->B95!=0.0) ? p->wavetime : p->simtime;

    if(p->B101==1 && tr<p->B102*p->wT)
    f = tr/(p->B102*p->wT) - (1.0/PI)*sin(PI*(tr/(p->B102*p->wT)));

    if(p->B101==2 && tr<p->B102)
    f = tr/(p->B102) - (1.0/PI)*sin(PI*(tr/(p->B102)));

    return f;
}

// the frame of the waves (iowave xgen, ygen): along B 105_1 from (B 105_2, B 105_3)
void ioflow_amr_block::wave_frame(lexer *p, double x, double y, double &xg, double &yg) const
{
    const double g = p->B105_1*(PI/180.0);
    xg = (x-p->B105_2)*cos(g) + (y-p->B105_3)*sin(g);
    yg = -(x-p->B105_2)*sin(g) + (y-p->B105_3)*cos(g);
}

// the non-hydrostatic pressure of the incoming waves at the nodes of the -i edge column (i,j) (iowave::nhflow_edge_pressure):
// p(z) = rho int_z^eta dw/dt dz, the Dirichlet value of the Poisson equation at the edge
void ioflow_amr_block::edge_pressure(lexer *p, fdm_nhf *d, ghostcell *pgc, double xs, double ys)
{
    if(d->Pbc==nullptr)
    p->Darray(d->Pbc,p->imax*p->jmax*(p->kmax+2));

    const double rw = wramp(p);
    const double dt = p->dt>0.0 ? p->dt : 1.0e-3;
    const double wt = p->wavetime;
    std::vector<double> a(p->knoz+1), zz(p->knoz+1);

    for(k=0; k<=p->knoz; ++k)
    {
        zz[k] = p->ZN[KP]*d->WL(i,j) + d->bed(i,j);
        a[k] = pwave->wave_w(p,pgc,xs,ys,zz[k]-p->phimean);
    }

    p->wavetime = wt - dt;
    for(k=0; k<=p->knoz; ++k)
    a[k] = rw*(a[k] - pwave->wave_w(p,pgc,xs,ys,zz[k]-p->phimean))/dt;
    p->wavetime = wt;

    double pn = 0.0;
    k = p->knoz;
    d->Pbc[FIm1JK] = 0.0;
    for(k=p->knoz-1; k>=0; --k)
    {
        pn += p->W1*0.5*(a[k]+a[k+1])*(zz[k+1]-zz[k]);
        d->Pbc[FIm1JK] = pn;
    }
}

// the weights of the columns: the distance from the face along the grid line through the cell centres (the cells of
// the block: i zi0 .. zi0+zni-1, j zj0 .. zj0+znj-1 of the grid, its faces at the nodes zi0 and zi0+zni)
void ioflow_amr_block::zones_build(lexer *p)
{
    zbuilt = true;
    const size_t n2 = size_t(p->imax)*size_t(p->jmax);
    zcol.clear();
    vector<double> xg(n2,-1.0), rb(n2,0.0);
    vector<int> cb(n2,0), zg(n2,-1), zb(n2,-1);

    for(int f=0; f<4; ++f)
    {
        if(rz[f]<0)
        continue;

        const bc_zone &z = rzone[rz[f]];
        const double w = z.d;

        // the lines across the face: along i from an i face, along j from a j face
        const bool jf = (f>=2);
        const int nline = jf ? zni : znj, nlen = jf ? znj : zni;

        for(int l=0; l<nline; ++l)
        {
            int a, c;
            if(!jf)
            {
                j = zj0+l;
                const int in = (f==0) ? zi0 : zi0+zni;
                a = p->nij(in,j);
                c = p->nij(in,j+1);
            }
            else
            {
                i = zi0+l;
                const int jn = (f==2) ? zj0 : zj0+znj;
                a = p->nij(i,jn);
                c = p->nij(i+1,jn);
            }
            double xo = 0.5*(p->XN2D[a] + p->XN2D[c]);
            double yo = 0.5*(p->YN2D[a] + p->YN2D[c]);
            double s = 0.0;

            for(int m=0; m<nlen; ++m)
            {
                if(!jf)
                i = (f==0) ? zi0+m : zi0+zni-1-m;
                else
                j = (f==2) ? zj0+m : zj0+znj-1-m;
                const double xc = p->XC2D[IJ], yc = p->YC2D[IJ];
                const double ds = sqrt((xc-xo)*(xc-xo) + (yc-yo)*(yc-yo));
                s += ds;
                xo = xc;
                yo = yc;

                if(s>=w)
                break;

                // the zone reaches beyond the other face of the block (a G 55 piece): it ends there
                if(m==nlen-1 && s+ds<w && !ztold)
                {
                    ztold = true;
                    std::cout<<"!!! zone "<<z.id<<" (B 522) is "<<w<<" m wide, its block (G 55: its piece) ends "<<s
                             <<" m from its face: the zone ends there !!!"<<std::endl;
                }

                if(z.method==bc_method::relax)
                {
                    // the closest zone sets the weight and the target (bc_zone_set::relax_weight)
                    const double x = MAX(1.0-s/w,0.0);
                    if(x>xg[IJ])
                    {
                        xg[IJ] = x;
                        zg[IJ] = rz[f];
                    }
                }
                else
                {
                    // overlapping beach zones: the mean of their weights (bc_zone_set::beach_weight)
                    const double x = MAX((w-s)/(w*z.fac),0.0);
                    rb[IJ] += zone_r(x);
                    ++cb[IJ];
                    zb[IJ] = rz[f];
                }
            }
        }
    }

    zgen0.clear(); zbch0.clear(); zgz.clear(); zbz.clear();
    for(i=zi0; i<zi0+zni; ++i)
    for(j=zj0; j<zj0+znj; ++j)
    if(zg[IJ]>=0 || zb[IJ]>=0)
    {
        zcol.push_back(IJ);
        zgen0.push_back(zg[IJ]>=0 ? zone_r(xg[IJ]) : 1.0);
        zbch0.push_back(zb[IJ]>=0 ? rb[IJ]/double(cb[IJ]) : 1.0);
        zgz.push_back(zg[IJ]);
        zbz.push_back(zb[IJ]);
    }
    zgen = zgen0;
    zbch = zbch0;
    zfac = 1.0;
}

// B 95 (with B 98 2 on level 0, iowave::relax_stage_factor): the weights to the power dt / dt_ref
void ioflow_amr_block::zone_stage(lexer *p)
{
    if(!zbuilt)
    zones_build(p);

    double dtref = p->B95;
    if(p->B95<0.0)
    dtref = p->dt;

    const double f = p->dt/(dtref>1.0e-20 ? dtref : 1.0e-20);
    if(f==zfac)
    return;

    zfac = f;
    for(size_t q=0; q<zcol.size(); ++q)
    {
        const double rw = zgen0[q], rb = zbch0[q];
        zgen[q] = (rw>0.0 && rw<1.0) ? pow(rw,f) : rw;
        zbch[q] = (rb>0.0 && rb<1.0) ? pow(rb,f) : rb;
    }
}

void ioflow_amr_block::wavegen_stage_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double t)
{
    if(rzone.empty() || p->B95==0.0 || !gen_waves)
    return;

    zone_stage(p);
    tstage = t;
}

// the time of the targets: B 95 the stage, else the start of the step
double ioflow_amr_block::zone_time(lexer *p) const
{
    return (p->B95!=0.0 && gen_waves && tstage>-1.0e19) ? tstage : p->simtime;
}

// every column (i, j set) of the zones: f(q, r_gen, r_beach, zone_gen, zone_beach)
template<class F>
void ioflow_amr_block::zone_columns(lexer *p, F &&fn)
{
    if(!zbuilt)
    zones_build(p);

    for(size_t q=0; q<zcol.size(); ++q)
    {
        i = zcol[q]/p->jmax + p->imin;
        j = zcol[q]%p->jmax + p->jmin;
        fn(q,zgen[q],zbch[q],zgz[q],zbz[q]);
    }
}

void ioflow_amr_block::zone_cells(lexer *p, const std::function<void(int,int)> &f)
{
    if(rzone.empty())
    return;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        if(zg>=0 || zb>=0)
        f(i,j);
    });
}

// the turbulence in the zones (k, epsilon/omega, the eddy viscosities): towards 0 with the weight of the zone, as
// iowave::turb_relax_nhflow
void ioflow_amr_block::turb_relax_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *F)
{
    if(rzone.empty())
    return;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        for(k=0; k<p->knoz; ++k)
        if(p->flag4[IJK]>0)
        {
            if(zg>=0)
            F[IJK] = rg*F[IJK];
            if(zb>=0)
            F[IJK] = rb*F[IJK];
        }
    });
}

void ioflow_amr_block::WL_relax(lexer *p, ghostcell *pgc, slice &WL, slice &depth)
{
    if(rzone.empty())
    return;

    const bool waves = gen_waves && wave_on(p,pgc);
    wave_run l0;
    if(waves)
    {
    wave_in(p,l0);
    p->wavetime = zone_time(p);
    }
    const double rmp = waves ? wramp(p) : 0.0;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        const double x = p->XC2D[IJ], y = p->YC2D[IJ];

        // relaxation zone: still water depth + background + the ramped waves
        if(zg>=0 && rg<1.0 && p->wet[IJ]==1 && p->deep[IJ]==1)
        {
            const bc_zone &z = rzone[zg];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            double ew = 0.0;
            if(waves)
            {
                double xg,yg;
                wave_frame(p,x,y,xg,yg);
                pwave->select_sources(z.sources.empty() ? nullptr : &z.sources);
                ew = pwave->wave_eta(p,pgc,xg,yg);
                pwave->select_sources(nullptr);
            }
            const double eb = (b>=0) ? bgs.eta(b,x,y) : 0.0;
            WL(i,j) = (1.0-rg)*(depth(i,j) + eb + rmp*ew) + rg*WL(i,j);
        }

        // beach zone: still water depth + background
        if(zb>=0 && rb<1.0 && p->wet[IJ]==1 && !zgen_only)
        {
            const bc_zone &z = rzone[zb];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            const double eb = (b>=0) ? bgs.eta(b,x,y) : 0.0;
            WL(i,j) = (1.0-rb)*(depth(i,j) + eb) + rb*WL(i,j);
        }
    });

    if(waves)
    wave_out(p,l0);
}

// the targets of the relaxation zones as iowave::U_relax, V_relax, W_relax (the inflow current B 60 aside)
void ioflow_amr_block::U_relax(lexer *p, ghostcell *pgc, double *U, double *UH)
{
    if(rzone.empty())
    return;

    const bool waves = gen_waves && wave_on(p,pgc);
    wave_run l0;
    if(waves)
    {
    wave_in(p,l0);
    p->wavetime = zone_time(p);
    }
    const double rmp = waves ? wramp(p) : 0.0;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        const double x = p->XC2D[IJ], y = p->YC2D[IJ];
        const double h0 = p->wd - p->bed[IJ];

        if(zg>=0 && rg<1.0 && p->wet[IJ]==1 && p->deep[IJ]==1)
        {
            const bc_zone &z = rzone[zg];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            double xg=0.0, yg=0.0, ew=0.0;
            if(waves)
            {
                wave_frame(p,x,y,xg,yg);
                pwave->select_sources(z.sources.empty() ? nullptr : &z.sources);
                ew = pwave->wave_eta(p,pgc,xg,yg);
            }
            double ub=0.0, vb=0.0, eb=0.0;
            if(b>=0)
            {
                bgs.vel(b,h0,x,y,ub,vb);
                eb = bgs.eta(b,x,y);
            }
            for(k=0; k<p->knoz; ++k)
            if(p->flag4[IJK]>0)
            {
                const double uw = waves ? pwave->wave_u(p,pgc,xg,yg,p->ZSP[IJK]-p->phimean) : 0.0;
                if(b<0)
                {
                    U[IJK]  = (1.0-rg)*rmp*uw + rg*U[IJK];
                    UH[IJK] = (1.0-rg)*rmp*(ew + h0)*uw + rg*UH[IJK];
                }
                else
                {
                    const double fp = bgs.profiled(b) ? bg_prof(p,b,k,h0+eb) : 1.0;
                    const double ut = ub*fp + rmp*uw;
                    const double ht = h0 + eb + rmp*ew;
                    U[IJK]  = (1.0-rg)*ut + rg*U[IJK];
                    UH[IJK] = (1.0-rg)*ht*ut + rg*UH[IJK];
                }
            }
            if(waves)
            pwave->select_sources(nullptr);
        }

        if(zb>=0 && rb<1.0)
        {
            const bc_zone &z = rzone[zb];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            double ub=0.0, vb=0.0, eb=0.0;
            if(b>=0)
            {
                bgs.vel(b,h0,x,y,ub,vb);
                eb = bgs.eta(b,x,y);
            }
            for(k=0; k<p->knoz; ++k)
            if(p->flag4[IJK]>0)
            {
                if(b<0)
                {
                    U[IJK]  = rb*U[IJK];
                    UH[IJK] = rb*UH[IJK];
                }
                else
                {
                    const double fp = bgs.profiled(b) ? bg_prof(p,b,k,h0+eb) : 1.0;
                    U[IJK]  = rb*U[IJK]  + (1.0-rb)*ub*fp;
                    UH[IJK] = rb*UH[IJK] + (1.0-rb)*ub*fp*(h0+eb);
                }
            }
        }
    });

    if(waves)
    wave_out(p,l0);
}

void ioflow_amr_block::V_relax(lexer *p, ghostcell *pgc, double *V, double *VH)
{
    if(rzone.empty() || p->j_dir==0)
    return;

    const bool waves = gen_waves && wave_on(p,pgc);
    wave_run l0;
    if(waves)
    {
    wave_in(p,l0);
    p->wavetime = zone_time(p);
    }
    const double rmp = waves ? wramp(p) : 0.0;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        const double x = p->XC2D[IJ], y = p->YC2D[IJ];
        const double h0 = p->wd - p->bed[IJ];

        if(zg>=0 && rg<1.0 && p->wet[IJ]==1 && p->deep[IJ]==1)
        {
            const bc_zone &z = rzone[zg];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            double xg=0.0, yg=0.0, ew=0.0;
            if(waves)
            {
                wave_frame(p,x,y,xg,yg);
                pwave->select_sources(z.sources.empty() ? nullptr : &z.sources);
                ew = pwave->wave_eta(p,pgc,xg,yg);
            }
            double ub=0.0, vb=0.0, eb=0.0;
            if(b>=0)
            {
                bgs.vel(b,h0,x,y,ub,vb);
                eb = bgs.eta(b,x,y);
            }
            for(k=0; k<p->knoz; ++k)
            if(p->flag4[IJK]>0)
            {
                const double vw = waves ? pwave->wave_v(p,pgc,xg,yg,p->ZSP[IJK]-p->phimean) : 0.0;
                if(b<0)
                {
                    V[IJK]  = (1.0-rg)*rmp*vw + rg*V[IJK];
                    VH[IJK] = (1.0-rg)*rmp*(ew + h0)*vw + rg*VH[IJK];
                }
                else
                {
                    const double fp = bgs.profiled(b) ? bg_prof(p,b,k,h0+eb) : 1.0;
                    const double vt = vb*fp + rmp*vw;
                    const double ht = h0 + eb + rmp*ew;
                    V[IJK]  = (1.0-rg)*vt + rg*V[IJK];
                    VH[IJK] = (1.0-rg)*ht*vt + rg*VH[IJK];
                }
            }
            if(waves)
            pwave->select_sources(nullptr);
        }

        if(zb>=0 && rb<1.0)
        {
            const bc_zone &z = rzone[zb];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            double ub=0.0, vb=0.0, eb=0.0;
            if(b>=0)
            {
                bgs.vel(b,h0,x,y,ub,vb);
                eb = bgs.eta(b,x,y);
            }
            for(k=0; k<p->knoz; ++k)
            if(p->flag4[IJK]>0)
            {
                if(b<0)
                {
                    V[IJK]  = rb*V[IJK];
                    VH[IJK] = rb*VH[IJK];
                }
                else
                {
                    const double fp = bgs.profiled(b) ? bg_prof(p,b,k,h0+eb) : 1.0;
                    V[IJK]  = rb*V[IJK]  + (1.0-rb)*vb*fp;
                    VH[IJK] = rb*VH[IJK] + (1.0-rb)*vb*fp*(h0+eb);
                }
            }
        }
    });

    if(waves)
    wave_out(p,l0);
}

void ioflow_amr_block::W_relax(lexer *p, ghostcell *pgc, double *W, double *WH)
{
    if(rzone.empty())
    return;

    const bool waves = gen_waves && wave_on(p,pgc);
    wave_run l0;
    if(waves)
    {
    wave_in(p,l0);
    p->wavetime = zone_time(p);
    }
    const double rmp = waves ? wramp(p) : 0.0;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        const double x = p->XC2D[IJ], y = p->YC2D[IJ];
        const double h0 = p->wd - p->bed[IJ];

        if(zg>=0 && rg<1.0 && p->wet[IJ]==1 && p->deep[IJ]==1)
        {
            const bc_zone &z = rzone[zg];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            double xg=0.0, yg=0.0, ew=0.0;
            if(waves)
            {
                wave_frame(p,x,y,xg,yg);
                pwave->select_sources(z.sources.empty() ? nullptr : &z.sources);
                ew = pwave->wave_eta(p,pgc,xg,yg);
            }
            const double eb = (b>=0) ? bgs.eta(b,x,y) : 0.0;
            for(k=0; k<p->knoz; ++k)
            if(p->flag4[IJK]>0)
            {
                const double ww = waves ? pwave->wave_w(p,pgc,xg,yg,p->ZSP[IJK]-p->phimean) : 0.0;
                const double ht = (b<0) ? ew + h0 : h0 + eb + rmp*ew;
                W[IJK]  = (1.0-rg)*rmp*ww + rg*W[IJK];
                WH[IJK] = (b<0) ? (1.0-rg)*rmp*ht*ww + rg*WH[IJK] : (1.0-rg)*ht*rmp*ww + rg*WH[IJK];
            }
            if(waves)
            pwave->select_sources(nullptr);
        }

        if(zb>=0 && rb<1.0)
        for(k=0; k<p->knoz; ++k)
        if(p->flag4[IJK]>0)
        {
            W[IJK]  = rb*W[IJK];
            WH[IJK] = rb*WH[IJK];
        }
    });

    if(waves)
    wave_out(p,l0);
}

void ioflow_amr_block::P_relax(lexer *p, ghostcell *pgc, double *P)
{
    if(rzone.empty())
    return;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        if(zb>=0 && rb<1.0)
        for(k=0; k<=p->knoz; ++k)
        P[FIJK] = rb*P[FIJK];
    });
}

// the zones at the start, as iowave::ini_nhflow and driver_ini_nhflow: the velocities relaxed once, the water level of
// the relaxation zones 20 times (their surface is the target then; the beach zones start from the initial state), eta
// from the new water level
void ioflow_amr_block::zones_ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(rzone.empty())
    return;

    U_relax(p,pgc,d->U,d->UH);
    V_relax(p,pgc,d->V,d->VH);
    W_relax(p,pgc,d->W,d->WH);

    zgen_only = true;
    for(int qn=0; qn<20; ++qn)
    WL_relax(p,pgc,d->WL,d->depth);
    zgen_only = false;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        if(zg>=0 && rg<1.0 && p->wet[IJ]==1 && p->deep[IJ]==1)
        {
            d->eta(i,j) = d->WL(i,j) - d->depth(i,j);
            d->eta_n(i,j) = d->eta(i,j);
            p->WL[IJ] = d->WL(i,j);
        }
    });
}

// ------------------------------------------------------------------------------------------------- SFLOW
// The zones of an SFLOW block as iowave_sflow_relax.cpp on the native grid: a relaxation zone relaxes eta towards
// its background + the ramped waves (also dry columns), the depth-averaged velocities of the wet columns towards the
// background current + the ramped waves averaged over B 160 + 1 levels from the bed to the wave surface
// (iowave::wavegen_2D_precalc), W towards the waves; a beach zone eta (wet columns) towards its background, U, V
// towards the background current, W and the non-hydrostatic pressure towards zero.  The targets at the start of the
// step (p->simtime).

// the waves of a column at (x, y): eta and the depth-averaged u, v, w (not ramped)
void ioflow_amr_block::wave_col2D(lexer *p, ghostcell *pgc, const bc_zone &z, double eb, double &ew, double &uw, double &vw, double &ww)
{
    double xg,yg;
    wave_frame(p,p->XC2D[IJ],p->YC2D[IJ],xg,yg);
    pwave->select_sources(z.sources.empty() ? nullptr : &z.sources);
    ew = pwave->wave_eta(p,pgc,xg,yg);
    const double dz = (eb + ew + p->wd - zb2->bed(i,j))/double(p->B160);
    double su=0.0, sv=0.0, sw=0.0;
    double zl = -p->wd;
    for(int qn=0; qn<=p->B160; ++qn)
    {
        su += pwave->wave_u(p,pgc,xg,yg,zl);
        sv += pwave->wave_v(p,pgc,xg,yg,zl);
        sw += pwave->wave_w(p,pgc,xg,yg,zl);
        zl += dz;
    }
    uw = su/double(p->B160+1);
    vw = (p->j_dir==1) ? sv/double(p->B160+1) : 0.0;
    ww = sw/double(p->B160+1);
    pwave->select_sources(nullptr);
}

void ioflow_amr_block::eta_relax(lexer *p, ghostcell *pgc, slice &f)
{
    if(rzone.empty() || zb2==nullptr)
    return;

    const bool waves = gen_waves && wave_on(p,pgc);
    wave_run l0;
    if(waves)
    {
    wave_in(p,l0);
    p->wavetime = p->simtime;
    }
    const double rmp = waves ? wramp(p) : 0.0;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        const double x = p->XC2D[IJ], y = p->YC2D[IJ];

        if(zg>=0 && rg<1.0)
        {
            const bc_zone &z = rzone[zg];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            const double eb = (b>=0) ? bgs.eta(b,x,y) : 0.0;
            double ew = 0.0;
            if(waves)
            {
                double xg,yg;
                wave_frame(p,x,y,xg,yg);
                pwave->select_sources(z.sources.empty() ? nullptr : &z.sources);
                ew = pwave->wave_eta(p,pgc,xg,yg);
                pwave->select_sources(nullptr);
            }
            f(i,j) = (1.0-rg)*(eb + rmp*ew) + rg*f(i,j);
        }

        if(zb>=0 && rb<1.0 && p->wet[IJ]==1)
        {
            const bc_zone &z = rzone[zb];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            const double eb = (b>=0) ? bgs.eta(b,x,y) : 0.0;
            f(i,j) = rb*f(i,j) + (1.0-rb)*eb;
        }
    });

    if(waves)
    wave_out(p,l0);
}

// comp 0 u, 1 v, 2 w
void ioflow_amr_block::uvw_relax2D(lexer *p, ghostcell *pgc, slice &U, slice &UH, slice &WL, int comp)
{
    if(rzone.empty() || zb2==nullptr || (comp==1 && p->j_dir==0))
    return;

    const bool waves = gen_waves && wave_on(p,pgc);
    wave_run l0;
    if(waves)
    {
    wave_in(p,l0);
    p->wavetime = p->simtime;
    }
    const double rmp = waves ? wramp(p) : 0.0;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        const double x = p->XC2D[IJ], y = p->YC2D[IJ];
        const double h0 = zb2->depth(i,j);

        if(zg>=0 && rg<1.0 && p->wet[IJ]==1)
        {
            const bc_zone &z = rzone[zg];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            double ub=0.0, vb=0.0;
            if(b>=0 && comp<2)
            bgs.vel(b,h0,x,y,ub,vb);
            double ew=0.0, uw=0.0, vw=0.0, ww=0.0;
            if(waves)
            wave_col2D(p,pgc,z,(b>=0) ? bgs.eta(b,x,y) : 0.0,ew,uw,vw,ww);
            const double tg = (comp==0) ? ub + rmp*uw : (comp==1) ? vb + rmp*vw : rmp*ww;
            U(i,j)  = (1.0-rg)*tg + rg*U(i,j);
            UH(i,j) = (1.0-rg)*tg*WL(i,j) + rg*UH(i,j);
        }

        if(zb>=0 && rb<1.0)
        {
            const bc_zone &z = rzone[zb];
            const int b = z.bg>0 ? bgs.index(z.bg) : -1;
            double ub=0.0, vb=0.0;
            if(b>=0 && comp<2)
            bgs.vel(b,h0,x,y,ub,vb);
            const double tg = (comp==0) ? ub : (comp==1) ? vb : 0.0;
            U(i,j)  = rb*U(i,j) + (1.0-rb)*tg;
            UH(i,j) = rb*UH(i,j) + (1.0-rb)*tg*WL(i,j);
        }
    });

    if(waves)
    wave_out(p,l0);
}

void ioflow_amr_block::um_relax(lexer *p, ghostcell *pgc, slice &U, slice &UH, slice &WL)
{
    uvw_relax2D(p,pgc,U,UH,WL,0);
}

void ioflow_amr_block::vm_relax(lexer *p, ghostcell *pgc, slice &V, slice &VH, slice &WL)
{
    uvw_relax2D(p,pgc,V,VH,WL,1);
}

void ioflow_amr_block::wm_relax(lexer *p, ghostcell *pgc, slice &W, slice &WH, slice &WL)
{
    uvw_relax2D(p,pgc,W,WH,WL,2);
}

void ioflow_amr_block::pm_relax(lexer *p, ghostcell *pgc, slice &f)
{
    if(rzone.empty())
    return;

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        if(zb>=0 && rb<1.0)
        f(i,j) = rb*f(i,j);
    });
}

// the zones at the start, as sflow_f::ini: U, V relaxed once with the water depth of the initial state, then eta
// (sflow_f::ini_fsf_2) and the water depth from it
void ioflow_amr_block::zones_ini2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
    zb2 = b;
    if(rzone.empty())
    return;

    um_relax(p,pgc,b->U,b->UH,b->WL);
    vm_relax(p,pgc,b->V,b->VH,b->WL);
    eta_relax(p,pgc,b->eta);

    zone_columns(p,[&](size_t q, double rg, double rb, int zg, int zb)
    {
        b->eta_n(i,j) = b->eta(i,j);
        b->WL(i,j) = MAX(b->eta(i,j) + b->depth(i,j),0.0);
    });
}
