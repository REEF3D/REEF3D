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
#include"nhflow_tracer_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<cmath>
#include<algorithm>

// ---------------------------------------------------------------------------
// Density from the tracers and its pressure gradient (A 620; P12b), hydrostatic and Boussinesq:
//
//     p = g rho0 (eta - z) + p',   p'(z) = g int_z^eta (rho - rho0) dz'
//     d(UH)/dt = ... - g WL grad(eta) - WL (1/rho0) grad_z p'      (the first term in the eta system)
//
// The baroclinic term is explicit and acts on the layer face discharges (face predictor, si_assemble) and on
// the cell momentum (cell_update), not in the RK stages. Its density is that of a predictor of the tracers
// with the old part (1 - theta) Q^n of the transport (nhflow_momentum_SI::start, nhflow_amr::si_step: the
// tracers saved, moved, the density computed, the tracers restored); the step after the momentum moves them
// with the whole transport T = theta Q^{n+1} + (1 - theta) Q^n as without density. Forward-backward for the
// internal waves, neutral up to omega dt = 2; with the density of t^n instead the transport T amplifies them
// by 1 + (1 - theta) (omega dt)^2 per step.
//
// The gradient along the sigma layers (per face between the cells L and R, normal n from L to R):
//
//     (1/rho0) [ (P_R - P_L) + g/2 ((rho'_R + rho'_L)(z_R - z_L) - 1/5 (...)) ] / dn
//
// A 622 1 (default): Shchepetkin and McWilliams (2003), the density Jacobian with the cubic corrections of
// the harmonic means of the differences (ROMS prsgrd32): P from the surface with the same corrections in the
// vertical, the horizontal ones from the columns LL, L, R, RR (a missing outer column: its difference
// continued). A 622 2: the same without the corrections (the standard two-term form).
// Both are the transpose of the tracer advection with centred face values (the energy of the internal
// waves conserved by the space discretisation). An integral at constant z (Stelling and van Kester 1994) was
// tried: no force from a level stratification, but not of that form, and the stratified fluid at rest over a
// steep ridge was unstable with it.
//
// Equation of state (A 620): 1 linear rho = rho0 (1 - alpha (T - T0) + beta (S - S0)) (A 621);
// 2 UNESCO 1981 at one atmosphere (Millero and Poisson 1981). Without temperature (A 600 1) T = T0.
// ---------------------------------------------------------------------------

double nhflow_momentum_SI::eos(lexer *p, double S, double T)
{
    if(p->A620==2)
    {
        const double T2=T*T, T3=T2*T, T4=T3*T, T5=T4*T;
        const double rw = 999.842594 + 6.793952e-2*T - 9.095290e-3*T2 + 1.001685e-4*T3 - 1.120083e-6*T4 + 6.536332e-9*T5;
        const double A = 8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T2 - 8.2467e-7*T3 + 5.3875e-9*T4;
        const double B = -5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T2;
        const double C = 4.8314e-4;
        const double Sp = S>0.0 ? S : 0.0;
        return rw + A*Sp + B*Sp*sqrt(Sp) + C*Sp*Sp;
    }

    return p->A621_rho0*(1.0 - p->A621_alpha*(T - p->A621_t0) + p->A621_beta*(S - p->A621_s0));
}

// the density of all cells (ghost cells included) from the tracers
void nhflow_momentum_SI::density(lexer *p, fdm_nhf *d)
{
    const int K = p->knoz;

    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
        for(int kk=0; kk<K; ++kk)
        {
            const int c = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
            const double S = d->SALT!=nullptr ? d->SALT[c] : p->A621_s0;
            const double T = d->TEMP!=nullptr ? d->TEMP[c] : p->A621_t0;
            d->RO[c] = eos(p,S,T);
        }

        const int c0 = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax - p->kmin;
        d->RO[c0-1] = d->RO[c0];
        d->RO[c0+K] = d->RO[c0+K-1];
    }
}

// the speed of the first internal mode of column (i,j), 1/2 sqrt(g' h) with g' from the density range of the
// column (two-layer bound; the time step, nhflow_timestep_SI)
double nhflow_momentum_SI::internal_speed(lexer *p, fdm_nhf *d, int ii, int jj)
{
    if(p->A620==0)
    return 0.0;

    double lo=1.0e30, hi=-1.0e30;

    for(int kk=0; kk<p->knoz; ++kk)
    {
    const double r = d->RO[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin];
    lo = std::min(lo,r);
    hi = std::max(hi,r);
    }

    const double gp = fabs(p->W22)*std::max(hi-lo,0.0)/p->A621_rho0;
    return 0.5*sqrt(gp*std::max(d->WL(ii,jj),0.0));
}

// the pressure p' = g int_z^eta (rho - rho0) dz' of column (a,b) at its layer centres: the top half layer with
// the density continued linearly, below it the trapezoid between the centres with (A 622 1) the cubic
// corrections of the harmonic means of the vertical differences (Shchepetkin and McWilliams 2003)
void nhflow_momentum_SI::baroclinic_column(lexer *p, fdm_nhf *d, int a, int b, double *P)
{
    const int K = p->knoz;
    const int mg = increment::marge;
    const double g = fabs(p->W22), rho0 = p->A621_rho0;
    const double h = d->WL(a,b), zb = d->bed(a,b);
    const bool sm = p->A622!=2;

    bcr.resize(K);
    bcz.resize(K);
    bcd.resize(K);
    bce.resize(K);
    vector<double> &r = bcr, &z = bcz, &dr = bcd, &dz = bce;

    for(int k=0; k<K; ++k)
    {
    r[k] = d->RO[(a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + k-p->kmin] - rho0;
    z[k] = zb + p->ZP[k+mg]*h;
    }

    const double eta = zb + h;

    if(K==1)
    {
    P[0] = g*r[0]*(eta - z[0]);
    return;
    }

    // the derivatives at the centres: harmonic means of the differences above and below (zero at an extremum)
    for(int k=0; k<K; ++k)
    {
        const int ka = k<K-1 ? k : K-2, kb = k>0 ? k-1 : 0;
        const double ru = r[ka+1] - r[ka], rd = r[kb+1] - r[kb];
        const double zu = z[ka+1] - z[ka], zd = z[kb+1] - z[kb];
        const double c = 2.0*ru*rd;

        dr[k] = c>1.0e-14 ? c/(ru + rd) : 0.0;
        dz[k] = (zu + zd)>0.0 ? 2.0*zu*zd/(zu + zd) : 0.0;
    }

    const double c2 = 0.5*(r[K-1] - r[K-2])*(eta - z[K-1])/(z[K-1] - z[K-2]);
    P[K-1] = g*(r[K-1] + c2)*(eta - z[K-1]);

    const double f12 = 1.0/12.0;

    for(int k=K-2; k>=0; --k)
    {
        double v = (r[k+1] + r[k])*(z[k+1] - z[k]);

        if(sm)
        v -= 0.2*((dr[k+1] - dr[k])*(z[k+1] - z[k] - f12*(dz[k+1] + dz[k])) - (dz[k+1] - dz[k])*(r[k+1] - r[k] - f12*(dr[k+1] + dr[k])));

        P[k] = P[k+1] + 0.5*g*v;
    }
}

// the normal component (L -> R) of (1/rho0) grad_z p' on the face (dir, fi, fj), per layer
void nhflow_momentum_SI::baroclinic_face(lexer *p, fdm_nhf *d, int dir, int fi, int fj, double *bn)
{
    const int K = p->knoz;
    const int mg = increment::marge;
    const double g = fabs(p->W22), rho0 = p->A621_rho0;
    const bool sm = p->A622!=2;
    const int li = fi, lj = fj, ri = dir==1 ? fi+1 : fi, rj = dir==2 ? fj+1 : fj;
    const int lli = dir==1 ? li-1 : li, llj = dir==2 ? lj-1 : lj;
    const int rri = dir==1 ? ri+1 : ri, rrj = dir==2 ? rj+1 : rj;
    const double dn = fdn(p,dir,fi,fj);

    auto c3 = [&](int a, int b, int k) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + k-p->kmin; };
    auto zc = [&](int a, int b, int k) { return d->bed(a,b) + p->ZP[k+mg]*d->WL(a,b); };

    // the outer columns: fluid (the tracer mask, also deep in the halo and around a block) and wet
    auto colwet = [&](int a, int b)
    {
        return a>=p->imin && a<p->imin+p->imax && b>=p->jmin && b<p->jmin+p->jmax
               && ptrc->fluid(p,a,b) && d->WL(a,b)>p->A544;
    };

    const bool hl = sm && colwet(lli,llj), hr = sm && colwet(rri,rrj);

    bpl.resize(K);
    bpr.resize(K);
    baroclinic_column(p,d,li,lj,bpl.data());
    baroclinic_column(p,d,ri,rj,bpr.data());

    auto harm = [](double a, double b) { const double c = 2.0*a*b; return c>1.0e-14 ? c/(a+b) : 0.0; };
    const double f12 = 1.0/12.0;

    for(int k=0; k<K; ++k)
    {
        const double rL = d->RO[c3(li,lj,k)] - rho0, rR = d->RO[c3(ri,rj,k)] - rho0;
        const double fc = rR - rL, az = zc(ri,rj,k) - zc(li,lj,k);

        double v = (rR + rL)*az;

        if(sm)
        {
            const double fl = hl ? rL - (d->RO[c3(lli,llj,k)] - rho0) : fc;
            const double zl = hl ? zc(li,lj,k) - zc(lli,llj,k) : az;
            const double fr = hr ? (d->RO[c3(rri,rrj,k)] - rho0) - rR : fc;
            const double zr = hr ? zc(rri,rrj,k) - zc(ri,rj,k) : az;

            const double dRl = harm(fl,fc), dRr = harm(fc,fr);
            const double dZl = harm(zl,az), dZr = harm(az,zr);

            v -= 0.2*((dRr - dRl)*(az - f12*(dZr + dZl)) - (dZr - dZl)*(fc - f12*(dRr + dRl)));
        }

        bn[k] = ((bpr[k] - bpl[k]) + 0.5*g*v)/(dn*rho0);
    }
}

// the baroclinic acceleration (1/rho0) grad_z p' per layer: its normal component on the open faces between
// two wet cells (BFx, BFy: the layer face discharges) and the cell vector from the faces of the cell (Bx,
// By: the cell momentum; the mean of the two faces per index direction, the transpose of the face discharges
// as the mean of the cells)
void nhflow_momentum_SI::baroclinic(lexer *p, fdm_nhf *d)
{
    const int K = p->knoz;
    const int N = p->imax*p->jmax*(p->kmax+2);
    vector<double> bf(K);

    auto c3 = [&](int a, int b, int k) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + k-p->kmin; };
    auto wetc = [&](int a, int b) { return p->flagslice4[(a-p->imin)*p->jmax + (b-p->jmin)]>0 && d->WL(a,b)>p->A544; };

    for(int q=0; q<N; ++q)
    BFx[q] = BFy[q] = Bx[q] = By[q] = 0.0;

    // faces (dir 1: ii from -1 so that every own cell has both x faces)
    for(int ii=-1; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    if(Hx(ii,jj)>0.0 && wetc(ii,jj) && wetc(ii+1,jj))
    {
        baroclinic_face(p,d,1,ii,jj,bf.data());
        for(int k=0; k<K; ++k)
        BFx[c3(ii,jj,k)] = bf[k];
    }

    if(p->j_dir==1)
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=-1; jj<p->knoy; ++jj)
    if(Hy(ii,jj)>0.0 && wetc(ii,jj) && wetc(ii,jj+1))
    {
        baroclinic_face(p,d,2,ii,jj,bf.data());
        for(int k=0; k<K; ++k)
        BFy[c3(ii,jj,k)] = bf[k];
    }

    // cells
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        if(!wetc(ii,jj))
        continue;

        const bool oe = Hx(ii,jj)>0.0 && wetc(ii+1,jj);
        const bool ow = Hx(ii-1,jj)>0.0 && wetc(ii-1,jj);
        const bool on = p->j_dir==1 && Hy(ii,jj)>0.0 && wetc(ii,jj+1);
        const bool os = p->j_dir==1 && Hy(ii,jj-1)>0.0 && wetc(ii,jj-1);

        const int na = (oe?1:0) + (ow?1:0), nb = (on?1:0) + (os?1:0);

        if(na==0 && nb==0)
        continue;

        // the mean normals per index direction (curvilinear: unit normals of the faces)
        double ax=1.0, ay=0.0, bx=0.0, by=1.0;

        if(p->geo_curv==1)
        {
            double nx,ny;
            ax = ay = bx = by = 0.0;
            if(oe) {fnorm(p,1,ii,jj,nx,ny); ax+=nx; ay+=ny;}
            if(ow) {fnorm(p,1,ii-1,jj,nx,ny); ax+=nx; ay+=ny;}
            if(on) {fnorm(p,2,ii,jj,nx,ny); bx+=nx; by+=ny;}
            if(os) {fnorm(p,2,ii,jj-1,nx,ny); bx+=nx; by+=ny;}

            const double la = sqrt(ax*ax+ay*ay), lb = sqrt(bx*bx+by*by);
            if(la>0.0) {ax/=la; ay/=la;}
            if(lb>0.0) {bx/=lb; by/=lb;}
        }

        const double det = ax*by - ay*bx;

        for(int k=0; k<K; ++k)
        {
            const double ga = na>0 ? ((oe?BFx[c3(ii,jj,k)]:0.0) + (ow?BFx[c3(ii-1,jj,k)]:0.0))/double(na) : 0.0;
            const double gb = nb>0 ? ((on?BFy[c3(ii,jj,k)]:0.0) + (os?BFy[c3(ii,jj-1,k)]:0.0))/double(nb) : 0.0;
            double gx=0.0, gy=0.0;

            if(na>0 && nb>0 && fabs(det)>1.0e-3)
            {
            gx = (ga*by - ay*gb)/det;
            gy = (ax*gb - ga*bx)/det;
            }

            else
            if(na>0)
            {
            gx = ga*ax;
            gy = ga*ay;
            }

            else
            {
            gx = gb*bx;
            gy = gb*by;
            }

            Bx[c3(ii,jj,k)] = gx;
            By[c3(ii,jj,k)] = p->j_dir==1 ? gy : 0.0;
        }
    }
}
