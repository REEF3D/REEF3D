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

// ---------------------------------------------------------------------------
// NHFLOW semi-implicit step (A 510 11/12/13) with the non-hydrostatic pressure (A 520 > 0), after the multilayer solver
// of Basilisk (Popinet 2020, layered/nh.h with the implicit free surface layered/implicit.h): the non-hydrostatic
// pressure phi = p/rho and the free surface eta^{n+1} come out of one coupled linear system per step.
//
// Layers k = 0 .. K-1 of the sigma grid, thickness h_k = WL dsig_k. phi lives on the layer interfaces: phi_k on the
// lower interface of layer k (phi_0 on the bed), phi_K = 0 at the surface (the K unknowns of a column, as the explicit
// NHFLOW, d->P). Per layer (layer-integrated, kinematic):
//
//   d_t (h u)_k = ... - g h_k grad(eta) - grad(h phi)_k + [phi grad(z)]_k      (the last two: pg_k, on the faces)
//   d_t (h w)_k = ...                                    - [phi]_k
//   div(h u)_k - [u . grad(z)]_k + [w]_k = 0                                    (incompressibility of layer k)
//
// with the Keller box scheme in the vertical: phi_k of a layer = (phi_{k-1/2} + phi_{k+1/2})/2, the interface
// velocities from w_{k+1/2} = 2 w_k - w_{k-1/2}, w_{-1/2} = u_0 . grad(z_b) (Basilisk takes 0 there: it holds on a
// flat bed only). The face term pg_k is Basilisk's hpg(): the conservative form of -grad(h phi) + [phi grad z] with the
// face differences of the interface heights (slopes limited to tan 30 deg, max_slope), times H_f/(mean h) of the face.
// [u . grad(z)]_k with the interface slopes of the cell (centred) and u_k on the upper, u_{k-1} on the lower interface;
// with A 505 1 (tl2) u^{n+1} estimated from the predictor and the acceleration of the step before (nh_prepare).
//
// Time levels: phi fully implicit, it acts on the velocities over dt; the constraint holds for the velocities of
// t^{n+1}, the layer discharges Q_k^{n+1} and w^{n+1} = w* - dt [phi]_k/h_k (Basilisk: those of t^{n+theta}; with them
// the vertical velocity of the Beji & Battjes bar developed an odd-even mode in the vertical at dt 0.01, theta 0.5,
// w^{n+1} = 2 w^{n+theta} - w*). The mass row moves WL with the transport of the step, T_k = theta Q_k^{n+1} +
// (1-theta) Q_k^n (eta_update). The hydrostatic part is the semi-implicit step of nhflow_momentum_SI
// (face_update: Q_k^{n+1} = Q~_k - dt g H_f (theta G^{n+1} + (1-theta) G^n) r_k), so eta^{n+1} and phi^{n+1} are
// coupled through the face discharges: the mass row of a column (eta_update) and the K constraint rows.
// Layer thicknesses of the face terms of the pressure: WL extrapolated to t^{n+1/2} (nh_levels); the geometry of the
// constraint (interface slopes, the depth that turns WH into w and the Keller box layer thickness) at t^{n+1}, WL
// extrapolated (nh_wl1), as the constraint holds for the velocities of t^{n+1}.
//
// After the solve: Q_k += dt pg_k/dsig_k on the faces (face_update), UH_k += the mean of that of the open faces of a
// cell (as the gravity term), WH_k -= dt [phi]_k/dsig_k (cell_update).
//
// Non-hydrostatic cells: wet, deeper than 10 A 544, p->deep, wet neighbours; phi = 0 elsewhere. Open boundary faces
// carry no pressure term.
//
// Solution (nh_solve): the operator is the residual F(x) of the full step with x = (eta, phi_0 .. phi_{K-1}) per
// column; A x = F(x) - F(0), A assembled directly from the derivatives of the face discharges (nh_dfaces, nh_assemble;
// probing F with three or five colours, as before, took about half of the time of an NH step and gives the same
// entries to round-off, environment NHCHECK). One cell in y on one rank: band LU of the interleaved unknowns (bandwidth
// 2 (K+1)). Otherwise (2D, MPI): FGMRES(30) preconditioned with the multigrid of nhflow_si_nhmg, whose relaxation
// solves the K+1 unknowns of a column at once, as Basilisk's relax() of layered/nh.h; the coarse levels by aggregation,
// K-cycle (grid-independent iteration counts: 10 - 13 for 16^2 .. 64^2 columns at |r| <= 1e-6 |b|; one pre- and
// post-sweep, the initial guess extrapolated in time, |r| <= A 506 |b| with A 506 = 1e-5 by default: 6 - 7).
// ---------------------------------------------------------------------------

#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sflow_momentum_SI_nh.h"
#include"nhflow_si_nhmg.h"
#include<vector>
#include<memory>
#include<algorithm>
#include<cmath>

namespace
{
// Basilisk's default limit of the interface slopes (hydro.h, tan 30 deg)
const double nh_max_slope = 0.577350269189626;

inline double nh_lim(double dz, double dx)
{
    const double s = dz/dx;
    return dx*(fabs(s)<nh_max_slope ? s : (s>0.0 ? nh_max_slope : -nh_max_slope));
}
}

void nhflow_momentum_SI::nh_ini(lexer *p)
{
    nh = (p->A520>0);

    if(!nh)
    return;

    nh_brk = 0.0;

    const int N = p->imax*p->jmax*(p->kmax+2);
    double **arr[] = {&NQx,&NQy,&WHS,&PHI,&PHIo,&NSX,&NSY};
    for(double **a : arr)
    p->Darray(*a,N);

    nWLh.assign((size_t)p->imax*p->jmax,0.0);
    nWLo.assign((size_t)p->imax*p->jmax,-1.0);
    nAct.assign((size_t)p->imax*p->jmax,0);
}

void nhflow_momentum_SI::nh_free()
{
    double **arr[] = {&NQx,&NQy,&WHS,&PHI,&PHIo,&NSX,&NSY};
    for(double **a : arr)
    if(*a!=nullptr)
    {
    delete [] *a;
    *a = nullptr;
    }
}

void nhflow_momentum_SI::nh_begin(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // time levels of the nonlinear terms of the step (as the semi-implicit SFLOW, nh_levels there), with the
    // non-hydrostatic pressure (A 520 > 0) and with the one-step predictor (A 510 11): with t^n the waves on the Beji &
    // Battjes bar steepened too much behind the bar, the more the larger dt
    //  - the face depths of the gravity term (H_f, face_depth) at t^{n+1/2}, extrapolated from t^n and t^{n-1}
    //  - omega (the vertical advection of the column operator) at t^{n+1/2}
    // The advecting transport of the predictor is that of the estimate pass of step_single (tl2)
    const double r = p->dt/MAX(p->dt_old,1.0e-20);
    const size_t ns = (size_t)p->imax*p->jmax;


    // the face depth between two wet non-shoreline cells second order in space, as Basilisk's face heights (hydro.h, a
    // variant of BCG): WL of the upwind cell extrapolated to the face and to t^{n+1/2}, WL_f = WL_up + a (1 - a c) g dx/2,
    // c = |u| dt/dx, g the central gradient of WL of the upwind cell
    {
    // wet deep cells; MPI: from the neighbouring ranks (flagslice4 marks only the first ghost column of a rank boundary)
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
    const int q=(ii-p->imin)*p->jmax + (jj-p->jmin);
    tS(ii,jj) = (p->flagslice4[q]>0 && p->wet[q]==1 && p->deep[q]==1) ? 1.0 : 0.0;
    }

    exchange4(p,pgc,tS);

    auto fsw = [&](int a, int c) { return tS(a,c)>0.5; };

        SLICELOOP1
        {
        if(Hx(i,j)>0.0 && fsw(i,j) && fsw(i+1,j) && fsw(i-1,j) && fsw(i+2,j))
        {
        const double u = qx(i,j)/Hx(i,j);
        const double dx = p->DXP[IP];
        const double c = u*p->dt/dx;
        const double a = c>=0.0 ? 1.0 : -1.0;
        const int iu = c>=0.0 ? i : i+1;
        const double g = 0.5*(d->WL(iu+1,j) - d->WL(iu-1,j));
        Hx(i,j) = MAX(d->WL(iu,j) + a*(1.0 - a*c)*0.5*g, p->A544);
        }
        }

        if(p->j_dir==1)
        SLICELOOP2
        if(Hy(i,j)>0.0 && fsw(i,j) && fsw(i,j+1) && fsw(i,j-1) && fsw(i,j+2))
        {
        const double u = qy(i,j)/Hy(i,j);
        const double dy = p->DYP[JP];
        const double c = u*p->dt/dy;
        const double a = c>=0.0 ? 1.0 : -1.0;
        const int ju = c>=0.0 ? j : j+1;
        const double g = 0.5*(d->WL(i,ju+1) - d->WL(i,ju-1));
        Hy(i,j) = MAX(d->WL(i,ju) + a*(1.0 - a*c)*0.5*g, p->A544);
        }
    }

    // omega of the step before (from its transport, t^{n-1+theta}) to t^{n+1/2}, for the Crank-Nicolson vertical
    // advection of the column operator (finish() computes it new at the end of the step): (1 - theta) dt_old + dt/2
    // ahead, e = (1 - theta) + r/2 (written so that theta 1/2 gives 0.5 (1 + r) exactly)
    const double e = 0.5*(1.0+r) + (0.5-theta);

    if(onestep || nh || tl2)
    {
        const int NF = p->imax*p->jmax*p->kmaxF;

        if(std::isnan(nOmP[0]))
        nOmP.assign(d->omegaF,d->omegaF+NF);

        for(int q=0; q<NF; ++q)
        {
        const double o = d->omegaF[q];
        d->omegaF[q] = o + e*(o - nOmP[q]);
        nOmP[q] = o;
        }
    }

    SLICELOOP1
    {
    const double h = Hx(i,j);
        if(h>0.0 && nHxo[IJ]>0.0)
        Hx(i,j) = MAX(h + 0.5*r*(h - nHxo[IJ]), p->A544);
    nHxo[IJ] = h;
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    const double h = Hy(i,j);
        if(h>0.0 && nHyo[IJ]>0.0)
        Hy(i,j) = MAX(h + 0.5*r*(h - nHyo[IJ]), p->A544);
    nHyo[IJ] = h;
    }

    exchange1(p,pgc,Hx);
    exchange2(p,pgc,Hy);
}

void nhflow_momentum_SI::nh_levels(lexer *p, fdm_nhf *d)
{
    // the layer thicknesses of the non-hydrostatic terms at t^{n+1/2}, extrapolated from t^n and t^{n-1}
    const double r = p->dt/MAX(p->dt_old,1.0e-20);

    SLICELOOP4
    {
    const double wo = nWLo[IJ];
    const bool ok = wo>0.0 && WL0(i,j)>p->A544;
    nWLh[IJ] = ok ? MAX(WL0(i,j) + 0.5*r*(WL0(i,j) - wo), p->A544) : WL0(i,j);
    nWLo[IJ] = WL0(i,j);
    }
}

double nhflow_momentum_SI::nh_wl1(lexer *p, int a, int c)
{
    // the depth of t^{n+1} for the constraint: WL^n + r (WL^n - WL^{n-1}) = 2 WLh - WL^n (ghost columns: as WLh)
    const int q = (a-p->imin)*p->jmax + (c-p->jmin);
    return MAX(2.0*nWLh[q] - WL0(a,c),p->A544);
}

void nhflow_momentum_SI::nh_prepare(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // non-hydrostatic cells, w* (the predictor with the implicit vertical advection, as cell_update), the slopes of
    // the interfaces of the cells (centred, limited) for [u . grad z] and the bed condition
    const int K = p->knoz;
    auto in = [&](int a, int c) { return p->flagslice4[(a-p->imin)*p->jmax + (c-p->jmin)]>0; };
    auto wt = [&](int a, int c) { const int q = (a-p->imin)*p->jmax + (c-p->jmin); return p->flagslice4[q]>0 && p->wet[q]==1; };
    auto dry = [&](int a, int c) { return in(a,c) && !wt(a,c); };

    nh_levels(p,d);

    // the layer thickness of the nh terms of the ghost cells: their WL (walls: mirrored)
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    if(!in(ii,jj))
    nWLh[(ii-p->imin)*p->jmax + (jj-p->jmin)] = MAX(d->WL(ii,jj),0.0);

    // MPI: nWLh of the columns of the neighbouring ranks (nh_faces, the slopes)
    if(p->mpi_size>1)
    {
        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
        tS(ii,jj) = nWLh[(ii-p->imin)*p->jmax + (jj-p->jmin)];

        exchange4(p,pgc,tS);

        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
        nWLh[(ii-p->imin)*p->jmax + (jj-p->jmin)] = tS(ii,jj);
    }

    SLICELOOP4
    {
    const bool a = p->wet[IJ]==1 && WL0(i,j)>10.0*p->A544 && p->deep[IJ]==1
                && !dry(i-1,j) && !dry(i+1,j) && (p->j_dir==0 || (!dry(i,j-1) && !dry(i,j+1)));
    nAct[IJ] = a ? 1 : 0;
    }

    // w*: the vertical advection of the step (column operator without viscosity), as cell_update
    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        for(k=0;k<K;++k)
        {
        wnu[k] = 0.0;
        wb[k] = WHp[IJK];
        }

        col_coef(p,d,WL0(i,j),i,j,i,j,wnu.data(),wcu.data(),wcl.data(),wcd.data());

        if(thv<1.0)
        {
        for(k=0;k<K;++k)
        wx[k] = WH0[IJK];
        col_adv_rhs(p,wx.data(),wb.data());
        }

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),0.0,wb.data(),wx.data());

        for(k=0;k<K;++k)
        WHS[IJK] = wx[k];
        }

        else
        for(k=0;k<K;++k)
        WHS[IJK] = 0.0;
    }

    // tl2: the velocities of [u . grad z] and of the bed condition at t^{n+1}, as the constraint holds for them: the cell
    // predictor plus the acceleration of the implicit part of the step before, UH* + dt WL^{n+1} acc (O(dt^2)); with the
    // predictor alone (u*, without the gravity and pressure of the step) the nonlinear non-hydrostatic step was first
    // order (standing wave a/h 0.1, A 510 11: order 1.35 at T/1280, the error of the first quarter period halved with dt)
    if(tl2)
    {
    const int N = p->imax*p->jmax*(p->kmax+2);
    nUZx.assign(N,0.0);
    nUZy.assign(N,0.0);

        SLICELOOP4
        if(p->wet[IJ]==1)
        {
        const double H1 = nh_wl1(p,i,j);
            for(k=0;k<K;++k)
            {
            nUZx[IJK] = UHD[IJK] + p->dt*H1*accU[IJK];
            nUZy[IJK] = VHD[IJK] + p->dt*H1*accV[IJK];
            }
        }
    }

    // interface heights z_m = z_b + sig_m WL1, m = 0 .. K; centred slopes of the cells (one-sided next to solids). The
    // constraint holds at t^{n+1}, so its geometry is that of t^{n+1}: WL1 = WL^n + r (WL^n - WL^{n-1}) = 2 WLh - WL^n
    // (with t^{n+1/2} the Stokes-5 wave ran 0.15 % slower)
    auto zif = [&](int a, int c, int m) { return d->bed(a,c) + p->ZN[m+marge]*nh_wl1(p,a,c); };

    SLICELOOP4
    for(int m=0; m<=K; ++m)
    {
    k = m;
    const int ie = in(i+1,j) ? i+1 : i, iw = in(i-1,j) ? i-1 : i;
    const double dx = p->XP[ie+marge] - p->XP[iw+marge];
    NSX[IJK] = dx>0.0 ? nh_lim(zif(ie,j,m) - zif(iw,j,m),dx)/dx : 0.0;
    NSY[IJK] = 0.0;

        if(p->j_dir==1)
        {
        const int jn = in(i,j+1) ? j+1 : j, js = in(i,j-1) ? j-1 : j;
        const double dy = p->YP[jn+marge] - p->YP[js+marge];
        NSY[IJK] = dy>0.0 ? nh_lim(zif(i,jn,m) - zif(i,js,m),dy)/dy : 0.0;
        }
    }
}

void nhflow_momentum_SI::nh_faces(lexer *p, fdm_nhf *d)
{
    // the layer discharges of phi on the faces, NQ_k = dt pg_k/dsig_k (per unit sigma, as Qx), pg_k Basilisk's hpg
    const double dt = p->dt;
    const int K = p->knoz;
    auto in = [&](int a, int c) { return p->flagslice4[(a-p->imin)*p->jmax + (c-p->jmin)]>0; };
    auto hh = [&](int a, int c) { return nWLh[(a-p->imin)*p->jmax + (c-p->jmin)]; };
    auto zif = [&](int a, int c, int m) { return d->bed(a,c) + p->ZN[m+marge]*hh(a,c); };

    // face of cells (iL,jL) - (iR,jR), distance dx, face depth Hf; NQ of its layers
    auto face = [&](int iL, int jL, int iR, int jR, double dx, double Hf, double *NQ)
    {
        const int il=i, jl=j;

        if(!(in(iL,jL) && in(iR,jR)) || Hf<=0.0)
        {
        for(k=0;k<K;++k)
        NQ[IJK] = 0.0;
        return;
        }

        const double hL = hh(iL,jL), hR = hh(iR,jR);
        const double c = (hL+hR)>1.0e-20 ? Hf/(dx*(hL + hR)) : 0.0;

        for(k=0;k<K;++k)
        {
        const double s  = p->DZN[KP];
        const double sl = nh_lim(zif(iR,jR,k) - zif(iL,jL,k),dx);
        const double su = nh_lim(zif(iR,jR,k+1) - zif(iL,jL,k+1),dx);

        i=iL; j=jL;
        const double fL = PHI[IJK], uL = (k<K-1) ? PHI[IJKp1] : 0.0;
        i=iR; j=jR;
        const double fR = PHI[IJK], uR = (k<K-1) ? PHI[IJKp1] : 0.0;
        i=il; j=jl;

        // (h_L - s_l) phi_L - (h_R + s_l) phi_R + (h_L + s_u) phi_L^up - (h_R - s_u) phi_R^up, the layer h = WL dsig
        const double br = (s*hL - sl)*fL - (s*hR + sl)*fR + (s*hL + su)*uL - (s*hR - su)*uR;

        NQ[IJK] = dt*br*c/s;
        }
    };

    SLICELOOP1
    face(i,j,i+1,j,p->DXP[IP],Hx(i,j),NQx);

    if(p->j_dir==1)
    SLICELOOP2
    face(i,j,i,j+1,p->DYP[JP],Hy(i,j),NQy);

    // open boundary faces: no pressure term
    for(const bc_face &f : bcf)
    {
    i=f.fi;
    j=f.fj;
    double *NQ = (f.cs==1 || f.cs==4) ? NQx : NQy;
        for(k=0;k<K;++k)
        NQ[IJK] = 0.0;
    }
}

void nhflow_momentum_SI::nh_residual(lexer *p, fdm_nhf *d, ghostcell *pgc, const std::vector<double> &x, std::vector<double> &y, bool zero)
{
    // F(x): the mass rows (eta_update) and the constraint rows of the step with x = (eta, phi_0 .. phi_{K-1}) per column
    const double dt = p->dt;
    const int K = p->knoz, B = K+1;
    const bool yd = (p->j_dir==1);

    int n=0;
    SLICELOOP4
    {
    etaI(i,j) = x[(size_t)B*n];

        for(k=0;k<K;++k)
        PHI[IJK] = nAct[IJ] ? x[(size_t)B*n+1+k] : 0.0;
    ++n;
    }

    exchange4(p,pgc,etaI);
    exchangeC(p,pgc,PHI);

    // x = 0: phi = 0, the phi faces zero (nh_faces gives 0 on every face, and so the exchange)
    if(zero)
    {
    const int NF = p->imax*p->jmax*(p->kmax+2);
    std::fill(NQx,NQx+NF,0.0);
    std::fill(NQy,NQy+NF,0.0);
    }

    else
    {
    nh_faces(p,d);
    exchangeF(p,pgc,NQx);
    exchangeF(p,pgc,NQy);
    }

    face_update(p,d,pgc);

    n=0;
    SLICELOOP4
    {
    // mass row: as eta_update, without the clipping at WL = 0
    double s=0.0;

    for(k=0;k<K;++k)
    s += p->DZN[KP]*layer_div(p,Tx,Ty);

    double w = WL0(i,j) - dt*s;

    if(mass_source)
    w += dt*wsrc(i,j);

    y[(size_t)B*n] = etaI(i,j) + d->depth(i,j) - w;

        if(!nAct[IJ])
        {
        for(k=0;k<K;++k)
        y[(size_t)B*n+1+k] = x[(size_t)B*n+1+k];
        }

        else
        {
        // the depth of t^{n+1} (extrapolated): w^{n+1} = (WH* - dt [phi]/dsig)/WL^{n+1} as nh_cell and velcalc form it
        const double H0 = nh_wl1(p,i,j), Hh = H0;
        const double iH0 = 1.0/H0;

        // the bed: w_b = u_0 . grad(z_b)
        k=0;
        const double *UU = tl2 ? nUZx.data() : UHD, *VV = tl2 ? nUZy.data() : VHD;
        double wlo = UU[IJK]*iH0*NSX[IJK] + (yd ? VV[IJK]*iH0*NSY[IJK] : 0.0);

        for(k=0;k<K;++k)
        {
        const double ds = p->DZN[KP];
        const double h = Hh*ds;
        const double ph_up = (k<K-1) ? PHI[IJKp1] : 0.0;

        // divergence of the layer discharges of t^{n+1}
        double D = ds*layer_div(p,Qx,Qy);

        // [u . grad z]_k: u_k on the upper interface (k+1), u_{k-1} on the lower one, u_0 on the bed
        const double uk = UU[IJK]*iH0, vk = yd ? VV[IJK]*iH0 : 0.0;
        double uz = uk*NSX[IJKp1] + vk*NSY[IJKp1];

            if(k>0)
            {
            const double um = UU[IJKm1]*iH0, vm = yd ? VV[IJKm1]*iH0 : 0.0;
            uz -= um*NSX[IJK] + vm*NSY[IJK];
            }

            else
            uz -= uk*NSX[IJK] + vk*NSY[IJK];

        D -= uz;

        // [w]_k, Keller box, with w^{n+1}
        const double wk = WHS[IJK]*iH0 - dt*(ph_up - PHI[IJK])/h;
        const double wup = 2.0*wk - wlo;
        D += wup - wlo;
        wlo = wup;

        // the scaling of the row (Basilisk's 2 h/(theta dt) with theta 1)
        y[(size_t)B*n+1+k] = D*2.0*h/dt;
        }
        }
    ++n;
    }
}

void nhflow_momentum_SI::nh_dfaces(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // the derivatives of the layer discharges of the faces of the step (face_update with nh_faces) with respect to the
    // unknowns: Q_k = QA_k - dt g H_f (theta G^{n+1} + (1-theta) G^n) R_k - K_f G^{n+1} + NQ_k, G^{n+1} = (eta_R - eta_L)/dn
    // + c, NQ_k = dt c_f/dsig_k ((h_L dsig_k - s_l) phi_L,k - (h_R dsig_k + s_l) phi_R,k + (h_L dsig_k + s_u) phi_L,k+1
    // - (h_R dsig_k - s_u) phi_R,k+1); open boundary faces as face_update (riemann/flather: inflow or outflow as at
    // eta = 0, as F(0))
    const double dt = p->dt;
    const int K = p->knoz;
    const int N = p->imax*p->jmax*(p->kmax+2);

    for(std::vector<double> *v : {&dQxL,&dQxR,&dQyL,&dQyR,&pQxL0,&pQxR0,&pQxL1,&pQxR1,&pQyL0,&pQyR0,&pQyL1,&pQyR1})
    v->assign(N,0.0);

    auto in = [&](int a, int c) { return p->flagslice4[(a-p->imin)*p->jmax + (c-p->jmin)]>0; };
    auto hh = [&](int a, int c) { return nWLh[(a-p->imin)*p->jmax + (c-p->jmin)]; };
    auto zif = [&](int a, int c, int m) { return d->bed(a,c) + p->ZN[m+marge]*hh(a,c); };

    // nAct with the columns of the neighbouring ranks
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    tS(ii,jj) = 0.0;

    SLICELOOP4
    tS(i,j) = nAct[IJ] ? 1.0 : 0.0;

    exchange4(p,pgc,tS);

    nActG.assign((size_t)p->imax*p->jmax,0.0);
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    nActG[(ii-p->imin)*p->jmax + (jj-p->jmin)] = tS(ii,jj);

    auto act = [&](int a, int c) { return nActG[(a-p->imin)*p->jmax + (c-p->jmin)]>0.5; };

    // face of the cells (iL,jL) - (iR,jR) with the face index (iL,jL)
    auto face = [&](int iL, int jL, int iR, int jR, double dn, double dxp, double Hf, double Kf, const double *R,
                    double *dL, double *dR, double *pL0, double *pR0, double *pL1, double *pR1)
    {
        const int il=i, jl=j;
        i=iL; j=jL;

        for(k=0;k<K;++k)
        {
        const double a = (dt*g*Hf*theta*R[IJK] + Kf)/dn;
        dL[IJK] = a;
        dR[IJK] = -a;
        }

        if(nh && in(iL,jL) && in(iR,jR))
        {
        const double hL = hh(iL,jL), hR = hh(iR,jR);
        const double c = (hL+hR)>1.0e-20 ? Hf/(dxp*(hL + hR)) : 0.0;

            for(k=0;k<K;++k)
            {
            const double sg = p->DZN[KP];
            const double sl = nh_lim(zif(iR,jR,k) - zif(iL,jL,k),dxp);
            const double su = nh_lim(zif(iR,jR,k+1) - zif(iL,jL,k+1),dxp);
            const double f = dt*c/sg;

            pL0[IJK] = act(iL,jL) ? f*(sg*hL - sl) : 0.0;
            pR0[IJK] = act(iR,jR) ? -f*(sg*hR + sl) : 0.0;
            pL1[IJK] = (k<K-1 && act(iL,jL)) ? f*(sg*hL + su) : 0.0;
            pR1[IJK] = (k<K-1 && act(iR,jR)) ? -f*(sg*hR - su) : 0.0;
            }
        }

        i=il; j=jl;
    };

    SLICELOOP1
    if(Hx(i,j)>0.0)
    face(i,j,i+1,j,fdn(p,1,i,j),p->DXP[IP],Hx(i,j),kx(i,j),Rx,dQxL.data(),dQxR.data(),pQxL0.data(),pQxR0.data(),pQxL1.data(),pQxR1.data());

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0)
    face(i,j,i,j+1,fdn(p,2,i,j),p->DYP[JP],Hy(i,j),ky(i,j),Ry,dQyL.data(),dQyR.data(),pQyL0.data(),pQyR0.data(),pQyL1.data(),pQyR1.data());

    // open boundary faces: no pressure term, the discharge as face_update
    for(const bc_face &f : bcf)
    {
        const bool xf = (f.cs==1 || f.cs==4);
        // the interior cell is the right one of the face for x- and y- (cs 1, 3), the left one for x+ and y+ (cs 4, 2)
        const bool intR = (f.cs==1 || f.cs==3);
        double *dL = xf ? dQxL.data() : dQyL.data();
        double *dR = xf ? dQxR.data() : dQyR.data();
        double *P4[4] = {xf ? pQxL0.data() : pQyL0.data(), xf ? pQxR0.data() : pQyR0.data(),
                         xf ? pQxL1.data() : pQyL1.data(), xf ? pQxR1.data() : pQyR1.data()};
        const double *R = xf ? Rx : Ry;
        const double Hf = xf ? Hx(f.fi,f.fj) : Hy(f.fi,f.fj);

        i=f.fi;
        j=f.fj;

        for(k=0;k<K;++k)
        {
        double a = 0.0;

            if(f.type==bc_level && Hf>0.0)
            a = -dt*g*Hf*theta*f.s/f.dxp*R[IJK];

            if(f.type==bc_riemann || f.type==bc_flather)
            {
            const bool inflow = f.A - f.B*eta0(f.i,f.j)>0.0;
            a = -f.s*f.B*(inflow ? f.Qg[k] : f.wo[k]);
            }

        dL[IJK] = intR ? 0.0 : a;
        dR[IJK] = intR ? a : 0.0;

            for(int q=0;q<4;++q)
            P4[q][IJK] = 0.0;
        }
    }

    for(std::vector<double> *v : {&dQxL,&dQxR,&dQyL,&dQyR,&pQxL0,&pQxR0,&pQxL1,&pQxR1,&pQyL0,&pQyR0,&pQyL1,&pQyR1})
    exchangeF(p,pgc,v->data());
}

void nhflow_momentum_SI::nh_assemble(lexer *p, const std::function<void(int,int,int,int,int,double)> &sink)
{
    // the rows of the coupled system of every column (as nh_residual, linear in x = (eta, phi_0 .. phi_{K-1})), as
    // blocks of the five columns: 0 the column, 1 +x, 2 -x, 3 +y, 4 -y; sink(i, j, block, row, unknown, value)
    const double dt = p->dt;
    const int K = p->knoz, B = K+1;
    const bool yd = (p->j_dir==1);
    std::vector<double> A((size_t)5*B*B);
    std::vector<double> wc(K), wl(K), dc(K);

    auto at = [&](int dd, int tr, int t) -> double& { return A[((size_t)dd*B + tr)*B + t]; };

    SLICELOOP4
    {
    std::fill(A.begin(),A.end(),0.0);

    // the derivatives of the layer discharge k of a face (face index fi,fj) into row tr with factor w: the left cell
    // of the face is block bL, the right one block bR
    auto addf = [&](int tr, double w, int fi, int fj, int bL, int bR, int kk, const std::vector<double> &dL,
                    const std::vector<double> &dR, const std::vector<double> &pL0, const std::vector<double> &pR0,
                    const std::vector<double> &pL1, const std::vector<double> &pR1)
    {
        const int il=i, jl=j, kl=k;
        i=fi; j=fj; k=kk;
        at(bL,tr,0) += w*dL[IJK];
        at(bR,tr,0) += w*dR[IJK];
        at(bL,tr,1+kk) += w*pL0[IJK];
        at(bR,tr,1+kk) += w*pR0[IJK];
            if(kk<K-1)
            {
            at(bL,tr,2+kk) += w*pL1[IJK];
            at(bR,tr,2+kk) += w*pR1[IJK];
            }
        i=il; j=jl; k=kl;
    };

    // a row with the divergence of the layer discharges k: east and west, north and south faces
    auto adddiv = [&](int tr, double w, int kk)
    {
        const double ax = w/p->DXN[IP];
        addf(tr, ax,i,  j,0,1,kk,dQxL,dQxR,pQxL0,pQxR0,pQxL1,pQxR1);
        addf(tr,-ax,i-1,j,2,0,kk,dQxL,dQxR,pQxL0,pQxR0,pQxL1,pQxR1);

            if(yd)
            {
            const double ay = w/p->DYN[JP];
            addf(tr, ay,i,j,  0,3,kk,dQyL,dQyR,pQyL0,pQyR0,pQyL1,pQyR1);
            addf(tr,-ay,i,j-1,4,0,kk,dQyL,dQyR,pQyL0,pQyR0,pQyL1,pQyR1);
            }
    };

    // mass row: eta + d - WL^n + dt sum_k dsig_k div(theta Q_k)
    at(0,0,0) += 1.0;
    for(int kk=0;kk<K;++kk)
    adddiv(0,dt*theta*p->DZN[kk+marge],kk);

        if(!nAct[IJ])
        {
        for(int kk=0;kk<K;++kk)
        at(0,1+kk,1+kk) = 1.0;
        }

        else
        {
        const double H1 = nh_wl1(p,i,j);
        std::fill(wl.begin(),wl.end(),0.0);

        // w-coefficients of the Keller box: w_k = WH*/H - dt (phi_{k+1} - phi_k)/h_k, D_k += 2 w_k - 2 w_lo,
        // w_lo' = 2 w_k - w_lo
            for(int kk=0;kk<K;++kk)
            {
            const double ds = p->DZN[kk+marge];
            const double h = H1*ds;
            const double sc = 2.0*h/dt;

            adddiv(1+kk,sc*ds,kk);

            std::fill(wc.begin(),wc.end(),0.0);
            wc[kk] += dt/h;
            if(kk<K-1)
            wc[kk+1] -= dt/h;

                for(int m=0;m<K;++m)
                {
                at(0,1+kk,1+m) += sc*(2.0*wc[m] - 2.0*wl[m]);
                wl[m] = 2.0*wc[m] - wl[m];
                }
            }
        }

        for(int dd=0;dd<5;++dd)
        {
        if(!yd && dd>2)
        continue;

            for(int tr=0;tr<B;++tr)
            for(int t=0;t<B;++t)
            {
            const double v = at(dd,tr,t);
            if(v!=0.0)
            sink(i,j,dd,tr,t,v);
            }
        }
    }
}

void nhflow_momentum_SI::nh_solve(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double starttime_nh = pgc->timer();
    const int K = p->knoz, B = K+1;

    nh_prepare(p,d,pgc);
    nh_fixed = true;

    int N=0;
    SLICELOOP4
    ++N;

    const int n2 = B*N;
    std::vector<double> &x = nsx, &F0 = nsF0, &y = nsy;
    x.assign(n2,0.0);
    F0.resize(n2);
    y.resize(n2);

    if(!nhlu)
    nhlu = std::make_shared<sflow_si_nh_lu>();

    auto &S = *nhlu;

    // F(0): the affine part; b = -F(0)
    nh_residual(p,d,pgc,x,F0,true);

    // ---- multigrid with the column relaxation (Basilisk's solver, nhflow_si_nhmg) inside FGMRES; one row of columns
    // on one rank: the band LU (direct, the cheapest there)
    if(p->j_dir==1 || p->mpi_size>1)
    {
    nh_solve_mg(p,d,pgc,F0,x);
    }

    else
    {
    std::vector<int> ci(N);
    {
    int n=0;
    SLICELOOP4
    ci[n++] = i;
    }

    std::vector<int> cell(p->imax,-1);
    for(int n=0; n<N; ++n)
    cell[ci[n]-p->imin] = n;

    S.kl = S.ku = 2*B-1;
    S.resize(n2);
    const bool probe = getenv("NHCHECK")!=nullptr;

    // direct assembly (nh_dfaces, nh_assemble); environment NHCHECK: also by probing nh_residual (3 B evaluations)
    // and the largest difference
    auto direct = [&](const std::function<void(int,int,double)> &put)
    {
        nh_dfaces(p,d,pgc);
        nh_assemble(p,[&](int ii, int jj, int dd, int tr, int t, double v)
        {
            const int n = cell[ii-p->imin];
            const int ij = ii + (dd==1 ? 1 : (dd==2 ? -1 : 0));
            if(n<0 || ij<p->imin || ij>=p->imin+p->imax)
            return;
            const int m = cell[ij-p->imin];
            if(m<0)
            return;
            const int rr = B*n+tr, cc = B*m+t;
            if(cc-rr>=-S.kl && cc-rr<=S.ku)
            put(rr,cc,v);
        });
    };

    if(!probe)
    direct([&](int rr, int cc, double v){ S.at(rr,cc) = v; });

    std::vector<double> e(probe ? n2 : 0);

    if(probe)
    for(int col=0; col<3; ++col)
    for(int t=0; t<B; ++t)
    {
        std::fill(e.begin(),e.end(),0.0);

        for(int n=0; n<N; ++n)
        if(((ci[n]%3)+3)%3==col)
        e[(size_t)B*n+t] = 1.0;

        nh_residual(p,d,pgc,e,y);

        for(int n=0; n<N; ++n)
        for(int tr=0; tr<B; ++tr)
        {
            const double v = y[(size_t)B*n+tr] - F0[(size_t)B*n+tr];

            if(v==0.0)
            continue;

            // the column of this colour within one cell
            for(int dd=-1; dd<=1; ++dd)
            {
            const int ii = ci[n]+dd;

            if(((ii%3)+3)%3!=col || ii<p->imin || ii>=p->imin+p->imax)
            continue;

            const int m = cell[ii-p->imin];

            if(m<0)
            continue;

            const int rr = B*n+tr, cc = B*m+t;

            if(cc-rr>=-S.kl && cc-rr<=S.ku)
            S.at(rr,cc) = v;
            }
        }
    }

    if(getenv("NHCHECK"))
    {
    double mx=0.0, dm=0.0;
    sflow_si_nh_lu D;
    D.kl = D.ku = S.kl;
    D.resize(n2);
    direct([&](int rr, int cc, double v){ D.at(rr,cc) = v; });
        for(int rr=0; rr<n2; ++rr)
        for(int cc=MAX(0,rr-S.kl); cc<=MIN(n2-1,rr+S.ku); ++cc)
        {
        mx = MAX(mx,fabs(S.at(rr,cc)));
        dm = MAX(dm,fabs(S.at(rr,cc) - D.at(rr,cc)));
        }
        if(p->mpirank==0)
        cout<<"NHCHECK band LU: max |probe - direct| "<<dm<<" of max |entry| "<<mx<<endl;
    }

    if(!S.factor())
    {
    if(p->mpirank==0)
    cout<<"!!! NHFLOW A 510 11/12/13 with A 520 "<<p->A520<<": the coupled eta-phi system is singular !!!"<<endl;
    exit(1);
    }

    for(int q=0; q<n2; ++q)
    x[q] = -F0[q];

    S.solve(x);
    nh_iter = 1;
    }

    // the solution: etaI, PHI, the face terms; the pressure for the output
    nh_residual(p,d,pgc,x,y);

    double res=0.0;
    for(int q=0; q<n2; ++q)
    res = MAX(res,fabs(y[q]));
    nh_res = res;

    FLOOP
    d->P[FIJK] = 0.0;

    SLICELOOP4
    for(k=0;k<K;++k)
    d->P[FIJK] = p->W1*PHI[IJK];

    nh_fixed = false;
    nh_time = pgc->timer() - starttime_nh;
}

void nhflow_momentum_SI::nh_cell(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // UH, VH: the mean of the face terms of the open faces (as the gravity term of cell_update); WH: -dt [phi]_k/dsig_k
    const double dt = p->dt;
    const int K = p->knoz;
    auto wn = [&](int a, int b) { return etaI(a,b) + d->depth(a,b) > p->A544; };

    SLICELOOP4
    if(p->wet[IJ]==1 && wn(i,j))
    {
    const double we = (Hx(i,j)>0.0 && wn(i+1,j)) ? 1.0 : 0.0;
    const double ww = (Hx(i-1,j)>0.0 && wn(i-1,j)) ? 1.0 : 0.0;
    const bool wz = tl2;
    const double wsx = MAX(we + ww + (wz ? wall(p,i,j,1) + wall(p,i,j,2) : 0.0),1.0);

    double ve=0.0, vw=0.0, wsy=1.0;

        if(p->j_dir==1)
        {
        ve = (Hy(i,j)>0.0 && wn(i,j+1)) ? 1.0 : 0.0;
        vw = (Hy(i,j-1)>0.0 && wn(i,j-1)) ? 1.0 : 0.0;
        wsy = MAX(ve + vw + (wz ? wall(p,i,j,3) + wall(p,i,j,4) : 0.0),1.0);
        }

        for(k=0;k<K;++k)
        {
        d->UH[IJK] += (ww*NQx[Im1JK] + we*NQx[IJK])/wsx;

        if(p->j_dir==1)
        d->VH[IJK] += (vw*NQy[IJm1K] + ve*NQy[IJK])/wsy;

            if(nAct[IJ])
            {
            const double ph_up = (k<K-1) ? PHI[IJKp1] : 0.0;
            d->WH[IJK] -= dt*(ph_up - PHI[IJK])/p->DZN[KP];
            }

            // Basilisk's breaking parameter: |w| <= b sqrt(g H)
            if(nh_brk>0.0)
            {
            const double H = MAX(etaI(i,j) + d->depth(i,j),0.0);
            const double wm = nh_brk*sqrt(g*H)*H;
            d->WH[IJK] = MAX(MIN(d->WH[IJK],wm),-wm);
            }
        }
    }
}

void nhflow_momentum_SI::nh_solve_mg(lexer *p, fdm_nhf *d, ghostcell *pgc, const std::vector<double> &F0, std::vector<double> &xr)
{
    // the coupled system A x = -F(0): A assembled from the face derivatives (block 5-point over the columns),
    // FGMRES(30) right-preconditioned with one K-cycle of nhflow_si_nhmg (single precision), stopping test below (A 506);
    // A 507 1: eta, then phi (two smaller systems, the same solver)
    const int K = p->knoz, B = K+1;
    const int nx = p->knox, ny = p->knoy;
    const bool yd = (p->j_dir==1);
    const int nc = yd ? 5 : 3;
    const size_t nm = (size_t)nx*ny*B;
    const int nbr[4] = {p->nb1, p->nb4, yd ? p->nb3 : -2, yd ? p->nb2 : -2};

    // the preconditioner is rebuilt at most every lagmax steps (below); 256^2 x 4: rebuilt every step 51.8 s, every 4th
    // 49.6 s, the iterations the same (5.1)
    const int lagmax = 4;

    if(!nhmg)
    {
    nhmg = std::make_shared<nhflow_si_nhmg>();
    nhmg->setup(pgc->mpi_comm,B,nx,ny,nbr);

        if(p->mpirank==0)
        cout<<"NHFLOW semi-implicit non-hydrostatic: multigrid of the coupled system with "<<nhmg->levels()<<" levels, "<<B<<" unknowns per column"<<endl;
    }

    nhflow_si_nhmg &M = *nhmg;

    // the cells of the residual (SLICELOOP4 order) in the box of the multigrid
    std::vector<int> cix, ciy;
    std::vector<int> cell((size_t)nx*ny,-1);
    int nact = 0;
    {
    int n=0;
    SLICELOOP4
    {
    cix.push_back(i);
    ciy.push_back(j);
    cell[(size_t)i*ny + j] = n;
    nact += nAct[IJ] ? 1 : 0;
    ++n;
    }
    }

    const int N = (int)cix.size();
    auto colour = [&](int ii, int jj) { const int I = ii + p->origin_i, J = yd ? jj + p->origin_j : 0; return (((I + 2*J)%nc) + nc)%nc; };
    const int dix[5] = {0, 1,-1, 0, 0};
    const int diy[5] = {0, 0, 0, 1,-1};

    // ---- assembly: directly from the face derivatives (nh_dfaces, nh_assemble); before: by probing nh_residual (B nc
    // evaluations of the full residual per step). Environment NHCHECK: both, and the largest difference of the entries
    // (round-off: 1e-13 of the largest entry)
    M.zero();
    const bool probe = getenv("NHCHECK")!=nullptr;
    std::vector<double> e(probe ? (size_t)N*B : 0), y(probe ? (size_t)N*B : 0);

    if(!probe)
    {
    nh_dfaces(p,d,pgc);
    nh_assemble(p,[&](int ii, int jj, int dd, int tr, int t, double v){ M.set(ii,jj,dd,tr,t,v); });
    }

    if(probe)
    {
    M.mark_all();

    for(int col=0; col<nc; ++col)
    for(int t=0; t<B; ++t)
    {
        std::fill(e.begin(),e.end(),0.0);

        for(int n=0; n<N; ++n)
        if(colour(cix[n],ciy[n])==col)
        e[(size_t)B*n+t] = 1.0;

        nh_residual(p,d,pgc,e,y);

        for(int n=0; n<N; ++n)
        for(int dd=0; dd<5; ++dd)
        {
            if(!yd && dd>2)
            continue;

            if(colour(cix[n]+dix[dd],ciy[n]+diy[dd])!=col)
            continue;

            double *a = M.block(cix[n],ciy[n],dd);
            for(int tr=0; tr<B; ++tr)
            a[tr*B+t] = y[(size_t)B*n+tr] - F0[(size_t)B*n+tr];
        }
    }
    }

    if(getenv("NHCHECK"))
    {
    nh_dfaces(p,d,pgc);
    double mx=0.0, dm=0.0;
    std::vector<double> Dg((size_t)nx*ny*5*B*B,0.0);
    nh_assemble(p,[&](int ii, int jj, int dd, int tr, int t, double v){ Dg[(((size_t)ii*ny + jj)*5 + dd)*B*B + tr*B + t] = v; });
        for(int ii=0; ii<nx; ++ii)
        for(int jj=0; jj<ny; ++jj)
        if(cell[(size_t)ii*ny+jj]>=0)
        for(int dd=0; dd<5; ++dd)
        for(int q=0; q<B*B; ++q)
        {
        const double a = M.block(ii,jj,dd)[q], b = Dg[(((size_t)ii*ny + jj)*5 + dd)*B*B + q];
        mx = MAX(mx,fabs(a));
        dm = MAX(dm,fabs(a-b));
        }
    mx = pgc->globalmax(mx);
    dm = pgc->globalmax(dm);
        if(p->mpirank==0)
        cout<<"NHCHECK multigrid blocks: max |probe - direct| "<<dm<<" of max |entry| "<<mx<<endl;
    }

    // solid columns of the box: identity
    for(int ii=0; ii<nx; ++ii)
    for(int jj=0; jj<ny; ++jj)
    if(cell[(size_t)ii*ny+jj]<0)
    for(int r=0; r<B; ++r)
    M.set(ii,jj,0,r,r,1.0);

    // rows in metres: the mass row is the error of eta; a constraint row k (D_k 2 h_k/dt) times dt^2/2 is D_k h_k dt, the
    // change of the layer thickness in the step that the error of the divergence makes. Block Gauss-Seidel and the
    // coarse levels do not change with the same scaling of every column; FGMRES and the stopping test measure metres
    const double rs = 0.5*p->dt*p->dt;

    for(int n=0; n<N; ++n)
    {
    i = cix[n]; j = ciy[n];

        if(nAct[IJ])
        for(int dd=0; dd<5; ++dd)
        {
        if(!yd && dd>2)
        continue;

        double *a = M.block(i,j,dd);
            for(int tr=1; tr<B; ++tr)
            for(int t=0; t<B; ++t)
            a[tr*B+t] *= rs;
        }
    }

    // ---- the right-hand side and the initial guess on the box vectors
    std::vector<double> &b = nmb, &x = nmx, &w = nmw;
    b.assign(nm,0.0);
    x.assign(nm,0.0);
    w.resize(nm);
    const bool guess = true;
    const bool hist = nEtaP.size()==(size_t)p->imax*p->jmax;

    for(int n=0; n<N; ++n)
    {
    const size_t q = ((size_t)cix[n]*ny + ciy[n])*B;
    i = cix[n]; j = ciy[n];
        for(int t=0; t<B; ++t)
        b[q+t] = -F0[(size_t)B*n+t]*(t>0 && nAct[IJ] ? rs : 1.0);

    // initial guess: eta^n and phi of the step before, both extrapolated in time (r = dt/dt_old)
    x[q] = eta0(i,j);
        for(k=0;k<K;++k)
        x[q+1+k] = nAct[IJ] ? PHI[IJK] : 0.0;

        if(guess && hist)
        {
        const double r = p->dt/MAX(p->dt_old,1.0e-20);
        x[q] += r*(eta0(i,j) - nEtaP[IJ]);
            for(k=0;k<K;++k)
            if(nAct[IJ])
            x[q+1+k] += r*(PHI[IJK] - nPhiP[IJK]);
        }
    }

    // the levels of this step for the guess of the next one
    {
    const int NC = p->imax*p->jmax*(p->kmax+2);
    nEtaP.resize((size_t)p->imax*p->jmax);
    nPhiP.resize(NC);
    SLICELOOP4
    {
    nEtaP[IJ] = eta0(i,j);
        for(k=0;k<K;++k)
        nPhiP[IJK] = PHI[IJK];
    }
    }

    // ---- FGMRES(30) of a system of the multigrid MM (n unknowns, Krylov vectors V, Z kept between the steps)
    // stopping test (A 506, default 1e-4): max |r| <= max(A 506 max |r0|, 1e-12 max(|b|, |A x0|)), r0 = b - A x0 the
    // residual of the guess, rows in metres (above). Relative to the step's own starting residual, not to |b|: b holds
    // eta^n, so the |b| test of v7b loosened with a tide or a still water level away from the initial one (0.4 m in
    // 0.45 m water: 8.0 -> 5.5 iterations, the error of the solve 7.8e-11 -> 9.8e-9 m). The floor is the round-off of
    // b - A x (state at rest: no iterations). FGMRES minimises the 2-norm: the inner target of a cycle is the 2-norm of
    // its starting residual times the reduction the max-norm still needs; the max-norm is checked on the true residual
    // of every restart. 2D standing wave 64 x 64 x 4, T/40 .. T/1280, error against T/5120: 1e-4 within 1.5 % of the
    // converged solve, 6.0 iterations (|b| test 1e-5: 8.0); 1e-3 up to 27 % at T/1280, 5.0 iterations
    struct fg_result { int it=0, cyc=0; double r0max=0.0, rmax=0.0, tolmax=0.0; };

    auto fgmres = [&](nhflow_si_nhmg &MM, const std::vector<double> &bb, std::vector<double> &xx,
                      std::vector<std::vector<double>> &V, std::vector<std::vector<double>> &Z, std::vector<double> &wv) -> fg_result
    {
        const size_t n = bb.size();
        const int m = 30, maxit = 300;
        const double tol = p->A506;
        fg_result R;

        auto dot = [&](const std::vector<double> &u, const std::vector<double> &v)
        {
            double s = 0.0;
            for(size_t q=0; q<n; ++q)
            s += u[q]*v[q];
            return pgc->globalsum(s);
        };

        // the Krylov vectors, kept from step to step and allocated when an iteration first needs them (allocating all
        // 61 of them every step cost about as much as two iterations)
        if(V.empty() || V[0].size()!=n)
        {
        V.assign(m+1,std::vector<double>());
        Z.assign(m,std::vector<double>());
        }

        auto need = [&](std::vector<double> &v) { if(v.size()!=n) v.assign(n,0.0); };
        need(V[0]);
        std::vector<double> H((m+1)*m,0.0), cs(m,0.0), sn(m,0.0), g(m+1,0.0), hv(m+2);
        wv.resize(n);

        while(true)
        {
            MM.matvec(xx.data(),wv.data());

            double mx[3] = {0.0, 0.0, 0.0};
            for(size_t q=0; q<n; ++q)
            {
            V[0][q] = bb[q] - wv[q];
            mx[0] = MAX(mx[0],fabs(V[0][q]));
            mx[1] = MAX(mx[1],fabs(bb[q]));
            mx[2] = MAX(mx[2],fabs(wv[q]));
            }

            pgc->globalmax(mx,3);
            R.rmax = mx[0];
            const double beta = sqrt(dot(V[0],V[0]));

            if(R.cyc==0)
            {
            R.r0max = R.rmax;
            R.tolmax = MAX(tol*R.r0max,1.0e-12*MAX(mx[1],mx[2]));
            }

            if(R.rmax<=R.tolmax || R.it>=maxit || beta<=0.0)
            break;

            ++R.cyc;
            const double target = 0.9*beta*R.tolmax/R.rmax;

            for(size_t q=0; q<n; ++q)
            V[0][q] /= beta;

            std::fill(g.begin(),g.end(),0.0);
            g[0] = beta;

            int jj=0;
            for(; jj<m && R.it<maxit; ++jj, ++R.it)
            {
                need(Z[jj]);
                need(V[jj+1]);
                MM.vcycle(V[jj].data(),Z[jj].data());
                MM.matvec(Z[jj].data(),V[jj+1].data());

                // classical Gram-Schmidt: all products in one pass and one reduction, the update in a second pass with
                // the new norm (before: modified Gram-Schmidt, a pass and a reduction per vector); again when the
                // norm drops below 1/sqrt(2) of the one before (DGKS)
                std::vector<double> &wj = V[jj+1];
                double h1 = 0.0;

                for(int pass=0; pass<2; ++pass)
                {
                    std::fill(hv.begin(),hv.begin()+jj+2,0.0);

                    for(size_t q=0; q<n; ++q)
                    {
                    const double wq = wj[q];
                    for(int l=0; l<=jj; ++l)
                    hv[l] += wq*V[l][q];
                    hv[jj+1] += wq*wq;
                    }

                    MPI_Allreduce(MPI_IN_PLACE,hv.data(),jj+2,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

                    double nn = 0.0;
                    for(size_t q=0; q<n; ++q)
                    {
                    double s = wj[q];
                    for(int l=0; l<=jj; ++l)
                    s -= hv[l]*V[l][q];
                    wj[q] = s;
                    nn += s*s;
                    }

                    nn = pgc->globalsum(nn);

                    for(int l=0; l<=jj; ++l)
                    H[l*m+jj] = (pass==0 ? 0.0 : H[l*m+jj]) + hv[l];

                    h1 = sqrt(MAX(nn,0.0));

                    if(h1*h1 > 0.5*hv[jj+1])
                    break;
                }

                H[(jj+1)*m+jj] = h1;

                if(h1>0.0)
                for(size_t q=0; q<n; ++q)
                V[jj+1][q] /= h1;

                for(int l=0; l<jj; ++l)
                {
                const double t0 = cs[l]*H[l*m+jj] + sn[l]*H[(l+1)*m+jj];
                H[(l+1)*m+jj] = -sn[l]*H[l*m+jj] + cs[l]*H[(l+1)*m+jj];
                H[l*m+jj] = t0;
                }

                const double a0 = H[jj*m+jj], a1 = H[(jj+1)*m+jj];
                const double rr = sqrt(a0*a0 + a1*a1);
                cs[jj] = rr>0.0 ? a0/rr : 1.0;
                sn[jj] = rr>0.0 ? a1/rr : 0.0;
                H[jj*m+jj] = rr;
                H[(jj+1)*m+jj] = 0.0;
                g[jj+1] = -sn[jj]*g[jj];
                g[jj]   =  cs[jj]*g[jj];

                if(fabs(g[jj+1])<=target)
                {
                ++jj; ++R.it;
                break;
                }
            }

            std::vector<double> yk(jj,0.0);
            for(int l=jj-1; l>=0; --l)
            {
            double s = g[l];
            for(int q=l+1; q<jj; ++q)
            s -= H[l*m+q]*yk[q];
            yk[l] = s/H[l*m+l];
            }

            for(size_t q=0; q<n; ++q)
            {
            double s = 0.0;
            for(int l=0; l<jj; ++l)
            s += yk[l]*Z[l][q];
            xx[q] += s;
            }
        }

        return R;
    };

    // ---- the preconditioner of an earlier step: rebuilt after lagmax steps, when dt changed by more than
    // 10 %, when the active columns changed, or when the iterations grew by more than one over those of the rebuild
    const bool refresh = lagmax<=1 || nh_lagc+1>=lagmax || fabs(p->dt/MAX(nh_lagdt,1.0e-20) - 1.0)>0.1
                         || pgc->globalisum(nact)!=nh_lagnact || nh_iter>nh_lagit+1;

    const int nsplit = p->A507;
    fg_result R, R2;

    if(nsplit==0)
    {
    M.finalize(refresh);
    R = fgmres(M,b,x,fgV,fgZ,nmv);
    }

    // ---- A 507 1: eta with phi of the guess (the eta system: the mass rows and eta), then phi with the new eta (the
    // constraint rows and phi): one block Gauss-Seidel pass over the coupled system, each part solved to A 506
    else
    {
    if(!nhme)
    {
    nhme = std::make_shared<nhflow_si_nhmg>();
    nhme->setup(pgc->mpi_comm,1,nx,ny,nbr);
    nhmp = std::make_shared<nhflow_si_nhmg>();
    nhmp->setup(pgc->mpi_comm,K,nx,ny,nbr);
    }

    nhflow_si_nhmg &Me = *nhme, &Mp = *nhmp;
    Me.zero();
    Mp.zero();

    for(int ii=0; ii<nx; ++ii)
    for(int jj=0; jj<ny; ++jj)
    for(int dd=0; dd<5; ++dd)
    {
        const double *a = M.block(ii,jj,dd);

        if(a[0]!=0.0)
        Me.set(ii,jj,dd,0,0,a[0]);

        for(int r=1; r<B; ++r)
        for(int c=1; c<B; ++c)
        if(a[r*B+c]!=0.0)
        Mp.set(ii,jj,dd,r-1,c-1,a[r*B+c]);
    }

    M.pattern();
    Me.finalize(refresh);
    Mp.finalize(refresh);

    const size_t nc1 = (size_t)nx*ny;
    std::vector<double> xs(nm,0.0), be(nc1), xe(nc1), bp(nc1*K), xp(nc1*K);

    // eta: b_eta - A_eta,phi phi
    for(size_t c=0; c<nc1; ++c)
    for(int t=1; t<B; ++t)
    xs[c*B+t] = x[c*B+t];

    M.matvec(xs.data(),w.data());

    for(size_t c=0; c<nc1; ++c)
    {
    be[c] = b[c*B] - w[c*B];
    xe[c] = x[c*B];
    }

    std::vector<double> wve, wvp;
    R = fgmres(Me,be,xe,fgVe,fgZe,wve);

    // phi: b_phi - A_phi,eta eta
    std::fill(xs.begin(),xs.end(),0.0);
    for(size_t c=0; c<nc1; ++c)
    xs[c*B] = xe[c];

    M.matvec(xs.data(),w.data());

    for(size_t c=0; c<nc1; ++c)
    for(int t=1; t<B; ++t)
    {
    bp[c*K+t-1] = b[c*B+t] - w[c*B+t];
    xp[c*K+t-1] = x[c*B+t];
    }

    R2 = fgmres(Mp,bp,xp,fgVp,fgZp,wvp);

    for(size_t c=0; c<nc1; ++c)
    {
    x[c*B] = xe[c];
        for(int t=1; t<B; ++t)
        x[c*B+t] = xp[c*K+t-1];
    }
    }

    if(refresh)
    {
    nh_lagc = 0;
    nh_lagdt = p->dt;
    nh_lagnact = pgc->globalisum(nact);
    nh_lagit = R.it + R2.it;
    }

    else
    ++nh_lagc;

    nh_iter = R.it + R2.it;

    if(getenv("NHLOG"))
    {
    // the residual of the coupled system (A 507 1: what the pass leaves)
    M.matvec(x.data(),w.data());
    double rc = 0.0;
    for(size_t q=0; q<nm; ++q)
    rc = MAX(rc,fabs(b[q]-w[q]));
    rc = pgc->globalmax(rc);

        if(p->mpirank==0)
        cout<<"NHLOG "<<p->count<<" it "<<R.it<<" + "<<R2.it<<" cycles "<<R.cyc<<" r0 "<<R.r0max<<" r "<<R.rmax<<" target "<<R.tolmax
            <<" coupled r "<<rc<<(refresh ? " refresh" : "")<<endl;
    }

    if((R.rmax>R.tolmax || R2.rmax>R2.tolmax) && p->mpirank==0)
    cout<<"!!! NHFLOW semi-implicit non-hydrostatic: FGMRES not converged in 300 iterations, max |r| "<<MAX(R.rmax,R2.rmax)<<" m, target "<<R.tolmax<<" m !!!"<<endl;

    for(int n=0; n<N; ++n)
    {
    const size_t q = ((size_t)cix[n]*ny + ciy[n])*B;
        for(int t=0; t<B; ++t)
        xr[(size_t)B*n+t] = x[q+t];
    }
}
