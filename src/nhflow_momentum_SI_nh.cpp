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
// column; A x = F(x) - F(0), its columns probed with three colours (one cell in y) or five ((i + 2 j) mod 5, the global
// indices: the colours agree across ranks). One cell in y on one rank: band LU of the interleaved unknowns (bandwidth
// 2 (K+1)). Otherwise (2D, MPI): FGMRES(30) preconditioned with the multigrid of nhflow_si_nhmg, whose relaxation
// solves the K+1 unknowns of a column at once, as Basilisk's relax() of layered/nh.h; the coarse levels by aggregation,
// K-cycle (grid-independent iteration counts: 10 - 13 for 16^2 .. 64^2 columns at |r| <= 1e-6 |b|).
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

void nhflow_momentum_SI::nh_residual(lexer *p, fdm_nhf *d, ghostcell *pgc, const std::vector<double> &x, std::vector<double> &y)
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

    nh_faces(p,d);
    exchangeF(p,pgc,NQx);
    exchangeF(p,pgc,NQy);

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

void nhflow_momentum_SI::nh_solve(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double starttime_nh = pgc->timer();
    const int K = p->knoz, B = K+1;

    nh_prepare(p,d,pgc);

    int N=0;
    SLICELOOP4
    ++N;

    const int n2 = B*N;
    std::vector<double> x(n2), F0(n2), y(n2), e(n2);

    if(!nhlu)
    nhlu = std::make_shared<sflow_si_nh_lu>();

    auto &S = *nhlu;

    // F(0): the affine part; b = -F(0)
    std::fill(x.begin(),x.end(),0.0);
    nh_residual(p,d,pgc,x,F0);

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
    // the coupled system A x = -F(0): A assembled by probing (block 5-point over the columns, five colours, three with
    // one cell in y), FGMRES(30) right-preconditioned with one K-cycle of nhflow_si_nhmg, |r| <= 1e-6 |b|
    const int K = p->knoz, B = K+1;
    const int nx = p->knox, ny = p->knoy;
    const bool yd = (p->j_dir==1);
    const int nc = yd ? 5 : 3;
    const size_t nm = (size_t)nx*ny*B;

    if(!nhmg)
    {
    nhmg = std::make_shared<nhflow_si_nhmg>();
    const int nbr[4] = {p->nb1, p->nb4, yd ? p->nb3 : -2, yd ? p->nb2 : -2};
    nhmg->setup(pgc->mpi_comm,B,nx,ny,nbr);

        if(p->mpirank==0)
        cout<<"NHFLOW semi-implicit non-hydrostatic: multigrid of the coupled system with "<<nhmg->levels()<<" levels, "<<B<<" unknowns per column"<<endl;
    }

    nhflow_si_nhmg &M = *nhmg;

    // the cells of the residual (SLICELOOP4 order) in the box of the multigrid
    std::vector<int> cix, ciy;
    std::vector<int> cell((size_t)nx*ny,-1);
    {
    int n=0;
    SLICELOOP4
    {
    cix.push_back(i);
    ciy.push_back(j);
    cell[(size_t)i*ny + j] = n;
    ++n;
    }
    }

    const int N = (int)cix.size();
    auto colour = [&](int ii, int jj) { const int I = ii + p->origin_i, J = yd ? jj + p->origin_j : 0; return (((I + 2*J)%nc) + nc)%nc; };
    const int dix[5] = {0, 1,-1, 0, 0};
    const int diy[5] = {0, 0, 0, 1,-1};

    // ---- assembly
    M.zero();
    std::vector<double> e((size_t)N*B), y((size_t)N*B);

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

    // solid columns of the box: identity
    for(int ii=0; ii<nx; ++ii)
    for(int jj=0; jj<ny; ++jj)
    if(cell[(size_t)ii*ny+jj]<0)
    {
    double *a = M.block(ii,jj,0);
    for(int r=0; r<B; ++r)
    a[r*B+r] = 1.0;
    }

    M.finalize();

    // ---- FGMRES on the box vectors
    std::vector<double> b(nm,0.0), x(nm,0.0), w(nm), r0(nm);

    for(int n=0; n<N; ++n)
    {
    const size_t q = ((size_t)cix[n]*ny + ciy[n])*B;
    i = cix[n]; j = ciy[n];
        for(int t=0; t<B; ++t)
        b[q+t] = -F0[(size_t)B*n+t];

    // initial guess: eta^n and phi of the step before
    x[q] = eta0(i,j);
        for(k=0;k<K;++k)
        x[q+1+k] = nAct[IJ] ? PHI[IJK] : 0.0;
    }

    auto dot = [&](const std::vector<double> &u, const std::vector<double> &v)
    {
        double s = 0.0;
        for(size_t q=0; q<nm; ++q)
        s += u[q]*v[q];
        return pgc->globalsum(s);
    };

    const int m = 30, maxit = 300;
    const double tol = 1.0e-6;
    std::vector<std::vector<double>> V(m+1,std::vector<double>(nm)), Z(m,std::vector<double>(nm));
    std::vector<double> H((m+1)*m,0.0), cs(m,0.0), sn(m,0.0), g(m+1,0.0);

    const double bn = sqrt(dot(b,b));
    int it = 0;

    while(it<maxit)
    {
        M.matvec(x.data(),w.data());

        for(size_t q=0; q<nm; ++q)
        V[0][q] = b[q] - w[q];

        const double beta = sqrt(dot(V[0],V[0]));

        if(beta<=tol*MAX(bn,1.0e-30))
        break;

        for(size_t q=0; q<nm; ++q)
        V[0][q] /= beta;

        std::fill(g.begin(),g.end(),0.0);
        g[0] = beta;

        int jj=0;
        for(; jj<m && it<maxit; ++jj, ++it)
        {
            M.vcycle(V[jj].data(),Z[jj].data());
            M.matvec(Z[jj].data(),V[jj+1].data());

            for(int l=0; l<=jj; ++l)
            {
            const double h = dot(V[jj+1],V[l]);
            H[l*m+jj] = h;
            for(size_t q=0; q<nm; ++q)
            V[jj+1][q] -= h*V[l][q];
            }

            const double h1 = sqrt(dot(V[jj+1],V[jj+1]));
            H[(jj+1)*m+jj] = h1;

            if(h1>0.0)
            for(size_t q=0; q<nm; ++q)
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

            if(fabs(g[jj+1])<=tol*MAX(bn,1.0e-30))
            {
            ++jj; ++it;
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

        for(int l=0; l<jj; ++l)
        for(size_t q=0; q<nm; ++q)
        x[q] += yk[l]*Z[l][q];

        if(jj>0 && fabs(g[jj])<=tol*MAX(bn,1.0e-30))
        break;
    }

    nh_iter = it;

    if(it>=maxit && p->mpirank==0)
    cout<<"!!! NHFLOW semi-implicit non-hydrostatic: FGMRES not converged in "<<maxit<<" iterations, |r|/|b| "<<fabs(g[0])/MAX(bn,1.0e-30)<<" !!!"<<endl;

    for(int n=0; n<N; ++n)
    {
    const size_t q = ((size_t)cix[n]*ny + ciy[n])*B;
        for(int t=0; t<B; ++t)
        xr[(size_t)B*n+t] = x[q+t];
    }
}
