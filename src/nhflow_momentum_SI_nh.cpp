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
// of Popinet 2020 (non-hydrostatic, with the implicit free surface): the non-hydrostatic
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
// velocities from w_{k+1/2} = 2 w_k - w_{k-1/2}, w_{-1/2} = u_0 . grad(z_b) (Popinet 2020 takes 0 there: it holds
// on a flat bed only). The face term pg_k: the conservative form of -grad(h phi) + [phi grad z] with the
// face differences of the interface heights (slopes limited to tan 30 deg, max_slope), times H_f/(mean h) of the face.
// [u . grad(z)]_k with the interface slopes of the cell (centred) and u_k on the upper, u_{k-1} on the lower interface;
// with tl2 u^{n+1} estimated from the predictor and the acceleration of the step before (nh_prepare).
//
// Time levels: phi fully implicit, it acts on the velocities over dt; the constraint holds for the velocities of
// t^{n+1}, the layer discharges Q_k^{n+1} and w^{n+1} = w* - dt [phi]_k/h_k (Popinet 2020: those of t^{n+theta}; with them
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
// [u . grad z]_k with the interface slopes of the geometry of t^{n+1} from the solve itself: z_m = z_b + sig_m (WL^n +
// eta^{n+1} - eta^n), centred (nh_prepare: the part of t^n, NSX/NSY, and the cells whose eta enters; nh_residual,
// nh_assemble: the change of eta, implicit). Before, the geometry was extrapolated to t^{n+1} (2 WLh - WL^n): in a
// current the waves of 5 - 6 cells grew by 6 - 8 % per step from a gravity Courant number of about 1 (2D flume, h 0.109
// m, U 0.42 m/s, dx 1 cm, dt 0.012 s; one layer or ten). A linear (von Neumann) analysis of the one-layer step (Keller
// box, the two passes of tl2, the multilayer flux, theta) gives the same growth and locates it in this term: the
// extrapolation of an oscillation of omega_abs dt ~ 1 over a step. With eta^{n+1} the step is stable there up to the
// advective limit of the predictor (Courant 0.6 with SSP-RK2 and WENO5, 0.84 with SSP-RK3 tested), for theta 0.55
// and 0.5. A slope that already exceeds the limit at t^n (tan 30 deg, nh_lim) stays at the limit, without eta^{n+1}, so
// that the rows stay linear.
//
// Solution (nh_solve): the operator is the residual F(x) of the full step with x = (eta, phi_0 .. phi_{K-1}) per
// column; A x = F(x) - F(0), A assembled directly from the derivatives of the face discharges (nh_dfaces, nh_assemble;
// checked against probing F with three or five colours with the environment variable NHCHECK). One cell in y on one rank: band LU of the interleaved unknowns (bandwidth
// 2 (K+1)). Otherwise (2D, MPI): FGMRES(30) preconditioned with the multigrid of nhflow_si_nhmg, whose relaxation
// solves the K+1 unknowns of a column at once, as in Popinet 2020; the coarse levels by aggregation,
// K-cycle (grid-independent iteration counts: 10 - 13 for 16^2 .. 64^2 columns at |r| <= 1e-6 |b|; one pre- and
// post-sweep, the initial guess extrapolated in time, |r| <= A 506 |b| with A 506 = 1e-5 by default: 6 - 7).
// ---------------------------------------------------------------------------

#include"nhflow_momentum_SI.h"
#include"gcsl_yskip.h"
#include"sflow_curv.h"
#include"si_eta2D.h"
#include"nhflow_drywall.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"slide2D.h"
#include"ghostcell.h"
#include"sflow_momentum_SI_nh.h"
#include"nhflow_si_nhmg.h"
#include<vector>
#include<memory>
#include<algorithm>
#include<cmath>

namespace
{
// the limit of the interface slopes (tan 30 deg, as Popinet 2020)
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

    p->Darray_lazy(NSFx,lN);
    p->Darray_lazy(NSFy,lN);

    nh_brk = p->A558;   // A 558: the breaking limiter |w| <= A 558 sqrt(g H) of Popinet 2020 (nh_cell)

    const int N = lN;
    double **arr[] = {&NQx,&NQy,&WHS,&PHI,&PHIo,&NSX,&NSY};
    for(double **a : arr)
    p->Darray_lazy(*a,N);

    nWLh.assign((size_t)p->imax*p->jmax,0.0);
    nWLo.assign((size_t)p->imax*p->jmax,-1.0);
    nAct.assign((size_t)p->imax*p->jmax,0);
    nPw.assign((size_t)p->imax*p->jmax,1.0);

    // one curvilinear grid (nhflow_momentum_SI_nhcurv.cpp)
    nhcv = (p->geo_curv==1);
    nhco = nhcv && p->A539==1;

    if(nhcv)
    {
    nCL1.assign((size_t)p->imax*p->jmax,0.0);
    nCL2.assign((size_t)p->imax*p->jmax,0.0);
    nCT.assign((size_t)4*p->imax*p->jmax,0.0);
    }
}

void nhflow_momentum_SI::nh_free()
{
    double **arr[] = {&NQx,&NQy,&WHS,&PHI,&PHIo,&NSX,&NSY,&NSFx,&NSFy};
    for(double **a : arr)
    if(*a!=nullptr)
    {
    free(*a);
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

    // the face depth of the mass flux (hMx, hMy; si_assemble) between two wet non-shoreline cells, second order in
    // space (BCG): WL of the upwind cell extrapolated to the face and to t^{n+1/2}, WL_f = WL_up + a (1 - a c) g dx/2,
    // c = |u| dt/dx, g the slope of WL of the upwind cell, limited: van Leer, Superbee or van Albada (A 514 1, 2, 3),
    // minmod otherwise. The face depth of the gravity term (Hx, Hy) stays the mean of the two cells (mlflux_faces)
    {
    auto glim = [&](double sl, double sr)
    {
        if(sl*sr<=0.0)
        return 0.0;

        const double a = fabs(sl), b = fabs(sr), sg = sl>0.0 ? 1.0 : -1.0;

        if(p->A514==1)
        return 2.0*sl*sr/(sl + sr);

        if(p->A514==2)
        return sg*MAX(MIN(2.0*a,b),MIN(a,2.0*b));

        if(p->A514==3)
        return sl*sr*(sl + sr)/(sl*sl + sr*sr);

        return sg*MIN(a,b);
    };

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
        if(Hx(i,j)>0.0 && fsw(i,j) && fsw(i+1,j) && fsw(i-1,j) && fsw(i+2,j))
        {
        const double u = qx(i,j)/Hx(i,j);
        const double dx = p->geo_curv==1 ? MAX(sflow_curv::dnorm(p,1,i,j),1.0e-12) : p->DXP[IP];   // curvilinear: along the normal
        const double c = u*p->dt/dx;
        const double a = c>=0.0 ? 1.0 : -1.0;
        const int iu = c>=0.0 ? i : i+1;
        const double g = glim(d->WL(iu,j) - d->WL(iu-1,j),d->WL(iu+1,j) - d->WL(iu,j));
        hMx[IJ] = MAX(d->WL(iu,j) + a*(1.0 - a*c)*0.5*g, p->A544);
        }

        if(p->j_dir==1)
        SLICELOOP2
        if(Hy(i,j)>0.0 && fsw(i,j) && fsw(i,j+1) && fsw(i,j-1) && fsw(i,j+2))
        {
        const double u = qy(i,j)/Hy(i,j);
        const double dy = p->geo_curv==1 ? MAX(sflow_curv::dnorm(p,2,i,j),1.0e-12) : p->DYP[JP];
        const double c = u*p->dt/dy;
        const double a = c>=0.0 ? 1.0 : -1.0;
        const int ju = c>=0.0 ? j : j+1;
        const double g = glim(d->WL(i,ju) - d->WL(i,ju-1),d->WL(i,ju+1) - d->WL(i,ju));
        hMy[IJ] = MAX(d->WL(i,ju) + a*(1.0 - a*c)*0.5*g, p->A544);
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
    nWLh[(ii-p->imin)*p->jmax + (jj-p->jmin)] = MAX(d->WL(ii,gcsl_yrow(p,jj)),0.0);

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
    // dry land next to deep water (A 541, nhflow_drywall.h): the column stays non-hydrostatic, its dry faces carry
    // no discharge
    const bool a = p->wet[IJ]==1 && WL0(i,j)>10.0*p->A544 && p->deep[IJ]==1
                && (nhflow_drywall::cliff(p,d,i,j,WL0(i,j))
                    || (!dry(i-1,j) && !dry(i+1,j) && (p->j_dir==0 || (!dry(i,j-1) && !dry(i,j+1)))));
    nAct[IJ] = (a && !df_col(p,i,j)) ? 1 : 0;
    }

    // w*: the vertical advection of the step (column operator without viscosity), as cell_update
    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        for(k=0;k<K;++k)
        {
        wnu[k] = 0.0;
        wb[k] = WHp[PIJK];
        }

        if(ld)
        ld_cell(p,2);
        col_coef(p,d,WL0(i,j),i,j,i,j,wnu.data(),wcu.data(),wcl.data(),wcd.data());

        if(thv<1.0)
        {
        for(k=0;k<K;++k)
        wx[k] = WH0[PIJK];
        col_adv_rhs(p,wx.data(),wb.data());
        }

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),0.0,wb.data(),wx.data());

        for(k=0;k<K;++k)
        WHS[PIJK] = wx[k];
        }

        else
        for(k=0;k<K;++k)
        WHS[PIJK] = 0.0;
    }

    // tl2: the velocities of [u . grad z] and of the bed condition at t^{n+1}, as the constraint holds for them: the cell
    // predictor plus the acceleration of the implicit part of the step before, UH* + dt WL^{n+1} acc (O(dt^2)); with the
    // predictor alone (u*, without the gravity and pressure of the step) the nonlinear non-hydrostatic step was first
    // order (standing wave a/h 0.1, A 510 11: order 1.35 at T/1280, the error of the first quarter period halved with dt)
    if(tl2)
    {
    nUZx.assign(lN,0.0);
    nUZy.assign(lN,0.0);

        SLICELOOP4
        if(p->wet[IJ]==1)
        {
        const double H1 = nh_wl1(p,i,j);
            for(k=0;k<K;++k)
            {
            nUZx[PIJK] = UHD[PIJK] + p->dt*H1*accU[PIJK];
            nUZy[PIJK] = VHD[PIJK] + p->dt*H1*accV[PIJK];
            }
        }
    }

    // interface heights z_m = z_b + sig_m WL, m = 0 .. K; centred slopes of the cells (one-sided next to solids). The
    // constraint holds at t^{n+1}, so its geometry is that of t^{n+1}: WL^n + (eta^{n+1} - eta^n), the change of eta
    // from the solve (nh_residual, nh_assemble; see the top of the file). Here the part of t^n (NSX, NSY) and the
    // cells whose eta enters (nSx, nSy: the east and west, north and south cell, -9 none: an open boundary ghost or the
    // column itself one-sided on both sides; NSFx, NSFy 0 where the slope stays at its limit)
    auto zif = [&](int a, int c, int m) { return d->bed(a,c) + p->ZN[m+marge]*WL0(a,c); };

    nSx.assign((size_t)2*p->imax*p->jmax,-9);
    nSy.assign((size_t)2*p->imax*p->jmax,-9);
    nSdx.assign((size_t)p->imax*p->jmax,0.0);
    nSdy.assign((size_t)p->imax*p->jmax,0.0);

    // a side with an open boundary face (inflow, outflow, open edge): centred as in the interior, the ghost column
    // with the boundary state of the step (WLb of bc_faces: the inflow level, the outflow level, the open edge) on the
    // bed of the boundary cell. One-sided, the slope of the boundary column was twice that of its faces (the face
    // towards the ghost has the level of the boundary, as the face discharge): [u . grad z] of the constraint drove a
    // downward w in the column behind a fixed-discharge inflow, which the phi of the step held up as a steady jump of
    // the level at the inflow (abutment flume, 16 l/s, dx 1 cm: +21 mm in the first column, a lee wave behind it)
    wall_sides(p);

    auto zob = [&](int a, int c, int m) { return d->bed(i,j) + p->ZN[m+marge]*MAX(WLb(a,c),0.0); };

    SLICELOOP4
    for(int m=0; m<=K; ++m)
    {
    k = m;
    // (a column in front of dry land, A 541: one-sided next to the dry cells too, their level is that of the land)
    const bool cl = nhflow_drywall::cliff(p,d,i,j,WL0(i,j));
    auto ok = [&](int a, int c) { return cl ? wt(a,c) : in(a,c); };
    const unsigned char ws = wSide[(i-p->imin)*p->jmax + (j-p->jmin)];

    const bool oe = !ok(i+1,j) && (ws&2), ow = !ok(i-1,j) && (ws&1);
    const int ie = (ok(i+1,j) || oe) ? i+1 : i, iw = (ok(i-1,j) || ow) ? i-1 : i;
    // curvilinear grid: the centre distances along the index directions, the velocity transform (once per cell)
    if(nhcv && m==0)
    {
    const bool on_ = !ok(i,j+1) && (ws&8), os_ = !ok(i,j-1) && (ws&4);
    nh_curv_metric(p,ie,iw,oe,ow,(ok(i,j+1) || on_) ? j+1 : j,(ok(i,j-1) || os_) ? j-1 : j,on_,os_);
    }

    const double dx = nhcv ? nCL1[(i-p->imin)*p->jmax + (j-p->jmin)] : p->XP[ie+marge] - p->XP[iw+marge];
    const double ze = oe ? zob(ie,j,m) : zif(ie,j,m), zw = ow ? zob(iw,j,m) : zif(iw,j,m);
    NSX[PIJK] = dx>0.0 ? nh_lim(ze - zw,dx)/dx : 0.0;
    NSFx[PIJK] = (dx>0.0 && fabs(ze - zw)<nh_max_slope*dx) ? 1.0 : 0.0;
    NSY[PIJK] = 0.0;
    NSFy[PIJK] = 0.0;

        if(m==0)
        {
        const size_t q = (size_t)((i-p->imin)*p->jmax + (j-p->jmin));
        nSx[2*q]   = oe ? -9 : ie-i;
        nSx[2*q+1] = ow ? -9 : iw-i;
        nSdx[q] = dx;
        }

        if(p->j_dir==1)
        {
        const bool on = !ok(i,j+1) && (ws&8), os = !ok(i,j-1) && (ws&4);
        const int jn = (ok(i,j+1) || on) ? j+1 : j, js = (ok(i,j-1) || os) ? j-1 : j;
        const double dy = nhcv ? nCL2[(i-p->imin)*p->jmax + (j-p->jmin)] : p->YP[jn+marge] - p->YP[js+marge];
        const double zn = on ? zob(i,jn,m) : zif(i,jn,m), zs = os ? zob(i,js,m) : zif(i,js,m);
        NSY[PIJK] = dy>0.0 ? nh_lim(zn - zs,dy)/dy : 0.0;
        NSFy[PIJK] = (dy>0.0 && fabs(zn - zs)<nh_max_slope*dy) ? 1.0 : 0.0;

            if(m==0)
            {
            const size_t q = (size_t)((i-p->imin)*p->jmax + (j-p->jmin));
            nSy[2*q]   = on ? -9 : jn-j;
            nSy[2*q+1] = os ? -9 : js-j;
            nSdy[q] = dy;
            }
        }
    }

    // curvilinear grid: the velocities of [u . grad z] along the index directions
    if(nhcv)
    nh_curv_vel(p);
}

void nhflow_momentum_SI::nh_faces(lexer *p, fdm_nhf *d)
{
    // the layer discharges of phi on the faces, NQ_k = dt pg_k/dsig_k (per unit sigma, as Qx), pg_k the face term of the header
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
        NQ[PIJK] = 0.0;
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
        const double fL = PHI[PIJK], uL = (k<K-1) ? PHI[PIJKp1] : 0.0;
        i=iR; j=jR;
        const double fR = PHI[PIJK], uR = (k<K-1) ? PHI[PIJKp1] : 0.0;
        i=il; j=jl;

        // (h_L - s_l) phi_L - (h_R + s_l) phi_R + (h_L + s_u) phi_L^up - (h_R - s_u) phi_R^up, the layer h = WL dsig
        const double br = (s*hL - sl)*fL - (s*hR + sl)*fR + (s*hL + su)*uL - (s*hR - su)*uR;

        NQ[PIJK] = dt*br*c/s;

        // VRANS: the inertia factor and the drag of the face layer
        if(vr)
        NQ[PIJK] *= (NQ==NQx ? vrSx : vrSy)[PIJK];
        }
    };

    SLICELOOP1
    face(i,j,i+1,j,nhcv ? fdn(p,1,i,j) : p->DXP[IP],Hx(i,j),NQx);

    if(p->j_dir==1)
    SLICELOOP2
    face(i,j,i,j+1,nhcv ? fdn(p,2,i,j) : p->DYP[JP],Hy(i,j),NQy);

    // open boundary faces: no pressure term
    for(const bc_face &f : bcf)
    {
    i=f.fi;
    j=f.fj;
    double *NQ = (f.cs==1 || f.cs==4) ? NQx : NQy;
        for(k=0;k<K;++k)
        NQ[PIJK] = 0.0;
    }
}

void nhflow_momentum_SI::nh_residual(lexer *p, fdm_nhf *d, ghostcell *pgc, const std::vector<double> &x, std::vector<double> &y, bool zero, bool rows)
{
    // F(x): the mass rows (eta_update) and the constraint rows of the step with x = (eta, phi_0 .. phi_{K-1}) per column;
    // rows false: only the state of x (etaI, PHI, the face terms and discharges), for the solution
    const double dt = p->dt;
    const int K = p->knoz, B = K+1;
    const bool yd = (p->j_dir==1);

    int n=0;
    SLICELOOP4
    {
    etaI(i,j) = x[(size_t)B*n];

        for(k=0;k<K;++k)
        PHI[PIJK] = nAct[IJ] ? x[(size_t)B*n+1+k] : 0.0;
    ++n;
    }

    exchange4(p,pgc,etaI);
    pexchC(p,pgc,PHI);

    // x = 0: phi = 0, the phi faces zero (nh_faces gives 0 on every face, and so the exchange)
    if(zero)
    {
    std::fill(NQx,NQx+lN,0.0);
    std::fill(NQy,NQy+lN,0.0);
    }

    else
    {
    nh_faces(p,d);
    pexchF(p,pgc,NQx);
    pexchF(p,pgc,NQy);
    }

    if(ho || nhco)
    ho_phi_add(p);

    face_update(p,d,pgc);

    if(!rows)
    return;

    n=0;
    SLICELOOP4
    {
    // mass row: as eta_update, without the clipping at WL = 0
    double s=0.0;

    for(k=0;k<K;++k)
    s += p->DZN[KP]*layer_div_p(p,Tx,Ty);

    double w = WL0(i,j) - dt*s;

    if(mass_source)
    w += dt*wsrc(i,j);

    y[(size_t)B*n] = nPw[IJ]*(etaI(i,j) + d->depth(i,j)) - w;

    // VRANS: the storage n_s
    if(vr)
    y[(size_t)B*n] += (d->PORS(i,j) - 1.0)*(nPw[IJ]*(etaI(i,j) + d->depth(i,j)) - WL0(i,j));

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

        // the interface slopes of the geometry of t^{n+1}: NSX + sig_m (de_e - de_w)/dx, de = eta^{n+1} - eta^n
        const size_t qc = (size_t)((i-p->imin)*p->jmax + (j-p->jmin));
        auto de = [&](int di, int dj) { return etaI(i+di,j+dj) - eta0(i+di,j+dj); };
        const double gx = nSdx[qc]>0.0 ? ((nSx[2*qc]>-9 ? de(nSx[2*qc],0) : 0.0) - (nSx[2*qc+1]>-9 ? de(nSx[2*qc+1],0) : 0.0))/nSdx[qc] : 0.0;
        const double gy = (yd && nSdy[qc]>0.0) ? ((nSy[2*qc]>-9 ? de(0,nSy[2*qc]) : 0.0) - (nSy[2*qc+1]>-9 ? de(0,nSy[2*qc+1]) : 0.0))/nSdy[qc] : 0.0;
        auto SX = [&](int m) { k=m; return NSX[PIJK] + NSFx[PIJK]*p->ZN[m+marge]*gx; };
        auto SY = [&](int m) { k=m; return NSY[PIJK] + NSFy[PIJK]*p->ZN[m+marge]*gy; };

        // the bed: w_b = u_0 . grad(z_b)
        k=0;
        const double *UU = (tl2 || nhcv) ? nUZx.data() : UHD, *VV = (tl2 || nhcv) ? nUZy.data() : VHD;
        double wlo = UU[PIJK]*iH0*SX(0) + (yd ? VV[PIJK]*iH0*SY(0) : 0.0);

        // the slide layer (W 150): the bed moves, w_b = dz_b/dt + u_0 . grad(z_b)
        if(d->slide!=nullptr)
        wlo += d->slide->wb(i,j);
        k=0;

        for(k=0;k<K;++k)
        {
        const double ds = p->DZN[KP];
        const double h = Hh*ds;
        const double ph_up = (k<K-1) ? PHI[PIJKp1] : 0.0;

        // divergence of the layer discharges of t^{n+1}
        double D = ds*layer_div_p(p,Qx,Qy);

        // [u . grad z]_k: u_k on the upper interface (k+1), u_{k-1} on the lower one, u_0 on the bed
        const int kc = k;
        const double uk = UU[PIJK]*iH0, vk = yd ? VV[PIJK]*iH0 : 0.0;
        const double um = kc>0 ? UU[PIJKm1]*iH0 : uk, vm = kc>0 ? (yd ? VV[PIJKm1]*iH0 : 0.0) : vk;
        double uz = uk*SX(kc+1) + (yd ? vk*SY(kc+1) : 0.0);
        uz -= um*SX(kc) + (yd ? vm*SY(kc) : 0.0);
        k=kc;

        D -= uz;

        // [w]_k, Keller box, with w^{n+1}; 6DOF: in a cell of a body the body velocity (fb_wcell)
        const int fbc = fbon ? fb_wcell(p) : 0;
        const double dtw = vr ? dt*vr_sw(p) : dt;    // VRANS: the inertia factor and the drag of w
        const double wk = fbc>0 ? fbW[IJK] : WHS[PIJK]*iH0 - dtw*(ph_up - PHI[PIJK])/h;
        const double wup = 2.0*wk - wlo;
        D += wup - wlo;
        wlo = wup;

        // the scaling of the row (2 h/(theta dt) with theta 1, as Popinet 2020); 6DOF: a cell of a body above another
        // one, phi_k = phi_{k-1}
        y[(size_t)B*n+1+k] = fbc==2 ? PHI[PIJK] - PHI[PIJKm1] : D*2.0*h/dt;
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

    // the layers of the columns only (index dqi: no vertical ghost cells, a third of a 3D array with 4 layers)
    const int NC = p->imax*p->jmax*K;

    // one cell in y (j_dir 0): only the row j = 0 is written and read (x faces), the others stay 0 from the first
    // step; zeroing whole arrays was 5 % of a 1D non-hydrostatic step
    const bool row0 = p->j_dir==0 && dQxL.size()==(size_t)NC;

    for(std::vector<double> *v : {&dQxL,&dQxR,&dQyL,&dQyR,&pQxL0,&pQxR0,&pQxL1,&pQxR1,&pQyL0,&pQyR0,&pQyL1,&pQyR1})
    {
        if(!row0)
        v->assign(NC,0.0);

        else
        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        std::fill(v->begin()+dqi(p,ii,0,0),v->begin()+dqi(p,ii,0,0)+K,0.0);
    }

    auto ci = [&]() { return dqi(p,i,j,k); };

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
        const double a = (dt*g*Hf*theta*R[PIJK] + Kf)/dn;
        dL[ci()] = a;
        dR[ci()] = -a;
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
            const double f = vr ? dt*c/sg*(iR!=iL ? vrSx : vrSy)[PIJK] : dt*c/sg;

            pL0[ci()] = act(iL,jL) ? f*(sg*hL - sl) : 0.0;
            pR0[ci()] = act(iR,jR) ? -f*(sg*hR + sl) : 0.0;
            pL1[ci()] = (k<K-1 && act(iL,jL)) ? f*(sg*hL + su) : 0.0;
            pR1[ci()] = (k<K-1 && act(iR,jR)) ? -f*(sg*hR - su) : 0.0;
            }
        }

        // 6DOF: a layer blocked by a body carries the discharge of the body motion (face_update), independent of the
        // unknowns
        if(fbon)
        for(k=0;k<K;++k)
        if(fb_blocked(p,iR!=iL ? 1 : 2))
        {
        dL[ci()] = dR[ci()] = 0.0;
        pL0[ci()] = pR0[ci()] = pL1[ci()] = pR1[ci()] = 0.0;
        }

        i=il; j=jl;
    };

    SLICELOOP1
    if(Hx(i,j)>0.0)
    face(i,j,i+1,j,fdn(p,1,i,j),nhcv ? fdn(p,1,i,j) : p->DXP[IP],Hx(i,j),kx(i,j),Rx,dQxL.data(),dQxR.data(),pQxL0.data(),pQxR0.data(),pQxL1.data(),pQxR1.data());

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0)
    face(i,j,i,j+1,fdn(p,2,i,j),nhcv ? fdn(p,2,i,j) : p->DYP[JP],Hy(i,j),ky(i,j),Ry,dQyL.data(),dQyR.data(),pQyL0.data(),pQyR0.data(),pQyL1.data(),pQyR1.data());

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
            a = -dt*g*Hf*theta*f.s/f.dxp*R[PIJK];

            if(f.type==bc_riemann || f.type==bc_flather)
            {
            const bool inflow = f.A - f.B*eta0(f.i,f.j)>0.0;
            a = -f.s*f.B*(inflow ? f.Qg[k] : f.wo[k]);
            }

        dL[ci()] = intR ? 0.0 : a;
        dR[ci()] = intR ? a : 0.0;

            for(int q=0;q<4;++q)
            P4[q][ci()] = 0.0;
        }
    }

    // MPI: the faces of the neighbouring ranks through a 3D array (exchangeF)
    if(p->mpi_size>1)
    {
    const int N = p->imax*p->jmax*(p->kmax+2);
    std::vector<double> scr(N,0.0);

        for(std::vector<double> *v : {&dQxL,&dQxR,&dQyL,&dQyR,&pQxL0,&pQxR0,&pQxL1,&pQxR1,&pQyL0,&pQyR0,&pQyL1,&pQyR1})
        {
        for(i=p->imin; i<p->imin+p->imax; ++i)
        for(j=p->jmin; j<p->jmin+p->jmax; ++j)
        for(k=0; k<K; ++k)
        scr[IJK] = (*v)[ci()];

        exchangeF(p,pgc,scr.data());

        for(i=p->imin; i<p->imin+p->imax; ++i)
        for(j=p->jmin; j<p->jmin+p->jmax; ++j)
        for(k=0; k<K; ++k)
        (*v)[ci()] = scr[IJK];
        }
    }
}

void nhflow_momentum_SI::nh_assemble(lexer *p, const std::function<void(int,int,int,int,int,double)> &sink)
{
    // entry by entry: the nonzero entries of the blocks of every column (nh_assemble_col)
    const int B = p->knoz+1;
    const int nd = (p->j_dir==1) ? 5 : 3;

    nh_assemble_col(p,[&](int ii, int jj, const double *A)
    {
        for(int dd=0;dd<nd;++dd)
        for(int tr=0;tr<B;++tr)
        for(int t=0;t<B;++t)
        {
        const double v = A[((size_t)dd*B + tr)*B + t];
        if(v!=0.0)
        sink(ii,jj,dd,tr,t,v);
        }
    });
}

void nhflow_momentum_SI::nh_assemble_col(lexer *p, const std::function<void(int,int,const double*)> &sink)
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
        const size_t q = dqi(p,fi,fj,kk);
        at(bL,tr,0) += w*dL[q];
        at(bR,tr,0) += w*dR[q];
        at(bL,tr,1+kk) += w*pL0[q];
        at(bR,tr,1+kk) += w*pR0[q];
            if(kk<K-1)
            {
            at(bL,tr,2+kk) += w*pL1[q];
            at(bR,tr,2+kk) += w*pR1[q];
            }
    };

    // a row with the divergence of the layer discharges k: east and west, north and south faces
    auto adddiv = [&](int tr, double w, int kk)
    {
        // (curvilinear grid: the face lengths over the cell area, as div)
        const double ax = nhcv ? w*flen(p,1,i,j)/carea(p,i,j) : w/p->DXN[IP];
        const double axw = nhcv ? w*flen(p,1,i-1,j)/carea(p,i,j) : ax;
        addf(tr, ax,i,  j,0,1,kk,dQxL,dQxR,pQxL0,pQxR0,pQxL1,pQxR1);
        addf(tr,-axw,i-1,j,2,0,kk,dQxL,dQxR,pQxL0,pQxR0,pQxL1,pQxR1);

            if(yd)
            {
            const double ay = nhcv ? w*flen(p,2,i,j)/carea(p,i,j) : w/p->DYN[JP];
            const double ays = nhcv ? w*flen(p,2,i,j-1)/carea(p,i,j) : ay;
            addf(tr, ay,i,j,  0,3,kk,dQyL,dQyR,pQyL0,pQyR0,pQyL1,pQyR1);
            addf(tr,-ays,i,j-1,4,0,kk,dQyL,dQyR,pQyL0,pQyR0,pQyL1,pQyR1);
            }
    };

    // mass row: eta + d - WL^n + dt sum_k dsig_k div(theta Q_k)
    at(0,0,0) += nPw[IJ];
    if(vr)
    at(0,0,0) += ((*psi->stor)(i,j) - 1.0)*nPw[IJ];    // VRANS: the storage n_s (d->PORS)
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

                // -[u . grad z]_k: the change of eta of the cells of the interface slopes (implicit geometry of t^{n+1})
                {
                const double *UU = (tl2 || nhcv) ? nUZx.data() : UHD, *VV = (tl2 || nhcv) ? nUZy.data() : VHD;
                const double iH = 1.0/H1;
                const size_t qc = (size_t)((i-p->imin)*p->jmax + (j-p->jmin));
                k=kk;
                const double uk = UU[PIJK]*iH, vk = yd ? VV[PIJK]*iH : 0.0;
                const double um = kk>0 ? UU[PIJKm1]*iH : uk, vm = kk>0 ? (yd ? VV[PIJKm1]*iH : 0.0) : vk;
                const double su = p->ZN[kk+1+marge], sl = p->ZN[kk+marge];
                auto blk = [&](int dir, int off) { return off==0 ? 0 : (dir==1 ? (off>0 ? 1 : 2) : (off>0 ? 3 : 4)); };

                    if(nSdx[qc]>0.0)
                    {
                    k=kk+1; const double fu = NSFx[PIJK];
                    k=kk;   const double fl = NSFx[PIJK];
                    const double cx = (uk*su*fu - um*sl*fl)/nSdx[qc];
                    if(nSx[2*qc]>-9)   at(blk(1,nSx[2*qc]),1+kk,0)   -= sc*cx;
                    if(nSx[2*qc+1]>-9) at(blk(1,nSx[2*qc+1]),1+kk,0) += sc*cx;
                    }

                    if(yd && nSdy[qc]>0.0)
                    {
                    k=kk+1; const double fu = NSFy[PIJK];
                    k=kk;   const double fl = NSFy[PIJK];
                    const double cy = (vk*su*fu - vm*sl*fl)/nSdy[qc];
                    if(nSy[2*qc]>-9)   at(blk(2,nSy[2*qc]),1+kk,0)   -= sc*cy;
                    if(nSy[2*qc+1]>-9) at(blk(2,nSy[2*qc+1]),1+kk,0) += sc*cy;
                    }
                k=kk;
                }

            // 6DOF: w of a cell of a body is the body velocity, independent of phi (fb_wcell)
            k=kk;
            const int fbc = fbon ? fb_wcell(p) : 0;

            std::fill(wc.begin(),wc.end(),0.0);

            if(fbc==0)
            {
            const double dtw = vr ? dt*vr_sw(p) : dt;    // VRANS (k = kk)
            wc[kk] += dtw/h;
            if(kk<K-1)
            wc[kk+1] -= dtw/h;
            }

                for(int m=0;m<K;++m)
                {
                at(0,1+kk,1+m) += sc*(2.0*wc[m] - 2.0*wl[m]);
                wl[m] = 2.0*wc[m] - wl[m];
                }

                // a cell of a body above another one: phi_k - phi_{k-1}
                if(fbc==2)
                {
                for(int dd=0;dd<5;++dd)
                for(int t=0;t<B;++t)
                at(dd,1+kk,t) = 0.0;

                at(0,1+kk,1+kk) = 1.0;
                at(0,1+kk,kk) = -1.0;
                }
            }
        }

    sink(i,j,A.data());
    }
}

void nhflow_momentum_SI::nh_solve(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double starttime_nh = pgc->timer();
    const int K = p->knoz, B = K+1;

    nh_prepare(p,d,pgc);
    nh_fixed = true;

    // fourth order: the explicit part of the phi face terms of this step
    if(ho)
    ho_phi_corr(p,d,pgc);

    // curvilinear grid with A 539 1: the non-orthogonal part of the phi face terms
    if(nhco)
    nh_curv_phi_corr(p,d,pgc);

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

    // wetting and drying (A 537 1): the mass row of a column is pw (eta + d) - WL^{n+1}(fluxes), pw = 1 where the column
    // is wet at the iterate (eta + d > 0) or has no open face, 0 otherwise (a dry column: its eta from the balance of
    // its faces, its water level 0), as the hydrostatic eta solve (si_eta2D::solve_wetdry, Casulli 2009); the coupled
    // system solved again with the new pattern until it no longer changes (at most 20 times). A 537 0: pw = 1, once.
    // With the second-order time levels (tl2) the first pass starts from the wet pattern of the hydrostatic estimate of
    // this step (its eta solve iterates the pattern as well, si_eta2D::solve_wetdry): from all wet, a step with a dry
    // column next to a wet one took a second pass (conical island: 1.6 passes per step, 0.04 from the estimate; 4 layers
    // 531 s -> 406 s, A 537 0: 419 s, the same solution). A further pass then starts FGMRES from the solution of the pass
    // before (the preconditioner rebuilt for the changed mass rows, nh_solve_mg); where the pattern still changes after
    // 4 passes, the passes start again from all wet, the start of Casulli's iteration. Without tl2: from all wet, every pass from the
    // guess of the step, the preconditioner rebuilt for a new pattern
    const int itwd = p->A537==1 ? 20 : 1;
    const bool est = itwd>1 && tl2;
    nPw.assign((size_t)p->imax*p->jmax,1.0);

    if(est)
    SLICELOOP4
    nPw[IJ] = psi->wd_pattern(i,j)>0.5 ? 1.0 : 0.0;

    nh_wdpass = false;
    nh_warm = false;
    bool restarted = !est;
    int wdit = 0;

    for(int itw=0; itw<itwd; ++itw)
    {
    ++wdit;

    // F(0): the affine part; b = -F(0)
    std::fill(x.begin(),x.end(),0.0);
    nh_residual(p,d,pgc,x,F0,true);

    // ---- multigrid with the column relaxation (nhflow_si_nhmg) inside FGMRES; one row of columns
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

        if(itwd==1)
        break;

    // the new wet pattern from eta of x
    int change = 0;
    {
    int n=0;
    SLICELOOP4
    {
    const bool coupled = Hx(i,j)>0.0 || Hx(i-1,j)>0.0 || (p->j_dir==1 && (Hy(i,j)>0.0 || Hy(i,j-1)>0.0));
    const double pn = (x[(size_t)B*n] + d->depth(i,j)>0.0 || !coupled) ? 1.0 : 0.0;

        if(pn!=nPw[IJ])
        {
        nPw[IJ] = pn;
        ++change;
        }
    ++n;
    }
    }

    change = pgc->globalisum(change);

        if(change==0)
        break;

        // the pattern of the step before did not settle: from all wet
        if(itw==3 && !restarted)
        {
        std::fill(nPw.begin(),nPw.end(),1.0);
        restarted = true;
        }

    nh_wdpass = true;
    nh_warm = est;
    }

    nh_wdpass = false;
    nh_warm = false;

    // the iterations of the wet pattern and the dry columns for the log (as the hydrostatic solve)
    nh_wdit = wdit;
    nh_dry = 0;
    SLICELOOP4
    if(nPw[IJ]<0.5)
    ++nh_dry;
    nh_dry = pgc->globalisum(nh_dry);

    // the solution: etaI, PHI, the face terms and discharges (face_update: si_end_faces skips its own call where the
    // open boundaries do not depend on nh_fixed, no riemann/flather faces); the pressure for the output
    nh_residual(p,d,pgc,x,y,false,false);

    faces_set = true;
    for(const bc_face &f : bcf)
    if(f.type==bc_riemann || f.type==bc_flather)
    faces_set = false;

    // the same on all ranks: face_update exchanges the face discharges (with an absorbing edge B 99 3-5 on one rank
    // only, the ranks went different ways and waited for each other)
    faces_set = pgc->globalimin(faces_set ? 1 : 0)==1;

    FLOOP
    d->P[FIJK] = 0.0;

    SLICELOOP4
    for(k=0;k<K;++k)
    d->P[FIJK] = p->W1*PHI[PIJK];

    nh_fixed = false;
    nh_time = pgc->timer() - starttime_nh;
}

void nhflow_momentum_SI::nh_cell_faces(lexer *p, fdm_nhf *d, double *cx, double *cy)
{
    // the phi face terms of the cell (i,j), for UH and VH: the mean of those of its open faces, as the gravity term of
    // cell_update; fourth order: (-1, 9, 9, -1)/16 of four faces
    const int K = p->knoz;
    auto wn = [&](int a, int b) { return etaI(a,b) + d->depth(a,b) > p->A544; };

    // curvilinear grid: the Cartesian vector from the face terms along the normals
    if(nhcv)
    {
    nh_cell_faces_curv(p,d,cx,cy);
    return;
    }

    const bool h4x = ho && ho_cell(p,d,1,i,j);
    const bool h4y = ho && ho_cell(p,d,2,i,j);
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
    cx[k] = h4x ? (9.0*(NQx[PIm1JK] + NQx[PIJK]) - NQx[PIm2JK] - NQx[PIp1JK])/16.0 : (ww*NQx[PIm1JK] + we*NQx[PIJK])/wsx;

    cy[k] = 0.0;

    if(h4y)
    cy[k] = (9.0*(NQy[PIJm1K] + NQy[PIJK]) - NQy[PIJm2K] - NQy[PIJp1K])/16.0;

    else
    if(p->j_dir==1)
    cy[k] = (vw*NQy[PIJm1K] + ve*NQy[PIJK])/wsy;
    }
}

void nhflow_momentum_SI::nh_cell(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // WH: -dt [phi]_k/dsig_k (UH, VH: the face terms of phi are in the column solve of cell_update, nh_cell_faces)
    const double dt = p->dt;
    const int K = p->knoz;
    auto wn = [&](int a, int b) { return etaI(a,b) + d->depth(a,b) > p->A544; };

    SLICELOOP4
    if(p->wet[IJ]==1 && wn(i,j))
    {
        for(k=0;k<K;++k)
        {
            if(nAct[IJ])
            {
            const double ph_up = (k<K-1) ? PHI[PIJKp1] : 0.0;
            d->WH[IJK] -= (vr ? dt*vr_sw(p) : dt)*(ph_up - PHI[PIJK])/p->DZN[KP];
            }

            // the breaking parameter of Popinet 2020: |w| <= b sqrt(g H)
            if(nh_brk>0.0)
            {
            const double H = MAX(etaI(i,j) + d->depth(i,j),0.0);
            const double wm = nh_brk*sqrt(g*H)*H;
            d->WH[IJK] = MAX(MIN(d->WH[IJK],wm),-wm);
            }
        }
    }
}

std::vector<char> nhflow_momentum_SI::nh_mask(int B, bool yd)
{
    // the entries of the blocks of the coupled system that nh_assemble can set (the multigrid stores only these): the
    // column in full; a neighbour: the mass row every unknown (the discharges of all layers of the face), the constraint
    // row of layer k eta, phi_k and phi_k+1 (addf); no y neighbours with one cell in y
    std::vector<char> mk((size_t)5*B*B,0);

    for(int dd=0; dd<5; ++dd)
    {
        if(!yd && dd>2)
        continue;

        for(int r=0; r<B; ++r)
        for(int c=0; c<B; ++c)
        if(dd==0 || r==0 || c==0 || c==r || c==r+1)
        mk[(size_t)dd*B*B + r*B + c] = 1;
    }

    return mk;
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
    const std::vector<char> mk = nh_mask(B,yd);
    nhmg->setup(pgc->mpi_comm,B,nx,ny,nbr,&mk);

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

    // rows in metres (below): the constraint rows of an active column times rs, here with the assembly of the column
    const double rs = 0.5*p->dt*p->dt;

    if(!probe)
    {
    nh_dfaces(p,d,pgc);
    nh_assemble_col(p,[&](int ii, int jj, const double *A)
    {
        const int qq = (ii-p->imin)*p->jmax + (jj-p->jmin);
        M.set_column(ii,jj,nc==5 ? 5 : 3,A,1,nAct[qq] ? rs : 1.0);
    });
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
        const double a = M.get(ii,jj,dd,q), b = Dg[(((size_t)ii*ny + jj)*5 + dd)*B*B + q];
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
    if(probe)
    for(int n=0; n<N; ++n)
    {
    i = cix[n]; j = ciy[n];

        if(nAct[IJ])
        for(int dd=0; dd<5; ++dd)
        {
        if(!yd && dd>2)
        continue;

        M.scale_rows(i,j,dd,1,rs);
        }
    }

    // ---- the right-hand side and the initial guess on the box vectors
    std::vector<double> &b = nmb, &x = nmx, &w = nmw;
    b.assign(nm,0.0);
    w.resize(nm);
    const bool warm = nh_warm && x.size()==nm;
    if(!warm)
    x.assign(nm,0.0);
    const bool guess = true;
    const bool hist = nEtaP.size()==(size_t)p->imax*p->jmax;

    for(int n=0; n<N; ++n)
    {
    const size_t q = ((size_t)cix[n]*ny + ciy[n])*B;
    i = cix[n]; j = ciy[n];
        for(int t=0; t<B; ++t)
        b[q+t] = -F0[(size_t)B*n+t]*(t>0 && nAct[IJ] ? rs : 1.0);

    // a further pass of the wet pattern: the solution of the pass before (x as it is)
    if(warm)
    continue;

    // initial guess: eta^n and phi of the step before, both extrapolated in time (r = dt/dt_old)
    x[q] = eta0(i,j);
        for(k=0;k<K;++k)
        x[q+1+k] = nAct[IJ] ? PHI[PIJK] : 0.0;

        if(guess && hist)
        {
        const double r = p->dt/MAX(p->dt_old,1.0e-20);
        x[q] += r*(eta0(i,j) - nEtaP[IJ]);
            for(k=0;k<K;++k)
            if(nAct[IJ])
            x[q+1+k] += r*(PHI[PIJK] - nPhiP[IJK]);
        }
    }

    // the levels of this step for the guess of the next one (once per step: a further pass of the wet pattern found
    // them set, and with them the extrapolation of its own guess gone)
    if(!warm)
    {
    const int NC = p->imax*p->jmax*(p->kmax+2);
    nEtaP.resize((size_t)p->imax*p->jmax);
    nPhiP.resize(NC);
    SLICELOOP4
    {
    nEtaP[IJ] = eta0(i,j);
        for(k=0;k<K;++k)
        nPhiP[IJK] = PHI[PIJK];
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
    // reb, itreb: the preconditioner rebuilt inside the solve (stall, below) and the iterations before that
    struct fg_result { int it=0, cyc=0, reb=0, itreb=0; double r0max=0.0, rmax=0.0, tolmax=0.0; };

    // tabs > 0: at most this stopping residual (a further pass of the wet pattern: that of the first pass, or A 506 of its
    // own starting residual where that is smaller)
    auto fgmres = [&](nhflow_si_nhmg &MM, const std::vector<double> &bb, std::vector<double> &xx,
                      std::vector<std::vector<double>> &V, std::vector<std::vector<double>> &Z, std::vector<double> &wv,
                      double tabs=0.0) -> fg_result
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
        std::vector<std::pair<int,double>> hst;     // iterations and max |r| at every restart (stall test)
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
            R.tolmax = MAX(tabs>0.0 ? MIN(tabs,tol*R.r0max) : tol*R.r0max,1.0e-12*MAX(mx[1],mx[2]));
            }

            if(R.rmax<=R.tolmax || R.it>=maxit || beta<=0.0)
            break;

            // stall: the max-norm residual not reduced tenfold over the last 2 m iterations (two restart cycles; a
            // solve takes about 5): the preconditioner of an earlier step no longer fits this operator, rebuilt from
            // it, once per solve. Measured against the residual 2 m iterations back, not r0: the stalls of the stale
            // preconditioner reduce r0 1e4 times in a few iterations and then stay at 2-7 mm (harbour smoke test)
            hst.push_back({R.it,R.rmax});

            if(!R.reb && R.it>=2*m)
            {
            double rb = -1.0;
            for(const auto &h : hst)
            if(h.first<=R.it-2*m)
            rb = h.second;

                if(rb>0.0 && R.rmax>0.1*rb)
                {
                MM.finalize(true);
                R.reb = 1;
                R.itreb = R.it;
                }
            }

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
    // 10 %, when the active columns changed, or when the iterations of the first pass grew by more than one over those of
    // the rebuild. The iterations of a further pass of the wet pattern from the solution of the pass before (nh_warm)
    // set no rebuild; a further pass from the guess of the step rebuilds it.
    // Also rebuilt when a column changed its wet pattern (pw of the mass row, A 537 1) or became active or inactive since
    // the rebuild, a further pass of the pattern included: a dry mass row (pw = 0) has only the face terms, about
    // dt^2 g A 544/dx^2 of a wet one (2e-5 at dt 0.31 s, dx 10 m, A 544 0.01), so a preconditioner of the other pattern
    // inverts that row 1e5 times too large or too small (harbour 10 m, 7 columns re-wetted under a preconditioner two
    // steps old: 300 iterations instead of 4)
    const size_t ncol = (size_t)p->imax*p->jmax;
    auto colpat = [&](int q) { return (char)((nPw[q]>0.5 ? 1 : 0) + (nAct[q] ? 2 : 0)); };
    int pchg = nh_lagpat.size()!=ncol ? 1 : 0;

    if(pchg==0)
    SLICELOOP4
    if(nh_lagpat[IJ]!=colpat(IJ))
    ++pchg;

    pchg = pgc->globalisum(pchg);

    const bool refresh = lagmax<=1 || nh_lagc+1>=lagmax || fabs(p->dt/MAX(nh_lagdt,1.0e-20) - 1.0)>0.1
                         || pgc->globalisum(nact)!=nh_lagnact || nh_iter1>nh_lagit+1 || (nh_wdpass && !nh_warm)
                         || pchg>0;

    const int nsplit = p->A507;
    fg_result R, R2;

    if(nsplit==0)
    {
    M.finalize(refresh);
    R = fgmres(M,b,x,fgV,fgZ,nmv,warm ? nh_tolabs : 0.0);

        if(!warm)
        nh_tolabs = R.tolmax;
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
    // the phi system: the constraint rows and phi of the mask of the coupled one
    const std::vector<char> mk = nh_mask(B,yd);
    std::vector<char> mkp((size_t)5*K*K,0);
    for(int dd=0; dd<5; ++dd)
    for(int r=1; r<B; ++r)
    for(int c=1; c<B; ++c)
    mkp[(size_t)dd*K*K + (r-1)*K + c-1] = mk[(size_t)dd*B*B + r*B + c];
    nhmp->setup(pgc->mpi_comm,K,nx,ny,nbr,&mkp);
    }

    nhflow_si_nhmg &Me = *nhme, &Mp = *nhmp;
    Me.zero();
    Mp.zero();

    for(int ii=0; ii<nx; ++ii)
    for(int jj=0; jj<ny; ++jj)
    for(int dd=0; dd<5; ++dd)
    {
        if(M.get(ii,jj,dd,0)!=0.0)
        Me.set(ii,jj,dd,0,0,M.get(ii,jj,dd,0));

        for(int r=1; r<B; ++r)
        for(int c=1; c<B; ++c)
        if(M.get(ii,jj,dd,r*B+c)!=0.0)
        Mp.set(ii,jj,dd,r-1,c-1,M.get(ii,jj,dd,r*B+c));
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

    // the preconditioner now in use: rebuilt at the start (refresh) or inside FGMRES (stall); the iterations of the
    // rebuild those after it
    nh_reb = R.reb + R2.reb;
    const bool rebuilt = refresh || nh_reb>0;

    if(rebuilt)
    {
    nh_lagpat.resize(ncol);
    for(size_t q=0; q<ncol; ++q)
    nh_lagpat[q] = colpat((int)q);
    }

    if(rebuilt && !warm)
    {
    nh_lagc = 0;
    nh_lagdt = p->dt;
    nh_lagnact = pgc->globalisum(nact);
    nh_lagit = R.it - R.itreb + R2.it - R2.itreb;
    }

    else if(!warm)
    ++nh_lagc;

    nh_iter = R.it + R2.it;

    if(!warm)
    nh_iter1 = nh_iter;

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
            <<" coupled r "<<rc<<(refresh ? " refresh" : "")<<(pchg>0 ? " pattern "+std::to_string(pchg) : "")
            <<(nh_reb>0 ? " stall-rebuild at "+std::to_string(R.itreb+R2.itreb) : "")<<endl;
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

// ---------------------------------------------------------------------------------------------------------------------
// fourth-order horizontal operators (nhflow_momentum_SI_ho.cpp): the non-hydrostatic pressure
// ---------------------------------------------------------------------------------------------------------------------

// the explicit part of the phi face terms (NQcx, NQcy): the face term of nh_faces with -delta^2 phi/12 along the face
// normal in place of phi, phi extrapolated to this step from the two steps before (PHI, PHIo)
void nhflow_momentum_SI::ho_phi_corr(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double dt = p->dt;
    const int K = p->knoz;
    const int N = lN;

    if(NQcx==nullptr)
    {
    p->Darray_lazy(NQcx,N);
    p->Darray_lazy(NQcy,N);
    }

    hoPS.resize(N);

    const double r = ho_phi ? p->dt/MAX(p->dt_old,1.0e-20) : 0.0;

    for(int q=0; q<N; ++q)
    {
    hoPS[q] = PHI[q] + r*(PHI[q] - PHIo[q]);
    PHIo[q] = PHI[q];
    NQcx[q] = 0.0;
    NQcy[q] = 0.0;
    }

    ho_phi = true;

    auto hh = [&](int a, int c) { return nWLh[(a-p->imin)*p->jmax + (c-p->jmin)]; };
    auto zif = [&](int a, int c, int m) { return d->bed(a,c) + p->ZN[m+marge]*hh(a,c); };
    auto ps = [&](int a, int c, int kk)
    {
        return kk<K ? hoPS[PIX(a,c,kk)] : 0.0;
    };

    // face of the cells (iL,jL) - (iR,jR) along (di,dj)
    auto face = [&](int iL, int jL, int iR, int jR, int di, int dj, double dx, double Hf, double *NQ)
    {
        const double hL = hh(iL,jL), hR = hh(iR,jR);
        const double c = (hL+hR)>1.0e-20 ? Hf/(dx*(hL + hR)) : 0.0;
        auto cc = [&](int a, int b, int kk) { return -(ps(a+di,b+dj,kk) - 2.0*ps(a,b,kk) + ps(a-di,b-dj,kk))/12.0; };

        for(k=0;k<K;++k)
        {
        const double s  = p->DZN[KP];
        const double sl = nh_lim(zif(iR,jR,k) - zif(iL,jL,k),dx);
        const double su = nh_lim(zif(iR,jR,k+1) - zif(iL,jL,k+1),dx);

        const double fL = cc(iL,jL,k), uL = cc(iL,jL,k+1);
        const double fR = cc(iR,jR,k), uR = cc(iR,jR,k+1);

        const double br = (s*hL - sl)*fL - (s*hR + sl)*fR + (s*hL + su)*uL - (s*hR - su)*uR;

        NQ[PIJK] = dt*br*c/s;
        }
    };

    SLICELOOP1
    if(Hx(i,j)>0.0 && ho_face(p,1,i,j))
    face(i,j,i+1,j,1,0,p->DXP[IP],Hx(i,j),NQcx);

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0 && ho_face(p,2,i,j))
    face(i,j,i,j+1,0,1,p->DYP[JP],Hy(i,j),NQcy);

    // open boundary faces: no pressure term
    for(const bc_face &f : bcf)
    {
    i=f.fi;
    j=f.fj;
    double *NQ = (f.cs==1 || f.cs==4) ? NQcx : NQcy;
        for(k=0;k<K;++k)
        NQ[PIJK] = 0.0;
    }

    pexchF(p,pgc,NQcx);
    pexchF(p,pgc,NQcy);
}

// the phi face terms of the step with their explicit part (nh_residual, after nh_faces and its exchange)
void nhflow_momentum_SI::ho_phi_add(lexer *p)
{
    const int N = lN;

    for(int q=0; q<N; ++q)
    {
    NQx[q] += NQcx[q];
    NQy[q] += NQcy[q];
    }
}
