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
// SFLOW semi-implicit step (A 226 1) with the quadratic non-hydrostatic pressure (A 220 2/3).
//
// The pressure of sflow_pjm_quad (Jeschke et al. 2017; A 220 3 with the improved dispersion, alpha A 224) is not
// projected after the step, but solved together with eta^{n+1} in one linear system, so that neither sqrt(g h)
// nor the dispersive wave speed limits the time step (as the theta = 1/2 coupled solve of Basilisk's layered
// solver, Popinet 2020). Per cell the unknowns eta^{n+1} and P^{n+1} = p/rho, P fully implicit; per face (x and y):
//
//  face  q_f = q^E_f - a_f/(theta dt) (eta_R - eta_L)                                      (hydrostatic SI)
//              - dt alpha_f [ (H_R P_R - H_L P_L) - qb_f (d_R - d_L) ] / dxp_f              (pressure, bed pressure)
//              + dt alpha_f g B [ theta (L_R - L_L)^{n+1} + (1-theta) (L_R - L_L)^n ] / dxp_f  (A 220 3)
//        qb_f = cb (P_L + P_R)/2 + alpha (H phi)_f / 4,  L = h^3 lap(eta),  B = (alpha - 1)/3
//  mass  eta_i + theta dt div(q) = WL^n - d - (1-theta) dt div(q^n)                    (the row of si_eta2D)
//  incompressibility (new velocities)
//        div(q) + U* (2 d_x - H_x) + V* (2 d_y - H_y) + 2 W = 0,  W = (WH* + dt (cw P + H phi/4))/H,
//        U*, V*: the cell velocities extrapolated to t^{n+1} (the slope term implicit, with the face velocities,
//        needed about 100 Krylov iterations on the Beji & Battjes bar instead of 4, the errors there the same)
//
// compact face differences throughout (no checkerboard mode, the dispersion of the staggered grid). The cell
// discharges take the mean of the non-hydrostatic face terms (as the gravity term, cell_update). phi, the bed
// acceleration Dw_b/Dt, is explicit: the acceleration of the previous step (an estimate of U^{n+1} with P^n as in
// sflow_pjm_quad is unstable for large steps with A 220 3). WH*: WH^n advected upwind with q^n (Adams-Bashforth 2).
// Time levels of the nonlinear terms (nh_levels): the depths of the gravity and A 220 3 terms at t^{n+1/2}, the depth
// of the pressure terms and of the constraint at t^{n+1}, the advecting discharge at t^{n+1/2}, extrapolated from
// t^n and t^{n-1}.
//
// Solution of the coupled system (nh_solve). The operator is applied matrix-free (apply: face discharges from
// per-face coefficients, then the rows), so that the same code serves both solvers:
//  - one cell in y, one rank: band LU (the band of the interleaved unknowns, read from the operator by probing with
//    five colours), exact;
//  - otherwise (2D, MPI): right-preconditioned FGMRES(30), relative tolerance 1e-6, initial guess extrapolated in
//    time. The preconditioner eliminates eta locally in the non-hydrostatic cells (without the slope and A 220 3
//    terms the constraint gives eta_i = s_i + g_i P_i) and solves the remaining 5-point equation for P (eta in the
//    hydrostatic cells) with one V-cycle of REEFMG2D (own hierarchy). A 220 3: the h^3 lap(eta) term makes the
//    reduced operator fourth order, approximated as M5 (I - rho lap): a second V-cycle with I - div(rho grad).
//    4-10 iterations per step (1D bar, 2D shoal).
//
// Wetting and drying (A 229 1): Newton iteration on the wet pattern of the mass rows (Casulli 2009, as si_eta2D).
// Breaking (A 246, A 242): sflow_eta::breaking before the predictor with (eta^n - eta^{n-1})/dt; the pressure stays
// on in breaking cells, the dissipation is the eddy viscosity (A 246 4: crest tracking and Kennedy et al. 2000 as
// NHFLOW A 550 4, sflow_breaking_kb; A 246 2: A 250) through the implicit diffusion of the predictor.
// Time step (sflow_timestep_SI): the advective limit is halved (WH with AB2), it acts in run-down and bores.
//
// Limits (check()): one grid, Cartesian. Non-hydrostatic cells: wet,
// deeper than 10 A 244, not next to the shoreline (p->deep), wet neighbours; P = 0 elsewhere. Open boundary faces carry no pressure term.
// ---------------------------------------------------------------------------

#include"sflow_momentum_SI.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"si_eta2D.h"
#include"sflow_fsf.h"
#include"solver2D.h"
#include"reefmg2D.h"
#include"matrix2D.h"
#include"vec2D.h"
#include<vector>
#include<memory>
#include<algorithm>

// work space of the coupled solve: band LU, FGMRES vectors, the face map of the open boundaries
struct sflow_si_nh_lu
{
    // band LU (row-major, each row holds the columns r-kl .. r+ku+kl)
    int n=0, kl=5, ku=5, W=0;
    std::vector<double> AB;
    std::vector<int> piv;

    double& at(int r, int c) { return AB[(size_t)r*W + (c - r + kl)]; }

    void resize(int nn) { n=nn; W=2*kl+ku+1; AB.assign((size_t)n*W,0.0); piv.assign(n,0); }

    bool factor()
    {
        for(int k=0; k<n; ++k)
        {
            int pr=k; double pm=fabs(at(k,k));
            for(int r=k+1; r<=std::min(n-1,k+kl); ++r)
            if(fabs(at(r,k))>pm) { pm=fabs(at(r,k)); pr=r; }

            if(pm==0.0)
            return false;

            piv[k]=pr;
            const int cmax = std::min(n-1,k+ku+kl);

            if(pr!=k)
            for(int c=k; c<=cmax; ++c)
            std::swap(at(k,c),at(pr,c));

            const double d = at(k,k);

            for(int r=k+1; r<=std::min(n-1,k+kl); ++r)
            {
                const double l = at(r,k)/d;
                at(r,k) = l;

                if(l!=0.0)
                for(int c=k+1; c<=cmax; ++c)
                at(r,c) -= l*at(k,c);
            }
        }
        return true;
    }

    void solve(std::vector<double> &x)
    {
        for(int k=0; k<n; ++k)
        {
            if(piv[k]!=k)
            std::swap(x[k],x[piv[k]]);

            for(int r=k+1; r<=std::min(n-1,k+kl); ++r)
            x[r] -= at(r,k)*x[k];
        }

        for(int k=n-1; k>=0; --k)
        {
            double s = x[k];
            for(int c=k+1; c<=std::min(n-1,k+ku+kl); ++c)
            s -= at(k,c)*x[c];
            x[k] = s/at(k,k);
        }
    }

    // FGMRES
    std::vector<std::vector<double>> V, Z;
    std::vector<double> H, cs, sn, g;

    // open boundary face of x-face / y-face (i,j): index into bcf, -1 none
    std::vector<int> bfx, bfy;
    int nbcf=-1;

    // preconditioner: M5 (the 5-point Schur operator) and A2 = I - div(rho grad) (the A 220 3 term), one V-cycle each
    std::unique_ptr<reefmg2D> mg1, mg2;
    std::unique_ptr<matrix2D> M1, M2;
    std::unique_ptr<vec2D> r1, r2;
};

void sflow_momentum_SI::nh_ini(lexer *p)
{
    nh = (p->A220==2 || p->A220==3);

    cb = 1.5;
    adisp = (p->A220==3) ? MAX(p->A224,1.0) : 1.0;
    cwn = cb/adisp;
    Bdisp = (adisp-1.0)/3.0;

    SLICELOOP4
    {
    phiN(i,j) = 0.0;
    WHp(i,j) = 0.0;
    Lc0(i,j) = 0.0;
    Unm1(i,j) = 0.0;
    Vnm1(i,j) = 0.0;
    WLo(i,j) = -1.0;
    }

    SLICELOOP1
    {
    qxN(i,j) = 0.0;
    Hxo(i,j) = 0.0;
    qxo(i,j) = 0.0;
    }

    SLICELOOP2
    {
    qyN(i,j) = 0.0;
    Hyo(i,j) = 0.0;
    qyo(i,j) = 0.0;
    }

    // the Krylov solver also on a 1D grid with one rank (development: nh_krylov = true compares it with the band LU)
    nh_krylov = false;
}

void sflow_momentum_SI::nh_levels(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // time levels of the nonlinear terms, extrapolated from t^n and t^{n-1}: with the depth of t^n the waves on the
    // Beji & Battjes bar lagged in proportion to dt (0.05 s at dt 0.02 s after the bar crest) and the released
    // harmonics grew by 5 %:
    //  - face and cell depth of the gravity and A 220 3 terms at t^{n+1/2} (theta 1/2)
    //  - depth of the pressure terms and of the constraint at t^{n+1} (P^{n+1} and the velocities of t^{n+1})
    //  - discharge of the advection (predictor) at t^{n+1/2}
    // the continuity equation keeps WL^n and q^n: the volume stays exact
    const double r = p->dt/MAX(p->dt_old,1.0e-20);
    const double e1 = 1.0 + 0.5*r, e0 = 0.5*r;

    SLICELOOP1
    {
    const double h = Hx(i,j);

        if(h>0.0 && Hxo(i,j)>0.0)
        Hx(i,j) = MAX(e1*h - e0*Hxo(i,j), p->A244);

    Hxo(i,j) = h;

    qxa(i,j) = nh_first ? qxn(i,j) : e1*qxn(i,j) - e0*qxo(i,j);
    qxo(i,j) = qxn(i,j);
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    const double h = Hy(i,j);

        if(h>0.0 && Hyo(i,j)>0.0)
        Hy(i,j) = MAX(e1*h - e0*Hyo(i,j), p->A244);

    Hyo(i,j) = h;

    qya(i,j) = nh_first ? qyn(i,j) : e1*qyn(i,j) - e0*qyo(i,j);
    qyo(i,j) = qyn(i,j);
    }

    SLICELOOP4
    {
    const bool ok = WLo(i,j)>0.0 && WL0(i,j)>p->A244;
    WLh(i,j) = ok ? MAX(e1*WL0(i,j) - e0*WLo(i,j), p->A244) : WL0(i,j);
    WLq(i,j) = ok ? MAX((1.0+r)*WL0(i,j) - r*WLo(i,j), p->A244) : WL0(i,j);
    WLo(i,j) = WL0(i,j);
    }

    // wave breaking (A 246, A 242; sflow_eta::breaking): the rise of eta over the last step, (eta^n - eta^{n-1})/dt;
    // the non-hydrostatic pressure stays on in breaking cells (nh_active), the dissipation is the eddy viscosity of
    // A 246 4 (crest tracking and Kennedy et al. 2000, as NHFLOW A 550 4) or A 246 2 (A 250); the A 246 1/2 flags are
    // reset after the step (sflow_f, A 248 0) as for the explicit stages
    if(nh_first)
    SLICELOOP4
    nEb(i,j) = b->eta(i,j);

    pfsf->breaking(p,b,pgc,b->eta,nEb,p->dt_old/MAX(p->dt,1.0e-20));

    SLICELOOP4
    nEb(i,j) = b->eta(i,j);

    exchange1(p,pgc,Hx);
    exchange1(p,pgc,qxa);
    exchange2(p,pgc,Hy);
    exchange2(p,pgc,qya);
    exchange4(p,pgc,WL0);
    exchange4(p,pgc,WLh);
    exchange4(p,pgc,WLq);
}

int sflow_momentum_SI::nh_active(lexer *p, fdm2D *b, int ii, int jj)
{
    // wet, deep, the neighbours in the domain wet (as sflow_pjm_quad::active without its breaking switch; a wall is no dry cell)
    auto id = [&](int a, int c) { return (a-p->imin)*p->jmax + (c-p->jmin); };
    auto dry = [&](int a, int c) { return p->flagslice4[id(a,c)]>0 && p->wet[id(a,c)]==0; };

    // breaking cells keep the pressure: switched off, the edges of the breaking regions radiated spurious waves and
    // drove the waves up on the Mase & Kirby slope; the breaking dissipation is the eddy viscosity (A 246 4 or 2)
    if(p->wet[id(ii,jj)]==0 || WL0(ii,jj)<=10.0*p->A244)
    return 0;

    // shallow cells next to the shoreline (sflow_eta::wetdrydeep: dry cells within two cells): hydrostatic
    if(p->deep[id(ii,jj)]==0)
    return 0;

    if(dry(ii-1,jj) || dry(ii+1,jj))
    return 0;

    if(p->j_dir==1 && (dry(ii,jj-1) || dry(ii,jj+1)))
    return 0;

    return 1;
}

void sflow_momentum_SI::nh_predictor(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // h^3 lap(eta^n) (A 220 3), WH predictor, bed acceleration phi
    const double dt = p->dt;
    const bool yd = (p->j_dir==1);
    auto id = [&](int a, int c) { return (a-p->imin)*p->jmax + (c-p->jmin); };
    auto in = [&](int a, int c) { return p->flagslice4[id(a,c)]>0; };

    // L = h^3 lap(eta), compact, mirrored at the domain edges (b->eta: eta^n with the halo of the last step)
    SLICELOOP4
    {
    Lc0(i,j) = 0.0;

        if(Bdisp>0.0 && nh_active(p,b,i,j)==1)
        {
        const double h = WLh(i,j);
        double lap = ((in(i+1,j) ? (b->eta(i+1,j) - b->eta(i,j))/p->DXP[IP] : 0.0)
                    - (in(i-1,j) ? (b->eta(i,j) - b->eta(i-1,j))/p->DXP[IM1] : 0.0))/p->DXN[IP];

        if(yd)
        lap += ((in(i,j+1) ? (b->eta(i,j+1) - b->eta(i,j))/p->DYP[JP] : 0.0)
              - (in(i,j-1) ? (b->eta(i,j) - b->eta(i,j-1))/p->DYP[JM1] : 0.0))/p->DYN[JP];

        Lc0(i,j) = h*h*h*lap;
        }
    }

    exchange4(p,pgc,Lc0);

    // W^n advected upwind with q^n (conservative), Adams-Bashforth 2 in time (the forward step gave 1.5 % too large
    // harmonics on the lee side of a bar, dt 0.02 s)
    auto Wn = [&](int a, int c) { return (in(a,c) && WL0(a,c)>p->A244) ? b->WH(a,c)/WL0(a,c) : 0.0; };
    auto fWx = [&](int a, int c) { const double q = qxn(a,c); return q>=0.0 ? q*Wn(a,c) : q*Wn(a+1,c); };
    auto fWy = [&](int a, int c) { const double q = qyn(a,c); return q>=0.0 ? q*Wn(a,c) : q*Wn(a,c+1); };

    const double r = dt/MAX(p->dt_old,1.0e-20);
    const bool ab2 = !nh_first;

    SLICELOOP4
    {
    double A = 0.0;

        if(p->wet[IJ]==1)
        {
        A = (fWx(i,j) - fWx(i-1,j))/p->DXN[IP];

        if(yd)
        A += (fWy(i,j) - fWy(i,j-1))/p->DYN[JP];
        }

        if(p->wet[IJ]==1)
        WHp(i,j) = b->WH(i,j) - dt*(ab2 ? (1.0+0.5*r)*A - 0.5*r*AWo(i,j) : A);

        else
        WHp(i,j) = 0.0;

    AWo(i,j) = A;
    }

    // U^n, V^n of the cells; the acceleration of the previous step, (U^n - U^{n-1})/dt, gives phi: the estimate of
    // U^{n+1} with the explicit terms of the step (sflow_pjm_quad, P^n) is not balanced by the pressure for large
    // time steps, the h^3 lap(eta) term of A 220 3 then drives a grid mode on a sloping bed (Beji & Battjes bar,
    // dt 0.005-0.02 s); solving again with phi from U^{n+1} changed nothing there
    SLICELOOP4
    {
    const bool w = (p->wet[IJ]==1 && WL0(i,j)>p->A244);
    nU(i,j) = w ? UH0(i,j)/WL0(i,j) : 0.0;
    nV(i,j) = (w && yd) ? VH0(i,j)/WL0(i,j) : 0.0;

        if(nh_first)
        {
        Unm1(i,j) = nU(i,j);
        Vnm1(i,j) = nV(i,j);
        }
    }

    exchange4(p,pgc,nU);
    exchange4(p,pgc,nV);

    nh_first = false;

    // bed acceleration phi = Dw_b/Dt, w_b = -u.grad(d) (sflow_pjm_quad::quad_calc)
    auto Ue = [&](int a, int c) { return in(a,c) ? nU(a,c) : b->U(a,c); };
    auto Ve = [&](int a, int c) { return in(a,c) ? nV(a,c) : b->V(a,c)*p->y_dir; };
    auto dep = [&](int a, int c, int a0, int c0) { return in(a,c) ? b->depth(a,c) : b->depth(a0,c0); };

    SLICELOOP4
    {
    phiN(i,j) = 0.0;

        if(nh_active(p,b,i,j)==1)
        {
        const double dx = p->DXP[IP] + p->DXP[IM1];
        const double dE = dep(i+1,j,i,j), dW = dep(i-1,j,i,j), dC = b->depth(i,j);

        const double dddx  = (dE - dW)/dx;
        const double dddxx = ((dE - dC)/p->DXP[IP] - (dC - dW)/p->DXP[IM1])/p->DXN[IP];

        const double u = nU(i,j), v = nV(i,j);
        double dudt = (u - Unm1(i,j))/p->dt_old + u*(Ue(i+1,j) - Ue(i-1,j))/dx;

            if(yd)
            {
            const double dy = p->DYP[JP] + p->DYP[JM1];
            const double dN = dep(i,j+1,i,j), dS = dep(i,j-1,i,j);
            const double dddy  = (dN - dS)/dy;
            const double dddyy = ((dN - dC)/p->DYP[JP] - (dC - dS)/p->DYP[JM1])/p->DYN[JP];

            double dddxy = 0.0;
            if(in(i+1,j+1) && in(i+1,j-1) && in(i-1,j+1) && in(i-1,j-1))
            dddxy = (b->depth(i+1,j+1) - b->depth(i+1,j-1) - b->depth(i-1,j+1) + b->depth(i-1,j-1))/(dx*dy);

            dudt += v*(Ue(i,j+1) - Ue(i,j-1))/dy;
            const double dvdt = (v - Vnm1(i,j))/p->dt_old + u*(Ve(i+1,j) - Ve(i-1,j))/dx + v*(Ve(i,j+1) - Ve(i,j-1))/dy;

            phiN(i,j) = -dvdt*dddy - v*v*dddyy - 2.0*u*v*dddxy;
            }

        phiN(i,j) += -dudt*dddx - u*u*dddxx;
        }
    }

    exchange4(p,pgc,phiN);

    // U^{n+1}, V^{n+1} extrapolated, for the slope term of the constraint
    SLICELOOP4
    {
    nUs(i,j) = (1.0+r)*nU(i,j) - r*Unm1(i,j);
    nVs(i,j) = (1.0+r)*nV(i,j) - r*Vnm1(i,j);
    Unm1(i,j) = nU(i,j);
    Vnm1(i,j) = nV(i,j);
    }
}

void sflow_momentum_SI::nh_solve(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // the coupled eta-P system, unknowns 2n (eta) and 2n+1 (P) of cell n (SLICELOOP4 order)
    const double starttime_nh = pgc->timer();

    const double dt = p->dt;
    const double tdt = theta*dt;
    const bool yd = (p->j_dir==1);

    auto id = [&](int a, int c) { return (a-p->imin)*p->jmax + (c-p->jmin); };
    auto in = [&](int a, int c) { return p->flagslice4[id(a,c)]>0; };

    if(!nhlu)
    nhlu = std::make_shared<sflow_si_nh_lu>();

    auto &S = *nhlu;

    // open boundary faces of x and y faces
    if(S.nbcf!=(int)bcf.size() || (int)S.bfx.size()!=p->imax*p->jmax)
    {
    S.bfx.assign(p->imax*p->jmax,-1);
    S.bfy.assign(p->imax*p->jmax,-1);

        for(int k=0; k<(int)bcf.size(); ++k)
        {
        const bc_face &f = bcf[k];
        if(f.cs==1 || f.cs==4) S.bfx[id(f.fi,f.fj)] = k;
        else                   S.bfy[id(f.fi,f.fj)] = k;
        }

    S.nbcf = (int)bcf.size();
    }

    // cells, non-hydrostatic cells
    int N=0;
    nh_cells_active = 0;

    SLICELOOP4
    {
    nA(i,j) = (double)nh_active(p,b,i,j);
    nh_cells_active += (int)nA(i,j);
    ++N;
    }

    exchange4(p,pgc,nA);
    auto ac = [&](int a, int c) { return in(a,c) && nA(a,c)>0.5; };

    // ---- face coefficients: q_NH = -R P_R + L P_L + B (L_R - L_L) + c0
    auto fcoef = [&](int a, int c, int a1, int c1, double dxp, double alf, double Hf,
                     double &cR, double &cL, double &cB, double &c0)
    {
        cR = cL = cB = c0 = 0.0;

        if(!(in(a,c) && in(a1,c1)) || Hf<=0.0)
        return;

        const double dd = b->depth(a1,c1) - b->depth(a,c);
        const double HL = WLq(a,c), HR = WLq(a1,c1);
        const double s = dt*alf/dxp;

        cR = s*(HR - 0.5*cb*dd);
        cL = s*(HL + 0.5*cb*dd);
        c0 = s*0.125*adisp*(HL*phiN(a,c) + HR*phiN(a1,c1))*dd;

        // A 220 3 term only between two non-hydrostatic cells: L = 0 in a hydrostatic (breaking, shallow) cell, the
        // jump of L at the edge of a breaking region drove the waves up on the Mase & Kirby slope (Hm0 x5 after 90 s)
        if(Bdisp>0.0 && ac(a,c) && ac(a1,c1))
        {
        const double sb = s*g*Bdisp;
        cB = sb*theta;
        c0 += sb*(1.0-theta)*(Lc0(a1,c1) - Lc0(a,c));
        }
    };

    SLICELOOP1
    fcoef(i,j,i+1,j,p->DXP[IP],alx(i,j),Hx(i,j),nxR(i,j),nxL(i,j),nxB(i,j),nx0(i,j));

    if(yd)
    SLICELOOP2
    fcoef(i,j,i,j+1,p->DYP[JP],aly(i,j),Hy(i,j),nyR(i,j),nyL(i,j),nyB(i,j),ny0(i,j));

    exchange1(p,pgc,nxR); exchange1(p,pgc,nxL); exchange1(p,pgc,nxB); exchange1(p,pgc,nx0);
    exchange2(p,pgc,nyR); exchange2(p,pgc,nyL); exchange2(p,pgc,nyB); exchange2(p,pgc,ny0);

    // ---- cell coefficients: slope of the constraint, U (2 d_x - H_x) + V (2 d_y - H_y)
    SLICELOOP4
    {
    const double H = WLq(i,j);
    const double dx = p->DXP[IP] + p->DXP[IM1];
    nsx(i,j) = (2.0*((in(i+1,j) ? b->depth(i+1,j) : b->depth(i,j)) - (in(i-1,j) ? b->depth(i-1,j) : b->depth(i,j)))
                    - ((in(i+1,j) ? WLq(i+1,j) : H) - (in(i-1,j) ? WLq(i-1,j) : H)))/dx;
    nsy(i,j) = 0.0;

        if(yd)
        {
        const double dy = p->DYP[JP] + p->DYP[JM1];
        nsy(i,j) = (2.0*((in(i,j+1) ? b->depth(i,j+1) : b->depth(i,j)) - (in(i,j-1) ? b->depth(i,j-1) : b->depth(i,j)))
                        - ((in(i,j+1) ? WLq(i,j+1) : H) - (in(i,j-1) ? WLq(i,j-1) : H)))/dy;
        }
    }

    // ---- the operator: x (eta, P) -> y (mass row, constraint row); with c the affine terms
    std::vector<double> xv(2*N), yv(2*N);

    // face coefficients of the discharge of the step: q = EL E_L + ER E_R + PL P_L + PR P_R + B (L_R - L_L)
    // + k1 (mode 1: the constants of the mass row; the q^E of the faces are in rhsW) + k2 (mode 2: all constants);
    // open boundary faces in the same form (the cell inside on its side), closed faces zero
    auto kface = [&](int dir, int a, int c, slice &EL, slice &ER, slice &PL, slice &PR, slice &K1, slice &K2)
    {
        const int a1 = dir==1 ? a+1 : a, c1 = dir==1 ? c : c+1;
        const double Hf = dir==1 ? Hx(a,c) : Hy(a,c);
        const double af = dir==1 ? ax(a,c) : ay(a,c);
        const double qE = dir==1 ? qxE(a,c) : qyE(a,c);

        EL(a,c) = ER(a,c) = PL(a,c) = PR(a,c) = K1(a,c) = K2(a,c) = 0.0;

        if(in(a,c) && in(a1,c1))
        {
            if(Hf<=0.0)
            return;

        EL(a,c) =  af/tdt;
        ER(a,c) = -af/tdt;
        PL(a,c) =  dir==1 ? nxL(a,c) : nyL(a,c);
        PR(a,c) = -(dir==1 ? nxR(a,c) : nyR(a,c));
        K1(a,c) =  dir==1 ? nx0(a,c) : ny0(a,c);
        K2(a,c) =  K1(a,c) + qE;
        return;
        }

        const int k = dir==1 ? S.bfx[id(a,c)] : S.bfy[id(a,c)];

        if(k<0)
        return;

        const bc_face &f = bcf[k];
        const bool left = (f.i==a && f.j==c);      // the cell inside is on the left of the face
        slice &Ei = left ? EL : ER;

        if(f.type==bc_level && Hf>0.0)
        {
        Ei(a,c) = -af/tdt*f.s;
        K1(a,c) =  af/tdt*f.s*etaI(f.gi,f.gj);
        K2(a,c) =  K1(a,c) + qE;
        }

        if(f.type==bc_flux)
        K2(a,c) = f.s*f.qin;

        if(f.type==bc_riemann || f.type==bc_flather)
        {
        Ei(a,c) = -f.s*f.B;
        K2(a,c) =  f.s*f.A;
        }
    };

    SLICELOOP1
    kface(1,i,j,kxEL,kxER,kxPL,kxPR,kx1,kx2);

    if(yd)
    SLICELOOP2
    kface(2,i,j,kyEL,kyER,kyPL,kyPR,ky1,ky2);

    for(const bc_face &f : bcf)
    {
        if(f.cs==1 || f.cs==4)
        kface(1,f.fi,f.fj,kxEL,kxER,kxPL,kxPR,kx1,kx2);

        else
        kface(2,f.fi,f.fj,kyEL,kyER,kyPL,kyPR,ky1,ky2);
    }

    exchange1(p,pgc,kxEL); exchange1(p,pgc,kxER); exchange1(p,pgc,kxPL); exchange1(p,pgc,kxPR);
    exchange1(p,pgc,kx1); exchange1(p,pgc,kx2);
    exchange2(p,pgc,kyEL); exchange2(p,pgc,kyER); exchange2(p,pgc,kyPL); exchange2(p,pgc,kyPR);
    exchange2(p,pgc,ky1); exchange2(p,pgc,ky2);

    const bool bd = (Bdisp>0.0);

    // flat indices: cell n at slice index cq[n]; x neighbours +-J, y neighbours +-1; the east face of a cell has
    // its index, the west face q-J, the north face q, the south face q-1
    const int J = p->jmax;
    std::vector<int> cq(N);
    std::vector<double> idx(N), idy(N), lxe(N), lxw(N), lyn(N), lys(N), h3(N), w2(N), cst(N);

    {
    int n=0;
    SLICELOOP4
    {
    const int q = IJ;
    cq[n] = q;
    idx[n] = 1.0/p->DXN[IP];
    idy[n] = yd ? 1.0/p->DYN[JP] : 0.0;

    // h^3 lap of the A 220 3 term (non-hydrostatic cells, mirrored at the domain edges)
    const bool a = nA(i,j)>0.5 && bd;
    lxe[n] = (a && in(i+1,j)) ? 1.0/(p->DXP[IP]*p->DXN[IP]) : 0.0;
    lxw[n] = (a && in(i-1,j)) ? 1.0/(p->DXP[IM1]*p->DXN[IP]) : 0.0;
    lyn[n] = (a && yd && in(i,j+1)) ? 1.0/(p->DYP[JP]*p->DYN[JP]) : 0.0;
    lys[n] = (a && yd && in(i,j-1)) ? 1.0/(p->DYP[JM1]*p->DYN[JP]) : 0.0;
    h3[n] = a ? WLh(i,j)*WLh(i,j)*WLh(i,j) : 0.0;

    // constraint: 2 dt cw/H P, and the affine terms (W predictor, phi, slope term with U^{n+1} extrapolated:
    // implicit, the preconditioner below misses it - about 100 iterations on the Beji & Battjes bar instead of
    // 4, the errors on the bar the same)
    const double H = WLq(i,j);
    w2[n] = 2.0*dt*cwn/H;
    cst[n] = 2.0*(WHp(i,j) + dt*0.25*H*phiN(i,j))/H + nsx(i,j)*nUs(i,j) + (yd ? nsy(i,j)*nVs(i,j) : 0.0);
    ++n;
    }
    }

    double *E = nE.V, *Pp = nP.V, *Lp = nL.V;
    const double *XEL = kxEL.V, *XER = kxER.V, *XPL = kxPL.V, *XPR = kxPR.V, *XB = nxB.V, *X1 = kx1.V, *X2 = kx2.V;
    const double *YEL = kyEL.V, *YER = kyER.V, *YPL = kyPL.V, *YPR = kyPR.V, *YB = nyB.V, *Y1 = ky1.V, *Y2 = ky2.V;
    const double *Ac = nA.V;
    double *PW = nPW.V;

    // wet pattern of the mass rows (A 229 1: Newton iteration on WL = max(eta + d, 0), Casulli 2009), all wet first
    SLICELOOP4
    nPW(i,j) = 1.0;

    // x-face q (between q and q+J), y-face q (between q and q+1); mode 1: the constants of the mass row (the q^E of
    // the faces are in rhsW), mode 2: all constants
    auto qxf = [&](int q, int mode)
    {
        double v = XEL[q]*E[q] + XER[q]*E[q+J] + XPL[q]*Pp[q] + XPR[q]*Pp[q+J];
        if(bd) v += XB[q]*(Lp[q+J] - Lp[q]);
        return v + (mode==1 ? X1[q] : (mode==2 ? X2[q] : 0.0));
    };

    auto qyf = [&](int q, int mode)
    {
        double v = YEL[q]*E[q] + YER[q]*E[q+1] + YPL[q]*Pp[q] + YPR[q]*Pp[q+1];
        if(bd) v += YB[q]*(Lp[q+1] - Lp[q]);
        return v + (mode==1 ? Y1[q] : (mode==2 ? Y2[q] : 0.0));
    };

    auto lapL = [&]()
    {
        for(int n=0; n<N; ++n)
        {
        const int q = cq[n];
        Lp[q] = h3[n]*(lxe[n]*(E[q+J] - E[q]) - lxw[n]*(E[q] - E[q-J]) + lyn[n]*(E[q+1] - E[q]) - lys[n]*(E[q] - E[q-1]));
        }
    };

    auto apply = [&](const std::vector<double> &x, std::vector<double> &y, bool affine)
    {
        for(int n=0; n<N; ++n)
        {
        E[cq[n]]  = x[2*n];
        Pp[cq[n]] = x[2*n+1];
        }

        exchange4(p,pgc,nE);
        exchange4(p,pgc,nP);

        if(bd)
        {
        lapL();
        exchange4(p,pgc,nL);
        }

        const int m1 = affine ? 1 : 0, m2 = affine ? 2 : 0;

        for(int n=0; n<N; ++n)
        {
        const int q = cq[n];

        // mass row
        double dv = (qxf(q,m1) - qxf(q-J,m1))*idx[n];
        if(yd) dv += (qyf(q,m1) - qyf(q-1,m1))*idy[n];

        y[2*n] = PW[q]*E[q] + tdt*dv;

            // constraint, or P = 0
            if(Ac[q]<0.5)
            y[2*n+1] = Pp[q];

            else
            {
            double cc = (qxf(q,m2) - qxf(q-J,m2))*idx[n];
            if(yd) cc += (qyf(q,m2) - qyf(q-1,m2))*idy[n];

            cc += w2[n]*Pp[q];

            if(affine)
            cc += cst[n];

            y[2*n+1] = cc;
            }
        }
    };

    // initial guess: eta and P extrapolated from t^n and t^{n-1} (fewer Krylov iterations)
    {
    const double r = dt/MAX(p->dt_old,1.0e-20);
    int n=0;
    SLICELOOP4
    {
    const double P0 = b->press(i,j)/p->W1;
    xv[2*n]   = nh_x0 ? (1.0+r)*eta0(i,j) - r*nEo(i,j) : eta0(i,j);
    xv[2*n+1] = nA(i,j)>0.5 ? (nh_x0 ? (1.0+r)*P0 - r*nPo(i,j) : P0) : 0.0;
    nEo(i,j) = eta0(i,j);
    nPo(i,j) = P0;
    ++n;
    }
    nh_x0 = true;
    }

    std::vector<double> rv(2*N), zv(2*N);
    nh_wd_iter = 0;
    nh_wd_dry = 0;

    for(int wdit=0; wdit<(p->A229==1 ? 20 : 1); ++wdit)
    {
    ++nh_wd_iter;

    // right-hand side: b = (rhsW - d, 0) - A(0)
    std::fill(zv.begin(),zv.end(),0.0);
    apply(zv,yv,true);

    {
    int n=0;
    SLICELOOP4
    {
    rv[2*n]   = rhsW(i,j) - nPW(i,j)*b->depth(i,j) - yv[2*n];
    rv[2*n+1] = -yv[2*n+1];
    ++n;
    }
    }

    nh_iter = 0;

    if(!yd && p->mpi_size==1 && !nh_krylov)
    {
        // ---- band LU, the operator probed with five colours (couplings reach two cells)
        if(S.n!=2*N)
        S.resize(2*N);
        else
        std::fill(S.AB.begin(),S.AB.end(),0.0);

        std::vector<int> ci(N);
        {
        int n=0;
        SLICELOOP4
        ci[n++] = i;
        }

        std::vector<int> cell(p->imax*p->jmax,-1);
        {
        int n=0;
        SLICELOOP4
        cell[IJ] = n++;
        }

        std::vector<double> e(2*N), ye(2*N);

        for(int col=0; col<5; ++col)
        for(int t=0; t<2; ++t)
        {
            std::fill(e.begin(),e.end(),0.0);

            for(int n=0; n<N; ++n)
            if(((ci[n]%5)+5)%5==col)
            e[2*n+t] = 1.0;

            apply(e,ye,false);

            for(int n=0; n<N; ++n)
            for(int tr=0; tr<2; ++tr)
            {
                const double v = ye[2*n+tr];
                if(v==0.0)
                continue;

                // the column of this colour within two cells
                for(int d=-2; d<=2; ++d)
                {
                const int ii = ci[n]+d;
                if(((ii%5)+5)%5!=col || ii<p->imin || ii>=p->imin+p->imax)
                continue;

                const int m = cell[(ii-p->imin)*p->jmax + (0-p->jmin)];
                if(m<0)
                continue;

                const int rr = 2*n+tr, cc = 2*m+t;
                if(cc-rr>=-S.kl && cc-rr<=S.ku)
                S.at(rr,cc) = v;
                }
            }
        }

        if(!S.factor())
        {
        if(p->mpirank==0)
        cout<<"!!! SFLOW A 226 1 with A 220 "<<p->A220<<": the coupled eta-P system is singular !!!"<<endl;
        exit(1);
        }

        xv = rv;
        S.solve(xv);
        nh_iter = 1;
    }

    else
    {
        // ---- FGMRES, right preconditioner: eta eliminated locally, one V-cycle for P (eta in hydrostatic cells)
        // psi-coefficients of the cells: E = bE + aE psi, P = bP + aP psi
        SLICELOOP4
        {
            if(nA(i,j)>0.5)
            {
            nAE(i,j) = tdt*2.0*dt*cwn/WLq(i,j);
            nAP(i,j) = 1.0;
            }

            else
            {
            nAE(i,j) = 1.0;
            nAP(i,j) = 0.0;
            }
        }

        exchange4(p,pgc,nAE);
        exchange4(p,pgc,nAP);

        // the 5-point matrix of psi (rows in SLICELOOP4 order)
        if(!S.mg1)
        {
        S.mg1.reset(new reefmg2D(p,pgc));
        S.M1.reset(new matrix2D(p));
        S.r1.reset(new vec2D(p));

            if(!yd)
            S.mg1->set_exact_line(true);
        }

        if(Bdisp>0.0 && !S.mg2)
        {
        S.mg2.reset(new reefmg2D(p,pgc));
        S.M2.reset(new matrix2D(p));
        S.r2.reset(new vec2D(p));

            if(!yd)
            S.mg2->set_exact_line(true);
        }

        matrix2D &M = *S.M1;

        auto face_pre = [&](int dir, int a, int c, int a1, int c1, bool self_left, double dn, double &diag, double &off)
        {
            // contribution of face (a,c)-(a1,c1) to the row of the cell on side self_left (true: the left cell)
            diag = off = 0.0;
            const double Hf = dir==1 ? Hx(a,c) : Hy(a,c);
            const double af = dir==1 ? ax(a,c) : ay(a,c);

            if(in(a,c) && in(a1,c1))
            {
                if(Hf<=0.0)
                return;

            const double cR = dir==1 ? nxR(a,c) : nyR(a,c);
            const double cL = dir==1 ? nxL(a,c) : nyL(a,c);
            const int si = self_left ? a : a1, sj = self_left ? c : c1;
            const int oi = self_left ? a1 : a, oj = self_left ? c1 : c;
            const double cs = self_left ? cL : cR, co = self_left ? cR : cL;

            diag = (af/dn)*nAE(si,sj) + tdt*cs/dn*nAP(si,sj);
            off  = -(af/dn)*nAE(oi,oj) - tdt*co/dn*nAP(oi,oj);
            return;
            }

            // open boundary faces: implicit part on the diagonal
            const int k = dir==1 ? S.bfx[id(a,c)] : S.bfy[id(a,c)];
            if(k<0)
            return;

            const bc_face &f = bcf[k];
            const int si = f.i, sj = f.j;

            if(f.type==bc_level && Hf>0.0)
            diag = (af/dn)*nAE(si,sj);

            if(f.type==bc_riemann || f.type==bc_flather)
            diag = tdt*f.B/dn*nAE(si,sj);
        };

        {
        int n=0;
        SLICELOOP4
        {
        double d, o;
        M.p[n] = nPW(i,j)*nAE(i,j);
        M.n[n] = M.s[n] = M.w[n] = M.e[n] = 0.0;

        face_pre(1,i,j,i+1,j,true,p->DXN[IP],d,o);    M.p[n] += d; M.n[n] = o;
        face_pre(1,i-1,j,i,j,false,p->DXN[IP],d,o);   M.p[n] += d; M.s[n] = o;

            if(yd)
            {
            face_pre(2,i,j,i,j+1,true,p->DYN[JP],d,o);  M.p[n] += d; M.w[n] = o;
            face_pre(2,i,j-1,i,j,false,p->DYN[JP],d,o); M.p[n] += d; M.e[n] = o;
            }
        ++n;
        }
        }

        // A 220 3: the h^3 lap(eta) term makes the reduced operator fourth order, M = M5 + D lap^2; with
        // M ~ M5 (I - rho lap), rho = D/K (K: the 5-point coefficient of the face), the preconditioner is
        // A2^-1 M5^-1 with A2 = I - div(rho grad)
        if(Bdisp>0.0)
        {
            matrix2D &M2 = *S.M2;

            auto rho = [&](int dir, int a, int c, int a1, int c1)
            {
                const double Hf = dir==1 ? Hx(a,c) : Hy(a,c);

                if(!(ac(a,c) && ac(a1,c1)) || Hf<=0.0)
                return 0.0;

                const double af = dir==1 ? ax(a,c) : ay(a,c);
                const double cB = dir==1 ? nxB(a,c) : nyB(a,c);
                const double cc = 0.5*(dir==1 ? nxR(a,c) + nxL(a,c) : nyR(a,c) + nyL(a,c));
                const double aE = 0.5*(nAE(a,c) + nAE(a1,c1));
                const double h3 = 0.5*(pow(WLh(a,c),3.0) + pow(WLh(a1,c1),3.0));
                const double K = af*aE + tdt*cc;

                return K>0.0 ? tdt*cB*h3*aE/K : 0.0;
            };

            int n=0;
            SLICELOOP4
            {
            const double re = rho(1,i,j,i+1,j)/(p->DXP[IP]*p->DXN[IP]);
            const double rw = rho(1,i-1,j,i,j)/(p->DXP[IM1]*p->DXN[IP]);
            double rn = 0.0, rs = 0.0;

                if(yd)
                {
                rn = rho(2,i,j,i,j+1)/(p->DYP[JP]*p->DYN[JP]);
                rs = rho(2,i,j-1,i,j)/(p->DYP[JM1]*p->DYN[JP]);
                }

            M2.p[n] = 1.0 + re + rw + rn + rs;
            M2.n[n] = -re;
            M2.s[n] = -rw;
            M2.w[n] = -rn;
            M2.e[n] = -rs;
            ++n;
            }
        }

        bool rebuild = true;

        // z = M^-1 v for the coupled system
        auto precond = [&](const std::vector<double> &v, std::vector<double> &z)
        {
            int n=0;
            SLICELOOP4
            {
            const double rm = v[2*n], rc = v[2*n+1];

                if(nA(i,j)>0.5)
                {
                nBE(i,j) = rm - tdt*rc;
                nBP(i,j) = 0.0;
                }

                else
                {
                nBE(i,j) = 0.0;
                nBP(i,j) = rc;
                }
            ++n;
            }

            // right-hand side: rm minus the mass row of (bE, bP), without the A 220 3 term (nL = 0)
            for(int n=0; n<N; ++n)
            {
            const int q = cq[n];
            E[q]  = nBE.V[q];
            Pp[q] = nBP.V[q];
            Lp[q] = 0.0;
            }

            exchange4(p,pgc,nE);
            exchange4(p,pgc,nP);
            exchange4(p,pgc,nL);

            for(int n=0; n<N; ++n)
            {
            const int q = cq[n];
            double dv = (qxf(q,0) - qxf(q-J,0))*idx[n];
            if(yd) dv += (qyf(q,0) - qyf(q-1,0))*idy[n];

            S.r1->V[n] = v[2*n] - (E[q] + tdt*dv);
            nZ.V[q] = 0.0;
            }

            S.mg1->vcycle(p,pgc,nZ,M,*S.r1,rebuild);

            if(Bdisp>0.0)
            {
            n=0;
            SLICELOOP4
            S.r2->V[n++] = nZ(i,j);

            S.mg2->vcycle(p,pgc,nZ,*S.M2,*S.r2,rebuild);
            }

            n=0;
            SLICELOOP4
            {
            z[2*n]   = nBE(i,j) + nAE(i,j)*nZ(i,j);
            z[2*n+1] = nBP(i,j) + nAP(i,j)*nZ(i,j);
            ++n;
            }
        };

        auto dot = [&](const std::vector<double> &a, const std::vector<double> &c)
        {
            double s = 0.0;
            for(int k=0; k<2*N; ++k)
            s += a[k]*c[k];
            return pgc->globalsum(s);
        };

        const int m = 30, maxit = 300;
        // relative tolerance: 1e-6 and 1e-8 gave the same gauge records (8 digits) on the 2D shoal case
        const double tol = 1.0e-6;

        S.V.assign(m+1,std::vector<double>(2*N));
        S.Z.assign(m,std::vector<double>(2*N));
        S.H.assign((m+1)*m,0.0);
        S.cs.assign(m,0.0); S.sn.assign(m,0.0); S.g.assign(m+1,0.0);

        const double bn = sqrt(dot(rv,rv));
        int it = 0;

        while(it<maxit)
        {
            apply(xv,yv,false);

            for(int k=0; k<2*N; ++k)
            S.V[0][k] = rv[k] - yv[k];

            double beta = sqrt(dot(S.V[0],S.V[0]));

            if(beta<=tol*MAX(bn,1.0e-30))
            break;

            for(int k=0; k<2*N; ++k)
            S.V[0][k] /= beta;

            std::fill(S.g.begin(),S.g.end(),0.0);
            S.g[0] = beta;

            int jj=0;
            for(; jj<m && it<maxit; ++jj, ++it)
            {
                precond(S.V[jj],S.Z[jj]);
                apply(S.Z[jj],S.V[jj+1],false);

                for(int l=0; l<=jj; ++l)
                {
                const double h = dot(S.V[jj+1],S.V[l]);
                S.H[l*m+jj] = h;
                for(int k=0; k<2*N; ++k)
                S.V[jj+1][k] -= h*S.V[l][k];
                }

                const double h1 = sqrt(dot(S.V[jj+1],S.V[jj+1]));
                S.H[(jj+1)*m+jj] = h1;

                if(h1>0.0)
                for(int k=0; k<2*N; ++k)
                S.V[jj+1][k] /= h1;

                for(int l=0; l<jj; ++l)
                {
                const double t0 = S.cs[l]*S.H[l*m+jj] + S.sn[l]*S.H[(l+1)*m+jj];
                S.H[(l+1)*m+jj] = -S.sn[l]*S.H[l*m+jj] + S.cs[l]*S.H[(l+1)*m+jj];
                S.H[l*m+jj] = t0;
                }

                const double a0 = S.H[jj*m+jj], a1 = S.H[(jj+1)*m+jj];
                const double rr = sqrt(a0*a0 + a1*a1);
                S.cs[jj] = rr>0.0 ? a0/rr : 1.0;
                S.sn[jj] = rr>0.0 ? a1/rr : 0.0;
                S.H[jj*m+jj] = rr;
                S.H[(jj+1)*m+jj] = 0.0;
                S.g[jj+1] = -S.sn[jj]*S.g[jj];
                S.g[jj]   =  S.cs[jj]*S.g[jj];

                if(fabs(S.g[jj+1])<=tol*MAX(bn,1.0e-30))
                {
                ++jj; ++it;
                break;
                }
            }

            // x += Z y,  H y = g
            std::vector<double> yk(jj,0.0);
            for(int l=jj-1; l>=0; --l)
            {
            double s = S.g[l];
            for(int q=l+1; q<jj; ++q)
            s -= S.H[l*m+q]*yk[q];
            yk[l] = s/S.H[l*m+l];
            }

            for(int l=0; l<jj; ++l)
            for(int k=0; k<2*N; ++k)
            xv[k] += yk[l]*S.Z[l][k];

            if(jj>0 && fabs(S.g[jj])<=tol*MAX(bn,1.0e-30))
            break;
        }

        nh_iter = it;
    }

        if(p->A229!=1)
        break;

        // new wet pattern (si_eta2D::wd_update): dry where eta + d <= 0 and the cell has an open face
        int change=0;
        {
        int n=0;
        SLICELOOP4
        {
            const bool coupled = ax(i,j)>0.0 || ax(i-1,j)>0.0 || bdiag(i,j)>0.0
                              || (yd && (ay(i,j)>0.0 || ay(i,j-1)>0.0));
            const double pn = (xv[2*n] + b->depth(i,j)>0.0 || !coupled) ? 1.0 : 0.0;

            if(pn!=nPW(i,j))
            {
            nPW(i,j) = pn;
            ++change;
            }
        ++n;
        }
        }

        change = pgc->globalisum(change);

        if(change==0)
        break;
    }

    SLICELOOP4
    if(nPW(i,j)<0.5)
    ++nh_wd_dry;

    nh_wd_dry = pgc->globalisum(nh_wd_dry);

    // ---- the solution; nE, nP, nL of the solution (halo included) for the face terms
    apply(xv,yv,false);

    {
    int n=0;
    SLICELOOP4
    {
    etaI(i,j) = xv[2*n];
    b->press(i,j) = nA(i,j)>0.5 ? p->W1*xv[2*n+1] : 0.0;
    ++n;
    }
    }

    auto nhq = [&](int dir, int a, int c)
    {
        const int a1 = dir==1 ? a+1 : a, c1 = dir==1 ? c : c+1;
        const double Hf = dir==1 ? Hx(a,c) : Hy(a,c);

        if(!(in(a,c) && in(a1,c1)) || Hf<=0.0)
        return 0.0;

        const double cR = dir==1 ? nxR(a,c) : nyR(a,c);
        const double cL = dir==1 ? nxL(a,c) : nyL(a,c);
        const double cB = dir==1 ? nxB(a,c) : nyB(a,c);
        const double c0 = dir==1 ? nx0(a,c) : ny0(a,c);

        double q = -cR*nP(a1,c1) + cL*nP(a,c) + c0;

        if(Bdisp>0.0)
        q += cB*(nL(a1,c1) - nL(a,c));

        return q;
    };

    SLICELOOP1
    qxN(i,j) = nhq(1,i,j);

    if(yd)
    SLICELOOP2
    qyN(i,j) = nhq(2,i,j);

    exchange1(p,pgc,qxN);
    exchange2(p,pgc,qyN);
    exchange4(p,pgc,etaI);      // face_update and cell_update read eta^{n+1} of the halo
    pgc->gcsl_start4(p,b->press,40);

    nh_time = pgc->timer() - starttime_nh;
    psi->iter = nh_iter;
    psi->time = nh_time;
}

void sflow_momentum_SI::nh_cell(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // UH, VH: the mean of the non-hydrostatic face terms of the open faces (as the gravity term); WH^{n+1}
    const double dt = p->dt;
    auto id = [&](int a, int c) { return (a-p->imin)*p->jmax + (c-p->jmin); };
    auto in = [&](int a, int c) { return p->flagslice4[id(a,c)]>0; };
    auto wn = [&](int a, int c) { return etaI(a,c) + b->depth(a,c) > p->A244; };

    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        const double we = (Hx(i,j)>0.0 && in(i+1,j) && wn(i+1,j)) ? 1.0 : 0.0;
        const double ww = (Hx(i-1,j)>0.0 && in(i-1,j) && wn(i-1,j)) ? 1.0 : 0.0;

        b->UH(i,j) += (ww*qxN(i-1,j) + we*qxN(i,j))/MAX(we+ww,1.0);

            if(p->j_dir==1)
            {
            const double wn_ = (Hy(i,j)>0.0 && in(i,j+1) && wn(i,j+1)) ? 1.0 : 0.0;
            const double ws_ = (Hy(i,j-1)>0.0 && in(i,j-1) && wn(i,j-1)) ? 1.0 : 0.0;

            b->VH(i,j) += (ws_*qyN(i,j-1) + wn_*qyN(i,j))/MAX(wn_+ws_,1.0);
            }

        b->WH(i,j) = WHp(i,j);

        if(nA(i,j)>0.5)
        b->WH(i,j) += dt*(cwn*b->press(i,j)/p->W1 + 0.25*WLq(i,j)*phiN(i,j));
        }

        else
        b->WH(i,j) = 0.0;
    }
}
