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
// SFLOW semi-implicit step (A 226 1) with the quadratic non-hydrostatic pressure (A 220 2/3): prototype.
//
// The pressure of sflow_pjm_quad (Jeschke et al. 2017; A 220 3 with the improved dispersion, alpha A 224) is not
// projected after the step, but solved together with eta^{n+1} in one linear system, so that neither sqrt(g h)
// nor the dispersive wave speed limits the time step (as the theta = 1/2 coupled solve of Basilisk's layered
// solver, Popinet 2020). Per cell the unknowns eta^{n+1} and P^{n+1} = p/rho, P fully implicit:
//
//  face  q_f = q^E_f - a_f/(theta dt) (eta_R - eta_L)                                      (hydrostatic SI)
//              - dt alpha_f [ (H_R P_R - H_L P_L) - qb_f (d_R - d_L) ] / dxp_f              (pressure, bed pressure)
//              + dt alpha_f g B [ theta (L_R - L_L)^{n+1} + (1-theta) (L_R - L_L)^n ] / dxp_f  (A 220 3)
//        qb_f = cb (P_L + P_R)/2 + alpha (H phi)_f / 4,  L = h^3 lap(eta),  B = (alpha - 1)/3
//  mass  (1 + bdiag) eta_i - div(a grad eta) + theta dt div(q_NH) = rhsW - d               (row of si_eta2D)
//  incompressibility (new velocities)
//        div(q) + U (2 d_x - H_x) + 2 W = 0,  U the mean of q_f/H_f, W = (WH* + dt (cw P + H phi/4))/H
//
// compact face differences throughout (no checkerboard mode, the dispersion of the staggered grid). The cell
// discharges take the mean of the non-hydrostatic face terms (as the gravity term, cell_update). phi, the bed
// acceleration Dw_b/Dt, is explicit: the acceleration of the previous step (an estimate of U^{n+1} with P^n as in
// sflow_pjm_quad is unstable for large steps with A 220 3). WH*: WH^n advected upwind with q^n. H of the step: the
// water depth of t^n.
//
// Limits of the prototype (check()): one grid, 1D (one cell in y), Cartesian, serial, no wetting and drying
// iteration; sparse direct solve (Eigen SparseLU). Non-hydrostatic cells: wet, deeper than 10 A 244, not
// breaking, wet neighbours; P = 0 elsewhere. Open boundary faces carry no pressure term.
// ---------------------------------------------------------------------------

#include"sflow_momentum_SI.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"si_eta2D.h"
#include<Eigen/Sparse>
#include<Eigen/SparseLU>
#include<vector>
#include<memory>

// sparse LU of the coupled system: the symbolic analysis is kept while the pattern does not change
struct sflow_si_nh_lu
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::NaturalOrdering<int>> lu;
    std::vector<int> outer, inner;
    bool analysed = false;
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
    }

    SLICELOOP1
    qxN(i,j) = 0.0;
}

int sflow_momentum_SI::nh_active(lexer *p, fdm2D *b, int ii, int jj)
{
    // wet, deep, not breaking, the neighbours in the domain wet (as sflow_pjm_quad::active, a wall is no dry cell)
    auto id = [&](int a, int c) { return (a-p->imin)*p->jmax + (c-p->jmin); };

    if(p->wet[id(ii,jj)]==0 || WL0(ii,jj)<=10.0*p->A244 || b->breaking(ii,jj)>0)
    return 0;

    if(p->flagslice4[id(ii-1,jj)]>0 && p->wet[id(ii-1,jj)]==0)
    return 0;

    if(p->flagslice4[id(ii+1,jj)]>0 && p->wet[id(ii+1,jj)]==0)
    return 0;

    return 1;
}

void sflow_momentum_SI::nh_predictor(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // h^3 lap(eta^n) (A 220 3), WH predictor, bed acceleration phi
    const double dt = p->dt;
    auto id = [&](int a, int c) { return (a-p->imin)*p->jmax + (c-p->jmin); };
    auto in = [&](int a, int c) { return p->flagslice4[id(a,c)]>0; };

    // L = h^3 lap(eta), compact, mirrored at the domain edges
    SLICELOOP4
    {
    Lc0(i,j) = 0.0;

        if(Bdisp>0.0 && nh_active(p,b,i,j)==1)
        {
        const double h = WL0(i,j);
        const double ge = in(i+1,j) ? (eta0(i+1,j) - eta0(i,j))/p->DXP[IP] : 0.0;
        const double gw = in(i-1,j) ? (eta0(i,j) - eta0(i-1,j))/p->DXP[IM1] : 0.0;
        Lc0(i,j) = h*h*h*(ge - gw)/p->DXN[IP];
        }
    }

    // W^n advected upwind with q^n (conservative, first order)
    auto Wn = [&](int a, int c) { return (in(a,c) && WL0(a,c)>p->A244) ? b->WH(a,c)/WL0(a,c) : 0.0; };
    auto fW = [&](int a, int c) { const double q = qxn(a,c); return q>=0.0 ? q*Wn(a,c) : q*Wn(a+1,c); };

    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        WHp(i,j) = b->WH(i,j) - dt*(fW(i,j) - fW(i-1,j))/p->DXN[IP];

        else
        WHp(i,j) = 0.0;
    }

    // U^n of the cells; the acceleration of the previous step, (U^n - U^{n-1})/dt, gives phi: the estimate of
    // U^{n+1} with the explicit terms of the step (sflow_pjm_quad, P^n) is not balanced by the pressure for large
    // time steps, the h^3 lap(eta) term of A 220 3 then drives a grid mode on a sloping bed (Beji & Battjes bar,
    // dt 0.005-0.02 s); solving again with phi from U^{n+1} changed nothing there
    SLICELOOP4
    {
    tU(i,j) = (p->wet[IJ]==1 && WL0(i,j)>p->A244) ? UH0(i,j)/WL0(i,j) : 0.0;

        if(nh_first)
        Unm1(i,j) = tU(i,j);
    }

    nh_first = false;

    // bed acceleration phi = Dw_b/Dt, w_b = -u d_x (sflow_pjm_quad::quad_calc)
    auto Ue = [&](int a, int c) { return in(a,c) ? tU(a,c) : b->U(a,c); };
    auto dep = [&](int a, int c, int a0) { return in(a,c) ? b->depth(a,c) : b->depth(a0,c); };

    SLICELOOP4
    {
    phiN(i,j) = 0.0;

        if(nh_active(p,b,i,j)==1)
        {
        const double dx = p->DXP[IP] + p->DXP[IM1];
        const double dE = dep(i+1,j,i), dW = dep(i-1,j,i), dC = b->depth(i,j);

        const double dddx  = (dE - dW)/dx;
        const double dddxx = ((dE - dC)/p->DXP[IP] - (dC - dW)/p->DXP[IM1])/p->DXN[IP];

        const double u = tU(i,j);

        const double dudt = (u - Unm1(i,j))/p->dt_old + u*(Ue(i+1,j) - Ue(i-1,j))/dx;

        phiN(i,j) = -dudt*dddx - u*u*dddxx;
        }
    }

    SLICELOOP4
    Unm1(i,j) = tU(i,j);
}

void sflow_momentum_SI::nh_solve(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // the coupled eta-P system, unknowns 2n (eta) and 2n+1 (P) of cell n (SLICELOOP4 order)
    const double starttime_nh = pgc->timer();

    const double dt = p->dt;
    const double tdt = theta*dt;

    auto id = [&](int a, int c) { return (a-p->imin)*p->jmax + (c-p->jmin); };
    auto in = [&](int a, int c) { return p->flagslice4[id(a,c)]>0; };

    std::vector<int> cell(p->imax*p->jmax,-1);
    std::vector<char> act(p->imax*p->jmax,0);

    int N=0;
    nh_cells_active = 0;

    SLICELOOP4
    {
    cell[IJ] = N++;
    act[IJ] = (char)nh_active(p,b,i,j);
    nh_cells_active += act[IJ];
    }

    auto ce = [&](int a, int c) { return in(a,c) ? cell[id(a,c)] : -1; };
    auto ac = [&](int a, int c) { return in(a,c) && act[id(a,c)]==1; };

    // linear form of a face quantity: sum c_k x_k + c0
    struct lin
    {
        int k[16]; double c[16]; int n=0; double c0=0.0;
        void add(int kk, double cc) { if(kk>=0) { k[n]=kk; c[n]=cc; ++n; } }
    };

    // L_c^{n+1} = h^3 lap(eta^{n+1}) of cell (a,c), times s, into f
    auto addL = [&](lin &f, int a, int c, double s)
    {
        if(Bdisp<=0.0 || !ac(a,c))
        return;

        const double h3 = WL0(a,c)*WL0(a,c)*WL0(a,c)/p->DXN[a+marge];
        const int n0 = cell[id(a,c)];

        if(in(a+1,c))
        {
        const double w = h3/p->DXP[a+marge];
        f.add(2*ce(a+1,c), s*w);
        f.add(2*n0, -s*w);
        }

        if(in(a-1,c))
        {
        const double w = h3/p->DXP[a-1+marge];
        f.add(2*n0, -s*w);
        f.add(2*ce(a-1,c), s*w);
        }
    };

    // non-hydrostatic part of the discharge of the interior face a+1/2
    auto fNH = [&](int a, int c)
    {
        lin f;

        if(!(in(a,c) && in(a+1,c)) || Hx(a,c)<=0.0)
        return f;

        const double dxp = p->DXP[a+marge];
        const double dd = b->depth(a+1,c) - b->depth(a,c);
        const double HL = WL0(a,c), HR = WL0(a+1,c);
        const double s = dt*alx(a,c)/dxp;
        const int L = ce(a,c), R = ce(a+1,c);

        // pressure and bed pressure (P of a hydrostatic cell is 0: its row is P = 0)
        f.add(2*R+1, -s*(HR - 0.5*cb*dd));
        f.add(2*L+1,  s*(HL + 0.5*cb*dd));
        f.c0 += s*0.125*adisp*(HL*phiN(a,c) + HR*phiN(a+1,c))*dd;

        // dispersion term (A 220 3), theta-weighted
        if(Bdisp>0.0)
        {
        const double sb = s*g*Bdisp;
        addL(f,a+1,c, sb*theta);
        addL(f,a,c, -sb*theta);
        f.c0 += sb*(1.0-theta)*(Lc0(a+1,c) - Lc0(a,c));
        }

        return f;
    };

    // full discharge of the face a+1/2 at t^{n+1} (interior and open boundary faces, walls 0)
    auto fQ = [&](int a, int c)
    {
        lin f;

        if(in(a,c) && in(a+1,c))
        {
            if(Hx(a,c)<=0.0)
            return f;

        f = fNH(a,c);
        f.c0 += qxE(a,c);
        f.add(2*ce(a+1,c), -ax(a,c)/tdt);
        f.add(2*ce(a,c), ax(a,c)/tdt);
        return f;
        }

        for(const bc_face &bf : bcf)
        if((bf.cs==1 || bf.cs==4) && bf.fi==a && bf.fj==c)
        {
            const int n0 = cell[id(bf.i,bf.j)];

            if(bf.type==bc_level && Hx(a,c)>0.0)
            {
            f.c0 = qxE(a,c) + ax(a,c)/tdt*bf.s*etaI(bf.gi,bf.gj);
            f.add(2*n0, -ax(a,c)/tdt*bf.s);
            }

            if(bf.type==bc_flux)
            f.c0 = bf.s*bf.qin;

            if(bf.type==bc_riemann || bf.type==bc_flather)
            {
            f.c0 = bf.s*bf.A;
            f.add(2*n0, -bf.s*bf.B);
            }
        }

        return f;
    };

    std::vector<Eigen::Triplet<double>> T;
    T.reserve(32*N);
    Eigen::VectorXd r(2*N), x(2*N);

    auto row = [&](int rr, const lin &f, double s)
    {
        for(int q=0; q<f.n; ++q)
        T.emplace_back(rr,f.k[q],s*f.c[q]);

        r[rr] -= s*f.c0;
    };

    auto put = [&](int rr, int cc, double v) { T.emplace_back(rr,cc,v); };

    auto build = [&]()
    {
    SLICELOOP4
    {
    const int n0 = cell[IJ];
    const int re = 2*n0, rp = 2*n0+1;

    r[re] = 0.0;
    r[rp] = 0.0;

    // mass: the row of si_eta2D (wet pattern 1), plus theta dt div(q_NH)
    const double cn = ax(i,j)/p->DXN[IP];
    const double cs = ax(i-1,j)/p->DXN[IP];

    put(re,re,1.0 + bdiag(i,j) + cn + cs);
    r[re] += rhsW(i,j) - b->depth(i,j);

        if(cn!=0.0)
        {
        if(in(i+1,j))
        put(re,2*ce(i+1,j),-cn);
        else
        r[re] += cn*etaI(i+1,j);
        }

        if(cs!=0.0)
        {
        if(in(i-1,j))
        put(re,2*ce(i-1,j),-cs);
        else
        r[re] += cs*etaI(i-1,j);
        }

    row(re, fNH(i,j),    tdt/p->DXN[IP]);
    row(re, fNH(i-1,j), -tdt/p->DXN[IP]);

    // incompressibility of the new velocities, or P = 0
        if(act[IJ]==0)
        {
        put(rp,rp,1.0);
        continue;
        }

    const double H = WL0(i,j);
    const lin qe = fQ(i,j), qw = fQ(i-1,j);

    row(rp, qe,  1.0/p->DXN[IP]);
    row(rp, qw, -1.0/p->DXN[IP]);

    // U (2 d_x - H_x), U the mean of the open face velocities
    const double dx = p->DXP[IP] + p->DXP[IM1];
    const double dE = in(i+1,j) ? b->depth(i+1,j) : b->depth(i,j);
    const double dW = in(i-1,j) ? b->depth(i-1,j) : b->depth(i,j);
    const double HE = in(i+1,j) ? WL0(i+1,j) : H;
    const double HW = in(i-1,j) ? WL0(i-1,j) : H;
    const double sl = (2.0*(dE - dW) - (HE - HW))/dx;

    const double fe = Hx(i,j)>0.0 ? 1.0 : 0.0;
    const double fw = Hx(i-1,j)>0.0 ? 1.0 : 0.0;

        if(fe+fw>0.0 && fabs(sl)>0.0)
        {
        if(fe>0.0)
        row(rp, qe, sl/(fe+fw)/Hx(i,j));

        if(fw>0.0)
        row(rp, qw, sl/(fe+fw)/Hx(i-1,j));
        }

    // 2 W^{n+1}
    put(rp,rp,2.0*dt*cwn/H);
    r[rp] -= 2.0*(WHp(i,j) + dt*0.25*H*phiN(i,j))/H;
    }
    };

    build();

    Eigen::SparseMatrix<double> A(2*N,2*N);
    A.setFromTriplets(T.begin(),T.end());

    // unknowns interleaved per cell: a band matrix, natural ordering
    A.makeCompressed();

    if(!nhlu)
    nhlu = std::make_shared<sflow_si_nh_lu>();

    auto &S = *nhlu;
    const int no = A.outerSize()+1, ni = A.nonZeros();

    bool same = S.analysed && (int)S.outer.size()==no && (int)S.inner.size()==ni
             && std::equal(S.outer.begin(),S.outer.end(),A.outerIndexPtr())
             && std::equal(S.inner.begin(),S.inner.end(),A.innerIndexPtr());

    if(!same)
    {
    S.lu.analyzePattern(A);
    S.outer.assign(A.outerIndexPtr(),A.outerIndexPtr()+no);
    S.inner.assign(A.innerIndexPtr(),A.innerIndexPtr()+ni);
    S.analysed = true;
    }

    S.lu.factorize(A);
    auto &lu = S.lu;

    if(lu.info()!=Eigen::Success)
    {
    if(p->mpirank==0)
    cout<<"!!! SFLOW A 226 1 with A 220 "<<p->A220<<": the coupled eta-P system is singular !!!"<<endl;
    exit(1);
    }

    x = lu.solve(r);

    auto val = [&](const lin &f) { double q = f.c0; for(int k=0; k<f.n; ++k) q += f.c[k]*x[f.k[k]]; return q; };

    SLICELOOP4
    {
    etaI(i,j) = x[2*cell[IJ]];
    b->press(i,j) = act[IJ]==1 ? p->W1*x[2*cell[IJ]+1] : 0.0;
    }

    // the non-hydrostatic face discharges of the solution
    SLICELOOP1
    {
    qxN(i,j) = 0.0;

        if(in(i,j) && in(i+1,j) && Hx(i,j)>0.0)
        qxN(i,j) = val(fNH(i,j));
    }

    exchange1(p,pgc,qxN);
    pgc->gcsl_start4(p,b->press,40);

    nh_time = pgc->timer() - starttime_nh;
    psi->iter = 1;
    psi->time = nh_time;
}

void sflow_momentum_SI::nh_cell(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // UH: the mean of the non-hydrostatic face terms of the open faces (as the gravity term); WH^{n+1}
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
        const double ws = MAX(we+ww,1.0);

        b->UH(i,j) += (ww*qxN(i-1,j) + we*qxN(i,j))/ws;

        b->WH(i,j) = WHp(i,j);

        if(nh_active(p,b,i,j)==1)
        b->WH(i,j) += dt*(cwn*b->press(i,j)/p->W1 + 0.25*WL0(i,j)*phiN(i,j));
        }

        else
        b->WH(i,j) = 0.0;
    }
}
