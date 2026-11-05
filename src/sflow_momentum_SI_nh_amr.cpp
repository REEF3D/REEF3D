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
// SFLOW semi-implicit step with the coupled non-hydrostatic pressure (A 226 1, A 220 2/3) on the grids of the mesh
// refinement (G 1): the parts of the coupled system of one grid. sflow_amr::si_nh_solve (sflow_amr_si.cpp) runs the
// composite solve over the leaf cells of all grids with them:
//  nh_amr_cells    the cell terms of the constraint (2 dt cw/H, the affine part), h^3 of the A 220 3 term, wet pattern
//  nh_amr_guess    the initial guess (eta and P extrapolated in time) into nE, nP
//  nh_amr_lap      nL = h^3 lap(nE) on all cells of the grid (the cells around the grid hold the filled nE)
//  nh_amr_faces    the face discharges of (nE, nP, nL): nF1 with the constants of the mass row, nF2 with all
//                  constants (the constraint); the composite solve replaces those of the coarse faces next to a finer
//                  patch by the mean of the fine faces (the divergence of the discharges the face update takes)
//  nh_amr_row      the mass row and the constraint of a cell from the face discharges
//  nh_amr_leaf     the leaf cells of the solve as slice indices, nh_amr_rows: nh_amr_row of all of them
//  nh_amr_prec_*   the preconditioner of nh_solve on the own cells of the grid (eta eliminated locally, one V-cycle
//                  for psi, A 220 3 a second one), the cells around the grid 0 (Dirichlet)
//  nh_amr_finish   eta^{n+1}, P^{n+1} and the non-hydrostatic face discharges qxN, qyN of the solution
//  nh_amr_wetdry   the wet pattern of the mass rows (A 229 1) on the leaf cells
// The arithmetic of a row is that of nh_solve (one grid).
// ---------------------------------------------------------------------------

#include"sflow_momentum_SI.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"reefmg2D.h"
#include"matrix2D.h"
#include"vec2D.h"
#include<vector>
#include<memory>

#include"sflow_momentum_SI_nh.h"

void sflow_momentum_SI::nh_amr_cells(lexer *p, fdm2D *b)
{
    const double dt = p->dt;
    const bool yd = (p->j_dir==1);
    auto in = [&](int a, int c) { return p->flagslice4[(a-p->imin)*p->jmax + (c-p->jmin)]>0; };

    // the coefficients of lap(nE) of the cells (nh_amr_lap), fixed for the solve
    nh_lc.assign((size_t)4*p->imax*p->jmax,0.0);

    SLICELOOP4
    {
    const double H = WLq(i,j);
    const bool a = nA(i,j)>0.5 && Bdisp>0.0;

    nH3(i,j) = a ? WLh(i,j)*WLh(i,j)*WLh(i,j) : 0.0;
    nW2(i,j) = 2.0*dt*cwn/H;

        if(a)
        {
        double *c = &nh_lc[4*(size_t)(IJ)];
        c[0] = in(i+1,j) ? 1.0/(p->DXP[IP]*p->DXN[IP]) : 0.0;
        c[1] = in(i-1,j) ? 1.0/(p->DXP[IM1]*p->DXN[IP]) : 0.0;
        c[2] = (yd && in(i,j+1)) ? 1.0/(p->DYP[JP]*p->DYN[JP]) : 0.0;
        c[3] = (yd && in(i,j-1)) ? 1.0/(p->DYP[JM1]*p->DYN[JP]) : 0.0;
        }
    nCs(i,j) = 2.0*(WHp(i,j) + dt*0.25*H*phiN(i,j))/H + nsx(i,j)*nUs(i,j) + (yd ? nsy(i,j)*nVs(i,j) : 0.0);
    nPW(i,j) = 1.0;
    }
}

void sflow_momentum_SI::nh_amr_guess(lexer *p, fdm2D *b)
{
    // eta and P extrapolated from t^n and t^{n-1} (nh_solve)
    const double r = p->dt/MAX(p->dt_old,1.0e-20);

    SLICELOOP4
    {
    const double P0 = b->press(i,j)/p->W1;
    nE(i,j) = nh_x0 ? (1.0+r)*eta0(i,j) - r*nEo(i,j) : eta0(i,j);
    nP(i,j) = nA(i,j)>0.5 ? (nh_x0 ? (1.0+r)*P0 - r*nPo(i,j) : P0) : 0.0;
    nEo(i,j) = eta0(i,j);
    nPo(i,j) = P0;
    }

    nh_x0 = true;
}

void sflow_momentum_SI::nh_amr_lap(lexer *p, ghostcell *pgc_nh)
{
    const bool yd = (p->j_dir==1);
    const int J = p->jmax;
    const double *H3 = nH3.V, *E = nE.V, *lc = nh_lc.data();
    double *Lv = nL.V;

    SLICELOOP4
    {
    const int q = IJ;
    double v = 0.0;

        if(H3[q]>0.0)
        {
        const double *c = lc + 4*(size_t)q;
        double lap = c[0]*(E[q+J] - E[q]) - c[1]*(E[q] - E[q-J]);

            if(yd)
            lap += c[2]*(E[q+1] - E[q]) - c[3]*(E[q] - E[q-1]);

        v = H3[q]*lap;
        }

    Lv[q] = v;
    }

    // the partition halo of level 0 (the faces of the halo take the B term)
    exchange4(p,pgc_nh,nL);
}

void sflow_momentum_SI::nh_amr_faces(lexer *p, bool affine)
{
    // q = EL E_L + ER E_R + PL P_L + PR P_R + B (L_R - L_L) + k1 / k2 (kface of nh_assemble)
    const bool bd = (Bdisp>0.0);
    const int J = p->jmax;
    const double *E = nE.V, *Pp = nP.V, *Lp = nL.V;

    // all faces of the arrays (the faces of the partition halo with the exchanged coefficients and values, as
    // nh_solve computes them on the fly; the coefficients of faces nh_assemble does not set are 0)
    {
    const double *EL = kxEL.V, *ER = kxER.V, *PL = kxPL.V, *PR = kxPR.V, *B = nxB.V, *K1 = kx1.V, *K2 = kx2.V;
    double *F1 = nF1x.V, *F2 = nF2x.V;
    const int nq = (p->imax-1)*J;

    for(int q=0; q<nq; ++q)
    {
    const double v = EL[q]*E[q] + ER[q]*E[q+J] + PL[q]*Pp[q] + PR[q]*Pp[q+J]
                   + (bd ? B[q]*(Lp[q+J] - Lp[q]) : 0.0);
    F1[q] = v + (affine ? K1[q] : 0.0);
    F2[q] = v + (affine ? K2[q] : 0.0);
    }
    }

    if(p->j_dir==1)
    {
    const double *EL = kyEL.V, *ER = kyER.V, *PL = kyPL.V, *PR = kyPR.V, *B = nyB.V, *K1 = ky1.V, *K2 = ky2.V;
    double *F1 = nF1y.V, *F2 = nF2y.V;

    for(int a=0; a<p->imax; ++a)
    for(int q=a*J; q<a*J+J-1; ++q)
    {
    const double v = EL[q]*E[q] + ER[q]*E[q+1] + PL[q]*Pp[q] + PR[q]*Pp[q+1]
                   + (bd ? B[q]*(Lp[q+1] - Lp[q]) : 0.0);
    F1[q] = v + (affine ? K1[q] : 0.0);
    F2[q] = v + (affine ? K2[q] : 0.0);
    }
    }
}

void sflow_momentum_SI::nh_amr_leaf(lexer *p, const std::vector<int> &leaf)
{
    // the leaf cells (pairs i, j) of the solve as slice indices with their widths, for nh_amr_rows
    nh_lq.clear(); nh_ldx.clear(); nh_ldy.clear();

    for(size_t k=0; k<leaf.size(); k+=2)
    {
    const int ii = leaf[k], jj = leaf[k+1];
    nh_lq.push_back((ii-p->imin)*p->jmax + (jj-p->jmin));
    nh_ldx.push_back(p->DXN[ii+marge]);
    nh_ldy.push_back(p->j_dir==1 ? p->DYN[jj+marge] : 1.0);
    }
}

void sflow_momentum_SI::nh_amr_rows(lexer *p, bool affine, double *y)
{
    // the mass row y[2k] and the constraint y[2k+1] of the leaf cells of nh_amr_leaf (the rows of nh_amr_row)
    const double tdt = theta*p->dt;
    const bool yd = (p->j_dir==1);
    const int J = p->jmax, nl = (int)nh_lq.size();
    const double *F1x = nF1x.V, *F2x = nF2x.V, *F1y = nF1y.V, *F2y = nF2y.V;
    const double *PW = nPW.V, *E = nE.V, *Pp = nP.V, *A = nA.V, *W2 = nW2.V, *Cs = nCs.V;

    for(int k=0; k<nl; ++k)
    {
    const int q = nh_lq[k];

    double d1 = (F1x[q] - F1x[q-J])/nh_ldx[k];
    double d2 = (F2x[q] - F2x[q-J])/nh_ldx[k];

        if(yd)
        {
        d1 += (F1y[q] - F1y[q-1])/nh_ldy[k];
        d2 += (F2y[q] - F2y[q-1])/nh_ldy[k];
        }

    y[2*k] = PW[q]*E[q] + tdt*d1;

    if(A[q]<0.5)
    y[2*k+1] = Pp[q];

    else
    y[2*k+1] = d2 + W2[q]*Pp[q] + (affine ? Cs[q] : 0.0);
    }
}

void sflow_momentum_SI::nh_amr_row(lexer *p, fdm2D *b, int ii, int jj, bool affine, double &ym, double &yc)
{
    const double tdt = theta*p->dt;
    const int ip = ii+marge, jp = jj+marge;

    double d1 = (nF1x(ii,jj) - nF1x(ii-1,jj))/p->DXN[ip];
    double d2 = (nF2x(ii,jj) - nF2x(ii-1,jj))/p->DXN[ip];

    if(p->j_dir==1)
    {
    d1 += (nF1y(ii,jj) - nF1y(ii,jj-1))/p->DYN[jp];
    d2 += (nF2y(ii,jj) - nF2y(ii,jj-1))/p->DYN[jp];
    }

    ym = nPW(ii,jj)*nE(ii,jj) + tdt*d1;

    if(nA(ii,jj)<0.5)
    yc = nP(ii,jj);

    else
    yc = d2 + nW2(ii,jj)*nP(ii,jj) + (affine ? nCs(ii,jj) : 0.0);
}

void sflow_momentum_SI::nh_amr_prec_setup(lexer *p, ghostcell *pgc, const std::function<bool(int,int)> &own)
{
    // the matrices of the preconditioner of nh_solve; the rows of the cells the grid does not hold: identity
    const double dt = p->dt;
    const double tdt = theta*dt;
    const bool yd = (p->j_dir==1);
    auto id = [&](int a, int c) { return (a-p->imin)*p->jmax + (c-p->jmin); };
    auto in = [&](int a, int c) { return p->flagslice4[id(a,c)]>0; };
    auto ac = [&](int a, int c) { return in(a,c) && nA(a,c)>0.5; };
    auto &S = *nhlu;

    // the own cells of the solve and the cells in SLICELOOP4 order with their widths, for nh_amr_prec
    nh_own.assign((size_t)p->imax*p->jmax,0);
    nh_cq.clear(); nh_dxn.clear(); nh_dyn.clear();
    SLICELOOP4
    {
    nh_own[IJ] = own(i,j) ? 1 : 0;
    nh_cq.push_back(IJ);
    nh_dxn.push_back(p->DXN[IP]);
    nh_dyn.push_back(yd ? p->DYN[JP] : 1.0);
    }

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

    if(!S.mg1)
    {
    S.mg1.reset(patch ? new reefmg2D(p,pgc,true) : new reefmg2D(p,pgc));
    S.mg1->set_keep_matrix(true);
    S.M1.reset(new matrix2D(p));
    S.r1.reset(new vec2D(p));

        if(!yd)
        S.mg1->set_exact_line(true);
    }

    if(Bdisp>0.0 && !S.mg2)
    {
    S.mg2.reset(patch ? new reefmg2D(p,pgc,true) : new reefmg2D(p,pgc));
    S.mg2->set_keep_matrix(true);
    S.mg2->set_coarse_sweeps(2);
    S.M2.reset(new matrix2D(p));
    S.r2.reset(new vec2D(p));

        if(!yd)
        S.mg2->set_exact_line(true);
    }

    nh_rb1 = nh_rb2 = true;

    matrix2D &M = *S.M1;

    auto face_pre = [&](int dir, int a, int c, int a1, int c1, bool self_left, double dn, double &diag, double &off)
    {
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

        const int k = dir==1 ? S.bfx[id(a,c)] : S.bfy[id(a,c)];
        if(k<0)
        return;

        const bc_face &f = bcf[k];

        if(f.type==bc_level && Hf>0.0)
        diag = (af/dn)*nAE(f.i,f.j);

        if(f.type==bc_riemann || f.type==bc_flather)
        diag = tdt*f.B/dn*nAE(f.i,f.j);
    };

    {
    int n=0;
    SLICELOOP4
    {
    M.p[n] = 1.0;
    M.n[n] = M.s[n] = M.w[n] = M.e[n] = 0.0;

        if(own(i,j))
        {
        double d, o;
        M.p[n] = nPW(i,j)*nAE(i,j);

        face_pre(1,i,j,i+1,j,true,p->DXN[IP],d,o);    M.p[n] += d; M.n[n] = o;
        face_pre(1,i-1,j,i,j,false,p->DXN[IP],d,o);   M.p[n] += d; M.s[n] = o;

            if(yd)
            {
            face_pre(2,i,j,i,j+1,true,p->DYN[JP],d,o);  M.p[n] += d; M.w[n] = o;
            face_pre(2,i,j-1,i,j,false,p->DYN[JP],d,o); M.p[n] += d; M.e[n] = o;
            }
        }
    ++n;
    }
    }

    if(Bdisp>0.0)
    {
        matrix2D &M2 = *S.M2;

        // h^3 of the cells, computed once per cell
        S.h3c.assign((size_t)p->imax*p->jmax,-1.0);
        auto h3of = [&](int a, int c)
        {
            double &v = S.h3c[id(a,c)];
            if(v<0.0)
            v = pow(WLh(a,c),3.0);
            return v;
        };

        auto rho = [&](int dir, int a, int c, int a1, int c1)
        {
            const double Hf = dir==1 ? Hx(a,c) : Hy(a,c);

            if(!(ac(a,c) && ac(a1,c1)) || Hf<=0.0)
            return 0.0;

            const double af = dir==1 ? ax(a,c) : ay(a,c);
            const double cB = dir==1 ? nxB(a,c) : nyB(a,c);
            const double cc = 0.5*(dir==1 ? nxR(a,c) + nxL(a,c) : nyR(a,c) + nyL(a,c));
            const double aE = 0.5*(nAE(a,c) + nAE(a1,c1));
            const double h3 = 0.5*(h3of(a,c) + h3of(a1,c1));
            const double K = af*aE + tdt*cc;

            return K>0.0 ? tdt*cB*h3*aE/K : 0.0;
        };

        int n=0;
        SLICELOOP4
        {
        M2.p[n] = 1.0;
        M2.n[n] = M2.s[n] = M2.w[n] = M2.e[n] = 0.0;

            if(own(i,j))
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
            }
        ++n;
        }
    }
}

void sflow_momentum_SI::nh_amr_prec(lexer *p, ghostcell *pgc, const std::function<bool(int,int)> &own)
{
    // z = M^-1 r on the own cells (r: nRm, nRc; z: nZE, nZP), the preconditioner of nh_solve; the own cells and the
    // cell list of nh_amr_prec_setup
    const double tdt = theta*p->dt;
    const bool yd = (p->j_dir==1);
    auto &S = *nhlu;
    const char *ow = nh_own.data();
    const int N = (int)nh_cq.size(), J = p->jmax;
    const int *cq = nh_cq.data();
    const double *Rm = nRm.V, *Rc = nRc.V, *A = nA.V, *AE = nAE.V, *AP = nAP.V;
    double *BE = nBE.V, *BP = nBP.V, *Z = nZ.V, *ZE = nZE.V, *ZP = nZP.V;

    for(int n=0; n<N; ++n)
    {
    const int q = cq[n];
    const double rm = Rm[q], rc = Rc[q];
    BE[q] = 0.0;
    BP[q] = 0.0;

        if(ow[q])
        {
            if(A[q]>0.5)
            BE[q] = rm - tdt*rc;

            else
            BP[q] = rc;
        }
    }

    exchange4(p,pgc,nBE);
    exchange4(p,pgc,nBP);

    // the mass row of (bE, bP) without the A 220 3 term
    const double *XEL = kxEL.V, *XER = kxER.V, *XPL = kxPL.V, *XPR = kxPR.V;
    const double *YEL = kyEL.V, *YER = kyER.V, *YPL = kyPL.V, *YPR = kyPR.V;
    auto qx = [&](int q) { return XEL[q]*BE[q] + XER[q]*BE[q+J] + XPL[q]*BP[q] + XPR[q]*BP[q+J]; };
    auto qy = [&](int q) { return YEL[q]*BE[q] + YER[q]*BE[q+1] + YPL[q]*BP[q] + YPR[q]*BP[q+1]; };

    double *r1 = S.r1->V.data();
    for(int n=0; n<N; ++n)
    {
    const int q = cq[n];
    double r = 0.0;

        if(ow[q])
        {
        double dv = (qx(q) - qx(q-J))/nh_dxn[n];
        if(yd) dv += (qy(q) - qy(q-1))/nh_dyn[n];

        r = Rm[q] - (BE[q] + tdt*dv);
        }

    r1[n] = r;
    Z[q] = 0.0;
    }

    S.mg1->vcycle(p,pgc,nZ,*S.M1,*S.r1,nh_rb1);
    nh_rb1 = false;

    if(Bdisp>0.0)
    {
    double *r2 = S.r2->V.data();
    for(int n=0; n<N; ++n)
    r2[n] = ow[cq[n]] ? Z[cq[n]] : 0.0;

    S.mg2->vcycle(p,pgc,nZ,*S.M2,*S.r2,nh_rb2);
    nh_rb2 = false;
    }

    for(int n=0; n<N; ++n)
    {
    const int q = cq[n];
    const bool o = ow[q];
    ZE[q] = o ? BE[q] + AE[q]*Z[q] : 0.0;
    ZP[q] = o ? BP[q] + AP[q]*Z[q] : 0.0;
    }
}

void sflow_momentum_SI::nh_amr_finish(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // eta^{n+1}, P^{n+1} (nE, nP and nL of the solution, the cells around the grid filled) and the face terms
    auto in = [&](int a, int c) { return p->flagslice4[(a-p->imin)*p->jmax + (c-p->jmin)]>0; };

    SLICELOOP4
    {
    etaI(i,j) = nE(i,j);
    b->press(i,j) = nA(i,j)>0.5 ? p->W1*nP(i,j) : 0.0;
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

    if(p->j_dir==1)
    SLICELOOP2
    qyN(i,j) = nhq(2,i,j);

    exchange1(p,pgc,qxN);
    exchange2(p,pgc,qyN);
    exchange4(p,pgc,etaI);
    pgc->gcsl_start4(p,b->press,40);
}

int sflow_momentum_SI::nh_amr_wetdry(lexer *p, fdm2D *b, const std::function<bool(int,int)> &leaf)
{
    // new wet pattern of the leaf cells (nh_solve): dry where eta + d <= 0 and the cell has an open face
    const bool yd = (p->j_dir==1);
    int change = 0;

    SLICELOOP4
    if(leaf(i,j))
    {
        const bool coupled = ax(i,j)>0.0 || ax(i-1,j)>0.0 || bdiag(i,j)>0.0 || (yd && (ay(i,j)>0.0 || ay(i,j-1)>0.0));
        const double pn = (nE(i,j) + b->depth(i,j)>0.0 || !coupled) ? 1.0 : 0.0;

        if(pn!=nPW(i,j))
        {
        nPW(i,j) = pn;
        ++change;
        }
    }

    return change;
}
