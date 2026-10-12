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

#include"nhflow_curv.h"
#include"sflow_curv.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"slice.h"
#include<cmath>
#include<iostream>
#include<vector>

void nhflow_curv::check(lexer *p)
{
    int bad = 0;

    auto refuse = [&](bool cond, const char *what)
    {
        if(cond)
        {
            if(p->mpirank==0)
            cout<<"!!! NHFLOW on a curvilinear grid (G 50): "<<what<<" is not available yet !!!"<<endl;
            ++bad;
        }
    };

    refuse(p->j_dir!=1, "a 2D vertical slice (j_dir 0)");
    refuse(p->A520==3, "the non-hydrostatic scheme A 520 3");
    refuse(p->A511==2, "the HLLC flux A 511 2 (A 511 1: HLL)");
    refuse(p->geo_mb>0, "a multi-block grid (DIVEMesh R 20 1)");
    refuse(p->G1>0 && p->P18>0, "the LAGOON output of the mesh refinement (P 18)");
    refuse(p->X10>0, "6DOF bodies (X 10)");
    refuse(p->X330>0, "membranes (X 330)");
    refuse(p->S10>0, "sediment transport (S 10)");
    refuse(p->Z30>0, "FEM structures (Z 30)");
    refuse(p->A580>0 || p->A599>0 || p->A581>0 || p->A583>0 || p->A584>0 || p->A585>0 || p->A586>0
           || p->A587>0 || p->A588>0 || p->A589>0 || p->A590>0, "immersed solids (A 580-599)");
    refuse(p->A570>0, "wind forcing (A 570)");
    refuse(p->A750>0, "the SEASTATE coupling (A 750)");
    refuse(p->I11>0, "the potential flow initialisation (I 11)");
    refuse(p->A519==2, "the side-wall friction A 519 2");
    refuse(p->B98==2 && p->B96_1>0.0 && p->B108==0, "the generation zone of B 96 (an x range of the box): give it by B 108 in physical coordinates");
    refuse((p->B99==1 || p->B99==2) && p->B96_2>0.0 && p->B107==0, "the beach zone of B 96 (an x range of the box): give it by B 107 in physical coordinates");
    refuse(p->B97==1, "the relaxation to the inflow current in the numerical beach (B 97 1)");
    refuse(p->B98>=3 && fabs(p->W10)>0.0, "Dirichlet or active wave generation (B 98 3, 4) on a current (W 10)");

    if(bad>0)
    exit(1);
}

void nhflow_curv::normal_velocities(lexer *p, fdm_nhf *d)
{
    double nx,ny,len;

    ULOOP
    {
    normal(p,1,i,j,nx,ny,len);
    const double us = d->Us[IJK]*nx + d->Vs[IJK]*ny;
    const double un = d->Un[IJK]*nx + d->Vn[IJK]*ny;
    d->Us[IJK] = us;
    d->Un[IJK] = un;
    }

    VLOOP
    {
    normal(p,2,i,j,nx,ny,len);
    const double ve = d->Ue[IJK]*nx + d->Ve[IJK]*ny;
    const double vw = d->Uw[IJK]*nx + d->Vw[IJK]*ny;
    d->Ve[IJK] = ve;
    d->Vw[IJK] = vw;
    }
}

double nhflow_curv::bedslope(lexer *p, fdm_nhf *d, int comp, int ii, int jj)
{
    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
    const int qw = q - p->jmax, qs = q - 1;
    const double *S1 = (comp==1) ? p->SX1 : p->SY1;
    const double *S2 = (comp==1) ? p->SX2 : p->SY2;

    return fabs(p->W22)*d->eta(ii,jj)*(d->dfx(ii,jj)*S1[q] - d->dfx(ii-1,jj)*S1[qw]
                                      + d->dfy(ii,jj)*S2[q] - d->dfy(ii,jj-1)*S2[qs])/p->AREA2D[q];
}

void nhflow_curv::grad5(lexer *p, int ii, int jj, double fE, double fW, double fN, double fS, double &fx, double &fy)
{
    sflow_curv::grad5(p,ii,jj,fE,fW,fN,fS,fx,fy);
}

void nhflow_curv::grad_idx(lexer *p, int ii, int jj, double fxi, double feta, double &fx, double &fy)
{
    // grad5 takes differences over two cells
    sflow_curv::grad5(p,ii,jj,2.0*fxi,0.0,2.0*feta,0.0,fx,fy);
}

void nhflow_curv::grad_lim(lexer *p, int ii, int jj, double fC, double fE, double fW, double fN, double fS, double &fx, double &fy)
{
    auto X = [&](int a, int b) { return p->XC2D[(a-p->imin)*p->jmax + (b-p->jmin)]; };
    auto Y = [&](int a, int b) { return p->YC2D[(a-p->imin)*p->jmax + (b-p->jmin)]; };
    auto dist = [&](int a, int b, int c, int e) { const double dx = X(c,e)-X(a,b), dy = Y(c,e)-Y(a,b); return sqrt(dx*dx + dy*dy); };
    auto limiter = [](double v1, double v2)
    {
        double denom = fabs(v1) + fabs(v2);
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        return (v1*fabs(v2) + fabs(v1)*v2)/denom;
    };
    
    const double gi = limiter((fE-fC)/dist(ii,jj,ii+1,jj), (fC-fW)/dist(ii-1,jj,ii,jj));
    const double gj = limiter((fN-fC)/dist(ii,jj,ii,jj+1), (fC-fS)/dist(ii,jj-1,ii,jj));
    
    // unit vectors along the lines
    double e1x = X(ii+1,jj)-X(ii-1,jj), e1y = Y(ii+1,jj)-Y(ii-1,jj);
    double e2x = X(ii,jj+1)-X(ii,jj-1), e2y = Y(ii,jj+1)-Y(ii,jj-1);
    const double l1 = sqrt(e1x*e1x + e1y*e1y), l2 = sqrt(e2x*e2x + e2y*e2y);
    e1x/=l1; e1y/=l1; e2x/=l2; e2y/=l2;
    
    const double det = e1x*e2y - e1y*e2x;
    fx = (gi*e2y - gj*e1y)/det;
    fy = (gj*e1x - gi*e2x)/det;
}

void nhflow_curv::grad_lim(lexer *p, slice &f, int ii, int jj, double &fx, double &fy)
{
    grad_lim(p,ii,jj,f(ii,jj),f(ii+1,jj),f(ii-1,jj),f(ii,jj+1),f(ii,jj-1),fx,fy);
}

double nhflow_curv::lap(lexer *p, slice &f, int ii, int jj)
{
    double cE,cW,cN,cS;
    sflow_curv::lapcoef(p,ii,jj,cE,cW,cN,cS);
    const double fc = f(ii,jj);

    double r = cE*(f(ii+1,jj)-fc) + cW*(f(ii-1,jj)-fc) + cN*(f(ii,jj+1)-fc) + cS*(f(ii,jj-1)-fc);

    r += sflow_curv::kflux(p,f,ii,jj,[](int,int){ return true; });

    return r;
}

double nhflow_curv::width1(lexer *p, int ii, int jj)
{
    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
    const double l = 0.5*(p->LEN1[q] + p->LEN1[q-p->jmax]);
    return (l>0.0) ? p->AREA2D[q]/l : 1.0e20;
}

double nhflow_curv::width2(lexer *p, int ii, int jj)
{
    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
    const double l = 0.5*(p->LEN2[q] + p->LEN2[q-1]);
    return (l>0.0) ? p->AREA2D[q]/l : 1.0e20;
}

void nhflow_curv::kvec(lexer *p, int dir, int ii, int jj, double &kx, double &ky)
{
    sflow_curv::kvec(p,dir,ii,jj,kx,ky);
}

void nhflow_curv::lapcoef(lexer *p, int ii, int jj, double &cIp, double &cIm, double &cJp, double &cJm)
{
    sflow_curv::lapcoef(p,ii,jj,cIp,cIm,cJp,cJm);
}

void nhflow_curv::pgrad_layer(lexer *p, const double *P, int ii, int jj, int kk, double &px, double &py)
{
    auto L = [&](int a, int b)
    {
        const int q = (a-p->imin)*p->jmax*p->kmaxF + (b-p->jmin)*p->kmaxF + kk-p->kmin;
        return 0.5*(P[q+1] + P[q]);
    };
    grad5(p,ii,jj,L(ii+1,jj),L(ii-1,jj),L(ii,jj+1),L(ii,jj-1),px,py);
}

double nhflow_curv::divh_level(lexer *p, fdm_nhf *d, const double *U, const double *V, int ii, int jj, int kk)
{
    auto C = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    auto xc = [&](int a, int b) { return p->XC2D[(a-p->imin)*p->jmax + (b-p->jmin)]; };
    auto yc = [&](int a, int b) { return p->YC2D[(a-p->imin)*p->jmax + (b-p->jmin)]; };

    double faci = p->DZN[kk-1+marge]/(p->DZN[kk+marge]+p->DZN[kk-1+marge]);
    double facj = faci;

    // next to the bed: the slope of the bed along the index directions (Cartesian: Bx, By)
    if(kk==0)
    {
        double ex = xc(ii+1,jj)-xc(ii-1,jj), ey = yc(ii+1,jj)-yc(ii-1,jj);
        double l = sqrt(ex*ex+ey*ey);
        const double bi = (l>0.0) ? fabs(d->Bx(ii,jj)*ex + d->By(ii,jj)*ey)/l : 0.0;
        ex = xc(ii,jj+1)-xc(ii,jj-1); ey = yc(ii,jj+1)-yc(ii,jj-1);
        l = sqrt(ex*ex+ey*ey);
        const double bj = (l>0.0) ? fabs(d->Bx(ii,jj)*ex + d->By(ii,jj)*ey)/l : 0.0;

        faci = MAX((1.0 - p->A522*bi),0.0)*p->DZN[kk-1+marge]/(p->DZN[kk+marge]+p->DZN[kk-1+marge]);
        facj = MAX((1.0 - p->A522*bj),0.0)*p->DZN[kk-1+marge]/(p->DZN[kk+marge]+p->DZN[kk-1+marge]);
    }

    auto li = [&](const double *f, int a, int b) { return (1.0-faci)*f[C(a,b,kk)] + faci*f[C(a,b,kk-1)]; };
    auto lj = [&](const double *f, int a, int b) { return (1.0-facj)*f[C(a,b,kk)] + facj*f[C(a,b,kk-1)]; };

    double ux,uy,vx,vy;
    grad5(p,ii,jj,li(U,ii+1,jj),li(U,ii-1,jj),lj(U,ii,jj+1),lj(U,ii,jj-1),ux,uy);
    grad5(p,ii,jj,li(V,ii+1,jj),li(V,ii-1,jj),lj(V,ii,jj+1),lj(V,ii,jj-1),vx,vy);

    return ux + vy;
}

double nhflow_curv::sigma_cross(lexer *p, const double *P, int ii, int jj, int kk)
{
    auto q = [&](int a, int b)
    {
        const int f = (a-p->imin)*p->jmax*p->kmaxF + (b-p->jmin)*p->kmaxF + kk-p->kmin;
        return (P[f+1] - P[f-1])/(p->DZN[kk+marge]+p->DZN[kk-1+marge]);
    };
    double qx,qy;
    grad5(p,ii,jj,q(ii+1,jj),q(ii-1,jj),q(ii,jj+1),q(ii,jj-1),qx,qy);

    const int f = (ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF + kk-p->kmin;
    return 2.0*(p->sigx[f]*qx + p->sigy[f]*qy);
}

void nhflow_curv::pair_bc(lexer *p, double *u, double *v)
{
    auto C = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    double nx,ny,len;
    
    // river outflow (B 60, start4V outflow 2): the normal velocity of the ghost cells advected out and
    // kept positive (outwards), the tangential one zero gradient
    if(p->B60>=1 && p->open_xp==0)
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    for(int kk=0; kk<p->knoz; ++kk)
    if(p->flag4[C(ii,jj,kk)]>0 && p->flag4[C(ii+1,jj,kk)]<0)
    {
        normal(p,1,ii,jj,nx,ny,len);
        const int a = C(ii,jj,kk), am = C(ii-1,jj,kk);
        const double un = u[a]*nx + v[a]*ny, unm = u[am]*nx + v[am]*ny;
        const double ut = -u[a]*ny + v[a]*nx;
        
        for(int q=1; q<=3; ++q)
        {
            const double ung = (q==1) ? MAX(0.0,un - p->dt/width1(p,ii,jj)*p->Uo*(un-unm)) : MAX(0.0,un);
            const int g = C(ii+q,jj,kk);
            u[g] = ung*nx - ut*ny;
            v[g] = ung*ny + ut*nx;
        }
    }
    

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    for(int kk=0; kk<p->knoz; ++kk)
    if(p->flag4[C(ii,jj,kk)]>0)
    {
        // wall on the -j side
        if(p->flag4[C(ii,jj-1,kk)]<0 && (p->open_ym==0 || p->origin_j+jj>0))
        {
            normal(p,2,ii,jj-1,nx,ny,len);
            for(int q=0; q<3; ++q)
            {
                const int a = C(ii,jj+q,kk), g = C(ii,jj-1-q,kk);
                const double un = u[a]*nx + v[a]*ny;
                u[g] = u[a] - 2.0*un*nx;
                v[g] = v[a] - 2.0*un*ny;
            }
        }

        // wall on the +j side
        if(p->flag4[C(ii,jj+1,kk)]<0 && (p->open_yp==0 || p->origin_j+jj<p->gknoy-1))
        {
            normal(p,2,ii,jj,nx,ny,len);
            for(int q=0; q<3; ++q)
            {
                const int a = C(ii,jj-q,kk), g = C(ii,jj+1+q,kk);
                const double un = u[a]*nx + v[a]*ny;
                u[g] = u[a] - 2.0*un*nx;
                v[g] = v[a] - 2.0*un*ny;
            }
        }
    }
}

void nhflow_curv::bc_fluxes(lexer *p, fdm_nhf *d, int comp, double *Fx, double *Fy)
{
    // as ghostcell::start1V / start2V, the fluxes of the Cartesian components through the face normal:
    //   walls:              P n (momentum), 0 (continuity, left from start1V / start2V)
    //   in- / outflow, open edges: the ghost cell's state, q (u.n) + P n (momentum), UH.n (continuity)
    auto C = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    const double g = fabs(p->W22);
    auto P = [&](double eta, double df) { return 0.5*g*eta*eta + g*eta*df; };
    
    // flux of the ghost cell q through the face (dir, a, b): wall 0, open 1
    auto flux = [&](int dir, int fa, int fb, int q, double eta, double df, int open)
    {
        double nx,ny,len;
        normal(p,dir,fa,fb,nx,ny,len);
        
        if(comp==4)
        return open ? d->UH[q]*nx + d->VH[q]*ny : 0.0;
        
        const double nc = (comp==1) ? nx : ((comp==2) ? ny : 0.0);
        
        if(open==0)
        return P(eta,df)*nc;
        
        // summed in the order of start1V (identical to it on Cartesian nodes)
        const double un = d->U[q]*nx + d->V[q]*ny;
        const double *Q = (comp==1) ? d->UH : ((comp==2) ? d->VH : d->WH);
        return Q[q]*un + 0.5*g*eta*eta*nc + g*eta*df*nc;
    };
    
    int inflow = 0, outflow = 0;
    
    if(p->B60>=1 || p->B98>=3 || p->open_xm==1)
    inflow = 1;
    
    if(p->B99>=3 || p->B60>=1 || p->open_xp==1)
    outflow = 1;
    
    // i-faces
    for(int ii=0; ii<p->knox-p->ulast; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    for(int kk=0; kk<p->knoz; ++kk)
    if(p->flag1[C(ii,jj,kk)]>0)
    {
        if(p->flag1[C(ii-1,jj,kk)]<0 && (comp<4 || inflow==1))
        Fx[C(ii-1,jj,kk)] = flux(1,ii-1,jj,C(ii-1,jj,kk),d->eta(ii-1,jj),d->dfx(ii,jj),inflow);
        
        if(p->flag1[C(ii+1,jj,kk)]<0 && (comp<4 || outflow==1))
        Fx[C(ii+1,jj,kk)] = flux(1,ii+1,jj,C(ii+2,jj,kk),d->eta(ii+2,jj),d->dfx(ii+1,jj),outflow);
    }
    
    // j-faces (open: the y edges of iowave)
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy-p->vlast; ++jj)
    for(int kk=0; kk<p->knoz; ++kk)
    if(p->flag2[C(ii,jj,kk)]>0)
    {
        const int openm = (p->open_ym==1 && p->origin_j+jj==0) ? 1 : 0;
        const int openp = (p->open_yp==1 && p->origin_j+jj+1==p->gknoy-1) ? 1 : 0;
        
        if(p->flag2[C(ii,jj-1,kk)]<0 && (comp<4 || openm==1))
        Fy[C(ii,jj-1,kk)] = flux(2,ii,jj-1,C(ii,jj-1,kk),d->eta(ii,jj-1),d->dfy(ii,jj),openm);
        
        if(p->flag2[C(ii,jj+1,kk)]<0 && (comp<4 || openp==1))
        Fy[C(ii,jj+1,kk)] = flux(2,ii,jj+1,C(ii,jj+2,kk),d->eta(ii,jj+2),d->dfy(ii,jj+1),openp);
    }
}

void nhflow_curv::grad_dsig(lexer *p, const double *f, int ii, int jj, int kk, double &qx, double &qy)
{
    auto C = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    const double ds = p->DZP[kk+marge] + p->DZP[kk-1+marge];
    auto q = [&](int a, int b) { return (f[C(a,b,kk+1)] - f[C(a,b,kk-1)])/ds; };
    grad5(p,ii,jj,q(ii+1,jj),q(ii-1,jj),q(ii,jj+1),q(ii,jj-1),qx,qy);
}

// horizontal divergence of (UH, VH) in the cells and one ring around them (diffusion of u, v)
static std::vector<double> curv_diff_D;

void nhflow_curv::diff_div(lexer *p, const double *UH, const double *VH)
{
    auto C = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    std::vector<double> &D = curv_diff_D;
    D.assign(size_t(p->imax)*size_t(p->jmax)*size_t(p->kmax),0.0);
    
    for(int ii=-1; ii<=p->knox; ++ii)
    for(int jj=-1; jj<=p->knoy; ++jj)
    for(int kk=0; kk<p->knoz; ++kk)
    {
        double ux,uy,vx,vy;
        grad5(p,ii,jj,UH[C(ii+1,jj,kk)],UH[C(ii-1,jj,kk)],UH[C(ii,jj+1,kk)],UH[C(ii,jj-1,kk)],ux,uy);
        grad5(p,ii,jj,VH[C(ii+1,jj,kk)],VH[C(ii-1,jj,kk)],VH[C(ii,jj+1,kk)],VH[C(ii,jj-1,kk)],vx,vy);
        D[C(ii,jj,kk)] = ux + vy;
    }
}

double nhflow_curv::kflux3(lexer *p, const double *f, int ii, int jj, int kk)
{
    auto C = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    auto lay = [&](int a, int b) { return f[C(a,b,kk)]; };
    auto g = [&](int a, int b, double &gx, double &gy) { grad5(p,a,b,lay(a+1,b),lay(a-1,b),lay(a,b+1),lay(a,b-1),gx,gy); };
    auto ok = [&](int a, int b) { const int n = C(a,b,kk); return p->flag4[n]>0 && p->DF[n]>0; };
    double gx0,gy0,gx,gy;
    g(ii,jj,gx0,gy0);
    const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
    double sk = 0.0;
    
    for(int dd=0; dd<4; ++dd)
    {
        const int a = ii+di[dd], b = jj+dj[dd];
        if(!ok(a,b))
        continue;
        const int dir = (dd<2) ? 1 : 2;
        const int fi = (dd==1) ? ii-1 : ii, fj = (dd==3) ? jj-1 : jj;
        double kx,ky;
        kvec(p,dir,fi,fj,kx,ky);
        g(a,b,gx,gy);
        const double sg = (dd==0 || dd==2) ? 1.0 : -1.0;
        sk += sg*0.5*((gx0+gx)*kx + (gy0+gy)*ky);
    }
    
    return sk/p->AREA2D[(ii-p->imin)*p->jmax + (jj-p->jmin)];
}

double nhflow_curv::diff_explicit(lexer *p, int comp, int ii, int jj, int kk, const double *UH, const double *VH, const double *WH)
{
    auto C = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    const double *f = (comp==1) ? UH : ((comp==2) ? VH : WH);
    const int q2 = (ii-p->imin)*p->jmax + (jj-p->jmin);
    const int fl = (ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF + kk-p->kmin, flp = fl+1;
    double r = 0.0;
    
    // non-orthogonal part of the horizontal Laplacian (faces to fluid cells)
    r += kflux3(p,f,ii,jj,kk);
    
    // d/dx, d/dy of the horizontal divergence (the transpose and the second normal stress)
    double qx,qy;
    if(comp<3)
    {
    const std::vector<double> &D = curv_diff_D;
    double Dx,Dy;
    grad5(p,ii,jj,D[C(ii+1,jj,kk)],D[C(ii-1,jj,kk)],D[C(ii,jj+1,kk)],D[C(ii,jj-1,kk)],Dx,Dy);
    r += (comp==1) ? Dx : Dy;
    
    // vertical transpose stress (leading order, as the Cartesian kernels)
    grad_dsig(p,WH,ii,jj,kk,qx,qy);
    r += ((comp==1) ? qx : qy)*p->sigz[q2];
    }
    else
    {
    double ux,uy,vx,vy;
    grad_dsig(p,UH,ii,jj,kk,ux,uy);
    grad_dsig(p,VH,ii,jj,kk,vx,vy);
    r += (ux + vy)*p->sigz[q2];
    }
    
    // sigma cross terms 2 (sigx d/dx + sigy d/dy) df/dsigma
    grad_dsig(p,f,ii,jj,kk,qx,qy);
    r += 2.0*0.5*(p->sigx[fl]+p->sigx[flp])*qx + 2.0*0.5*(p->sigy[fl]+p->sigy[flp])*qy;
    
    return r;
}

void nhflow_curv::ediff(lexer *p, fdm_nhf *d, int comp, double *F, double *UH, double *VH, double *WH)
{
    auto C = [&](int a, int b, int c) { return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin; };
    const double *f = (comp==1) ? UH : ((comp==2) ? VH : WH);
    
    if(comp<3)
    diff_div(p,UH,VH);
    
    // neighbour values as nhflow_ediff (A 513 1: zero gradient towards solid or ghost cells)
    auto nb = [&](int c, int n) { return (p->A513==1 && (p->flag4[n]<0 || p->DF[n]<0)) ? f[c] : f[n]; };
    
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    for(int kk=0; kk<p->knoz; ++kk)
    {
        const int c = C(ii,jj,kk);
        const int q2 = (ii-p->imin)*p->jmax + (jj-p->jmin);
        
        if(p->flag4[c]<0 || p->wet[q2]==0 || p->DF[c]<0)
        continue;
        
        const double visc = d->VISC[c] + d->EV[c];
        const int fl = (ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF + kk-p->kmin, flp = fl+1;
        const double sz2 = (comp==3 ? 2.0 : 1.0)*pow(p->sigz[q2],2.0);
        const double s2b = pow(p->sigx[fl],2.0) + pow(p->sigy[fl],2.0) + sz2;
        const double s2t = pow(p->sigx[flp],2.0) + pow(p->sigy[flp],2.0) + sz2;
        
        const double fc = f[c];
        const double fkm = nb(c,C(ii,jj,kk-1)), fkp = nb(c,C(ii,jj,kk+1));
        
        // compact horizontal Laplacian
        double cE,cW,cN,cS;
        lapcoef(p,ii,jj,cE,cW,cN,cS);
        const double lh = cE*(nb(c,C(ii+1,jj,kk))-fc) + cW*(nb(c,C(ii-1,jj,kk))-fc)
                        + cN*(nb(c,C(ii,jj+1,kk))-fc) + cS*(nb(c,C(ii,jj-1,kk))-fc);
        
        F[c] += visc*lh
              + visc*(s2t*(fkp-fc)/p->DZP[kk+marge] - s2b*(fc-fkm)/p->DZP[kk-1+marge])/p->DZN[kk+marge]
              + visc*p->sigxx[fl]*(fkp-fkm)/(p->DZP[kk+marge]+p->DZP[kk-1+marge])
              + visc*diff_explicit(p,comp,ii,jj,kk,UH,VH,WH);
    }
}
