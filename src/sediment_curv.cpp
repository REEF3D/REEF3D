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

#include"sediment_curv.h"
#include"sediment_fdm.h"
#include"ghostcell.h"
#include"slice.h"
#include<iostream>

void sediment_curv::check(lexer *p)
{
    if(p->S10==0)
    return;

    int bad = 0;

    auto refuse = [&](bool cond, const char *what)
    {
        if(cond)
        {
            if(p->mpirank==0)
            cout<<"!!! sediment transport on a curvilinear grid (G 50): "<<what<<" is not available yet !!!"<<endl;
            ++bad;
        }
    };

    refuse(p->S31==3, "the advective Exner form S 31 3 (not conservative; S 31 1 or 2)");
    refuse(p->S32>=3, "the Exner discretisations S 32 3-5 (S 32 1 FOU or 2 CDS)");
    refuse(p->S10==2, "the VRANS sediment bed (S 10 2)");

    if(bad>0)
    exit(1);
}

void sediment_curv::face_normals(lexer *p, ghostcell *pgc, sediment_fdm *s, int wetcheck)
{
    double nx,ny,len;

    SLICELOOP1
    {
    normal(p,1,i,j,nx,ny,len);
    s->P(i,j) = 0.5*((s->ub(i,j)+s->ub(i+1,j))*nx + (s->vb(i,j)+s->vb(i+1,j))*ny);

        if(wetcheck==1 && (p->wet[IJ]==0 || p->wet[Ip1J]==0))
        s->P(i,j) = 0.0;
    }

    SLICELOOP2
    {
    normal(p,2,i,j,nx,ny,len);
    s->Q(i,j) = 0.5*((s->ub(i,j)+s->ub(i,j+1))*nx + (s->vb(i,j)+s->vb(i,j+1))*ny);

        if(wetcheck==1 && (p->wet[IJ]==0 || p->wet[IJp1]==0))
        s->Q(i,j) = 0.0;
    }

    pgc->gcsl_start1(p,s->P,10);
    pgc->gcsl_start2(p,s->Q,11);
}

void sediment_curv::facevec(lexer *p, sediment_fdm *s, int dir, int ii, int jj, double &u, double &v)
{
    double nx,ny,len;
    normal(p,dir,ii,jj,nx,ny,len);

    const int i2 = ii + (dir==1?1:0);
    const int j2 = jj + (dir==2?1:0);

    u = 0.5*(s->ub(ii,jj) + s->ub(i2,j2));
    v = 0.5*(s->vb(ii,jj) + s->vb(i2,j2));

    // the normal component of the face (boundary conditions of P, Q)
    const double un = (dir==1) ? s->P(ii,jj) : s->Q(ii,jj);
    const double dn = un - (u*nx + v*ny);

    u += dn*nx;
    v += dn*ny;
}

void sediment_curv::cellvel(lexer *p, sediment_fdm *s, int ii, int jj, double &u, double &v)
{
    double pu1,pv1,pu2,pv2,qu1,qv1,qu2,qv2;

    facevec(p,s,1,ii-1,jj,pu1,pv1);
    facevec(p,s,1,ii,jj,pu2,pv2);
    facevec(p,s,2,ii,jj-1,qu1,qv1);
    facevec(p,s,2,ii,jj,qu2,qv2);

    const double pu = 0.5*(pu1+pu2), pv = 0.5*(pv1+pv2);
    const double qu = 0.5*(qu1+qu2), qv = 0.5*(qv1+qv2);

    // grid vectors of the cell: between the midpoints of its opposite faces
    auto xn = [&](int a, int b) { return p->XN2D[p->nij(a,b)]; };
    auto yn = [&](int a, int b) { return p->YN2D[p->nij(a,b)]; };

    const double g1x = 0.5*(xn(ii+1,jj)+xn(ii+1,jj+1)) - 0.5*(xn(ii,jj)+xn(ii,jj+1));
    const double g1y = 0.5*(yn(ii+1,jj)+yn(ii+1,jj+1)) - 0.5*(yn(ii,jj)+yn(ii,jj+1));
    const double g2x = 0.5*(xn(ii,jj+1)+xn(ii+1,jj+1)) - 0.5*(xn(ii,jj)+xn(ii+1,jj));
    const double g2y = 0.5*(yn(ii,jj+1)+yn(ii+1,jj+1)) - 0.5*(yn(ii,jj)+yn(ii+1,jj));

    const double det = g1x*g2y - g1y*g2x;

    if(fabs(det)<1.0e-30)
    {
    u = 0.5*(pu+qu);
    v = 0.5*(pv+qv);
    return;
    }

    // contravariant components: a1 = (pu,pv).g^1, a2 = (qu,qv).g^2, g^1 = (g2y,-g2x)/det, g^2 = (-g1y,g1x)/det
    const double a1 = (pu*g2y - pv*g2x)/det;
    const double a2 = (-qu*g1y + qv*g1x)/det;

    u = a1*g1x + a2*g2x;
    v = a1*g1y + a2*g2y;
}

double sediment_curv::face_sign(lexer *p, sediment_fdm *s, int dir, int ii, int jj)
{
    double nx,ny,len;
    normal(p,dir,ii,jj,nx,ny,len);

    const double un = (dir==1) ? s->P(ii,jj) : s->Q(ii,jj);

    // tangential part: the mean of the four crossing faces
    double tu=0.0, tv=0.0, a, b;

    if(dir==1)
    {
    facevec(p,s,2,ii,jj,a,b);     tu+=a; tv+=b;
    facevec(p,s,2,ii+1,jj,a,b);   tu+=a; tv+=b;
    facevec(p,s,2,ii,jj-1,a,b);   tu+=a; tv+=b;
    facevec(p,s,2,ii+1,jj-1,a,b); tu+=a; tv+=b;
    }

    if(dir==2)
    {
    facevec(p,s,1,ii,jj,a,b);     tu+=a; tv+=b;
    facevec(p,s,1,ii,jj+1,a,b);   tu+=a; tv+=b;
    facevec(p,s,1,ii-1,jj,a,b);   tu+=a; tv+=b;
    facevec(p,s,1,ii-1,jj+1,a,b); tu+=a; tv+=b;
    }

    const double ut = 0.25*(-tu*ny + tv*nx);
    const double mag = sqrt(un*un + ut*ut);

    return fabs(mag)>1.0e-10 ? un/fabs(mag) : 0.0;
}

void sediment_curv::bedgrad(lexer *p, slice &f, int ii, int jj, int zero, double &fx, double &fy)
{
    auto xc = [&](int a, int b) { return p->XC2D[q(p,a,b)]; };
    auto yc = [&](int a, int b) { return p->YC2D[q(p,a,b)]; };
    auto df = [&](int a, int b) { return p->DFBED[q(p,a,b)]; };

    // stencil per direction as bedslope::slope_cds / bedload_direction_f: forward next to a closed
    // cell at -1, backward next to a closed cell at +1 (both closed: backward, or no gradient with zero=1)
    int im=ii-1, ip=ii+1, jm=jj-1, jp=jj+1;
    double si=1.0, sj=1.0;

    if(df(ii-1,jj)<0)
    {
    im = ii; ip = ii+1;
    }
    if(df(ii+1,jj)<0)
    {
    im = ii-1; ip = ii;
    }
    if(zero==1 && df(ii-1,jj)<0 && df(ii+1,jj)<0)
    {
    im = ii-1; ip = ii+1;
    si = 0.0;
    }

    if(df(ii,jj-1)<0)
    {
    jm = jj; jp = jj+1;
    }
    if(df(ii,jj+1)<0)
    {
    jm = jj-1; jp = jj;
    }
    if(zero==1 && df(ii,jj-1)<0 && df(ii,jj+1)<0)
    {
    jm = jj-1; jp = jj+1;
    sj = 0.0;
    }

    // metrics from the same stencil (central where both sides are closed)
    const double xe = xc(ip,jj)-xc(im,jj), ye = yc(ip,jj)-yc(im,jj);
    const double xn = xc(ii,jp)-xc(ii,jm), yn = yc(ii,jp)-yc(ii,jm);
    const double fe = si*(f(ip,jj)-f(im,jj));
    const double fn = sj*(f(ii,jp)-f(ii,jm));

    const double det = xe*yn - ye*xn;

    if(fabs(det)>0.0)
    {
    fx = (fe*yn - ye*fn)/det;
    fy = (xe*fn - fe*xn)/det;
    }

    else
    fx = fy = 0.0;
}

double sediment_curv::contra(lexer *p, int dir, int ii, int jj, double u, double v)
{
    auto xc = [&](int a, int b) { return p->XC2D[q(p,a,b)]; };
    auto yc = [&](int a, int b) { return p->YC2D[q(p,a,b)]; };

    const double xe = 0.5*(xc(ii+1,jj)-xc(ii-1,jj)), ye = 0.5*(yc(ii+1,jj)-yc(ii-1,jj));
    const double xn = 0.5*(xc(ii,jj+1)-xc(ii,jj-1)), yn = 0.5*(yc(ii,jj+1)-yc(ii,jj-1));
    const double det = xe*yn - ye*xn;

    if(fabs(det)<1.0e-30)
    return 0.0;

    if(dir==1)
    return (u*yn - v*xn)/det;

    return (-u*ye + v*xe)/det;
}

double sediment_curv::width(lexer *p, int dir, int ii, int jj)
{
    const int c = q(p,ii,jj);
    double l;

    if(dir==1)
    l = 0.5*(p->LEN1[c] + p->LEN1[q(p,ii-1,jj)]);
    else
    l = 0.5*(p->LEN2[c] + p->LEN2[q(p,ii,jj-1)]);

    return (l>0.0) ? p->AREA2D[c]/l : 1.0e20;
}
