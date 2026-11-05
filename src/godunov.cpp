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

#include"godunov.h"
#include"lexer.h"
#include"fdm.h"
#include"field.h"
#include<cmath>

godunov::godunov(lexer *p)
{
}

godunov::~godunov()
{
}

// monotonised central limiter of the two one-sided slopes
double godunov::mc(double dl, double dr)
{
    if(dl*dr<=0.0)
    return 0.0;

    const double s = (dl>0.0) ? 1.0 : -1.0;
    return s*std::min(std::min(2.0*fabs(dl),2.0*fabs(dr)),0.5*fabs(dl+dr));
}

// the centre of the control volume idx of the variable ipol in direction d (the velocities on the
// faces of the cells, the level set and the scalars in the cells)
double godunov::pos(lexer *p, int ipol, int d, int idx) const
{
    if(d==0)
    return (ipol==1) ? p->XN[idx+1+marge] : p->XP[idx+marge];

    if(d==1)
    return (ipol==2) ? p->YN[idx+1+marge] : p->YP[idx+marge];

    return (ipol==3) ? p->ZN[idx+1+marge] : p->ZP[idx+marge];
}

// the width of the control volume idx of the variable ipol in direction d
double godunov::width(lexer *p, int ipol, int d, int idx) const
{
    if(d==0)
    return (ipol==1) ? p->XP[idx+1+marge]-p->XP[idx+marge] : p->XN[idx+1+marge]-p->XN[idx+marge];

    if(d==1)
    return (ipol==2) ? p->YP[idx+1+marge]-p->YP[idx+marge] : p->YN[idx+1+marge]-p->YN[idx+marge];

    return (ipol==3) ? p->ZP[idx+1+marge]-p->ZP[idx+marge] : p->ZN[idx+1+marge]-p->ZN[idx+marge];
}

// the advecting velocities on the faces of the control volume (as flux_face_CDS2)
void godunov::vx(int ipol, int ii, int jj, int kk, field &u, double &v1, double &v2)
{
    if(ipol==1)
    {
    v1 = 0.5*(u(ii,jj,kk)+u(ii-1,jj,kk));
    v2 = 0.5*(u(ii,jj,kk)+u(ii+1,jj,kk));
    }

    if(ipol==2)
    {
    v1 = 0.5*(u(ii-1,jj,kk)+u(ii-1,jj+1,kk));
    v2 = 0.5*(u(ii,jj,kk)+u(ii,jj+1,kk));
    }

    if(ipol==3)
    {
    v1 = 0.5*(u(ii-1,jj,kk)+u(ii-1,jj,kk+1));
    v2 = 0.5*(u(ii,jj,kk)+u(ii,jj,kk+1));
    }

    if(ipol==4)
    {
    v1 = u(ii-1,jj,kk);
    v2 = u(ii,jj,kk);
    }
}

void godunov::vy(int ipol, int ii, int jj, int kk, field &v, double &v1, double &v2)
{
    if(ipol==1)
    {
    v1 = 0.5*(v(ii,jj-1,kk)+v(ii+1,jj-1,kk));
    v2 = 0.5*(v(ii,jj,kk)+v(ii+1,jj,kk));
    }

    if(ipol==2)
    {
    v1 = 0.5*(v(ii,jj,kk)+v(ii,jj-1,kk));
    v2 = 0.5*(v(ii,jj,kk)+v(ii,jj+1,kk));
    }

    if(ipol==3)
    {
    v1 = 0.5*(v(ii,jj-1,kk)+v(ii,jj-1,kk+1));
    v2 = 0.5*(v(ii,jj,kk)+v(ii,jj,kk+1));
    }

    if(ipol==4)
    {
    v1 = v(ii,jj-1,kk);
    v2 = v(ii,jj,kk);
    }
}

void godunov::vz(int ipol, int ii, int jj, int kk, field &w, double &v1, double &v2)
{
    if(ipol==1)
    {
    v1 = 0.5*(w(ii,jj,kk-1)+w(ii+1,jj,kk-1));
    v2 = 0.5*(w(ii,jj,kk)+w(ii+1,jj,kk));
    }

    if(ipol==2)
    {
    v1 = 0.5*(w(ii,jj,kk-1)+w(ii,jj+1,kk-1));
    v2 = 0.5*(w(ii,jj,kk)+w(ii,jj+1,kk));
    }

    if(ipol==3)
    {
    v1 = 0.5*(w(ii,jj,kk)+w(ii,jj,kk-1));
    v2 = 0.5*(w(ii,jj,kk)+w(ii,jj,kk+1));
    }

    if(ipol==4)
    {
    v1 = w(ii,jj,kk-1);
    v2 = w(ii,jj,kk);
    }
}

void godunov::start(lexer *p, fdm *a, field &f, int ipol, field &U, field &V, field &W)
{
    if(ipol==2 && p->j_dir==0)
    return;

    if(ipol<1 || ipol>4)
    return;

    const bool ydir = (p->j_dir==1);
    const double hdt = 0.5*p->dt;

    // the extended range: the control volumes of the loops and two layers around them
    i0 = -2; ni = p->knox+4;
    j0 = ydir ? -2 : 0; nj = ydir ? p->knoy+4 : 1;
    k0 = -2; nk = p->knoz+4;
    const size_t n = size_t(ni)*nj*nk;

    for(std::vector<double> *v : {&sf,&sx,&sy,&sz,&uc,&vc,&wc,&tx,&ty,&tz})
    v->assign(n,0.0);

    for(int ii=i0; ii<i0+ni; ++ii)
    for(int jj=j0; jj<j0+nj; ++jj)
    for(int kk=k0; kk<k0+nk; ++kk)
    {
        const int q = id(ii,jj,kk);
        const double fc = f(ii,jj,kk);
        double v1, v2;

        // x
        const double xl = (fc-f(ii-1,jj,kk))/(pos(p,ipol,0,ii)-pos(p,ipol,0,ii-1));
        const double xr = (f(ii+1,jj,kk)-fc)/(pos(p,ipol,0,ii+1)-pos(p,ipol,0,ii));
        vx(ipol,ii,jj,kk,U,v1,v2);
        uc[q] = 0.5*(v1+v2);
        sx[q] = mc(xl,xr);
        tx[q] = uc[q]*((uc[q]>0.0) ? xl : xr);

        // y
        if(ydir)
        {
        const double yl = (fc-f(ii,jj-1,kk))/(pos(p,ipol,1,jj)-pos(p,ipol,1,jj-1));
        const double yr = (f(ii,jj+1,kk)-fc)/(pos(p,ipol,1,jj+1)-pos(p,ipol,1,jj));
        vy(ipol,ii,jj,kk,V,v1,v2);
        vc[q] = 0.5*(v1+v2);
        sy[q] = mc(yl,yr);
        ty[q] = vc[q]*((vc[q]>0.0) ? yl : yr);
        }

        // z
        const double zl = (fc-f(ii,jj,kk-1))/(pos(p,ipol,2,kk)-pos(p,ipol,2,kk-1));
        const double zr = (f(ii,jj,kk+1)-fc)/(pos(p,ipol,2,kk+1)-pos(p,ipol,2,kk));
        vz(ipol,ii,jj,kk,W,v1,v2);
        wc[q] = 0.5*(v1+v2);
        sz[q] = mc(zl,zr);
        tz[q] = wc[q]*((wc[q]>0.0) ? zl : zr);

        sf[q] = fc;
    }

    // the state on the face between the control volumes qa (upwind side for a positive velocity)
    // and qb, in direction d, with the face velocity un
    auto state = [&](int d, int qa, int qb, double wa, double wb, double un) -> double
    {
        double fa, fb;

        if(d==0)
        {
        fa = sf[qa] + (0.5*wa - hdt*uc[qa])*sx[qa] - hdt*(ty[qa]+tz[qa]);
        fb = sf[qb] - (0.5*wb + hdt*uc[qb])*sx[qb] - hdt*(ty[qb]+tz[qb]);
        }
        else if(d==1)
        {
        fa = sf[qa] + (0.5*wa - hdt*vc[qa])*sy[qa] - hdt*(tx[qa]+tz[qa]);
        fb = sf[qb] - (0.5*wb + hdt*vc[qb])*sy[qb] - hdt*(tx[qb]+tz[qb]);
        }
        else
        {
        fa = sf[qa] + (0.5*wa - hdt*wc[qa])*sz[qa] - hdt*(tx[qa]+ty[qa]);
        fb = sf[qb] - (0.5*wb + hdt*wc[qb])*sz[qb] - hdt*(tx[qb]+ty[qb]);
        }

        if(un>0.0)
        return fa;

        if(un<0.0)
        return fb;

        return 0.5*(fa+fb);
    };

    auto term = [&](int ii, int jj, int kk) -> double
    {
        const int q = id(ii,jj,kk);
        double v1, v2, L = 0.0;

        // x
        {
        vx(ipol,ii,jj,kk,U,v1,v2);
        const double w = width(p,ipol,0,ii), wl = width(p,ipol,0,ii-1), wr = width(p,ipol,0,ii+1);
        const double fe = state(0,q,id(ii+1,jj,kk),w,wr,v2);
        const double fw = state(0,id(ii-1,jj,kk),q,wl,w,v1);
        L -= 0.5*(v1+v2)*(fe-fw)/w;
        }

        // y
        if(ydir)
        {
        vy(ipol,ii,jj,kk,V,v1,v2);
        const double w = width(p,ipol,1,jj), wl = width(p,ipol,1,jj-1), wr = width(p,ipol,1,jj+1);
        const double fn = state(1,q,id(ii,jj+1,kk),w,wr,v2);
        const double fs = state(1,id(ii,jj-1,kk),q,wl,w,v1);
        L -= 0.5*(v1+v2)*(fn-fs)/w;
        }

        // z
        {
        vz(ipol,ii,jj,kk,W,v1,v2);
        const double w = width(p,ipol,2,kk), wl = width(p,ipol,2,kk-1), wr = width(p,ipol,2,kk+1);
        const double ft = state(2,q,id(ii,jj,kk+1),w,wr,v2);
        const double fb = state(2,id(ii,jj,kk-1),q,wl,w,v1);
        L -= 0.5*(v1+v2)*(ft-fb)/w;
        }

        return L;
    };

    if(ipol==1)
    ULOOP
    a->F(i,j,k) += term(i,j,k);

    if(ipol==2)
    VLOOP
    a->G(i,j,k) += term(i,j,k);

    if(ipol==3)
    WLOOP
    a->H(i,j,k) += term(i,j,k);

    if(ipol==4)
    LOOP
    a->L(i,j,k) += term(i,j,k);
}
