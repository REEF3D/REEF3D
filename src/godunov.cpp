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
#include"interface_width.h"
#include"cfd_airband_dt.h"
#include<cmath>
#include<algorithm>

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

// N 59: the density of the control volume of the variable ipol, that of the projection on the face
// of the cells (density_f::roface: the Heaviside of the mean of the level set of the two cells)
double godunov::rho_cv(lexer *p, field &phi, int ipol, int ii, int jj, int kk) const
{
    const int aa = (ipol==1) ? 1 : 0;
    const int bb = (ipol==2) ? 1 : 0;
    const int cc = (ipol==3) ? 1 : 0;

    const double phival = 0.5*(phi(ii,jj,kk) + phi(ii+aa,jj+bb,kk+cc));
    const double H = heaviside_face(p,phi,phival,ii,jj,kk,aa,bb,cc);

    return p->W1*H + p->W3*(1.0-H);
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

        // N 56: the momentum first-order upwind in the air further than N56 cells from the free
        // surface (the slopes zero, the transverse terms kept)
        if(ipol<=3 && p->N56>0.0 && p->F30>0)
        {
            double ph = a->phi(ii,jj,kk);
            if(ipol==1) ph = 0.5*(ph+a->phi(ii+1,jj,kk));
            if(ipol==2) ph = 0.5*(ph+a->phi(ii,jj+1,kk));
            if(ipol==3) ph = 0.5*(ph+a->phi(ii,jj,kk+1));

            double h = std::min(width(p,ipol,0,ii),width(p,ipol,2,kk));
            if(ydir)
            h = std::min(h,width(p,ipol,1,jj));

            if(ph < -p->N56*h)
            sx[q]=sy[q]=sz[q]=0.0;
        }

        sf[q] = fc;
    }

    // N 59: the density of the control volumes and its reconstruction, as that of f (the densities
    // of the outermost layer are those of the layer inside it: the level set of its interface
    // width would be read beyond the ghost cells)
    const bool cons = (ipol<=3 && p->N59>0 && p->F30>0);

    if(cons)
    {
        field &ph = (cphi!=nullptr) ? *cphi : a->phi;

        for(std::vector<double> *v : {&rf,&rx,&ry,&rz,&rtx,&rty,&rtz})
        v->assign(n,0.0);

        for(int ii=i0; ii<i0+ni; ++ii)
        for(int jj=j0; jj<j0+nj; ++jj)
        for(int kk=k0; kk<k0+nk; ++kk)
        rf[id(ii,jj,kk)] = rho_cv(p,ph,ipol,std::min(ii,p->knox),ydir ? std::min(jj,p->knoy) : jj,std::min(kk,p->knoz));

        for(int ii=i0+1; ii<i0+ni-1; ++ii)
        for(int jj=(ydir ? j0+1 : j0); jj<(ydir ? j0+nj-1 : j0+nj); ++jj)
        for(int kk=k0+1; kk<k0+nk-1; ++kk)
        {
            const int q = id(ii,jj,kk);
            const double rc = rf[q];

            const double xl = (rc-rf[id(ii-1,jj,kk)])/(pos(p,ipol,0,ii)-pos(p,ipol,0,ii-1));
            const double xr = (rf[id(ii+1,jj,kk)]-rc)/(pos(p,ipol,0,ii+1)-pos(p,ipol,0,ii));
            rx[q] = mc(xl,xr);
            rtx[q] = uc[q]*((uc[q]>0.0) ? xl : xr);

            if(ydir)
            {
            const double yl = (rc-rf[id(ii,jj-1,kk)])/(pos(p,ipol,1,jj)-pos(p,ipol,1,jj-1));
            const double yr = (rf[id(ii,jj+1,kk)]-rc)/(pos(p,ipol,1,jj+1)-pos(p,ipol,1,jj));
            ry[q] = mc(yl,yr);
            rty[q] = vc[q]*((vc[q]>0.0) ? yl : yr);
            }

            const double zl = (rc-rf[id(ii,jj,kk-1)])/(pos(p,ipol,2,kk)-pos(p,ipol,2,kk-1));
            const double zr = (rf[id(ii,jj,kk+1)]-rc)/(pos(p,ipol,2,kk+1)-pos(p,ipol,2,kk));
            rz[q] = mc(zl,zr);
            rtz[q] = wc[q]*((wc[q]>0.0) ? zl : zr);
        }
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

    // N 59: the face state of the density, as state() for f
    auto rstate = [&](int d, int qa, int qb, double wa, double wb, double un) -> double
    {
        double fa, fb;

        if(d==0)
        {
        fa = rf[qa] + (0.5*wa - hdt*uc[qa])*rx[qa] - hdt*(rty[qa]+rtz[qa]);
        fb = rf[qb] - (0.5*wb + hdt*uc[qb])*rx[qb] - hdt*(rty[qb]+rtz[qb]);
        }
        else if(d==1)
        {
        fa = rf[qa] + (0.5*wa - hdt*vc[qa])*ry[qa] - hdt*(rtx[qa]+rtz[qa]);
        fb = rf[qb] - (0.5*wb + hdt*vc[qb])*ry[qb] - hdt*(rtx[qb]+rtz[qb]);
        }
        else
        {
        fa = rf[qa] + (0.5*wa - hdt*wc[qa])*rz[qa] - hdt*(rtx[qa]+rty[qa]);
        fb = rf[qb] - (0.5*wb + hdt*wc[qb])*rz[qb] - hdt*(rtx[qb]+rty[qb]);
        }

        if(un>0.0)
        return fa;

        if(un<0.0)
        return fb;

        return 0.5*(fa+fb);
    };

    // N 59 2, 3: the face density of the mass flux from the level set at the half step of its own
    // transport, (phi_start + phi)/2, where the level set has already been moved with the advecting
    // velocity of this update (cphi: the level set before that transport; N 40 5 with the level set
    // kept before its last transport, the full step of N 40 6); otherwise the Godunov states of N 59 1
    const bool hls = cons && p->N59>=2 && cphi!=nullptr;

    auto phh = [&](int ii, int jj, int kk) -> double
    {
        return 0.5*((*cphi)(ii,jj,kk) + a->phi(ii,jj,kk));
    };

    // the density on the upper face in direction d of the control volume (ii,jj,kk): in the direction
    // of the staggering the face is the centre of the cell above, in the others the edge of the four
    // cells around it (the mean of the level set there)
    auto hface = [&](int d, int ii, int jj, int kk) -> double
    {
        const int s = ipol-1;
        const int es[3] = {s==0, s==1, s==2};
        const int ed[3] = {d==0, d==1, d==2};
        double x, H;

        if(d==s)
        {
            const int ia = ii+es[0], ja = jj+es[1], ka = kk+es[2];
            x = phh(ia,ja,ka);
            H = heaviside_face(p,a->phi,x,ia,ja,ka,0,0,0);
        }
        else
        {
            x = 0.25*(phh(ii,jj,kk) + phh(ii+es[0],jj+es[1],kk+es[2])
                    + phh(ii+ed[0],jj+ed[1],kk+ed[2]) + phh(ii+es[0]+ed[0],jj+es[1]+ed[1],kk+es[2]+ed[2]));
            H = heaviside_face(p,a->phi,x,ii,jj,kk,es[0],es[1],es[2]);
        }

        return p->W1*H + p->W3*(1.0-H);
    };

    // N 59: the convection term (u* - u)/dt of the update u* = (rho u - dt div(m u))/(rho - dt div(m))
    const double dtu = (cdt>0.0) ? cdt : p->dt;
    const double rmin = 0.5*std::min(p->W1,p->W3);

    // N 59 3: the momentum flux m u of a face into the control volume q. The face density of N 59 2 is
    // the level set at the face, not the density of the fluid on the upwind side: where light fluid
    // flows in and the level set there says heavier (air pushing into the water, an interface that the
    // level set moved without a flux), m u carried the velocity of the light fluid with the heavy
    // density, and the update gave the control volume the velocity of the air (the range of the face
    // states holds it, and the control volume gains mass instead of losing it). The face density
    // above that of the upwind control volume (the larger of its density at the start of the transport
    // and of the current level set) carries the velocity of q itself; the mass flux stays that of
    // N 59 2. Outflow faces carry their face state as before.
    const bool hup = hls && p->N59==3;

    auto mflux = [&](double m, double vn, bool in, int qu, int iu, int ju, int ku, double us, double uq) -> double
    {
        if(!hup || !in)
        return m*us;

        const double ru = std::max(rf[qu], rho_cv(p,a->phi,ipol,std::min(iu,p->knox),
                                                  (p->j_dir==1) ? std::min(ju,p->knoy) : ju,std::min(ku,p->knoz)));
        const double r = m/vn;

        if(r<=ru)
        return m*us;

        return vn*(ru*us + (r-ru)*uq);
    };

    auto term_cons = [&](int ii, int jj, int kk) -> double
    {
        const int q = id(ii,jj,kk);
        double v1, v2, Dm = 0.0, Du = 0.0, Dv = 0.0;
        double umin = sf[q], umax = sf[q];

        // x
        {
        vx(ipol,ii,jj,kk,U,v1,v2);
        const int qe = id(ii+1,jj,kk), qw = id(ii-1,jj,kk);
        const double w = width(p,ipol,0,ii), wl = width(p,ipol,0,ii-1), wr = width(p,ipol,0,ii+1);
        const double m2 = v2*(hls ? hface(0,ii,jj,kk) : rstate(0,q,qe,w,wr,v2));
        const double m1 = v1*(hls ? hface(0,ii-1,jj,kk) : rstate(0,qw,q,wl,w,v1));
        Dm += (m2-m1)/w;
        const double u2 = state(0,q,qe,w,wr,v2), u1 = state(0,qw,q,wl,w,v1);
        umin = std::min(umin,std::min(u1,u2));  umax = std::max(umax,std::max(u1,u2));
        Du += (mflux(m2,v2,v2<0.0,qe,ii+1,jj,kk,u2,sf[q]) - mflux(m1,v1,v1>0.0,qw,ii-1,jj,kk,u1,sf[q]))/w;
        Dv += (v2-v1)/w;
        }

        // y
        if(ydir)
        {
        vy(ipol,ii,jj,kk,V,v1,v2);
        const int qn = id(ii,jj+1,kk), qs = id(ii,jj-1,kk);
        const double w = width(p,ipol,1,jj), wl = width(p,ipol,1,jj-1), wr = width(p,ipol,1,jj+1);
        const double m2 = v2*(hls ? hface(1,ii,jj,kk) : rstate(1,q,qn,w,wr,v2));
        const double m1 = v1*(hls ? hface(1,ii,jj-1,kk) : rstate(1,qs,q,wl,w,v1));
        Dm += (m2-m1)/w;
        const double u2 = state(1,q,qn,w,wr,v2), u1 = state(1,qs,q,wl,w,v1);
        umin = std::min(umin,std::min(u1,u2));  umax = std::max(umax,std::max(u1,u2));
        Du += (mflux(m2,v2,v2<0.0,qn,ii,jj+1,kk,u2,sf[q]) - mflux(m1,v1,v1>0.0,qs,ii,jj-1,kk,u1,sf[q]))/w;
        Dv += (v2-v1)/w;
        }

        // z
        {
        vz(ipol,ii,jj,kk,W,v1,v2);
        const int qt = id(ii,jj,kk+1), qb = id(ii,jj,kk-1);
        const double w = width(p,ipol,2,kk), wl = width(p,ipol,2,kk-1), wr = width(p,ipol,2,kk+1);
        const double m2 = v2*(hls ? hface(2,ii,jj,kk) : rstate(2,q,qt,w,wr,v2));
        const double m1 = v1*(hls ? hface(2,ii,jj,kk-1) : rstate(2,qb,q,wl,w,v1));
        Dm += (m2-m1)/w;
        const double u2 = state(2,q,qt,w,wr,v2), u1 = state(2,qb,q,wl,w,v1);
        umin = std::min(umin,std::min(u1,u2));  umax = std::max(umax,std::max(u1,u2));
        Du += (mflux(m2,v2,v2<0.0,qt,ii,jj,kk+1,u2,sf[q]) - mflux(m1,v1,v1>0.0,qb,ii,jj,kk-1,u1,sf[q]))/w;
        Dv += (v2-v1)/w;
        }

        // the advecting velocities of the control volume (means of the face velocities of the cells)
        // are not exactly divergence free; their divergence Dv is taken out of the mass and the
        // momentum fluxes alike, so that a constant density stays constant and a light control volume
        // next to diverging velocities keeps its mass
        const double rnew = rf[q] - dtu*(Dm - rf[q]*Dv);

        // the second-order update where it keeps the mass of the control volume and the velocity
        // within the range of the velocity of the control volume and its face states; a large
        // relative loss of mass with outflow states other than the velocity of the control volume
        // would extrapolate the velocity, u* = (rho u - m u_out)/(rho - m), and let it grow
        if(rnew >= rmin && rnew >= 0.5*rf[q])
        {
            const double un = (rf[q]*sf[q] - dtu*(Du - rf[q]*sf[q]*Dv))/rnew;

            if(un >= umin && un <= umax)
            return (un - sf[q])/dtu;
        }

        // otherwise (a light control volume next to a heavy one, whose reconstructed outflow density
        // can exceed its own, or a velocity out of range) the same update with the first-order
        // upwind states of the density and the velocity: u* = u + dt sum m_in (u_in - u)/rho*, a
        // convex combination of the velocities of the control volume and its inflow neighbours;
        // where even that fails, the advective form
        {
            double Dm1 = 0.0, Du1 = 0.0;

            auto up = [&](int qa, int qb, double un, double &r, double &u)
            {
                const int qu = (un>=0.0) ? qa : qb;
                r = rf[qu];
                u = sf[qu];
            };

            double r2, u2, r1, u1;

            vx(ipol,ii,jj,kk,U,v1,v2);
            up(q,id(ii+1,jj,kk),v2,r2,u2);   up(id(ii-1,jj,kk),q,v1,r1,u1);
            Dm1 += (v2*r2 - v1*r1)/width(p,ipol,0,ii);
            Du1 += (v2*r2*u2 - v1*r1*u1)/width(p,ipol,0,ii);

            if(ydir)
            {
            vy(ipol,ii,jj,kk,V,v1,v2);
            up(q,id(ii,jj+1,kk),v2,r2,u2);   up(id(ii,jj-1,kk),q,v1,r1,u1);
            Dm1 += (v2*r2 - v1*r1)/width(p,ipol,1,jj);
            Du1 += (v2*r2*u2 - v1*r1*u1)/width(p,ipol,1,jj);
            }

            vz(ipol,ii,jj,kk,W,v1,v2);
            up(q,id(ii,jj,kk+1),v2,r2,u2);   up(id(ii,jj,kk-1),q,v1,r1,u1);
            Dm1 += (v2*r2 - v1*r1)/width(p,ipol,2,kk);
            Du1 += (v2*r2*u2 - v1*r1*u1)/width(p,ipol,2,kk);

            const double rnew1 = rf[q] - dtu*(Dm1 - rf[q]*Dv);

            if(rnew1 >= rmin)
            return ((rf[q]*sf[q] - dtu*(Du1 - rf[q]*sf[q]*Dv))/rnew1 - sf[q])/dtu;
        }

        return term(ii,jj,kk);
    };

    // N 65 1 (with N 57): the momentum of the air beyond the band of N 57 semi-Lagrangian. The value at the
    // departure point x - dt U of the control volume (U its advecting velocity uc, vc, wc), multilinear
    // between the control volumes; the displacement at most cfd_sl_reach cells per direction (the time
    // step keeps it there, cfd_airband_dt.h, cfd_amr::timestep) and within the extended range. A
    // convex combination of the values around the departure point: stable at any Courant number, so
    // that the fast air away from the free surface (pushed out of a closing gap, the roller of a jump,
    // under a deck) no longer sets the time step. The water and the band keep the Godunov convection
    // (and N 59).
    const bool sl = (ipol<=3 && p->N65>0 && p->N57>0.0 && p->F30>0);

    auto air = [&](int ii, int jj, int kk) -> bool
    {
        const int aa = (ipol==1), bb = (ipol==2), cc = (ipol==3);
        const double ph = 0.5*(a->phi(ii,jj,kk) + a->phi(ii+aa,jj+bb,kk+cc));
        double h = std::min(width(p,ipol,0,ii),width(p,ipol,2,kk));
        if(ydir)
        h = std::min(h,width(p,ipol,1,jj));
        return ph < -p->N57*h;
    };

    auto sl_term = [&](int ii, int jj, int kk) -> double
    {
        const int q = id(ii,jj,kk);
        const double r = cfd_sl_reach;

        auto loc = [&](int d, int c0, int lo, int n, double vel, int &ib, double &t)
        {
            const double s = std::max(-r,std::min(r,dtu*vel/width(p,ipol,d,c0)));
            const double x = double(c0) - s;
            ib = int(std::floor(x));
            ib = std::max(lo,std::min(lo+n-2,ib));
            t = std::max(0.0,std::min(1.0,x-double(ib)));
        };

        int ib, jb = jj, kb;
        double ti, tj = 0.0, tk;
        loc(0,ii,i0,ni,uc[q],ib,ti);
        loc(2,kk,k0,nk,wc[q],kb,tk);
        if(ydir)
        loc(1,jj,j0,nj,vc[q],jb,tj);

        double val = 0.0;
        for(int a1=0; a1<2; ++a1)
        for(int b1=0; b1<(ydir ? 2 : 1); ++b1)
        for(int c1=0; c1<2; ++c1)
        {
            const double wgt = (a1 ? ti : 1.0-ti)*(ydir ? (b1 ? tj : 1.0-tj) : 1.0)*(c1 ? tk : 1.0-tk);
            val += wgt*sf[id(ib+a1,jb+b1,kb+c1)];
        }

        return (val - sf[q])/dtu;
    };

    if(cons)
    {
        if(ipol==1)
        ULOOP
        a->F(i,j,k) += (sl && air(i,j,k)) ? sl_term(i,j,k) : term_cons(i,j,k);

        if(ipol==2)
        VLOOP
        a->G(i,j,k) += (sl && air(i,j,k)) ? sl_term(i,j,k) : term_cons(i,j,k);

        if(ipol==3)
        WLOOP
        a->H(i,j,k) += (sl && air(i,j,k)) ? sl_term(i,j,k) : term_cons(i,j,k);

        return;
    }

    if(ipol==1)
    ULOOP
    a->F(i,j,k) += (sl && air(i,j,k)) ? sl_term(i,j,k) : term(i,j,k);

    if(ipol==2)
    VLOOP
    a->G(i,j,k) += (sl && air(i,j,k)) ? sl_term(i,j,k) : term(i,j,k);

    if(ipol==3)
    WLOOP
    a->H(i,j,k) += (sl && air(i,j,k)) ? sl_term(i,j,k) : term(i,j,k);

    if(ipol==4)
    LOOP
    a->L(i,j,k) += term(i,j,k);
}
