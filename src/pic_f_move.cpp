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

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

// Phase 0: the particles move with the CFD velocity of the new time level, frozen over the
// step (as the particles of PICIN move through the divergence-free grid velocity), with
// Ralston's third-order Runge-Kutta scheme (Chen et al. 2018):
//   k1 = u(x^n), k2 = u(x^n + dt/2 k1), k3 = u(x^n + 3dt/4 k2), x^n+1 = x^n + dt/9 (2 k1 + 3 k2 + 4 k3)
// Domain walls reflect the particles; at an inflow or outflow boundary (flag1) they leave the
// domain. A particle that ends in a solid (flag4 < 0, or the solid / topo level set < 0) slides
// back: the components of its step that take it there are undone.

void pic_f::move(lexer *p, fdm *a, ghostcell *pgc)
{
    const double dt = p->dt;
    lost = 0;

    size_t n=0;
    while(n<P.size())
    {
        pic_data &q = P[n];

        const double x0=q.x, y0=q.y, z0=q.z;
        double u1,v1,w1, u2,v2,w2, u3,v3,w3;

        // F 81 1: spray more than 2 cells away from the water with its own velocity (pic_f::g2p)
        if(single && ballistic(p,a,q))
        {
        q.x = x0 + dt*q.u;
        q.y = y0 + dt*q.v;
        q.z = z0 + dt*q.w;
        }
        else
        {
        velocity(p,a,x0,y0,z0,u1,v1,w1);
        velocity(p,a,x0+0.5*dt*u1, y0+0.5*dt*v1, z0+0.5*dt*w1, u2,v2,w2);
        velocity(p,a,x0+0.75*dt*u2, y0+0.75*dt*v2, z0+0.75*dt*w2, u3,v3,w3);

        q.x = x0 + dt/9.0*(2.0*u1 + 3.0*u2 + 4.0*u3);
        q.y = y0 + dt/9.0*(2.0*v1 + 3.0*v2 + 4.0*v3);
        q.z = z0 + dt/9.0*(2.0*w1 + 3.0*w2 + 4.0*w3);
        }

        if(!domain(p,q.x,q.y,q.z))
        {
            ++lost;
            P[n] = P.back();
            P.pop_back();
            continue;
        }

        slide(p,a,q,x0,y0,z0);

        // floating bodies (X 10): a body that moved onto the particle pushes it out
        if(p->X10>0)
        push_out(p,a,q);

        // the velocity at the new position (F 81 1, 2: the particles carry their own, pic_f::g2p)
        if(!flipmode)
        velocity(p,a,q.x,q.y,q.z,q.u,q.v,q.w);

        ++n;
    }

    lost = (long long)pgc->globalsum(double(lost));
    lost_total += lost;
}

// F 81 1, 2 with F 76 1: the particles move with the stages of the SSP-RK3 step of the momentum
// (momentum_rk::step_ssp), as the level set does inside the stages (N 40 3):
//     x1 = xn + dt u(n)(xn)                      before stage 1, s = 1
//     x2 = 3/4 xn + 1/4 (x1 + dt u1(x1))         before stage 2, s = 2
//     x(n+1) = 1/3 xn + 2/3 (x2 + dt u2(x2))     after the step, s = 3
// with u(n), u1, u2 the input velocities of the stages (extended into the air), and the free surface
// from the particles before stages 1 and 2, so that each stage has the surface of its time (as the
// level set in the stages of N 40 3). Moved once per step with the frozen u(n+1) (F 76 0), the
// sloshing tank at 2 cells per wave height and CFL 0.25 lost per period 3.2 % of its energy with
// F 81 2 and 4.0 % with F 81 1, less with a smaller time step; with the stages 0.5 % and 1.1 %
// (the level set 0.0 %). The velocity of the particles (pic_f::g2p) is that at their position at
// the start of the step (xn).
void pic_f::stage(lexer *p, fdm *a, ghostcell *pgc, int s, field &U, field &V, field &W)
{
    if(!rkstages)
    return;

    extend(p,a,pgc,U,*uw,0,3,single);
    extend(p,a,pgc,V,*vw,1,3,single);
    extend(p,a,pgc,W,*ww,2,3,single);

    const double dt = p->dt;
    const double c0 = s==1 ? 0.0 : (s==2 ? 0.75 : 1.0/3.0);   // weight of xn

    if(s==1)
    lost = 0;

    size_t n=0;
    while(n<P.size())
    {
        pic_data &q = P[n];

        if(s==1)
        {
        q.xn[0] = q.x;
        q.xn[1] = q.y;
        q.xn[2] = q.z;
        }

        const double x0=q.x, y0=q.y, z0=q.z;
        double up,vp,wp;

        // F 81 1: spray more than 2 cells away from the water with its own velocity (pic_f::g2p)
        if(single && ballistic(p,a,q))
        {
        up = q.u;
        vp = q.v;
        wp = q.w;
        }
        else
        velocity(p,a,x0,y0,z0,up,vp,wp);

        q.x = c0*q.xn[0] + (1.0-c0)*(x0 + dt*up);
        q.y = p->j_dir==1 ? c0*q.xn[1] + (1.0-c0)*(y0 + dt*vp) : y0;
        q.z = c0*q.xn[2] + (1.0-c0)*(z0 + dt*wp);

        // out through an open boundary: removed at the end of the step; in the stages before it
        // goes on beyond the boundary (in the ghost cells, with the velocity there), its position
        // of the step may come back, and the free surface of the stages leaves it out (phi_build,
        // as that of the end of the step without it); walls reflect it
        double xd=q.x, yd=q.y, zd=q.z;
        const bool inside = domain(p,xd,yd,zd);

        if(!inside && s==3)
        {
            ++lost;
            P[n] = P.back();
            P.pop_back();
            continue;
        }

        if(inside)
        {
        q.x = xd;
        q.y = yd;
        q.z = zd;

        slide(p,a,q,x0,y0,z0);

        if(p->X10>0)
        push_out(p,a,q);
        }

        ++n;
    }

    xchange(p,pgc);

    if(s<3)
    phi_build(p,a,pgc,true);

    if(s==3)
    {
    lost = (long long)pgc->globalsum(double(lost));
    lost_total += lost;
    staged = true;
    }
}

void pic_f::velocity(lexer *p, fdm *a, double xp, double yp, double zp, double &up, double &vp, double &wp)
{
    // F 81 1, 2: the velocity of the water extended into the air (pic_f_extend.cpp)
    field &fu = surface ? (field&)*uw : (field&)a->u;
    field &fv = surface ? (field&)*vw : (field&)a->v;
    field &fw = surface ? (field&)*ww : (field&)a->w;

    up = ipol(p,fu,1,xp,yp,zp);
    vp = p->j_dir==1 ? ipol(p,fv,2,xp,yp,zp) : 0.0;
    wp = ipol(p,fw,3,xp,yp,zp);
}

// the distance to the nearest solid surface (a->solid, a->topo) and its normal out of the solid;
// false without a solid within a cell
bool pic_f::solid_distance(lexer *p, fdm *a, double xp, double yp, double zp, double &f, double *n)
{
    auto level = [&](double x, double y, double z)
    {
        double v = 1.0e20;
        if(p->toporead>0)
        v = ipol(p,a->topo,4,x,y,z);
        if(p->solidread>0)
        v = std::min(v, ipol(p,a->solid,4,x,y,z));
        return v;
    };

    const int ci = std::min(std::max(cell(p->XN,p->knox,xp),0),p->knox-1);
    const double h = p->DXN[ci+marge];

    f = level(xp,yp,zp);

    if(f>h)
    return false;

    const double d = 0.5*h;
    n[0] = (level(xp+d,yp,zp) - level(xp-d,yp,zp))/(2.0*d);
    n[1] = p->j_dir==1 ? (level(xp,yp+d,zp) - level(xp,yp-d,zp))/(2.0*d) : 0.0;
    n[2] = (level(xp,yp,zp+d) - level(xp,yp,zp-d))/(2.0*d);
    const double g = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

    if(g<1.0e-12)
    return false;

    n[0]/=g;
    n[1]/=g;
    n[2]/=g;
    return true;
}

// a solid at the position: a blocked cell (flag4 < 0: topography, solid, object of the grid;
// the ghost cells of the neighbour subdomains have their flags) or the immersed solid,
// topography or floating body (X 10 1) (level sets < 0). The static loop indices are left
// alone (seed calls it in LOOP).
bool pic_f::blocked(lexer *p, fdm *a, double xp, double yp, double zp)
{
    // B 60: the ghost cells of the inflow (pic_f_inflow.cpp)
    if(inflowsrc && xp<xglo_s)
    return ipol(p,a->topo,4,xp,yp,zp)<0.0;

    const int ii = cell(p->XN,p->knox,xp);
    const int jj = p->j_dir==1 ? cell(p->YN,p->knoy,yp) : 0;
    const int kk = cell(p->ZN,p->knoz,zp);

    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    return true;

    if(p->solidread>0 && ipol(p,a->solid,4,xp,yp,zp)<0.0)
    return true;

    if(ipol(p,a->topo,4,xp,yp,zp)<0.0)
    return true;

    if(p->X10==1 && ipol(p,a->fb,4,xp,yp,zp)<0.0)
    return true;

    return false;
}

// the domain boundaries: false when the particle leaves through an inflow or outflow
// boundary (the u-flag of the boundary face); domain walls reflect it
bool pic_f::domain(lexer *p, double &xp, double &yp, double &zp)
{
    if(xp<xglo_s || xp>=xglo_e)
    {
        int jj = p->j_dir==1 ? cell(p->YN,p->knoy,yp) : 0;
        int kk = cell(p->ZN,p->knoz,zp);
        jj = std::max(0, std::min(jj, p->knoy-1));
        kk = std::max(0, std::min(kk, p->knoz-1));

        const int ib = xp<xglo_s ? -1 : p->knox;

        if((xp<xglo_s && xloc_s==xglo_s) || (xp>=xglo_e && xloc_e==xglo_e))
        {
            const int f = p->flag1[(ib-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin];

            // B 60: the particles of the inflow live in its three ghost layers (pic_f_inflow.cpp)
            if(inflowsrc && f==INFLOW_FLAG && xp<xglo_s)
            {
                if(xp<p->XN[marge-3])
                return false;

                auto clampw = [](double &c, double s, double e)
                {
                    const double eps = 1.0e-9*(e-s);
                    if(c<s)
                    c = 2.0*s - c;
                    if(c>=e)
                    c = 2.0*e - c;
                    c = std::max(s+eps, std::min(c, e-eps));
                };

                if(p->j_dir==1)
                clampw(yp,yglo_s,yglo_e);
                clampw(zp,zglo_s,zglo_e);
                return true;
            }

            // the relaxation zone of the wave generation (B 98 2) has a wall at its boundary: the
            // particles are reflected there (removed, the zone lost its water to the waves of a
            // steep wave: 70 % of the tank in 24 s for the cnoidal waves of the spilling breaker)
            if((f==INFLOW_FLAG && !relaxzones) || f==OUTFLOW_FLAG)
            return false;
        }
    }

    // reflection at the walls, and strictly inside the domain
    auto reflect = [](double &c, double s, double e)
    {
        const double eps = 1.0e-9*(e-s);

        if(c<s)
        c = 2.0*s - c;

        if(c>=e)
        c = 2.0*e - c;

        c = std::max(s+eps, std::min(c, e-eps));
    };

    reflect(xp,xglo_s,xglo_e);
    if(p->j_dir==1)
    reflect(yp,yglo_s,yglo_e);
    reflect(zp,zglo_s,zglo_e);

    return true;
}

// a particle that ended in a solid: the components of its step that take it there are
// undone, the vertical one first (the bed), then x and y, then two of them, then all
void pic_f::slide(lexer *p, fdm *a, pic_data &q, double x0, double y0, double z0)
{
    if(!blocked(p,a,q.x,q.y,q.z))
    return;

    const double xn=q.x, yn=q.y, zn=q.z;

    // which components are undone: bit 0 x, bit 1 y, bit 2 z
    const int tries[6] = {4, 1, 2, 5, 6, 3};

    for(int t=0; t<6; ++t)
    {
        const int m = tries[t];

        const double xc = (m&1) ? x0 : xn;
        const double yc = (m&2) ? y0 : yn;
        const double zc = (m&4) ? z0 : zn;

        if(!blocked(p,a,xc,yc,zc))
        {
            q.x=xc; q.y=yc; q.z=zc;
            return;
        }
    }

    q.x=x0; q.y=y0; q.z=z0;
}

// a particle inside a floating body (its level set fb < 0, the body moved onto it) goes out along
// the gradient of fb, to 0.1 cells outside
void pic_f::push_out(lexer *p, fdm *a, pic_data &q)
{
    const double f = ipol(p,a->fb,4,q.x,q.y,q.z);

    if(f>=0.0)
    return;

    const int ci = cell(p->XN,p->knox,q.x);
    const double h = p->DXN[ci+marge];
    const double d = 0.5*h;

    const double gx = (ipol(p,a->fb,4,q.x+d,q.y,q.z) - ipol(p,a->fb,4,q.x-d,q.y,q.z))/(2.0*d);
    const double gy = p->j_dir==1 ? (ipol(p,a->fb,4,q.x,q.y+d,q.z) - ipol(p,a->fb,4,q.x,q.y-d,q.z))/(2.0*d) : 0.0;
    const double gz = (ipol(p,a->fb,4,q.x,q.y,q.z+d) - ipol(p,a->fb,4,q.x,q.y,q.z-d))/(2.0*d);
    const double g = std::sqrt(gx*gx + gy*gy + gz*gz);

    if(g<1.0e-12)
    return;

    const double s = (-f + 0.1*h)/g;
    double x = q.x + s*gx/g;
    double y = q.y + s*gy/g;
    double z = q.z + s*gz/g;

    if(!domain(p,x,y,z))
    return;

    q.x = x;
    q.y = y;
    q.z = z;
}
