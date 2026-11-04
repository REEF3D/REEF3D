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

        velocity(p,a,x0,y0,z0,u1,v1,w1);
        velocity(p,a,x0+0.5*dt*u1, y0+0.5*dt*v1, z0+0.5*dt*w1, u2,v2,w2);
        velocity(p,a,x0+0.75*dt*u2, y0+0.75*dt*v2, z0+0.75*dt*w2, u3,v3,w3);

        q.x = x0 + dt/9.0*(2.0*u1 + 3.0*u2 + 4.0*u3);
        q.y = y0 + dt/9.0*(2.0*v1 + 3.0*v2 + 4.0*v3);
        q.z = z0 + dt/9.0*(2.0*w1 + 3.0*w2 + 4.0*w3);

        if(!domain(p,q.x,q.y,q.z))
        {
            ++lost;
            P[n] = P.back();
            P.pop_back();
            continue;
        }

        slide(p,a,q,x0,y0,z0);

        velocity(p,a,q.x,q.y,q.z,q.u,q.v,q.w);

        ++n;
    }

    lost = (long long)pgc->globalsum(double(lost));
    lost_total += lost;
}

void pic_f::velocity(lexer *p, fdm *a, double xp, double yp, double zp, double &up, double &vp, double &wp)
{
    up = ipol(p,a->u,1,xp,yp,zp);
    vp = p->j_dir==1 ? ipol(p,a->v,2,xp,yp,zp) : 0.0;
    wp = ipol(p,a->w,3,xp,yp,zp);
}

// a solid at the position: a blocked cell (flag4 < 0: topography, solid, object of the grid;
// the ghost cells of the neighbour subdomains have their flags) or the immersed solid,
// topography or floating body (X 10 1) (level sets < 0). The static loop indices are left
// alone (seed calls it in LOOP).
bool pic_f::blocked(lexer *p, fdm *a, double xp, double yp, double zp)
{
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

            if(f==INFLOW_FLAG || f==OUTFLOW_FLAG)
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
