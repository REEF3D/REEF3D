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

#include"fnpf_ice.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<algorithm>

void fnpf_ice::kinematics(fnpf_ice_floe &fl)
{
    quat_to_matrix(fl.q,fl.R);

    // bottom face: normal R*e_z, centre x - h/2*nb
    for(int q=0;q<3;++q)
    {
    fl.nb[q] = fl.R[q][2];
    fl.pb[q] = fl.x[q] - 0.5*fl.h*fl.nb[q];
    }

    const int nv = int(fl.bx.size());
    fl.wx.resize(nv);
    fl.wy.resize(nv);

    fl.bbox[0]=fl.bbox[2]= 1.0e20;
    fl.bbox[1]=fl.bbox[3]=-1.0e20;

    const double zb = (fl.type==0) ? -0.5*fl.h : 0.0;

    for(int q=0; q<nv; ++q)
    {
    fl.wx[q] = fl.x[0] + fl.R[0][0]*fl.bx[q] + fl.R[0][1]*fl.by[q] + fl.R[0][2]*zb;
    fl.wy[q] = fl.x[1] + fl.R[1][0]*fl.bx[q] + fl.R[1][1]*fl.by[q] + fl.R[1][2]*zb;

    fl.bbox[0] = MIN(fl.bbox[0],fl.wx[q]);
    fl.bbox[1] = MAX(fl.bbox[1],fl.wx[q]);
    fl.bbox[2] = MIN(fl.bbox[2],fl.wy[q]);
    fl.bbox[3] = MAX(fl.bbox[3],fl.wy[q]);
    }
}

double fnpf_ice::zbottom(const fnpf_ice_floe &fl, double xp, double yp) const
{
    const double nz = MAX(fl.nb[2],1.0e-6);
    return fl.pb[2] - (fl.nb[0]*(xp-fl.pb[0]) + fl.nb[1]*(yp-fl.pb[1]))/nz;
}

double fnpf_ice::zbottom_t(const fnpf_ice_floe &fl, double xp, double yp, double zp) const
{
    // Eulerian rate of the bottom plane height at fixed (x,y):
    // dz_b/dt = u_z - dz_b/dx*u_x - dz_b/dy*u_y,  u = v + w x r
    const double rx = xp-fl.x[0], ry = yp-fl.x[1], rz = zp-fl.x[2];
    const double ux = fl.v[0] + fl.w[1]*rz - fl.w[2]*ry;
    const double uy = fl.v[1] + fl.w[2]*rx - fl.w[0]*rz;
    const double uz = fl.v[2] + fl.w[0]*ry - fl.w[1]*rx;

    const double nz = MAX(fl.nb[2],1.0e-6);
    const double sx = -fl.nb[0]/nz;
    const double sy = -fl.nb[1]/nz;

    return uz - sx*ux - sy*uy;
}

double fnpf_ice::ramp(double d, double taper) const
{
    // smooth, area preserving edge taper: 0 at d = -taper, 1/2 on the edge, 1 at d = +taper
    if(taper<=0.0)
    return d>=0.0 ? 1.0 : 0.0;
    
    const double t = d/taper;
    
    if(t<=-1.0)
    return 0.0;
    
    if(t>=1.0)
    return 1.0;
    
    return 0.5 + 0.5*sin(0.5*PI*t);
}

double fnpf_ice::edge_distance(const fnpf_ice_floe &fl, double xp, double yp) const
{
    // signed distance to the convex footprint, positive inside (outside: exact away from corners)
    const int nv = int(fl.wx.size());
    double d = 1.0e20;
    
    for(int q=0; q<nv; ++q)
    {
    const int q2=(q+1)%nv;
    const double ex = fl.wx[q2]-fl.wx[q];
    const double ey = fl.wy[q2]-fl.wy[q];
    const double len = sqrt(ex*ex + ey*ey);
    
    if(len>0.0)
    d = MIN(d, (-ey*(xp-fl.wx[q]) + ex*(yp-fl.wy[q]))/len);
    }
    return d;
}

double fnpf_ice::coverage(lexer *p, const fnpf_ice_floe &fl, int ii, int jj) const
{
    const double taper = fl.tap;
    const double x0 = p->XN[ii+marge], x1 = p->XN[ii+1+marge];
    
    // 2D flume: taper along x
    if(is2D)
    {
    if(x1<=fl.bbox[0]-taper || x0>=fl.bbox[1]+taper)
    return 0.0;
    
    if(x0>=fl.bbox[0]+taper && x1<=fl.bbox[1]-taper)
    return 1.0;
    
    const int ns = 4*nsub;
    double sum=0.0;
    for(int a=0; a<ns; ++a)
    {
    const double xs = x0 + (a+0.5)*(x1-x0)/double(ns);
    sum += ramp(MIN(xs-fl.bbox[0], fl.bbox[1]-xs),taper);
    }
    return sum/double(ns);
    }
    
    const double y0 = p->YN[jj+marge], y1 = p->YN[jj+1+marge];
    
    if(x1<=fl.bbox[0]-taper || x0>=fl.bbox[1]+taper || y1<=fl.bbox[2]-taper || y0>=fl.bbox[3]+taper)
    return 0.0;
    
    // the distance is concave inside a convex polygon: its minimum over the cell sits on a corner
    if(MIN(MIN(edge_distance(fl,x0,y0),edge_distance(fl,x1,y0)),MIN(edge_distance(fl,x1,y1),edge_distance(fl,x0,y1))) >= taper)
    return 1.0;
    
    // product of the edge ramps: smooth (C1) also at the corners, where the distance function has a
    // crease; a 90 deg corner is rounded instead of pointed (area deficit is taken up by the rescaling)
    const int nv = int(fl.wx.size());
    double sum=0.0;
    for(int a=0; a<nsub; ++a)
    for(int b=0; b<nsub; ++b)
    {
    const double xs = x0 + (a+0.5)*(x1-x0)/double(nsub);
    const double ys = y0 + (b+0.5)*(y1-y0)/double(nsub);
    
    double prod=1.0;
    for(int q=0; q<nv && prod>0.0; ++q)
    {
    const int q2=(q+1)%nv;
    const double ex = fl.wx[q2]-fl.wx[q];
    const double ey = fl.wy[q2]-fl.wy[q];
    const double len = sqrt(ex*ex + ey*ey);
    if(len>0.0)
    prod *= ramp((-ey*(xs-fl.wx[q]) + ex*(ys-fl.wy[q]))/len, taper);
    }
    sum += prod;
    }
    
    return sum/double(nsub*nsub);
}

void fnpf_ice::footprint(lexer *p, fdm_fnpf *c)
{
    cell.clear();
    
    const double *xn = p->XN + marge;
    const double *yn = p->YN + marge;
    
    for(size_t f=0; f<floe.size(); ++f)
    {
        fnpf_ice_floe &fl = floe[f];
        
        if(fl.type!=0)
        continue;
        
        // taper at most a quarter of the floe width, small floes keep a sharper edge
        fl.tap = MIN(taper, 0.25*fl.wmin);
        const double tp = fl.tap;
        
        // local cell range of the bounding box
        const int is = MAX(0, int(upper_bound(xn, xn+p->knox+1, fl.bbox[0]-tp) - xn) - 1);
        const int ie = MIN(p->knox-1, int(lower_bound(xn, xn+p->knox+1, fl.bbox[1]+tp) - xn));
        
        int js=0, je=0;
        if(!is2D)
        {
        js = MAX(0, int(upper_bound(yn, yn+p->knoy+1, fl.bbox[2]-tp) - yn) - 1);
        je = MIN(p->knoy-1, int(lower_bound(yn, yn+p->knoy+1, fl.bbox[3]+tp) - yn));
        }
        
        if(is>ie || js>je)
        continue;
        
        for(i=is; i<=ie; ++i)
        for(j=js; j<=je; ++j)
        {
            PSLICECHECK4
            {
            const double phi = coverage(p,fl,i,j);
            
            if(phi>1.0e-12)
            {
            lidcell e;
            e.i = i;
            e.j = j;
            e.f = int(f);
            e.phi = phi;
            e.pl = 0.0;
            e.off = 0.0;
            e.area = is2D ? p->DXN[IP]*fl.width2D : p->DXN[IP]*p->DYN[JP];
            
            e.xc = p->XP[IP];
            e.yc = is2D ? fl.x[1] : p->YP[JP];
            cell.push_back(e);
            }
            }
        }
    }
    
    // per floe: sum phi*A, sum phi*A*r^2 over all ranks
    const int nv = 2;
    vector<double> sum(nv*floe.size(),0.0);
    
    for(auto &e : cell)
    {
    const fnpf_ice_floe &fl = floe[e.f];
    const double A = e.phi*e.area;
    const double rx = e.xc-fl.x[0], ry = e.yc-fl.x[1];
    sum[nv*e.f]   += A;
    sum[nv*e.f+1] += A*(rx*rx + ry*ry);
    }
    
    if(p->mpi_size>1 && !sum.empty())
    MPI_Allreduce(MPI_IN_PLACE, sum.data(), int(sum.size()), MPI_DOUBLE, MPI_SUM, comm);
    
    // The taper preserves the area of straight edges, not at corners: rescale phi so the footprint
    // carries the floe area exactly (weight = buoyancy at rest). Floes cut by the domain boundary or
    // over dry cells are left as they are.
    vector<double> scale(floe.size(),1.0);
    
    for(size_t f=0; f<floe.size(); ++f)
    {
        fnpf_ice_floe &fl = floe[f];
        
        if(fl.type!=0)
        continue;
        
        const double S0 = sum[nv*f];
        const double tp = fl.tap;
        
        const int inside = fl.bbox[0]-tp>p->global_xmin && fl.bbox[1]+tp<p->global_xmax
                        && (is2D || (fl.bbox[2]-tp>p->global_ymin && fl.bbox[3]+tp<p->global_ymax));
        
        if(inside && S0>0.0)
        scale[f] = MAX(0.5, MIN(2.0, fl.area/S0));
        
        // lid-spring frequencies of the footprint, for the time step
        const double Kz = fl.klid*scale[f]*S0;
        const double Kr = fl.klid*scale[f]*sum[nv*f+1];
        const double Imin = MAX(is2D ? fl.Ib[1][1] : MIN(fl.Ib[0][0],fl.Ib[1][1]), 1.0e-20);
        fl.omega = MAX(sqrt(Kz/MAX(fl.mass,1.0e-20)), sqrt(Kr/Imin));
    }
    
    for(auto &e : cell)
    e.phi *= scale[e.f];
    
    // equilibrium depression of the surface, sum over the floes sharing a cell
    SLICELOOP4
    dtot(i,j) = 0.0;
    
    for(auto &e : cell)
    {
    i=e.i;
    j=e.j;
    dtot(i,j) += e.phi*floe[e.f].rho*floe[e.f].h/rhow;
    }
    
    for(auto &e : cell)
    {
    i=e.i;
    j=e.j;
    e.off = floe[e.f].rho*floe[e.f].h/rhow - dtot(i,j);
    }
}

void fnpf_ice::stage_forces(lexer *p, fdm_fnpf *c, slice &K, slice &eta)
{
    // per floe: F[3], M[3], Awet
    const int nv = 7;
    vector<double> sum(nv*floe.size(),0.0);
    
    for(auto &e : cell)
    {
        i=e.i;
        j=e.j;
        
        WETDRY
        {
        const fnpf_ice_floe &fl = floe[e.f];
        double *s = &sum[nv*e.f];
        
        // bottom plane of the stage floe state at the cell centre; in the edge taper and in cells shared
        // by several floes the bottom is shifted by draft_f - sum_g phi_g*draft_g, so every floe is in
        // equilibrium with the depressed surface (single floe: (1-phi)*draft)
        const double zb  = zbottom(fl,e.xc,e.yc) + e.off;
        const double zbt = zbottom_t(fl,e.xc,e.yc,zb);
        
        const double arg = fl.pw + fl.klid*(wd + eta(i,j) - zb) + fl.clid*(etat(i,j) - zbt);
        
        e.pl = MAX(arg,0.0);
        
        if(arg<=0.0)
        continue;
        
        // lid pressure into the dynamic FSBC: dFi/dt = ... - p/rho_w
        K(i,j) -= e.phi*arg/rhow;
        
        const double A = e.phi*e.area;
        const double pA = arg*A;
        
        const double rx = e.xc - fl.x[0];
        const double ry = e.yc - fl.x[1];
        const double rz = zb - fl.x[2];
        
        // Lid pressure on the loaded surface z = eta: dF = p*A_h*(-deta/dx, -deta/dy, 1).
        // This is exactly the reaction of the pressure given to the fluid, so horizontal momentum is
        // conserved. Under the rigid part eta follows the floe bottom (slope force); in the edge taper
        // the surface gradient carries the side (waterline) pressure of the floe, which gives the
        // relative-elevation part of the mean wave drift force. The bottom tilt alone misses it and
        // lets floes drift up-wave in reflected wave fields.
        const double ex = (eta(i+1,j)-eta(i-1,j))/(p->XP[IP1]-p->XP[IM1]);
        const double ey = is2D ? 0.0 : (eta(i,j+1)-eta(i,j-1))/(p->YP[JP1]-p->YP[JM1]);
        double fx = -pA*ex;
        double fy = -pA*ey;
        double fz = pA;
        
        // ice-water skin drag on the relative horizontal velocity at the surface
        if(Cd>0.0)
        {
        const double us = c->Fx(i,j) - c->Ex(i,j)*c->Fz(i,j);
        const double vs = is2D ? 0.0 : c->Fy(i,j) - c->Ey(i,j)*c->Fz(i,j);
        
        const double upx = fl.v[0] + fl.w[1]*rz - fl.w[2]*ry;
        const double upy = fl.v[1] + fl.w[2]*rx - fl.w[0]*rz;
        
        const double du = us - upx;
        const double dv = vs - upy;
        const double mag = sqrt(du*du + dv*dv);
        
        fx += rhow*Cd*A*mag*du;
        fy += rhow*Cd*A*mag*dv;
        }
        
        s[0] += fx;
        s[1] += fy;
        s[2] += fz;
        s[3] += ry*fz - rz*fy;
        s[4] += rz*fx - rx*fz;
        s[5] += rx*fy - ry*fx;
        s[6] += A;
        }
    }
    
    if(p->mpi_size>1 && !sum.empty())
    MPI_Allreduce(MPI_IN_PLACE, sum.data(), int(sum.size()), MPI_DOUBLE, MPI_SUM, comm);
    
    for(size_t f=0; f<floe.size(); ++f)
    {
        fnpf_ice_floe &fl = floe[f];
        const double *s = &sum[nv*f];
        
        for(int q=0;q<3;++q)
        {
        fl.F[q] = s[q];
        fl.M[q] = s[3+q];
        }
        fl.Awet = s[6];
    }
}

void fnpf_ice::quat_to_matrix(const double *q, double (*R)[3])
{
    const double w=q[0], x=q[1], y=q[2], z=q[3];

    R[0][0] = 1.0 - 2.0*(y*y + z*z);
    R[0][1] = 2.0*(x*y - w*z);
    R[0][2] = 2.0*(x*z + w*y);
    R[1][0] = 2.0*(x*y + w*z);
    R[1][1] = 1.0 - 2.0*(x*x + z*z);
    R[1][2] = 2.0*(y*z - w*x);
    R[2][0] = 2.0*(x*z - w*y);
    R[2][1] = 2.0*(y*z + w*x);
    R[2][2] = 1.0 - 2.0*(x*x + y*y);
}

void fnpf_ice::quat_to_euler(const double *q, double &roll, double &pitch, double &yaw)
{
    const double w=q[0], x=q[1], y=q[2], z=q[3];

    roll  = atan2(2.0*(w*x + y*z), 1.0 - 2.0*(x*x + y*y));
    const double sp = 2.0*(w*y - z*x);
    pitch = fabs(sp)>=1.0 ? copysign(0.5*PI,sp) : asin(sp);
    yaw   = atan2(2.0*(w*z + x*y), 1.0 - 2.0*(y*y + z*z));
}

void fnpf_ice::quat_normalize(double *q)
{
    const double n = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    for(int a=0;a<4;++a)
    q[a]/=n;
}
