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

// External dynamics of a 6DOF body (6DOF_external.h): the state comes from the external solver, the
// loads of the fluid go to it.

#include"6DOF_obj.h"
#include"lexer.h"
#include<algorithm>
#include<cmath>

void sixdof_obj::external_set_state(const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Vector3d &v, const Eigen::Vector3d &w)
{
    // CoG, quaternion (R = E G^T, e = (e0, e1, e2, e3), sign kept continuous), momentum and the body
    // frame angular momentum; the velocities themselves are kept for update_fbvel
    Eigen::Quaterniond q(R);
    q.normalize();
    Eigen::Vector4d e(q.w(), q.x(), q.y(), q.z());
    
    if(e.dot(rb.e)<0.0)
    e = -e;
    
    rb.c = c;
    rb.e = e;
    rb.quat_matrices();
    rb.p = rb.mass*v;
    rb.h = rb.I*(rb.R.transpose()*w);
    rb.dc = v;
    
    ext_v_ = v;
    ext_w_ = w;
    ext_on_ = true;
}

void sixdof_obj::external_fluid_load(lexer *p, Eigen::Matrix<double,6,1> &F) const
{
    // the fluid load of the last force evaluation without the weight (Xe .. Ne include gravity)
    F(0) = Xe - p->W20*Mass_fb;
    F(1) = Ye - p->W21*Mass_fb;
    F(2) = Ze - p->W22*Mass_fb;
    F(3) = Ke;
    F(4) = Me;
    F(5) = Ne;
}

void sixdof_obj::external_properties(double &m, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Matrix3d &I) const
{
    m = Mass_fb;
    c = c_;
    R = R_;
    I = I_;
}

void sixdof_obj::external_strips_set(const Eigen::Vector3d &a, const Eigen::Vector3d &e, double L, int n)
{
    strip_a_ = a;
    strip_e_ = e.normalized();
    strip_L_ = L>0.0 ? L : 1.0;
    strip_n_ = n>0 ? n : 0;
    strip_F_.assign(3*strip_n_,0.0);
}

void sixdof_obj::external_strips_zero()
{
    std::fill(strip_F_.begin(),strip_F_.end(),0.0);
}

void sixdof_obj::external_strips_tri(const double a[3], const double b[3], const double c[3], const double pm[3],
                                     double nx, double ny, double nz)
{
    if(strip_n_<=0)
    return;
    
    // polygon vertices: position, pressure (linear through the edge-midpoint values), strip coordinate
    struct vtx {double x[3]; double p; double u;};
    const double scale = strip_n_/strip_L_;
    auto make = [&](const double *x, double pv)
    {
        vtx v;
        for(int i=0; i<3; ++i)
        v.x[i] = x[i];
        v.p = pv;
        v.u = ((x[0] - strip_a_(0))*strip_e_(0) + (x[1] - strip_a_(1))*strip_e_(1) + (x[2] - strip_a_(2))*strip_e_(2))*scale;
        return v;
    };
    const vtx tri[3] = {make(a, pm[0] + pm[2] - pm[1]), make(b, pm[0] + pm[1] - pm[2]), make(c, pm[1] + pm[2] - pm[0])};
    
    auto strip_of = [&](double u) {int k = (int)std::floor(u); return k<0 ? 0 : (k>=strip_n_ ? strip_n_-1 : k);};
    const int k0 = strip_of(std::min(tri[0].u, std::min(tri[1].u, tri[2].u)));
    const int k1 = strip_of(std::max(tri[0].u, std::max(tri[1].u, tri[2].u)));
    
    // integral of p over a convex polygon (fan of triangles, exact for linear p)
    auto integral = [](const vtx *v, int nv)
    {
        double I = 0.0;
        for(int s=1; s<nv-1; ++s)
        {
            const double e1[3] = {v[s].x[0] - v[0].x[0], v[s].x[1] - v[0].x[1], v[s].x[2] - v[0].x[2]};
            const double e2[3] = {v[s+1].x[0] - v[0].x[0], v[s+1].x[1] - v[0].x[1], v[s+1].x[2] - v[0].x[2]};
            const double cx = e1[1]*e2[2] - e1[2]*e2[1], cy = e1[2]*e2[0] - e1[0]*e2[2], cz = e1[0]*e2[1] - e1[1]*e2[0];
            I += 0.5*std::sqrt(cx*cx + cy*cy + cz*cz)*(v[0].p + v[s].p + v[s+1].p)/3.0;
        }
        return I;
    };
    // keep the part with sgn*(u - lim) >= 0 (Sutherland-Hodgman on the strip coordinate)
    auto clip = [](const vtx *in, int nin, vtx *out, double lim, double sgn)
    {
        int nout = 0;
        for(int q=0; q<nin; ++q)
        {
            const vtx &P = in[q], &Q = in[(q+1)%nin];
            const bool inP = sgn*(P.u - lim) >= 0.0, inQ = sgn*(Q.u - lim) >= 0.0;
            if(inP)
            out[nout++] = P;
            if(inP != inQ)
            {
                const double t = (lim - P.u)/(Q.u - P.u);
                vtx &R = out[nout++];
                for(int i=0; i<3; ++i)
                R.x[i] = P.x[i] + t*(Q.x[i] - P.x[i]);
                R.p = P.p + t*(Q.p - P.p);
                R.u = lim;
            }
        }
        return nout;
    };
    
    for(int k=k0; k<=k1; ++k)
    {
        double I;
        if(k0==k1)
        I = integral(tri, 3);
        else
        {
            vtx v1[5], v2[6];
            int n1 = 3;
            const vtx *cur = tri;
            if(k>0)
            {
                n1 = clip(tri, 3, v1, (double)k, 1.0);
                cur = v1;
            }
            int n2 = n1;
            if(k<strip_n_-1)
            {
                n2 = clip(cur, n1, v2, k + 1.0, -1.0);
                cur = v2;
            }
            if(n2<3)
            continue;
            I = integral(cur, n2);
        }
        strip_F_[3*k]   -= I*nx;
        strip_F_[3*k+1] -= I*ny;
        strip_F_[3*k+2] -= I*nz;
    }
}
