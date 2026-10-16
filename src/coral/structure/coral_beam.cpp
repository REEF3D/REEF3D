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

#include"coral/structure/coral_beam.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include<cmath>
#include<stdexcept>

namespace coral
{

beam_section beam_section::tube(double D, double t, double E, double G, double rho)
{
    beam_section s;
    const double Di = D - 2.0*t;
    const double A = 0.25*M_PI*(D*D - Di*Di);
    const double I = M_PI/64.0*(D*D*D*D - Di*Di*Di*Di);
    s.EA = E*A;
    s.EIy = s.EIz = E*I;
    s.GJ = G*2.0*I;
    s.GAy = s.GAz = 0.5*G*A;
    s.m = rho*A;
    s.rIx = rho*2.0*I;
    s.rIy = s.rIz = rho*I;
    s.D = D;
    return s;
}

beam_section beam_section::blend(const beam_section &a, const beam_section &b, double s)
{
    beam_section r;
    auto mix = [s](double x, double y) {return (1.0 - s)*x + s*y;};
    r.EA = mix(a.EA, b.EA); r.EIy = mix(a.EIy, b.EIy); r.EIz = mix(a.EIz, b.EIz); r.GJ = mix(a.GJ, b.GJ);
    r.GAy = mix(a.GAy, b.GAy); r.GAz = mix(a.GAz, b.GAz); r.m = mix(a.m, b.m);
    r.rIx = mix(a.rIx, b.rIx); r.rIy = mix(a.rIy, b.rIy); r.rIz = mix(a.rIz, b.rIz); r.D = mix(a.D, b.D);
    return r;
}

beam_element::beam_element(const model &m, const beam_end &a, const beam_end &b, const beam_section &s,
                           const Vec3 &orient, double beta, const std::string &name)
    : name_(name), ba_(a.body), bb_(b.body), s_(s), beta_(beta)
{
    const Vec3 d = b.p - a.p;
    L0_ = d.norm();
    if(L0_ <= 0.0)
        throw std::runtime_error("coral::beam_element " + name + ": zero length");
    const Vec3 ex = d/L0_;
    Vec3 ey = orient - orient.dot(ex)*ex;
    if(ey.norm() < 1.0e-6)
    {
        const Vec3 alt = std::fabs(ex(2)) < 0.9 ? Vec3::UnitZ() : Vec3::UnitX();
        ey = alt - alt.dot(ex)*ex;
    }
    ey.normalize();
    Mat3 E0;
    E0.col(0) = ex;
    E0.col(1) = ey;
    E0.col(2) = ex.cross(ey);
    const rigid_body &A = m.body(a.body), &B = m.body(b.body);
    sa_ = A.to_body(a.p);
    sb_ = B.to_body(b.p);
    Qa_ = A.R.transpose()*E0;
    Qb_ = B.R.transpose()*E0;
}

void beam_element::ends(const model &m, Vec3 &pa, Mat3 &Ta, Vec3 &pb, Mat3 &Tb) const
{
    const rigid_body &A = m.body(ba_), &B = m.body(bb_);
    pa = A.point(sa_);
    pb = B.point(sb_);
    Ta = A.R*Qa_;
    Tb = B.R*Qb_;
}

static Mat3 element_frame(const Vec3 &pa, const Mat3 &Ta, const Vec3 &pb, const Mat3 &Tb, double &L)
{
    const Vec3 d = pb - pa;
    L = d.norm();
    const Vec3 e1 = d/L;
    const Vec3 q = 0.5*(Ta.col(1) + Tb.col(1));
    Vec3 e3 = e1.cross(q);
    e3.normalize();
    Mat3 E;
    E.col(0) = e1;
    E.col(1) = e3.cross(e1);
    E.col(2) = e3;
    return E;
}

Eigen::Matrix<double,7,1> beam_element::deformation(const Vec3 &pa, const Mat3 &Ta, const Vec3 &pb, const Mat3 &Tb) const
{
    double L;
    const Mat3 E = element_frame(pa, Ta, pb, Tb, L);
    Eigen::Matrix<double,7,1> d;
    d(0) = L - L0_;
    d.segment<3>(1) = log_so3(E.transpose()*Ta);
    d.segment<3>(4) = log_so3(E.transpose()*Tb);
    return d;
}

Eigen::Matrix<double,7,1> beam_element::deformation(const model &m) const
{
    Vec3 pa, pb; Mat3 Ta, Tb;
    ends(m, pa, Ta, pb, Tb);
    return deformation(pa, Ta, pb, Tb);
}

Mat3 beam_element::frame(const model &m, Vec3 &pa, Vec3 &pb) const
{
    Mat3 Ta, Tb;
    ends(m, pa, Ta, pb, Tb);
    double L;
    return element_frame(pa, Ta, pb, Tb, L);
}

Eigen::Matrix<double,7,1> beam_element::local_forces(const Eigen::Matrix<double,7,1> &d) const
{
    const double L = L0_;
    Eigen::Matrix<double,7,1> f;
    f(0) = s_.EA*d(0)/L;
    const double T = s_.GJ*(d(4) - d(1))/L;
    f(1) = -T;
    f(4) = T;
    // bending about local y (rotations ty) and z (rotations tz)
    const double Py = s_.GAz > 0.0 ? 12.0*s_.EIy/(s_.GAz*L*L) : 0.0;
    const double Pz = s_.GAy > 0.0 ? 12.0*s_.EIz/(s_.GAy*L*L) : 0.0;
    const double ky = s_.EIy/(L*(1.0 + Py)), kz = s_.EIz/(L*(1.0 + Pz));
    f(2) = ky*((4.0 + Py)*d(2) + (2.0 - Py)*d(5));
    f(5) = ky*((2.0 - Py)*d(2) + (4.0 + Py)*d(5));
    f(3) = kz*((4.0 + Pz)*d(3) + (2.0 - Pz)*d(6));
    f(6) = kz*((2.0 - Pz)*d(3) + (4.0 + Pz)*d(6));
    return f;
}

double beam_element::potential(const model &m) const
{
    const Eigen::Matrix<double,7,1> d = deformation(m);
    return 0.5*d.dot(local_forces(d));
}

void beam_element::eval(const model &m, double, VecX &f) const
{
    Vec3 pa, pb; Mat3 Ta, Tb;
    ends(m, pa, Ta, pb, Tb);
    const Eigen::Matrix<double,7,1> d = deformation(pa, Ta, pb, Tb);
    // B = dd/dq by central differences: q = (dx_a, dth_a, dx_b, dth_b), rotations in the body frames
    const rigid_body &A = m.body(ba_), &B = m.body(bb_);
    Eigen::Matrix<double,7,12> Bm = Eigen::Matrix<double,7,12>::Zero();
    const double ex = 1.0e-6*L0_, er = 1.0e-6;
    for(int k=0; k<12; ++k)
    {
        const bool onA = k < 6;
        const int body = onA ? ba_ : bb_;
        if(body < 0)
            continue;
        const rigid_body &Bd = onA ? A : B;
        const int c = k%6;
        Eigen::Matrix<double,7,1> dp, dm;
        for(int sgn=-1; sgn<=1; sgn+=2)
        {
            Vec3 x = Bd.x;
            Mat3 R = Bd.R;
            if(c < 3)
                x(c) += sgn*ex;
            else
                R = R*exp_so3(sgn*er*Vec3::Unit(c-3));
            const Vec3 p = x + R*(onA ? sa_ : sb_);
            const Mat3 T = R*(onA ? Qa_ : Qb_);
            const Eigen::Matrix<double,7,1> dd = onA ? deformation(p, T, pb, Tb) : deformation(pa, Ta, p, T);
            if(sgn < 0) dm = dd; else dp = dd;
        }
        Bm.col(k) = (dp - dm)/(2.0*(c < 3 ? ex : er));
    }
    Eigen::Matrix<double,7,1> fl = local_forces(d);
    if(beta_ > 0.0)
    {
        Eigen::Matrix<double,12,1> qd = Eigen::Matrix<double,12,1>::Zero();
        if(ba_ >= 0) {qd.segment<3>(0) = A.v; qd.segment<3>(3) = A.w;}
        if(bb_ >= 0) {qd.segment<3>(6) = B.v; qd.segment<3>(9) = B.w;}
        fl += beta_*local_forces(Bm*qd);
    }
    const Eigen::Matrix<double,12,1> g = -Bm.transpose()*fl;
    if(ba_ >= 0)
    {
        f.segment<3>(6*ba_) += g.segment<3>(0);
        f.segment<3>(6*ba_+3) += g.segment<3>(3);
    }
    if(bb_ >= 0)
    {
        f.segment<3>(6*bb_) += g.segment<3>(6);
        f.segment<3>(6*bb_+3) += g.segment<3>(9);
    }
}

}
