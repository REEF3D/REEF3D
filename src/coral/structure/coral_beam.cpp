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

namespace
{
    // inverse of the left Jacobian of SO(3): exp(psi + Jl^-1(psi) eta) = exp(eta) exp(psi) to first order
    Mat3 jl_inv(const Vec3 &psi)
    {
        const double th = psi.norm();
        const Mat3 S = skew(psi);
        double c;
        if(th < 1.0e-4)
            c = 1.0/12.0 + th*th/720.0;
        else
            c = 1.0/(th*th) - (1.0 + std::cos(th))/(2.0*th*std::sin(th));
        return Mat3::Identity() - 0.5*S + c*S*S;
    }
}

void beam_element::kinematics(const Vec3 &xa, const Mat3 &Ra, const Vec3 &xb, const Mat3 &Rb, Vec7 &d, Mat7x12 &B) const
{
    const Vec3 pa = xa + Ra*sa_, pb = xb + Rb*sb_;
    const Mat3 Ta = Ra*Qa_, Tb = Rb*Qb_;
    double L;
    const Mat3 E = element_frame(pa, Ta, pb, Tb, L);
    const Vec3 psa = log_so3(E.transpose()*Ta), psb = log_so3(E.transpose()*Tb);
    d(0) = L - L0_;
    d.segment<3>(1) = psa;
    d.segment<3>(4) = psb;
    // variations of the chord dd, the end spins wa, wb (global) in the 12 coordinates
    Eigen::Matrix<double,3,12> Dd = Eigen::Matrix<double,3,12>::Zero(), Wa = Dd, Wb = Dd;
    Dd.block<3,3>(0,0) = -Mat3::Identity();
    Dd.block<3,3>(0,3) = Ra*skew(sa_);
    Dd.block<3,3>(0,6) = Mat3::Identity();
    Dd.block<3,3>(0,9) = -Rb*skew(sb_);
    Wa.block<3,3>(0,3) = Ra;
    Wb.block<3,3>(0,9) = Rb;
    const Vec3 e1 = E.col(0), e2 = E.col(1), e3 = E.col(2);
    const Vec3 ya = Ta.col(1), yb = Tb.col(1), q = 0.5*(ya + yb);
    const double nn = e1.cross(q).norm();
    const Mat3 P = Mat3::Identity() - e1*e1.transpose();
    // spin of the element frame: e3 component from the chord in e2, e2 component from the chord in e3,
    // e1 component from the rotation of e3 = e1 x q / |e1 x q| about e1
    const Eigen::Matrix<double,1,12> w3 = e2.transpose()*Dd/L;
    const Eigen::Matrix<double,1,12> w2 = -e3.transpose()*Dd/L;
    // (d e3 = w x e3, so w.e1 = -e2.d e3 = -e2.d n/|n|, d n = d e1 x q + e1 x d q)
    const Eigen::Matrix<double,1,12> w1 = -((P*q.cross(e2)).transpose()*Dd/L - 0.5*ya.cross(e3).transpose()*Wa
                                           - 0.5*yb.cross(e3).transpose()*Wb)/nn;
    const Eigen::Matrix<double,3,12> WE = e1*w1 + e2*w2 + e3*w3;
    B.row(0) = e1.transpose()*Dd;
    B.middleRows<3>(1) = jl_inv(psa)*E.transpose()*(Wa - WE);
    B.middleRows<3>(4) = jl_inv(psb)*E.transpose()*(Wb - WE);
    if(ba_ < 0) B.leftCols<6>().setZero();
    if(bb_ < 0) B.rightCols<6>().setZero();
}

Eigen::Matrix<double,7,7> beam_element::local_stiffness() const
{
    Eigen::Matrix<double,7,7> K;
    for(int k=0; k<7; ++k)
        K.col(k) = local_forces(Vec7::Unit(k));
    return K;
}

beam_element::Vec12 beam_element::element_force(const Vec3 &xa, const Mat3 &Ra, const Vec3 &xb, const Mat3 &Rb, const Vec12 &qd) const
{
    Vec7 d;
    Mat7x12 B;
    kinematics(xa, Ra, xb, Rb, d, B);
    Vec7 fl = local_forces(d);
    if(beta_ > 0.0)
        fl += beta_*local_forces(B*qd);
    return -B.transpose()*fl;
}

beam_element::Vec12 beam_element::rates(const model &m) const
{
    Vec12 qd = Vec12::Zero();
    if(ba_ >= 0) {qd.segment<3>(0) = m.body(ba_).v; qd.segment<3>(3) = m.body(ba_).w;}
    if(bb_ >= 0) {qd.segment<3>(6) = m.body(bb_).v; qd.segment<3>(9) = m.body(bb_).w;}
    return qd;
}

Eigen::Matrix<double,7,12> beam_element::strain_matrix(const model &m) const
{
    const rigid_body &A = m.body(ba_), &B = m.body(bb_);
    Vec7 d;
    Mat7x12 Bm;
    kinematics(A.x, A.R, B.x, B.R, d, Bm);
    return Bm;
}

void beam_element::eval(const model &m, double, VecX &f) const
{
    const rigid_body &A = m.body(ba_), &B = m.body(bb_);
    const Vec12 g = element_force(A.x, A.R, B.x, B.R, rates(m));
    if(ba_ >= 0)
        f.segment<6>(6*ba_) += g.head<6>();
    if(bb_ >= 0)
        f.segment<6>(6*bb_) += g.tail<6>();
}

void beam_element::tangent(const model &m, double, double cK, double cC, double, std::vector<triplet> &T) const
{
    const rigid_body &A = m.body(ba_), &Bd = m.body(bb_);
    const Vec12 qd = rates(m);
    Eigen::Matrix<double,12,12> K = Eigen::Matrix<double,12,12>::Zero();
    if(cK != 0.0)
    {
        // dF/dq by central differences of the element force (analytic B), velocities fixed
        const double ex = 1.0e-5*L0_, er = 1.0e-5;
        for(int k=0; k<12; ++k)
        {
            const bool onA = k < 6;
            if((onA ? ba_ : bb_) < 0)
                continue;
            const int c = k%6;
            Vec12 gp, gm;
            for(int sgn=-1; sgn<=1; sgn+=2)
            {
                Vec3 xa = A.x, xb = Bd.x;
                Mat3 Ra = A.R, Rb = Bd.R;
                Vec3 &x = onA ? xa : xb;
                Mat3 &R = onA ? Ra : Rb;
                if(c < 3) x(c) += sgn*ex;
                else R = R*exp_so3(sgn*er*Vec3::Unit(c-3));
                (sgn < 0 ? gm : gp) = element_force(xa, Ra, xb, Rb, qd);
            }
            K.col(k) = cK*(gp - gm)/(2.0*(c < 3 ? ex : er));
        }
    }
    if(cC != 0.0 && beta_ > 0.0)
    {
        Vec7 d;
        Mat7x12 B;
        kinematics(A.x, A.R, Bd.x, Bd.R, d, B);
        K -= cC*beta_*B.transpose()*local_stiffness()*B;
    }
    auto dof = [&](int k) {return k < 6 ? (ba_ >= 0 ? 6*ba_ + k : -1) : (bb_ >= 0 ? 6*bb_ + k - 6 : -1);};
    for(int j=0; j<12; ++j)
    {
        const int cj = dof(j);
        if(cj < 0) continue;
        for(int i=0; i<12; ++i)
        {
            const int ri = dof(i);
            if(ri >= 0 && K(i, j) != 0.0)
                T.emplace_back(ri, cj, K(i, j));
        }
    }
}

}
