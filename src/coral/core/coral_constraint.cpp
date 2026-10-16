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

#include"coral/core/coral_constraint.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include<cmath>
#include<stdexcept>

namespace coral
{

joint::joint(const model &m, joint_type type, int a, int b, const Vec3 &P, const Vec3 &e_in, const std::string &name)
    : name_(name), type_(type)
{
    if(a==b)
        throw std::runtime_error("coral::joint: both ends on the same body");
    body_a = a;
    body_b = b;
    const rigid_body &A = m.body(a);
    const rigid_body &B = m.body(b);

    Vec3 e = e_in;
    if(e.norm() < 1.0e-12)
        e = Vec3::UnitZ();
    e.normalize();
    // two normals to the axis
    Vec3 n1 = std::fabs(e(0)) < 0.9 ? e.cross(Vec3::UnitX()) : e.cross(Vec3::UnitY());
    n1.normalize();
    const Vec3 n2 = e.cross(n1);

    sa_ = A.R.transpose()*(P - A.x);
    sb_ = B.R.transpose()*(P - B.x);
    ea_ = A.R.transpose()*e;
    eb_ = B.R.transpose()*e;
    pa_ = A.R.transpose()*n1;
    pb_ = B.R.transpose()*n1;
    const Vec3 n2a = A.R.transpose()*n2, n2b = B.R.transpose()*n2;
    const Vec3 n1b = pb_;
    const Vec3 n1a = pa_;

    auto lock_orientation = [&]()
    {
        const Vec3 e1a = A.R.transpose()*Vec3::UnitX(), e2a = A.R.transpose()*Vec3::UnitY();
        const Vec3 e2b = B.R.transpose()*Vec3::UnitY(), e3b = B.R.transpose()*Vec3::UnitZ();
        add_perp(e1a, e2b);
        add_perp(e1a, e3b);
        add_perp(e2a, e3b);
    };

    switch(type)
    {
        case joint_type::fixed:
            npoint_ = 1;
            lock_orientation();
            break;
        case joint_type::spherical:
            npoint_ = 1;
            break;
        case joint_type::revolute:
            npoint_ = 1;
            add_perp(ea_, n1b);
            add_perp(ea_, n2b);
            break;
        case joint_type::prismatic:
            lock_orientation();
            add_dist(n1a);
            add_dist(n2a);
            break;
        case joint_type::cylindrical:
            add_perp(ea_, n1b);
            add_perp(ea_, n2b);
            add_dist(n1a);
            add_dist(n2a);
            break;
        case joint_type::planar:
            add_dist(ea_);
            add_perp(ea_, n1b);
            add_perp(ea_, n2b);
            break;
    }
    (void)n2a;
    d0_ = 0.0;    // the joint point is the same on both bodies at creation
}

void joint::eval(const model &m, Eigen::Ref<VecX> g, Eigen::Ref<MatX> G) const
{
    const rigid_body &A = m.body(body_a);
    const rigid_body &B = m.body(body_b);
    const int ca = 6*body_a, cb = 6*body_b;
    const bool ha = body_a>=0, hb = body_b>=0;

    const Vec3 pA = A.x + A.R*sa_;
    const Vec3 pB = B.x + B.R*sb_;
    const Vec3 d = pB - pA;
    const Mat3 RaSa = A.R*skew(sa_);
    const Mat3 RbSb = B.R*skew(sb_);

    int row = 0;
    if(npoint_)
    {
        g.segment<3>(0) = pA - pB;
        if(ha)
        {
            G.block<3,3>(0,ca) = Mat3::Identity();
            G.block<3,3>(0,ca+3) = -RaSa;
        }
        if(hb)
        {
            G.block<3,3>(0,cb) = -Mat3::Identity();
            G.block<3,3>(0,cb+3) = RbSb;
        }
        row = 3;
    }
    for(const prim &p : prims_)
    {
        const Vec3 ug = A.R*p.u;
        if(p.kind==0)
        {
            const Vec3 wg = B.R*p.w;
            g(row) = ug.dot(wg);
            if(ha) G.block<1,3>(row,ca+3) = -wg.transpose()*A.R*skew(p.u);
            if(hb) G.block<1,3>(row,cb+3) = -ug.transpose()*B.R*skew(p.w);
        }
        else
        {
            g(row) = d.dot(ug);
            if(ha)
            {
                G.block<1,3>(row,ca) = -ug.transpose();
                G.block<1,3>(row,ca+3) = ug.transpose()*RaSa - d.transpose()*A.R*skew(p.u);
            }
            if(hb)
            {
                G.block<1,3>(row,cb) = ug.transpose();
                G.block<1,3>(row,cb+3) = -ug.transpose()*RbSb;
            }
        }
        ++row;
    }
}

Vec3 joint::axis(const model &m) const
{
    return m.body(body_a).R*ea_;
}

double joint::coordinate(const model &m) const
{
    const rigid_body &A = m.body(body_a);
    const rigid_body &B = m.body(body_b);
    return (B.point(sb_) - A.point(sa_)).dot(A.R*ea_) - d0_;
}

double joint::rate(const model &m) const
{
    const rigid_body &A = m.body(body_a);
    const rigid_body &B = m.body(body_b);
    const Vec3 d = B.point(sb_) - A.point(sa_);
    const Vec3 dd = B.point_velocity(sb_) - A.point_velocity(sa_);
    return dd.dot(A.R*ea_) + d.dot(A.R*A.w.cross(ea_));
}

double joint::rotation_coordinate(const model &m) const
{
    const rigid_body &A = m.body(body_a);
    const rigid_body &B = m.body(body_b);
    const Vec3 u = A.R*pa_, w = B.R*pb_, e = A.R*ea_;
    return std::atan2(u.cross(w).dot(e), u.dot(w));
}

double joint::rotation_rate(const model &m) const
{
    const rigid_body &A = m.body(body_a);
    const rigid_body &B = m.body(body_b);
    return (B.omega_global() - A.omega_global()).dot(A.R*ea_);
}

void joint::add_axial_force(const model &m, double F, VecX &f) const
{
    const Vec3 e = axis(m);
    const Vec3 p = m.body(body_b).point(sb_);
    m.add_force(f, body_b, p, F*e);
    m.add_force(f, body_a, p, -F*e);
}

void joint::add_axial_torque(const model &m, double T, VecX &f) const
{
    const Vec3 e = axis(m);
    m.add_moment(f, body_b, T*e);
    m.add_moment(f, body_a, -T*e);
}

}
