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

#include"coral/structure/coral_frame.h"
#include"coral/core/coral_model.h"
#include<cmath>
#include<memory>
#include<stdexcept>

namespace coral
{

int frame_builder::node_at(const Vec3 &p, const std::string &name)
{
    for(size_t k=0; k<nodes_.size(); ++k)
        if((m_.bodies[nodes_[k]].x - p).norm() <= tol_*(1.0 + p.norm()))
            return nodes_[k];
    rigid_body b;
    b.name = name;
    b.x = p;
    const int i = m_.add_body(b);
    index_[i] = (int)nodes_.size();
    nodes_.push_back(i);
    mass_.push_back(0.0);
    inertia_.push_back(Mat3::Zero());
    return i;
}

std::vector<int> frame_builder::add_member(const member_spec &s)
{
    if(finished_)
        throw std::runtime_error("coral::frame_builder: member after finish()");
    const int n = std::max(1, s.n);
    std::vector<int> st(n+1);
    for(int i=0; i<=n; ++i)
    {
        const Vec3 p = s.p1 + (double)i/n*(s.p2 - s.p1);
        if(i==0 && s.body1 > -2) {st[i] = s.body1; continue;}
        if(i==n && s.body2 > -2) {st[i] = s.body2; continue;}
        st[i] = node_at(p, s.name + ":" + std::to_string(i));
    }
    const Vec3 dir = (s.p2 - s.p1).normalized();
    for(int e=0; e<n; ++e)
    {
        const Vec3 pa = s.p1 + (double)e/n*(s.p2 - s.p1), pb = s.p1 + (double)(e+1)/n*(s.p2 - s.p1);
        const beam_section sec = s.tapered ? beam_section::blend(s.s1, s.s2, (e + 0.5)/n) : s.s1;
        beam_end a, b;
        a.body = st[e]; a.p = pa;
        b.body = st[e+1]; b.p = pb;
        auto el = std::make_unique<beam_element>(m_, a, b, sec, s.orient, s.beta, s.name + ".e" + std::to_string(e));
        const double L = el->length0();
        // lumped mass: half element per end, rotary inertia of the section and of the half element
        // (a rod of length L/2 about its end: m L^2/24 about the transverse axes)
        const double mh = 0.5*sec.m*L;
        Mat3 Jl = Mat3::Zero();
        Jl(0,0) = 0.5*sec.rIx*L;
        Jl(1,1) = 0.5*sec.rIy*L + mh*L*L/24.0;
        Jl(2,2) = 0.5*sec.rIz*L + mh*L*L/24.0;
        // element axes in the global frame (the node bodies start with R = I)
        Vec3 ey = s.orient - s.orient.dot(dir)*dir;
        if(ey.norm() < 1.0e-6)
        {
            const Vec3 alt = std::fabs(dir(2)) < 0.9 ? Vec3::UnitZ() : Vec3::UnitX();
            ey = alt - alt.dot(dir)*dir;
        }
        ey.normalize();
        Mat3 E;
        E.col(0) = dir; E.col(1) = ey; E.col(2) = dir.cross(ey);
        const Mat3 Jg = E*Jl*E.transpose();
        for(int end=0; end<2; ++end)
        {
            const int nb = st[e+end];
            const auto it = index_.find(nb);
            if(it != index_.end())
            {
                mass_[it->second] += mh;
                inertia_[it->second] += Jg;
            }
            else if(nb >= 0)
                m_.add_point_mass(nb, end ? pb : pa, mh, Jg);   // end on a rigid body: mass at the end point
        }
        m_.add_load(std::move(el));
    }
    return st;
}

void frame_builder::finish()
{
    if(finished_)
        return;
    for(size_t k=0; k<nodes_.size(); ++k)
    {
        rigid_body &b = m_.bodies[nodes_[k]];
        b.mass += mass_[k];
        // inertia in the node body frame (R0 = I for nodes created here)
        b.inertia += b.R.transpose()*inertia_[k]*b.R;
    }
    finished_ = true;
}

}
