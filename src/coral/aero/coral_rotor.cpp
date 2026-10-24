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

#include"coral/aero/coral_rotor.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_model.h"
#include"coral/structure/coral_frame.h"
#include<cmath>
#include<stdexcept>

namespace coral
{

namespace
{
    Vec3 rotate(const Vec3 &v, const Vec3 &e, double a)
    {
        return v*std::cos(a) + e.cross(v)*std::sin(a) + e*e.dot(v)*(1.0 - std::cos(a));
    }
}

rotor_handles build_rotor(model &m, std::shared_ptr<const wind_inflow> wind, const rotor_spec &sp)
{
    if(sp.B < 1) throw std::runtime_error("rotor " + sp.name + ": at least one blade");
    if(sp.stations.size() < 2) throw std::runtime_error("rotor " + sp.name + ": at least two blade stations");
    if(!wind) throw std::runtime_error("rotor " + sp.name + ": no wind inflow");
    if(sp.blade_elements > 0 && sp.structure.frac.size() < 2)
        throw std::runtime_error("rotor " + sp.name + ": flexible blades need the structural table (structure)");
    for(const blade_station &bs : sp.stations)
        if(bs.airfoil < 0 || bs.airfoil >= (int)sp.airfoils.size())
            throw std::runtime_error("rotor " + sp.name + ": a blade station refers to airfoil " + std::to_string(bs.airfoil + 1) + " of "
                                     + std::to_string(sp.airfoils.size()));
    rotor_handles h;
    const Vec3 e = sp.axis.normalized();
    const double sg = sp.ccw ? -1.0 : 1.0;
    Vec3 up = Vec3::UnitZ() - Vec3::UnitZ().dot(e)*e;
    if(up.norm() < 1.0e-6) up = Vec3::UnitX() - Vec3::UnitX().dot(e)*e;
    up.normalize();
    const double cb = std::cos(sp.precone), sb = std::sin(sp.precone);
    std::vector<Vec3> span(sp.B), nrm(sp.B), tan(sp.B), pax(sp.B);
    for(int k=0; k<sp.B; ++k)
    {
        span[k] = cb*rotate(up, e, sg*(sp.azimuth0 + 2.0*M_PI*k/sp.B)) + sb*e;
        nrm[k] = (e - e.dot(span[k])*span[k]).normalized();
        tan[k] = sg*e.cross(span[k]).normalized();
        pax[k] = nrm[k].cross(tan[k]);          // pitch axis: positive towards feather
    }
    const double R_tip = sp.stations.back().r;
    const double L = R_tip - sp.R_hub;
    const bool flex = sp.blade_elements > 0;
    const bool pj = sp.pitch_joints;
    // the pitch in the geometry (pitch joints), else in the aerodynamics
    const double pitch0 = pj ? sp.pitch : 0.0;
    auto chord_at = [&](double r)
    {
        const auto &st = sp.stations;
        if(r <= st.front().r) return st.front().chord;
        if(r >= st.back().r) return st.back().chord;
        size_t q = 0;
        while(q + 2 < st.size() && st[q+1].r < r) ++q;
        return st[q].chord + (r - st[q].r)/(st[q+1].r - st[q].r)*(st[q+1].chord - st[q].chord);
    };
    auto chord_dir = [&](int k, double tw) {return Vec3(std::cos(tw + pitch0)*tan[k] - std::sin(tw + pitch0)*nrm[k]);};
    // a rigid body set turning with the rotor
    auto spin = [&](int b)
    {
        rigid_body &B = m.bodies[b];
        B.w = B.R.transpose()*(sp.omega0*e);
        B.v = sp.omega0*e.cross(B.x - sp.hub);
    };

    // the hub (with the rigid blades' masses unless they sit on pitch bearings)
    rigid_body hub;
    hub.name = sp.name + "_hub";
    hub.mass = sp.hub_mass > 0.0 ? sp.hub_mass : 1.0;
    const double J = sp.hub_J > 0.0 ? sp.hub_J : 1.0;
    const double Jt = sp.hub_Jt >= 0.0 ? sp.hub_Jt : 0.5*J;
    const Mat3 P = e*e.transpose();
    hub.inertia = J*P + Jt*(Mat3::Identity() - P);
    hub.x = sp.hub;
    hub.w = sp.omega0*e;
    h.hub = m.add_body(hub);
    h.mass = hub.mass;
    h.J = J;
    const blade_structure &bs = sp.structure;
    const int ns = (int)bs.frac.size();
    std::vector<double> Jp(sp.B, 0.0);        // blade inertia about the pitch axis (actuator)
    for(int k=0; k<sp.B; ++k)
        for(int i=0; i<ns; ++i)
        {
            const double w = 0.5*L*((i > 0 ? bs.frac[i] - bs.frac[i-1] : 0.0) + (i + 1 < ns ? bs.frac[i+1] - bs.frac[i] : 0.0));
            const double mi = w*bs.mass[i];
            const double r = sp.R_hub + bs.frac[i]*L;
            const Vec3 d = r*span[k];
            h.mass += mi;
            h.J += mi*(d - d.dot(e)*e).squaredNorm();
            const double c = chord_at(r);
            Jp[k] += mi*(0.04*0.04 + 0.2*0.2)*c*c;
        }

    // pitch bearings
    std::vector<int> root(sp.B, h.hub);
    if(pj)
    {
        for(int k=0; k<sp.B; ++k)
        {
            rigid_body pb;
            pb.name = sp.name + "_pitch" + std::to_string(k + 1);
            pb.x = sp.hub + sp.R_hub*span[k];
            const double mb = std::max(0.01*bs.total_mass(L), 1.0);
            pb.mass = mb;
            const double jp = std::max(flex ? 0.01*Jp[k] : Jp[k], 1.0);
            pb.inertia = jp*Mat3::Identity();
            root[k] = m.add_body(pb);
            spin(root[k]);
            h.pitch_bodies.push_back(root[k]);
            h.pitch_joints.push_back(&static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, h.hub, root[k], pb.x, pax[k],
                                                                                                  sp.name + "_pitch" + std::to_string(k + 1)))));
        }
    }

    // blade structure: point masses (rigid) or beams (flexible)
    std::vector<std::vector<int>> nodes(sp.B);
    std::vector<std::vector<double>> node_r(sp.B);
    if(flex)
    {
        frame_builder fb(m);
        const int n = sp.blade_elements;
        for(int k=0; k<sp.B; ++k)
        {
            nodes[k].push_back(root[k]);
            node_r[k].push_back(sp.R_hub);
            for(int q=0; q<n; ++q)
            {
                const double f1 = (double)q/n, f2 = (double)(q + 1)/n;
                auto section = [&](double f)
                {
                    auto at = [&](const std::vector<double> &v)
                    {
                        if(f <= bs.frac.front()) return v.front();
                        if(f >= bs.frac.back()) return v.back();
                        size_t i = 0;
                        while(i + 2 < bs.frac.size() && bs.frac[i+1] < f) ++i;
                        return v[i] + (f - bs.frac[i])/(bs.frac[i+1] - bs.frac[i])*(v[i+1] - v[i]);
                    };
                    beam_section s;
                    s.m = at(bs.mass);
                    s.EIy = at(bs.EI_flap);
                    s.EIz = at(bs.EI_edge);
                    s.EA = sp.ea_per_mass*s.m;
                    s.GJ = sp.gj_factor*s.EIy;
                    const double c = chord_at(sp.R_hub + f*L);
                    s.rIy = s.m*0.04*0.04*c*c;
                    s.rIz = s.m*0.2*0.2*c*c;
                    s.rIx = s.rIy + s.rIz;
                    return s;
                };
                member_spec ms;
                ms.name = sp.name + "_b" + std::to_string(k + 1) + "e" + std::to_string(q);
                ms.p1 = sp.hub + (sp.R_hub + f1*L)*span[k];
                ms.p2 = sp.hub + (sp.R_hub + f2*L)*span[k];
                ms.n = 1;
                ms.s1 = section(f1); ms.s2 = section(f2); ms.tapered = true;
                // the chord direction at the element's middle: local y (EIy flapwise, about the chord)
                double tw;
                {
                    const double f = 0.5*(f1 + f2);
                    size_t i = 0;
                    while(i + 2 < bs.frac.size() && bs.frac[i+1] < f) ++i;
                    tw = f <= bs.frac.front() ? bs.twist.front() : f >= bs.frac.back() ? bs.twist.back()
                         : bs.twist[i] + (f - bs.frac[i])/(bs.frac[i+1] - bs.frac[i])*(bs.twist[i+1] - bs.twist[i]);
                }
                ms.orient = chord_dir(k, tw);
                ms.beta = sp.blade_beta;
                ms.body1 = q==0 ? root[k] : -2;
                ms.body2 = -2;
                const std::vector<int> st = fb.add_member(ms);
                nodes[k].push_back(st.back());
                node_r[k].push_back(sp.R_hub + f2*L);
            }
        }
        fb.finish();
        for(int k=0; k<sp.B; ++k)
            for(size_t q=1; q<nodes[k].size(); ++q)
                spin(nodes[k][q]);
        h.blade_nodes = nodes;
    }
    else if(ns > 0)
    {
        for(int k=0; k<sp.B; ++k)
            for(int i=0; i<ns; ++i)
            {
                const double w = 0.5*L*((i > 0 ? bs.frac[i] - bs.frac[i-1] : 0.0) + (i + 1 < ns ? bs.frac[i+1] - bs.frac[i] : 0.0));
                const double mi = w*bs.mass[i];
                if(mi > 0.0)
                    m.add_point_mass(root[k], sp.hub + (sp.R_hub + bs.frac[i]*L)*span[k], mi);
            }
    }

    // shaft and drivetrain
    h.shaft = &static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, sp.nacelle, h.hub, sp.hub, e,
                                                                            sp.name + "_shaft")));
    const double Jg = sp.J_gen*sp.dt.N*sp.dt.N;
    if(sp.dt.k > 0.0)
    {
        rigid_body g;
        g.name = sp.name + "_gen";
        g.mass = 1.0;
        const double ja = Jg > 0.0 ? Jg : 1.0;
        g.inertia = ja*P + 1.0e-3*ja*(Mat3::Identity() - P);
        g.x = sp.hub;
        g.w = sp.omega0*e;
        h.gen = m.add_body(g);
        h.gen_joint = &static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, sp.nacelle, h.gen, sp.hub, e,
                                                                                    sp.name + "_gen_shaft")));
    }
    else if(Jg > 0.0)
        m.add_point_mass(h.hub, sp.hub, 0.0, Jg*P);
    h.drive = &static_cast<drivetrain&>(m.add_load(std::make_unique<drivetrain>(*h.shaft, h.gen_joint, sp.dt, sp.name + "_dt")));
    if(pj)
    {
        double jp = 0.0;
        for(double x : Jp) jp += x/sp.B;
        jp = std::max(jp, 1.0);
        const double w = 2.0*M_PI*sp.pitch_freq;
        std::vector<const joint*> js(h.pitch_joints.begin(), h.pitch_joints.end());
        auto pa = std::make_unique<pitch_actuator>(js, jp*w*w, 2.0*sp.pitch_zeta*jp*w, sp.pitch, sp.name + "_pitch");
        if(!sp.pitch_t.empty()) pa->set_table(sp.pitch_t, sp.pitch_v);
        h.pitch = &static_cast<pitch_actuator&>(m.add_load(std::move(pa)));
    }

    // aerodynamic stations: on the hub, the pitch bearing or the nearest node
    std::vector<aero_station> st;
    std::vector<Vec3> roots;
    const int n = (int)sp.stations.size();
    for(int k=0; k<sp.B; ++k)
    {
        const Vec3 z = span[k];
        roots.push_back(sp.hub + sp.R_hub*z);
        for(int i=0; i<n; ++i)
        {
            const blade_station &bst = sp.stations[i];
            aero_station a;
            a.body = root[k];
            if(flex)
            {
                size_t best = 0;
                for(size_t q=1; q<node_r[k].size(); ++q)
                    if(std::fabs(node_r[k][q] - bst.r) < std::fabs(node_r[k][best] - bst.r)) best = q;
                a.body = nodes[k][best];
            }
            const Vec3 p = sp.hub + bst.r*z;
            const rigid_body &B = m.body(a.body);
            a.s = B.to_body(p);
            a.c0 = B.R.transpose()*chord_dir(k, bst.twist);
            a.z0 = B.R.transpose()*z;
            a.chord = bst.chord; a.airfoil = bst.airfoil; a.blade = k;
            a.w = 0.5*((i > 0 ? bst.r - sp.stations[i-1].r : 0.0) + (i + 1 < n ? sp.stations[i+1].r - bst.r : 0.0));
            st.push_back(a);
        }
    }
    rotor_aero_settings as = sp.aero;
    as.B = sp.B;
    as.sense = sg;
    as.R_tip = R_tip*cb;
    as.R_hub = sp.R_hub*cb;
    h.R_tip = as.R_tip;
    auto ra = std::make_unique<rotor_aero>(m, wind, sp.airfoils, h.hub, sp.hub, *h.shaft, st, roots, as, sp.name);
    if(!pj)
    {
        for(int k=0; k<sp.B; ++k) ra->set_pitch(k, sp.pitch);
        if(!sp.pitch_t.empty()) ra->set_pitch_table(sp.pitch_t, sp.pitch_v);
    }
    h.aero = &static_cast<rotor_aero&>(m.add_load(std::move(ra)));
    return h;
}

}
