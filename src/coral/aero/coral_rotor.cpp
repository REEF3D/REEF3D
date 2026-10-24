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
    rotor_handles h;
    const Vec3 e = sp.axis.normalized();
    const double sg = sp.ccw ? -1.0 : 1.0;
    Vec3 up = Vec3::UnitZ() - Vec3::UnitZ().dot(e)*e;
    if(up.norm() < 1.0e-6) up = Vec3::UnitX() - Vec3::UnitX().dot(e)*e;
    up.normalize();
    const double cb = std::cos(sp.precone), sb = std::sin(sp.precone);
    std::vector<Vec3> span(sp.B);
    for(int k=0; k<sp.B; ++k)
        span[k] = cb*rotate(up, e, sg*(sp.azimuth0 + 2.0*M_PI*k/sp.B)) + sb*e;
    const double R_tip = sp.stations.back().r;
    const double L = R_tip - sp.R_hub;

    // the rotor body: hub at H with the blade masses
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
    if(!sp.structure.frac.empty())
    {
        const blade_structure &bs = sp.structure;
        const int n = (int)bs.frac.size();
        for(int k=0; k<sp.B; ++k)
            for(int i=0; i<n; ++i)
            {
                const double w = 0.5*L*((i > 0 ? bs.frac[i] - bs.frac[i-1] : 0.0) + (i + 1 < n ? bs.frac[i+1] - bs.frac[i] : 0.0));
                const double mi = w*bs.mass[i];
                if(mi <= 0.0) continue;
                const Vec3 p = sp.hub + (sp.R_hub + bs.frac[i]*L)*span[k];
                m.add_point_mass(h.hub, p, mi);
                h.mass += mi;
                const Vec3 d = p - sp.hub;
                h.J += mi*(d - d.dot(e)*e).squaredNorm();
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

    // aerodynamic stations
    std::vector<aero_station> st;
    std::vector<Vec3> roots;
    const int n = (int)sp.stations.size();
    for(int k=0; k<sp.B; ++k)
    {
        const Vec3 z = span[k];
        const Vec3 nrm = (e - e.dot(z)*z).normalized();
        const Vec3 t = sg*e.cross(z).normalized();
        roots.push_back(sp.hub + sp.R_hub*z);
        for(int i=0; i<n; ++i)
        {
            const blade_station &bs = sp.stations[i];
            aero_station a;
            a.body = h.hub;
            const Vec3 p = sp.hub + bs.r*z;
            a.s = m.body(h.hub).to_body(p);
            const Mat3 &R = m.body(h.hub).R;
            a.n0 = R.transpose()*nrm; a.t0 = R.transpose()*t; a.z0 = R.transpose()*z;
            a.chord = bs.chord; a.twist = bs.twist; a.airfoil = bs.airfoil; a.blade = k;
            a.w = 0.5*((i > 0 ? bs.r - sp.stations[i-1].r : 0.0) + (i + 1 < n ? sp.stations[i+1].r - bs.r : 0.0));
            st.push_back(a);
        }
    }
    rotor_aero_settings as = sp.aero;
    as.B = sp.B;
    as.R_tip = R_tip*cb;
    as.R_hub = sp.R_hub*cb;
    h.R_tip = as.R_tip;
    for(const blade_station &bs : sp.stations)
        if(bs.airfoil < 0 || bs.airfoil >= (int)sp.airfoils.size())
            throw std::runtime_error("rotor " + sp.name + ": a blade station refers to airfoil " + std::to_string(bs.airfoil + 1) + " of "
                                     + std::to_string(sp.airfoils.size()));
    auto ra = std::make_unique<rotor_aero>(m, wind, sp.airfoils, h.hub, sp.hub, *h.shaft, st, roots, as, sp.name);
    for(int k=0; k<sp.B; ++k) ra->set_pitch(k, sp.pitch);
    if(!sp.pitch_t.empty()) ra->set_pitch_table(sp.pitch_t, sp.pitch_v);
    h.aero = &static_cast<rotor_aero&>(m.add_load(std::move(ra)));
    return h;
}

}
