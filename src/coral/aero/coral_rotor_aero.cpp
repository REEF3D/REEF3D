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

#include"coral/aero/coral_rotor_aero.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include<algorithm>
#include<cmath>

namespace coral
{

rotor_aero::rotor_aero(const model &m, std::shared_ptr<const wind_inflow> wind, std::vector<airfoil> af, int hub, const Vec3 &H,
                       const joint &shaft, std::vector<aero_station> st, const std::vector<Vec3> &roots,
                       const rotor_aero_settings &s, const std::string &name)
    : name_(name), wind_(std::move(wind)), af_(std::move(af)), hub_(hub), shaft_(&shaft), st_(std::move(st)), set_(s)
{
    sH_ = m.body(hub).to_body(H);
    for(const Vec3 &r : roots)
        roots_.push_back(m.body(hub).to_body(r));
    tb_ = m.body(s.tower_body).to_body(s.tower_base);
    tt_ = m.body(s.tower_body).to_body(s.tower_top);
    pitch_.assign(std::max(s.B, 1), 0.0);
}

void rotor_aero::coupled_bodies(std::vector<int> &b) const
{
    b.push_back(hub_);
    if(shaft_->body_a >= 0) b.push_back(shaft_->body_a);
    for(const aero_station &s : st_)
        if(std::find(b.begin(), b.end(), s.body)==b.end()) b.push_back(s.body);
    if(set_.tower && set_.tower_body >= 0 && std::find(b.begin(), b.end(), set_.tower_body)==b.end()) b.push_back(set_.tower_body);
}

double rotor_aero::pitch(int blade, double t) const
{
    if(pt_.empty())
        return pitch_[blade];
    if(t <= pt_.front()) return pv_.front();
    if(t >= pt_.back()) return pv_.back();
    const size_t k = std::upper_bound(pt_.begin(), pt_.end(), t) - pt_.begin() - 1;
    return pv_[k] + (t - pt_[k])/(pt_[k+1] - pt_[k])*(pv_[k+1] - pv_[k]);
}

Vec3 rotor_aero::hub_centre(const model &m) const {return m.body(hub_).point(sH_);}

Vec3 rotor_aero::tower_influence(const model &m, const Vec3 &P, const Vec3 &u) const
{
    if(!set_.tower)
        return Vec3::Zero();
    const rigid_body &tb = m.body(set_.tower_body);
    const Vec3 b = tb.point(tb_), top = tb.point(tt_);
    const Vec3 L = top - b;
    const double len = L.norm();
    if(len <= 0.0) return Vec3::Zero();
    const Vec3 a = L/len;
    const double z = (P - b).dot(a)/len;
    if(z < 0.0 || z > 1.0) return Vec3::Zero();
    const Vec3 un = u - u.dot(a)*a;
    const double U = un.norm();
    if(U <= 0.0) return Vec3::Zero();
    const Vec3 ex = un/U, ey = a.cross(ex);
    const double rt = 0.5*(set_.D_base + z*(set_.D_top - set_.D_base));
    if(rt <= 0.0) return Vec3::Zero();
    const Vec3 d = P - (b + z*len*a);
    const double x = d.dot(ex)/rt, y = d.dot(ey)/rt;
    const double q = x*x + y*y;
    if(q <= 1.0) return Vec3::Zero();
    return U*(-(x*x - y*y)/(q*q))*ex + U*(-2.0*x*y/(q*q))*ey;
}

rotor_aero::hub_state rotor_aero::hub(const model &m, double t) const
{
    hub_state h;
    const rigid_body &hb = m.body(hub_);
    h.H = hb.point(sH_);
    h.e = shaft_->axis(m);
    h.gam = 0.0;
    h.cross = Vec3::Zero();
    if(set_.skew)
    {
        const Vec3 Vh = wind_->velocity(h.H, t) - hb.point_velocity(sH_);
        const double vn = Vh.dot(h.e);
        const Vec3 c = Vh - vn*h.e;
        if(c.norm() > 1.0e-9*std::max(Vh.norm(), 1.0e-12))
        {
            h.gam = std::atan2(c.norm(), vn);
            h.cross = c.normalized();
        }
    }
    return h;
}

void rotor_aero::station(const model &m, const rigid_body &b, const aero_station &s, const hub_state &h, double t,
                         station_load &o, const Vec3 &dV) const
{
    const Vec3 &e = h.e;
    o.P = b.point(s.s);
    const Vec3 z = b.R*s.z0, c = b.R*s.c0;
    o.n = (e - e.dot(z)*z).normalized();
    o.t = set_.sense*e.cross(z).normalized();
    o.theta = std::atan2(-c.dot(o.n), c.dot(o.t)) + pitch(s.blade, t);
    Vec3 u = wind_->velocity(o.P, t);
    u += tower_influence(m, o.P, u);
    const Vec3 V = u - b.point_velocity(s.s) + dV;
    const double Vx = V.dot(o.n), Vy = -V.dot(o.t);
    const Vec3 d = o.P - h.H;
    const Vec3 dr = d - d.dot(e)*e;
    bem_element el;
    el.r = dr.norm(); el.chord = s.chord; el.theta = o.theta;
    el.R = set_.R_tip; el.R_hub = set_.R_hub; el.B = set_.B;
    el.af = (s.airfoil >= 0 && s.airfoil < (int)af_.size()) ? &af_[s.airfoil] : nullptr;
    o.bem = bem_solve(el, Vx, Vy, set_.rho, set_.bem);
    if(h.gam > 0.0 && o.bem.region != 0 && el.r > 0.0)
    {
        const double chi = (0.6*o.bem.a + 1.0)*h.gam;
        const double cpsi = dr.dot(h.cross)/el.r;
        const double K = 15.0*M_PI/32.0*(el.r/set_.R_tip)*std::tan(0.5*chi)*cpsi;
        o.bem = bem_solve(el, Vx, Vy, set_.rho, set_.bem, K);
    }
    o.F = s.w*(o.bem.fn*o.n + o.bem.ft*o.t);
    o.M = s.w*o.bem.mz*o.t.cross(o.n);
}

void rotor_aero::station_loads(const model &m, double t, std::vector<station_load> &out) const
{
    out.resize(st_.size());
    const hub_state h = hub(m, t);
    for(size_t i=0; i<st_.size(); ++i)
        station(m, m.body(st_[i].body), st_[i], h, t, out[i]);
}

void rotor_aero::eval(const model &m, double t, VecX &f) const
{
    std::vector<station_load> L;
    station_loads(m, t, L);
    for(size_t i=0; i<st_.size(); ++i)
    {
        m.add_force(f, st_[i].body, L[i].P, L[i].F);
        m.add_moment(f, st_[i].body, L[i].M);
    }
}

void rotor_aero::tangent(const model &m, double t, double cK, double cC, double, std::vector<triplet> &T) const
{
    if(cC==0.0 && cK==0.0)
        return;
    const hub_state h = hub(m, t);
    // the generalized force of one station on its body: [F (global); M about the CoG (body frame)]
    auto gen = [&](const rigid_body &b, const aero_station &s)
    {
        station_load o;
        station(m, b, s, h, t, o);
        Vec6 g;
        g.head<3>() = o.F;
        g.tail<3>() = b.R.transpose()*((o.P - b.x).cross(o.F) + o.M);
        return g;
    };
    for(const aero_station &s : st_)
    {
        if(s.body < 0)
            continue;
        const rigid_body &b0 = m.body(s.body);
        Mat6 Kv = Mat6::Zero(), Kq = Mat6::Zero();
        const double ev = 1.0e-5*std::max(1.0, b0.point_velocity(s.s).norm()), er = 1.0e-7;
        // dF/dV and dM/dV of the station's force and moment by its relative inflow V = u - v_P
        station_load o0;
        station(m, b0, s, h, t, o0);
        Mat3 dF, dM;
        for(int j=0; j<3; ++j)
        {
            Vec3 d = Vec3::Zero(); d(j) = ev;
            station_load o1;
            station(m, b0, s, h, t, o1, d);
            dF.col(j) = (o1.F - o0.F)/ev;
            dM.col(j) = (o1.M - o0.M)/ev;
        }
        // v_P = v + R (w x s): dV/dv = -I, dV/dw = R [s]x
        const Vec3 r = o0.P - b0.x;
        auto skew = [](const Vec3 &a) {Mat3 S; S << 0, -a(2), a(1), a(2), 0, -a(0), -a(1), a(0), 0; return S;};
        const Mat3 dVdw = b0.R*skew(s.s);
        const Mat3 Rt = b0.R.transpose();
        const Mat3 dGm = Rt*(skew(r)*dF + dM);               // moment (body frame) by V
        Kv.block<3,3>(0, 0) = -dF;
        Kv.block<3,3>(3, 0) = -dGm;
        Kv.block<3,3>(0, 3) = dF*dVdw;
        Kv.block<3,3>(3, 3) = dGm*dVdw;
        // rotation about the span: pitch and torsion change the angle of attack (the other rotations, which tilt the
        // blade against the shaft axis of the nacelle, are left out)
        {
            rigid_body rp = b0, rm = b0;
            rp.R = b0.R*exp_so3(er*s.z0); rm.R = b0.R*exp_so3(-er*s.z0);
            const Vec6 G = (gen(rp, s) - gen(rm, s))/(2.0*er);
            Kq.rightCols<3>() = G*s.z0.transpose();
        }
        const int o = 6*s.body;
        for(int i=0; i<6; ++i)
            for(int j=0; j<6; ++j)
            {
                const double v = cK*Kq(i, j) + cC*Kv(i, j);
                if(v != 0.0)
                    T.emplace_back(o + i, o + j, v);
            }
    }
}

rotor_aero::rotor_values rotor_aero::values(const model &m, double t, const std::vector<station_load> &L) const
{
    rotor_values v{};
    const Vec3 H = hub_centre(m), e = shaft_->axis(m);
    for(const station_load &s : L)
    {
        v.T += s.F.dot(e);
        v.Q += ((s.P - H).cross(s.F) + s.M).dot(e);
    }
    v.omega = shaft_->rotation_rate(m);
    v.P = v.Q*v.omega;
    const Vec3 uh = wind_->velocity(H, t);
    v.U = uh.norm();
    const double A = M_PI*set_.R_tip*set_.R_tip;
    if(v.U > 0.0)
    {
        v.Cp = v.P/(0.5*set_.rho*A*v.U*v.U*v.U);
        v.Ct = v.T/(0.5*set_.rho*A*v.U*v.U);
        v.lambda = v.omega*set_.R_tip/v.U;
        v.yaw = std::acos(std::clamp(uh.dot(e)/v.U, -1.0, 1.0));
    }
    v.azimuth = shaft_->rotation_coordinate(m);
    return v;
}

void rotor_aero::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"azimuth", "rpm", "T", "Q", "P", "Cp", "Ct", "lambda", "U", "yaw", "pitch"})
        c.push_back(name_ + "_" + x);
    for(int k=0; k<set_.B; ++k)
    {
        c.push_back(name_ + "_b" + std::to_string(k + 1) + "_Moop");
        c.push_back(name_ + "_b" + std::to_string(k + 1) + "_Mip");
    }
}

void rotor_aero::channel_values(const model &m, double t, std::vector<double> &v) const
{
    std::vector<station_load> L;
    station_loads(m, t, L);
    const rotor_values r = values(m, t, L);
    const double deg = 180.0/M_PI;
    double az = r.azimuth*deg;
    if(az < 0.0) az += 360.0;
    for(double x : {az, r.omega*30.0/M_PI, r.T, r.Q, r.P, r.Cp, r.Ct, r.lambda, r.U, r.yaw*deg, pitch(0, t)*deg})
        v.push_back(x);
    const rigid_body &hb = m.body(hub_);
    for(int k=0; k<set_.B; ++k)
    {
        const Vec3 root = k < (int)roots_.size() ? hb.point(roots_[k]) : hub_centre(m);
        Vec3 M = Vec3::Zero(), n = Vec3::Zero(), td = Vec3::Zero(), z = Vec3::Zero();
        bool first = true;
        for(size_t i=0; i<st_.size(); ++i)
        {
            if(st_[i].blade != k) continue;
            M += (L[i].P - root).cross(L[i].F) + L[i].M;
            if(first)
            {
                n = L[i].n; td = L[i].t; z = m.body(st_[i].body).R*st_[i].z0;
                first = false;
            }
        }
        v.push_back(M.dot(z.cross(n)));     // out of the plane: positive for a downwind thrust
        v.push_back(M.dot(z.cross(td)));    // in the plane: positive for a driving force
    }
}

}
