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
    // the stations at the tip and at the hub (in the set-up geometry): no load with the tip and hub loss, whatever
    // rounding or deflection does to their radius later (a load switching on and off there stalls Newton)
    {
        const Vec3 e = shaft.axis(m);
        for(const aero_station &a : st_)
        {
            const Vec3 d = m.body(a.body).point(a.s) - H;
            const double r = (d - d.dot(e)*e).norm();
            edge_.push_back(r >= s.R_tip*(1.0 - 1.0e-9) ? 1 : (s.R_hub > 0.0 && r <= s.R_hub*(1.0 + 1.0e-9)) ? 2 : 0);
        }
    }
    if(s.dyn_stall.model > 0)
        for(const airfoil &a : af_)
            sp_.emplace_back(a, 0.0, s.dyn_stall);
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
                         station_load &o, const Vec3 &dV, int i, inflow_state *dio, stall_state *dso) const
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
    const int ie = int(&s - st_.data());
    if(ie >= 0 && ie < (int)edge_.size() && edge_[ie])
        el.r = edge_[ie]==1 ? set_.R_tip : set_.R_hub;
    el.af = (s.airfoil >= 0 && s.airfoil < (int)af_.size()) ? &af_[s.airfoil] : nullptr;
    o.bem = bem_solve(el, Vx, Vy, set_.rho, set_.bem);
    if(h.gam > 0.0 && o.bem.region != 0 && el.r > 0.0)
    {
        const double chi = (0.6*o.bem.a + 1.0)*h.gam;
        const double cpsi = dr.dot(h.cross)/el.r;
        const double K = 15.0*M_PI/32.0*(el.r/set_.R_tip)*std::tan(0.5*chi)*cpsi;
        o.bem = bem_solve(el, Vx, Vy, set_.rho, set_.bem, K);
    }
    if(i >= 0 && unsteady())
        unsteady_station(m, b, s, h, t, o, el, Vx, Vy, c, i, dio, dso);
    o.F = s.w*(o.bem.fn*o.n + o.bem.ft*o.t);
    o.M = s.w*o.bem.mz*o.t.cross(o.n);
}

void rotor_aero::unsteady_station(const model &m, const rigid_body &b, const aero_station &s, const hub_state &h, double t,
                                  station_load &o, const bem_element &el, double Vx, double Vy, const Vec3 &c, int i,
                                  inflow_state *dio, stall_state *dso) const
{
    // at the tip or the hub with their loss: no load, the states kept
    if(o.bem.region==0 && o.bem.F==0.0)
    {
        if(dio) *dio = i < (int)di_.size() ? di_[i] : inflow_state();
        if(dso) *dso = i < (int)ds_.size() ? ds_[i] : stall_state();
        return;
    }
    // the quasi-steady induced velocities (skewed wake included) and the dynamic ones
    const double Wqs[2] = {o.bem.a*Vx, o.bem.ap*Vy};
    const double dt = t - tc_;
    inflow_state di = inflow_steady(Wqs);
    double ux = Vx - Wqs[0], uy = Vy + Wqs[1];
    if(set_.dyn_inflow.on && started_)
    {
        di = inflow_step(di_[i], Wqs, tau1_, oye_tau2(tau1_, el.r/set_.R_tip), set_.dyn_inflow.k, dt);
        ux = Vx - di.W[0];
        uy = Vy + di.W[1];
    }
    if(dio) *dio = di;
    bem_result &r = o.bem;
    r.phi = std::atan2(ux, uy);
    r.alpha = r.phi - el.theta;
    r.W = std::hypot(ux, uy);
    airfoil_coef k = el.af ? el.af->coefficients(r.alpha, el.re) : airfoil_coef{};
    if(set_.dyn_stall.model > 0 && el.af)
    {
        // the section's nose-up angular velocity: its rotation about the span relative to the hub (pitch joints,
        // torsion) and the rate of a pitch table; the rotation of the hub itself leaves the chord's angle to the rotor
        // plane unchanged
        const Vec3 z = b.R*s.z0;
        const rigid_body &hb = m.body(hub_);
        const double wz = (b.R*b.w - hb.R*hb.w).dot(z);
        const Vec3 zc = z.cross(c);
        const double cn = c.dot(o.n), ct = c.dot(o.t);
        const double g = (-zc.dot(o.n)*ct + cn*zc.dot(o.t))/std::max(cn*cn + ct*ct, 1.0e-30);
        double pr = 0.0;
        if(!pt_.empty())
            pr = (pitch(s.blade, t + 1.0e-3) - pitch(s.blade, t - 1.0e-3))/2.0e-3;
        stall_input in;
        in.alpha = r.alpha; in.U = r.W; in.chord = el.chord;
        in.omega = -(g*wz + pr);
        const stall_polar &sp = sp_[s.airfoil];
        stall_state ds;
        if(started_)
            k = stall_step(sp, set_.dyn_stall.model, ds_[i], in, dt, ds);
        else
            ds = stall_steady(sp, set_.dyn_stall.model, in);
        if(dso) *dso = ds;
    }
    r.cl = k.cl; r.cd = k.cd; r.cm = k.cm;
    r.cn = k.cl*std::cos(r.phi) + k.cd*std::sin(r.phi);
    r.ct = k.cl*std::sin(r.phi) - k.cd*std::cos(r.phi);
    const double q = 0.5*set_.rho*r.W*r.W;
    r.fn = q*el.chord*r.cn;
    r.ft = q*el.chord*r.ct;
    r.mz = q*el.chord*el.chord*r.cm;
}

void rotor_aero::advance(const model &m, double t)
{
    const hub_state h = hub(m, t);
    const size_t n = st_.size();
    std::vector<inflow_state> di(n);
    std::vector<stall_state> ds(n);
    double sa = 0.0, sw = 0.0;
    for(size_t i=0; i<n; ++i)
    {
        station_load o;
        station(m, m.body(st_[i].body), st_[i], h, t, o, Vec3::Zero(), (int)i, &di[i], &ds[i]);
        // the mean axial induction (area weighted) of tau1
        const double r = (o.P - h.H - (o.P - h.H).dot(h.e)*h.e).norm();
        if(o.bem.region != 0)
        {
            sa += st_[i].w*r*o.bem.a;
            sw += st_[i].w*r;
        }
    }
    di_ = std::move(di);
    ds_ = std::move(ds);
    if(set_.dyn_inflow.tau1 > 0.0)
        tau1_ = set_.dyn_inflow.tau1;
    else
        tau1_ = oye_tau1(sw > 0.0 ? sa/sw : 0.0, set_.R_tip, wind_->velocity(h.H, t).norm());
    tc_ = t;
}

void rotor_aero::initialize(const model &m, double t0)
{
    if(!unsteady())
        return;
    started_ = false;
    advance(m, t0);
    started_ = true;
}

void rotor_aero::commit(const model &m, double t, double)
{
    if(!unsteady() || !started_)
        return;
    advance(m, t);
}

void rotor_aero::save_state(std::vector<double> &v) const
{
    if(!unsteady())
        return;
    v.push_back(started_ ? 1.0 : 0.0);
    v.push_back(tc_);
    v.push_back(tau1_);
    v.push_back((double)di_.size());
    for(size_t i=0; i<di_.size(); ++i)
    {
        const inflow_state &a = di_[i];
        for(int j=0; j<2; ++j) {v.push_back(a.Wqs[j]); v.push_back(a.Wint[j]); v.push_back(a.W[j]);}
        const stall_state &b = i < ds_.size() ? ds_[i] : stall_state();
        for(int j=0; j<4; ++j) v.push_back(b.x[j]);
        for(double x : {b.a34, b.clp, b.fF, b.fst}) v.push_back(x);
    }
}

bool rotor_aero::load_state(const std::vector<double> &v, size_t &k)
{
    if(!unsteady())
        return true;
    if(k + 4 > v.size())
        return false;
    started_ = v[k] != 0.0;
    tc_ = v[k+1];
    tau1_ = v[k+2];
    const size_t n = (size_t)v[k+3];
    k += 4;
    if(n != 0 && n != st_.size())
        return false;
    if(k + 14*n > v.size())
        return false;
    di_.assign(n, inflow_state());
    ds_.assign(n, stall_state());
    for(size_t i=0; i<n; ++i)
    {
        inflow_state &a = di_[i];
        for(int j=0; j<2; ++j) {a.Wqs[j] = v[k++]; a.Wint[j] = v[k++]; a.W[j] = v[k++];}
        stall_state &b = ds_[i];
        for(int j=0; j<4; ++j) b.x[j] = v[k++];
        b.a34 = v[k++]; b.clp = v[k++]; b.fF = v[k++]; b.fst = v[k++];
    }
    return true;
}

void rotor_aero::station_loads(const model &m, double t, std::vector<station_load> &out) const
{
    out.resize(st_.size());
    const hub_state h = hub(m, t);
    for(size_t i=0; i<st_.size(); ++i)
        station(m, m.body(st_[i].body), st_[i], h, t, out[i], Vec3::Zero(), (int)i);
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
    auto gen = [&](const rigid_body &b, const aero_station &s, int i)
    {
        station_load o;
        station(m, b, s, h, t, o, Vec3::Zero(), i);
        Vec6 g;
        g.head<3>() = o.F;
        g.tail<3>() = b.R.transpose()*((o.P - b.x).cross(o.F) + o.M);
        return g;
    };
    for(size_t is=0; is<st_.size(); ++is)
    {
        const aero_station &s = st_[is];
        const int i = (int)is;
        if(s.body < 0)
            continue;
        const rigid_body &b0 = m.body(s.body);
        Mat6 Kv = Mat6::Zero(), Kq = Mat6::Zero();
        const double ev = 1.0e-5*std::max(1.0, b0.point_velocity(s.s).norm()), er = 1.0e-7;
        // dF/dV and dM/dV of the station's force and moment by its relative inflow V = u - v_P
        station_load o0;
        station(m, b0, s, h, t, o0, Vec3::Zero(), i);
        Mat3 dF, dM;
        for(int j=0; j<3; ++j)
        {
            Vec3 d = Vec3::Zero(); d(j) = ev;
            station_load o1;
            station(m, b0, s, h, t, o1, d, i);
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
            const Vec6 G = (gen(rp, s, i) - gen(rm, s, i))/(2.0*er);
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
    if(set_.dyn_inflow.on)
        c.push_back(name_ + "_tau1");
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
    if(set_.dyn_inflow.on)
        v.push_back(tau1_);
}

}
