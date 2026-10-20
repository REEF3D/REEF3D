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

#include"coral/mooring/coral_mooring.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>
#include<stdexcept>

namespace coral
{

static double gmag(const model &m) {return m.gravity.norm() > 0.0 ? m.gravity.norm() : 9.81;}

double default_seabed(const mooring_end &a, const mooring_end &b)
{
    const bool ga = a.body < 0 && a.point < 0, gb = b.body < 0 && b.point < 0;
    if(ga && gb) return std::min(a.p(2), b.p(2));
    if(ga) return a.p(2);
    if(gb) return b.p(2);
    return -1.0e30;
}

// ---------------------------------------------------------------- quasi-static

mooring_qs::mooring_qs(const model &m, const mooring_end &a, const mooring_end &b, const mooring_props &p, const std::string &name,
                       double z_seabed)
    : name_(name), ba_(a.body), bb_(b.body), pr_(p)
{
    // the seabed portion: at a ground end on (or below) the seabed, the lower ground end of two
    const double zbed = std::isnan(z_seabed) ? default_seabed(a, b) : z_seabed;
    const double tol = 1.0e-6*std::max(p.L, 1.0);
    const bool a_bed = a.body < 0 && a.p(2) <= zbed + tol, b_bed = b.body < 0 && b.p(2) <= zbed + tol;
    if(p.CB >= 0.0 && (a_bed || b_bed))
        bed_end_ = (a_bed && (!b_bed || a.p(2) <= b.p(2))) ? 0 : 1;
    sa_ = m.body(a.body).to_body(a.p);
    sb_ = m.body(b.body).to_body(b.p);
    cl_.L = p.L;
    cl_.EA = p.EA;
    cl_.CB = p.CB;
    cl_.w = (p.mass - p.rho*0.25*M_PI*p.d*p.d)*gmag(m);
    if(cl_.w <= 0.0)
        throw std::runtime_error("coral::mooring_qs: the line must sink (mass > rho A)");
}

void mooring_qs::forces(const model &m, Vec3 &Fa, Vec3 &Fb) const
{
    const Vec3 pa = m.body(ba_).point(sa_), pb = m.body(bb_).point(sb_);
    const bool a_low = pa(2) <= pb(2);
    const Vec3 lo = a_low ? pa : pb, up = a_low ? pb : pa;
    Vec3 dh = up - lo;
    const double zF = dh(2);
    dh(2) = 0.0;
    const double xF = dh.norm();
    const Vec3 eh = xF > 0.0 ? Vec3(dh/xF) : Vec3::UnitX();
    // seabed portion only from the ground end on the seabed while it is the lower end
    catenary_line cl = cl_;
    if(bed_end_ < 0 || (bed_end_==0) != a_low)
        cl.CB = -1.0;
    double H = H_, V = V_;
    if(!catenary_solve(cl, xF, zF, H, V))
    {
        H = -1.0;
        if(!catenary_solve(cl, xF, zF, H, V))
            throw std::runtime_error("coral::mooring_qs: no catenary solution for line " + name_);
    }
    H_ = H;
    V_ = V;
    double HA, VA;
    catenary_anchor(cl, H, V, HA, VA);
    const Vec3 Fup = -H*eh - V*Vec3::UnitZ();
    const Vec3 Flo = HA*eh + VA*Vec3::UnitZ();
    Fa = a_low ? Flo : Fup;
    Fb = a_low ? Fup : Flo;
}

void mooring_qs::eval(const model &m, double, VecX &f) const
{
    Vec3 Fa, Fb;
    forces(m, Fa, Fb);
    m.add_force(f, ba_, m.body(ba_).point(sa_), Fa);
    m.add_force(f, bb_, m.body(bb_).point(sb_), Fb);
}

void mooring_qs::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"Fa_x", "Fa_y", "Fa_z", "Fb_x", "Fb_y", "Fb_z"})
        c.push_back(name_ + "_" + x);
}

void mooring_qs::channel_values(const model &m, double, std::vector<double> &v) const
{
    Vec3 Fa, Fb;
    forces(m, Fa, Fb);
    for(int i=0; i<3; ++i) v.push_back(Fa(i));
    for(int i=0; i<3; ++i) v.push_back(Fb(i));
}

// ---------------------------------------------------------------- mooring system: dynamics
// (the set-up, outputs and checkpoints in coral_mooring_setup.cpp, the static equilibrium in
// coral_mooring_static.cpp: next to the sub-step loop here they change the compiler's inlining of it)

const mooring_system::node_kin &mooring_system::kinematics(int slot, double t) const
{
    node_kin &k = kin_[slot];
    if(k.gen==gen_ && k.t==t)
        return k;
    const int n = (int)r_.size();
    k.t = t;
    k.gen = gen_;
    k.eta.resize(n);
    k.u.resize(n);
    k.du.resize(n);
    for(int i=0; i<n; ++i)
    {
        k.eta[i] = waves_->elevation_at(this, i, r_[i], t);
        double pr;
        waves_->kinematics_at(this, i, r_[i], t, k.u[i], k.du[i], pr);
    }
    return k;
}

void mooring_system::water(int i, const Vec3 &x, double t, double &eta, Vec3 &u, Vec3 &du) const
{
    switch(kmode_)
    {
    case 1:
    {
        eta = waves_->elevation_at(this, i, x, t);
        double pr;
        waves_->kinematics_at(this, i, x, t, u, du, pr);
        return;
    }
    case 2:
    {
        const double f = kt1_ > kt0_ ? (t - kt0_)/(kt1_ - kt0_) : 1.0;
        const node_kin &a = kin_[0], &b = kin_[1];
        eta = a.eta[i] + f*(b.eta[i] - a.eta[i]);
        u = a.u[i] + f*(b.u[i] - a.u[i]);
        du = a.du[i] + f*(b.du[i] - a.du[i]);
        return;
    }
    case 3:
        eta = kin_[0].eta[i];
        u.setZero();
        du.setZero();
        return;
    default:
        eta = 0.0;
        u.setZero();
        du.setZero();
    }
}

void mooring_system::line_accel(const line_t &L, const state &s, double t, std::vector<Vec3> &a, double extra) const
{
    const int N = L.N, o = L.off;
    const mooring_props &p = L.p;
    const double A = 0.25*M_PI*p.d*p.d;
    const double l_ = L.l, zbed_ = L.zbed;
    const Vec3 &g = g_;
    auto R = [&](int i) -> const Vec3& {return s.r[o+i];};
    auto V = [&](int i) -> const Vec3& {return s.v[o+i];};
    std::vector<Vec3> &T_ = L.T;
    // segment tensions (no compression)
    T_.resize(N);
    for(int j=0; j<N; ++j)
    {
        const Vec3 d = R(j+1) - R(j);
        const double Ls = d.norm();
        double Tj = 0.0;
        if(Ls > l_)
        {
            const Vec3 e = d/Ls;
            const double sr = (V(j+1) - V(j)).dot(e)/l_;
            Tj = std::max(0.0, p.EA*(Ls/l_ - 1.0) + p.BA*sr)/Ls;
        }
        else if(creg_ > 0.0 && Ls > 0.0)
            Tj = creg_*p.EA*(Ls/l_ - 1.0)/Ls;      // static solve only: a compressed segment pushes a little
        T_[j] = Tj*d;
    }
    // external forces on a node of length ll (weight, hydrodynamics, seabed); q tangent; added mass
    // alpha (normal) and beta (tangential)
    auto external = [&](int node, const Vec3 &x, const Vec3 &vx, const Vec3 &q, double ll, double &an, double &at)
    {
        Vec3 F = p.mass*ll*g;
        an = at = 0.0;
        double eta;
        Vec3 u, aw;
        water(o + node, x, t, eta, u, aw);
        // submerged fraction of the node's piece of line: linear over its vertical extent (at least d), so
        // a node crossing the free surface takes its buoyancy, drag and added mass gradually
        const double hz = std::max(p.d, std::fabs(q(2))*ll);
        const double fs = std::clamp((eta - x(2))/hz + 0.5, 0.0, 1.0);
        if(fs > 0.0)
        {
            F -= fs*p.rho*A*ll*g;
            const Vec3 uw = current_ + u;
            const Vec3 vr = uw - vx;
            const double vq = vr.dot(q);
            const Vec3 vrt = vq*q, vrn = vr - vrt;
            F += fs*(0.5*p.rho*p.d*ll*p.Cdn*vrn.norm()*vrn + 0.5*p.rho*M_PI*p.d*ll*p.Cdt*std::fabs(vq)*vrt);
            if(kmode_==1 || kmode_==2)
            {
                const Vec3 awt = aw.dot(q)*q, awn = aw - awt;
                F += fs*p.rho*A*ll*((1.0 + p.Can)*awn + (1.0 + p.Cat)*awt);
            }
            an = fs*p.rho*A*ll*p.Can;
            at = fs*p.rho*A*ll*p.Cat;
        }
        // seabed: penalty and damping normal to it, Coulomb friction (regularized) with mu across and mu_axial
        // along the line (elliptic limit: F = -Fz (mua^2 va ea + mu^2 vt et) / |(mua va, mu vt)|, the
        // isotropic Coulomb law for mua = mu)
        if(x(2) < zbed_)
        {
            const double Fz = std::max(0.0, (p.kb*(zbed_ - x(2)) - p.cb*vx(2))*p.d*ll);
            F(2) += Fz;
            const double mua = p.mu_axial >= 0.0 ? p.mu_axial : p.mu;
            const Vec3 vh(vx(0), vx(1), 0.0);
            const double vn = vh.norm();
            if(Fz > 0.0 && vn > 0.0 && (p.mu > 0.0 || mua > 0.0))
            {
                Vec3 e(q(0), q(1), 0.0);
                const double en = e.norm();
                e = en > 1.0e-6 ? Vec3(e/en) : Vec3(vh/vn);
                const Vec3 va = vh.dot(e)*e, vt = vh - va;
                const double wn = std::sqrt(mua*mua*va.squaredNorm() + p.mu*p.mu*vt.squaredNorm());
                if(wn > 0.0)
                    F -= (Fz*std::tanh(vn/p.v_eps)/wn)*(mua*mua*va + p.mu*p.mu*vt);
            }
        }
        return F;
    };
    const double mi = p.mass*l_;
    double an, at;
    for(int i=1; i<N; ++i)
    {
        Vec3 q = R(i+1) - R(i-1);
        const double qn = q.norm();
        q = qn > 0.0 ? Vec3(q/qn) : Vec3::UnitZ();
        const Vec3 F = T_[i] - T_[i-1] + external(i, R(i), V(i), q, l_, an, at) - extra*mi*V(i);
        // (m I + an (I - q q^T) + at q q^T)^-1 F
        const double al = mi + an, be = at - an;
        const double fq = F.dot(q);
        a[o+i] = (F - (be/(al + be))*fq*q)/al;
    }
    // end forces: tension of the end segment, weight and hydrodynamics of the end half-segment; its mass
    // and added mass (carried by the attached body or point)
    Vec3 qa = R(1) - R(0), qb = R(N) - R(N-1);
    qa = qa.norm() > 0.0 ? Vec3(qa/qa.norm()) : Vec3::UnitZ();
    qb = qb.norm() > 0.0 ? Vec3(qb/qb.norm()) : Vec3::UnitZ();
    const double me = 0.5*mi;
    L.F[0] = T_[0] + external(0, R(0), V(0), qa, 0.5*l_, an, at);
    L.M[0] = me*Mat3::Identity() + an*(Mat3::Identity() - qa*qa.transpose()) + at*qa*qa.transpose();
    L.F[1] = -T_[N-1] + external(N, R(N), V(N), qb, 0.5*l_, an, at);
    L.M[1] = me*Mat3::Identity() + an*(Mat3::Identity() - qb*qb.transpose()) + at*qb*qb.transpose();
    a[o].setZero();
    a[o+N].setZero();
}

void mooring_system::accel(const state &s, double t, std::vector<Vec3> &a, double extra) const
{
    a.resize(r_.size());
    for(const line_t &L : lines_)
        line_accel(L, s, t, a, extra);
    // free points: own weight, buoyancy, drag, added mass, seabed; the attached line ends
    for(const point_t &P : points_)
    {
        const mooring_point_props &p = P.p;
        const int i = P.idx;
        const Vec3 &x = s.r[i], &vx = s.v[i];
        double eta;
        Vec3 u, aw;
        water(i, x, t, eta, u, aw);
        const double rho = lines_.empty() ? 1025.0 : lines_[P.ends.empty() ? 0 : P.ends[0].first].p.rho;
        const double rv = std::cbrt(0.75*p.volume/M_PI);                  // volume-equivalent sphere
        const double hz = std::max(2.0*rv, 0.1);
        const double fs = std::clamp((eta - x(2))/hz + 0.5, 0.0, 1.0);
        Vec3 F = p.mass*g_ - fs*rho*p.volume*g_;
        Mat3 M = (p.mass + fs*rho*p.volume*p.Ca)*Mat3::Identity();
        if(fs > 0.0)
        {
            const Vec3 vr = current_ + u - vx;
            F += fs*0.5*rho*p.CdA*vr.norm()*vr;
            if(kmode_==1 || kmode_==2)
                F += fs*rho*p.volume*(1.0 + p.Ca)*aw;
        }
        // seabed: contact of the sphere of the contact radius
        const double zb = P.zbed;
        const double rc = p.radius >= 0.0 ? p.radius : std::max(rv, 0.05);
        if(x(2) - rc < zb)
        {
            const double Ac = M_PI*rc*rc;
            const double Fz = std::max(0.0, (p.kb*(zb - x(2) + rc) - p.cb*vx(2))*Ac);
            F(2) += Fz;
            const Vec3 vh(vx(0), vx(1), 0.0);
            const double vn = vh.norm();
            if(Fz > 0.0 && vn > 0.0 && p.mu > 0.0)
                F -= (p.mu*Fz*std::tanh(vn/p.v_eps)/vn)*vh;
        }
        for(const auto &e : P.ends)
        {
            const line_t &L = lines_[e.first];
            F += L.F[e.second];
            M += L.M[e.second];
        }
        F -= extra*M*vx;
        a[i] = M.inverse()*F;
    }
}

void mooring_system::point_ends(state &s) const
{
    for(const line_t &L : lines_)
        for(int e=0; e<2; ++e)
            if(L.pt[e] >= 0)
            {
                const int i = L.off + (e==0 ? 0 : L.N), j = points_[L.pt[e]].idx;
                s.r[i] = s.r[j];
                s.v[i] = s.v[j];
            }
}

void mooring_system::integrate(state &s, const std::vector<end_motion> &em, double t0, double t1, double extra, bool at_rest) const
{
    const double T = t1 - t0;
    // the water: at rest (settling), exact, or from the converged nodes at t0 and t1
    if(!waves_)
        kmode_ = 0;
    else if(at_rest)
    {
        kinematics(0, t0);
        kmode_ = 3;
    }
    else if(exact_kin_)
        kmode_ = 1;
    else
    {
        kinematics(0, t0);
        kinematics(1, t1);
        kt0_ = t0;
        kt1_ = t1;
        kmode_ = 2;
    }
    const int n = std::max(1, (int)std::ceil(T/dt_ - 1.0e-9));
    const double h = T/n;
    std::vector<Vec3> &k1 = w_[0], &k2 = w_[1], &k3 = w_[2], &k4 = w_[3];
    state &s2 = ws_;
    s2 = s;
    const int nn = (int)s.r.size();
    std::vector<Vec3> rs(nn), vs(nn);
    // the body ends on cubic Hermite curves through the converged and the trial positions and velocities;
    // the line ends on points take the points' state
    auto ends = [&](state &x, double tt)
    {
        for(size_t k=0; k<att_.size(); ++k)
        {
            const end_motion &e = em[k];
            Vec3 &r = x.r[att_[k].idx], &v = x.v[att_[k].idx];
            if(T <= 0.0)
            {
                r = e.r1;
                v = e.v1;
                continue;
            }
            const double y = (tt - t0)/T, y2 = y*y, y3 = y2*y;
            r = (2.0*y3 - 3.0*y2 + 1.0)*e.r0 + T*(y3 - 2.0*y2 + y)*e.v0 + (3.0*y2 - 2.0*y3)*e.r1 + T*(y3 - y2)*e.v1;
            v = (6.0*(y2 - y)/T)*(e.r0 - e.r1) + (3.0*y2 - 4.0*y + 1.0)*e.v0 + (3.0*y2 - 2.0*y)*e.v1;
        }
        point_ends(x);
    };
    for(int k=0; k<n; ++k)
    {
        const double tk = t0 + k*h;
        // RK4 on (r, v) of the inner nodes and the points; the stage velocities are the stage states' v
        ends(s, tk);
        accel(s, tk, k1, extra);
        for(int i : free_) {s2.r[i] = s.r[i] + 0.5*h*s.v[i]; s2.v[i] = s.v[i] + 0.5*h*k1[i];}
        rs = s.v;
        ends(s2, tk + 0.5*h);
        accel(s2, tk + 0.5*h, k2, extra);
        vs = s2.v;
        for(int i : free_) {s2.r[i] = s.r[i] + 0.5*h*vs[i]; s2.v[i] = s.v[i] + 0.5*h*k2[i];}
        for(int i : free_) rs[i] += 2.0*vs[i];
        point_ends(s2);
        accel(s2, tk + 0.5*h, k3, extra);
        vs = s2.v;
        for(int i : free_) {s2.r[i] = s.r[i] + h*vs[i]; s2.v[i] = s.v[i] + h*k3[i];}
        for(int i : free_) rs[i] += 2.0*vs[i];
        ends(s2, tk + h);
        accel(s2, tk + h, k4, extra);
        for(int i : free_)
        {
            s.r[i] += h/6.0*(rs[i] + s2.v[i]);
            s.v[i] += h/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        }
    }
    for(size_t k=0; k<att_.size(); ++k)
    {
        s.r[att_[k].idx] = em[k].r1;
        s.v[att_[k].idx] = em[k].v1;
    }
    point_ends(s);
    accel(s, t1, k1, extra);
}

void mooring_system::eval(const model &m, double t, VecX &f) const
{
    const size_t na = att_.size();
    bool hit = cvalid_ && t==ct_;
    std::vector<Vec3> r1(na), v1(na);
    for(size_t k=0; k<na; ++k)
    {
        const rigid_body &B = m.body(att_[k].body);
        r1[k] = B.point(att_[k].s);
        v1[k] = B.point_velocity(att_[k].s);
        hit = hit && r1[k]==cr_[k] && v1[k]==cv_[k];
    }
    if(!hit)
    {
        std::vector<end_motion> em(na);
        for(size_t k=0; k<na; ++k)
            em[k] = end_motion{att_[k].r, att_[k].v, r1[k], v1[k]};
        state s{r_, v_};
        integrate(s, em, t_, t, 0.0);
        cr_ = r1;
        cv_ = v1;
        cF_.resize(na);
        cM_.resize(na);
        for(size_t k=0; k<na; ++k)
        {
            const line_t &L = lines_[att_[k].line];
            cF_[k] = L.F[att_[k].end];
            cM_[k] = L.M[att_[k].end];
        }
        ct_ = t;
        cvalid_ = true;
    }
    // the body carries the end node: minus its mass and added mass times the attachment's acceleration
    for(size_t k=0; k<na; ++k)
    {
        const attach_t &a = att_[k];
        m.add_force(f, a.body, r1[k], cF_[k] - cM_[k]*m.body(a.body).point_acceleration(a.s));
    }
}

void mooring_system::kinematics_points(const model&, const void *source, std::vector<Vec3> &pts) const
{
    // the converged nodes of the lines (0 .. N each), then the points
    if(waves_.get() != source)
        return;
    pts.insert(pts.end(), r_.begin(), r_.end());
}

void mooring_system::commit(const model &m, double t, double)
{
    const size_t na = att_.size();
    std::vector<end_motion> em(na);
    for(size_t k=0; k<na; ++k)
    {
        const rigid_body &B = m.body(att_[k].body);
        em[k] = end_motion{att_[k].r, att_[k].v, B.point(att_[k].s), B.point_velocity(att_[k].s)};
    }
    state s{r_, v_};
    integrate(s, em, t_, t, 0.0);
    for(size_t k=0; k<na; ++k)
    {
        attach_t &a = att_[k];
        const line_t &L = lines_[a.line];
        a.F = L.F[a.end] - L.M[a.end]*m.body(a.body).point_acceleration(a.s);
        a.r = em[k].r1;
        a.v = em[k].v1;
    }
    r_ = s.r;
    v_ = s.v;
    t_ = t;
    keep_converged();
    // failures: a line lets go at its time or when an end segment reaches its breaking tension
    for(size_t k=0; k<lines_.size(); ++k)
    {
        const line_t &L = lines_[k];
        if(!L.has_break || L.broken >= 0)
            continue;
        int e = -1;
        if(t >= L.brk.t)
            e = L.brk.end >= 0 ? L.brk.end : 1;
        else
        {
            const double Ta = line_tension((int)k, 0), Tb = line_tension((int)k, L.N - 1);
            if(L.brk.end != 1 && Ta >= L.brk.T) e = 0;
            if(L.brk.end != 0 && Tb >= L.brk.T && (e < 0 || Tb > Ta)) e = 1;
        }
        if(e >= 0)
            break_line((int)k, e, t);
    }
    ++gen_;
    cvalid_ = false;
}

// ---- restart from a checkpoint: the start values of the catenary iteration

void mooring_qs::save_state(std::vector<double> &v) const
{
    state_out o{v};
    o.d(H_);
    o.d(V_);
}

bool mooring_qs::load_state(const std::vector<double> &v, size_t &k)
{
    state_in in{v, k};
    in.d(H_);
    in.d(V_);
    k = in.k;
    return in.ok;
}

}
