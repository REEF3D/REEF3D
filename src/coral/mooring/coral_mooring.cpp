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
    if(a.body < 0 && b.body < 0) return std::min(a.p(2), b.p(2));
    if(a.body < 0) return a.p(2);
    if(b.body < 0) return b.p(2);
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

// ---------------------------------------------------------------- lumped mass

mooring_line::mooring_line(const model &m, const mooring_end &a, const mooring_end &b, const mooring_props &p, int nseg,
                           std::shared_ptr<const wave_kinematics> waves, const Vec3 &current, double z_seabed, double dt_line,
                           const std::string &name)
    : name_(name), ba_(a.body), bb_(b.body), pr_(p), nseg_(std::max(nseg, 2)), waves_(std::move(waves)),
      current_(current), zbed_(z_seabed)
{
    sa_ = m.body(a.body).to_body(a.p);
    sb_ = m.body(b.body).to_body(b.p);
    g_ = m.gravity;
    l_ = p.L/nseg_;
    const double A = 0.25*M_PI*p.d*p.d;
    const double mn = (p.mass + p.rho*A*std::max(p.Can, p.Cat))*l_;
    const double mn0 = p.mass*l_;
    // explicit RK4: dt from the axial and seabed frequencies and the damping
    double dtc = 0.5*std::sqrt(mn0/(p.EA/l_));
    if(p.kb > 0.0) dtc = std::min(dtc, 0.5*std::sqrt(mn0/(p.kb*p.d*l_)));
    const double cmax = 2.0*p.BA/l_ + p.cb*p.d*l_;
    if(cmax > 0.0) dtc = std::min(dtc, 0.5*mn0/cmax);
    (void)mn;
    dt_ = dt_line > 0.0 ? dt_line : dtc;

    // initial shape: catenary between the lower and the upper end (seabed at the lower end), straight
    // without gravity or for a line that floats
    const Vec3 pa = a.p, pb = b.p;
    const bool a_low = pa(2) <= pb(2);
    const Vec3 lo = a_low ? pa : pb, up = a_low ? pb : pa;
    Vec3 dh = up - lo;
    const double zF = dh(2);
    dh(2) = 0.0;
    const double xF = dh.norm();
    const Vec3 eh = xF > 0.0 ? Vec3(dh/xF) : Vec3::UnitX();
    catenary_line cl;
    cl.L = p.L;
    cl.EA = p.EA;
    cl.CB = (zbed_ > lo(2) - 1.0e-6*p.L) ? std::max(p.CB, 0.0) : -1.0;
    cl.w = (p.mass - p.rho*A)*m.gravity.norm();     // no gravity: straight
    double H = -1.0, V = 0.0;
    const bool ok = cl.w > 0.0 && xF > 1.0e-6*p.L && catenary_solve(cl, xF, zF, H, V);
    r_.assign(nseg_+1, Vec3::Zero());
    v_.assign(nseg_+1, Vec3::Zero());
    // the profile's horizontal extent is scaled to the span (exact for a taut catenary; a slack line
    // is spread over the span with slack segments)
    double xend = xF, zend = zF;
    if(ok)
        catenary_profile(cl, H, V, p.L, xend, zend);
    const double sx = (ok && xend > 0.0) ? xF/xend : 1.0;
    for(int i=0; i<=nseg_; ++i)
    {
        const double s = i*l_;
        Vec3 q;
        if(ok)
        {
            double x, z;
            catenary_profile(cl, H, V, s, x, z);
            q = lo + sx*x*eh + z*Vec3::UnitZ();
        }
        else
            q = lo + (up - lo)*(s/p.L);
        r_[a_low ? i : nseg_ - i] = q;
    }
    r_[0] = pa;
    r_[nseg_] = pb;
    ra_ = pa;
    rb_ = pb;
}

const mooring_line::node_kin &mooring_line::kinematics(int slot, double t) const
{
    node_kin &k = kin_[slot];
    if(k.gen==gen_ && k.t==t)
        return k;
    const int n = nseg_ + 1;
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

void mooring_line::water(int i, const Vec3 &x, double t, double &eta, Vec3 &u, Vec3 &du) const
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

double mooring_line::tension(int j) const
{
    return j >= 0 && j < (int)Tc_.size() ? Tc_[j].norm() : 0.0;
}

void mooring_line::keep_converged()
{
    Tc_ = T_;
    bool ok = Fa_.allFinite() && Fb_.allFinite();
    for(int i=0; i<=nseg_ && ok; ++i)
        ok = r_[i].allFinite() && v_[i].allFinite();
    if(!ok)
        throw std::runtime_error("coral::mooring_line " + name_ + ": the line state is not finite (the explicit "
                                 "scheme became unstable): set a smaller dt_line (now " + std::to_string(dt_) + " s)");
}

void mooring_line::accel(const state &s, const Vec3 &ra, const Vec3 &rb, const Vec3 &va, const Vec3 &vb, double t,
                         std::vector<Vec3> &a, Vec3 &Fa, Vec3 &Fb, double extra) const
{
    const int N = nseg_;
    const mooring_props &p = pr_;
    const double A = 0.25*M_PI*p.d*p.d;
    const Vec3 &g = g_;
    auto R = [&](int i) -> const Vec3& {return i==0 ? ra : (i==N ? rb : s.r[i]);};
    auto V = [&](int i) -> const Vec3& {return i==0 ? va : (i==N ? vb : s.v[i]);};
    // segment tensions (no compression)
    T_.resize(N);
    for(int j=0; j<N; ++j)
    {
        const Vec3 d = R(j+1) - R(j);
        const double L = d.norm();
        double Tj = 0.0;
        if(L > l_)
        {
            const Vec3 e = d/L;
            const double sr = (V(j+1) - V(j)).dot(e)/l_;
            Tj = std::max(0.0, p.EA*(L/l_ - 1.0) + p.BA*sr)/L;
        }
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
        water(node, x, t, eta, u, aw);
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
    a.resize(N+1);
    a[0].setZero();
    a[N].setZero();
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
        a[i] = (F - (be/(al + be))*fq*q)/al;
    }
    // end forces: tension of the end segment, weight and hydrodynamics of the end half-segment; its mass
    // and added mass (carried by the attached body, end_inertia)
    Vec3 qa = R(1) - R(0), qb = R(N) - R(N-1);
    qa = qa.norm() > 0.0 ? Vec3(qa/qa.norm()) : Vec3::UnitZ();
    qb = qb.norm() > 0.0 ? Vec3(qb/qb.norm()) : Vec3::UnitZ();
    const double me = 0.5*mi;
    Fa = T_[0] + external(0, ra, va, qa, 0.5*l_, an, at);
    Ma_ = me*Mat3::Identity() + an*(Mat3::Identity() - qa*qa.transpose()) + at*qa*qa.transpose();
    Fb = -T_[N-1] + external(N, rb, vb, qb, 0.5*l_, an, at);
    Mb_ = me*Mat3::Identity() + an*(Mat3::Identity() - qb*qb.transpose()) + at*qb*qb.transpose();
}

void mooring_line::end_inertia(Vec3 &Fa, Vec3 &Fb, const Vec3 &aa, const Vec3 &ab) const
{
    Fa -= Ma_*aa;
    Fb -= Mb_*ab;
}

void mooring_line::integrate(state &s, const end_motion &ea, const end_motion &eb, double t0, double t1, Vec3 &Fa, Vec3 &Fb,
                             double extra, bool at_rest) const
{
    const int N = nseg_;
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
    std::vector<Vec3> rs(N+1), vs(N+1);
    // the ends on cubic Hermite curves through the converged and the trial positions and velocities
    auto hermite = [&](const end_motion &e, double tt, Vec3 &r, Vec3 &v)
    {
        if(T <= 0.0)
        {
            r = e.r1;
            v = e.v1;
            return;
        }
        const double x = (tt - t0)/T, x2 = x*x, x3 = x2*x;
        r = (2.0*x3 - 3.0*x2 + 1.0)*e.r0 + T*(x3 - 2.0*x2 + x)*e.v0 + (3.0*x2 - 2.0*x3)*e.r1 + T*(x3 - x2)*e.v1;
        v = (6.0*(x2 - x)/T)*(e.r0 - e.r1) + (3.0*x2 - 4.0*x + 1.0)*e.v0 + (3.0*x2 - 2.0*x)*e.v1;
    };
    Vec3 ra, rb, va, vb;
    auto ends = [&](double tt)
    {
        hermite(ea, tt, ra, va);
        hermite(eb, tt, rb, vb);
    };
    for(int k=0; k<n; ++k)
    {
        const double tk = t0 + k*h;
        // RK4 on (r, v) of the inner nodes; the stage velocities are the stage states' v
        ends(tk);
        accel(s, ra, rb, va, vb, tk, k1, Fa, Fb, extra);
        for(int i=1; i<N; ++i) {s2.r[i] = s.r[i] + 0.5*h*s.v[i]; s2.v[i] = s.v[i] + 0.5*h*k1[i];}
        rs = s.v;
        ends(tk + 0.5*h);
        accel(s2, ra, rb, va, vb, tk + 0.5*h, k2, Fa, Fb, extra);
        vs = s2.v;
        for(int i=1; i<N; ++i) {s2.r[i] = s.r[i] + 0.5*h*vs[i]; s2.v[i] = s.v[i] + 0.5*h*k2[i];}
        for(int i=1; i<N; ++i) rs[i] += 2.0*vs[i];
        accel(s2, ra, rb, va, vb, tk + 0.5*h, k3, Fa, Fb, extra);
        vs = s2.v;
        for(int i=1; i<N; ++i) {s2.r[i] = s.r[i] + h*vs[i]; s2.v[i] = s.v[i] + h*k3[i];}
        for(int i=1; i<N; ++i) rs[i] += 2.0*vs[i];
        ends(tk + h);
        accel(s2, ra, rb, va, vb, tk + h, k4, Fa, Fb, extra);
        for(int i=1; i<N; ++i)
        {
            s.r[i] += h/6.0*(rs[i] + s2.v[i]);
            s.v[i] += h/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        }
    }
    s.r[0] = ea.r1; s.r[N] = eb.r1;
    s.v[0] = ea.v1; s.v[N] = eb.v1;
    accel(s, ea.r1, eb.r1, ea.v1, eb.v1, t1, k1, Fa, Fb, extra);
}

void mooring_line::initialize(const model &m, double t0)
{
    if(static_start_ && !started_)
        equilibrium(m, t0);
    started_ = true;
    t_ = t0;
    const rigid_body &A = m.body(ba_), &B = m.body(bb_);
    ra_ = A.point(sa_);
    rb_ = B.point(sb_);
    va_ = A.point_velocity(sa_);
    vb_ = B.point_velocity(sb_);
    r_[0] = ra_;
    r_[nseg_] = rb_;
    v_[0] = va_;
    v_[nseg_] = vb_;
    ++gen_;
    state s{r_, v_};
    std::vector<Vec3> a;
    if(waves_ && !exact_kin_)
    {
        kinematics(0, t0);
        kinematics(1, t0);
        kt0_ = kt1_ = t0;
        kmode_ = 2;
    }
    else
        kmode_ = waves_ ? 1 : 0;
    accel(s, ra_, rb_, va_, vb_, t0, a, Fa_, Fb_, 0.0);
    end_inertia(Fa_, Fb_, A.point_acceleration(sa_), B.point_acceleration(sb_));
    keep_converged();
    cvalid_ = false;
}

void mooring_line::settle(const model &m, double T)
{
    ra_ = m.body(ba_).point(sa_);
    rb_ = m.body(bb_).point(sb_);
    r_[0] = ra_;
    r_[nseg_] = rb_;
    ++gen_;
    state s{r_, v_};
    // critical-type damping of the slowest pendulum mode of the line, about 1/s; the water at rest (the
    // surface at the start time t_)
    const end_motion ea{ra_, Vec3::Zero(), ra_, Vec3::Zero()}, eb{rb_, Vec3::Zero(), rb_, Vec3::Zero()};
    integrate(s, ea, eb, t_, t_ + T, Fa_, Fb_, 2.0, true);
    for(auto &x : s.v) x.setZero();
    va_.setZero();
    vb_.setZero();
    r_ = s.r;
    v_ = s.v;
    keep_converged();
    ++gen_;
    cvalid_ = false;
}

void mooring_line::eval(const model &m, double t, VecX &f) const
{
    const rigid_body &A = m.body(ba_), &B = m.body(bb_);
    const Vec3 ra1 = A.point(sa_), rb1 = B.point(sb_);
    const Vec3 va1 = A.point_velocity(sa_), vb1 = B.point_velocity(sb_);
    if(!(cvalid_ && t==ct_ && ra1==cra_ && rb1==crb_ && va1==cva_ && vb1==cvb_))
    {
        state s{r_, v_};
        integrate(s, end_motion{ra_, va_, ra1, va1}, end_motion{rb_, vb_, rb1, vb1}, t_, t, cFa_, cFb_, 0.0);
        cMa_ = Ma_;
        cMb_ = Mb_;
        cra_ = ra1;
        crb_ = rb1;
        cva_ = va1;
        cvb_ = vb1;
        ct_ = t;
        cvalid_ = true;
    }
    // the body carries the end node: minus its mass and added mass times the attachment's acceleration
    m.add_force(f, ba_, ra1, cFa_ - cMa_*A.point_acceleration(sa_));
    m.add_force(f, bb_, rb1, cFb_ - cMb_*B.point_acceleration(sb_));
}

void mooring_line::kinematics_points(const model&, const void *source, std::vector<Vec3> &pts) const
{
    // the converged nodes 0 .. N
    if(waves_.get() != source)
        return;
    pts.insert(pts.end(), r_.begin(), r_.end());
}

void mooring_line::commit(const model &m, double t, double)
{
    const rigid_body &A = m.body(ba_), &B = m.body(bb_);
    const Vec3 ra1 = A.point(sa_), rb1 = B.point(sb_);
    const Vec3 va1 = A.point_velocity(sa_), vb1 = B.point_velocity(sb_);
    state s{r_, v_};
    integrate(s, end_motion{ra_, va_, ra1, va1}, end_motion{rb_, vb_, rb1, vb1}, t_, t, Fa_, Fb_, 0.0);
    end_inertia(Fa_, Fb_, A.point_acceleration(sa_), B.point_acceleration(sb_));
    r_ = s.r;
    v_ = s.v;
    ra_ = ra1;
    rb_ = rb1;
    va_ = va1;
    vb_ = vb1;
    t_ = t;
    keep_converged();
    ++gen_;
    cvalid_ = false;
}

void mooring_line::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"Fa_x", "Fa_y", "Fa_z", "Fb_x", "Fb_y", "Fb_z", "Tmax", "Ta", "Tb"})
        c.push_back(name_ + "_" + x);
}

void mooring_line::channel_values(const model&, double, std::vector<double> &v) const
{
    for(int i=0; i<3; ++i) v.push_back(Fa_(i));
    for(int i=0; i<3; ++i) v.push_back(Fb_(i));
    double tm = 0.0;
    for(int j=0; j<nseg_; ++j) tm = std::max(tm, tension(j));
    v.push_back(tm);
    v.push_back(tension(0));
    v.push_back(tension(nseg_ - 1));
}

// ---- restart from a checkpoint: the start values of the catenary iteration; the converged line

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

void mooring_line::save_state(std::vector<double> &v) const
{
    state_out o{v};
    o.list(r_, [&](const Vec3 &x){o.m(x);});
    o.list(v_, [&](const Vec3 &x){o.m(x);});
    o.list(Tc_, [&](const Vec3 &x){o.m(x);});
    o.m(ra_);
    o.m(rb_);
    o.d(t_);
    o.m(Fa_);
    o.m(Fb_);
    o.m(va_);
    o.m(vb_);
}

bool mooring_line::load_state(const std::vector<double> &v, size_t &k)
{
    state_in in{v, k};
    in.list(r_, [&](Vec3 &x){in.m(x);});
    in.list(v_, [&](Vec3 &x){in.m(x);});
    in.list(Tc_, [&](Vec3 &x){in.m(x);});
    in.m(ra_);
    in.m(rb_);
    in.d(t_);
    in.m(Fa_);
    in.m(Fb_);
    in.m(va_);
    in.m(vb_);
    started_ = true;
    ++gen_;
    cvalid_ = false;
    k = in.k;
    return in.ok;
}

}
