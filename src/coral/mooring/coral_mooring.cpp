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

// ---------------------------------------------------------------- quasi-static

mooring_qs::mooring_qs(const model &m, const mooring_end &a, const mooring_end &b, const mooring_props &p, const std::string &name)
    : name_(name), ba_(a.body), bb_(b.body), pr_(p)
{
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
    double H = H_, V = V_;
    if(!catenary_solve(cl_, xF, zF, H, V))
    {
        H = -1.0;
        if(!catenary_solve(cl_, xF, zF, H, V))
            throw std::runtime_error("coral::mooring_qs: no catenary solution for line " + name_);
    }
    H_ = H;
    V_ = V;
    double HA, VA;
    catenary_anchor(cl_, H, V, HA, VA);
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
                           std::shared_ptr<const wave_field> waves, const Vec3 &current, double z_seabed, double dt_line,
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

    // initial shape: catenary between the lower and the upper end (seabed at the lower end)
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
    cl.w = (p.mass - p.rho*A)*gmag(m);
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

double mooring_line::tension(int j) const
{
    const double L = (r_[j+1] - r_[j]).norm();
    return L > l_ ? pr_.EA*(L/l_ - 1.0) : 0.0;
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
    auto external = [&](const Vec3 &x, const Vec3 &vx, const Vec3 &q, double ll, double &an, double &at)
    {
        Vec3 F = p.mass*ll*g;
        an = at = 0.0;
        const double eta = waves_ ? waves_->elevation(x(0), x(1), t) : 0.0;
        if(x(2) > eta)
            return F;
        F -= p.rho*A*ll*g;
        Vec3 uw = current_, aw = Vec3::Zero();
        if(waves_)
        {
            Vec3 u, du; double pr;
            waves_->kinematics(x, t, u, du, pr);
            uw += u;
            aw = du;
        }
        const Vec3 vr = uw - vx;
        const double vq = vr.dot(q);
        const Vec3 vrt = vq*q, vrn = vr - vrt;
        F += 0.5*p.rho*p.d*ll*p.Cdn*vrn.norm()*vrn + 0.5*p.rho*M_PI*p.d*ll*p.Cdt*std::fabs(vq)*vrt;
        if(waves_)
        {
            const Vec3 awt = aw.dot(q)*q, awn = aw - awt;
            F += p.rho*A*ll*((1.0 + p.Can)*awn + (1.0 + p.Cat)*awt);
        }
        an = p.rho*A*ll*p.Can;
        at = p.rho*A*ll*p.Cat;
        // seabed
        if(x(2) < zbed_)
        {
            const double Fz = std::max(0.0, (p.kb*(zbed_ - x(2)) - p.cb*vx(2))*p.d*ll);
            F(2) += Fz;
            if(p.mu > 0.0)
            {
                const double vn = std::sqrt(vx(0)*vx(0) + vx(1)*vx(1));
                if(vn > 0.0)
                {
                    const double c = p.mu*Fz*std::tanh(vn/p.v_eps)/vn;
                    F(0) -= c*vx(0);
                    F(1) -= c*vx(1);
                }
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
        const Vec3 F = T_[i] - T_[i-1] + external(R(i), V(i), q, l_, an, at) - extra*mi*V(i);
        // (m I + an (I - q q^T) + at q q^T)^-1 F
        const double al = mi + an, be = at - an;
        const double fq = F.dot(q);
        a[i] = (F - (be/(al + be))*fq*q)/al;
    }
    // end forces: tension of the end segment, weight and hydrodynamics of the end half-segment
    Vec3 qa = R(1) - R(0), qb = R(N) - R(N-1);
    qa = qa.norm() > 0.0 ? Vec3(qa/qa.norm()) : Vec3::UnitZ();
    qb = qb.norm() > 0.0 ? Vec3(qb/qb.norm()) : Vec3::UnitZ();
    Fa = T_[0] + external(ra, va, qa, 0.5*l_, an, at);
    Fb = -T_[N-1] + external(rb, vb, qb, 0.5*l_, an, at);
}

void mooring_line::integrate(state &s, const Vec3 &ra0, const Vec3 &ra1, const Vec3 &rb0, const Vec3 &rb1,
                             double t0, double t1, Vec3 &Fa, Vec3 &Fb, double extra) const
{
    const int N = nseg_;
    const double T = t1 - t0;
    const int n = std::max(1, (int)std::ceil(T/dt_ - 1.0e-9));
    const double h = T/n;
    const Vec3 va = T > 0.0 ? Vec3((ra1 - ra0)/T) : Vec3::Zero();
    const Vec3 vb = T > 0.0 ? Vec3((rb1 - rb0)/T) : Vec3::Zero();
    std::vector<Vec3> &k1 = w_[0], &k2 = w_[1], &k3 = w_[2], &k4 = w_[3];
    state &s2 = ws_;
    s2 = s;
    std::vector<Vec3> rs(N+1), vs(N+1);
    auto ends = [&](double tt, Vec3 &ra, Vec3 &rb)
    {
        const double f = T > 0.0 ? (tt - t0)/T : 1.0;
        ra = ra0 + f*(ra1 - ra0);
        rb = rb0 + f*(rb1 - rb0);
    };
    Vec3 ra, rb;
    for(int k=0; k<n; ++k)
    {
        const double tk = t0 + k*h;
        // RK4 on (r, v) of the inner nodes; the stage velocities are the stage states' v
        ends(tk, ra, rb);
        accel(s, ra, rb, va, vb, tk, k1, Fa, Fb, extra);
        for(int i=1; i<N; ++i) {s2.r[i] = s.r[i] + 0.5*h*s.v[i]; s2.v[i] = s.v[i] + 0.5*h*k1[i];}
        rs = s.v;
        ends(tk + 0.5*h, ra, rb);
        accel(s2, ra, rb, va, vb, tk + 0.5*h, k2, Fa, Fb, extra);
        vs = s2.v;
        for(int i=1; i<N; ++i) {s2.r[i] = s.r[i] + 0.5*h*vs[i]; s2.v[i] = s.v[i] + 0.5*h*k2[i];}
        for(int i=1; i<N; ++i) rs[i] += 2.0*vs[i];
        accel(s2, ra, rb, va, vb, tk + 0.5*h, k3, Fa, Fb, extra);
        vs = s2.v;
        for(int i=1; i<N; ++i) {s2.r[i] = s.r[i] + h*vs[i]; s2.v[i] = s.v[i] + h*k3[i];}
        for(int i=1; i<N; ++i) rs[i] += 2.0*vs[i];
        ends(tk + h, ra, rb);
        accel(s2, ra, rb, va, vb, tk + h, k4, Fa, Fb, extra);
        for(int i=1; i<N; ++i)
        {
            s.r[i] += h/6.0*(rs[i] + s2.v[i]);
            s.v[i] += h/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        }
    }
    s.r[0] = ra1; s.r[N] = rb1;
    s.v[0] = va; s.v[N] = vb;
    accel(s, ra1, rb1, va, vb, t1, k1, Fa, Fb, extra);
}

void mooring_line::initialize(const model &m, double t0)
{
    t_ = t0;
    ra_ = m.body(ba_).point(sa_);
    rb_ = m.body(bb_).point(sb_);
    r_[0] = ra_;
    r_[nseg_] = rb_;
    state s{r_, v_};
    std::vector<Vec3> a;
    accel(s, ra_, rb_, Vec3::Zero(), Vec3::Zero(), t0, a, Fa_, Fb_, 0.0);
    cvalid_ = false;
}

void mooring_line::settle(const model &m, double T)
{
    ra_ = m.body(ba_).point(sa_);
    rb_ = m.body(bb_).point(sb_);
    state s{r_, v_};
    // critical-type damping of the slowest pendulum mode of the line, about 1/s
    integrate(s, ra_, ra_, rb_, rb_, 0.0, T, Fa_, Fb_, 2.0);
    for(auto &x : s.v) x.setZero();
    r_ = s.r;
    v_ = s.v;
    cvalid_ = false;
}

void mooring_line::eval(const model &m, double t, VecX &f) const
{
    const Vec3 ra1 = m.body(ba_).point(sa_), rb1 = m.body(bb_).point(sb_);
    if(!(cvalid_ && t==ct_ && ra1==cra_ && rb1==crb_))
    {
        state s{r_, v_};
        integrate(s, ra_, ra1, rb_, rb1, t_, t, cFa_, cFb_, 0.0);
        cra_ = ra1;
        crb_ = rb1;
        ct_ = t;
        cvalid_ = true;
    }
    m.add_force(f, ba_, ra1, cFa_);
    m.add_force(f, bb_, rb1, cFb_);
}

void mooring_line::commit(const model &m, double t, double)
{
    const Vec3 ra1 = m.body(ba_).point(sa_), rb1 = m.body(bb_).point(sb_);
    state s{r_, v_};
    integrate(s, ra_, ra1, rb_, rb1, t_, t, Fa_, Fb_, 0.0);
    r_ = s.r;
    v_ = s.v;
    ra_ = ra1;
    rb_ = rb1;
    t_ = t;
    cvalid_ = false;
}

void mooring_line::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"Fa_x", "Fa_y", "Fa_z", "Fb_x", "Fb_y", "Fb_z", "Tmax"})
        c.push_back(name_ + "_" + x);
}

void mooring_line::channel_values(const model&, double, std::vector<double> &v) const
{
    for(int i=0; i<3; ++i) v.push_back(Fa_(i));
    for(int i=0; i<3; ++i) v.push_back(Fb_(i));
    double tm = 0.0;
    for(int j=0; j<nseg_; ++j) tm = std::max(tm, tension(j));
    v.push_back(tm);
}

}
