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

// REEF3D::CORAL: mooring systems: set-up (lines, points, start shapes, sub-step), start, settling, outputs and
// checkpoints (the dynamics in coral_mooring.cpp, the static equilibrium in coral_mooring_static.cpp)

#include"coral/mooring/coral_mooring.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>
#include<stdexcept>

namespace coral
{

mooring_system::mooring_system(const model &m, std::shared_ptr<const wave_kinematics> waves, const Vec3 &current,
                               const std::string &name)
    : name_(name), waves_(std::move(waves)), current_(current), g_(m.gravity)
{
}

int mooring_system::find_point(const std::string &n) const
{
    for(size_t j=0; j<points_.size(); ++j)
        if(points_[j].name==n)
            return (int)j;
    return -1;
}

// the node layout: the nodes of the lines (N + 1 each) in their order, then the points; the attachments
// of the body ends; the nodes the scheme advances
void mooring_system::layout(const std::vector<Vec3> &line_r, const std::vector<Vec3> &point_r)
{
    int n = 0;
    for(line_t &L : lines_)
    {
        L.off = n;
        n += L.N + 1;
    }
    for(point_t &P : points_)
        P.idx = n++;
    r_ = line_r;
    r_.insert(r_.end(), point_r.begin(), point_r.end());
    v_.assign(r_.size(), Vec3::Zero());
    free_.clear();
    for(const line_t &L : lines_)
        for(int i=1; i<L.N; ++i)
            free_.push_back(L.off + i);
    for(const point_t &P : points_)
        free_.push_back(P.idx);
    for(attach_t &a : att_)
        a.idx = lines_[a.line].off + (a.end==0 ? 0 : lines_[a.line].N);
    for(point_t &P : points_)
        P.ends.clear();
    for(size_t k=0; k<lines_.size(); ++k)
        for(int e=0; e<2; ++e)
            if(lines_[k].pt[e] >= 0)
                points_[lines_[k].pt[e]].ends.push_back({(int)k, e});
    ++gen_;
    cvalid_ = false;
}

int mooring_system::add_point(const std::string &name, const Vec3 &p, const mooring_point_props &pp, double z_seabed)
{
    std::vector<Vec3> lr(r_.begin(), r_.begin() + (r_.size() - points_.size())), pr;
    for(const point_t &P : points_)
        pr.push_back(r_[P.idx]);
    point_t P;
    P.name = name;
    P.p = pp;
    P.zbed = z_seabed;
    points_.push_back(P);
    pr.push_back(p);
    layout(lr, pr);
    update_dt();
    return (int)points_.size() - 1;
}

// the elastic catenary of line L between pa and pb (seabed portion at the lower end if it lies on the line's
// seabed, else hanging free): forces of the line on its ends and the nodes along it; straight without
// gravity or for a line that floats (false: no catenary, the forces zero)
bool mooring_system::line_catenary(const line_t &L, const Vec3 &pa, const Vec3 &pb, double gm, Vec3 &Fa, Vec3 &Fb,
                                   std::vector<Vec3> *nodes) const
{
    const mooring_props &p = L.p;
    const int N = L.N;
    const double l_ = L.l;
    const double A = 0.25*M_PI*p.d*p.d;
    // a floating line (buoyancy above its weight) is the catenary of its net buoyancy upside down
    const double wn = (p.mass - p.rho*A)*gm;
    if(wn < 0.0 && gm > 0.0)
    {
        line_t M = L;
        M.p.mass = 2.0*p.rho*A - p.mass;    // the mirrored line weighs what the line lifts
        M.zbed = -1.0e30;
        const Vec3 fa(pa(0), pa(1), -pa(2)), fb(pb(0), pb(1), -pb(2));
        const bool ok = line_catenary(M, fa, fb, gm, Fa, Fb, nodes);
        Fa(2) = -Fa(2);
        Fb(2) = -Fb(2);
        if(nodes)
            for(Vec3 &x : *nodes) x(2) = -x(2);
        return ok;
    }
    // initial shape: catenary between the lower and the upper end (seabed at the lower end), straight
    // without gravity
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
    cl.CB = (L.zbed > lo(2) - 1.0e-6*p.L) ? std::max(p.CB, 0.0) : -1.0;
    cl.w = wn;     // no gravity: straight
    double H = -1.0, V = 0.0;
    const bool ok = cl.w > 0.0 && xF > 1.0e-6*p.L && catenary_solve(cl, xF, zF, H, V);
    Fa.setZero();
    Fb.setZero();
    if(ok)
    {
        double HA, VA;
        catenary_anchor(cl, H, V, HA, VA);
        const Vec3 Fup = -H*eh - V*Vec3::UnitZ(), Flo = HA*eh + VA*Vec3::UnitZ();
        Fa = a_low ? Flo : Fup;
        Fb = a_low ? Fup : Flo;
    }
    if(!nodes)
        return ok;
    std::vector<Vec3> &rl = *nodes;
    rl.assign(N+1, Vec3::Zero());
    // the profile's horizontal extent is scaled to the span (exact for a taut catenary; a slack line
    // is spread over the span with slack segments)
    double xend = xF, zend = zF;
    if(ok)
        catenary_profile(cl, H, V, p.L, xend, zend);
    const double sx = (ok && xend > 0.0) ? xF/xend : 1.0;
    for(int i=0; i<=N; ++i)
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
        rl[a_low ? i : N - i] = q;
    }
    rl[0] = pa;
    rl[N] = pb;
    return ok;
}

std::vector<Vec3> mooring_junctions(const Vec3 &pa, const Vec3 &pb, const std::vector<mooring_props> &sec, double z_seabed, double g)
{
    // one catenary of the whole line: the mean weight per length, the series axial stiffness
    double Lt = 0.0, wl = 0.0, c = 0.0;
    for(const mooring_props &p : sec)
    {
        Lt += p.L;
        wl += (p.mass - p.rho*0.25*M_PI*p.d*p.d)*g*p.L;
        c += p.L/p.EA;
    }
    std::vector<double> sa;          // arc length of the junctions from end a
    double acc = 0.0;
    for(size_t k=0; k+1<sec.size(); ++k)
    {
        acc += sec[k].L;
        sa.push_back(acc);
    }
    const bool a_low = pa(2) <= pb(2);
    const Vec3 lo = a_low ? pa : pb, up = a_low ? pb : pa;
    Vec3 dh = up - lo;
    const double zF = dh(2);
    dh(2) = 0.0;
    const double xF = dh.norm();
    const Vec3 eh = xF > 0.0 ? Vec3(dh/xF) : Vec3::UnitX();
    catenary_line cl;
    cl.L = Lt;
    cl.EA = Lt/c;
    cl.CB = (z_seabed > lo(2) - 1.0e-6*Lt) ? 0.0 : -1.0;
    cl.w = wl/Lt;
    double H = -1.0, V = 0.0;
    const bool ok = cl.w > 0.0 && xF > 1.0e-6*Lt && catenary_solve(cl, xF, zF, H, V);
    double xend = xF, zend = zF;
    if(ok)
        catenary_profile(cl, H, V, Lt, xend, zend);
    const double sx = (ok && xend > 0.0) ? xF/xend : 1.0;
    std::vector<Vec3> J;
    for(double s : sa)
    {
        const double sl = a_low ? s : Lt - s;
        if(ok)
        {
            double x, z;
            catenary_profile(cl, H, V, sl, x, z);
            J.push_back(lo + sx*x*eh + z*Vec3::UnitZ());
        }
        else
            J.push_back(pa + (pb - pa)*(s/Lt));
    }
    return J;
}

int mooring_system::add_line(const model &m, const std::string &name, const mooring_end &a, const mooring_end &b,
                             const mooring_props &p, int nseg, double z_seabed)
{
    line_t L;
    L.name = name;
    L.p = p;
    L.N = std::max(nseg, 2);
    L.l = p.L/L.N;
    L.zbed = z_seabed;
    const double l_ = L.l;
    const double mn0 = p.mass*l_;
    // explicit RK4: dt from the axial and seabed frequencies and the damping
    double dtc = 0.5*std::sqrt(mn0/(p.EA/l_));
    if(p.kb > 0.0) dtc = std::min(dtc, 0.5*std::sqrt(mn0/(p.kb*p.d*l_)));
    const double cmax = 2.0*p.BA/l_ + p.cb*p.d*l_;
    if(cmax > 0.0) dtc = std::min(dtc, 0.5*mn0/cmax);
    L.dt = dtc;
    // the ends: global positions
    const mooring_end *ends[2] = {&a, &b};
    Vec3 pe[2];
    for(int e=0; e<2; ++e)
    {
        const mooring_end &E = *ends[e];
        if(E.point >= 0)
        {
            if(E.point >= (int)points_.size())
                throw std::runtime_error("coral::mooring_system: line " + name + " ends on an unknown point");
            L.pt[e] = E.point;
            pe[e] = r_[points_[E.point].idx];
        }
        else
            pe[e] = E.p;
    }
    std::vector<Vec3> rl;
    Vec3 Fa, Fb;
    line_catenary(L, pe[0], pe[1], m.gravity.norm(), Fa, Fb, &rl);
    // the new layout: this line's nodes after the other lines'
    std::vector<Vec3> lr(r_.begin(), r_.begin() + (r_.size() - points_.size())), pr;
    for(const point_t &P : points_)
        pr.push_back(r_[P.idx]);
    lr.insert(lr.end(), rl.begin(), rl.end());
    lines_.push_back(L);
    const int k = (int)lines_.size() - 1;
    for(int e=0; e<2; ++e)
        if(lines_[k].pt[e] < 0)
        {
            attach_t at;
            at.body = ends[e]->body;
            at.s = m.body(at.body).to_body(pe[e]);
            at.line = k;
            at.end = e;
            at.r = pe[e];
            lines_[k].att[e] = (int)att_.size();
            att_.push_back(at);
        }
    layout(lr, pr);
    update_dt();
    return k;
}

void mooring_system::set_dt_line(double dt)
{
    dt_set_ = dt;
    update_dt();
}

void mooring_system::update_dt()
{
    double dt = 1.0e300;
    for(const line_t &L : lines_)
        dt = std::min(dt, L.dt);
    // points on the seabed: the contact frequency of the point with the end nodes of its lines
    for(const point_t &P : points_)
    {
        const mooring_point_props &p = P.p;
        double M = p.mass;
        for(const auto &e : P.ends)
            M += 0.5*lines_[e.first].p.mass*lines_[e.first].l;
        const double rv = std::cbrt(0.75*p.volume/M_PI);
        const double rc = p.radius >= 0.0 ? p.radius : std::max(rv, 0.05);
        const double Ac = M_PI*rc*rc;
        if(p.kb > 0.0 && M > 0.0)
            dt = std::min(dt, 0.5*std::sqrt(M/(p.kb*Ac)));
        if(p.cb > 0.0 && M > 0.0)
            dt = std::min(dt, 0.5*M/(p.cb*Ac));
    }
    dt_ = dt_set_ > 0.0 ? dt_set_ : dt;
}

void mooring_system::coupled_bodies(std::vector<int> &b) const
{
    for(const attach_t &a : att_)
        b.push_back(a.body);
}

std::vector<Vec3> mooring_system::line_nodes(int k) const
{
    const line_t &L = lines_[k];
    return std::vector<Vec3>(r_.begin() + L.off, r_.begin() + L.off + L.N + 1);
}

double mooring_system::line_tension(int k, int j) const
{
    const line_t &L = lines_[k];
    return j >= 0 && j < (int)L.Tc.size() ? L.Tc[j].norm() : 0.0;
}

void mooring_system::keep_converged()
{
    bool ok = true;
    for(line_t &L : lines_)
    {
        L.Tc = L.T;
        for(int e=0; e<2; ++e)
        {
            L.Fc[e] = L.att[e] >= 0 ? att_[L.att[e]].F : L.F[e];
            ok = ok && L.Fc[e].allFinite();
        }
    }
    for(size_t i=0; i<r_.size() && ok; ++i)
        ok = r_[i].allFinite() && v_[i].allFinite();
    if(!ok)
    {
        std::string n = name_;
        for(const line_t &L : lines_)
        {
            bool f = L.Fc[0].allFinite() && L.Fc[1].allFinite();
            for(int i=0; i<=L.N && f; ++i)
                f = r_[L.off+i].allFinite() && v_[L.off+i].allFinite();
            if(!f) {n = L.name; break;}
        }
        throw std::runtime_error("coral::mooring_line " + n + ": the line state is not finite (the explicit "
                                 "scheme became unstable): set a smaller dt_line (now " + std::to_string(dt_) + " s)");
    }
}

void mooring_system::initialize(const model &m, double t0)
{
    if(static_start_ && !started_)
        equilibrium(m, t0);
    started_ = true;
    t_ = t0;
    for(attach_t &a : att_)
    {
        const rigid_body &B = m.body(a.body);
        a.r = B.point(a.s);
        a.v = B.point_velocity(a.s);
        r_[a.idx] = a.r;
        v_[a.idx] = a.v;
    }
    {
        state p{r_, v_};
        point_ends(p);
        r_ = p.r;
        v_ = p.v;
    }
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
    accel(s, t0, a, 0.0);
    for(attach_t &at : att_)
    {
        const line_t &L = lines_[at.line];
        at.F = L.F[at.end] - L.M[at.end]*m.body(at.body).point_acceleration(at.s);
    }
    keep_converged();
    cvalid_ = false;
}

void mooring_system::settle(const model &m, double T)
{
    std::vector<end_motion> em;
    for(attach_t &a : att_)
    {
        a.r = m.body(a.body).point(a.s);
        r_[a.idx] = a.r;
        em.push_back(end_motion{a.r, Vec3::Zero(), a.r, Vec3::Zero()});
    }
    ++gen_;
    state s{r_, v_};
    // critical-type damping of the slowest pendulum mode of the line, about 1/s; the water at rest (the
    // surface at the start time t_)
    integrate(s, em, t_, t_ + T, 2.0, true);
    for(auto &x : s.v) x.setZero();
    for(attach_t &a : att_)
    {
        a.v.setZero();
        a.F = lines_[a.line].F[a.end];
    }
    r_ = s.r;
    v_ = s.v;
    keep_converged();
    ++gen_;
    cvalid_ = false;
}

void mooring_system::channel_names(std::vector<std::string> &c) const
{
    for(const line_t &L : lines_)
        for(const char *x : {"Fa_x", "Fa_y", "Fa_z", "Fb_x", "Fb_y", "Fb_z", "Tmax", "Ta", "Tb"})
            c.push_back(L.name + "_" + x);
    for(const point_t &P : points_)
        for(const char *x : {"x", "y", "z"})
            c.push_back(P.name + "_" + x);
}

void mooring_system::channel_values(const model&, double, std::vector<double> &v) const
{
    for(size_t k=0; k<lines_.size(); ++k)
    {
        const line_t &L = lines_[k];
        for(int i=0; i<3; ++i) v.push_back(L.Fc[0](i));
        for(int i=0; i<3; ++i) v.push_back(L.Fc[1](i));
        double tm = 0.0;
        for(int j=0; j<L.N; ++j) tm = std::max(tm, line_tension((int)k, j));
        v.push_back(tm);
        v.push_back(line_tension((int)k, 0));
        v.push_back(line_tension((int)k, L.N - 1));
    }
    for(const point_t &P : points_)
        for(int i=0; i<3; ++i) v.push_back(r_[P.idx](i));
}

// ---- restart from a checkpoint: the converged nodes, tensions and end forces, the attachments

void mooring_system::save_state(std::vector<double> &v) const
{
    state_out o{v};
    o.list(r_, [&](const Vec3 &x){o.m(x);});
    o.list(v_, [&](const Vec3 &x){o.m(x);});
    o.d(t_);
    o.list(lines_, [&](const line_t &L){o.list(L.Tc, [&](const Vec3 &x){o.m(x);}); o.m(L.Fc[0]); o.m(L.Fc[1]);});
    o.list(att_, [&](const attach_t &a){o.m(a.r); o.m(a.v); o.m(a.F);});
}

bool mooring_system::load_state(const std::vector<double> &v, size_t &k)
{
    state_in in{v, k};
    std::vector<Vec3> r, vv;
    in.list(r, [&](Vec3 &x){in.m(x);});
    in.list(vv, [&](Vec3 &x){in.m(x);});
    in.d(t_);
    if(r.size() != r_.size() || vv.size() != v_.size())
        in.ok = false;
    else
    {
        r_ = r;
        v_ = vv;
    }
    size_t nl = 0;
    in.i(nl);
    if(nl != lines_.size())
        in.ok = false;
    for(size_t q=0; q<nl && in.ok; ++q)
    {
        line_t &L = lines_[q];
        in.list(L.Tc, [&](Vec3 &x){in.m(x);});
        in.m(L.Fc[0]);
        in.m(L.Fc[1]);
    }
    size_t na = 0;
    in.i(na);
    if(na != att_.size())
        in.ok = false;
    for(size_t q=0; q<na && in.ok; ++q)
    {
        in.m(att_[q].r);
        in.m(att_[q].v);
        in.m(att_[q].F);
    }
    started_ = true;
    ++gen_;
    cvalid_ = false;
    k = in.k;
    return in.ok;
}

// ---------------------------------------------------------------- a single line

mooring_line::mooring_line(const model &m, const mooring_end &a, const mooring_end &b, const mooring_props &p, int nseg,
                           std::shared_ptr<const wave_kinematics> waves, const Vec3 &current, double z_seabed, double dt_line,
                           const std::string &name)
    : mooring_system(m, std::move(waves), current, name)
{
    add_line(m, name, a, b, p, nseg, z_seabed);
    set_dt_line(dt_line);
}

}
