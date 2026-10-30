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
#include<limits>
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
    cl_.EA = p.EA_secant();
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

void mooring_qs::shape(const model &m, int n, std::vector<Vec3> &x, std::vector<double> &T) const
{
    // as forces(), with copies of the start values
    const Vec3 pa = m.body(ba_).point(sa_), pb = m.body(bb_).point(sb_);
    const bool a_low = pa(2) <= pb(2);
    const Vec3 lo = a_low ? pa : pb, up = a_low ? pb : pa;
    Vec3 dh = up - lo;
    const double zF = dh(2);
    dh(2) = 0.0;
    const double xF = dh.norm();
    const Vec3 eh = xF > 0.0 ? Vec3(dh/xF) : Vec3::UnitX();
    catenary_line cl = cl_;
    if(bed_end_ < 0 || (bed_end_==0) != a_low)
        cl.CB = -1.0;
    double H = H_, V = V_;
    x.clear();
    T.clear();
    if(!catenary_solve(cl, xF, zF, H, V))
    {
        H = -1.0;
        if(!catenary_solve(cl, xF, zF, H, V))
        {
            x = {pa, pb};
            T = {0.0, 0.0};
            return;
        }
    }
    n = std::max(n, 1);
    const double L = cl.L, w = cl.w;
    const double LB = (cl.CB < 0.0 && H > 0.0) ? 0.0 : std::max(0.0, L - V/w);
    for(int i=0; i<=n; ++i)
    {
        const double s = L*double(i)/double(n);
        double xs, zs;
        catenary_profile(cl, H, V, s, xs, zs);
        x.push_back(lo + xs*eh + zs*Vec3::UnitZ());
        // vertical tension: that of the fairlead less the weight of the line above s; on the seabed none, the
        // horizontal one reduced by the friction
        double Vs = V - w*(L - s), Hs = std::max(H, 0.0);
        if(Vs < 0.0 || s < LB)
        {
            Vs = 0.0;
            if(cl.CB > 0.0 && H > 0.0)
                Hs = std::max(0.0, H - cl.CB*w*(LB - s));
        }
        T.push_back(std::sqrt(Hs*Hs + Vs*Vs));
    }
    if(!a_low)
    {
        std::reverse(x.begin(), x.end());
        std::reverse(T.begin(), T.end());
    }
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
    if(fbed_)
    {
        k.bed.resize(n);
        for(int i=0; i<n; ++i)
            k.bed[i] = waves_->bed_at(this, i, r_[i]);
    }
    return k;
}

bool mooring_system::seabed_at(int i, const Vec3 &x, bool surf, bool fbed, double &zb, Vec3 &n) const
{
    if(fbed && waves_)
    {
        const node_kin &k = kin_[0];
        const double b = (k.gen==gen_ && (int)k.bed.size() > i) ? k.bed[i] : waves_->bed_at(this, i, x);
        if(std::isfinite(b))
        {
            zb = b;
            n = Vec3::UnitZ();
            return true;
        }
    }
    if(surf && sb_)
    {
        zb = sb_->z(x(0), x(1));
        n = sb_->normal(x(0), x(1));
        return true;
    }
    return false;
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

void mooring_system::line_accel(const line_t &L, const state &s, double t, std::vector<Vec3> &a, double extra, bool upd) const
{
    const int N = L.N, o = L.off;
    const mooring_props &p = L.p;
    const double A = 0.25*M_PI*p.d*p.d;
    const double l_ = L.l, zbed_ = L.zbed;
    const Vec3 &g = g_;
    auto R = [&](int i) -> const Vec3& {return s.r[o+i];};
    auto V = [&](int i) -> const Vec3& {return s.v[o+i];};
    std::vector<Vec3> &T_ = L.T;
    // segment tensions (no compression; with BA_slack the damping acts in every segment and may push)
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
            const double Te = p.elastic(Ls/l_ - 1.0);
            Tj = p.BA_slack ? (Te + p.BA*sr)/Ls : std::max(0.0, Te + p.BA*sr)/Ls;
        }
        else if(creg_ > 0.0 && Ls > 0.0)
            Tj = creg_*p.EA_secant()*(Ls/l_ - 1.0)/Ls;      // static solve only: a compressed segment pushes a little
        else if(p.BA_slack && Ls > 0.0)
            Tj = p.BA*(V(j+1) - V(j)).dot(d)/(Ls*Ls*l_);
        T_[j] = Tj*d;
    }
    // external forces on a node of length ll (weight, hydrodynamics, seabed); q tangent; added mass
    // alpha (normal) and beta (tangential); side: where the node's piece of line lies along q (0: centred on an
    // inner node, +1: from the node along q (end a), -1: against q (end b))
    const bool fric = p.mu > 0.0 || p.mu_axial > 0.0;
    auto external = [&](int node, const Vec3 &x, const Vec3 &vx, const Vec3 &q, double ll, double side, double &an, double &at)
    {
        Vec3 F = p.mass*ll*g;
        an = at = 0.0;
        double eta;
        Vec3 u, aw;
        water(o + node, x, t, eta, u, aw);
        // submerged fraction of the node's piece of line: linear over its vertical extent (at least d), so
        // a node crossing the free surface takes its buoyancy, drag and added mass gradually; the piece of an end
        // node reaches from the node into the line
        const double hz = std::max(p.d, std::fabs(q(2))*ll);
        const double zc = x(2) + 0.5*side*ll*q(2);
        const double fs = std::clamp((eta - zc)/hz + 0.5, 0.0, 1.0);
        if(fs > 0.0)
        {
            F -= fs*p.rho*A*ll*g;
            const Vec3 uw = current_at(x(2)) + u;
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
        // seabed: penalty and damping normal to it, friction (seabed_friction); a surface or the fluid's bed:
        // along its normal n, the penetration (zb - z) n_z, the friction in its tangent plane
        double zs;
        Vec3 nb;
        bool contact = false, fr = false;
        if((L.surf || L.fbed) && seabed_at(o + node, x, L.surf, L.fbed, zs, nb))
        {
            if(x(2) < zs)
            {
                const double Fn = std::max(0.0, (p.kb*(zs - x(2))*nb(2) - p.cb*vx.dot(nb))*p.d*ll);
                F += Fn*nb;
                if(fric)
                {
                    const Vec3 Ff = seabed_friction(L, o + node, x, vx, q, nb, Fn, ll, upd, fr, contact);
                    if(fr) F += Ff;
                }
            }
        }
        else if(x(2) < zbed_)
        {
            const double Fz = std::max(0.0, (p.kb*(zbed_ - x(2)) - p.cb*vx(2))*p.d*ll);
            F(2) += Fz;
            if(fric)
            {
                const Vec3 Ff = seabed_friction(L, o + node, x, vx, q, Vec3::UnitZ(), Fz, ll, upd, fr, contact);
                if(fr) F += Ff;
            }
        }
        // off the seabed: no anchor
        if(fric && upd && !contact && sc_)
            lose_anchor(o + node);
        return F;
    };
    // bending: energy EI/l (1 - cos phi) at each inner node (phi: the turning angle), EI/(l/2) (1 - cos phi)
    // between a clamped end's direction and its segment; forces -dE/dx, moments on the clamping bodies
    const bool bend = p.EI > 0.0;
    if(bend)
    {
        std::vector<Vec3> &B = L.B;
        B.assign(N+1, Vec3::Zero());
        const double k = p.EI/l_;
        for(int i=1; i<N; ++i)
        {
            const Vec3 a = R(i) - R(i-1), b = R(i+1) - R(i);
            const double la = a.norm(), lb = b.norm();
            if(la <= 0.0 || lb <= 0.0)
                continue;
            const Vec3 ta = a/la, tb = b/lb;
            const double c = ta.dot(tb);
            const Vec3 dca = (tb - c*ta)/la, dcb = (ta - c*tb)/lb;
            B[i-1] -= k*dca;
            B[i] += k*(dca - dcb);
            B[i+1] += k*dcb;
        }
        for(int e=0; e<2; ++e)
        {
            L.Mb[e].setZero();
            if(!L.clamp[e])
                continue;
            const int i0 = e==0 ? 0 : N, i1 = e==0 ? 1 : N-1;
            const Vec3 a = R(i1) - R(i0);
            const double la = a.norm();
            if(la <= 0.0)
                continue;
            const Vec3 ta = a/la, &d = L.cd[e];
            const double c = d.dot(ta), k0 = 2.0*k;
            const Vec3 dc = (d - c*ta)/la;
            B[i1] += k0*dc;
            B[i0] -= k0*dc;
            L.Mb[e] = k0*d.cross(ta);
        }
    }
    const double mi = p.mass*l_;
    double an, at;
    for(int i=1; i<N; ++i)
    {
        Vec3 q = R(i+1) - R(i-1);
        const double qn = q.norm();
        q = qn > 0.0 ? Vec3(q/qn) : Vec3::UnitZ();
        Vec3 F = T_[i] - T_[i-1] + external(i, R(i), V(i), q, l_, 0.0, an, at) - extra*mi*V(i);
        if(bend)
            F += L.B[i];
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
    L.F[0] = T_[0] + external(0, R(0), V(0), qa, 0.5*l_, 1.0, an, at);
    L.M[0] = me*Mat3::Identity() + an*(Mat3::Identity() - qa*qa.transpose()) + at*qa*qa.transpose();
    L.F[1] = -T_[N-1] + external(N, R(N), V(N), qb, 0.5*l_, -1.0, an, at);
    L.M[1] = me*Mat3::Identity() + an*(Mat3::Identity() - qb*qb.transpose()) + at*qb*qb.transpose();
    if(bend)
    {
        L.F[0] += L.B[0];
        L.F[1] += L.B[N];
    }
    a[o].setZero();
    a[o+N].setZero();
}

void mooring_system::accel(const state &s, double t, std::vector<Vec3> &a, double extra, bool upd) const
{
    a.resize(r_.size());
    for(const line_t &L : lines_)
        line_accel(L, s, t, a, extra, upd);
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
            const Vec3 vr = current_at(x(2)) + u - vx;
            F += fs*0.5*rho*p.CdA*vr.norm()*vr;
            if(kmode_==1 || kmode_==2)
                F += fs*rho*p.volume*(1.0 + p.Ca)*aw;
        }
        // seabed: contact of the sphere of the contact radius (along the normal of a surface or the fluid's bed);
        // a point without a radius or volume (a junction, a bridle node) has none of its own, its line ends
        // touch the seabed like the line's nodes
        const double rc = p.radius >= 0.0 ? p.radius : rv;
        double zb = P.zbed;
        Vec3 nb = Vec3::UnitZ();
        if(rc > 0.0 && (P.surf || P.fbed))
            seabed_at(i, x, P.surf, P.fbed, zb, nb);
        bool contact = false;
        if(rc > 0.0 && x(2) - rc/nb(2) < zb)
        {
            const double Ac = M_PI*rc*rc;
            const double Fz = std::max(0.0, (p.kb*((zb - x(2))*nb(2) + rc) - p.cb*vx.dot(nb))*Ac);
            F += Fz*nb;
            const Vec3 vh = vx - vx.dot(nb)*nb;
            const double vn = vh.norm();
            if(p.v_eps > 0.0)
            {
                if(Fz > 0.0 && vn > 0.0 && p.mu > 0.0)
                    F -= (p.mu*Fz*std::tanh(vn/p.v_eps)/vn)*vh;
            }
            else if(Fz > 0.0 && p.mu > 0.0 && sc_)
            {
                F += stick_friction(i, x, vx, nb, Vec3::Zero(), Fz, Ac, p.ks >= 0.0 ? p.ks : p.kb, p.cs >= 0.0 ? p.cs : p.cb,
                                    p.mu, p.mu, p.mu_static, upd);
                contact = true;
            }
        }
        if(upd && !contact && sc_)
            lose_anchor(i);
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

void mooring_system::integrate(state &s, const std::vector<end_motion> &em, double t0, double t1, double extra, bool at_rest,
                               stick_t *fs) const
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
    const int nn = (int)s.r.size();
    // stick-slip friction: anchors for all nodes (new ones without), none while settling (in the water at rest)
    sc_ = nullptr;
    if(fs && stick_ && !at_rest)
    {
        fs->a.resize(nn, Vec3::Constant(std::numeric_limits<double>::quiet_NaN()));
        fs->slip.resize(nn, 0);
        sc_ = fs;
    }
    state &s2 = ws_;
    s2 = s;
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
        // clamp directions: between the converged and the trial ones
        for(size_t k=0; k<att_.size(); ++k)
            if(clamped_)
            {
                const line_t &L = lines_[att_[k].line];
                if(!L.clamp[att_[k].end])
                    continue;
                const end_motion &e = em[k];
                const double y = T > 0.0 ? (tt - t0)/T : 1.0;
                const Vec3 d = e.d0 + y*(e.d1 - e.d0);
                L.cd[att_[k].end] = d.norm() > 0.0 ? Vec3(d.normalized()) : e.d1;
            }
        point_ends(x);
    };
    for(int k=0; k<n; ++k)
    {
        const double tk = t0 + k*h;
        // RK4 on (r, v) of the inner nodes and the points; the stage velocities are the stage states' v
        ends(s, tk);
        accel(s, tk, k1, extra, true);
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
        if(clamped_ && lines_[att_[k].line].clamp[att_[k].end])
            lines_[att_[k].line].cd[att_[k].end] = em[k].d1;
    }
    point_ends(s);
    accel(s, t1, k1, extra, true);
    sc_ = nullptr;
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
    if(clamped_)
        for(size_t k=0; k<na; ++k)
            hit = hit && (att_[k].ds.isZero() || m.body(att_[k].body).R*att_[k].ds==cd_[k]);
    if(!hit)
    {
        std::vector<end_motion> em(na);
        cd_.resize(na);
        for(size_t k=0; k<na; ++k)
        {
            em[k] = end_motion{att_[k].r, att_[k].v, r1[k], v1[k]};
            if(clamped_ && !att_[k].ds.isZero())
            {
                em[k].d0 = att_[k].d;
                em[k].d1 = cd_[k] = m.body(att_[k].body).R*att_[k].ds;
            }
        }
        state s{r_, v_};
        stick_t f = st_;
        integrate(s, em, t_, t, 0.0, false, &f);
        cs_ = std::move(s);
        cst_ = std::move(f);
        cw_.resize(lines_.size());
        for(size_t k=0; k<lines_.size(); ++k)
        {
            const line_t &L = lines_[k];
            line_work &w = cw_[k];
            w.T = L.T;
            for(int e=0; e<2; ++e) {w.F[e] = L.F[e]; w.M[e] = L.M[e]; w.cd[e] = L.cd[e];}
        }
        cr_ = r1;
        cv_ = v1;
        cF_.resize(na);
        cM_.resize(na);
        cMb_.resize(na);
        for(size_t k=0; k<na; ++k)
        {
            const line_t &L = lines_[att_[k].line];
            cF_[k] = L.F[att_[k].end];
            cM_[k] = L.M[att_[k].end];
            cMb_[k] = L.Mb[att_[k].end];
        }
        ct_ = t;
        cvalid_ = true;
    }
    // the body carries the end node: minus its mass and added mass times the attachment's acceleration
    for(size_t k=0; k<na; ++k)
    {
        const attach_t &a = att_[k];
        m.add_force(f, a.body, r1[k], cF_[k] - cM_[k]*m.body(a.body).point_acceleration(a.s));
        if(clamped_ && !a.ds.isZero())
            m.add_moment(f, a.body, cMb_[k]);
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
        if(clamped_ && !att_[k].ds.isZero())
        {
            em[k].d0 = att_[k].d;
            em[k].d1 = B.R*att_[k].ds;
        }
    }
    // the converged state is the last trial that eval() integrated, up to the last Newton increment (the time
    // integrator stops on an increment below its tolerance, 1e-12 relative, without evaluating the loads there): take
    // that integration instead of integrating the step again (positions and velocities within 1e-10 relative, the
    // end forces of the lines then within 1e-10 EA/l of a fairlead)
    bool hit = cvalid_ && t==ct_ && cw_.size()==lines_.size();
    auto near = [](const Vec3 &a, const Vec3 &b) {return (a - b).cwiseAbs().maxCoeff() <= 1.0e-10*(1.0 + a.cwiseAbs().maxCoeff());};
    for(size_t k=0; k<na && hit; ++k)
        hit = near(em[k].r1, cr_[k]) && near(em[k].v1, cv_[k]) && (!clamped_ || att_[k].ds.isZero() || near(em[k].d1, cd_[k]));
    state s;
    if(hit)
    {
        s = std::move(cs_);
        st_ = std::move(cst_);
        for(size_t k=0; k<lines_.size(); ++k)
        {
            const line_t &L = lines_[k];
            const line_work &w = cw_[k];
            L.T = w.T;
            for(int e=0; e<2; ++e) {L.F[e] = w.F[e]; L.M[e] = w.M[e]; L.cd[e] = w.cd[e];}
        }
        ++commit_hits_;
    }
    else
    {
        s = state{r_, v_};
        integrate(s, em, t_, t, 0.0, false, &st_);
    }
    ++commits_;
    for(size_t k=0; k<na; ++k)
    {
        attach_t &a = att_[k];
        const line_t &L = lines_[a.line];
        a.F = L.F[a.end] - L.M[a.end]*m.body(a.body).point_acceleration(a.s);
        a.r = em[k].r1;
        a.v = em[k].v1;
        a.d = em[k].d1;
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
