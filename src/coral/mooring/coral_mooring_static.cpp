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

// REEF3D::CORAL: static equilibrium of a mooring system (mooring_system::equilibrium), in its own file: next
// to the time integration in coral_mooring.cpp it changed the compiler's inlining of the sub-step loop
// (5 times slower).

#include"coral/mooring/coral_mooring.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_linsolve.h"
#include<algorithm>
#include<cmath>
#include<limits>
#include<set>

namespace coral
{

// round-off of the accelerations of line L: 1000 machine epsilons of the coordinates (their largest
// magnitude) times the axial stiffness per node mass, EA / (l m l)
double mooring_system::roundoff(const line_t &L) const
{
    double X = 1.0;
    for(int i=0; i<=L.N; ++i)
        X = std::max(X, r_[L.off+i].cwiseAbs().maxCoeff());
    const double m = std::max(L.p.mass*L.l, 1.0e-300);
    return 1.0e3*std::numeric_limits<double>::epsilon()*X*L.p.EA_max()/(L.l*m);
}

bool mooring_system::equilibrium(const model &m, double t)
{
    for(attach_t &a : att_)
    {
        a.r = m.body(a.body).point(a.s);
        a.v.setZero();
        r_[a.idx] = a.r;
        if(!a.ds.isZero())
            lines_[a.line].cd[a.end] = a.d = m.body(a.body).R*a.ds;
    }
    for(auto &x : v_) x.setZero();
    state s0{r_, v_};
    point_ends(s0);
    r_ = s0.r;
    ++gen_;
    // the water at rest at t (the surface of the nodes' start positions), the current acts
    if(waves_)
    {
        kinematics(0, t);
        kmode_ = 3;
    }
    else
        kmode_ = 0;
    // lines on their own (block-tridiagonal Newton); if that fails (stiff lines whose catenary start has
    // slack segments) or with points, the whole system (compression stage, dynamic relaxation)
    bool flat = true;
    for(const line_t &L : lines_)
        flat = flat && !L.surf && !L.fbed && !(L.p.EI > 0.0);
    bool ok = points_.empty() && flat && equilibrium_lines(m, t);
    if(!ok)
        ok = equilibrium_system(m, t);
    for(auto &x : v_) x.setZero();
    // frictionless equilibrium (friction is indeterminate at rest): the nodes stick where they lie from the start
    st_ = stick_t();
    ++gen_;
    cvalid_ = false;
    return ok;
}

// each line on its own between its fixed ends: block-tridiagonal Jacobian
bool mooring_system::equilibrium_lines(const model &m, double t)
{
    bool all = true;
    for(const line_t &Ln : lines_)
    {
        const int N = Ln.N, o = Ln.off;
        const mooring_props &pr_ = Ln.p;
        const double l_ = Ln.l, zbed_ = Ln.zbed;
        state s{r_, v_}, st = s;
        std::vector<Vec3> a, ap;
        auto resid = [&](const state &x, std::vector<Vec3> &out)
        {
            out.resize(r_.size());
            line_accel(Ln, x, t, out, 0.0);
            double n = 0.0;
            for(int i=1; i<N; ++i)
            {
                if(!out[o+i].allFinite() || !x.r[o+i].allFinite())
                    return 1.0e300;
                n += out[o+i].squaredNorm();
            }
            return n;
        };
        auto rmax = [&](const std::vector<Vec3> &out)
        {
            double n = 0.0;
            for(int i=1; i<N; ++i)
                n = std::max(n, out[o+i].norm());
            return n;
        };
        // nodes on the seabed start at the depth the penalty lets them sink (contact active: Newton sees its
        // stiffness instead of the kink at the seabed)
        if(pr_.kb > 0.0 && zbed_ > -1.0e29)
        {
            const double wsub = (pr_.mass - pr_.rho*0.25*M_PI*pr_.d*pr_.d)*m.gravity.norm();
            const double dz = std::max(wsub, 0.0)/(pr_.kb*pr_.d);
            for(int i=1; i<N; ++i)
                if(s.r[o+i](2) <= zbed_ + 1.0e-9*std::max(1.0, std::fabs(zbed_)))
                    s.r[o+i](2) = zbed_ - dz;
        }
        const double tol = 1.0e-9*std::max(1.0, m.gravity.norm());
        const double eps = 1.0e-7*std::max(l_, 1.0e-3);
        std::vector<Mat3> L(N + 1, Mat3::Zero()), D(N + 1, Mat3::Zero()), U(N + 1, Mat3::Zero()), Dp(N + 1, Mat3::Zero());
        std::vector<Vec3> b(N + 1), dx(N + 1);
        double r = resid(s, a), lam = -1.0;
        bool ok = rmax(a) <= tol;
        for(int it=0; it<200 && !ok; ++it)
        {
            // Jacobian d a_i / d x_j (j = i-1, i, i+1): nodes of one colour (i mod 3) perturbed together
            for(int c=0; c<3; ++c)
                for(int d=0; d<3; ++d)
                {
                    st = s;
                    for(int j=1; j<N; ++j)
                        if(j%3==c) st.r[o+j](d) += eps;
                    resid(st, ap);
                    for(int i=1; i<N; ++i)
                    {
                        const Vec3 col = (ap[o+i] - a[o+i])/eps;
                        if((i-1)%3==c && i-1 >= 1) L[i].col(d) = col;
                        if(i%3==c) D[i].col(d) = col;
                        if((i+1)%3==c && i+1 <= N-1) U[i].col(d) = col;
                    }
                }
            // scale of the Jacobian: its largest diagonal entry, at least the axial stiffness per node mass
            double dmax = pr_.mass > 0.0 ? pr_.EA_max()/(l_*pr_.mass*l_) : 1.0;
            for(int i=1; i<N; ++i)
                dmax = std::max(dmax, D[i].cwiseAbs().maxCoeff());
            if(lam < 0.0)
                lam = 1.0e-10*dmax;
            // damped Newton (Levenberg-Marquardt type): (J - lam I) dx = -a, block-tridiagonal (Thomas); J = da/dx is
            // negative (semi-)definite, a large lam moves the nodes along their forces; the step must lower |a|
            bool moved = false;
            for(int k=0; k<30 && !moved; ++k)
            {
                for(int i=1; i<N; ++i)
                {
                    Mat3 Di = D[i] - lam*Mat3::Identity();
                    Vec3 bi = -a[o+i];
                    if(i > 1)
                    {
                        const Mat3 W = L[i]*Dp[i-1].inverse();
                        Di -= W*U[i-1];
                        bi -= W*b[i-1];
                    }
                    Dp[i] = Di;
                    b[i] = bi;
                }
                for(int i=N-1; i>=1; --i)
                    dx[i] = Dp[i].inverse()*(b[i] - (i < N-1 ? Vec3(U[i]*dx[i+1]) : Vec3::Zero()));
                st = s;
                for(int i=1; i<N; ++i)
                    st.r[o+i] += dx[i];
                const double rn = resid(st, ap);
                if(rn < r)
                {
                    s = st;
                    a = ap;
                    r = rn;
                    lam = std::max(0.1*lam, 1.0e-14*dmax);
                    moved = true;
                }
                else
                    lam = std::max(10.0*lam, 1.0e-10*dmax);
            }
            if(!moved)
                break;
            ok = rmax(a) <= tol;
        }
        // stiff lines: the accelerations cannot fall below their round-off (the axial stiffness per node mass
        // times the rounding of the coordinates); there the stagnated Newton iteration has converged
        ok = ok || rmax(a) <= roundoff(Ln);
        // not converged (e.g. a straight line without tension, no transverse stiffness): keep the start shape
        if(ok)
            for(int i=1; i<N; ++i)
                r_[o+i] = s.r[o+i];
        all = all && ok;
    }
    return all;
}

// the points where the catenaries of their lines (and their own weight, buoyancy and seabed contact) are in
// balance (Newton with backtracking on the point positions, finite-difference Jacobian), the nodes of all
// lines then on the catenaries between their ends: the start of the Newton iteration of the whole system
bool mooring_system::point_statics(const model &m)
{
    const int np = (int)points_.size();
    const double gm = m.gravity.norm();
    std::vector<Vec3> P(np);
    for(int j=0; j<np; ++j)
        P[j] = r_[points_[j].idx];
    auto end_pos = [&](const line_t &L, int e, const std::vector<Vec3> &Q) -> Vec3
    {
        return L.pt[e] >= 0 ? Q[L.pt[e]] : att_[L.att[e]].r;
    };
    auto resid = [&](const std::vector<Vec3> &Q, VecX &R)
    {
        R.setZero(3*np);
        for(const line_t &L : lines_)
        {
            if(L.pt[0] < 0 && L.pt[1] < 0)
                continue;
            Vec3 Fa, Fb;
            line_catenary(L, end_pos(L, 0, Q), end_pos(L, 1, Q), gm, Fa, Fb, nullptr);
            if(L.pt[0] >= 0) R.segment<3>(3*L.pt[0]) += Fa;
            if(L.pt[1] >= 0) R.segment<3>(3*L.pt[1]) += Fb;
        }
        for(int j=0; j<np; ++j)
        {
            const point_t &Pt = points_[j];
            const mooring_point_props &p = Pt.p;
            const double rho = lines_[Pt.ends.empty() ? 0 : Pt.ends[0].first].p.rho;
            double eta;
            Vec3 u, du;
            water(Pt.idx, Q[j], 0.0, eta, u, du);
            const double rv = std::cbrt(0.75*p.volume/M_PI), hz = std::max(2.0*rv, 0.1);
            const double fs = std::clamp((eta - Q[j](2))/hz + 0.5, 0.0, 1.0);
            Vec3 F = (p.mass - fs*rho*p.volume)*m.gravity;
            const double rc = p.radius >= 0.0 ? p.radius : rv;
            double zb = Pt.zbed;
            Vec3 nb = Vec3::UnitZ();
            if(rc > 0.0 && (Pt.surf || Pt.fbed))
                seabed_at(Pt.idx, Q[j], Pt.surf, Pt.fbed, zb, nb);
            if(rc > 0.0 && Q[j](2) - rc/nb(2) < zb)
                F += p.kb*((zb - Q[j](2))*nb(2) + rc)*M_PI*rc*rc*nb;
            // the seabed contact of the end nodes of its lines (their half-segments), as in the dynamics: it
            // holds a point without a sphere of its own (a junction on the seabed) at the line's depth
            for(const auto &e : Pt.ends)
            {
                const line_t &L = lines_[e.first];
                if(!(L.p.kb > 0.0) || !(L.zbed > -1.0e29))
                    continue;
                double zl = L.zbed;
                Vec3 nl = Vec3::UnitZ();
                if(L.surf || L.fbed)
                    seabed_at(L.off + (e.second==0 ? 0 : L.N), Q[j], L.surf, L.fbed, zl, nl);
                const double pen = (zl - Q[j](2))*nl(2);
                if(pen > 0.0)
                    F += L.p.kb*pen*L.p.d*0.5*L.l*nl;
            }
            R.segment<3>(3*j) += F;
        }
    };
    VecX R, Rp;
    resid(P, R);
    double scale = 1.0;
    for(const line_t &L : lines_)
        scale = std::max(scale, std::fabs((L.p.mass - L.p.rho*0.25*M_PI*L.p.d*L.p.d)*gm*L.p.L));
    for(const point_t &Pt : points_)
        scale = std::max(scale, std::fabs(Pt.p.mass*gm));
    bool ok = R.norm() < 1.0e-8*scale;
    for(int it=0; it<100 && !ok; ++it)
    {
        MatX J(3*np, 3*np);
        const double eps = 1.0e-6;
        for(int c=0; c<3*np; ++c)
        {
            std::vector<Vec3> Q = P;
            Q[c/3](c%3) += eps;
            resid(Q, Rp);
            J.col(c) = (Rp - R)/eps;
        }
        VecX dx = J.fullPivLu().solve(-R);
        if(!dx.allFinite())
            break;
        // backtracking: at most 10 % of the shortest line per step, the residual must fall
        double lmin = 1.0e300;
        for(const line_t &L : lines_) lmin = std::min(lmin, L.p.L);
        double f = std::min(1.0, 0.1*lmin/std::max(dx.lpNorm<Eigen::Infinity>(), 1.0e-300));
        bool moved = false;
        for(int k=0; k<30; ++k, f *= 0.5)
        {
            std::vector<Vec3> Q = P;
            for(int j=0; j<np; ++j)
                Q[j] += f*dx.segment<3>(3*j);
            resid(Q, Rp);
            if(Rp.allFinite() && Rp.norm() < R.norm())
            {
                P = Q;
                R = Rp;
                moved = true;
                break;
            }
        }
        if(!moved)
            break;
        ok = R.norm() < 1.0e-8*scale;
    }
    // the points there, the nodes of all lines on the catenaries between their ends
    std::vector<Vec3> rl;
    Vec3 Fa, Fb;
    for(int j=0; j<np; ++j)
        r_[points_[j].idx] = P[j];
    for(const line_t &L : lines_)
    {
        line_catenary(L, end_pos(L, 0, P), end_pos(L, 1, P), gm, Fa, Fb, &rl);
        for(int i=0; i<=L.N; ++i)
            r_[L.off+i] = rl[i];
    }
    return ok;
}

// lines and points together: sparse Jacobian of the free nodes by coloured finite differences
bool mooring_system::equilibrium_system(const model &m, double t)
{
    // start: the points on the catenaries of their lines (not with bending: the start shapes stay)
    bool bend = false;
    for(const line_t &L : lines_)
        bend = bend || L.p.EI > 0.0;
    if(!bend)
        point_statics(m);
    const int nf = (int)free_.size();
    std::vector<int> pos(r_.size(), -1);           // node -> unknown
    for(int k=0; k<nf; ++k)
        pos[free_[k]] = k;
    // rows(j): the unknowns whose acceleration depends on unknown j
    std::vector<std::set<int>> rows(nf);
    auto link = [&](int i, int j) {if(pos[i] >= 0 && pos[j] >= 0) {rows[pos[i]].insert(pos[j]); rows[pos[j]].insert(pos[i]);}};
    std::vector<int> node_of(r_.size());           // the node itself, or the point a line end sits on
    for(size_t i=0; i<r_.size(); ++i)
        node_of[i] = (int)i;
    for(const line_t &L : lines_)
        for(int e=0; e<2; ++e)
            if(L.pt[e] >= 0)
                node_of[L.off + (e==0 ? 0 : L.N)] = points_[L.pt[e]].idx;
    for(int k=0; k<nf; ++k)
        rows[k].insert(k);
    for(const line_t &L : lines_)
        for(int i=0; i<L.N; ++i)
        {
            // neighbours along the line (the ends on points stand for the points); the end forces depend on
            // the two nodes of the end segment and act on the point
            const int a = node_of[L.off + i], b = node_of[L.off + i + 1];
            link(a, b);
            // bending: the turning angle at a node joins its two neighbours
            if(L.p.EI > 0.0 && i+2 <= L.N)
                link(a, node_of[L.off + i + 2]);
        }
    // seabed contact: start at the penalty depth (see equilibrium_lines)
    for(const line_t &L : lines_)
    {
        const mooring_props &p = L.p;
        if(p.kb > 0.0 && L.zbed > -1.0e29)
        {
            const double wsub = (p.mass - p.rho*0.25*M_PI*p.d*p.d)*m.gravity.norm();
            const double dz = std::max(wsub, 0.0)/(p.kb*p.d);
            for(int i=1; i<L.N; ++i)
            {
                double zb = L.zbed;
                Vec3 nb;
                if(L.surf || L.fbed) seabed_at(L.off+i, r_[L.off+i], L.surf, L.fbed, zb, nb);
                if(r_[L.off+i](2) <= zb + 1.0e-9*std::max(1.0, std::fabs(zb)))
                    r_[L.off+i](2) = zb - dz;
            }
        }
    }
    // colours: unknowns of one colour have disjoint rows
    std::vector<int> color(nf, -1);
    std::vector<std::set<int>> used;
    int ncol = 0;
    for(int j=0; j<nf; ++j)
    {
        int c = 0;
        for(; c<ncol; ++c)
        {
            bool clash = false;
            for(int i : rows[j])
                if(used[c].count(i)) {clash = true; break;}
            if(!clash) break;
        }
        if(c==ncol) {used.emplace_back(); ++ncol;}
        color[j] = c;
        used[c].insert(rows[j].begin(), rows[j].end());
    }
    state s{r_, v_}, st = s;
    std::vector<Vec3> a, ap;
    auto resid = [&](state &x, std::vector<Vec3> &out)
    {
        point_ends(x);
        accel(x, t, out, 0.0);
        double n = 0.0;
        for(int i : free_)
        {
            if(!out[i].allFinite() || !x.r[i].allFinite())
                return 1.0e300;
            n += out[i].squaredNorm();
        }
        return n;
    };
    auto rmax = [&](const std::vector<Vec3> &out)
    {
        double n = 0.0;
        for(int i : free_)
            n = std::max(n, out[i].norm());
        return n;
    };
    double lmin = 1.0e300;
    for(const line_t &L : lines_)
        lmin = std::min(lmin, L.l);
    const double tol = 1.0e-9*std::max(1.0, m.gravity.norm());
    const double eps = 1.0e-7*std::max(lmin, 1.0e-3);
    std::vector<Eigen::Triplet<double>> J, Jl;
    linear_solver ls;
    auto newton = [&](bool watchdog) -> bool
    {
    bool ok = false;
    // two stages: slack segments resist compression with 1 % of EA (the start shapes have chords shorter
    // than the segments: without it the points float on slack segments, the Jacobian is singular), then
    // the line law itself from there
    for(int stage=0; stage<2; ++stage)
    {
    creg_ = stage==0 ? 0.01 : 0.0;
    double r = resid(s, a), lam = -1.0;
    ok = rmax(a) <= tol;
    // watchdog (if asked): up to 5 full Newton steps in a row may raise the residual (stiff axial springs that
    // the linearized step stretches, e.g. a bending beam starting straight); the best state is kept
    state sbest = s;
    std::vector<Vec3> abest = a;
    double rbest = r;
    int nonmono = 0;
    for(int it=0; it<300 && !ok; ++it)
    {
        J.clear();
        double dmax = 1.0;
        for(const line_t &L : lines_)
            if(L.p.mass > 0.0) dmax = std::max(dmax, L.p.EA_max()/(L.l*L.p.mass*L.l));
        for(int c=0; c<ncol; ++c)
            for(int d=0; d<3; ++d)
            {
                st = s;
                for(int j=0; j<nf; ++j)
                    if(color[j]==c) st.r[free_[j]](d) += eps;
                resid(st, ap);
                for(int j=0; j<nf; ++j)
                {
                    if(color[j] != c) continue;
                    for(int i : rows[j])
                    {
                        const Vec3 col = (ap[free_[i]] - a[free_[i]])/eps;
                        for(int q=0; q<3; ++q)
                            if(col(q) != 0.0)
                            {
                                J.emplace_back(3*i+q, 3*j+d, col(q));
                                if(i==j && q==d) dmax = std::max(dmax, std::fabs(col(q)));
                            }
                    }
                }
            }
        if(lam < 0.0)
            lam = 1.0e-10*dmax;
        bool moved = false;
        for(int k=0; k<30 && !moved; ++k)
        {
            Jl = J;
            for(int q=0; q<3*nf; ++q)
                Jl.emplace_back(q, q, -lam);
            VecX rhs(3*nf);
            for(int i=0; i<nf; ++i)
                rhs.segment<3>(3*i) = -a[free_[i]];
            VecX dx;
            if(ls.compute(3*nf, Jl))
                dx = ls.solve(rhs);
            if(dx.size()==3*nf && dx.allFinite())
            {
                st = s;
                for(int i=0; i<nf; ++i)
                    st.r[free_[i]] += dx.segment<3>(3*i);
                const double rn = resid(st, ap);
                const bool watch = watchdog && k==0 && rn >= r && rn < 1.0e12*rbest && nonmono < 5;
                if(rn < r || watch)
                {
                    s = st;
                    a = ap;
                    r = rn;
                    lam = std::max(0.1*lam, 1.0e-14*dmax);
                    moved = true;
                    nonmono = watch ? nonmono + 1 : 0;
                    if(r < rbest)
                    {
                        sbest = s;
                        abest = a;
                        rbest = r;
                        nonmono = 0;
                    }
                    continue;
                }
            }
            lam = std::max(10.0*lam, 1.0e-10*dmax);
        }
        if(!moved)
            break;
        ok = rmax(a) <= tol;
    }
    if(r > rbest)
    {
        s = sbest;
        a = abest;
        r = rbest;
    }
    }
    creg_ = 0.0;
    double tro = 0.0;
    for(const line_t &L : lines_)
        tro = std::max(tro, roundoff(L));
    ok = ok || rmax(a) <= tro;
    return ok;
    };
    // dynamic relaxation with kinetic damping (the velocities are set to zero whenever the kinetic energy
    // of the free nodes falls), the body ends held, in the water at rest, for at most T_relax; then Newton
    // from there (a short one first: the catenary start shapes have slightly slack segments)
    std::vector<end_motion> em;
    for(const attach_t &at : att_)
    {
        em.push_back(end_motion{at.r, Vec3::Zero(), at.r, Vec3::Zero()});
        em.back().d0 = em.back().d1 = at.d;
    }
    std::vector<double> mw(r_.size(), 0.0);
    for(const line_t &L : lines_)
        for(int i=1; i<L.N; ++i)
            mw[L.off+i] = L.p.mass*L.l;
    for(const point_t &P : points_)
    {
        mw[P.idx] = P.p.mass;
        for(const auto &e : P.ends)
            mw[P.idx] += 0.5*lines_[e.first].p.mass*lines_[e.first].l;
    }
    auto relax = [&](double T_relax)
    {
        double ke0 = 0.0, tt = t;
        const double chunk = std::max(20.0*dt_, 0.01);
        const int n = (int)std::ceil(T_relax/chunk);
        for(int c=0; c<n; ++c)
        {
            integrate(s, em, tt, tt + chunk, 0.0, true);
            tt += chunk;
            double ke = 0.0;
            for(int i : free_)
                ke += mw[i]*s.v[i].squaredNorm();
            if(ke < ke0)
                for(int i : free_)
                    s.v[i].setZero();
            ke0 = ke;
            if(c%50==49)
            {
                state q{s.r, std::vector<Vec3>(s.v.size(), Vec3::Zero())};
                resid(q, a);
                if(rmax(a) < 1.0e-3*std::max(1.0, m.gravity.norm()))
                    break;
            }
        }
        for(auto &x : s.v) x.setZero();
    };
    // Newton without the watchdog first (it can throw a nearly converged system far off, e.g. a bridle node
    // on two taut legs), then with it from the best state; with bending the other way round
    auto newton2 = [&]() {return newton(bend) || newton(!bend);};
    bool ok = newton2();
    for(double T_relax : {2.0, 30.0})
    {
        if(ok)
            break;
        relax(T_relax);
        ok = newton2();
    }
    // not converged (neutral directions, e.g. a line lying slack on a seabed without friction): the relaxed
    // shape, closer to the equilibrium than the start
    point_ends(s);
    r_ = s.r;
    return ok;
}

}
