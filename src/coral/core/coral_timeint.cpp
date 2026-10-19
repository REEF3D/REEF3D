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

#include"coral/core/coral_timeint.h"
#include"coral/core/coral_state.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include"coral/core/coral_fdjac.h"
#include"coral/core/coral_linsolve.h"
#include<algorithm>
#include<cmath>

namespace coral
{

timeint::timeint(const timeint_settings &s) : set_(s)
{
    const double r = std::max(0.0, std::min(1.0, s.rho_inf));
    am_ = (2.0*r - 1.0)/(r + 1.0);
    af_ = r/(r + 1.0);
    gam_ = 0.5 + af_ - am_;
    beta_ = 0.25*(gam_ + 0.5)*(gam_ + 0.5);
    lu_.dense_max = s.dense_max;
}

void timeint::store(const model &m)
{
    const int nb = m.nbody(), nq = m.ncoord();
    xn_.resize(nb);
    Rn_.resize(nb);
    nun_.resize(6*nb + nq);
    for(int i=0; i<nb; ++i)
    {
        const rigid_body &b = m.bodies[i];
        xn_[i] = b.x;
        Rn_[i] = b.R;
        nun_.segment<3>(6*i) = b.v;
        nun_.segment<3>(6*i+3) = b.w;
    }
    cqn_ = m.cq;
    if(nq)
        nun_.tail(nq) = m.cv;
}

void timeint::restore(model &m) const
{
    for(int i=0; i<m.nbody(); ++i)
    {
        rigid_body &b = m.bodies[i];
        b.x = xn_[i];
        b.R = Rn_[i];
        b.v = nun_.segment<3>(6*i);
        b.w = nun_.segment<3>(6*i+3);
        b.vd = dnun_.segment<3>(6*i);
        b.wd = dnun_.segment<3>(6*i+3);
    }
    const int nq = m.ncoord();
    if(nq)
    {
        m.cq = cqn_;
        m.cv = nun_.tail(nq);
        m.ca = dnun_.tail(nq);
    }
}

// state at t_{n+1} from the trial dnu_{n+1}
void timeint::trial(model &m, const VecX &dnu, double h) const
{
    const VecX a1 = ((1.0 - af_)*dnu + af_*dnun_ - am_*an_)/(1.0 - am_);
    const VecX nu1 = nun_ + h*((1.0 - gam_)*an_ + gam_*a1);
    const VecX dq = h*nun_ + h*h*((0.5 - beta_)*an_ + beta_*a1);
    for(int i=0; i<m.nbody(); ++i)
    {
        rigid_body &b = m.bodies[i];
        b.x = xn_[i] + dq.segment<3>(6*i);
        b.R = Rn_[i]*exp_so3(dq.segment<3>(6*i+3));
        b.v = nu1.segment<3>(6*i);
        b.w = nu1.segment<3>(6*i+3);
        b.vd = dnu.segment<3>(6*i);
        b.wd = dnu.segment<3>(6*i+3);
    }
    const int nq = m.ncoord();
    if(nq)
    {
        m.cq = cqn_ + dq.tail(nq);
        m.cv = nu1.tail(nq);
        m.ca = dnu.tail(nq);
    }
}

void timeint::residual(model &m, double t1, double h, const VecX &u, VecX &r, bool fd_only)
{
    const int n = m.ndof(), nc = (int)u.size() - n;
    trial(m, u.head(n), h);
    m.forces(t1, f_, fd_only);
    m.inertia_forces(in_);
    m.constraints(g_, G_);
    r.resize(n + nc);
    r.head(n) = in_ - f_;
    if(nc)
    {
        r.head(n) += scale_*G_.transpose()*u.tail(nc);
        r.tail(nc) = scale_*g_/(beta_*h*h);
    }
}

// t0: positions and velocities given, G dnu = -c (c = dG/dt nu)
void timeint::residual_initial(model &m, double t0, const VecX &c, const VecX &u, VecX &r, bool fd_only)
{
    const int n = m.ndof(), nc = (int)u.size() - n;
    for(int i=0; i<m.nbody(); ++i)
    {
        m.bodies[i].vd = u.segment<3>(6*i);
        m.bodies[i].wd = u.segment<3>(6*i+3);
    }
    if(m.ncoord())
        m.ca = u.segment(6*m.nbody(), m.ncoord());
    m.forces(t0, f_, fd_only);
    m.inertia_forces(in_);
    m.constraints(g_, G_);
    r.resize(n + nc);
    r.head(n) = in_ - f_;
    if(nc)
    {
        r.head(n) += scale_*G_.transpose()*u.tail(nc);
        r.tail(nc) = scale_*(G_*u.head(n) + c);
    }
}

void timeint::build_jacobian(model &m, double t1, double h, const VecX &u, const VecX &r, bool initial, const VecX &c)
{
    const int n = m.ndof(), N = (int)u.size(), nc = N - n, nb = m.nbody();
    const bool tang = m.has_tangents();
    T_.clear();
    VecX up = u, rp, r0;
    auto res = [&](const VecX &uu, VecX &rr)
    {
        if(initial) residual_initial(m, t1, c, uu, rr, tang);
        else residual(m, t1, h, uu, rr, tang);
    };
    // finite differences of the loads without an analytic tangent (and of inertia and constraints)
    if(tang)
        res(u, r0);
    fd_jacobian(m, T_, tang ? r0 : r,
        [&](int col, double e) {up(col) = u(col) + e;},
        [&](VecX &rr) {res(up, rr);},
        [&](int col) {return set_.fd_eps*std::max(1.0, std::fabs(u(col)));},
        false, tang);
    // restore the work arrays of the unperturbed state (all loads)
    VecX rr;
    if(initial)
        residual_initial(m, t1, c, u, rr);
    else
        residual(m, t1, h, u, rr);
    // analytic tangents: d/d dnu of -F = -(cq dF/dq + cv dF/dnu + ca dF/dnu'), with dq/d dnu = beta h^2 c,
    // dnu/d dnu = gamma h c (c = (1 - alpha_f)/(1 - alpha_m)) and the loads see dnu itself as the
    // acceleration (ca = 1); rotation columns through the tangent operator of the trial increment
    // (R = R_n exp(dq): R exp(T(dq) d(dq)))
    if(tang)
    {
        const double cf = initial ? 0.0 : (1.0 - af_)/(1.0 - am_);
        const double cK = initial ? 0.0 : -beta_*h*h*cf, cC = initial ? 0.0 : -gam_*h*cf, cM = -1.0;
        // velocity and acceleration parts directly
        m.tangents(t1, 0.0, cC, cM, T_);
        // configuration part: rotation columns through the tangent operator of the trial increment
        if(cK != 0.0)
        {
            const VecX a1 = ((1.0 - af_)*u.head(n) + af_*dnun_ - am_*an_)/(1.0 - am_);
            const VecX dq = h*nun_ + h*h*((0.5 - beta_)*an_ + beta_*a1);
            std::vector<Mat3> Tr(nb);
            for(int i=0; i<nb; ++i)
                Tr[i] = tangent_so3(dq.segment<3>(6*i+3));
            std::vector<triplet> K;
            m.tangents(t1, cK, 0.0, 0.0, K);
            for(const auto &t : K)
            {
                const int col = t.col();
                if(col < 6*nb && col%6 >= 3)
                {
                    const int i = col/6, k = col%6 - 3;
                    for(int j=0; j<3; ++j)
                        if(Tr[i](k, j) != 0.0)
                            T_.emplace_back(t.row(), 6*i+3+j, t.value()*Tr[i](k, j));
                }
                else
                    T_.push_back(t);
            }
        }
    }
    // lambda columns are exact
    if(nc)
        for(int j=0; j<nc; ++j)
            for(int i=0; i<n; ++i)
                if(G_(j, i) != 0.0)
                    T_.emplace_back(i, n + j, scale_*G_(j, i));
    lu_.compute(N, T_);
    jac_ok_ = true;
    jac_h_ = h;
    ++stats.jacobians;
}

void timeint::initialize(model &m, double t0)
{
    m.initialize(t0);
    const int n = m.ndof(), nc = m.ncon();

    double msum = 0.0;
    for(const auto &b : m.bodies)
        msum += b.mass;
    scale_ = std::max(1.0, m.nbody() ? msum/m.nbody() : 1.0);

    store(m);
    dnun_.setZero(n);
    an_.setZero(n);
    lam_.setZero(nc);

    // c = d/ds [G(q o exp(s nu)) nu] at s = 0, by central differences
    VecX c = VecX::Zero(nc);
    if(nc)
    {
        const double vmax = nun_.cwiseAbs().maxCoeff();
        if(vmax > 0.0)
        {
            const double ds = 1.0e-5/vmax;
            VecX gp, gm;
            MatX Gp, Gm;
            for(int sgn=-1; sgn<=1; sgn+=2)
            {
                for(int i=0; i<m.nbody(); ++i)
                {
                    rigid_body &b = m.bodies[i];
                    b.x = xn_[i] + sgn*ds*b.v;
                    b.R = Rn_[i]*exp_so3(sgn*ds*b.w);
                }
                if(sgn<0) m.constraints(gm, Gm);
                else      m.constraints(gp, Gp);
            }
            for(int i=0; i<m.nbody(); ++i)
            {
                m.bodies[i].x = xn_[i];
                m.bodies[i].R = Rn_[i];
            }
            c = (Gp - Gm)*nun_/(2.0*ds);
        }
    }

    VecX u = VecX::Zero(n + nc), r, du;
    for(int it=0; it<set_.max_iter; ++it)
    {
        residual_initial(m, t0, c, u, r);
        build_jacobian(m, t0, 0.0, u, r, true, c);
        du = -lu_.solve(r);
        u += du;
        if(du.head(n).cwiseAbs().maxCoeff() <= set_.tol*(1.0 + u.head(n).cwiseAbs().maxCoeff()))
            break;
    }
    residual_initial(m, t0, c, u, r);
    dnun_ = u.head(n);
    an_ = dnun_;
    lam_ = scale_*u.tail(nc);
    if(m.ncoord())
        m.ca = dnun_.tail(m.ncoord());
    jac_ok_ = false;
    if(nc)
        stats.max_constraint = std::max(stats.max_constraint, g_.cwiseAbs().maxCoeff());
}

// Newton has converged when the increment of dnu changes positions and velocities by less than
// tol (relative to 1 + the largest position and velocity): the accelerations themselves carry the
// round-off of the positions amplified by 1/(beta h^2) through the constraints
bool timeint::converged(const model &m, double ddnu, double h) const
{
    const double c = (1.0 - af_)/(1.0 - am_);
    double xmax = 0.0, vmax = 0.0;
    for(const auto &b : m.bodies)
    {
        xmax = std::max(xmax, b.x.cwiseAbs().maxCoeff());
        vmax = std::max(vmax, std::max(b.v.cwiseAbs().maxCoeff(), b.w.cwiseAbs().maxCoeff()));
    }
    return beta_*h*h*c*ddnu <= set_.tol*(1.0 + xmax) && gam_*h*c*ddnu <= set_.tol*(1.0 + vmax);
}

bool timeint::stagnated(const model &m, double ddnu, double h) const
{
    const double c = (1.0 - af_)/(1.0 - am_);
    double xmax = 0.0;
    for(const auto &b : m.bodies)
        xmax = std::max(xmax, b.x.cwiseAbs().maxCoeff());
    return beta_*h*h*c*ddnu <= 1.0e-9*(1.0 + xmax);
}

bool timeint::step(model &m, double t, double h)
{
    const int n = m.ndof(), nc = m.ncon();
    const double t1 = t + h;
    const VecX c;
    VecX u(n + nc), r, du;

    if(jac_h_!=h)
        jac_ok_ = false;

    // attempt 0: modified Newton (Jacobian kept while it converges fast), attempt 1: full Newton
    for(int attempt=0; attempt<2; ++attempt)
    {
        const bool full = set_.full_newton || attempt==1;
        if(attempt==1)
        {
            if(set_.full_newton) break;
            restore(m);
            jac_ok_ = false;
            ++stats.retries;
        }
        u.head(n) = dnun_;
        if(nc) u.tail(nc) = lam_/scale_;

        bool conv = false;
        double dprev = 0.0;
        for(int it=0; it<set_.max_iter; ++it)
        {
            residual(m, t1, h, u, r);
            if(!jac_ok_ || full)
                build_jacobian(m, t1, h, u, r, false, c);
            du = -lu_.solve(r);
            if(!du.allFinite())
                break;
            u += du;
            ++stats.iterations;
            const double dn = m.increment_norm(du.head(n));
            if(converged(m, dn, h))
            {
                conv = true;
                break;
            }
            // stagnation at the round-off level of the loads (e.g. finite-difference element forces):
            // the increments no longer fall but change the positions by less than 1e-9 of their size
            if(it >= 3 && dn > 0.5*dprev && stagnated(m, dn, h))
            {
                conv = true;
                break;
            }
            // slow convergence with an old Jacobian: rebuild it in the next iteration
            if(it>0 && dn > 0.25*dprev)
                jac_ok_ = false;
            dprev = dn;
        }
        if(!conv)
            continue;

        trial(m, u.head(n), h);
        const VecX dnu1 = u.head(n);
        an_ = ((1.0 - af_)*dnu1 + af_*dnun_ - am_*an_)/(1.0 - am_);
        dnun_ = dnu1;
        lam_ = scale_*u.tail(nc);
        for(auto &b : m.bodies)
            b.R = orthonormalize(b.R);
        store(m);
        m.commit(t1, h);
        ++stats.steps;
        if(nc)
        {
            m.constraints(g_, G_);
            stats.max_constraint = std::max(stats.max_constraint, g_.cwiseAbs().maxCoeff());
        }
        return true;
    }
    restore(m);
    jac_ok_ = false;
    ++stats.failures;
    return false;
}

VecX timeint::reactions(const model &m) const
{
    VecX g;
    MatX G;
    m.constraints(g, G);
    if(lam_.size()==0)
        return VecX::Zero(m.ndof());
    return -G.transpose()*lam_;
}

void timeint::save_state(std::vector<double> &v) const
{
    state_out o{v};
    o.d(scale_);
    o.list(xn_, [&](const Vec3 &x){o.m(x);});
    o.list(Rn_, [&](const Mat3 &x){o.m(x);});
    o.m(cqn_);
    o.m(nun_);
    o.m(dnun_);
    o.m(an_);
    o.m(lam_);
    o.b(jac_ok_);
    o.d(jac_h_);
    o.list(T_, [&](const Eigen::Triplet<double> &t){o.i(t.row()); o.i(t.col()); o.d(t.value());});
    for(long x : {stats.steps, stats.iterations, stats.jacobians, stats.retries, stats.failures})
        o.i(x);
    o.d(stats.max_constraint);
}

bool timeint::load_state(const model &m, const std::vector<double> &v, size_t &k)
{
    state_in in{v, k};
    in.d(scale_);
    in.list(xn_, [&](Vec3 &x){in.m(x);});
    in.list(Rn_, [&](Mat3 &x){in.m(x);});
    in.m(cqn_);
    in.m(nun_);
    in.m(dnun_);
    in.m(an_);
    in.m(lam_);
    in.b(jac_ok_);
    in.d(jac_h_);
    in.list(T_, [&](Eigen::Triplet<double> &t){const double r = in.next(), c = in.next(), x = in.next(); t = Eigen::Triplet<double>(int(r), int(c), x);});
    in.i(stats.steps); in.i(stats.iterations); in.i(stats.jacobians); in.i(stats.retries); in.i(stats.failures);
    in.d(stats.max_constraint);
    k = in.k;
    if(!in.ok)
        return false;
    // the factorization of the Jacobian in use (the same matrix gives the same factors)
    if(jac_ok_)
        jac_ok_ = lu_.compute(m.ndof() + m.ncon(), T_);
    return true;
}

}
