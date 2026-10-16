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
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include"coral/core/coral_fdjac.h"
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
}

void timeint::store(const model &m)
{
    const int nb = m.nbody();
    xn_.resize(nb);
    Rn_.resize(nb);
    nun_.resize(6*nb);
    for(int i=0; i<nb; ++i)
    {
        const rigid_body &b = m.bodies[i];
        xn_[i] = b.x;
        Rn_[i] = b.R;
        nun_.segment<3>(6*i) = b.v;
        nun_.segment<3>(6*i+3) = b.w;
    }
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
}

void timeint::residual(model &m, double t1, double h, const VecX &u, VecX &r)
{
    const int n = m.ndof(), nc = (int)u.size() - n;
    trial(m, u.head(n), h);
    m.forces(t1, f_);
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
void timeint::residual_initial(model &m, double t0, const VecX &c, const VecX &u, VecX &r)
{
    const int n = m.ndof(), nc = (int)u.size() - n;
    for(int i=0; i<m.nbody(); ++i)
    {
        m.bodies[i].vd = u.segment<3>(6*i);
        m.bodies[i].wd = u.segment<3>(6*i+3);
    }
    m.forces(t0, f_);
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
    const int n = m.ndof(), N = (int)u.size(), nc = N - n;
    J_.setZero(N, N);
    VecX up = u, rp;
    fd_jacobian(m, J_, r,
        [&](int col, double e) {up(col) = u(col) + e;},
        [&](VecX &rr) {if(initial) residual_initial(m, t1, c, up, rr); else residual(m, t1, h, up, rr);},
        [&](int col) {return set_.fd_eps*std::max(1.0, std::fabs(u(col)));},
        false);
    (void)n;
    // lambda columns are exact; restore the work arrays of the unperturbed state
    if(initial)
        residual_initial(m, t1, c, u, rp);
    else
        residual(m, t1, h, u, rp);
    if(nc)
        J_.block(0, n, n, nc) = scale_*G_.transpose();
    lu_.compute(J_);
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
            const double dn = du.head(n).cwiseAbs().maxCoeff();
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

}
