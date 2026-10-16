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

#include"coral/core/coral_static.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include"coral/core/coral_fdjac.h"
#include"coral/core/coral_linsolve.h"
#include<algorithm>
#include<cmath>
#include<vector>

namespace coral
{

namespace
{
    using config = model::config;
    config save(const model &m) {return m.save_config();}
    void set(model &m, const config &c, const VecX &dq, double s) {m.set_config(c, dq, s);}

    // residual for the configuration in the model and the scaled multipliers mu
    void residual(model &m, double t, double scale, const VecX &mu, VecX &r, VecX &f, VecX &g, MatX &G, bool fd_only=false)
    {
        const int n = m.ndof(), nc = m.ncon();
        m.forces(t, f, fd_only);
        m.constraints(g, G);
        r.resize(n + nc);
        r.head(n) = -f;
        if(nc)
        {
            r.head(n) += scale*G.transpose()*mu;
            r.tail(nc) = scale*g;
        }
    }
}

// one Newton solve from the configuration in the model; damped: Deuflhard's natural monotonicity
// test on every step, else full Newton steps (stopped when the corrections blow up)
static bool newton(model &m, double t, const static_settings &s, bool damped, double scale, VecX &mu, static_result &res)
{
    const int n = m.ndof(), nc = m.ncon(), N = n + nc;
    VecX r, f, g;
    MatX G;
    std::vector<triplet> T;
    linear_solver lu;
    lu.dense_max = s.dense_max;
    config c = save(m);
    const bool tang = m.has_tangents();
    const VecX zero = VecX::Zero(n);
    residual(m, t, scale, mu, r, f, g, G);
    double d0 = -1.0;
    for(int it=0; it<s.max_iter; ++it)
    {
        // Jacobian: analytic tangents of the loads that have one, the configuration columns of the rest
        // by central differences, multiplier columns exact
        T.clear();
        VecX dq = zero;
        const int nb = m.nbody();
        fd_jacobian(m, T, r,
            [&](int col, double e) {dq(col) = e;},
            [&](VecX &rr) {set(m, c, dq, 1.0); VecX ff, gg; MatX GG; residual(m, t, scale, mu, rr, ff, gg, GG, tang);},
            [&](int col) {return s.fd_eps*(col < 6*nb ? ((col%6) < 3 ? std::max(1.0, std::fabs(c.x[col/6](col%6))) : 1.0)
                                                      : std::max(1.0, std::fabs(c.cq(col - 6*nb))));},
            true, tang);
        set(m, c, zero, 1.0);
        residual(m, t, scale, mu, r, f, g, G);
        if(tang)
            m.tangents(t, -1.0, 0.0, 0.0, T);
        for(int j=0; j<nc; ++j)
            for(int i=0; i<n; ++i)
                if(G(j, i) != 0.0)
                    T.emplace_back(i, n + j, scale*G(j, i));
        lu.compute(N, T);
        const VecX d = -lu.solve(r);
        if(!d.allFinite())
            return false;
        const double dn = d.norm();
        if(d0 < 0.0) d0 = dn;
        if(!damped && dn > 1.0e3*std::max(d0, 1.0e-12))
            return false;
        double a = 1.0;
        if(damped)
        {
            bool ok = false;
            for(int k=0; k<30 && !ok; ++k)
            {
                set(m, c, d.head(n), a);
                VecX r1;
                residual(m, t, scale, VecX(mu + a*d.tail(nc)), r1, f, g, G);
                const VecX db = lu.solve(r1);
                ok = db.allFinite() && db.norm() <= (1.0 - 0.25*a)*dn;
                if(!ok) a *= 0.5;
            }
            if(!ok)
                return false;
        }
        set(m, c, d.head(n), a);
        mu += a*d.tail(nc);
        for(auto &b : m.bodies) b.R = orthonormalize(b.R);
        c = save(m);
        residual(m, t, scale, mu, r, f, g, G);
        ++res.iterations;
        double xmax = 0.0;
        for(const auto &b : m.bodies) xmax = std::max(xmax, b.x.cwiseAbs().maxCoeff());
        if(a*m.increment_norm(d.head(n)) <= s.tol*(1.0 + xmax))
            return true;
    }
    return false;
}

static_result static_solve(model &m, double t, const static_settings &s)
{
    static_result res;
    const int n = m.ndof(), nc = m.ncon();
    m.zero_rates();
    double msum = 0.0;
    for(const auto &b : m.bodies) msum += b.mass;
    const double scale = std::max(1.0, m.nbody() ? msum/m.nbody() : 1.0);
    const config c0 = save(m);
    VecX mu = VecX::Zero(nc);
    // full Newton first (fast for structures: large axial forces after a first step do not matter),
    // the damped variant from the start configuration if it fails
    res.converged = newton(m, t, s, false, scale, mu, res);
    if(!res.converged)
    {
        set(m, c0, VecX::Zero(n), 1.0);
        mu.setZero();
        res.converged = newton(m, t, s, true, scale, mu, res);
    }
    VecX r, f, g;
    MatX G;
    residual(m, t, scale, mu, r, f, g, G);
    res.residual = n ? r.head(n).cwiseAbs().maxCoeff() : 0.0;
    return res;
}

}
