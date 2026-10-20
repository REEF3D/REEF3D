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

// REEF3D::CORAL: static equilibrium of the dynamic lumped-mass line (mooring_line::equilibrium), in its own
// file: next to the time integration in coral_mooring.cpp it changed the compiler's inlining of the
// sub-step loop (5 times slower).

#include"coral/mooring/coral_mooring.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>

namespace coral
{

bool mooring_line::equilibrium(const model &m, double t)
{
    const int N = nseg_;
    ra_ = m.body(ba_).point(sa_);
    rb_ = m.body(bb_).point(sb_);
    r_[0] = ra_;
    r_[N] = rb_;
    for(auto &x : v_) x.setZero();
    va_.setZero();
    vb_.setZero();
    ++gen_;
    // the water at rest at t (the surface of the nodes' start positions), the current acts
    if(waves_)
    {
        kinematics(0, t);
        kmode_ = 3;
    }
    else
        kmode_ = 0;
    const Vec3 z0 = Vec3::Zero();
    Vec3 Fa, Fb;
    state s{r_, v_}, st = s;
    std::vector<Vec3> a, ap;
    auto resid = [&](const state &x, std::vector<Vec3> &out)
    {
        accel(x, ra_, rb_, z0, z0, t, out, Fa, Fb, 0.0);
        double n = 0.0;
        for(int i=1; i<N; ++i)
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
        for(int i=1; i<N; ++i)
            n = std::max(n, out[i].norm());
        return n;
    };
    // nodes on the seabed start at the depth the penalty lets them sink (contact active: Newton sees its
    // stiffness instead of the kink at the seabed)
    if(pr_.kb > 0.0 && zbed_ > -1.0e29)
    {
        const double wsub = (pr_.mass - pr_.rho*0.25*M_PI*pr_.d*pr_.d)*m.gravity.norm();
        const double dz = std::max(wsub, 0.0)/(pr_.kb*pr_.d);
        for(int i=1; i<N; ++i)
            if(s.r[i](2) <= zbed_ + 1.0e-9*std::max(1.0, std::fabs(zbed_)))
                s.r[i](2) = zbed_ - dz;
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
                    if(j%3==c) st.r[j](d) += eps;
                resid(st, ap);
                for(int i=1; i<N; ++i)
                {
                    const Vec3 col = (ap[i] - a[i])/eps;
                    if((i-1)%3==c && i-1 >= 1) L[i].col(d) = col;
                    if(i%3==c) D[i].col(d) = col;
                    if((i+1)%3==c && i+1 <= N-1) U[i].col(d) = col;
                }
            }
        // scale of the Jacobian: its largest diagonal entry, at least the axial stiffness per node mass
        double dmax = pr_.mass > 0.0 ? pr_.EA/(l_*pr_.mass*l_) : 1.0;
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
                Vec3 bi = -a[i];
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
                st.r[i] += dx[i];
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
    // not converged (e.g. a straight line without tension, no transverse stiffness): keep the start shape
    if(ok)
        r_ = s.r;
    for(auto &x : v_) x.setZero();
    ++gen_;
    cvalid_ = false;
    return ok;
}

}
