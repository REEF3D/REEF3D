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


// REEF3D::CORAL: seabed friction of the dynamic mooring lines and points (mooring_system::seabed_friction,
// stick_friction), in its own file: next to the sub-step loop in coral_mooring.cpp, or in the larger set-up file,
// the compiler's inlining of the small vector operations changes (several times slower).

#include"coral/mooring/coral_mooring.h"
#include<cmath>
#include<limits>

namespace coral
{

// seabed friction of a node of line L in contact (normal force Fn > 0 along n, tangent q of the line): mu across
// and mu_axial along the line (elliptic limit). Stick-slip (stick_friction) or, with v_eps > 0, the regularized
// Coulomb law F = -Fn tanh(v / v_eps) (mua^2 va ea + mu^2 vt et) / |(mua va, mu vt)| (the isotropic law for
// mua = mu). has: a force (else none); contact: the node has a stick-slip anchor
Vec3 mooring_system::seabed_friction(const line_t &L, int i, const Vec3 &x, const Vec3 &vx, const Vec3 &q, const Vec3 &nb,
                                     double Fn, double ll, bool upd, bool &has, bool &contact) const
{
    const mooring_props &p = L.p;
    const double mua = p.mu_axial >= 0.0 ? p.mu_axial : p.mu;
    has = false;
    if(!(Fn > 0.0) || !(p.mu > 0.0 || mua > 0.0))
        return Vec3::Zero();
    if(p.v_eps > 0.0)
    {
        const Vec3 vh = vx - vx.dot(nb)*nb;
        const double vn = vh.norm();
        if(!(vn > 0.0))
            return Vec3::Zero();
        Vec3 e = q - q.dot(nb)*nb;
        const double en = e.norm();
        e = en > 1.0e-6 ? Vec3(e/en) : Vec3(vh/vn);
        const Vec3 va = vh.dot(e)*e, vt = vh - va;
        const double wn = std::sqrt(mua*mua*va.squaredNorm() + p.mu*p.mu*vt.squaredNorm());
        if(!(wn > 0.0))
            return Vec3::Zero();
        has = true;
        return -((Fn*std::tanh(vn/p.v_eps)/wn)*(mua*mua*va + p.mu*p.mu*vt));
    }
    if(!sc_)
        return Vec3::Zero();
    Vec3 e = q - q.dot(nb)*nb;
    const double en = e.norm();
    e = en > 1.0e-6 ? Vec3(e/en) : Vec3(Vec3::Zero());
    has = contact = true;
    return stick_friction(i, x, vx, nb, e, Fn, p.d*ll, p.ks >= 0.0 ? p.ks : p.kb, p.cs >= 0.0 ? p.cs : p.cb, p.mu, mua, p.mu_static, upd);
}

void mooring_system::lose_anchor(int i) const
{
    if(sc_ && !std::isnan(sc_->a[i](0)))
    {
        sc_->a[i] = Vec3::Constant(std::numeric_limits<double>::quiet_NaN());
        sc_->slip[i] = 0;
    }
}

// stick-slip seabed friction (mooring_props::mu): the tangential force -ks A d - cs A v of the shear spring from the
// node's anchor (sc_; d: the tangential distance to it) while it lies within the elliptic limit
// |(Fa/mua, Ft/mu)| <= mus Fn (sticking) or Fn (sliding), else on the kinetic limit (its closest point). upd (the state ends a
// sub-step): a node without an anchor gets one where it is; a node beyond the limit slides, its anchor follows so
// that the spring holds the kinetic limit; a sliding node sticks again once the spring force falls within it
Vec3 mooring_system::stick_friction(int i, const Vec3 &x, const Vec3 &vx, const Vec3 &n, const Vec3 &e, double Fn, double A,
                                    double ks, double cs, double mu, double mua, double mus, bool upd) const
{
    Vec3 &an = sc_->a[i];
    signed char &sl = sc_->slip[i];
    const bool none = std::isnan(an(0));
    Vec3 d = none ? Vec3(Vec3::Zero()) : Vec3(x - an);
    d -= d.dot(n)*n;
    const Vec3 vt = vx - vx.dot(n)*n;
    Vec3 F = -ks*A*d - cs*A*vt;
    // along the line (e; zero: isotropic, all with mu) and across it; nothing in a direction without friction
    double fa = 0.0;
    Vec3 ft = F;
    if(e.squaredNorm() > 0.0)
    {
        fa = F.dot(e);
        ft = F - fa*e;
    }
    if(!(mua > 0.0)) fa = 0.0;
    if(!(mu > 0.0)) ft.setZero();
    F = fa*e + ft;
    const double ft2 = ft.squaredNorm();
    const double w = std::sqrt((mua > 0.0 ? fa*fa/(mua*mua) : 0.0) + (mu > 0.0 ? ft2/(mu*mu) : 0.0));
    const double lim = (sl && !none ? 1.0 : mus)*Fn;
    const bool slide = w > lim;
    if(slide)
    {
        if(mua > 0.0 && mu > 0.0 && mua != mu && fa != 0.0 && ft2 > 0.0)
        {
            // the closest point of the kinetic ellipse (associated sliding: the velocity along the normal of the
            // ellipse, as the Coulomb law of different coefficients along and across): F = (fa / (1 + lam/mua^2),
            // ft / (1 + lam/mu^2)), lam >= 0 from Newton on the ellipse equation (convex, decreasing)
            const double a2 = mua*mua, b2 = mu*mu, fa2 = fa*fa, F2 = Fn*Fn;
            double lam = 0.0;
            for(int it=0; it<60; ++it)
            {
                const double pa = 1.0 + lam/a2, pb = 1.0 + lam/b2;
                const double g = fa2/(a2*pa*pa) + ft2/(b2*pb*pb) - F2;
                if(g <= 1.0e-14*F2)
                    break;
                const double dg = -2.0*fa2/(a2*a2*pa*pa*pa) - 2.0*ft2/(b2*b2*pb*pb*pb);
                lam -= g/dg;
            }
            F = (fa/(1.0 + lam/a2))*e + ft/(1.0 + lam/b2);
        }
        else
            F *= Fn/w;      // isotropic (or one direction): onto the circle
    }
    if(upd)
    {
        if(none)
        {
            an = x;
            sl = 0;
        }
        else if(slide)
        {
            an = ks*A > 0.0 ? Vec3(x + F/(ks*A)) : x;
            sl = 1;
        }
        else
            sl = 0;
    }
    return F;
}

}
